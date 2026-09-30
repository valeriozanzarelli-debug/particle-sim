//! Thermal environment: the black-body photons of a room at temperature T,
//! drawn as a faint haze around the bench.
//!
//! PURELY VISUAL for now. The haze does not interact with the photon under
//! study; coupling the environment to the physics (thermal background
//! clicks on the screen, decoherence of the grains' phases) comes later.
//!
//! What follows Planck:
//!   • count ∝ T³ — the number density of black-body photons,
//!     n = (2ζ(3)/π²)(k_B T / ħc)³, sampled as HAZE_MAX_POINTS · (T/T_MAX)³;
//!     so at room temperature only a handful of points are drawn;
//!   • colour from the black-body spectrum (peak λ by Wien, λ_max = b/T):
//!     below the Draper point (~800 K) the peak is deep infrared, drawn as a
//!     dim dark red; above it the approximate sRGB colour of a black body.
//! The haze photons fly at c (same speed as the scene's light) and follow
//! the time bar: they freeze when paused and fly backwards when rewinding.

use bevy::prelude::*;
use rand::{rngs::StdRng, Rng, SeedableRng};

use crate::config::*;
use crate::render::{cloud_mesh, point_mesh, swap_mesh, Cloud, PARKED};
use crate::timeline::base_speed;
use crate::LabState;

/// The box the haze fills (it wraps around at the faces).
const BOX_MIN: Vec3 = Vec3::new(-90.0, -50.0, -90.0);
const BOX_MAX: Vec3 = Vec3::new(90.0, 50.0, 190.0);

#[derive(Resource)]
pub struct Haze {
    /// (start position, unit direction) of each haze photon.
    seeds:  Vec<(Vec3, Vec3)>,
    /// Distance the haze photons have flown (signed: rewinding flies back).
    flown:  f32,
    drawn:  usize,
}

pub fn setup_haze(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let mut rng = StdRng::seed_from_u64(0x7EA7);
    let size = BOX_MAX - BOX_MIN;
    let seeds = (0..HAZE_MAX_POINTS).map(|_| {
        let p = BOX_MIN + Vec3::new(rng.gen::<f32>(), rng.gen::<f32>(), rng.gen::<f32>()) * size;
        let z: f32 = rng.gen_range(-1.0..1.0);
        let a: f32 = rng.gen_range(0.0..std::f32::consts::TAU);
        let r = (1.0 - z * z).sqrt();
        (p, Vec3::new(r * a.cos(), r * a.sin(), z))
    }).collect();
    let mesh = meshes.add(point_mesh(1));
    let mat = materials.add(StandardMaterial { base_color: Color::WHITE, unlit: true, ..default() });
    commands.spawn((Cloud::Haze, Mesh3d(mesh), MeshMaterial3d(mat), Transform::default()));
    commands.insert_resource(Haze { seeds, flown: 0.0, drawn: usize::MAX });
}

/// Number of haze points drawn at temperature `t` (∝ T³).
pub fn haze_count(t: f32) -> usize {
    let x = (t / T_MAX).clamp(0.0, 1.0);
    (HAZE_MAX_POINTS as f32 * x * x * x).round() as usize
}

/// Wien peak wavelength in µm (None at 0 K).
pub fn wien_peak_um(t: f32) -> Option<f32> {
    (t > 0.5).then(|| WIEN_UM_K / t)
}

/// Approximate sRGB colour and brightness of a black body at `t` kelvin.
fn thermal_color(t: f32) -> [f32; 4] {
    // Tanner Helland's fit of the black-body locus (valid ≥ 1000 K).
    let x = t.max(1000.0) / 100.0;
    let r = if x <= 66.0 { 255.0 } else { 329.699 * (x - 60.0).powf(-0.133_205) };
    let g = if x <= 66.0 { 99.4708 * x.ln() - 161.1196 } else { 288.122 * (x - 60.0).powf(-0.075_515) };
    let b = if x >= 66.0 { 255.0 } else if x <= 19.0 { 0.0 } else { 138.518 * (x - 10.0).ln() - 305.045 };
    // Invisible (infrared) below the Draper point: dim; brighter as it glows.
    let s = ((t - 700.0) / (T_MAX - 700.0)).clamp(0.0, 1.0);
    let bright = 0.22 + 0.6 * s * s * (3.0 - 2.0 * s);
    let c = |v: f32| (v / 255.0).clamp(0.0, 1.0) * bright;
    [c(r), c(g), c(b), 1.0]
}

pub fn update_haze(
    time: Res<Time>,
    state: Res<LabState>,
    haze: Option<ResMut<Haze>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut clouds: Query<(&Cloud, &mut Mesh3d)>,
) {
    let Some(mut haze) = haze else { return };
    let n = haze_count(state.temperature);
    let v = base_speed(&state.sim.lab) * state.time.signed_factor();
    let moved = v != 0.0 && n > 0;
    if moved {
        let d = v * time.delta_secs().min(MAX_DT);
        haze.flown = (haze.flown + d).rem_euclid(1.0e5);
    }
    if !moved && n == haze.drawn { return; }
    haze.drawn = n;

    let size = BOX_MAX - BOX_MIN;
    let col = thermal_color(state.temperature);
    let mut pos = Vec::with_capacity(n.max(1));
    for &(p0, dir) in &haze.seeds[..n] {
        let p = p0 + dir * haze.flown - BOX_MIN;
        let w = Vec3::new(p.x.rem_euclid(size.x), p.y.rem_euclid(size.y), p.z.rem_euclid(size.z));
        pos.push((BOX_MIN + w).to_array());
    }
    if pos.is_empty() { pos.push(PARKED); }
    let cols = vec![col; pos.len()];
    swap_mesh(Cloud::Haze, cloud_mesh(pos, cols), &mut meshes, &mut clouds);
}

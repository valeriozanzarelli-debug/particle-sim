//! Drawing the photon (classic and 2-D pilot rules): grains as points
//! (colour = phase), E-field arrows, optical elements. The base 3-D rule,
//! the clicks (round dots) and the photon ball of the 3-D rule are drawn by
//! `render3d.rs`.

use bevy::{
    prelude::*,
    render::{
        mesh::PrimitiveTopology,
        render_asset::RenderAssetUsages,
    },
};

use crate::config::*;
use crate::physics::{Experiment, TAU};
use crate::LabState;

/// Which point cloud an entity draws.
#[derive(Component, Clone, Copy, PartialEq, Eq)]
pub enum Cloud { Grains, Haze }

/// Give the cloud a new mesh and remove the old one.
///
/// Bevy 0.15 never frees the GPU allocation of a mesh that is *modified*
/// (only of one that is *removed*), so editing a mesh in place every frame
/// leaks its whole size each time. Swapping the asset frees the old one.
pub(crate) fn swap_mesh(
    which: Cloud,
    mesh: Mesh,
    meshes: &mut Assets<Mesh>,
    clouds: &mut Query<(&Cloud, &mut Mesh3d)>,
) {
    for (cloud, mut mesh3d) in clouds.iter_mut() {
        if *cloud == which {
            let old = mesh3d.0.id();
            mesh3d.0 = meshes.add(mesh);
            meshes.remove(old);
            return;
        }
    }
}

/// A point cloud with the given positions and colours.
pub(crate) fn cloud_mesh(pos: Vec<[f32; 3]>, col: Vec<[f32; 4]>) -> Mesh {
    let n = pos.len();
    let mut mesh = Mesh::new(PrimitiveTopology::PointList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR,    col);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL,   vec![[0.0f32, 1.0, 0.0]; n]);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0,     vec![[0.0f32; 2]; n]);
    mesh
}

/// What the clouds last drew, to skip the upload when nothing changed.
#[derive(Resource, Default)]
pub struct Drawn {
    grains: Option<(u64, f32, usize)>,
}

/// Where dead grains are parked (outside any view).
pub(crate) const PARKED: [f32; 3] = [1.0e7, 1.0e7, 1.0e7];

pub(crate) fn point_mesh(n: usize) -> Mesh {
    cloud_mesh(vec![PARKED; n], vec![[0.0; 4]; n])
}

pub fn setup_render(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let mat = materials.add(StandardMaterial { base_color: Color::WHITE, unlit: true, ..default() });
    commands.spawn((Cloud::Grains, Mesh3d(meshes.add(point_mesh(1))), MeshMaterial3d(mat), Transform::default()));
    commands.init_resource::<Drawn>();
}

/// Phase colour wheel: the colour of a grain is the angle of its arrow.
pub(crate) fn phase_color(angle: f32) -> [f32; 4] {
    let t = angle.rem_euclid(TAU) / TAU;
    // Cyclic palette: cyan â†’ blue â†’ violet â†’ magenta â†’ gold â†’ cyan
    const STOPS: [[f32; 3]; 5] = [
        [0.20, 0.95, 1.00],
        [0.25, 0.45, 1.00],
        [0.70, 0.30, 1.00],
        [1.00, 0.35, 0.70],
        [1.00, 0.85, 0.30],
    ];
    let x = t * 5.0;
    let i = (x.floor() as usize) % 5;
    let s = x - x.floor();
    let (a, b) = (STOPS[i], STOPS[(i + 1) % 5]);
    [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s, 1.0]
}

pub fn update_meshes(
    state: Res<LabState>,
    drawn: Option<ResMut<Drawn>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut clouds: Query<(&Cloud, &mut Mesh3d)>,
) {
    let Some(mut drawn) = drawn else { return };
    let sim = &state.sim;
    let k = sim.lab.k();

    // â”€â”€ grains â”€â”€
    // All grains are simulated; at most MAX_GRAINS_DRAWN are drawn (every
    // `stride`-th one), because the cloud is uploaded again every frame.
    let total: usize = sim.photons.iter().map(|p| p.grains.len()).sum();
    let first_id = sim.photons.first().map_or(u64::MAX, |p| p.id);
    let key = (first_id, sim.clock, total);
    if drawn.grains != Some(key) {
        drawn.grains = Some(key);
        let stride = total.div_ceil(MAX_GRAINS_DRAWN).max(1);
        let n = total.div_ceil(stride).max(1);
        let mut pos: Vec<[f32; 3]> = Vec::with_capacity(n);
        let mut col: Vec<[f32; 4]> = Vec::with_capacity(n);
        for ph in &sim.photons {
            for g in ph.grains.iter().step_by(stride).filter(|g| g.alive) {
                pos.push(g.pos.to_array());
                let mut c = phase_color(g.phase(k) - k * sim.clock);
                if !ph.pol_defined {
                    // Entangled, unmeasured: phase yes, polarization none â†’ pale
                    for v in c.iter_mut().take(3) { *v = 0.55 + 0.45 * *v * 0.5; }
                }
                col.push(c);
            }
        }
        if pos.is_empty() {
            pos.push(PARKED);
            col.push([0.0; 4]);
        }
        swap_mesh(Cloud::Grains, cloud_mesh(pos, col), &mut meshes, &mut clouds);
    }
}

/// Optical bench and E-field arrows.
pub fn draw_bench(state: Res<LabState>, vis: Option<Res<crate::render3d::Visuals>>, mut gizmos: Gizmos) {
    let sim = &state.sim;
    let lab = &sim.lab;
    let dim   = Color::srgba(0.55, 0.6, 0.7, 0.6);
    let gold  = Color::srgb(1.0, 0.8, 0.3);
    let green = Color::srgb(0.6, 1.0, 0.5);

    let rect = |g: &mut Gizmos, z: f32, x0: f32, x1: f32, y0: f32, y1: f32, c: Color| {
        let p = [Vec3::new(x0, y0, z), Vec3::new(x1, y0, z), Vec3::new(x1, y1, z), Vec3::new(x0, y1, z)];
        for i in 0..4 { g.line(p[i], p[(i + 1) % 4], c); }
    };
    let polarizer = |g: &mut Gizmos, z: f32, angle: f32, c: Color| {
        let r = 10.0;
        let iso = Isometry3d::new(Vec3::new(0.0, 0.0, z), Quat::IDENTITY);
        g.circle(iso, r, c);
        let a = Vec3::new(angle.cos(), angle.sin(), 0.0) * r;
        g.line(Vec3::new(0.0, 0.0, z) - a, Vec3::new(0.0, 0.0, z) + a, c);
    };

    // Source
    if lab.experiment.only_3d() {
        // the source is the lit aperture / the beam around the disk
    } else if lab.experiment != Experiment::Bell {
        let zs = if lab.base3d() { lab.source_z } else { lab.source_z - 12.0 };
        gizmos.circle(Isometry3d::new(Vec3::new(0.0, 0.0, zs), Quat::IDENTITY), if lab.base3d() { lab.beam_radius(lab.source_z) } else { 2.0 }, gold);
    } else {
        gizmos.sphere(Isometry3d::from_translation(Vec3::ZERO), 1.5, gold);
    }

    match lab.experiment {
        Experiment::FreeBeam => {}
        Experiment::DoubleSlit if lab.base3d() => {
            // 3-D: slits of finite length; outline of the wall and of each opening
            let (c, h) = lab.slits();
            let z = lab.wall_z;
            let hy = 0.5 * lab.slit_height;
            rect(&mut gizmos, z, -40.0, 40.0, -30.0, 30.0, dim);
            for (s, f) in lab.slit_filters().iter().enumerate() {
                rect(&mut gizmos, z, c[s] - h, c[s] + h, -hy, hy, gold.with_alpha(0.8));
                if let Some(f) = f {
                    let a = f[1].re.atan2(f[0].re);
                    let d = Vec3::new(a.cos(), a.sin(), 0.0) * 3.0;
                    let p = Vec3::new(c[s], 0.0, z - 0.5);
                    gizmos.line(p - d, p + d, gold);
                }
            }
        }
        Experiment::DoubleSlit => {
            let (c, h) = lab.slits();
            let (x0, x1, yh) = (-40.0, 40.0, 25.0);
            let z = lab.wall_z;
            // Wall pieces between/around the slits
            rect(&mut gizmos, z, x0, c[0] - h, -yh, yh, dim);
            rect(&mut gizmos, z, c[0] + h, c[1] - h, -yh, yh, dim);
            rect(&mut gizmos, z, c[1] + h, x1, -yh, yh, dim);
            // Which-path polarizers on the slits
            for (s, f) in lab.slit_filters().iter().enumerate() {
                if let Some(f) = f {
                    let a = f[1].re.atan2(f[0].re);
                    let d = Vec3::new(a.cos(), a.sin(), 0.0) * 3.0;
                    let p = Vec3::new(c[s], 0.0, z - 0.5);
                    gizmos.line(p - d, p + d, gold);
                }
            }
        }
        Experiment::Polarizers => {
            for p in lab.polarizers.iter().filter(|p| p.enabled) {
                polarizer(&mut gizmos, p.z, p.angle, gold);
            }
        }
        Experiment::Airy => {
            // opaque screen with a circular hole (the hole is the source)
            let iso = Isometry3d::new(Vec3::new(0.0, 0.0, lab.wall_z), Quat::IDENTITY);
            gizmos.circle(iso, 0.5 * lab.aperture_d, gold);
            rect(&mut gizmos, lab.wall_z, -20.0, 20.0, -20.0, 20.0, dim);
        }
        Experiment::Arago => {
            let iso = Isometry3d::new(Vec3::new(0.0, 0.0, lab.wall_z), Quat::IDENTITY);
            gizmos.circle(iso, lab.disk_r, gold);
            gizmos.circle(iso, 0.6 * lab.disk_r, gold.with_alpha(0.5));
            gizmos.circle(iso, 0.25 * lab.disk_r, gold.with_alpha(0.35));
        }
        Experiment::Bell => {
            let last = sim.stats.last_bell;
            for (side, d) in [(0usize, lab.bell.dist_a), (1, lab.bell.dist_b)] {
                let z = if side == 0 { d } else { -d };
                let angle = match (last, lab.bell.random_chsh) {
                    (Some((a, _, b, _)), true) => if side == 0 { a } else { b },
                    _ => if side == 0 { lab.bell.angle_a } else { lab.bell.angle_b },
                };
                let c = if side == 0 { Color::srgb(0.4, 0.9, 1.0) } else { Color::srgb(1.0, 0.6, 0.3) };
                polarizer(&mut gizmos, z, angle, c);
                let zd = if side == 0 { d + 15.0 } else { -d - 15.0 };
                rect(&mut gizmos, zd, -6.0, 6.0, -6.0, 6.0, dim);
            }
        }
    }

    // Detection screen
    if lab.experiment != Experiment::Bell {
        let yh = if lab.base3d() { lab.screen_half } else { 25.0 };
        rect(&mut gizmos, lab.screen_z, -lab.screen_half, lab.screen_half, -yh, yh, green.with_alpha(0.5));
    }

    // Base 3-D rule: arrows on a few drawn wave balls
    if state.show_arrows && lab.base3d() {
        if let Some(v) = vis.as_deref() {
            crate::render3d::arrows(v, sim.clock, crate::render3d::vis_scale(lab), &mut gizmos);
        }
    }

    // E-field arrows on a subset of grains: rotate + translate made visible
    if state.show_arrows {
        let k = lab.k();
        // pilot rule: dimmer, bluish arrows, so the white photon ball stands out
        let arrow = if lab.pilot() { Color::srgba(0.6, 0.75, 1.0, 0.3) } else { Color::srgba(1.0, 1.0, 1.0, 0.55) };
        for ph in sim.photons.iter().filter(|p| p.pol_defined) {
            let alive: Vec<_> = ph.grains.iter().filter(|g| g.alive).collect();
            let stride = (alive.len() / MAX_ARROWS).max(1);
            for g in alive.iter().step_by(stride) {
                let e = g.e_field(k, sim.clock) * ARROW_LEN;
                gizmos.line(g.pos, g.pos + e, arrow);
            }
        }
    }
}

/// Gizmos of the photon ball: drawn on top of the wave balls, thicker.
#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct PhotonBallGizmos;

pub fn photon_ball_gizmo_config() -> GizmoConfig {
    GizmoConfig { line_width: 3.0, depth_bias: -1.0, ..default() }
}

/// Pilot rule: the photon ball (white sphere + star) and its recent path.
/// Base 3-D rule: the trail of the photon ball (its sphere and glow are
/// drawn by `render3d.rs`).
pub fn draw_photon_ball(state: Res<LabState>, mut gizmos: Gizmos<PhotonBallGizmos>) {
    let white = Color::srgb(1.0, 1.0, 1.0);
    let clock = state.sim.clock as f64;
    for b in state.sim.photons.iter().filter_map(|p| p.ball.as_ref()) {
        let pts: Vec<Vec3> = b.trail.iter().filter(|q| q.0 <= clock)
            .map(|q| Vec3::new(q.1[0] as f32, q.1[1] as f32, q.1[2] as f32))
            .chain(std::iter::once(b.pos_at(clock))).collect();
        let n = pts.len().max(1) as f32;
        gizmos.linestrip_gradient(pts.into_iter().enumerate().map(|(i, p)| {
            let a = 0.15 + 0.85 * (i as f32 + 1.0) / n;
            (p, Color::srgba(1.0, 0.95, 0.8, a))
        }));
    }
    for pp in state.sim.photons.iter().filter_map(|p| p.pilot.as_ref()) {
        let n = pp.trail.len().max(1) as f32;
        gizmos.linestrip_gradient(pp.trail.iter().enumerate().map(|(i, p)| {
            let a = 0.2 + 0.8 * (i as f32 + 1.0) / n;
            (*p, Color::srgba(1.0, 1.0, 1.0, a))
        }));
        let c = pp.pos();
        gizmos.sphere(Isometry3d::from_translation(c), 1.6, white);
        gizmos.sphere(Isometry3d::from_translation(c), 0.6, white);
        for d in [Vec3::X, Vec3::Y, Vec3::Z] {
            gizmos.line(c - d * 4.5, c + d * 4.5, white);
        }
    }
}

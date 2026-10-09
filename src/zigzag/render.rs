//! Drawing the walks as clouds of round soft balls (own small shader, `balls.wgsl`).
//!
//! • light: one ball per lattice site where the electric energy (Re F)² is large: the crests;
//! • electron: one ball per site, colour = handedness (orange = right-handed R, blue = left-handed L),
//!   brightness = density; the surfers and ⟨r⟩ (white) leave trails: the Zitterbewegung loops;
//! • hydrogen: balls drawn with |ψ|² of the chosen state, colour = phase of the wave (it turns in time
//!   with the level's frequency, slowed down), red nucleus, surfers with trails.
//! The live walks change every frame: the mesh asset is REPLACED (Bevy 0.15 leaks the GPU copy of a
//! mesh modified in place; see photon/render.rs swap_mesh).

use bevy::{
    asset::load_internal_asset,
    pbr::{MaterialPipeline, MaterialPipelineKey, NotShadowCaster, NotShadowReceiver},
    prelude::*,
    render::{
        mesh::{Indices, MeshVertexAttribute, MeshVertexBufferLayoutRef, PrimitiveTopology},
        render_asset::RenderAssetUsages,
        render_resource::{AsBindGroup, RenderPipelineDescriptor, ShaderRef, SpecializedMeshPipelineError, VertexFormat},
        view::NoFrustumCulling,
    },
};
use rayon::prelude::*;

use crate::atom::{self, Level, ALPHA};
use crate::config::*;
use crate::{scene, Lab, Modo, Rng, STATI};

pub const BALL_SHADER: Handle<Shader> = Handle::weak_from_u128(0x7a1c_33e9_0b54_4c2d_9e61_5f08_c2b4_11d7);
pub const ATTR_PHASE: MeshVertexAttribute = MeshVertexAttribute::new("Zz_Phase", 0x7a1c_0001, VertexFormat::Float32x4);

#[allow(dead_code)]
mod params {
    use bevy::render::render_resource::ShaderType;
    #[derive(Clone, Copy, Debug, ShaderType)]
    pub struct BallParams {
        pub clock_phase: f32,
        pub size: f32,
        pub gain: f32,
        pub soft: f32,
    }
}
pub use params::BallParams;

#[derive(Asset, TypePath, AsBindGroup, Clone, Debug)]
pub struct BallMaterial {
    #[uniform(0)]
    pub params: BallParams,
}

impl Material for BallMaterial {
    fn vertex_shader() -> ShaderRef { BALL_SHADER.into() }
    fn fragment_shader() -> ShaderRef { BALL_SHADER.into() }
    fn alpha_mode(&self) -> AlphaMode { AlphaMode::Add }
    fn specialize(
        _pipeline: &MaterialPipeline<Self>,
        descriptor: &mut RenderPipelineDescriptor,
        layout: &MeshVertexBufferLayoutRef,
        _key: MaterialPipelineKey<Self>,
    ) -> Result<(), SpecializedMeshPipelineError> {
        let vl = layout.0.get_layout(&[
            Mesh::ATTRIBUTE_POSITION.at_shader_location(0),
            ATTR_PHASE.at_shader_location(1),
            Mesh::ATTRIBUTE_COLOR.at_shader_location(2),
            Mesh::ATTRIBUTE_UV_0.at_shader_location(3),
        ])?;
        descriptor.vertex.buffers = vec![vl];
        descriptor.primitive.cull_mode = None;
        Ok(())
    }
}

pub struct BallsPlugin;

impl Plugin for BallsPlugin {
    fn build(&self, app: &mut App) {
        load_internal_asset!(app, BALL_SHADER, "balls.wgsl", Shader::from_wgsl);
        app.add_plugins(MaterialPlugin::<BallMaterial> { prepass_enabled: false, shadows_enabled: false, ..default() })
            .add_systems(Startup, setup);
    }
}

/// The entity drawing all the balls.
#[derive(Component)]
pub struct Nuvola;

#[derive(Resource)]
pub struct Visuals {
    mat: Handle<BallMaterial>,
    /// what the current mesh shows: (mode, frame or cloud version)
    key: (u8, u64),
}

/// One ball: position, phase data (phase, use-phase flag, size factor, 0), colour.
pub type Ball = ([f32; 3], [f32; 4], [f32; 4]);

pub fn quads(items: &[Ball]) -> Mesh {
    let hidden = [([0.0f32, -1.0e6, 0.0], [0.0f32; 4], [0.0f32; 4])];
    let items = if items.is_empty() { &hidden[..] } else { items };
    let n = items.len();
    let mut pos = Vec::with_capacity(4 * n);
    let mut ph = Vec::with_capacity(4 * n);
    let mut col = Vec::with_capacity(4 * n);
    let mut uv = Vec::with_capacity(4 * n);
    let mut idx: Vec<u32> = Vec::with_capacity(6 * n);
    const CORNERS: [[f32; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    for (q, it) in items.iter().enumerate() {
        for c in CORNERS {
            pos.push(it.0);
            ph.push(it.1);
            col.push(it.2);
            uv.push(c);
        }
        let b = 4 * q as u32;
        idx.extend_from_slice(&[b, b + 1, b + 2, b, b + 2, b + 3]);
    }
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    mesh.insert_attribute(ATTR_PHASE, ph);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, col);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uv);
    mesh.insert_indices(Indices::U32(idx));
    mesh
}

fn setup(mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>, mut mats: ResMut<Assets<BallMaterial>>) {
    let mat = mats.add(BallMaterial { params: BallParams { clock_phase: 0.0, size: 0.6, gain: 1.0, soft: 3.0 } });
    commands.spawn((Nuvola, Mesh3d(meshes.add(quads(&[]))), MeshMaterial3d(mat.clone()), Transform::default(),
                    NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
    commands.insert_resource(Visuals { mat, key: (255, u64::MAX) });
}

const C_R: [f32; 3] = [1.0, 0.55, 0.12];
const C_L: [f32; 3] = [0.2, 0.55, 1.0];
const GOLD: [f32; 3] = [1.0, 0.86, 0.45];

/// Light: balls on the sites where the electric energy (Re F)² is near a crest; the near half of the
/// box (physics y < −1) is cut away so the crests inside show.
fn balls_light(lab: &Lab) -> Vec<Ball> {
    let l = &lab.light;
    let n = l.n;
    let e: Vec<f32> = l.f.par_iter().map(|p| p.iter().map(|c| (c.re * c.re) as f32).sum()).collect();
    let emax = e.iter().cloned().fold(0.0, f32::max).max(1e-30);
    let thr = 0.12 * emax;
    let half = n as f64 / 2.0 - 0.5;
    e.iter().enumerate().filter(|(s, &v)| v > thr && ((s / n) % n) as f64 - half >= -1.0).map(|(s, &v)| {
        let p = scene([(s / (n * n)) as f64 - half, ((s / n) % n) as f64 - half, (s % n) as f64 - half]);
        let a = 0.30 * (v / emax).powf(0.8);
        ([p.x, p.y, p.z], [0.0, 0.0, 1.0, 0.0], [GOLD[0], GOLD[1], GOLD[2], a])
    }).take(MAX_BALLS).collect()
}

/// Electron: balls coloured by handedness R (orange) / L (blue), brightness ∝ density. One ball every
/// second site per axis (the cloud stays readable).
fn balls_electron(lab: &Lab) -> Vec<Ball> {
    let d = &lab.dirac;
    let n = d.n;
    let rl: Vec<(f32, f32)> = d.psi.par_iter().map(|p| ((p[0].norm2() + p[1].norm2()) as f32, (p[2].norm2() + p[3].norm2()) as f32)).collect();
    let rmax = rl.iter().map(|v| v.0 + v.1).fold(0.0, f32::max).max(1e-30);
    let thr = 0.06 * rmax;
    let half = n as f64 / 2.0 - 0.5;
    rl.iter().enumerate().filter(|(s, v)| {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        // every second site, and the near half cut away (physics y < 0) to see inside
        i % 2 == 0 && j % 2 == 0 && k % 2 == 0 && v.0 + v.1 > thr && (j as f64) >= half
    }).map(|(s, &(l, r))| {
        let p = scene([(s / (n * n)) as f64 - half, ((s / n) % n) as f64 - half, (s % n) as f64 - half]);
        let h = r / (l + r);
        let c = [C_L[0] + (C_R[0] - C_L[0]) * h, C_L[1] + (C_R[1] - C_L[1]) * h, C_L[2] + (C_R[2] - C_L[2]) * h];
        let a = 0.20 * ((l + r) / rmax).powf(0.5);
        ([p.x, p.y, p.z], [0.0, 0.0, 1.0, 0.0], [c[0], c[1], c[2], a])
    }).take(MAX_BALLS).collect()
}

/// Points drawn with |ψ|² of a hydrogen state (positions in λ_C) and the phase used for the colour.
pub fn sample_cloud(lv: &Level, kappa: i32, two_m: i32, n: usize, rng: &mut Rng) -> Vec<([f64; 3], f32)> {
    // radial distribution g² + f² on the grid → inverse cdf
    let w: Vec<f64> = lv.gf.iter().map(|v| v[0] * v[0] + v[1] * v[1]).collect();
    let mut cdf = Vec::with_capacity(w.len());
    let mut s = 0.0;
    for x in &w { s += x; cdf.push(s); }
    let mut out = Vec::with_capacity(n);
    let mut guard = 0usize;
    while out.len() < n && guard < 400 * n {
        guard += 1;
        let u = rng.uni() * s;
        let i = cdf.partition_point(|&c| c < u).min(w.len() - 1);
        let r = (i as f64 + rng.uni()) * lv.a;
        // uniform direction, accepted with the angular density of |ψ|² at this radius
        let ct = 2.0 * rng.uni() - 1.0;
        let ph = 2.0 * std::f64::consts::PI * rng.uni();
        let st = (1.0 - ct * ct).max(0.0).sqrt();
        let p = [r * st * ph.cos(), r * st * ph.sin(), r * ct];
        let sp = atom::spinor(lv, kappa, two_m, p);
        let dens: f64 = sp.iter().map(|c| c.norm2()).sum::<f64>() * r * r;
        let rad = (w[i] / lv.a.max(1e-300)).max(1e-300);
        let ang = dens / (rad * lv.a) * 4.0 * std::f64::consts::PI;
        if rng.uni() * 3.2 < ang {
            let phase = if sp[0].norm2() >= sp[1].norm2() { sp[0].arg() } else { sp[1].arg() };
            out.push((p, phase as f32));
        }
    }
    out
}

fn balls_atom(lab: &mut Lab) -> Vec<Ball> {
    let guard = lab.hydro.lock().unwrap();
    let Some(h) = guard.as_ref() else { return vec![]; };
    let st = lab.h_stato;
    let lv = &h.states[st];
    let mut rng = Rng::new(77 + st as u64);
    let pts = sample_cloud(lv, STATI[st].1, STATI[st].2, N_CLOUD, &mut rng);
    let sc = (ALPHA * A0_SCENE as f64) as f32;
    pts.iter().map(|(p, ph)| {
        let q = scene(*p) * sc;
        ([q.x, q.y, q.z], [*ph, 1.0, 0.55, 0.0], [1.0, 1.0, 1.0, 0.16])
    }).collect()
}

/// Rebuild the balls (every frame for the live walks, once per state for the atom), set the clock phase.
pub fn aggiorna(
    mut lab: ResMut<Lab>,
    mut vis: ResMut<Visuals>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut mats: ResMut<Assets<BallMaterial>>,
    mut q: Query<&mut Mesh3d, With<Nuvola>>,
) {
    let have_atom = lab.hydro.lock().map(|g| g.is_some()).unwrap_or(false);
    let key = match lab.modo {
        Modo::Luce => (0u8, lab.light.ticks),
        Modo::Elettrone => (1u8, lab.dirac.ticks + 1_000_000 * lab.cloud_version),
        Modo::Idrogeno => (2u8, lab.h_stato as u64 * 10 + have_atom as u64 + 100 * lab.cloud_version),
    };
    if vis.key != key {
        let items = match lab.modo {
            Modo::Luce => balls_light(&lab),
            Modo::Elettrone => balls_electron(&lab),
            Modo::Idrogeno => balls_atom(&mut lab),
        };
        for mut m in q.iter_mut() {
            let old = m.0.id();
            m.0 = meshes.add(quads(&items));
            meshes.remove(old);
        }
        vis.key = key;
    }
    // the phase of the atom's wave turns with the level's frequency (binding energy / ħ), slowed down
    let omega = {
        let g = lab.hydro.lock().unwrap();
        g.as_ref().map(|h| {
            let idx = [0usize, 1, 2, 3][lab.h_stato];
            -h.e_ext[idx] / atom::H_EVS * 2.0 * std::f64::consts::PI
        }).unwrap_or(0.0)
    };
    let (size, soft) = match lab.modo { Modo::Luce => (0.55, 3.0), Modo::Elettrone => (1.05, 3.0), Modo::Idrogeno => (0.42, 3.5) };
    let cp = (omega * lab.h_time).rem_euclid(2.0 * std::f64::consts::PI) as f32;
    if let Some(m) = mats.get_mut(&vis.mat) {
        m.params.clock_phase = cp;
        m.params.size = size * lab.ball_scale;
        m.params.gain = lab.gain;
        m.params.soft = soft;
    }
}

/// Lines and markers: box, surfers and their trails, ⟨r⟩, nucleus, spin axis.
pub fn disegna(lab: Res<Lab>, mut g: Gizmos) {
    let dim = Color::srgba(0.55, 0.6, 0.75, 0.18);
    match lab.modo {
        Modo::Luce | Modo::Elettrone => {
            let h = if lab.modo == Modo::Luce { N_LIGHT } else { N_ELECTRON } as f32 / 2.0;
            g.cuboid(Transform::from_scale(Vec3::splat(2.0 * h)), dim);
            if lab.modo == Modo::Luce {
                let dir = if lab.light_diag { Vec3::new(1.0, 1.0, -1.0).normalize() } else { Vec3::X };
                arrow(&mut g, -dir * 26.0, dir * 26.0, Color::srgba(1.0, 0.86, 0.45, 0.35));
            } else {
                let white = Color::srgba(1.0, 1.0, 1.0, 0.95);
                for s in &lab.surfers {
                    if s.trail.len() > 1 { g.linestrip(s.trail.iter().copied(), Color::srgba(0.85, 1.0, 0.75, 0.75)); }
                    if let Some(p) = s.trail.last() { g.sphere(Isometry3d::from_translation(*p), 0.35 * lab.ball_scale, Color::srgb(0.85, 1.0, 0.7)); }
                }
                if lab.mean_trail.len() > 1 { g.linestrip(lab.mean_trail.iter().copied(), white); }
                if let Some(p) = lab.mean_trail.last() { g.sphere(Isometry3d::from_translation(*p), 0.5 * lab.ball_scale, white); }
            }
        }
        Modo::Idrogeno => {
            let red = Color::srgb(1.0, 0.25, 0.2);
            g.sphere(Isometry3d::IDENTITY, 0.45 * lab.ball_scale, red);
            g.sphere(Isometry3d::IDENTITY, 0.2 * lab.ball_scale, red);
            // spin axis (physics z = screen up)
            arrow(&mut g, Vec3::new(0.0, -A0_SCENE * 5.0, 0.0), Vec3::new(0.0, A0_SCENE * 5.0, 0.0), Color::srgba(0.7, 0.75, 0.9, 0.3));
            for s in &lab.h_surfers {
                if s.trail.len() > 1 { g.linestrip(s.trail.iter().copied(), Color::srgba(1.0, 1.0, 1.0, 0.7)); }
                if let Some(p) = s.trail.last() { g.sphere(Isometry3d::from_translation(*p), 0.3 * lab.ball_scale, Color::WHITE); }
            }
        }
    }
}

/// A thin arrow with a small head.
fn arrow(g: &mut Gizmos, a: Vec3, b: Vec3, c: Color) {
    g.line(a, b, c);
    let d = (b - a).normalize_or_zero();
    let side = if d.y.abs() < 0.9 { d.cross(Vec3::Y).normalize() } else { d.cross(Vec3::X).normalize() };
    let up = side.cross(d);
    for s in [side, -side, up, -up] { g.line(b, b - d * 2.0 + s * 0.8, c); }
}

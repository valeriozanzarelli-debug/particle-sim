//! Drawing: the electron cloud and the photons as round soft balls (the shader of the photon laboratory,
//! src/photon/balls.wgsl), the nucleus, the surfers and the photon balls (gizmos).
//!
//! Screen axes: the photon comes along +Z; its polarization (the atom's dipole axis, physics z) is Y.
//! Physics (x, y, z) in bohr → screen (x, z, y) × ATOM_SCALE.
//! Meshes are REPLACED, never modified in place (Bevy 0.15 leaks a modified mesh; see `swap_mesh` in
//! src/photon/render.rs).

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
use std::f64::consts::{PI, TAU};

use crate::config::*;
use crate::dmc::Rng;
use crate::fisica::{Fase, Scontro};
use crate::Lab;

pub const BALL_SHADER: Handle<Shader> = Handle::weak_from_u128(0x2c71_9d0e_44a3_4b5f_8e21_6a0c_33d5_7f12);
pub const ATTR_DIR: MeshVertexAttribute = MeshVertexAttribute::new("Ball_Dir", 0x5b3d_0001, VertexFormat::Float32x4);
pub const ATTR_LIFE: MeshVertexAttribute = MeshVertexAttribute::new("Ball_Life", 0x5b3d_0002, VertexFormat::Float32x4);

// The ShaderType derive emits a per-field `check` fn that rustc reports as unused.
#[allow(dead_code)]
mod params {
    use bevy::render::render_resource::ShaderType;
    /// Same layout as `BallParams` in src/photon/balls.wgsl.
    #[derive(Clone, Copy, Debug, ShaderType)]
    pub struct BallParams {
        pub clock: f32,
        pub size: f32,
        pub lim_x: f32,
        pub lim_y: f32,
        pub gain: f32,
        pub soft: f32,
        pub star: f32,
        pub additive: f32,
        pub use_transform: f32,
        pub pad0: f32,
        pub pad1: f32,
        pub pad2: f32,
    }
}
use params::BallParams;

#[derive(Asset, TypePath, AsBindGroup, Clone, Debug)]
pub struct BallMaterial {
    #[uniform(0)]
    pub params: BallParams,
    pub additive: bool,
}

impl Material for BallMaterial {
    fn vertex_shader() -> ShaderRef { BALL_SHADER.into() }
    fn fragment_shader() -> ShaderRef { BALL_SHADER.into() }
    fn alpha_mode(&self) -> AlphaMode { if self.additive { AlphaMode::Add } else { AlphaMode::Blend } }
    fn specialize(
        _pipeline: &MaterialPipeline<Self>,
        descriptor: &mut RenderPipelineDescriptor,
        layout: &MeshVertexBufferLayoutRef,
        _key: MaterialPipelineKey<Self>,
    ) -> Result<(), SpecializedMeshPipelineError> {
        let vl = layout.0.get_layout(&[
            Mesh::ATTRIBUTE_POSITION.at_shader_location(0),
            ATTR_DIR.at_shader_location(1),
            ATTR_LIFE.at_shader_location(2),
            Mesh::ATTRIBUTE_COLOR.at_shader_location(3),
            Mesh::ATTRIBUTE_UV_0.at_shader_location(4),
        ])?;
        descriptor.vertex.buffers = vec![vl];
        descriptor.primitive.cull_mode = None;
        Ok(())
    }
}

pub struct BallsPlugin;

impl Plugin for BallsPlugin {
    fn build(&self, app: &mut App) {
        load_internal_asset!(app, BALL_SHADER, "../photon/balls.wgsl", Shader::from_wgsl);
        app.add_plugins(MaterialPlugin::<BallMaterial> { prepass_enabled: false, shadows_enabled: false, ..default() });
    }
}

/// Which sprite set an entity draws.
#[derive(Component, Clone, Copy, PartialEq, Eq)]
pub enum Sprites { Nuvola, FotoneIn, FotoneOut }

#[derive(Resource)]
pub struct Visuals {
    mats: [Handle<BallMaterial>; 3],
    key_in: u64,
    key_out: u64,
}

/// One sprite: (position, dir + depth, life = (t_emit, z_kill, size, _), colour).
type Item = ([f32; 3], [f32; 4], [f32; 4], [f32; 4]);

fn sprite_mesh(items: &[Item]) -> Mesh {
    let hidden = [([0.0f32; 3], [0.0f32; 4], [1.0e30f32, 0.0, 1.0, 0.0], [0.0f32; 4])];
    let items = if items.is_empty() { &hidden[..] } else { items };
    let n = items.len();
    let (mut pos, mut dir, mut life, mut col, mut uv) =
        (Vec::with_capacity(4 * n), Vec::with_capacity(4 * n), Vec::with_capacity(4 * n), Vec::with_capacity(4 * n), Vec::with_capacity(4 * n));
    let mut idx: Vec<u32> = Vec::with_capacity(6 * n);
    const CORNERS: [[f32; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    for (q, it) in items.iter().enumerate() {
        for c in CORNERS {
            pos.push(it.0);
            dir.push(it.1);
            life.push(it.2);
            col.push(it.3);
            uv.push(c);
        }
        let b = 4 * q as u32;
        idx.extend_from_slice(&[b, b + 1, b + 2, b, b + 2, b + 3]);
    }
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    mesh.insert_attribute(ATTR_DIR, dir);
    mesh.insert_attribute(ATTR_LIFE, life);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, col);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uv);
    mesh.insert_indices(Indices::U32(idx));
    mesh
}

/// Replace the mesh of a sprite set (never modify a mesh in place).
fn swap(which: Sprites, mesh: Mesh, meshes: &mut Assets<Mesh>, q: &mut Query<(&Sprites, &mut Mesh3d)>) {
    for (s, mut m) in q.iter_mut() {
        if *s == which {
            let old = m.0.id();
            m.0 = meshes.add(mesh);
            meshes.remove(old);
            return;
        }
    }
}

fn params(size: f32, gain: f32, soft: f32, additive: bool) -> BallParams {
    BallParams { clock: 0.0, size, lim_x: 1.0e9, lim_y: 1.0e9, gain, soft, star: 0.0, additive: if additive { 1.0 } else { 0.0 },
                 use_transform: 0.0, pad0: 0.0, pad1: 0.0, pad2: 0.0 }
}

pub fn setup_render(mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>, mut mats: ResMut<Assets<BallMaterial>>) {
    let nuvola = mats.add(BallMaterial { params: params(0.32, 0.55, 3.0, true), additive: true });
    let fin = mats.add(BallMaterial { params: params(0.55, 0.9, 3.0, false), additive: false });
    let fout = mats.add(BallMaterial { params: params(0.55, 1.0, 3.0, false), additive: false });
    for (which, mat) in [(Sprites::Nuvola, &nuvola), (Sprites::FotoneIn, &fin), (Sprites::FotoneOut, &fout)] {
        commands.spawn((which, Mesh3d(meshes.add(sprite_mesh(&[]))), MeshMaterial3d(mat.clone()), Transform::default(),
            Visibility::Visible, NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
    }
    commands.insert_resource(Visuals { mats: [nuvola, fin, fout], key_in: 0, key_out: 0 });
}

/// Physics (bohr) → screen.
pub fn schermo(p: [f64; 3]) -> Vec3 {
    let s = ATOM_SCALE as f32;
    Vec3::new(p[0] as f32 * s, p[2] as f32 * s, p[1] as f32 * s)
}

/// Cyclic phase palette of the photon laboratory (src/photon/render.rs).
fn phase_color(angle: f64, alpha: f32) -> [f32; 4] {
    let t = (angle.rem_euclid(TAU) / TAU) as f32;
    const STOPS: [[f32; 3]; 5] = [[0.20, 0.95, 1.00], [0.25, 0.45, 1.00], [0.70, 0.30, 1.00], [1.00, 0.35, 0.70], [1.00, 0.85, 0.30]];
    let x = t * 5.0;
    let i = (x.floor() as usize) % 5;
    let s = x - x.floor();
    let (a, b) = (STOPS[i], STOPS[(i + 1) % 5]);
    [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s, alpha]
}

/// Wave balls of the incoming photon: a Gaussian beam (waist w0 = λ at the atom) with a Gaussian envelope
/// of 1 λ along z, sampled ∝ |E|² at the start; each ball flies straight along its ray at c; colour = phase
/// of the field where it was born (the crests move with the packet). If the photon is absorbed or
/// scattered, its balls stop at the atom's plane (z = 0).
pub fn onde_in(f: &crate::fisica::FotoneIn, seed: u64) -> Vec<Item> {
    let mut rng = Rng::new(seed);
    let lam = f.lam_draw;
    let k = TAU / lam;
    let w0 = lam;
    let zr = PI * w0 * w0 / lam;
    let kill = if f.preso { 0.0 } else { 1.0e9 };
    let mut out = Vec::with_capacity(N_WAVE);
    for _ in 0..N_WAVE {
        let dz = Scontro::sigma_z(lam) * rng.gauss();
        let z = Z_START + dz;
        let w = w0 * (1.0 + (z / zr).powi(2)).sqrt();
        let (x, y) = (0.5 * w * rng.gauss(), 0.5 * w * rng.gauss());
        // ray of the Gaussian beam through (x, y, z): transverse slope ρ z/(z² + z_R²)
        let sl = z / (z * z + zr * zr);
        let d = Vec3::new((x * sl) as f32, (y * sl) as f32, 1.0).normalize();
        let env = (-0.5 * (dz / Scontro::sigma_z(lam)).powi(2)).exp();
        // brightness ∝ the field's energy density cos²(kz − ωt) × envelope: the crests show as bands
        let a = (env.sqrt() * (k * z).cos().powi(2)) as f32;
        out.push(([x as f32, y as f32, z as f32], [d.x, d.y, d.z, 0.0], [0.0, kill, 1.0, 0.0], phase_color(k * z, a)));
    }
    out
}

/// Wave balls of the outgoing photon: a dipole wave, directions ∝ sin²θ about Y, radial envelope of 1 λ;
/// each ball leaves the atom at c (born when the light front reaches it: `depth`).
pub fn onde_out(o: &crate::fisica::FotoneOut, seed: u64) -> Vec<Item> {
    let mut rng = Rng::new(seed);
    let lam = o.lam_draw;
    let k = TAU / lam;
    // a thinner envelope (0.6 λ) than the incoming packet, so that the shells (crests) read as rings
    let sz = 0.6 * Scontro::sigma_z(lam);
    let mut out = Vec::with_capacity(N_WAVE);
    for _ in 0..N_WAVE {
        let n = crate::fisica::dir_dipolo(&mut rng);
        let dr = sz * rng.gauss();
        let depth = 3.0 * Scontro::sigma_z(lam) + dr;
        let env = (-0.5 * (dr / sz).powi(2)).exp();
        out.push(([0.0; 3], [n[0] as f32, n[1] as f32, n[2] as f32, depth.max(0.0) as f32], [0.0, 1.0e9, 1.0, 0.0],
                  phase_color(-k * depth, (env.sqrt() * (k * depth).cos().powi(2)) as f32)));
    }
    out
}

/// Colours of the electrons: spin up cyan, spin down orange.
fn colore(spin: usize) -> [f32; 4] {
    if spin == 0 { [0.30, 0.80, 1.0, 1.0] } else { [1.0, 0.55, 0.20, 1.0] }
}

#[allow(clippy::too_many_arguments)]
pub fn aggiorna(
    mut lab: ResMut<Lab>,
    mut vis: ResMut<Visuals>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut mats: ResMut<Assets<BallMaterial>>,
    mut q: Query<(&Sprites, &mut Mesh3d)>,
) {
    // the cloud: rebuilt every frame (static sprites)
    let pts = lab.sc.punti();
    let items: Vec<Item> = pts.iter().filter(|(_, s)| lab.mostra_e[*s]).map(|(p, s)| {
        let v = schermo(*p);
        ([v.x, v.y, v.z], [0.0; 4], [-1.0e9, 1.0e9, 1.0, 0.0], colore(*s))
    }).collect();
    swap(Sprites::Nuvola, sprite_mesh(&items), &mut meshes, &mut q);
    let sc = &lab.sc;
    // the photons: rebuilt only when a new one appears; the clock moves them
    let id_in = sc.fin.as_ref().map_or(0, |f| f.id);
    if id_in != vis.key_in {
        vis.key_in = id_in;
        let it = sc.fin.as_ref().map_or(vec![], |f| onde_in(f, 0xF0F0 ^ f.id));
        swap(Sprites::FotoneIn, sprite_mesh(&it), &mut meshes, &mut q);
    }
    let id_out = sc.fout.as_ref().map_or(0, |f| f.id);
    if id_out != vis.key_out {
        vis.key_out = id_out;
        let it = sc.fout.as_ref().map_or(vec![], |f| onde_out(f, 0x0F0F ^ f.id));
        swap(Sprites::FotoneOut, sprite_mesh(&it), &mut meshes, &mut q);
    }
    if let Some(f) = &sc.fin {
        if let Some(m) = mats.get_mut(&vis.mats[1]) { m.params.clock = (C_SCR * (sc.t - f.t0)) as f32; }
    }
    if let Some(o) = &sc.fout {
        let st = Scontro::sigma_z(o.lam_draw) / C_SCR;
        if let Some(m) = mats.get_mut(&vis.mats[2]) { m.params.clock = (C_SCR * (sc.t - o.t0 + 3.0 * st)) as f32; }
    }
}

/// Nucleus, surfer, photon balls, polarization.
pub fn disegna(lab: Res<Lab>, mut g: Gizmos) {
    let sc = &lab.sc;
    let gold = Color::srgb(1.0, 0.82, 0.3);
    let white = Color::srgb(1.0, 1.0, 1.0);
    g.sphere(Isometry3d::IDENTITY, 0.35, gold);
    g.sphere(Isometry3d::IDENTITY, 0.18, gold);
    if lab.mostra_surfer {
        let a = schermo([sc.surfer[0], sc.surfer[1], sc.surfer[2]]);
        let b = schermo([sc.surfer[3], sc.surfer[4], sc.surfer[5]]);
        for c in [a, b] {
            g.sphere(Isometry3d::from_translation(c), 0.45, white);
            for d in [Vec3::X, Vec3::Y, Vec3::Z] { g.line(c - d * 1.1, c + d * 1.1, white); }
        }
        g.line(a, b, Color::srgba(1.0, 1.0, 1.0, 0.35));
    }
    // the photon balls (surfers of the photons)
    if let Some(f) = &sc.fin {
        let lam = f.lam_draw;
        let zr = PI * lam;
        let z = Z_START + C_SCR * (sc.t - f.t0);
        let sl = |zz: f64| (1.0 + (zz / zr).powi(2)).sqrt() / (1.0 + (Z_START / zr).powi(2)).sqrt();
        // a ray of the Gaussian beam: transverse position scales with the beam width w(z)
        let p = Vec3::new((f.surf0[0] * sl(z)) as f32, (f.surf0[1] * sl(z)) as f32, z as f32);
        if !(f.preso && z > 0.0) {
            g.sphere(Isometry3d::from_translation(p), 0.9, white);
            g.sphere(Isometry3d::from_translation(p), 1.6, Color::srgba(1.0, 0.95, 0.8, 0.4));
            // polarization: E along Y
            let e = 2.5 * Vec3::Y;
            g.arrow(p - e, p + e, Color::srgba(1.0, 0.95, 0.6, 0.6));
        }
    }
    if let Some(o) = &sc.fout {
        let r = C_SCR * (sc.t - o.t0);
        if r > 0.0 {
            let p = r as f32 * Vec3::new(o.dir[0] as f32, o.dir[1] as f32, o.dir[2] as f32);
            g.sphere(Isometry3d::from_translation(p), 0.9, white);
            g.sphere(Isometry3d::from_translation(p), 1.6, Color::srgba(1.0, 0.95, 0.8, 0.4));
        }
    }
    // the instantaneous dipole of the drawn cloud (×4): it swings up and down during the superposition
    if sc.b2 > 1e-3 && sc.b2 < 1.0 - 1e-3 {
        let d = (4.0 * ATOM_SCALE * sc.dipolo) as f32;
        if d.abs() > 0.2 { g.arrow(Vec3::ZERO, Vec3::new(0.0, d, 0.0), Color::srgb(1.0, 0.9, 0.3)); }
    }
    // the polarization axis of the atom (faint), only while something happens
    if sc.fase == Fase::Volo || sc.fase == Fase::Emissione {
        g.line(Vec3::new(0.0, -6.0, 0.0), Vec3::new(0.0, 6.0, 0.0), Color::srgba(1.0, 0.95, 0.6, 0.15));
    }
}

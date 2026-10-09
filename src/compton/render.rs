//! Drawing the collision as a cloud of round soft balls (own small shader, `balls.wgsl`).
//!
//! • incoming photon wave: balls flying at c along the beam, bright on the crests (spacing λ);
//!   it passes on almost unchanged (the scattered fraction is ~10⁻¹¹ … 10⁻¹⁹ per photon);
//! • scattered photon shell: balls flying out of the electron packet, directions drawn ∝ Klein–
//!   Nishina, colour = wavelength shift (white-gold = none … magenta-red = 2 h/mₑc), crests spaced
//!   by λ'(θ): the bands get wider backwards;
//! • recoil electron partial waves: for each drawn shell direction one electron ball at the
//!   conjugate direction and speed (the pairs of the sum of products); together they form the
//!   recoil "cone"; the electron packet at rest stays (most of it: nothing happened);
//! • the two surfers (sphere + glow + trail) and, after the event, the line of the correlated pair.
//! The shell and the recoil cloud are drawn MUCH brighter than their real weight (stated in the panel).
//! The mesh of a shot is built once; only the clock uniform changes per frame. A new shot REPLACES
//! the mesh asset (Bevy 0.15 leaks the GPU copy of a mesh modified in place; see photon/render.rs).

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

use crate::config::*;
use crate::fisica::*;
use crate::{Lab, Vista};

pub const BALL_SHADER: Handle<Shader> = Handle::weak_from_u128(0x2c7e_91d4_5a03_4f6b_b8e2_10c9_d3a7_6e15);
pub const ATTR_VEL: MeshVertexAttribute = MeshVertexAttribute::new("Ball_Vel", 0x6c0d_0001, VertexFormat::Float32x4);
pub const ATTR_LIFE: MeshVertexAttribute = MeshVertexAttribute::new("Ball_Life", 0x6c0d_0002, VertexFormat::Float32x4);

// The ShaderType derive emits a per-field `check` fn that rustc reports as unused.
#[allow(dead_code)]
mod params {
    use bevy::render::render_resource::ShaderType;
    #[derive(Clone, Copy, Debug, ShaderType)]
    pub struct BallParams {
        pub clock: f32,
        pub size: f32,
        pub gain: f32,
        pub soft: f32,
        pub star: f32,
        pub use_transform: f32,
        pub pad0: f32,
        pub pad1: f32,
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
            ATTR_VEL.at_shader_location(1),
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
        load_internal_asset!(app, BALL_SHADER, "balls.wgsl", Shader::from_wgsl);
        app.add_plugins(MaterialPlugin::<BallMaterial> { prepass_enabled: false, shadows_enabled: false, ..default() })
            .add_systems(Startup, setup);
    }
}

/// The entity drawing all the wave balls of the shot.
#[derive(Component)]
pub struct Onde;

/// The glow of a surfer (0 = photon, 1 = electron).
#[derive(Component)]
pub struct Glow(pub usize);

/// One drawn ball: origin, velocity, (t0, t1, size factor), colour.
type Item = ([f32; 3], [f32; 4], [f32; 4], [f32; 4]);

/// A correlated pair of partial waves drawn as a faint line: photon and electron ball.
#[derive(Clone, Copy, Debug)]
pub struct Coppia {
    pub g_origin: Vec3,
    pub g_vel: Vec3,
    pub e_origin: Vec3,
    pub e_vel: Vec3,
    pub t0: f32,
}

#[derive(Resource)]
pub struct Visuals {
    waves: Handle<BallMaterial>,
    glow: Handle<BallMaterial>,
    key: Option<(u64, u64, bool, u8)>,
    pub coppie: Vec<Coppia>,
}

fn quads(items: &[Item]) -> Mesh {
    let hidden = [([0.0f32; 3], [0.0f32; 4], [1.0e30f32, 0.0, 1.0, 0.0], [0.0f32; 4])];
    let items = if items.is_empty() { &hidden[..] } else { items };
    let n = items.len();
    let mut pos = Vec::with_capacity(4 * n);
    let mut vel = Vec::with_capacity(4 * n);
    let mut life = Vec::with_capacity(4 * n);
    let mut col = Vec::with_capacity(4 * n);
    let mut uv = Vec::with_capacity(4 * n);
    let mut idx: Vec<u32> = Vec::with_capacity(6 * n);
    const CORNERS: [[f32; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    for (q, it) in items.iter().enumerate() {
        for c in CORNERS {
            pos.push(it.0);
            vel.push(it.1);
            life.push(it.2);
            col.push(it.3);
            uv.push(c);
        }
        let b = 4 * q as u32;
        idx.extend_from_slice(&[b, b + 1, b + 2, b, b + 2, b + 3]);
    }
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    mesh.insert_attribute(ATTR_VEL, vel);
    mesh.insert_attribute(ATTR_LIFE, life);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, col);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uv);
    mesh.insert_indices(Indices::U32(idx));
    mesh
}

fn setup(mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>, mut mats: ResMut<Assets<BallMaterial>>) {
    let p = |size: f32, gain: f32, soft: f32, star: f32, tr: f32| BallParams {
        clock: 0.0, size, gain, soft, star, use_transform: tr, pad0: 0.0, pad1: 0.0,
    };
    let waves = mats.add(BallMaterial { params: p(0.55, 1.0, 3.0, 0.0, 0.0) });
    let glow = mats.add(BallMaterial { params: p(4.5, 1.0, 9.0, 1.0, 1.0) });
    commands.spawn((Onde, Mesh3d(meshes.add(quads(&[]))), MeshMaterial3d(waves.clone()), Transform::default(),
                    NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
    for (i, c) in [[1.0, 0.97, 0.85, 0.9], [0.45, 0.95, 1.0, 0.9]].into_iter().enumerate() {
        let item = ([0.0; 3], [0.0; 4], [-1.0e9, 1.0e9, 1.0, 0.0], c);
        commands.spawn((Glow(i), Mesh3d(meshes.add(quads(&[item]))), MeshMaterial3d(glow.clone()),
                        Transform::default(), NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
    }
    commands.insert_resource(Visuals { waves, glow, key: None, coppie: vec![] });
}

/// Colour of a scattered direction by its wavelength shift: x = (1 − cos θ)/2 = Δλ / (2 h/mₑc).
pub fn colore_shift(x: f32) -> [f32; 3] {
    let stops = [[1.0, 0.95, 0.72], [1.0, 0.62, 0.18], [1.0, 0.25, 0.25], [0.85, 0.2, 0.95]];
    let t = x.clamp(0.0, 1.0) * 3.0;
    let i = (t.floor() as usize).min(2);
    let f = t - i as f32;
    let (a, b) = (stops[i], stops[i + 1]);
    [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]
}

const C_IN: [f32; 3] = [1.0, 0.93, 0.72];
const C_E: [f32; 3] = [0.35, 0.85, 1.0];

/// Smooth window over the packet (−½ … ½ of its length).
fn finestra(x: f32) -> f32 { if x.abs() >= 0.5 { 0.0 } else { (std::f32::consts::PI * (x + 0.5)).sin().powi(2) } }

/// Position of a drawn ball along its packet, in wavelengths: the balls sit on the wave crests
/// (one crest per λ, a little spread), so the crests show as thin bands.
fn cresta(rng: &mut Rng) -> f32 {
    let m = (rng.uni() * N_CYCLES as f64).floor() as f32;
    m - 0.5 * (N_CYCLES - 1.0) + 0.07 * rng.gauss() as f32
}

/// All the wave balls of one shot.
fn costruisci(lab: &Lab) -> (Vec<Item>, Vec<Coppia>) {
    let kn = &lab.kn;
    let c = lab.colpo.clone();
    let mut rng = Rng::new(0xC0FFEE ^ c.id.wrapping_mul(0x9E37_79B9));
    let mut items = Vec::with_capacity(N_IN + N_SHELL + N_E_REST + N_E_RECOIL);
    let gain_w = lab.gain_onda;

    let sez = lab.vista == Vista::Sezione;
    // thickness of the section slab (λ) in the section view
    let slab = |rng: &mut Rng| 0.25 * rng.gauss() as f32;
    let pi = std::f64::consts::PI;

    // incoming photon wave (passes on)
    for _ in 0..N_IN {
        let (x, y) = loop {
            let (x, y) = (rng.gauss() as f32, rng.gauss() as f32);
            if x * x + y * y < 6.25 { break (x * S_GAMMA_DRAW, y * S_GAMMA_DRAW); }
        };
        let x = if sez { slab(&mut rng) } else { x };
        let z = cresta(&mut rng);
        let a = if sez { 0.06 } else { 0.22 } * finestra(0.8 * z / N_CYCLES);
        items.push(([x, y, Z_START + z], [0.0, 0.0, C_SCENE, 0.0], [0.0, 1.0e9, 1.0, 0.0], [C_IN[0], C_IN[1], C_IN[2], a]));
    }
    // electron packet at rest (what stays: almost all of it)
    for _ in 0..N_E_REST {
        let mut p = Vec3::new(rng.gauss() as f32, rng.gauss() as f32, rng.gauss() as f32) * S_E_DRAW;
        if sez { p.x = slab(&mut rng); }
        items.push(([p.x, p.y, p.z], [0.0; 4], [-1.0e9, 1.0e9, 0.8, 0.0], [C_E[0], C_E[1], C_E[2], 0.35]));
    }
    // section view: directions in the plane of the event, θ drawn ∝ dσ/dΩ(θ, φ) (no sin θ: a cut)
    let in_plane = |th: f64, ph: f64| kn_domega(kn.k, th, ph, lab.pol);
    let fmax = (0..=400).map(|i| {
        let th = pi * i as f64 / 400.0;
        in_plane(th, c.rot).max(in_plane(th, c.rot + pi))
    }).fold(0.0, f64::max) * 1.05;
    // scattered shell + recoil electron partial waves: the pairs of the sum of products
    let mut coppie = Vec::new();
    for j in 0..N_SHELL {
        let (theta, phi) = if sez {
            loop {
                let th = pi * rng.uni();
                let ph = c.rot + if rng.uni() < 0.5 { 0.0 } else { pi };
                if rng.uni() * fmax < in_plane(th, ph) { break (th, ph); }
            }
        } else {
            let th = kn.mu_of_rank(rng.uni()).clamp(-1.0, 1.0).acos();
            (th, kn.phi_of_rank(th, rng.uni(), lab.pol))
        };
        let n = c.ruota_dir(theta, phi);
        let r = (lam_out(kn.lam_pm, theta) / kn.lam_pm) as f32; // λ'/λ: crest spacing on screen
        let s = r * cresta(&mut rng);
        let mut jit = Vec3::new(rng.gauss() as f32, rng.gauss() as f32, rng.gauss() as f32) * (0.12 * S_E_DRAW);
        if sez { jit.x = slab(&mut rng); }
        let t0 = c.t_w + s / C_SCENE;
        let col = colore_shift(0.5 * (1.0 - theta.cos() as f32));
        let a = gain_w * if sez { 0.5 } else { 0.30 } * finestra(0.8 * s / (N_CYCLES * r));
        let o = jit;
        items.push(([o.x, o.y, o.z], [n.x * C_SCENE, n.y * C_SCENE, n.z * C_SCENE, 0.0], [t0, 1.0e9, 1.0, 0.0],
                    [col[0], col[1], col[2], a]));
        if j < N_E_RECOIL {
            let ev = cinematica(kn.e_kev, theta, phi);
            let ne = c.ruota(ev.pe).normalize_or_zero();
            let ve = ne * ev.beta_e as f32 * C_SCENE;
            let mut je = Vec3::new(rng.gauss() as f32, rng.gauss() as f32, rng.gauss() as f32) * (0.3 * S_E_DRAW);
            if sez { je.x = slab(&mut rng); }
            let oe = je;
            let ae = gain_w * 0.22 * finestra(0.8 * s / (N_CYCLES * r));
            items.push(([oe.x, oe.y, oe.z], [ve.x, ve.y, ve.z, 0.0], [t0, 1.0e9, 0.9, 0.0], [C_E[0], C_E[1], C_E[2], ae]));
            if coppie.len() < N_PAIR_LINES && (lab.vista != Vista::Meta || n.x >= 0.0) {
                coppie.push(Coppia { g_origin: o, g_vel: n * C_SCENE, e_origin: oe, e_vel: ve, t0 });
            }
        }
    }
    if lab.vista == Vista::Meta {
        // cut-away: keep the half of the waves on the far side of the beam plane (screen X ≥ 0)
        let keep = |o: [f32; 3], v: [f32; 4]| if v[0].abs() > 1e-6 { v[0] >= 0.0 } else { o[0] >= 0.0 };
        items.retain(|it| keep(it.0, it.1));
    }
    (items, coppie)
}

/// Update the clock, rebuild the balls when a new shot starts, move the glows.
pub fn aggiorna(
    lab: Res<Lab>,
    mut vis: ResMut<Visuals>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut mats: ResMut<Assets<BallMaterial>>,
    mut onde: Query<&mut Mesh3d, With<Onde>>,
    mut glows: Query<(&Glow, &mut Transform, &mut Visibility)>,
) {
    let key = (lab.colpo.id, lab.e_kev.to_bits() ^ lab.gain_onda.to_bits() as u64, lab.pol == Pol::LineareX, lab.vista as u8);
    if vis.key != Some(key) {
        let (items, coppie) = costruisci(&lab);
        for mut m in onde.iter_mut() {
            let old = m.0.id();
            m.0 = meshes.add(quads(&items));
            meshes.remove(old);
        }
        vis.coppie = coppie;
        vis.key = Some(key);
    }
    if let Some(m) = mats.get_mut(&vis.waves) {
        m.params.clock = lab.clock;
        m.params.size = 0.55 * lab.ball_scale;
    }
    let gs = 4.0 * lab.ball_scale;
    if mats.get(&vis.glow).is_some_and(|m| m.params.size != gs) {
        mats.get_mut(&vis.glow).unwrap().params.size = gs;
    }
    let t = lab.clock;
    for (g, mut tr, mut v) in glows.iter_mut() {
        let p = if g.0 == 0 { lab.colpo.fotone(t) } else { lab.colpo.elettrone(t) };
        tr.translation = p;
        *v = if lab.mostra_surfisti { Visibility::Visible } else { Visibility::Hidden };
    }
}

/// Lines: trails of the surfers, the correlated pair, faint pair lines between partial waves,
/// the recoil-electron fronts, the polarization arrow.
pub fn disegna(lab: Res<Lab>, vis: Res<Visuals>, mut g: Gizmos) {
    let t = lab.clock;
    let c = &lab.colpo;
    let white = Color::srgb(1.0, 0.97, 0.88);
    let cyan = Color::srgb(0.45, 0.9, 1.0);
    // beam axis, faint
    g.line(Vec3::new(0.0, 0.0, Z_START - 10.0), Vec3::new(0.0, 0.0, R_END + 10.0), Color::srgba(0.6, 0.6, 0.7, 0.12));
    // polarization of the incoming photon (physics x → screen up)
    if lab.pol == Pol::LineareX && t < c.t_w + 1.0 {
        let z = Z_START + C_SCENE * t;
        let col = Color::srgba(1.0, 0.9, 0.6, 0.7);
        let e = c.ruota([1.0, 0.0, 0.0]);
        let (a, b) = (Vec3::Z * z - e * 3.0 * S_GAMMA_DRAW, Vec3::Z * z + e * 3.0 * S_GAMMA_DRAW);
        g.line(a, b, col);
        g.line(b, b - e * 1.5 - Vec3::Z, col);
        g.line(b, b - e * 1.5 + Vec3::Z, col);
    }
    if lab.mostra_surfisti {
        if lab.trail_g.len() > 1 { g.linestrip(lab.trail_g.iter().copied(), Color::srgba(1.0, 0.95, 0.8, 0.8)); }
        if lab.trail_e.len() > 1 { g.linestrip(lab.trail_e.iter().copied(), Color::srgba(0.45, 0.9, 1.0, 0.8)); }
        let pg = c.fotone(t);
        let pe = c.elettrone(t);
        g.sphere(Isometry3d::from_translation(pg), 1.3 * lab.ball_scale, white);
        g.sphere(Isometry3d::from_translation(pg), 0.6 * lab.ball_scale, white);
        g.sphere(Isometry3d::from_translation(pe), 1.3 * lab.ball_scale, cyan);
        g.sphere(Isometry3d::from_translation(pe), 0.6 * lab.ball_scale, cyan);
        // the correlated pair: after the event the two surfers are entangled
        if c.evento().is_some() && t > c.t_c {
            g.line(pg, pe, Color::srgba(1.0, 0.45, 0.9, 0.9));
            let mid = 0.5 * (pg + pe);
            g.sphere(Isometry3d::from_translation(mid), 0.35 * lab.ball_scale, Color::srgba(1.0, 0.45, 0.9, 0.9));
        }
    }
    // a few correlated pairs of partial waves (photon direction ↔ electron direction)
    if lab.mostra_coppie {
        for p in &vis.coppie {
            if t < p.t0 { continue; }
            let gp = p.g_origin + p.g_vel * (t - p.t0);
            let ep = p.e_origin + p.e_vel * (t - p.t0);
            g.line(gp, ep, Color::srgba(1.0, 0.5, 0.9, 0.12));
        }
    }
    // the front of the recoil electron waves in four azimuth planes: the "cone"
    if t > c.t_w {
        let dt = t - c.t_w;
        let nk = if lab.vista == Vista::Sezione { 2 } else { 4 };
        for k in 0..nk {
            let phi = c.rot + k as f64 * 2.0 * std::f64::consts::PI / nk as f64;
            let pts = (0..=60).map(|i| {
                let th = std::f64::consts::PI * i as f64 / 60.0;
                let ev = cinematica(lab.kn.e_kev, th, phi);
                c.ruota(ev.pe).normalize_or_zero() * ev.beta_e as f32 * C_SCENE * dt
            });
            g.linestrip(pts, Color::srgba(0.45, 0.9, 1.0, 0.35));
        }
    }
}

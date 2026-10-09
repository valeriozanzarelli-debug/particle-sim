//! Drawing the base 3-D rule as a CLOUD OF BALLS.
//!
//! • Wave balls: round soft sprites (camera-facing quads, additive blending,
//!   `balls.wgsl`), size in world units, colour = phase, brightness = the
//!   ball's amplitude weight (√(|a|/max|a|) × obliquity cos θ). The GPU moves
//!   them: each ball flies at c from its origin, appears when the light front
//!   reaches its object and dies at the next plane, so the mesh is built once
//!   per photon and only the light-clock changes per frame.
//! • What is drawn is a RANDOM SUBSET of the rule's balls, fixed per photon
//!   (seed = laboratory seed × photon id): origins drawn ∝ |a| among the
//!   rule's origins, directions uniform in (sin θx, sin θy) over the unit disc
//!   (= the density of the rule's Nyquist direction grid). The rule itself
//!   uses 10⁷–10¹² balls per photon (see the panel), far too many to draw.
//! • PURELY VISUAL (the photon ball never reads drawn balls):
//!   – each drawn ball is placed a random distance (0 … `DEPTH`) behind the
//!     point its ray has reached, so that the balls of neighbouring origins
//!     do not line up on the exact wave front (the "streaks" of the old 2-D
//!     view); its colour is the phase where it is drawn;
//!   – by default the drawn direction is spread uniformly over its Nyquist
//!     cell: far from a small object the rule's discrete directions would
//!     show as separate clumps (a 6 λ aperture has only 24 directions per
//!     axis). The panel can show the exact grid directions instead.
//! • The photon ball: white sphere + soft glow + its trail; clicks: round dots.

use bevy::{
    asset::load_internal_asset,
    pbr::{MaterialPipeline, MaterialPipelineKey, NotShadowCaster, NotShadowReceiver},
    prelude::*,
    render::{
        mesh::{Indices, MeshVertexAttribute, MeshVertexBufferLayoutRef, PrimitiveTopology},
        render_asset::RenderAssetUsages,
        render_resource::{
            AsBindGroup, RenderPipelineDescriptor, ShaderRef, SpecializedMeshPipelineError, VertexFormat,
        },
        view::NoFrustumCulling,
    },
};

use crate::base3d::Bench3;
use crate::config::*;
use crate::physics::{transverse_basis, Experiment, Lab};
use crate::pilot::{hash2, uniform};
use crate::render::phase_color;
use crate::LabState;

pub const BALL_SHADER: Handle<Shader> = Handle::weak_from_u128(0x6b1f_3c2a_9d47_4e8b_a1c3_5f20_77e4_d901);
pub const ATTR_DIR: MeshVertexAttribute = MeshVertexAttribute::new("Ball_Dir", 0x5b3d_0001, VertexFormat::Float32x4);
pub const ATTR_LIFE: MeshVertexAttribute = MeshVertexAttribute::new("Ball_Life", 0x5b3d_0002, VertexFormat::Float32x4);

// The ShaderType derive emits a per-field `check` fn that rustc reports as
// unused: keep the struct in its own module with the lint allowed there.
#[allow(dead_code)]
mod params {
    use bevy::render::render_resource::ShaderType;

    #[derive(Clone, Copy, Debug, ShaderType)]
    pub struct BallParams {
        pub clock: f32,
        pub size: f32,
        pub lim_x: f32,
        pub lim_y: f32,
        pub gain: f32,
        pub soft: f32,
        /// > 0.5: draw a star (the photon's glow)
        pub star: f32,
        /// > 0.5: additive blending (else alpha blending)
        pub additive: f32,
        /// > 0.5: positions are relative to the entity's transform (the glow);
        /// else absolute (the entity's transform only sets its draw order)
        pub use_transform: f32,
        pub pad0: f32,
        pub pad1: f32,
        pub pad2: f32,
    }
}
pub use params::BallParams;

#[derive(Asset, TypePath, AsBindGroup, Clone, Debug)]
pub struct BallMaterial {
    #[uniform(0)]
    pub params: BallParams,
    /// Additive (clicks, glow: density shows as brightness) or alpha-blended
    /// (wave balls: each ball keeps its own phase colour where they overlap).
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
        load_internal_asset!(app, BALL_SHADER, "balls.wgsl", Shader::from_wgsl);
        app.add_plugins(MaterialPlugin::<BallMaterial> { prepass_enabled: false, shadows_enabled: false, ..default() })
            .add_systems(Startup, setup_visuals);
    }
}

/// One drawn wave ball.
#[derive(Clone, Copy, Debug)]
pub struct VisBall {
    pub origin: Vec3,
    pub dir: Vec3,
    pub depth: f32,
    pub t_emit: f32,
    pub z_kill: f32,
    /// Phase of the ball where it is drawn, minus k·clock (constant in time).
    pub phase: f32,
    pub bright: f32,
    /// Unit Jones vector of its origin (E-field arrows).
    pub jones: [f32; 4],
}

/// Which sprite set an entity draws.
#[derive(Component, Clone, Copy, PartialEq, Eq)]
pub enum Sprites { Waves, Clicks, Glow }

#[derive(Component)]
pub struct PhotonSphere;

#[derive(Resource)]
pub struct Visuals {
    waves: Handle<BallMaterial>,
    clicks: Handle<BallMaterial>,
    glow: Handle<BallMaterial>,
    wave_key: Option<(u64, String, usize, bool)>,
    click_key: Option<(usize, Option<[f32; 3]>, u32)>,
    /// The drawn subset of the current photon.
    pub subset: Vec<VisBall>,
}

/// Size of the drawing relative to the 200 λ benches (Arago's bench is 12 λ).
pub fn vis_scale(lab: &Lab) -> f32 {
    ((lab.screen_z - lab.source_z.min(lab.wall_z)) / 200.0).clamp(0.12, 1.0)
}

/// Orbit camera (focus, radius, yaw, pitch) that frames the bench of `lab`:
/// three-quarter side view, the light going from left to right.
pub fn camera_for(lab: &Lab) -> (Vec3, f32, f32, f32) {
    let z0 = lab.source_z.min(lab.wall_z);
    let len = lab.screen_z - z0;
    (Vec3::new(0.0, 0.0, z0 + 0.5 * len), (1.05 * len + 8.0).max(48.0), -2.0, 0.32)
}

/// Visual depth of the cloud behind the wave front (λ): see the module doc.
const DEPTH: f32 = 3.0;

/// The drawn subset of the rule's wave balls for one photon.
pub fn wave_subset(b: &Bench3, seed: u64, per_fan: usize, scale: f32, exact_dirs: bool) -> Vec<VisBall> {
    let k = b.lev.k as f32;
    let eta = b.rule.eta;
    let mut out = Vec::with_capacity(per_fan * b.fans.len());
    for (m, f) in b.fans.iter().enumerate() {
        if f.l0.is_empty() { continue; }
        let amp: Vec<f64> = f.l0.iter().map(|o| (o.a[0].norm2() + o.a[1].norm2()).sqrt()).collect();
        let amax = amp.iter().cloned().fold(0.0, f64::max).max(1e-300);
        let mut acc = 0.0;
        let cdf: Vec<f64> = amp.iter().map(|a| { acc += a; acc }).collect();
        // the rule's direction grid: d(sin θ) = λ/(η W) per axis, inside the unit circle
        let dsx = b.lam / (eta * f.w[0]);
        let dsy = b.lam / (eta * f.w[1]);
        for n in 0..per_fan {
            let key = hash2(seed ^ 0x0BA1_15D0, ((m as u64) << 32) | n as u64);
            let i = cdf.partition_point(|&c| c < uniform(key, 0) * acc).min(f.l0.len() - 1);
            let o = &f.l0[i];
            let (mut sx, mut sy) = (0.0, 0.0);
            for t in 0..16 {
                let r = uniform(key, 1 + 2 * t).sqrt();
                let ph = std::f64::consts::TAU * uniform(key, 2 + 2 * t);
                sx = r * ph.cos();
                sy = r * ph.sin();
                if exact_dirs {
                    sx = ((sx / dsx).floor() + 0.5) * dsx;
                    sy = ((sy / dsy).floor() + 0.5) * dsy;
                }
                if sx * sx + sy * sy < 1.0 { break; }
                sx = 0.0;
                sy = 0.0;
            }
            let cz = (1.0 - sx * sx - sy * sy).max(0.0).sqrt();
            let depth = DEPTH * scale * uniform(key, 40) as f32;
            let p = if o.a[0].norm2() >= o.a[1].norm2() { 0 } else { 1 };
            let t_emit = f.t_emit as f32;
            let norm = amp[i].max(1e-300);
            out.push(VisBall {
                origin: Vec3::new(o.p[0] as f32, o.p[1] as f32, f.z0 as f32),
                dir: Vec3::new(sx as f32, sy as f32, cz as f32),
                depth,
                t_emit,
                z_kill: f.z_kill as f32,
                phase: o.a[p].arg() as f32 - k * (t_emit + depth),
                bright: ((amp[i] / amax).sqrt() * (0.25 + 0.75 * cz)) as f32,
                jones: [(o.a[0].re / norm) as f32, (o.a[0].im / norm) as f32, (o.a[1].re / norm) as f32, (o.a[1].im / norm) as f32],
            });
        }
    }
    out
}

/// Quads (4 vertices) for sprites: (origin, dir+depth, life, colour).
fn sprite_mesh(items: &[([f32; 3], [f32; 4], [f32; 4], [f32; 4])]) -> Mesh {
    let hidden = [([0.0f32; 3], [0.0f32; 4], [1.0e30f32, 0.0, 1.0, 0.0], [0.0f32; 4])];
    let items = if items.is_empty() { &hidden[..] } else { items };
    let n = items.len();
    let mut pos = Vec::with_capacity(4 * n);
    let mut dir = Vec::with_capacity(4 * n);
    let mut life = Vec::with_capacity(4 * n);
    let mut col = Vec::with_capacity(4 * n);
    let mut uv = Vec::with_capacity(4 * n);
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

/// Replace the mesh of a sprite set (never modify a mesh in place: Bevy 0.15
/// leaks the GPU copy of a modified mesh; see `render::swap_mesh`).
fn swap_sprites(which: Sprites, mesh: Mesh, meshes: &mut Assets<Mesh>, q: &mut Query<(&Sprites, &mut Mesh3d)>) {
    for (s, mut m) in q.iter_mut() {
        if *s == which {
            let old = m.0.id();
            m.0 = meshes.add(mesh);
            meshes.remove(old);
            return;
        }
    }
}

fn params(size: f32, lim: f32, gain: f32, soft: f32, additive: bool) -> BallParams {
    BallParams {
        clock: 0.0, size, lim_x: lim, lim_y: lim, gain, soft, star: 0.0, additive: if additive { 1.0 } else { 0.0 },
        use_transform: 0.0, pad0: 0.0, pad1: 0.0, pad2: 0.0,
    }
}

fn setup_visuals(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut mats: ResMut<Assets<BallMaterial>>,
    mut std_mats: ResMut<Assets<StandardMaterial>>,
) {
    let waves = mats.add(BallMaterial { params: params(0.5, 100.0, 0.9, 3.0, false), additive: false });
    let clicks = mats.add(BallMaterial { params: params(0.4, 1.0e9, 0.9, 2.0, true), additive: true });
    let mut glow = params(4.0, 1.0e9, 1.0, 9.0, true);
    glow.star = 1.0;
    glow.use_transform = 1.0;
    let glow = mats.add(BallMaterial { params: glow, additive: true });
    for (which, mat) in [(Sprites::Waves, &waves), (Sprites::Clicks, &clicks)] {
        commands.spawn((which, Mesh3d(meshes.add(sprite_mesh(&[]))), MeshMaterial3d(mat.clone()), Transform::default(),
            Visibility::Visible, NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
    }
    let glow_quad = sprite_mesh(&[([0.0; 3], [0.0; 4], [-1.0e9, 1.0e9, 1.0, 0.0], [1.0, 0.85, 0.55, 1.0])]);
    commands.spawn((Sprites::Glow, Mesh3d(meshes.add(glow_quad)), MeshMaterial3d(glow.clone()), Transform::default(),
        Visibility::Hidden, NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
    let white = std_mats.add(StandardMaterial { base_color: Color::srgb(1.0, 0.97, 0.88), unlit: true, ..default() });
    commands.spawn((PhotonSphere, Mesh3d(meshes.add(Sphere::new(1.0).mesh().uv(24, 16))), MeshMaterial3d(white),
        Transform::default(), Visibility::Hidden, NotShadowCaster, NotShadowReceiver));
    commands.insert_resource(Visuals { waves, clicks, glow, wave_key: None, click_key: None, subset: Vec::new() });
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn update_visuals(
    state: Res<LabState>,
    vis: Option<ResMut<Visuals>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut mats: ResMut<Assets<BallMaterial>>,
    mut sprites: Query<(&Sprites, &mut Mesh3d)>,
    mut glow: Query<(&Sprites, &mut Transform, &mut Visibility), Without<PhotonSphere>>,
    mut sphere: Query<(&mut Transform, &mut Visibility), (With<PhotonSphere>, Without<Sprites>)>,
    camera: Query<&GlobalTransform, With<Camera3d>>,
) {
    let Some(mut vis) = vis else { return };
    let sim = &state.sim;
    let lab = &sim.lab;
    let scale = vis_scale(lab);
    let user = state.ball_scale;
    let lim = if lab.experiment == Experiment::Arago { 40.0 } else { 100.0 };
    if let Some(m) = mats.get_mut(&vis.waves) {
        m.params.clock = sim.clock;
        m.params.size = 0.5 * scale * user;
        m.params.lim_x = lim;
        m.params.lim_y = lim;
    }
    let click_size = 0.38 * scale.max(0.3) * user;
    if mats.get(&vis.clicks).is_some_and(|m| m.params.size != click_size) {
        mats.get_mut(&vis.clicks).unwrap().params.size = click_size;
    }
    let glow_size = 6.0 * scale.max(0.25) * user;
    if mats.get(&vis.glow).is_some_and(|m| m.params.size != glow_size) {
        mats.get_mut(&vis.glow).unwrap().params.size = glow_size;
    }

    // ── wave balls of the photon in flight (base 3-D rule) ──
    let ball = sim.photons.first().and_then(|p| p.ball.as_ref().map(|b| (p.id, b)));
    let key = ball.map(|(id, b)| (id, b.bench.key.clone(), state.balls_per_fan, state.exact_dirs));
    if vis.wave_key != key {
        vis.subset = match ball {
            Some((id, b)) => wave_subset(&b.bench, hash2(lab.seed, id), state.balls_per_fan, scale, state.exact_dirs),
            None => Vec::new(),
        };
        let items: Vec<_> = vis.subset.iter().map(|v| {
            let c = phase_color(v.phase);
            ([v.origin.x, v.origin.y, v.origin.z], [v.dir.x, v.dir.y, v.dir.z, v.depth],
             [v.t_emit, v.z_kill, 1.0, 0.0], [c[0], c[1], c[2], v.bright])
        }).collect();
        swap_sprites(Sprites::Waves, sprite_mesh(&items), &mut meshes, &mut sprites);
        vis.wave_key = key;
    }

    // ── clicks: round dots (every rule) ──
    let clicks = &sim.stats.clicks;
    let ck = (clicks.len(), clicks.last().map(|c| c.to_array()), (scale * 1000.0) as u32);
    if vis.click_key != Some(ck) {
        vis.click_key = Some(ck);
        let start = clicks.len().saturating_sub(MAX_CLICKS_DRAWN);
        let n = clicks.len() - start;
        let items: Vec<_> = clicks[start..].iter().enumerate().map(|(i, c)| {
            // the newest clicks are a little brighter
            let fresh = if n - i <= 3 { 1.0 } else { 0.55 };
            ([c.x, c.y, c.z], [0.0; 4], [-1.0e9, 1.0e9, 1.0, 0.0], [0.6, 1.0, 0.45, fresh])
        }).collect();
        swap_sprites(Sprites::Clicks, sprite_mesh(&items), &mut meshes, &mut sprites);
    }

    // Draw order of the transparent sprites (sorted back to front by their
    // entity's position): the wave balls and the clicks are placed (as
    // entities; their vertices are absolute) far behind everything, so the
    // photon's glow is always drawn over the cloud.
    let far = camera.iter().next().map_or(Vec3::new(0.0, 0.0, 1.0e6), |c| c.translation() + c.forward() * 1.0e6);
    for (s, mut t, _) in glow.iter_mut() {
        if *s != Sprites::Glow { t.translation = far; }
    }

    // ── the photon ball: sphere + glow ──
    let at = ball.filter(|(_, b)| b.end.is_none() || sim.clock as f64 <= b.t).map(|(_, b)| b.pos_at(sim.clock as f64));
    for (s, mut t, mut v) in glow.iter_mut() {
        if *s != Sprites::Glow { continue; }
        match at {
            Some(p) => { t.translation = p; *v = Visibility::Visible; }
            None => *v = Visibility::Hidden,
        }
    }
    for (mut t, mut v) in sphere.iter_mut() {
        match at {
            Some(p) => {
                t.translation = p;
                t.scale = Vec3::splat(1.1 * scale.max(0.25) * user);
                *v = Visibility::Visible;
            }
            None => *v = Visibility::Hidden,
        }
    }
}

/// E-field arrows on a few drawn balls (off by default; panel toggle).
pub fn arrows(vis: &Visuals, clock: f32, scale: f32, gizmos: &mut Gizmos) {
    let col = Color::srgba(0.7, 0.85, 1.0, 0.5);
    let mut shown = 0;
    for v in &vis.subset {
        if shown >= 600 { break; }
        let s = clock - v.t_emit - v.depth;
        if s < 0.0 { continue; }
        let p = v.origin + v.dir * s;
        if p.z > v.z_kill || p.x.abs() > 100.0 || p.y.abs() > 100.0 { continue; }
        // E = Re[J e^{iφ}]: a ball rides its crest, so its arrow keeps its angle
        let (e1, e2) = transverse_basis(v.dir);
        let (sn, c) = v.phase.sin_cos();
        let ex = v.jones[0] * c - v.jones[1] * sn;
        let ey = v.jones[2] * c - v.jones[3] * sn;
        let e = (e1 * ex + e2 * ey) * 2.0 * scale;
        gizmos.line(p, p + e, col);
        shown += 1;
    }
}

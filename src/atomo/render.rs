//! Drawing: EVERYTHING is a ball (no lines, no wireframes). Own small shader `balls.wgsl`: one camera-facing
//! quad per ball, additive glow.
//!
//! Layers:
//!   • the electron cloud: the electrons of the last `SLOT` ticks of the rule (≈ 90 000 balls), in a ring of
//!     meshes: each new tick replaces ONE mesh, so only ~9 000 balls are uploaded per frame. Colour = shell,
//!     found from the radial distribution of the balls themselves (minima between maxima), nothing labelled;
//!   • the nucleus: a glow ball at the atom scale; at the femtometre scale the B = 1 knot (proton) as 40 000
//!     balls drawn with its baryon density (colour = density);
//!   • the real electron: the zig-zag Dirac cloud (orange = right-handed, blue = left-handed) and ⟨r⟩;
//!   • photons (phase colours, as in the photon laboratory), photo-electrons, the thermal haze.
//! Meshes are REPLACED, never modified in place (Bevy 0.15 leaks a modified mesh; see `swap_mesh` in
//! src/photon/render.rs). Physics (x, y, z) in bohr → screen (x, z, −y), relative to the camera focus
//! (computed in f64, so the femtometre view stays exact).

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
use std::collections::HashMap;
use std::f64::consts::TAU;

use crate::camera::Cam;
use crate::carbonio::Rng;
use crate::mondo::Fotone;
use crate::{costanti as k, Lab};

pub const BALL_SHADER: Handle<Shader> = Handle::weak_from_u128(0x4c2a_91d3_70be_4f18_a2c6_33e1_0d97_5b21);
pub const ATTR_EXTRA: MeshVertexAttribute = MeshVertexAttribute::new("At_Extra", 0x4c2a_0001, VertexFormat::Float32x4);

#[allow(dead_code)]
mod params {
    use bevy::render::render_resource::ShaderType;
    #[derive(Clone, Copy, Debug, ShaderType)]
    pub struct BallParams {
        pub size: f32,
        pub gain: f32,
        pub soft: f32,
        pub core: f32,
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
            ATTR_EXTRA.at_shader_location(1),
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

/// Ticks of the rule shown together (each one a mesh).
pub const SLOT: usize = 10;

#[derive(Component, Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Strato { Nuvola(usize), Nucleo, Nodo, Dirac, Fotoni, Elettroni, Foschia }

/// One ball: position, size factor, colour (alpha = brightness).
pub type Ball = ([f32; 3], f32, [f32; 4]);

pub fn quads(items: &[Ball]) -> Mesh {
    let hidden = [([0.0f32, -1.0e9, 0.0], 0.0f32, [0.0f32; 4])];
    let items = if items.is_empty() { &hidden[..] } else { items };
    let n = items.len();
    let mut pos = Vec::with_capacity(4 * n);
    let mut ex = Vec::with_capacity(4 * n);
    let mut col = Vec::with_capacity(4 * n);
    let mut uv = Vec::with_capacity(4 * n);
    let mut idx: Vec<u32> = Vec::with_capacity(6 * n);
    const CORNERS: [[f32; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    for (q, it) in items.iter().enumerate() {
        for c in CORNERS {
            pos.push(it.0);
            ex.push([it.1, 0.0, 0.0, 0.0]);
            col.push(it.2);
            uv.push(c);
        }
        let b = 4 * q as u32;
        idx.extend_from_slice(&[b, b + 1, b + 2, b, b + 2, b + 3]);
    }
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    mesh.insert_attribute(ATTR_EXTRA, ex);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, col);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uv);
    mesh.insert_indices(Indices::U32(idx));
    mesh
}

/// Physics (bohr) → screen.
#[inline]
pub fn schermo(p: [f64; 3]) -> Vec3 { Vec3::new(p[0] as f32, p[2] as f32, -p[1] as f32) }
#[inline]
fn schermo32(p: [f32; 3]) -> [f32; 3] { [p[0], p[2], -p[1]] }

#[derive(Resource)]
pub struct Visuals {
    mats: HashMap<Strato, Handle<BallMaterial>>,
    /// what each mesh shows (tick of a cloud slot, version of the others)
    chiave: HashMap<Strato, u64>,
    /// balls of each photon in its own frame, by id
    fotoni: HashMap<u64, Vec<Ball>>,
    pub palline_disegnate: usize,
}

fn setup(mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>, mut mats: ResMut<Assets<BallMaterial>>) {
    let mut m: HashMap<Strato, Handle<BallMaterial>> = HashMap::new();
    let mut strati: Vec<Strato> = (0..SLOT).map(Strato::Nuvola).collect();
    strati.extend([Strato::Nucleo, Strato::Nodo, Strato::Dirac, Strato::Fotoni, Strato::Elettroni, Strato::Foschia]);
    for s in strati {
        // all the cloud slots share one material
        let h = if let Strato::Nuvola(i) = s {
            if i == 0 { mats.add(BallMaterial { params: BallParams { size: 0.03, gain: 1.0, soft: 3.0, core: 0.0 } }) }
            else { m[&Strato::Nuvola(0)].clone() }
        } else {
            mats.add(BallMaterial { params: BallParams { size: 0.03, gain: 1.0, soft: 3.0, core: 0.5 } })
        };
        commands.spawn((s, Mesh3d(meshes.add(quads(&[]))), MeshMaterial3d(h.clone()), Transform::default(),
                        NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
        m.insert(s, h);
    }
    commands.insert_resource(Visuals { mats: m, chiave: HashMap::new(), fotoni: HashMap::new(), palline_disegnate: 0 });
}

// ───────────────────────── colours ─────────────────────────
const K_SHELL: [f32; 3] = [1.0, 0.58, 0.22];
const L_SHELL: [f32; 3] = [0.30, 0.72, 1.0];
const M_SHELL: [f32; 3] = [0.78, 0.42, 1.0];

/// Shell colour from the radial distance only, with the boundaries found by the balls.
pub fn colore_guscio(r: f32, confini: &[f64]) -> [f32; 3] {
    let cols = [K_SHELL, L_SHELL, M_SHELL];
    let mut c = cols[0];
    for (i, &b) in confini.iter().enumerate().take(2) {
        let w = 0.12 * b as f32 + 0.03;
        let t = ((r - b as f32) / w * 0.5 + 0.5).clamp(0.0, 1.0);
        let t = t * t * (3.0 - 2.0 * t);
        let n = cols[i + 1];
        c = [c[0] + (n[0] - c[0]) * t, c[1] + (n[1] - c[1]) * t, c[2] + (n[2] - c[2]) * t];
    }
    c
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

/// Black-body colour at temperature T (rough RGB of the Planck spectrum, as in the photon haze).
pub fn colore_corpo_nero(t: f64) -> [f32; 3] {
    let x = (t / 1000.0) as f32;
    if x < 0.8 { [0.35, 0.05, 0.03] } else {
        let r = 1.0;
        let g = (0.25 * (x - 0.8)).clamp(0.0, 1.0).powf(0.6);
        let b = (0.12 * (x - 1.6)).clamp(0.0, 1.0);
        [r, 0.15 + 0.85 * g, 0.05 + 0.95 * b]
    }
}

/// Drawn wavelength of a photon (bohr): NOT to scale (a 166 nm photon is 3100 bohr long, the atom 4).
pub fn lambda_disegno(lambda_nm: f64) -> f64 { 3.0 * (lambda_nm / 165.7).sqrt() }

fn base(d: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let a = if d[2].abs() < 0.9 { [0.0, 0.0, 1.0] } else { [1.0, 0.0, 0.0] };
    let mut e1 = [d[1] * a[2] - d[2] * a[1], d[2] * a[0] - d[0] * a[2], d[0] * a[1] - d[1] * a[0]];
    let n = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
    for v in e1.iter_mut() { *v /= n; }
    let e2 = [d[1] * e1[2] - d[2] * e1[1], d[2] * e1[0] - d[0] * e1[2], d[0] * e1[1] - d[1] * e1[0]];
    (e1, e2)
}

const N_ONDA: usize = 5000;

/// Balls of an incoming photon in its own frame (centre of the packet at 0): a Gaussian beam with an
/// envelope of 1 λ along its direction; colour = phase, brightness ∝ cos² (the crests).
fn palline_fotone_in(f: &Fotone) -> Vec<Ball> {
    let mut rng = Rng::new(f.id * 7919 + 13);
    let lam = lambda_disegno(f.lambda_nm);
    let kk = TAU / lam;
    let (e1, e2) = base(f.dir);
    (0..N_ONDA).map(|_| {
        let z = lam * rng.gauss();
        let (x, y) = (0.6 * lam * rng.gauss(), 0.6 * lam * rng.gauss());
        let env = (-0.5 * (z / lam).powi(2)).exp();
        let a = (0.5 * env.sqrt() * (kk * z).cos().powi(2)) as f32;
        let p = [0, 1, 2].map(|q| z * f.dir[q] + x * e1[q] + y * e2[q]);
        (schermo(p).to_array(), 1.0, phase_color(kk * z, a))
    }).collect()
}

/// Balls of an outgoing photon: a dipole wave (directions ∝ sin²θ about its axis), a shell of 1 λ.
fn palline_fotone_out(f: &Fotone) -> Vec<Ball> {
    let mut rng = Rng::new(f.id * 7919 + 17);
    let lam = lambda_disegno(f.lambda_nm);
    let kk = TAU / lam;
    let sz = 0.6 * lam;
    let mut out = Vec::with_capacity(N_ONDA);
    while out.len() < N_ONDA {
        let v = [rng.gauss(), rng.gauss(), rng.gauss()];
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let u = [v[0] / n, v[1] / n, v[2] / n];
        let c = u[0] * f.dir[0] + u[1] * f.dir[1] + u[2] * f.dir[2];
        if rng.uni() > 1.0 - c * c { continue; }
        let dr = sz * rng.gauss();
        let r = (f.raggio + dr).max(0.0);
        let env = (-0.5 * (dr / sz).powi(2)).exp();
        let a = (0.6 * env.sqrt() * (kk * dr).cos().powi(2)) as f32;
        out.push((schermo([u[0] * r, u[1] * r, u[2] * r]).to_array(), 1.0, phase_color(-kk * dr, a)));
    }
    out
}

fn sostituisci(meshes: &mut Assets<Mesh>, m3: &mut Mesh3d, mesh: Mesh) {
    let old = m3.0.id();
    m3.0 = meshes.add(mesh);
    meshes.remove(old);
}

/// Fade factor: 1 inside [a, b] of the camera distance (log scale), 0 outside, smooth over a factor `w`.
fn finestra(r: f64, lo: f64, hi: f64, w: f64) -> f32 {
    let x = r.ln();
    let up = ((x - (lo / w).ln()) / w.ln()).clamp(0.0, 1.0);
    let down = (((hi * w).ln() - x) / w.ln()).clamp(0.0, 1.0);
    (up * down) as f32
}

#[allow(clippy::too_many_arguments)]
pub fn aggiorna(
    mut lab: ResMut<Lab>,
    cam: Res<Cam>,
    mut vis: ResMut<Visuals>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut mats: ResMut<Assets<BallMaterial>>,
    mut q: Query<(&Strato, &mut Mesh3d, &mut Transform)>,
) {
    let lab = &mut *lab;
    let vis = &mut *vis;
    let p = lab.storia.frame(lab.vista).clone();
    let focus = cam.focus;
    let rel = |a: [f64; 3]| schermo([a[0] - focus[0], a[1] - focus[1], a[2] - focus[2]]);
    let r_cam = cam.r();
    let nucleo = p.nucleo;
    // the snapshots to show and their slots
    let foto = lab.storia.foto_fino_a(p.tick, SLOT);
    let mut wanted: HashMap<usize, usize> = HashMap::new();
    for (i, f) in foto.iter().enumerate() { wanted.insert((f.tick % SLOT as u64) as usize, i); }
    let mut contate = 0usize;
    let confini = lab.confini.clone();
    let alfa_nuvola = 0.11f32;
    let gen_ver = lab.versione_dirac;
    for (s, mut m3, mut tr) in q.iter_mut() {
        match *s {
            Strato::Nuvola(i) => {
                tr.translation = rel(nucleo);
                match wanted.get(&i) {
                    Some(&j) => {
                        let f = foto[j];
                        contate += f.pos.len();
                        let key = f.tick * 8 + confini.len() as u64 + 1000 * lab.versione_colori;
                        if vis.chiave.get(s) != Some(&key) {
                            let items: Vec<Ball> = f.pos.iter().map(|&q| {
                                let r = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
                                let c = colore_guscio(r, &confini);
                                (schermo32(q), 1.0, [c[0], c[1], c[2], alfa_nuvola])
                            }).collect();
                            sostituisci(&mut meshes, &mut m3, quads(&items));
                            vis.chiave.insert(*s, key);
                        }
                    }
                    None => {
                        if vis.chiave.get(s) != Some(&u64::MAX) {
                            sostituisci(&mut meshes, &mut m3, quads(&[]));
                            vis.chiave.insert(*s, u64::MAX);
                        }
                    }
                }
            }
            Strato::Nucleo => {
                tr.translation = rel(nucleo);
                if vis.chiave.get(s).is_none() {
                    sostituisci(&mut meshes, &mut m3, quads(&[([0.0; 3], 1.0, [1.0, 0.35, 0.15, 1.0])]));
                    vis.chiave.insert(*s, 1);
                }
            }
            Strato::Nodo => {
                tr.translation = rel(nucleo);
                if vis.chiave.get(s).is_none() {
                    let fm = 1e-15 / k::bohr_m();
                    let items: Vec<Ball> = lab.nodo_palline.iter().map(|(q, d)| {
                        let d = d.clamp(0.0, 1.0);
                        // colour = density: deep red (tail) → orange → white-gold (core)
                        let c = [1.0, 0.25 + 0.7 * d.powf(0.5), 0.10 + 0.75 * d * d];
                        (schermo([q[0] * fm, q[1] * fm, q[2] * fm]).to_array(), 1.0, [c[0], c[1], c[2], 0.10 + 0.25 * d])
                    }).collect();
                    contate += items.len();
                    sostituisci(&mut meshes, &mut m3, quads(&items));
                    vis.chiave.insert(*s, 1);
                }
            }
            Strato::Dirac => {
                tr.translation = rel([nucleo[0] + lab.punto_elettrone[0], nucleo[1] + lab.punto_elettrone[1], nucleo[2] + lab.punto_elettrone[2]]);
                let key = if lab.mostra_elettrone { gen_ver } else { u64::MAX };
                if vis.chiave.get(s) != Some(&key) {
                    let items = if lab.mostra_elettrone { lab.palline_dirac() } else { vec![] };
                    sostituisci(&mut meshes, &mut m3, quads(&items));
                    vis.chiave.insert(*s, key);
                }
            }
            Strato::Fotoni => {
                tr.translation = rel([0.0; 3]);
                let mut items = Vec::new();
                for f in &p.fotoni {
                    let off = rel(f.pos) - rel([0.0; 3]);
                    let o = [off.x, off.y, off.z];
                    if f.uscente {
                        for b in palline_fotone_out(f) { items.push(([b.0[0] + o[0], b.0[1] + o[1], b.0[2] + o[2]], b.1, b.2)); }
                    } else {
                        let local = vis.fotoni.entry(f.id).or_insert_with(|| palline_fotone_in(f));
                        for b in local.iter() { items.push(([b.0[0] + o[0], b.0[1] + o[1], b.0[2] + o[2]], b.1, b.2)); }
                    }
                }
                let ids: Vec<u64> = p.fotoni.iter().map(|f| f.id).collect();
                vis.fotoni.retain(|id, _| ids.contains(id));
                contate += items.len();
                let key = if items.is_empty() { u64::MAX } else { lab.frame_count };
                if vis.chiave.get(s) != Some(&key) {
                    sostituisci(&mut meshes, &mut m3, quads(&items));
                    vis.chiave.insert(*s, key);
                }
            }
            Strato::Elettroni => {
                tr.translation = rel([0.0; 3]);
                let mut items = Vec::new();
                for (j, e) in p.elettroni.iter().enumerate() {
                    let mut rng = Rng::new(j as u64 + 5);
                    let o = rel(e.pos) - rel([0.0; 3]);
                    for _ in 0..600 {
                        let d = [0.35 * rng.gauss(), 0.35 * rng.gauss(), 0.35 * rng.gauss()];
                        items.push(([o.x + d[0] as f32, o.y + d[1] as f32, o.z + d[2] as f32], 1.0, [0.55, 0.85, 1.0, 0.25]));
                    }
                }
                let key = if items.is_empty() { u64::MAX } else { lab.frame_count };
                if vis.chiave.get(s) != Some(&key) {
                    sostituisci(&mut meshes, &mut m3, quads(&items));
                    vis.chiave.insert(*s, key);
                }
            }
            Strato::Foschia => {
                tr.translation = Vec3::ZERO;
                tr.scale = Vec3::splat((3.0 * r_cam) as f32);
                let key = (lab.t_kelvin / 25.0).round() as u64;
                if vis.chiave.get(s) != Some(&key) {
                    // thermal photons: visual only (number ∝ T³ like the black-body photon density; colour = Planck)
                    let n = (900.0 * (lab.t_kelvin / 3000.0).powi(3)).round() as usize;
                    let c = colore_corpo_nero(lab.t_kelvin);
                    let mut rng = Rng::new(99);
                    let items: Vec<Ball> = (0..n).map(|_| {
                        let v = [rng.gauss(), rng.gauss(), rng.gauss()];
                        let nn = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                        let rr = rng.uni().powf(1.0 / 3.0) / nn;
                        ([(v[0] * rr) as f32, (v[1] * rr) as f32, (v[2] * rr) as f32], 1.0, [c[0], c[1], c[2], 0.5])
                    }).collect();
                    sostituisci(&mut meshes, &mut m3, quads(&items));
                    vis.chiave.insert(*s, key);
                }
            }
        }
    }
    vis.palline_disegnate = contate;
    lab.palline_disegnate = contate;
    // ── sizes and brightness from the zoom ──
    let fm = 1e-15 / k::bohr_m();
    let set = |mats: &mut Assets<BallMaterial>, s: Strato, size: f64, gain: f32, soft: f32, core: f32| {
        if let Some(m) = mats.get_mut(&vis.mats[&s]) {
            m.params = BallParams { size: size as f32, gain, soft, core };
        }
    };
    // cloud: balls of ~0.5 % of the view, never larger than 0.035 bohr
    let s_nuv = (0.0045 * r_cam).min(0.035) * lab.grandezza as f64;
    let g_nuv = lab.luce * (1.0 - 0.75 * finestra(r_cam, 1e-9, 0.02, 2.0)) * if lab.mostra_nuvola { 1.0 } else { 0.0 };
    set(&mut mats, Strato::Nuvola(0), s_nuv, g_nuv, 3.2, 0.25);
    // nucleus glow: a few pixels at the atom scale; it gives way to the knot below ~1000 fm
    let g_nuc = finestra(r_cam, 3000.0 * fm, 1e3, 3.0);
    set(&mut mats, Strato::Nucleo, 0.012 * r_cam, 1.4 * g_nuc, 2.5, 1.0);
    let g_nodo = if lab.mostra_protone { finestra(r_cam, 1e-9, 600.0 * fm, 3.0) } else { 0.0 };
    set(&mut mats, Strato::Nodo, (0.004 * r_cam).min(0.06 * fm), 1.2 * g_nodo, 3.0, 0.3);
    // the Dirac cloud: λ_C = α a₀ per lattice unit
    let g_dir = if lab.mostra_elettrone { finestra(r_cam, 1e-9, 0.6, 3.0) } else { 0.0 };
    set(&mut mats, Strato::Dirac, (0.004 * r_cam).min(0.25 * k::alpha_fs()), 1.0 * g_dir, 3.0, 0.3);
    set(&mut mats, Strato::Fotoni, (0.004 * r_cam).min(0.07), 1.0 * finestra(r_cam, 0.05, 1e4, 3.0), 3.0, 0.3);
    set(&mut mats, Strato::Elettroni, (0.004 * r_cam).min(0.05), 1.0, 3.0, 0.5);
    set(&mut mats, Strato::Foschia, 0.006 * r_cam, 0.6, 3.0, 0.4);
}

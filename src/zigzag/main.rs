//! Zig-zag — ONE rule for light, the electron and the hydrogen atom (docs/ZIGZAG.md).
//!
//! Every wave ball moves at c and turns its internal state (a quantum walk on a cubic lattice):
//!   (a) light: spin-1 walkers carrying F = E + icB, no flip;
//!   (b) free electron: left- and right-handed spin-½ walkers mixed by the mass flip θ_m = m c a/ħ;
//!       the surfer v = c ψ†αψ/ψ†ψ draws the Zitterbewegung loops;
//!   (c) hydrogen: the same electron walker in the Coulomb phase of a point proton (radial form,
//!       levels computed live by "infinite listening", three lattices extrapolated to a → 0).
//!
//! Run the viewer:   cargo run --bin zigzag  [-- --modo luce|elettrone|idrogeno --stato 0..3 --pannello]
//! Tests:            cargo test --bin zigzag

use std::sync::{Arc, Mutex};

use bevy::{
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::{EguiContexts, EguiPlugin};

mod atom;
mod config;
mod cplx;
mod render;
mod ui;
mod walk;
#[path = "../rendering/camera.rs"]
mod camera;
#[cfg(test)]
mod tests;

use atom::{Hydrogen, ALPHA, LAMBDA_C_M, T_UNIT_S};
use camera::{orbit_camera, OrbitState};
use config::*;
use cplx::C;
use walk::{surfer_velocity, Dirac3D, Light3D, Walk1D};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Modo {
    Luce,
    Elettrone,
    Idrogeno,
}

/// The hydrogen states drawn in 3-D (κ, 2 m_j, index in `Hydrogen::states`).
pub const STATI: [(&str, i32, i32); 4] = [("1s1/2", -1, 1), ("2s1/2", -1, 1), ("2p1/2", 1, 1), ("2p3/2", -2, 1)];

pub struct Surfer {
    pub p: [f64; 3],
    pub trail: Vec<Vec3>,
}

/// Small deterministic random generator (xorshift64*).
pub struct Rng(u64);
impl Rng {
    pub fn new(seed: u64) -> Self { Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1) }
    pub fn uni(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn gauss(&mut self) -> f64 {
        let (u, v) = (self.uni().max(1e-300), self.uni());
        (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
    }
}

#[derive(Resource)]
pub struct Lab {
    pub modo: Modo,
    pub paused: bool,
    /// walk ticks per frame (fractional speeds accumulate)
    pub speed: f32,
    acc: f32,
    pub panel_open: bool,
    pub ball_scale: f32,
    pub gain: f32,
    // (a) light
    pub light: Light3D,
    pub light_diag: bool,
    pub light_start: [f64; 3],
    pub light_speed: f64,
    // (b) electron
    pub dirac: Dirac3D,
    /// false: purely right-handed start (zig-zag along the spin); true: mixed spinor (Zitterbewegung circle)
    pub cerchio: bool,
    pub surfers: Vec<Surfer>,
    pub mean_trail: Vec<Vec3>,
    pub zbw: Vec<[f64; 3]>,
    pub zbw_radius: f64,
    pub zbw_omega: f64,
    pub checker: Vec<[u8; 4]>,
    // (c) hydrogen
    pub hydro: Arc<Mutex<Option<Hydrogen>>>,
    pub h_stato: usize,
    pub h_surfers: Vec<Surfer>,
    /// physical time of the atom (s) and physical seconds per screen second
    pub h_time: f64,
    pub t_vis: f64,
    pub cloud_version: u64,
    pub frame: u64,
    rng: Rng,
}

/// Lattice index → scene position (physics z = screen up).
pub fn scene(p: [f64; 3]) -> Vec3 { Vec3::new(p[0] as f32, p[2] as f32, -p[1] as f32) }

impl Lab {
    fn new(modo: Modo) -> Self {
        let hydro = Arc::new(Mutex::new(None));
        let h2 = hydro.clone();
        std::thread::spawn(move || {
            let h = Hydrogen::compute();
            *h2.lock().unwrap() = Some(h);
        });
        let mut lab = Lab {
            modo, paused: false, speed: 1.0, acc: 0.0, panel_open: false, ball_scale: 1.0, gain: 1.0,
            light: Light3D::new(N_LIGHT), light_diag: false, light_start: [0.0; 3], light_speed: 0.0,
            dirac: Dirac3D::new(N_ELECTRON, A_ELECTRON, A_ELECTRON), cerchio: false, surfers: vec![], mean_trail: vec![], zbw: vec![],
            zbw_radius: 0.0, zbw_omega: 0.0, checker: checkerboard(),
            hydro, h_stato: 0, h_surfers: vec![], h_time: 0.0, t_vis: T_VIS_ATOM, cloud_version: 0, frame: 0,
            rng: Rng::new(20261009),
        };
        lab.reset_light();
        lab.reset_electron();
        lab
    }

    pub fn reset_light(&mut self) {
        let n = N_LIGHT;
        let k = 2.0 * std::f64::consts::PI / LAMBDA_CELLS;
        let (dir, e): ([f64; 3], [f64; 3]) = if self.light_diag {
            let s = 1.0 / 3f64.sqrt();
            ([s, s, s], [std::f64::consts::FRAC_1_SQRT_2, -std::f64::consts::FRAC_1_SQRT_2, 0.0])
        } else {
            ([1.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        };
        // k̂ × ê
        let b = [dir[1] * e[2] - dir[2] * e[1], dir[2] * e[0] - dir[0] * e[2], dir[0] * e[1] - dir[1] * e[0]];
        let c0 = [-14.0 * dir[0], -14.0 * dir[1], -14.0 * dir[2]];
        let half = n as f64 / 2.0 - 0.5;
        let mut l = Light3D::new(n);
        for s in 0..n * n * n {
            let x = [(s / (n * n)) as f64 - half, ((s / n) % n) as f64 - half, (s % n) as f64 - half];
            let d = [x[0] - c0[0], x[1] - c0[1], x[2] - c0[2]];
            // envelope: longer along the motion (σ∥ = 1.4 σ) than across (σ⊥ = σ)
            let par = d[0] * dir[0] + d[1] * dir[1] + d[2] * dir[2];
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - par * par;
            let env = (-r2 / (4.0 * SIGMA_LIGHT * SIGMA_LIGHT) - par * par / (4.0 * 1.96 * SIGMA_LIGHT * SIGMA_LIGHT)).exp();
            if env < 1e-6 { continue; }
            let ph = (k * (d[0] * dir[0] + d[1] * dir[1] + d[2] * dir[2])).cos() * env;
            // real fields: E = ê cos, cB = k̂×ê cos  → F = E + i cB
            l.f[s] = [C::new(e[0] * ph, b[0] * ph), C::new(e[1] * ph, b[1] * ph), C::new(e[2] * ph, b[2] * ph)];
        }
        self.light = l;
        self.light_start = centroid_light(&self.light);
        self.light_speed = 0.0;
    }

    pub fn reset_electron(&mut self) {
        let n = N_ELECTRON;
        let mut d = Dirac3D::new(n, A_ELECTRON, A_ELECTRON);
        // zig-zag start: a purely RIGHT-handed spin-up ball cloud (R₀ = 1): the flip turns it left-handed and
        // back at 2mc²/ħ, and the surfers run up and down the spin axis at up to c.
        // circle start: (L₀, L₁, R₀, R₁) = (0, i, 1, 0)/√2 → ⟨r⟩ runs on a circle of radius ħ/2mc (x–z plane).
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let sp = if self.cerchio { [C::ZERO, C::new(0.0, s), C::real(s), C::ZERO] } else { [C::ZERO, C::ZERO, C::ONE, C::ZERO] };
        let mut norm = 0.0;
        for i in 0..n * n * n {
            let x = [d.x(i / (n * n)), d.x((i / n) % n), d.x(i % n)];
            let r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
            let g = (-r2 / (4.0 * SIGMA_ELECTRON * SIGMA_ELECTRON)).exp();
            d.psi[i] = [sp[0] * g, sp[1] * g, sp[2] * g, sp[3] * g];
            norm += g * g;
        }
        let k = 1.0 / norm.sqrt();
        for p in d.psi.iter_mut() { for c in p.iter_mut() { *c = c.scale(k); } }
        self.dirac = d;
        // surfers born with |ψ|² (a Gaussian of std σ/√2 per axis), kept within one σ for a clear picture
        self.surfers.clear();
        while self.surfers.len() < N_SURFERS {
            let p = [self.rng.gauss(), self.rng.gauss(), self.rng.gauss()].map(|v| v * SIGMA_ELECTRON / 2f64.sqrt());
            if (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() < SIGMA_ELECTRON {
                self.surfers.push(Surfer { p, trail: vec![] });
            }
        }
        self.mean_trail.clear();
        self.zbw.clear();
        self.cloud_version += 1;
    }

    pub fn reset_atom_surfers(&mut self) {
        self.h_surfers.clear();
        self.h_time = 0.0;
        let guard = self.hydro.lock().unwrap();
        let Some(h) = guard.as_ref() else { return; };
        let lv = &h.states[self.h_stato];
        let pts = render::sample_cloud(lv, STATI[self.h_stato].1, STATI[self.h_stato].2, N_H_SURFERS, &mut self.rng);
        drop(guard);
        self.h_surfers = pts.into_iter().map(|(p, _)| Surfer { p, trail: vec![] }).collect();
    }

    fn step_light(&mut self) {
        self.light.tick();
        let c = centroid_light(&self.light);
        let t = self.light.ticks as f64;
        let dir = if self.light_diag { [1.0 / 3f64.sqrt(); 3] } else { [1.0, 0.0, 0.0] };
        let n = N_LIGHT as f64;
        let d: Vec<f64> = (0..3).map(|a| (c[a] - self.light_start[a] + n / 2.0).rem_euclid(n) - n / 2.0).collect();
        self.light_speed = (d[0] * dir[0] + d[1] * dir[1] + d[2] * dir[2]) / t;
        if t >= 34.0 { self.reset_light(); }
    }

    fn step_electron(&mut self) {
        let a = self.dirac.a;
        let v1: Vec<[f64; 3]> = self.surfers.iter().map(|s| surfer_velocity(&self.dirac.sample(s.p))).collect();
        self.dirac.tick();
        for (s, v) in self.surfers.iter_mut().zip(v1) {
            let q = [s.p[0] + v[0] * a, s.p[1] + v[1] * a, s.p[2] + v[2] * a];
            let v2 = surfer_velocity(&self.dirac.sample(q));
            for k in 0..3 { s.p[k] += 0.5 * (v[k] + v2[k]) * a; }
            let pc = scene([s.p[0] / a, s.p[1] / a, s.p[2] / a]);
            s.trail.push(pc);
            if s.trail.len() > 160 { s.trail.remove(0); }
        }
        // ⟨r⟩ of the wave: the Zitterbewegung circle
        let m = mean_dirac(&self.dirac);
        self.zbw.push(m);
        self.mean_trail.push(scene([m[0] / a, m[1] / a, m[2] / a]));
        if self.mean_trail.len() > 400 { self.mean_trail.remove(0); }
        // radius and frequency from the last ~2 periods (period π ħ/mc² = π/a ticks)
        let per = (std::f64::consts::PI / a).round() as usize;
        if self.zbw.len() > 2 * per + 2 {
            let w = &self.zbw[self.zbw.len() - 2 * per..];
            let span = |k: usize| {
                let (lo, hi) = w.iter().fold((f64::MAX, f64::MIN), |(lo, hi), v| (lo.min(v[k]), hi.max(v[k])));
                0.5 * (hi - lo)
            };
            self.zbw_radius = if self.cerchio { 0.5 * (span(0) + span(2)) } else { span(2) };
            // frequency from the angle swept by (⟨x⟩, ⟨z⟩) around its mean
            let cx = w.iter().map(|v| v[0]).sum::<f64>() / w.len() as f64;
            let cz = w.iter().map(|v| v[2]).sum::<f64>() / w.len() as f64;
            let mut ang = 0.0;
            for p in w.windows(2) {
                let a0 = (p[0][2] - cz).atan2(p[0][0] - cx);
                let a1 = (p[1][2] - cz).atan2(p[1][0] - cx);
                let mut d = a1 - a0;
                while d > std::f64::consts::PI { d -= 2.0 * std::f64::consts::PI; }
                while d < -std::f64::consts::PI { d += 2.0 * std::f64::consts::PI; }
                ang += d;
            }
            self.zbw_omega = ang.abs() / ((w.len() - 1) as f64 * a);
            if !self.cerchio {
                // linear zig-zag along z: frequency from the crossings of the mean
                let mut cross = vec![];
                for (q, p) in w.windows(2).enumerate() {
                    if (p[0][2] - cz).signum() != (p[1][2] - cz).signum() { cross.push(q as f64); }
                }
                if cross.len() >= 3 {
                    let half_per = (cross[cross.len() - 1] - cross[0]) / (cross.len() - 1) as f64 * a;
                    self.zbw_omega = std::f64::consts::PI / half_per;
                }
            }
        }
        if self.dirac.ticks >= 80 { self.reset_electron(); }
    }

    fn step_atom(&mut self, dt: f64) {
        let guard = self.hydro.lock().unwrap();
        let Some(h) = guard.as_ref() else { return; };
        let lv = &h.states[self.h_stato];
        let (kappa, tm, _) = (STATI[self.h_stato].1, STATI[self.h_stato].2, 0);
        let dts = dt * self.t_vis;
        self.h_time += dts;
        // surfers: v (units of c) → displacement in λ_C: v · dt / t_unit (RK2)
        let k = dts / T_UNIT_S;
        for s in self.h_surfers.iter_mut() {
            let v1 = atom::velocity_dirac(&atom::spinor(lv, kappa, tm, s.p));
            let q = [s.p[0] + v1[0] * k, s.p[1] + v1[1] * k, s.p[2] + v1[2] * k];
            let v2 = atom::velocity_dirac(&atom::spinor(lv, kappa, tm, q));
            for i in 0..3 { s.p[i] += 0.5 * (v1[i] + v2[i]) * k; }
            let sc = (ALPHA * A0_SCENE as f64) as f32;
            s.trail.push(scene(s.p) * sc);
            if s.trail.len() > 240 { s.trail.remove(0); }
        }
    }

    fn avanza(&mut self, dt: f32) {
        self.frame += 1;
        if self.paused { return; }
        match self.modo {
            Modo::Luce | Modo::Elettrone => {
                self.acc += self.speed;
                while self.acc >= 1.0 {
                    self.acc -= 1.0;
                    if self.modo == Modo::Luce { self.step_light(); } else { self.step_electron(); }
                }
            }
            Modo::Idrogeno => {
                if self.h_surfers.is_empty() { self.reset_atom_surfers(); }
                self.step_atom((dt * self.speed) as f64);
            }
        }
    }
}

/// Centre of the light energy |F|² as a circular mean on the periodic box (the tails that wrap around
/// do not pull it back), in cells from the box centre.
pub fn centroid_light(l: &Light3D) -> [f64; 3] {
    let n = l.n;
    let half = n as f64 / 2.0 - 0.5;
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let mut s = [C::ZERO; 3];
    for (i, p) in l.f.iter().enumerate() {
        let w: f64 = p.iter().map(|c| c.norm2()).sum();
        if w == 0.0 { continue; }
        let x = [(i / (n * n)) as f64 - half, ((i / n) % n) as f64 - half, (i % n) as f64 - half];
        for a in 0..3 { s[a] += C::expi(k * x[a]) * w; }
    }
    [s[0].arg() / k, s[1].arg() / k, s[2].arg() / k]
}

pub fn mean_dirac(d: &Dirac3D) -> [f64; 3] {
    let n = d.n;
    let mut s = [0.0; 4];
    for (i, p) in d.psi.iter().enumerate() {
        let w: f64 = p.iter().map(|c| c.norm2()).sum();
        s[0] += w * d.x(i / (n * n));
        s[1] += w * d.x((i / n) % n);
        s[2] += w * d.x(i % n);
        s[3] += w;
    }
    [s[0] / s[3], s[1] / s[3], s[2] / s[3]]
}

/// Feynman checkerboard: a 1-D walk started as ONE right-handed ball, θ = 0.08; space-time picture
/// (orange = right-handed, blue = left-handed, brightness = amplitude).
pub fn checkerboard() -> Vec<[u8; 4]> {
    let n = CHECK_W;
    let mut w = Walk1D { theta: 0.08, r: vec![C::ZERO; n], l: vec![C::ZERO; n] };
    w.r[n / 2] = C::ONE;
    let mut img = vec![[0u8, 0, 0, 255]; n * CHECK_T];
    for t in 0..CHECK_T {
        for i in 0..n {
            let (r, l) = (w.r[i].abs(), w.l[i].abs());
            let f = |x: f64| (255.0 * (x * 3.0).sqrt().min(1.0)) as u8;
            let row = CHECK_T - 1 - t;
            img[row * n + i] = [f(r).max(f(l) / 3), (f(r) as f32 * 0.55 + f(l) as f32 * 0.6) as u8, f(l), 255];
        }
        w.tick();
    }
    img
}

/// True while the pointer is over the HUD, so the camera ignores it.
#[derive(Resource, Default)]
struct PointerOverHud(bool);

fn arg<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).map(|s| s.as_str())
}

pub fn camera_for(modo: Modo, stato: usize) -> OrbitState {
    match modo {
        Modo::Luce => OrbitState { yaw: 0.35, pitch: 0.35, radius: 95.0, focus: Vec3::ZERO },
        Modo::Elettrone => OrbitState { yaw: 0.25, pitch: 0.18, radius: 105.0, focus: Vec3::ZERO },
        Modo::Idrogeno => OrbitState { yaw: 0.5, pitch: 0.3, radius: if stato == 0 { 95.0 } else { 230.0 }, focus: Vec3::ZERO },
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let modo = match arg(&args, "--modo") {
        Some("luce") => Modo::Luce,
        Some("idrogeno") => Modo::Idrogeno,
        _ => Modo::Elettrone,
    };
    let mut lab = Lab::new(modo);
    lab.panel_open = args.iter().any(|a| a == "--pannello");
    if let Some(s) = arg(&args, "--stato").and_then(|s| s.parse::<usize>().ok()) { lab.h_stato = s.min(STATI.len() - 1); }
    if args.iter().any(|a| a == "--diagonale") { lab.light_diag = true; lab.reset_light(); }
    if args.iter().any(|a| a == "--cerchio") { lab.cerchio = true; lab.reset_electron(); }
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Zig-zag — una regola per luce, elettrone e atomo".into(),
                resolution: WindowResolution::new(1600.0, 900.0),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        .add_plugins(render::BallsPlugin)
        .insert_resource(ClearColor(Color::srgb(0.012, 0.016, 0.03)))
        .insert_resource(camera_for(modo, lab.h_stato))
        .insert_resource(lab)
        .init_resource::<PointerOverHud>()
        .add_systems(Startup, setup_scene)
        .add_systems(Update, (
            run_lab,
            render::aggiorna.after(run_lab),
            render::disegna.after(run_lab),
            ui::draw_hud,
            hud_pointer.after(ui::draw_hud),
            orbit_camera.run_if(|p: Res<PointerOverHud>| !p.0),
        ))
        .run();
}

fn setup_scene(mut commands: Commands, mut store: ResMut<GizmoConfigStore>) {
    let (cfg, _) = store.config_mut::<DefaultGizmoConfigGroup>();
    cfg.line_width = 2.5;
    commands.spawn((Camera3d::default(), Transform::from_xyz(0.0, 40.0, 120.0).looking_at(Vec3::ZERO, Vec3::Y)));
}

fn hud_pointer(mut ctx: EguiContexts, mut over: ResMut<PointerOverHud>) {
    let c = ctx.ctx_mut();
    over.0 = c.wants_pointer_input() || c.is_pointer_over_area();
}

fn run_lab(time: Res<Time>, mut lab: ResMut<Lab>) {
    lab.avanza(time.delta_secs());
}

/// λ_C in metres (for the panel).
pub fn lambda_c_m() -> f64 { LAMBDA_C_M }

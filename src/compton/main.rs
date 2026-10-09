//! Compton — one photon hits one FREE electron at rest (docs/COMPTON.md).
//!
//! Before: a photon wave-ball packet and an electron wave packet with their surfers.
//! The outgoing two-particle wave is a sum of products over photon directions θ (photon partial
//! wave ⊗ recoil-electron partial wave), the kinematics fixed by energy–momentum conservation,
//! the amplitudes from QED's Klein–Nishina formula (an INPUT, not derived from the balls).
//! Only a tiny fraction of the incoming wave scatters; the rest passes straight on.
//! The surfers pick ONE correlated pair (or nothing happens) by their rank in the flux.
//!
//! Run the viewer:   cargo run --bin compton  [-- --kev 17.48 --pol u|x --seed 1 --ferma-a 9.0 --pannello --statistica --velocita 1]
//! Measures / CSV:   cargo run --bin compton -- --misura          (validation table)
//!                   cargo run --bin compton -- --eventi 20000 --kev 662 --seed 1 --pol u --out file.csv
//! Tests:            cargo test --bin compton

use bevy::{
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::{EguiContexts, EguiPlugin};

mod config;
mod fisica;
mod misura;
mod render;
mod scena;
mod ui;
#[path = "../rendering/camera.rs"]
mod camera;
#[cfg(test)]
mod tests;

use camera::{orbit_camera, OrbitState};
use config::*;
use fisica::{Kn, Pacchetti, Pol, Rng};
use scena::{Colpo, Statistica};

/// How the waves are drawn.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Vista {
    /// a thin slice through the plane of the event (beam, photon, electron): crests as rings
    Sezione,
    /// the far half of the 3-D waves (cut-away)
    Meta,
    /// all of the 3-D waves
    Intera,
}

/// The running viewer.
#[derive(Resource)]
pub struct Lab {
    pub e_kev: f64,
    pub kn: Kn,
    pub pk: Pacchetti,
    pub pol: Pol,
    /// show only shots whose surfers are born inside the scattering disk (conditioned)
    pub solo_diffusi: bool,
    pub colpo: Colpo,
    /// screen time of the shot (s)
    pub clock: f32,
    pub speed: f32,
    pub paused: bool,
    pub auto: bool,
    pub trail_g: Vec<Vec3>,
    pub trail_e: Vec<Vec3>,
    pub n_colpi: u64,
    pub n_diffusi: u64,
    pub panel_open: bool,
    pub mostra_surfisti: bool,
    pub mostra_coppie: bool,
    /// visual brightness factor of the scattered shell and of the recoil cloud
    pub gain_onda: f32,
    pub ball_scale: f32,
    pub vista: Vista,
    pub stat: Option<Statistica>,
    /// `--ferma-a T`: pause when the clock reaches T (for unattended screenshots)
    pub ferma_a: Option<f32>,
    rng: Rng,
    next_id: u64,
}

pub fn pacchetti(kn: &Kn) -> Pacchetti {
    Pacchetti { s_gamma_m: S_GAMMA_LAMBDA * kn.lam_pm * 1e-12, s_e_m: S_E_M }
}

impl Lab {
    fn new(e_kev: f64, pol: Pol, seed: u64) -> Self {
        let kn = Kn::new(e_kev);
        let pk = pacchetti(&kn);
        let mut rng = Rng::new(seed);
        let colpo = Colpo::nuovo(1, &mut rng, &kn, &pk, pol, true);
        Lab { e_kev, kn, pk, pol, solo_diffusi: true, colpo, clock: 0.0, speed: 1.0, paused: false, auto: true,
              trail_g: vec![], trail_e: vec![], n_colpi: 1, n_diffusi: 1, panel_open: false,
              mostra_surfisti: true, mostra_coppie: true, gain_onda: 1.0, ball_scale: 1.0, vista: Vista::Sezione, stat: None,
              ferma_a: None, rng, next_id: 2 }
    }

    pub fn nuovo_colpo(&mut self) {
        self.colpo = Colpo::nuovo(self.next_id, &mut self.rng, &self.kn, &self.pk, self.pol, self.solo_diffusi);
        self.next_id += 1;
        self.n_colpi += 1;
        if self.colpo.evento().is_some() { self.n_diffusi += 1; }
        self.clock = 0.0;
        self.trail_g.clear();
        self.trail_e.clear();
    }

    pub fn cambia_energia(&mut self, e_kev: f64) {
        if (e_kev - self.e_kev).abs() < 1e-12 * e_kev { return; }
        self.e_kev = e_kev;
        self.kn = Kn::new(e_kev);
        self.pk = pacchetti(&self.kn);
        self.stat = None;
        self.n_colpi = 0;
        self.n_diffusi = 0;
        self.nuovo_colpo();
    }

    fn avanza(&mut self, dt: f32) {
        if !self.paused {
            self.clock += dt * self.speed;
            if let Some(t) = self.ferma_a {
                if self.clock >= t { self.clock = t; self.paused = true; self.ferma_a = None; }
            }
        }
        if self.clock > self.colpo.t_end + PAUSE_END && self.auto && !self.paused {
            self.nuovo_colpo();
        }
        let t = self.clock.min(self.colpo.t_end);
        let (pg, pe) = (self.colpo.fotone(t), self.colpo.elettrone(t));
        if self.trail_g.last().is_none_or(|p| p.distance(pg) > 0.4) { self.trail_g.push(pg); }
        if self.trail_e.last().is_none_or(|p| p.distance(pe) > 0.1) { self.trail_e.push(pe); }
    }
}

/// True while the pointer is over the HUD, so the camera ignores it.
#[derive(Resource, Default)]
struct PointerOverHud(bool);

fn arg<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).map(|s| s.as_str())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let e_kev: f64 = arg(&args, "--kev").and_then(|s| s.parse().ok()).unwrap_or(17.48);
    let seed: u64 = arg(&args, "--seed").and_then(|s| s.parse().ok()).unwrap_or(1);
    let pol = match arg(&args, "--pol") { Some("x") => Pol::LineareX, _ => Pol::NonPolarizzata };
    if args.iter().any(|a| a == "--misura") {
        misura::tabella();
        return;
    }
    if let Some(n) = arg(&args, "--eventi").and_then(|s| s.parse().ok()) {
        let out = arg(&args, "--out").unwrap_or("lab/out/compton_eventi.csv");
        misura::csv(e_kev, pol, n, seed, out);
        return;
    }
    let mut lab = Lab::new(e_kev, pol, seed);
    lab.ferma_a = arg(&args, "--ferma-a").and_then(|s| s.parse().ok());
    lab.panel_open = args.iter().any(|a| a == "--pannello");
    lab.vista = match arg(&args, "--vista") { Some("meta") => Vista::Meta, Some("3d") => Vista::Intera, _ => Vista::Sezione };
    if let Some(v) = arg(&args, "--velocita").and_then(|s| s.parse().ok()) { lab.speed = v; }
    if args.iter().any(|a| a == "--statistica") { lab.stat = Some(scena::statistica(e_kev, pol, 20_000, 2026)); }
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Compton — fotone contro elettrone libero".into(),
                resolution: WindowResolution::new(1600.0, 900.0),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        .add_plugins(render::BallsPlugin)
        .insert_resource(ClearColor(Color::srgb(0.012, 0.016, 0.03)))
        .insert_resource(lab)
        .insert_resource(OrbitState { yaw: -std::f32::consts::FRAC_PI_2, pitch: 0.06, radius: CAM_RADIUS_DEFAULT, focus: Vec3::new(0.0, 0.0, 8.0) })
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

fn setup_scene(mut commands: Commands) {
    commands.spawn((Camera3d::default(), Transform::from_xyz(-120.0, 50.0, 60.0).looking_at(Vec3::ZERO, Vec3::Y)));
}

fn hud_pointer(mut ctx: EguiContexts, mut over: ResMut<PointerOverHud>) {
    let c = ctx.ctx_mut();
    over.0 = c.wants_pointer_input() || c.is_pointer_over_area();
}

fn run_lab(time: Res<Time>, mut lab: ResMut<Lab>) {
    lab.avanza(time.delta_secs());
}

//! Molecules from the electron-ball rule — level 1 → 2 of the zoom levels (atoms → molecules).
//!
//! Every electron "wave ball" does, each tick, only: a random step (D = ħ/2mₑ) and duplicate-or-vanish
//! by the local energy (diffusion Monte Carlo written as a rule per ball; see `dmc.rs`).
//! Inputs: ħ, e, mₑ, ε0, c, nuclear charges and masses. Nothing else.
//!
//! Run the viewer:     cargo run --release --bin molecole   [--sistema h2|h2piu|he|2p] [--giro]
//! Run the measures:   cargo run --release --bin molecole -- --misura idrogeno|eccitati|nodo2s|nodo3p|campo|riga|h2piu|h2|elio
//!                     (all of them: sh lab/molecole_misura.sh; numbers -> assets/molecole/risultati.txt)
//! Tests:              cargo test --release --bin molecole

use bevy::{
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::{EguiContexts, EguiPlugin};
use std::time::Instant;

mod config;
mod costanti;
mod dmc;
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
use misura::Risultati;
use scena::Scelta;

/// The running viewer: one system, its balls, the measured numbers.
#[derive(Resource)]
pub struct Lab {
    pub scelta: Scelta,
    pub dmc: dmc::Dmc,
    pub res: Risultati,
    pub mostra_e: [bool; 2],
    pub mostra_surfer: bool,
    /// the surfer: one ball (positions of all its electrons), at rest
    pub surfer: Vec<[f64; 3]>,
    pub vibra: bool,
    pub t_vib: f64,
    /// running average of the live energy estimate
    pub e_media: f64,
    pub panel_open: bool,
    /// `--giro`: visit every system for 3 s (with vibration where possible), for unattended checks
    pub giro: Option<f64>,
    seed: u64,
    ticks_since_switch: u64,
}

impl Lab {
    fn new(s: Scelta) -> Self {
        let res = Risultati::load(misura::FILE);
        let dmc = scena::nuova(s, &res, BALLS_LIVE, 1);
        Lab { scelta: s, dmc, res, mostra_e: [true, true], mostra_surfer: true, surfer: vec![], vibra: false,
              t_vib: 0.0, e_media: 0.0, panel_open: false, giro: None, seed: 1, ticks_since_switch: 0 }
    }

    pub fn cambia(&mut self, s: Scelta) {
        self.seed += 1;
        self.scelta = s;
        self.dmc = scena::nuova(s, &self.res, BALLS_LIVE, self.seed);
        self.surfer.clear();
        self.vibra = false;
        self.e_media = 0.0;
        self.ticks_since_switch = 0;
    }

    /// A new surfer: one ball drawn at random from the cloud.
    pub fn nuovo_surfer(&mut self) {
        let n = self.dmc.n();
        if n == 0 { return; }
        self.seed += 1;
        let i = (dmc::Rng::new(self.seed).next_u64() % n as u64) as usize;
        let ne = self.dmc.sys.ne();
        self.surfer = (0..ne).map(|e| self.dmc.elettrone(e).nth(i).unwrap()).collect();
    }

    /// Run the rule for (at most) the frame budget.
    fn avanza(&mut self, dt: f64) {
        if let Some(t) = self.giro.as_mut() {
            *t += dt;
            if *t > 3.0 {
                *t = 0.0;
                let all = Scelta::TUTTE;
                let i = all.iter().position(|s| *s == self.scelta).unwrap_or(0);
                self.cambia(all[(i + 1) % all.len()]);
                self.vibra = self.scelta.biatomica();
                self.panel_open = !self.panel_open;
            }
        }
        if self.vibra && self.scelta.biatomica() {
            if let Some(w) = self.scelta.omega(&self.res) {
                // amplitude of the zero-point motion: √(ħ/μω) (atomic units)
                let mu = (costanti::proton_me() + 1.0) / 2.0;
                let w_au = w / costanti::hartree_cm1();
                let amp = (1.0 / (mu * w_au)).sqrt();
                self.t_vib += dt;
                let r = self.scelta.r_eq(&self.res) + amp * (std::f64::consts::TAU * self.t_vib / VIB_PERIOD_S).sin();
                self.dmc.sys = self.scelta.sistema(&self.res, r);
                self.dmc.ricalcola();
            }
        }
        let t = Instant::now();
        for _ in 0..MAX_STEPS_PER_FRAME {
            let (em, _) = self.dmc.passo();
            self.ticks_since_switch += 1;
            let k = if self.ticks_since_switch < 50 { 1.0 } else { 0.01 };
            self.e_media += k * (em - self.e_media);
            if t.elapsed().as_secs_f64() > STEP_BUDGET_S { break; }
        }
        if self.surfer.is_empty() && self.ticks_since_switch > 300 { self.nuovo_surfer(); }
    }
}

/// True while the pointer is over the HUD, so the camera ignores it.
#[derive(Resource, Default)]
struct PointerOverHud(bool);

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Some(i) = args.iter().position(|a| a == "--misura") {
        misura::main(args.get(i + 1).map(|s| s.as_str()).unwrap_or("tutto"));
        return;
    }
    let start = match args.iter().position(|a| a == "--sistema").and_then(|i| args.get(i + 1)).map(|s| s.as_str()) {
        Some("h2") => Scelta::H2,
        Some("h2piu") => Scelta::H2piu,
        Some("he") => Scelta::He,
        Some("2p") => Scelta::H2p,
        _ => Scelta::H1s,
    };
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Molecole — palline-onda degli elettroni".into(),
                resolution: WindowResolution::new(1600.0, 900.0),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        .insert_resource(ClearColor(Color::srgb(0.015, 0.02, 0.035)))
        .insert_resource({
            let mut lab = Lab::new(start);
            if args.iter().any(|a| a == "--giro") { lab.giro = Some(0.0); }
            lab
        })
        .insert_resource(OrbitState { yaw: 0.6, pitch: 0.3, radius: CAM_RADIUS_DEFAULT, focus: Vec3::ZERO })
        .init_resource::<PointerOverHud>()
        .add_systems(Startup, (setup_scene, render::setup_render))
        .add_systems(Update, (
            run_lab,
            render::aggiorna_nuvola.after(run_lab),
            render::disegna.after(run_lab),
            ui::draw_hud,
            hud_pointer.after(ui::draw_hud),
            orbit_camera.run_if(|p: Res<PointerOverHud>| !p.0),
        ))
        .run();
}

fn setup_scene(mut commands: Commands) {
    commands.spawn((Camera3d::default(), Transform::from_xyz(8.0, 4.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y)));
}

fn hud_pointer(mut ctx: EguiContexts, mut over: ResMut<PointerOverHud>) {
    let c = ctx.ctx_mut();
    over.0 = c.wants_pointer_input() || c.is_pointer_over_area();
}

fn run_lab(time: Res<Time>, mut lab: ResMut<Lab>) {
    lab.avanza(time.delta_secs_f64());
}

//! SCONTRO — one photon meets one helium atom (2 electrons), everything computed by the electron balls.
//!
//! A separate small binary (not a mode of `molecole`): it needs the photon's ball renderer
//! (src/photon/balls.wgsl) next to the molecule level (src/molecole: the ball rule dmc.rs, the helium
//! measurements elio.rs and their results file), and keeping it apart leaves both viewers untouched.
//!
//! Run:   cargo run --bin scontro   [--seme N] [--lambda nm]   (default: the computed resonance, 58.437 nm)
//! Tests: cargo test --bin scontro
//! Needs the numbers of `molecole --misura elio1p` and `--misura eliocampo` (assets/molecole/risultati.txt).

use bevy::{
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::{EguiContexts, EguiPlugin};

#[allow(dead_code)]
#[path = "../molecole/costanti.rs"]
mod costanti;
#[allow(dead_code)]
#[path = "../molecole/dmc.rs"]
mod dmc;
#[allow(dead_code)]
#[path = "../molecole/misura.rs"]
mod misura;
#[allow(dead_code)]
#[path = "../molecole/elio.rs"]
mod elio;
#[path = "../rendering/camera.rs"]
mod camera;
mod config;
mod fisica;
mod render;
mod ui;
#[cfg(test)]
mod tests;

use camera::{orbit_camera, OrbitState};
use config::*;

/// Default seed of the event's random numbers (photon taken or not, waiting time, emission direction).
/// With seed 1 the first photon is taken (P = 0.30 per photon): the opening shows a whole cycle.
pub const SEME: u64 = 1;

#[derive(Resource)]
pub struct Lab {
    pub sc: fisica::Scontro,
    pub res: misura::Risultati,
    pub mostra_e: [bool; 2],
    pub mostra_surfer: bool,
    pub panel_open: bool,
}

#[derive(Resource, Default)]
struct PointerOverHud(bool);

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seme = args.iter().position(|a| a == "--seme").and_then(|i| args.get(i + 1)).and_then(|s| s.parse().ok()).unwrap_or(SEME);
    let res = misura::Risultati::load(misura::FILE);
    let mut sc = fisica::Scontro::new(&res, seme);
    if let Some(l) = args.iter().position(|a| a == "--lambda").and_then(|i| args.get(i + 1)).and_then(|s| s.parse::<f64>().ok()) {
        sc.lambda_nm = l;
    }
    eprintln!("lambda = {:.4} nm: sigma = {:.3e} m^2, P per photon = {:.3e} (line computed at {:.4} nm)",
              sc.lambda_nm, sc.luce.sigma(sc.lambda_nm), sc.prob(), sc.luce.lambda_nm);
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Scontro — un fotone e un atomo di elio".into(),
                resolution: WindowResolution::new(1600.0, 900.0),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        .add_plugins(render::BallsPlugin)
        .insert_resource(ClearColor(Color::srgb(0.012, 0.015, 0.03)))
        .insert_resource(Lab { sc, res, mostra_e: [true, true], mostra_surfer: true, panel_open: false })
        .insert_resource(OrbitState { yaw: -std::f32::consts::FRAC_PI_2, pitch: 0.25, radius: CAM_RADIUS_DEFAULT, focus: Vec3::new(0.0, 0.0, 4.0) })
        .init_resource::<PointerOverHud>()
        .add_systems(Startup, (setup_scene, render::setup_render))
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
    commands.spawn((Camera3d::default(), Transform::from_xyz(-80.0, 20.0, 30.0).looking_at(Vec3::ZERO, Vec3::Y)));
}

fn hud_pointer(mut ctx: EguiContexts, mut over: ResMut<PointerOverHud>) {
    let c = ctx.ctx_mut();
    over.0 = c.wants_pointer_input() || c.is_pointer_over_area();
}

fn run_lab(time: Res<Time>, mut lab: ResMut<Lab>) {
    // frames longer than 0.1 s (window dragged, start-up) do not make the event jump
    let dt = time.delta_secs_f64().min(0.1);
    lab.sc.avanza(dt);
}

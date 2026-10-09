//! Photon laboratory — one photon made of tens of thousands of grains.
//!
//! Each grain only translates at c and rotates its phase (2π per λ);
//! grains see each other by comparing phases where they cross. From that:
//! propagation, diffraction, polarization, single-photon interference and,
//! for pairs, real polarization entanglement. See `physics.rs`.
//!
//! Three rules move the photon (hidden panel, or `--regola base|pilota|classica`):
//!   • base 3-D (default for free flight, double slit, polarizers, Airy disk,
//!     Arago spot): the base rule v2 of docs/REGOLA_BASE.md in 3-D
//!     (`base3d.rs`) — wave balls + one photon ball guided by the 3-D
//!     current, drawn as a cloud of balls (`render3d.rs`);
//!   • pilot: the 2-D "two roles" rule of `pilot.rs`;
//!   • classic: grains + photon-level Born decisions (`physics.rs`).
//!
//! Run:  cargo run --bin fotone [-- --esperimento fenditure --regola base --fps]
//! Test: cargo test --bin fotone -- --nocapture --test-threads=1

use bevy::{
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::{EguiContexts, EguiPlugin};

mod base3d;
mod complex;
mod config;
mod entangle;
mod gpu;
mod haze;
mod physics;
mod pilot;
mod render;
mod render3d;
mod timeline;
mod ui;
#[path = "../rendering/camera.rs"]
mod camera;
#[cfg(test)]
mod validation;
#[cfg(test)]
mod pilot_validation;
#[cfg(test)]
mod base3d_validation;

use camera::{orbit_camera, OrbitState};
use physics::{Experiment, Lab, Rule, Sim};
use timeline::TimeCtl;

/// The running laboratory.
#[derive(Resource)]
pub struct LabState {
    pub sim: Sim,
    /// Emit a new photon as soon as the previous one is gone.
    pub auto: bool,
    pub show_arrows: bool,
    /// Reversible time: mode, speed factor, history (see `timeline.rs`).
    pub time: TimeCtl,
    /// Temperature of the environment, K (visual only for now: `haze.rs`).
    pub temperature: f32,
    /// The settings panel (Tab / ☰), hidden by default.
    pub panel_open: bool,
    /// Base 3-D rule: wave balls drawn per fan (a random subset of the rule's).
    pub balls_per_fan: usize,
    /// Size of the drawn balls (× the default).
    pub ball_scale: f32,
    /// Draw the wave balls on the rule's exact direction grid (else spread
    /// over each direction cell: visual only, see `render3d.rs`).
    pub exact_dirs: bool,
    /// The camera should frame the bench again (experiment switched).
    pub camera_dirty: bool,
}

impl LabState {
    pub fn new(e: Experiment) -> Self {
        Self::from_lab(Lab::default_for(e))
    }

    pub fn from_lab(lab: Lab) -> Self {
        let sim = Sim::new(lab);
        let time = TimeCtl::new(&sim);
        Self {
            sim, auto: true, show_arrows: false, time, temperature: config::T_DEFAULT, panel_open: false,
            balls_per_fan: config::BALLS_PER_FAN, ball_scale: 1.0, exact_dirs: false, camera_dirty: true,
        }
    }

    /// Switch experiment, keeping the viewer's preferences.
    pub fn switch(&mut self, e: Experiment) {
        let mut lab = Lab::default_for(e);
        // keep the chosen rule (Bell has only the classic one)
        if e != Experiment::Bell && self.sim.lab.experiment != Experiment::Bell {
            lab.rule = self.sim.lab.rule;
        }
        let mut next = LabState::from_lab(lab);
        next.auto = self.auto;
        next.show_arrows = self.show_arrows;
        next.temperature = self.temperature;
        next.panel_open = self.panel_open;
        next.balls_per_fan = self.balls_per_fan;
        next.ball_scale = self.ball_scale;
        next.exact_dirs = self.exact_dirs;
        next.time.factor = self.time.factor;
        next.time.mode = self.time.mode;
        *self = next;
    }
}

/// True while the pointer is over the HUD, so the camera ignores it.
#[derive(Resource, Default)]
struct PointerOverHud(bool);

fn main() {
    let fps = std::env::args().any(|a| a == "--fps");
    let mut app = App::new();
    app
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title:        "Photon Lab".into(),
                resolution:   WindowResolution::new(1600.0, 900.0),
                position:     WindowPosition::Centered(MonitorSelection::Primary),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        .add_plugins(render3d::BallsPlugin)
        .insert_resource(ClearColor(Color::srgb(0.015, 0.02, 0.035)))
        .insert_resource(Script::from_args())
        .insert_resource(LabState::from_lab(lab_from_args()))
        .insert_resource(OrbitState { yaw: 1.15, pitch: 0.35, radius: config::CAM_RADIUS_DEFAULT, focus: Vec3::new(0.0, 0.0, 60.0) })
        .init_resource::<PointerOverHud>()
        .insert_gizmo_config(render::PhotonBallGizmos, render::photon_ball_gizmo_config())
        .add_systems(Startup, (setup_scene, render::setup_render, haze::setup_haze))
        .add_systems(Update, (
            run_lab,
            render::update_meshes.after(run_lab),
            render3d::update_visuals.after(run_lab),
            frame_camera.before(orbit_camera),
            render::draw_bench,
            render::draw_photon_ball,
            haze::update_haze.after(run_lab),
            ui::draw_hud,
            hud_pointer.after(ui::draw_hud),
            orbit_camera.run_if(|p: Res<PointerOverHud>| !p.0),
        ))
        ;
    if fps {
        // `--fps`: average frame rate on stderr every 5 s (diagnostics only)
        app.add_plugins(bevy::diagnostic::FrameTimeDiagnosticsPlugin)
            .add_systems(Update, report_fps);
    }
    app.run();
}

fn report_fps(diag: Res<bevy::diagnostic::DiagnosticsStore>, time: Res<Time>, mut last: Local<f32>, state: Res<LabState>) {
    let t = time.elapsed_secs();
    if t - *last < 5.0 { return; }
    *last = t;
    let fps = diag.get(&bevy::diagnostic::FrameTimeDiagnosticsPlugin::FPS).and_then(|d| d.average());
    let balls: usize = state.sim.photons.iter().map(|p| p.grains.len()).sum();
    let rule = if state.sim.lab.base3d() { Rule::Base3D } else { state.sim.lab.rule };
    eprintln!("[fps] t={t:.0}s fps={:.1} rule={rule:?} balls_in_flight={balls} emitted={} detected={} clock={:.1}",
        fps.unwrap_or(0.0), state.sim.stats.emitted, state.sim.stats.detected, state.sim.clock);
}

/// Command line: `--esperimento libero|fenditure|polarizzatori|entanglement|airy|arago`,
/// `--regola base|pilota|classica`, `--palline N`, `--polarizzazione h|v|d|a|c+|c-`,
/// `--fps` (frame rate on stderr).
fn lab_from_args() -> Lab {
    let args: Vec<String> = std::env::args().collect();
    let value = |name: &str| args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).cloned();
    let experiment = match value("--esperimento").as_deref() {
        Some("libero") => Experiment::FreeBeam,
        Some("polarizzatori") => Experiment::Polarizers,
        Some("entanglement") => Experiment::Bell,
        Some("airy") => Experiment::Airy,
        Some("arago") => Experiment::Arago,
        _ => Experiment::DoubleSlit,
    };
    let mut lab = Lab::default_for(experiment);
    match value("--regola").as_deref() {
        Some("classica") => lab = Lab::preset(experiment),   // the classic laboratory, unchanged
        Some("pilota") if experiment != Experiment::Bell => lab.rule = Rule::Pilot,
        Some("base") if experiment != Experiment::Bell => lab.rule = Rule::Base3D,
        _ => {}
    }
    if let Some(n) = value("--palline").and_then(|v| v.parse::<usize>().ok()) {
        lab.grains = n.clamp(100, 4_000_000);
    }
    if let Some(p) = value("--polarizzazione") {
        use physics::Polarization::*;
        lab.polarization = match p.as_str() {
            "h" => H, "v" => V, "d" => D, "a" => A, "c-" => CircMinus, _ => CircPlus,
        };
    }
    lab
}

fn setup_scene(mut commands: Commands) {
    commands.spawn((Camera3d::default(), Transform::from_xyz(200.0, 80.0, 60.0).looking_at(Vec3::new(0.0, 0.0, 60.0), Vec3::Y)));
}

/// Frame the bench when the experiment changes.
fn frame_camera(mut state: ResMut<LabState>, mut orbit: ResMut<OrbitState>) {
    if !state.camera_dirty { return; }
    state.camera_dirty = false;
    let (focus, radius, yaw, pitch) = render3d::camera_for(&state.sim.lab);
    orbit.focus = focus;
    orbit.radius = radius;
    orbit.yaw = yaw;
    orbit.pitch = pitch;
}

fn hud_pointer(mut ctx: EguiContexts, mut over: ResMut<PointerOverHud>) {
    let c = ctx.ctx_mut();
    over.0 = c.wants_pointer_input() || c.is_pointer_over_area();
}

/// Advance (or rewind) the laboratory by one frame.
fn run_lab(time: Res<Time>, mut state: ResMut<LabState>, mut script: ResMut<Script>) {
    script.before(&mut state);
    state.advance(time.delta_secs());
    script.after(&mut state);
}

/// Command-line script for demos and screenshots (none by default):
///   `--clic N`     run at ×1000 until N clicks, then go on at ×1;
///   `--ferma-a CT` pause when the photon in flight reaches light-clock CT
///                  (after the `--clic` phase, if any).
#[derive(Resource, Default)]
struct Script {
    clicks: Option<u64>,
    stop_at: Option<f32>,
}

impl Script {
    fn from_args() -> Self {
        let args: Vec<String> = std::env::args().collect();
        let value = |name: &str| args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).cloned();
        Self {
            clicks: value("--clic").and_then(|v| v.parse().ok()),
            stop_at: value("--ferma-a").and_then(|v| v.parse().ok()),
        }
    }
    fn before(&mut self, state: &mut LabState) {
        if let Some(n) = self.clicks {
            if state.sim.stats.detected >= n {
                state.time.factor = 1.0;
                self.clicks = None;
            } else {
                state.time.factor = config::SPEED_MAX;
            }
        }
    }
    fn after(&mut self, state: &mut LabState) {
        if self.clicks.is_some() { return; }
        if let Some(ct) = self.stop_at {
            if state.sim.in_flight() && state.sim.clock >= ct {
                state.time.mode = timeline::Mode::Pause;
                self.stop_at = None;
            }
        }
    }
}

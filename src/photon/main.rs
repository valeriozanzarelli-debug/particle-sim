//! Photon laboratory — one photon made of tens of thousands of grains.
//!
//! Each grain only translates at c and rotates its phase (2π per λ);
//! grains see each other by comparing phases where they cross. From that:
//! propagation, diffraction, polarization, single-photon interference and,
//! for pairs, real polarization entanglement. See `physics.rs`.
//!
//! Run:  cargo run --release --bin fotone
//! Test: cargo test --release --bin fotone -- --nocapture --test-threads=1

use bevy::{
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::{EguiContexts, EguiPlugin};

mod complex;
mod config;
mod entangle;
mod gpu;
mod haze;
mod physics;
mod render;
mod timeline;
mod ui;
#[path = "../rendering/camera.rs"]
mod camera;
#[cfg(test)]
mod validation;

use camera::{orbit_camera, OrbitState};
use physics::{Experiment, Lab, Sim};
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
}

impl LabState {
    pub fn new(e: Experiment) -> Self {
        let sim = Sim::new(Lab::preset(e));
        let time = TimeCtl::new(&sim);
        Self { sim, auto: true, show_arrows: true, time, temperature: config::T_DEFAULT, panel_open: false }
    }

    /// Switch experiment, keeping the viewer's preferences.
    pub fn switch(&mut self, e: Experiment) {
        let mut next = LabState::new(e);
        next.auto = self.auto;
        next.show_arrows = self.show_arrows;
        next.temperature = self.temperature;
        next.panel_open = self.panel_open;
        next.time.factor = self.time.factor;
        next.time.mode = self.time.mode;
        *self = next;
    }
}

/// True while the pointer is over the HUD, so the camera ignores it.
#[derive(Resource, Default)]
struct PointerOverHud(bool);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title:        "Photon Lab — a photon made of grains".into(),
                resolution:   WindowResolution::new(1600.0, 900.0),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        .insert_resource(ClearColor(Color::srgb(0.015, 0.02, 0.035)))
        .insert_resource(LabState::new(Experiment::DoubleSlit))
        .insert_resource(OrbitState { yaw: 1.15, pitch: 0.35, radius: config::CAM_RADIUS_DEFAULT, focus: Vec3::new(0.0, 0.0, 60.0) })
        .init_resource::<PointerOverHud>()
        .add_systems(Startup, (setup_scene, render::setup_render, haze::setup_haze))
        .add_systems(Update, (
            run_lab,
            render::update_meshes.after(run_lab),
            render::draw_bench,
            haze::update_haze.after(run_lab),
            ui::draw_hud,
            hud_pointer.after(ui::draw_hud),
            orbit_camera.run_if(|p: Res<PointerOverHud>| !p.0),
        ))
        .run();
}

fn setup_scene(mut commands: Commands) {
    commands.spawn((Camera3d::default(), Transform::from_xyz(200.0, 80.0, 60.0).looking_at(Vec3::new(0.0, 0.0, 60.0), Vec3::Y)));
}

fn hud_pointer(mut ctx: EguiContexts, mut over: ResMut<PointerOverHud>) {
    let c = ctx.ctx_mut();
    over.0 = c.wants_pointer_input() || c.is_pointer_over_area();
}

/// Advance (or rewind) the laboratory by one frame.
fn run_lab(time: Res<Time>, mut state: ResMut<LabState>) {
    state.advance(time.delta_secs());
}

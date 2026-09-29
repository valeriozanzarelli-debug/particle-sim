//! Emergent Particle Simulation — main entry point.
//!
//! Concept: millions of identical dimensionless particles inside a 3-D toroidal
//! cube. Each moves forward at constant speed and turns according to the 3-D
//! Primordial Particle System (PPS) rule — no imposed physics. Structure
//! (clusters, orbits, stable bound states) emerges purely from the local rule.
//!
//! Controls (HUD): α, β, R, speed, particle count, steps/frame, camera.
//!
//! Compute backend: NVIDIA GPU via wgpu compute shaders (Vulkan/DX12).
//! Falls back to CPU rayon if no compatible GPU is found.

use bevy::{
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::EguiPlugin;

mod config;
mod simulation;
mod rendering;

use config::CAM_RADIUS_DEFAULT;
use simulation::{
    gpu_pps::GpuPpsPlugin,
    pps::step_simulation,
    world::ParticleWorld,
};
use rendering::{
    camera::{orbit_camera, OrbitState},
    point_cloud::{setup_point_cloud, update_point_cloud},
    ui::{draw_hud, SimStats},
};

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title:        "Emergent Particle Simulation — PPS 3D".into(),
                resolution:   WindowResolution::new(1600.0, 900.0),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        // ── GPU compute plugin (tries NVIDIA first, falls back to CPU) ──────
        .add_plugins(GpuPpsPlugin)
        // ── Resources ──────────────────────────────────────────────────────
        .init_resource::<ParticleWorld>()
        .init_resource::<OrbitState>()
        .init_resource::<SimStats>()
        // ── Startup ─────────────────────────────────────────────────────────
        .add_systems(Startup, (setup_scene, setup_point_cloud).chain())
        // ── Update ──────────────────────────────────────────────────────────
        .add_systems(Update, (
            // CPU fallback: only runs when GpuCompute resource is NOT present
            step_simulation
                .run_if(not(resource_exists::<simulation::gpu_pps::GpuCompute>)),
            update_point_cloud,
            draw_box_gizmo,
            orbit_camera,
            draw_hud,
        ))
        .run();
}

// ─── One-time scene setup ──────────────────────────────────────────────────────

fn setup_scene(mut commands: Commands) {
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(0.0, 40.0, CAM_RADIUS_DEFAULT)
            .looking_at(Vec3::ZERO, Vec3::Y),
    ));

    commands.insert_resource(AmbientLight {
        color:      Color::WHITE,
        brightness: 300.0,
    });
}

// ─── Simulation boundary cube (gizmos, zero draw-call overhead) ───────────────

fn draw_box_gizmo(mut gizmos: Gizmos) {
    use config::BOX_HALF;
    let h = BOX_HALF;
    let color = bevy::color::palettes::css::DIM_GRAY;
    let corners = [
        Vec3::new(-h, -h, -h), Vec3::new( h, -h, -h),
        Vec3::new( h,  h, -h), Vec3::new(-h,  h, -h),
        Vec3::new(-h, -h,  h), Vec3::new( h, -h,  h),
        Vec3::new( h,  h,  h), Vec3::new(-h,  h,  h),
    ];
    let edges = [
        (0,1),(1,2),(2,3),(3,0),
        (4,5),(5,6),(6,7),(7,4),
        (0,4),(1,5),(2,6),(3,7),
    ];
    for (a, b) in edges {
        gizmos.line(corners[a], corners[b], color);
    }
}

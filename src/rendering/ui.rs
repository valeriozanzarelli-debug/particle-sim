//! egui HUD: controls for the emergent particle simulation.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};

use crate::config::*;
use crate::simulation::world::ParticleWorld;
use crate::simulation::gpu_pps::GpuCompute;
use crate::rendering::camera::OrbitState;
use crate::rendering::point_cloud::density_color;

/// How many ticks to cache for a rough FPS-per-simulation-step estimate.
#[derive(Resource, Default)]
pub struct SimStats {
    pub fps: f32,
}

pub fn draw_hud(
    mut ctx:  EguiContexts,
    mut world: ResMut<ParticleWorld>,
    mut orbit: ResMut<OrbitState>,
    mut stats: ResMut<SimStats>,
    time:      Res<Time>,
    gpu:       Option<Res<GpuCompute>>,
) {
    stats.fps = 1.0 / time.delta_secs().max(1e-5);

    egui::Window::new("⚛  Emergent Particles")
        .default_width(290.0)
        .resizable(false)
        .show(ctx.ctx_mut(), |ui| {
            // ── Compute backend badge ────────────────────────────────────
            ui.horizontal(|ui| {
                match &gpu {
                    Some(g) => {
                        ui.colored_label(
                            egui::Color32::from_rgb(100, 220, 100),
                            format!("GPU  {}", g.gpu_name()),
                        );
                    }
                    None => {
                        ui.colored_label(
                            egui::Color32::YELLOW,
                            "CPU (rayon)",
                        );
                    }
                }
            });

            ui.add_space(2.0);

            // ── Status ──────────────────────────────────────────────────────
            ui.horizontal(|ui| {
                ui.label("FPS:");
                ui.colored_label(fps_color(stats.fps), format!("{:.0}", stats.fps));
                ui.separator();
                ui.label(format!("Tick: {}", world.ticks));
            });
            ui.horizontal(|ui| {
                ui.label(format!("N: {}", world.n));
                ui.separator();
                ui.label(format!("⟨nbr⟩: {:.1}", world.mean_neighbors));
                ui.separator();
                ui.label(format!("max: {}", world.max_cluster));
            });

            ui.add_space(4.0);
            ui.separator();

            // ── Steps per frame (time scale) ─────────────────────────────
            ui.label("Steps / frame  (time speed)");
            let mut spf = world.steps_per_frame as f32;
            ui.add(egui::Slider::new(&mut spf, 0.0..=20.0).step_by(1.0).suffix("×"));
            world.steps_per_frame = spf.round() as u32;

            ui.separator();

            // ── PPS rule parameters ───────────────────────────────────────
            ui.collapsing("⚙  Rule parameters", |ui| {
                let mut alpha = world.alpha;
                let mut beta  = world.beta;
                let mut r     = world.r_neighbor;
                let mut speed = world.speed;

                ui.label("α  base turn (rad)");
                ui.add(egui::Slider::new(&mut alpha, 0.0..=1.0).step_by(0.005));
                ui.label("β  neighbour influence (rad/nbr)");
                ui.add(egui::Slider::new(&mut beta,  0.0..=0.05).step_by(0.0005));
                ui.label("R  detection radius");
                ui.add(egui::Slider::new(&mut r,     0.5..=12.0).step_by(0.1));
                ui.label("v  speed (sim/tick)");
                ui.add(egui::Slider::new(&mut speed, 0.05..=2.0).step_by(0.01));

                world.alpha      = alpha;
                world.beta       = beta;
                world.r_neighbor = r;
                world.speed      = speed;
            });

            ui.separator();

            // ── Particle count ────────────────────────────────────────────
            ui.collapsing("🔢  Particle count  (restart)", |ui| {
                ui.horizontal_wrapped(|ui| {
                    for &n in &[10_000usize, 50_000, 100_000, 200_000, 500_000] {
                        if ui.button(format!("{}K", n / 1000)).clicked() {
                            world.reset(n);
                        }
                    }
                });
                ui.label(egui::RichText::new(
                    "⚠  Large counts need a fast CPU.\n   rayon uses all cores.",
                ).small().weak());
            });

            ui.separator();

            // ── Reset button ──────────────────────────────────────────────
            if ui.button("↺  Reset positions (keep params)").clicked() {
                let n = world.n;
                world.reset(n);
            }

            ui.separator();

            // ── Camera controls hint ──────────────────────────────────────
            ui.collapsing("🎥  Camera", |ui| {
                ui.label("LMB drag  → rotate");
                ui.label("RMB drag  → pan");
                ui.label("Scroll    → zoom");
                ui.separator();
                let mut r = orbit.radius;
                ui.add(egui::Slider::new(&mut r, CAM_RADIUS_MIN..=CAM_RADIUS_MAX)
                    .text("radius"));
                orbit.radius = r;
                if ui.button("Reset camera").clicked() {
                    *orbit = OrbitState::default();
                }
            });

            ui.separator();

            // ── Colour legend ─────────────────────────────────────────────
            ui.collapsing("🎨  Colour = neighbour density", |ui| {
                let steps = 10usize;
                ui.horizontal(|ui| {
                    for i in 0..=steps {
                        let n_val = (i as f32 / steps as f32 * COLOR_SCALE_N) as u8;
                        let [r,g,b,_] = density_color(n_val);
                        let c = egui::Color32::from_rgb(
                            (r*255.0) as u8, (g*255.0) as u8, (b*255.0) as u8);
                        let (rect, _) = ui.allocate_exact_size(
                            egui::vec2(18.0, 14.0), egui::Sense::hover());
                        ui.painter().rect_filled(rect, 2.0, c);
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("0");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(format!("{:.0}+ nbrs", COLOR_SCALE_N));
                    });
                });
                ui.label(egui::RichText::new(
                    "Dense (white/red) = bound cluster.\nSparse (blue) = free particle.",
                ).small().weak());
            });
        });
}

fn fps_color(fps: f32) -> egui::Color32 {
    if fps >= 55.0 { egui::Color32::GREEN }
    else if fps >= 30.0 { egui::Color32::YELLOW }
    else { egui::Color32::RED }
}

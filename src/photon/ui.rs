//! HUD of the photon laboratory: controls, live plots, live validation.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::Color32;

use crate::config::NM_PER_UNIT;
use crate::entangle::{correlation_local, correlation_quantum};
use crate::physics::*;
use crate::LabState;

const DEG: f32 = std::f32::consts::PI / 180.0;
const CYAN:  Color32 = Color32::from_rgb(90, 220, 255);
const GREEN: Color32 = Color32::from_rgb(170, 255, 140);
const GOLD:  Color32 = Color32::from_rgb(255, 205, 80);
const GREY:  Color32 = Color32::from_rgb(150, 155, 170);
const ORANGE: Color32 = Color32::from_rgb(255, 150, 80);

/// Exact double-slit pattern, recomputed only when the geometry changes.
#[derive(Default)]
pub struct RefCache {
    key: String,
    profile: Vec<f32>,
    period: f32,
    visibility: f32,
}

pub fn draw_hud(mut ctx: EguiContexts, mut state: ResMut<LabState>, mut cache: Local<RefCache>) {
    let state = &mut *state;
    egui::Window::new("Photon Lab")
        .default_width(340.0)
        .default_pos([12.0, 12.0])
        .show(ctx.ctx_mut(), |ui| {
            egui::ScrollArea::vertical().max_height(860.0).show(ui, |ui| {
                experiment_picker(ui, state);
                ui.separator();
                controls(ui, state);
                ui.separator();
                photon_params(ui, state);
                ui.separator();
                match state.sim.lab.experiment {
                    Experiment::FreeBeam   => free_beam(ui, state),
                    Experiment::DoubleSlit => double_slit(ui, state, &mut cache),
                    Experiment::Polarizers => polarizers(ui, state),
                    Experiment::Bell       => bell(ui, state),
                }
                ui.separator();
                legend(ui);
            });
        });
}

fn experiment_picker(ui: &mut egui::Ui, state: &mut LabState) {
    let current = state.sim.lab.experiment;
    let mut pick = current;
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut pick, Experiment::FreeBeam,   "Free flight");
        ui.selectable_value(&mut pick, Experiment::DoubleSlit, "Double slit");
        ui.selectable_value(&mut pick, Experiment::Polarizers, "Polarizers");
        ui.selectable_value(&mut pick, Experiment::Bell,       "Entanglement");
    });
    if pick != current {
        let (auto, speed, arrows) = (state.auto, state.speed, state.show_arrows);
        *state = LabState::new(pick);
        state.auto = auto;
        state.speed = speed;
        state.show_arrows = arrows;
    }
}

fn controls(ui: &mut egui::Ui, state: &mut LabState) {
    ui.horizontal(|ui| {
        let busy = state.sim.in_flight();
        if ui.add_enabled(!busy, egui::Button::new("▶ Emit one photon")).clicked() {
            state.sim.launch();
        }
        ui.checkbox(&mut state.auto, "continuous");
        if ui.button("↺ Clear counts").clicked() {
            state.sim.reset_stats();
        }
    });
    ui.label("Slow motion: light moves this far per frame");
    ui.add(egui::Slider::new(&mut state.speed, 0.05..=20.0).logarithmic(true).suffix(" λ-units"));
    ui.label("Fast count: whole photons per frame (no animation)");
    ui.add(egui::Slider::new(&mut state.turbo, 0..=200));
    ui.checkbox(&mut state.show_arrows, "show the E-field arrows (the rotation)");

    let s = &state.sim.stats;
    ui.label(egui::RichText::new(format!(
        "emitted {}   detected {}   absorbed {}   off-screen {}",
        s.emitted, s.detected, s.absorbed_wall + s.absorbed_pol.iter().sum::<u64>(), s.missed,
    )).monospace().small());
}

fn photon_params(ui: &mut egui::Ui, state: &mut LabState) {
    let lab = &mut state.sim.lab;
    ui.collapsing("Photon", |ui| {
        ui.label("wavelength λ (applies to the next photon)");
        ui.add(egui::Slider::new(&mut lab.lambda, 0.5..=2.0).step_by(0.05));
        let nm = lab.lambda * NM_PER_UNIT;
        ui.label(egui::RichText::new(format!(
            "scale 1 unit = {NM_PER_UNIT:.0} nm  →  λ = {nm:.0} nm,  ν = {:.0} THz,\nE = hν = {:.2} eV,  p = h/λ = {:.2} eV/c",
            299_792.458 / nm, 1239.84 / nm, 1239.84 / nm,
        )).small().weak());
        ui.label("grains per photon");
        let mut g = lab.grains as f32;
        ui.add(egui::Slider::new(&mut g, 1_000.0..=120_000.0).logarithmic(true));
        lab.grains = g as usize;
        ui.label("beam waist w0");
        ui.add(egui::Slider::new(&mut lab.waist, 1.5..=12.0).step_by(0.1));
        if lab.experiment != Experiment::Bell {
            egui::ComboBox::from_label("polarization")
                .selected_text(lab.polarization.label())
                .show_ui(ui, |ui| {
                    for p in Polarization::ALL {
                        ui.selectable_value(&mut lab.polarization, p, p.label());
                    }
                });
        }
    });
}

// ─── Experiments ──────────────────────────────────────────────────────────────

fn free_beam(ui: &mut egui::Ui, state: &mut LabState) {
    let sim = &state.sim;
    let lab = &sim.lab;
    ui.heading("Free flight");
    ui.label(egui::RichText::new(
        "The photon flies alone. Every grain moves at c and turns its arrow\n\
         once per λ: with circular light the arrows form a helix.",
    ).small().weak());

    // Live width of the photon in flight vs Gaussian-beam theory
    if let Some(ph) = sim.photons.first() {
        let alive: Vec<_> = ph.grains.iter().filter(|g| g.alive).collect();
        if alive.len() > 100 {
            let n = alive.len() as f32;
            let mx = alive.iter().map(|g| g.pos.x).sum::<f32>() / n;
            let mz = alive.iter().map(|g| g.pos.z).sum::<f32>() / n;
            let sx = (alive.iter().map(|g| (g.pos.x - mx).powi(2)).sum::<f32>() / n).sqrt();
            check(ui, "beam radius w(z)", 2.0 * sx, lab.beam_radius(mz), 0.03,
                &format!("at z = {mz:.0}, z_R = {:.1}", lab.rayleigh()));
        }
    }
    let w = lab.beam_radius(lab.screen_z);
    let theory: Vec<f32> = (0..lab.bins).map(|b| (-2.0 * lab.bin_x(b).powi(2) / (w * w)).exp()).collect();
    let hist: Vec<f32> = sim.stats.hist.iter().map(|&h| h as f32).collect();
    ui.label(format!("screen: {} clicks vs Gaussian profile w = {w:.1}", sim.stats.detected));
    plot(ui, &[(&theory, GREY, false), (&hist, GREEN, true)]);
}

fn double_slit(ui: &mut egui::Ui, state: &mut LabState, cache: &mut RefCache) {
    ui.heading("Double slit");
    {
        let lab = &mut state.sim.lab;
        ui.horizontal(|ui| {
            ui.label("slit distance d");
            ui.add(egui::Slider::new(&mut lab.slit_sep, 3.0..=16.0).step_by(0.5));
        });
        ui.horizontal(|ui| {
            ui.label("slit width a");
            ui.add(egui::Slider::new(&mut lab.slit_width, 0.5..=4.0).step_by(0.1));
        });
        ui.horizontal(|ui| {
            ui.label("which-path polarizers:");
            ui.selectable_value(&mut lab.which_path, WhichPath::Off, "off");
            ui.selectable_value(&mut lab.which_path, WhichPath::Parallel, "∥ ∥");
            ui.selectable_value(&mut lab.which_path, WhichPath::Orthogonal, "∥ ⊥");
        });
    }
    let sim = &state.sim;
    let lab = &sim.lab;
    let key = format!("{:?}{}{}{}{}{:?}{:?}", lab.which_path, lab.lambda, lab.slit_sep, lab.slit_width, lab.waist, lab.polarization, lab.screen_z);
    if cache.key != key {
        cache.profile = double_slit_reference(lab);
        cache.period = fringe_spacing(lab, &cache.profile).unwrap_or(lab.lambda * (lab.screen_z - lab.wall_z) / lab.slit_sep);
        cache.visibility = visibility_at(lab, &cache.profile, cache.period);
        cache.key = key;
    }
    let n = sim.stats.photons_in_mean.max(1) as f64;
    let mean: Vec<f32> = sim.stats.mean_intensity.iter().map(|v| (v / n) as f32).collect();
    let hist: Vec<f32> = sim.stats.hist.iter().map(|&h| h as f32).collect();

    ui.label(egui::RichText::new(
        "grey = exact wave theory · cyan = wave from the grains · green = photon clicks",
    ).small().weak());
    plot(ui, &[(&cache.profile, GREY, false), (&mean, CYAN, false), (&hist, GREEN, true)]);

    let l = lab.screen_z - lab.wall_z;
    let expected = lab.lambda * l / lab.slit_sep;
    if sim.stats.photons_in_mean >= 5 {
        if let Some(sp) = fringe_spacing(lab, &mean) {
            check(ui, "fringe spacing", sp, expected, 0.03, "λ·L/d");
        }
        let corr = correlation(&mean, &cache.profile);
        row(ui, "match with exact pattern", format!("{:.3}", corr), corr > 0.97);
        let v = visibility_at(lab, &mean, cache.period);
        row(ui, "fringe visibility", format!("{v:.2}  (theory {:.2})", cache.visibility), (v - cache.visibility).abs() < 0.1);
    } else {
        ui.label(egui::RichText::new("need ≥ 5 photons through the slits for the checks").small().weak());
    }
    let finished = sim.stats.detected + sim.stats.absorbed_wall + sim.stats.missed;
    if finished > 0 {
        ui.label(format!(
            "{:.0}% of the photons passed the slits, the rest was absorbed by the wall — whole.",
            100.0 * sim.stats.detected as f32 / finished as f32,
        ));
    }
}

fn polarizers(ui: &mut egui::Ui, state: &mut LabState) {
    ui.heading("Polarizers — Malus");
    {
        let lab = &mut state.sim.lab;
        for (i, p) in lab.polarizers.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.checkbox(&mut p.enabled, format!("P{} at z={:.0}", i + 1, p.z));
                let mut deg = p.angle / DEG;
                ui.add(egui::Slider::new(&mut deg, 0.0..=180.0).step_by(1.0).suffix("°"));
                p.angle = deg * DEG;
            });
        }
    }
    let sim = &state.sim;
    let s = &sim.stats;
    let finished = s.detected + s.absorbed_pol.iter().sum::<u64>() + s.missed;
    let theory = malus_theory(&sim.lab);
    if finished > 0 {
        let m = s.detected as f32 / finished as f32;
        let tol = 3.0 * (theory * (1.0 - theory) / finished as f32).sqrt() + 1e-3;
        row(ui, "photons through", format!("{m:.3}  (theory {theory:.3}, {finished} photons)"), (m - theory).abs() <= tol);
        ui.label(format!("absorbed at P1 {}  P2 {}  P3 {}", s.absorbed_pol[0], s.absorbed_pol[1], s.absorbed_pol[2]));
    }
    ui.label(egui::RichText::new(
        "Each photon passes whole or not at all, with probability cos²θ.\n\
         Try P1 0°, P3 90°: nothing passes. Add P2 at 45°: a quarter passes.",
    ).small().weak());
}

fn bell(ui: &mut egui::Ui, state: &mut LabState) {
    ui.heading("Entangled pair — Bell test");
    {
        let b = &mut state.sim.lab.bell;
        ui.checkbox(&mut b.hidden_variables, "replace entanglement with hidden variables (classical)");
        ui.checkbox(&mut b.random_chsh, "random settings (CHSH: Alice 0°/45°, Bob 22.5°/67.5°)");
        if !b.random_chsh {
            let (mut a, mut bb) = (b.angle_a / DEG, b.angle_b / DEG);
            ui.horizontal(|ui| { ui.label("Alice"); ui.add(egui::Slider::new(&mut a, 0.0..=180.0).suffix("°")); });
            ui.horizontal(|ui| { ui.label("Bob  "); ui.add(egui::Slider::new(&mut bb, 0.0..=180.0).suffix("°")); });
            b.angle_a = a * DEG;
            b.angle_b = bb * DEG;
        }
    }
    let sim = &state.sim;
    let st = &sim.stats.bell;
    let cfg = sim.lab.bell;

    if let Some(pair_state) = sim.photons.first().map(|p| p.pol_defined) {
        ui.label(if pair_state {
            "state: measured — each photon now has its own polarization"
        } else {
            "state: entangled — neither photon has a polarization of its own"
        });
    }
    if let Some((a, pa, b, pb)) = st_last(sim) {
        ui.label(format!("last pair: Alice {:.1}° → {}   Bob {:.1}° → {}",
            a / DEG, if pa { "+" } else { "−" }, b / DEG, if pb { "+" } else { "−" }));
    }

    if cfg.random_chsh {
        egui::Grid::new("chsh").striped(true).show(ui, |ui| {
            ui.label("E(a, b)"); ui.label("measured"); ui.label("quantum"); ui.label("local"); ui.end_row();
            for i in 0..2 { for j in 0..2 {
                ui.label(format!("{:.1}°, {:.1}°", CHSH_A[i] / DEG, CHSH_B[j] / DEG));
                let e = BellStats::correlation(&st.chsh[i][j]);
                ui.label(e.map_or("—".into(), |e| format!("{e:+.3}")));
                ui.label(format!("{:+.3}", correlation_quantum(CHSH_A[i], CHSH_B[j])));
                ui.label(format!("{:+.3}", correlation_local(CHSH_A[i], CHSH_B[j])));
                ui.end_row();
            }}
        });
        if let Some(s) = st.chsh_s() {
            let mut var = 0.0;
            for i in 0..2 { for j in 0..2 {
                let n = st.chsh[i][j].iter().sum::<u64>() as f32;
                let e = BellStats::correlation(&st.chsh[i][j]).unwrap();
                var += (1.0 - e * e) / n;
            }}
            let sigma = var.sqrt();
            ui.label(egui::RichText::new(format!("S = {s:.3} ± {sigma:.3}")).size(20.0).color(if s > 2.0 + 3.0 * sigma { CYAN } else { GOLD }));
            ui.label(format!("quantum limit 2√2 = 2.828 · any local model ≤ 2 · {} pairs", st.pairs()));
            if s > 2.0 + 3.0 * sigma {
                ui.colored_label(GREEN, format!("✔ Bell violated by {:.0} σ: no local explanation", (s - 2.0) / sigma));
            }
            if let (Some(m0), Some(m1)) = (st.bob_plus_given_a(0), st.bob_plus_given_a(1)) {
                row(ui, "no signalling", format!("P(Bob +) = {m0:.3} | Alice 0°, {m1:.3} | Alice 45°"), (m0 - m1).abs() < 0.03);
            }
        } else {
            ui.label(egui::RichText::new("collecting pairs…").weak());
        }
    } else {
        let n = st.manual.iter().sum::<u64>();
        if let Some(e) = BellStats::correlation(&st.manual) {
            let t = if cfg.hidden_variables { correlation_local(cfg.angle_a, cfg.angle_b) } else { correlation_quantum(cfg.angle_a, cfg.angle_b) };
            let tol = 3.0 / (n as f32).sqrt();
            row(ui, "E(a, b)", format!("{e:+.3}  (theory {t:+.3}, {n} pairs)"), (e - t).abs() <= tol);
        }
        ui.label(egui::RichText::new("Clear the counts after changing the angles.").small().weak());
    }
    ui.label(egui::RichText::new(
        "Alice's splitter is closer: when her photon is measured, Bob's photon —\n\
         still flying the other way — gets its polarization at that instant\n\
         (its grains turn from pale to coloured, its arrows appear).",
    ).small().weak());
}

fn st_last(sim: &Sim) -> Option<(f32, bool, f32, bool)> { sim.stats.last_bell }

fn legend(ui: &mut egui::Ui) {
    ui.collapsing("How to read it", |ui| {
        ui.label(egui::RichText::new(
            "Points = the grains of ONE photon. Colour = the angle of each grain's\n\
             phase arrow: bands of colour are the wave crests, one per λ.\n\
             White arrows = electric field E: they rotate as the grains advance.\n\
             Pale grains = entangled, polarization not defined yet.\n\
             Green dots on the screen = where each photon was detected.\n\n\
             Grains never push each other: light does not interact with light.\n\
             They see each other only where they cross, by comparing phases.",
        ).small());
    });
}

// ─── Small widgets ────────────────────────────────────────────────────────────

fn row(ui: &mut egui::Ui, name: &str, value: String, ok: bool) {
    ui.horizontal(|ui| {
        ui.colored_label(if ok { GREEN } else { ORANGE }, if ok { "✔" } else { "…" });
        ui.label(format!("{name}: {value}"));
    });
}

fn check(ui: &mut egui::Ui, name: &str, measured: f32, theory: f32, rel_tol: f32, what: &str) {
    let ok = (measured / theory - 1.0).abs() <= rel_tol;
    row(ui, name, format!("{measured:.2}  ({what} {theory:.2})"), ok);
}

/// Several profiles over the screen, each scaled to its own maximum.
fn plot(ui: &mut egui::Ui, series: &[(&[f32], Color32, bool)]) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(320.0, 120.0), egui::Sense::hover());
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 3.0, Color32::from_rgb(12, 14, 24));
    p.line_segment([rect.center_top(), rect.center_bottom()], egui::Stroke::new(1.0, Color32::from_gray(40)));
    for (data, color, bars) in series {
        let max = data.iter().cloned().fold(0.0f32, f32::max);
        if max <= 0.0 { continue; }
        let n = data.len() as f32;
        let at = |i: usize, v: f32| egui::pos2(
            rect.left() + (i as f32 + 0.5) / n * rect.width(),
            rect.bottom() - 4.0 - v / max * (rect.height() - 8.0),
        );
        if *bars {
            for (i, &v) in data.iter().enumerate() {
                if v > 0.0 {
                    p.line_segment([at(i, 0.0), at(i, v)], egui::Stroke::new(1.0, color.gamma_multiply(0.7)));
                }
            }
        } else {
            let pts: Vec<egui::Pos2> = data.iter().enumerate().map(|(i, &v)| at(i, v)).collect();
            p.add(egui::Shape::line(pts, egui::Stroke::new(1.5, *color)));
        }
    }
}

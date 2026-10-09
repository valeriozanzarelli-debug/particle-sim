//! Panel of the zig-zag viewer (hidden by default: Tab or ☰). Italian, for a 20-year-old.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::{Color32, RichText};

use crate::atom::{Hydrogen, ALPHA, C_MS, LEVELS, MC2_EV, THETAS};
use crate::camera::OrbitState;
use crate::config::*;
use crate::{camera_for, lambda_c_m, Lab, Modo, STATI};

const CYAN: Color32 = Color32::from_rgb(90, 220, 255);
const GOLD: Color32 = Color32::from_rgb(255, 205, 80);
const GREY: Color32 = Color32::from_rgb(150, 155, 170);
const ORANGE: Color32 = Color32::from_rgb(255, 150, 40);
const BLUE: Color32 = Color32::from_rgb(70, 150, 255);
const GREEN: Color32 = Color32::from_rgb(140, 230, 140);

/// Measured values (cited, not fitted).
const FS_MEAS_HZ: f64 = 10.969_045e9;
const LYA_MEAS_NM: f64 = 121.567;
const HAL_MEAS_NM: f64 = 656.46;
const ION_MEAS_EV: f64 = 13.598_434_6;

fn stato(lab: &Lab) -> String {
    let p = if lab.paused { "  || pausa" } else { "" };
    match lab.modo {
        Modo::Luce => format!("LUCE · palline spin 1 (F = E + icB) · tick {} · centro del pacchetto a {:.3} c{p}", lab.light.ticks, lab.light_speed),
        Modo::Elettrone => format!(
            "ELETTRONE · arancio = destrorso R, blu = sinistrorso L · tick {} · tremolio di <r>: ampiezza {:.3} ħ/mc, {:.2} mc²/ħ{p}",
            lab.dirac.ticks, lab.zbw_radius, lab.zbw_omega),
        Modo::Idrogeno => {
            let ready = lab.hydro.lock().map(|g| g.is_some()).unwrap_or(false);
            if ready {
                format!("IDROGENO · stato {} · colore = fase dell'onda · t = {:.2e} s (1 s di schermo = {:.0e} s){p}",
                        STATI[lab.h_stato].0, lab.h_time, lab.t_vis)
            } else {
                "IDROGENO · le palline stanno ascoltando l'atomo (calcolo dei livelli)...".to_string()
            }
        }
    }
}

pub fn draw_hud(
    mut ctx: EguiContexts,
    mut lab: ResMut<Lab>,
    mut cam: ResMut<OrbitState>,
    keys: Res<ButtonInput<KeyCode>>,
    mut tex: Local<Option<egui::TextureHandle>>,
) {
    let lab = &mut *lab;
    let ctx = ctx.ctx_mut();
    if tex.is_none() {
        let px: Vec<u8> = lab.checker.iter().flatten().copied().collect();
        let img = egui::ColorImage::from_rgba_unmultiplied([CHECK_W, CHECK_T], &px);
        *tex = Some(ctx.load_texture("scacchiera", img, egui::TextureOptions::NEAREST));
    }
    let mut new_mode = None;
    if !ctx.wants_keyboard_input() {
        if keys.just_pressed(KeyCode::Tab) { lab.panel_open = !lab.panel_open; }
        if keys.just_pressed(KeyCode::Space) { lab.paused = !lab.paused; }
        if keys.just_pressed(KeyCode::Digit1) { new_mode = Some(Modo::Luce); }
        if keys.just_pressed(KeyCode::Digit2) { new_mode = Some(Modo::Elettrone); }
        if keys.just_pressed(KeyCode::Digit3) { new_mode = Some(Modo::Idrogeno); }
        if keys.just_pressed(KeyCode::KeyR) {
            match lab.modo { Modo::Luce => lab.reset_light(), Modo::Elettrone => lab.reset_electron(), Modo::Idrogeno => lab.reset_atom_surfers() }
        }
    }
    egui::Area::new(egui::Id::new("menu_button"))
        .anchor(egui::Align2::LEFT_TOP, [10.0, 10.0])
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                let b = egui::Button::new(RichText::new("☰").size(16.0)).fill(Color32::from_black_alpha(90)).selected(lab.panel_open);
                if ui.add(b).on_hover_text("pannello (Tab) · pausa (Spazio) · 1 luce · 2 elettrone · 3 idrogeno · R ricomincia").clicked() {
                    lab.panel_open = !lab.panel_open;
                }
                for (m, name) in [(Modo::Luce, "1 luce"), (Modo::Elettrone, "2 elettrone"), (Modo::Idrogeno, "3 idrogeno")] {
                    if ui.add(egui::Button::new(name).fill(Color32::from_black_alpha(90)).selected(lab.modo == m)).clicked() {
                        new_mode = Some(m);
                    }
                }
            });
        });
    if let Some(m) = new_mode {
        if m != lab.modo {
            lab.modo = m;
            *cam = camera_for(m, lab.h_stato);
            lab.speed = if m == Modo::Idrogeno { 1.0 } else { lab.speed.min(2.0) };
        }
    }
    egui::Area::new(egui::Id::new("status"))
        .anchor(egui::Align2::CENTER_BOTTOM, [0.0, -10.0])
        .show(ctx, |ui| {
            egui::Frame::none().fill(Color32::from_rgba_unmultiplied(10, 12, 20, 170)).rounding(6.0)
                .inner_margin(egui::Margin::symmetric(10.0, 5.0)).show(ui, |ui| {
                    ui.label(RichText::new(stato(lab)).monospace().color(CYAN));
                });
        });
    if !lab.panel_open { return; }
    let mut open = true;
    egui::Window::new("Zig-zag: una sola regola")
        .open(&mut open)
        .default_pos([12.0, 46.0])
        .default_width(500.0)
        .default_height(840.0)
        .resizable(true)
        .vscroll(true)
        .show(ctx, |ui| {
            ui.label(RichText::new("Ogni pallina-onda si muove alla velocità della luce e gira il suo stato interno.").color(GOLD));
            ui.label("Un tick: sposta la parte +1 di una cella avanti e la parte −1 indietro lungo x, y, z \
                      (ordine alternato), poi ruota: la massa mescola sinistrorso e destrorso di un angolo θ.");
            ui.separator();
            ui.horizontal(|ui| {
                ui.label("velocità");
                ui.add(egui::Slider::new(&mut lab.speed, 0.0..=4.0));
                if ui.button(if lab.paused { ">" } else { "||" }).clicked() { lab.paused = !lab.paused; }
                if ui.button("ricomincia").clicked() {
                    match lab.modo { Modo::Luce => lab.reset_light(), Modo::Elettrone => lab.reset_electron(), Modo::Idrogeno => lab.reset_atom_surfers() }
                }
            });
            ui.horizontal(|ui| {
                ui.label("palline");
                ui.add(egui::Slider::new(&mut lab.ball_scale, 0.3..=3.0).text("grandezza"));
                ui.add(egui::Slider::new(&mut lab.gain, 0.2..=4.0).text("luce"));
            });
            ui.separator();
            ingressi(ui, lab.modo);
            ui.separator();
            match lab.modo {
                Modo::Luce => pannello_luce(ui, lab),
                Modo::Elettrone => pannello_elettrone(ui, lab, tex.as_ref()),
                Modo::Idrogeno => pannello_idrogeno(ui, lab, &mut cam),
            }
        });
    if !open { lab.panel_open = false; }
}

fn ingressi(ui: &mut egui::Ui, modo: Modo) {
    ui.label(RichText::new("INGRESSI (messi a mano, non calcolati)").strong().color(GOLD));
    ui.label("• c, ħ: le unità (una cella per tick = c)");
    match modo {
        Modo::Luce => {
            ui.label("• spin 1, nessuna massa (θ = 0), nessuna carica");
        }
        Modo::Elettrone => {
            ui.label("• spin ½; massa dell'elettrone me = 0,511 MeV/c² -> angolo di flip θ = me c a/ħ");
            ui.label(format!("  qui a = {A_ELECTRON} ħ/mc = {:.3e} m, quindi θ = {A_ELECTRON} rad per tick", A_ELECTRON * lambda_c_m()));
        }
        Modo::Idrogeno => {
            ui.label("• spin ½, me (θ = a/λ_C), carica e -> α = 1/137,035999");
            ui.label("• protone: Z = 1, puntiforme; la sua massa entra solo nella massa ridotta μ = me m_p/(me+m_p)");
        }
    }
    ui.label(RichText::new("Il passo del reticolo a NON è fisica: ogni numero è dato per più a ed estrapolato a a -> 0.").color(GREY).small());
}

fn pannello_luce(ui: &mut egui::Ui, lab: &mut Lab) {
    ui.label(RichText::new("LUCE").strong().color(GOLD));
    ui.label("Le palline portano F = E + icB (vettore di Riemann–Silberstein). Si vedono le creste dell'energia elettrica (Re F)².");
    let mut diag = lab.light_diag;
    ui.horizontal(|ui| {
        ui.radio_value(&mut diag, false, "lungo l'asse x");
        ui.radio_value(&mut diag, true, "lungo la diagonale (1,1,1)");
    });
    if diag != lab.light_diag { lab.light_diag = diag; lab.reset_light(); }
    ui.label(format!("velocità del centro del pacchetto (misurata ora): {:.4} c", lab.light_speed));
    ui.label(RichText::new("Laboratorio (lab/zz_light.py, zz_slit*.py):").color(CYAN));
    for s in [
        "• lungo gli assi la regola è esatta: ω = c|k| per ogni λ;",
        "• in diagonale ω è più bassa: −3,5 % a λ = 8a, −0,86 % a 16a, −0,21 % a 32a (va come (a/λ)²);",
        "• trasversalità: energia nel modo fermo 4,3·10-³ a λ = 8a, 3,5·10-^4 a 16a;",
        "• spin 1: fase di Berry su un cono di 60° -> π (Tomita–Chiao 1986), errore 0,005 rad a λ = 64a;",
        "• doppia fenditura: la diffrazione è troppo debole di un fattore ka·cot(ka): a λ = 8a correlazione 0,60 (FALLITO), \
           a λ = 16a 0,98 (schermo a 40 λ).",
    ] { ui.label(s); }
}

fn pannello_elettrone(ui: &mut egui::Ui, lab: &mut Lab, tex: Option<&egui::TextureHandle>) {
    ui.label(RichText::new("ELETTRONE LIBERO").strong().color(GOLD));
    ui.horizontal(|ui| {
        ui.label(RichText::new("■ destrorso R").color(ORANGE));
        ui.label(RichText::new("■ sinistrorso L").color(BLUE));
        ui.label(RichText::new("— surfisti").color(GREEN));
        ui.label(RichText::new("— <r>").color(Color32::WHITE));
    });
    let mut cerchio = lab.cerchio;
    ui.horizontal(|ui| {
        ui.radio_value(&mut cerchio, false, "parte tutto destrorso (zig-zag)");
        ui.radio_value(&mut cerchio, true, "spinore misto (cerchio)");
    });
    if cerchio != lab.cerchio { lab.cerchio = cerchio; lab.reset_electron(); }
    ui.label("Ogni pallina corre a c; il flip della massa la fa girare da R a L e indietro: zig-zag. Un pacchetto fermo e stretto contiene energia positiva e negativa: la loro interferenza fa tremare <r> e i surfisti a 2mc²/ħ (Zitterbewegung, Schrödinger 1930): avanti e indietro lungo lo spin se parte tutto destrorso, in cerchio con lo spinore misto.");
    ui.label(format!("ampiezza del moto di <r> (ora): {:.3} ħ/mc = {:.3e} m   [ħ/2mc = 0,5 per un pacchetto largo]",
                     lab.zbw_radius, lab.zbw_radius * lambda_c_m()));
    ui.label(format!("frequenza (ora): {:.3} mc²/ħ   [2mc²/ħ = 1,553·10²¹ rad/s]", lab.zbw_omega));
    // dispersion check from the Bloch matrix of the walk (computed here)
    let th = A_ELECTRON;
    let mut s = String::from("dispersione della regola (θ = 0,25): ");
    for p in [0.0, 0.5, 1.0] {
        let e = crate::walk::walk_energy(p, th);
        s += &format!("p = {p}: E = {e:.5} (Dirac {:.5}); ", (1.0f64 + p * p).sqrt());
    }
    ui.label(s);
    ui.label(RichText::new("Laboratorio (lab/zz_free.py):").color(CYAN));
    for s in [
        "• E² = p²c² + m²c^4: energia a riposo esatta; errore 0,19 % fino a p = mc con θ = 0,1, scende come θ²;",
        "• massa 'cinetica' (1 + θ²/3) m: è l'errore del reticolo, sparisce per a -> 0;",
        "• Zitterbewegung 1-D: ω = 2,0095 mc²/ħ, ampiezza 0,491 ħ/mc (continuo 2,0095 e 0,4905);",
        "• spin ½: una rotazione di 90° del reticolo vuole lo spinore girato di 45° (metà angolo): 4 giri da 90° = −1;",
        "• g dalla corrente di Dirac: vedi docs/ZIGZAG.md.",
    ] { ui.label(s); }
    if let Some(t) = tex {
        ui.separator();
        ui.label(RichText::new("Scacchiera di Feynman (1-D): una sola pallina R al centro, il tempo va in su").color(CYAN));
        ui.image((t.id(), egui::vec2(CHECK_W as f32 * 1.4, CHECK_T as f32 * 1.4)));
        ui.label(RichText::new("arancio = R (va a destra a c), blu = L (va a sinistra a c); ogni flip è un angolo della scacchiera.").small());
    }
}

fn pannello_idrogeno(ui: &mut egui::Ui, lab: &mut Lab, cam: &mut OrbitState) {
    ui.label(RichText::new("IDROGENO").strong().color(GOLD));
    let guard = lab.hydro.lock().unwrap();
    let Some(h) = guard.as_ref() else {
        ui.label("calcolo dei livelli in corso (tre reticoli, nove livelli)...");
        return;
    };
    let h: &Hydrogen = h;
    let mut st = lab.h_stato;
    ui.horizontal(|ui| {
        for (i, s) in STATI.iter().enumerate() { ui.radio_value(&mut st, i, s.0); }
    });
    ui.horizontal(|ui| {
        ui.label("tempo: 1 s di schermo =");
        let mut e = lab.t_vis.log10();
        if ui.add(egui::Slider::new(&mut e, -18.0..=-15.0).text("10^x s")).changed() { lab.t_vis = 10f64.powf(e); }
    });
    ui.label(format!("livelli calcolati dalla regola in {:.1} s (θ = 0,2 · 0,1 · 0,05 -> a -> 0); nodi giusti: {}",
                     h.seconds, if h.nodes_ok { "sì" } else { "NO" }));
    egui::Grid::new("livelli").striped(true).show(ui, |ui| {
        ui.label("livello"); for t in THETAS { ui.label(format!("θ={t}")); } ui.label("a -> 0"); ui.label("Dirac (formula)"); ui.end_row();
        for (k, l) in LEVELS.iter().enumerate() {
            ui.label(l.0);
            for v in h.e_lat[k] { ui.label(format!("{v:.5}")); }
            ui.label(RichText::new(format!("{:.6}", h.e_ext[k])).color(GOLD));
            ui.label(format!("{:.6}", h.e_dirac[k]));
            ui.end_row();
        }
    });
    ui.label(RichText::new("energie di legame in eV, nucleo infinitamente pesante").small().color(GREY));
    let mu = Hydrogen::mu();
    let fs = h.fine_structure_hz();
    let lya = h.lyman_alpha_nm();
    let hal = h.balmer_alpha_nm();
    let ion = h.e_ext[0] * mu;
    let v1 = h.surfer_1s();
    egui::Grid::new("misure").striped(true).show(ui, |ui| {
        ui.label("grandezza"); ui.label("regola"); ui.label("misura"); ui.label("differenza"); ui.end_row();
        ui.label("ionizzazione 1s (massa ridotta)"); ui.label(format!("{ion:.5} eV")); ui.label(format!("{ION_MEAS_EV:.5} eV"));
        ui.label(format!("{:+.1e} (reticolo sugli s; Lamb 1s +2e-6)", ion / ION_MEAS_EV - 1.0)); ui.end_row();
        ui.label("struttura fine 2p3/2 − 2p1/2"); ui.label(format!("{:.4} GHz", fs / 1e9)); ui.label("10,9690 GHz");
        ui.label(format!("{:+.1} MHz (QED)", (fs - FS_MEAS_HZ) / 1e6)); ui.end_row();
        ui.label("Lyman-α"); ui.label(format!("{lya:.4} nm")); ui.label(format!("{LYA_MEAS_NM} nm")); ui.label(format!("{:+.0} ppm (reticolo sull'1s)", (lya / LYA_MEAS_NM - 1.0) * 1e6)); ui.end_row();
        ui.label("Balmer-α"); ui.label(format!("{hal:.4} nm")); ui.label(format!("{HAL_MEAS_NM} nm")); ui.label(format!("{:+.0} ppm", (hal / HAL_MEAS_NM - 1.0) * 1e6)); ui.end_row();
        ui.label("surfista 1s, |v| all'equatore"); ui.label(format!("{:.4} αc", v1)); ui.label("αc = 2,188·10^6 m/s (Dirac)");
        ui.label(format!("{:.3e} m/s", v1 * ALPHA * C_MS)); ui.end_row();
    });
    ui.label(RichText::new(format!("riposo: mec² = {MC2_EV:.2} eV; i livelli sono le frequenze a cui la regola 'suona' (vedi lab/zz_listen.py)")).small().color(GREY));
    ui.label(RichText::new("Che cosa NON spiega").color(CYAN));
    for s in [
        "• perché la massa dell'elettrone e α hanno quei valori (sono ingressi);",
        "• lo spostamento di Lamb (1057,8 MHz) e g − 2: QED, il vuoto quantistico;",
        "• il protone (è solo una carica puntiforme), il rinculo oltre la massa ridotta.",
    ] { ui.label(s); }
    if st != lab.h_stato {
        drop(guard);
        lab.h_stato = st;
        lab.reset_atom_surfers();
        *cam = camera_for(Modo::Idrogeno, st);
        return;
    }
}

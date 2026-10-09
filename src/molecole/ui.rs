//! Panel of the molecule viewer (hidden by default: Tab or ☰): pick the system, computed vs measured
//! numbers, the light of hydrogen, and what the cloud is (and is not).

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::Color32;

use crate::costanti as k;
use crate::scena::Scelta;
use crate::Lab;

const CYAN: Color32 = Color32::from_rgb(90, 220, 255);
const GOLD: Color32 = Color32::from_rgb(255, 205, 80);
const GREY: Color32 = Color32::from_rgb(150, 155, 170);
const GREEN: Color32 = Color32::from_rgb(170, 255, 140);

pub fn draw_hud(mut ctx: EguiContexts, mut lab: ResMut<Lab>, keys: Res<ButtonInput<KeyCode>>) {
    let lab = &mut *lab;
    let ctx = ctx.ctx_mut();
    if keys.just_pressed(KeyCode::Tab) && !ctx.wants_keyboard_input() {
        lab.panel_open = !lab.panel_open;
    }
    egui::Area::new(egui::Id::new("menu_button"))
        .anchor(egui::Align2::LEFT_TOP, [10.0, 10.0])
        .show(ctx, |ui| {
            let b = egui::Button::new(egui::RichText::new("☰").size(16.0))
                .fill(Color32::from_black_alpha(90))
                .selected(lab.panel_open);
            if ui.add(b).on_hover_text("pannello (Tab)").clicked() {
                lab.panel_open = !lab.panel_open;
            }
        });
    // one discreet line at the bottom: what is on screen
    egui::Area::new(egui::Id::new("status"))
        .anchor(egui::Align2::CENTER_BOTTOM, [0.0, -10.0])
        .show(ctx, |ui| {
            egui::Frame::none().fill(Color32::from_rgba_unmultiplied(10, 12, 20, 170)).rounding(6.0)
                .inner_margin(egui::Margin::symmetric(10.0, 5.0)).show(ui, |ui| {
                    ui.label(egui::RichText::new(format!(
                        "{}   ·   {} palline   ·   E (dal vivo) = {:.4} hartree = {:.3} eV",
                        lab.scelta.nome(), lab.dmc.n(), lab.e_media, lab.e_media * k::hartree_ev()))
                        .monospace().color(CYAN));
                });
        });
    if lab.panel_open {
        let mut open = true;
        egui::Window::new("Molecole dalle palline")
            .open(&mut open)
            .default_width(400.0)
            .default_pos([12.0, 48.0])
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().max_height(780.0).show(ui, |ui| panel(ui, lab));
            });
        if !open { lab.panel_open = false; }
    }
}

fn panel(ui: &mut egui::Ui, lab: &mut Lab) {
    let mut pick = lab.scelta;
    ui.horizontal_wrapped(|ui| {
        for s in Scelta::TUTTE { ui.selectable_value(&mut pick, s, s.nome()); }
    });
    if pick != lab.scelta { lab.cambia(pick); }
    ui.separator();
    ui.horizontal(|ui| {
        if lab.scelta.elettroni() == 2 {
            ui.checkbox(&mut lab.mostra_e[0], egui::RichText::new("elettrone spin su").color(CYAN));
            ui.checkbox(&mut lab.mostra_e[1], egui::RichText::new("elettrone spin giù").color(Color32::from_rgb(255, 160, 70)));
        }
        ui.checkbox(&mut lab.mostra_surfer, "surfer");
        if ui.button("nuovo surfer").clicked() { lab.nuovo_surfer(); }
    });
    if lab.scelta.biatomica() {
        ui.checkbox(&mut lab.vibra, "vibrazione dei nuclei (frequenza calcolata, rallentata)");
    }
    ui.label(egui::RichText::new(format!("tick τ = {} a.u. di tempo ({:.2e} s) · E media dal vivo {:.5} hartree",
        lab.dmc.tau, lab.dmc.tau * k::au_time_s(), lab.e_media)).small().color(GREY));
    ui.separator();
    tabella(ui, lab);
    ui.separator();
    spiegazione(ui, lab);
}

/// computed ± error | measured | difference
fn riga(ui: &mut egui::Ui, nome: &str, v: Option<f64>, e: Option<f64>, exp: Option<f64>, fmt: usize) {
    ui.label(nome);
    match v {
        Some(v) => {
            let s = match e { Some(e) if e.is_finite() => format!("{v:.fmt$} ± {e:.fmt$}"), _ => format!("{v:.fmt$}") };
            ui.label(egui::RichText::new(s).color(CYAN));
        }
        None => { ui.label(egui::RichText::new("— (non ancora misurato)").color(GREY)); }
    }
    match exp { Some(x) => { ui.label(format!("{x:.fmt$}")); } None => { ui.label("—"); } }
    match (v, exp) {
        (Some(v), Some(x)) => { ui.label(egui::RichText::new(format!("{:+.2} %", 100.0 * (v - x) / x.abs())).color(GOLD)); }
        _ => { ui.label(""); }
    }
    ui.end_row();
}

fn tabella(ui: &mut egui::Ui, lab: &Lab) {
    let r = &lab.res;
    let ev = k::hartree_ev();
    let f = |key: &str| r.f(key);
    let fe = |key: &str, s: f64| r.f(key).map(|x| x * s);
    ui.label(egui::RichText::new("Calcolato con le palline  vs  misurato").strong());
    egui::Grid::new("tab").striped(true).num_columns(4).show(ui, |ui| {
        ui.label(egui::RichText::new("grandezza").weak());
        ui.label(egui::RichText::new("calcolato").weak());
        ui.label(egui::RichText::new("riferimento").weak());
        ui.label(egui::RichText::new("diff.").weak());
        ui.end_row();
        match lab.scelta {
            Scelta::H1s | Scelta::H2p | Scelta::H2s | Scelta::H3d | Scelta::H3p => {
                riga(ui, "E 1s (eV), senza guida", fe("h_E_plain", ev), fe("h_E_plain_err", ev), Some(-13.605693), 4);
                riga(ui, "E 1s (eV), con guida", fe("h_E_guided", ev), fe("h_E_guided_err", ev), Some(-13.605693), 5);
                for (st, n) in [("2p", 2.0), ("2s", 2.0), ("3d", 3.0), ("3p", 3.0)] {
                    riga(ui, &format!("E {st} (eV)"), fe(&format!("h{st}_E"), ev), fe(&format!("h{st}_E_err"), ev),
                         Some(-13.605693 / (n * n)), 5);
                }
                riga(ui, "nodo 2s (a0)", f("h2s_node"), f("h2s_node_err"), Some(2.0), 3);
                riga(ui, "nodo 3p (a0)", f("h3p_node"), f("h3p_node_err"), Some(6.0), 3);
            }
            Scelta::H2piu => {
                riga(ui, "distanza Re (Å)", f("h2p_mol_re"), f("h2p_mol_re_err"), Some(1.052), 4);
                riga(ui, "ωe (cm-1)", f("h2p_mol_we"), f("h2p_mol_we_err"), Some(2321.7), 1);
                riga(ui, "ν(0->1) (cm-1)", f("h2p_mol_nu01"), f("h2p_mol_nu01_err"), Some(2191.13), 1);
                riga(ui, "D0 (eV)", f("h2p_mol_d0"), f("h2p_mol_d0_err"), Some(2.6507), 4);
            }
            Scelta::H2 => {
                riga(ui, "distanza Re (Å)", f("h2_mol_re"), f("h2_mol_re_err"), Some(0.74144), 4);
                riga(ui, "ωe (cm-1)", f("h2_mol_we"), f("h2_mol_we_err"), Some(4401.21), 1);
                riga(ui, "ν(0->1) (cm-1)", f("h2_mol_nu01"), f("h2_mol_nu01_err"), Some(4161.17), 1);
                riga(ui, "D0 (eV)", f("h2_mol_d0"), f("h2_mol_d0_err"), Some(4.47813), 4);
            }
            Scelta::He => {
                riga(ui, "E (eV), nucleo fermo", fe("he_E", ev), fe("he_E_err", ev), Some(-79.0052), 4);
                riga(ui, "E (eV), massa del nucleo", f("he_E_ev_mu"), fe("he_E_err", ev), Some(-79.0052), 4);
                riga(ui, "E 1s2p ¹P (Eₕ), nodo fisso", f("he1p_E"), f("he1p_E_err"), Some(-2.123843), 5);
                riga(ui, "riga 1¹S->2¹P λ (nm)", f("he1p_lambda_nm"), f("he1p_lambda_err"), Some(58.4334), 3);
                riga(ui, "forza d'oscillatore f", f("he1p_f"), f("he1p_f_err"), Some(0.2762), 4);
                riga(ui, "A(2¹P->1¹S) (1e9/s)", f("he1p_A").map(|x| x / 1e9), f("he1p_A_err").map(|x| x / 1e9), Some(1.799), 3);
                riga(ui, "polarizzabilità α (a0³)", f("he_alpha"), f("he_alpha_err"), Some(1.3832), 3);
            }
        }
    });
    ui.label(egui::RichText::new(match lab.scelta {
        Scelta::H2piu | Scelta::H2 => "Riferimenti: Huber & Herzberg 1979; D0(H2): Liu et al. 2009; ν(0->1): Dickenson et al. 2013.",
        Scelta::He => "Riferimenti: NIST ASD (energie di ionizzazione 24.5874 + 54.4178 eV; riga 58.4334 nm, f, A);                        E(2¹P) esatta non relativistica (Drake); α esatta 1.3832. Lo scontro con un fotone: cargo run --bin scontro",
        _ => "Livelli: soluzione esatta −13.6057 eV / n² (nucleo fermo). La misura vera è nelle righe di luce qui sotto (NIST).",
    }).small().color(GREY));
    if matches!(lab.scelta, Scelta::H2piu | Scelta::H2) {
        let key = if lab.scelta == Scelta::H2 { "h2_mol_line_um" } else { "h2p_mol_line_um" };
        if let Some(l) = r.f(key) {
            ui.label(egui::RichText::new(format!("La sua riga di vibrazione (0->1): λ = {l:.3} µm, infrarosso (misura: {:.3} µm).",
                if lab.scelta == Scelta::H2 { 1e4 / 4161.17 } else { 1e4 / 2191.13 })).color(GREEN));
        }
    }
    if !matches!(lab.scelta, Scelta::H2piu | Scelta::H2 | Scelta::He) { righe_luce(ui, lab); }
}

/// Wavelength (nm) → approximate sRGB of the colour seen; None outside 380–780 nm.
pub fn colore_lambda(nm: f64) -> Option<Color32> {
    if !(380.0..=780.0).contains(&nm) { return None; }
    let (r, g, b) = match nm {
        l if l < 440.0 => (-(l - 440.0) / 60.0, 0.0, 1.0),
        l if l < 490.0 => (0.0, (l - 440.0) / 50.0, 1.0),
        l if l < 510.0 => (0.0, 1.0, -(l - 510.0) / 20.0),
        l if l < 580.0 => ((l - 510.0) / 70.0, 1.0, 0.0),
        l if l < 645.0 => (1.0, -(l - 645.0) / 65.0, 0.0),
        _ => (1.0, 0.0, 0.0),
    };
    let fade = if nm > 700.0 { 0.3 + 0.7 * (780.0 - nm) / 80.0 } else if nm < 420.0 { 0.3 + 0.7 * (nm - 380.0) / 40.0 } else { 1.0 };
    let c = |x: f64| ((x * fade).clamp(0.0, 1.0).powf(0.8) * 255.0) as u8;
    Some(Color32::from_rgb(c(r), c(g), c(b)))
}

fn righe_luce(ui: &mut egui::Ui, lab: &Lab) {
    let r = &lab.res;
    ui.separator();
    ui.label(egui::RichText::new("La luce che l'idrogeno emette (ΔE = hc/λ, nel vuoto)").strong());
    egui::Grid::new("righe").striped(true).num_columns(4).show(ui, |ui| {
        for (nome, key) in [("Lyman-α  2p->1s", "lyman_alpha"), ("Lyman-β  3p->1s", "lyman_beta"), ("Balmer-α  3d->2p", "balmer_alpha")] {
            let lam = r.f(&format!("line_{key}_nm"));
            let err = r.f(&format!("line_{key}_err"));
            let nist = r.f(&format!("line_{key}_nist"));
            let (rect, _) = ui.allocate_exact_size(egui::vec2(18.0, 14.0), egui::Sense::hover());
            match lam.and_then(colore_lambda) {
                Some(c) => { ui.painter().rect_filled(rect, 2.0, c); }
                None => { ui.painter().rect_stroke(rect, 2.0, egui::Stroke::new(1.0, Color32::from_rgb(140, 90, 220))); }
            }
            ui.label(nome);
            match (lam, err) {
                (Some(l), Some(e)) => { ui.label(egui::RichText::new(format!("{l:.3} ± {e:.3} nm")).color(CYAN)); }
                _ => { ui.label("—"); }
            }
            let visib = lam.map_or("", |l| if l < 380.0 { "UV, invisibile" } else { "visibile" });
            ui.label(format!("NIST {:.3} nm · {visib}", nist.unwrap_or(f64::NAN)));
            ui.end_row();
        }
    });
    ui.separator();
    ui.label(egui::RichText::new("Regole per la luce dell'atomo di idrogeno").strong());
    egui::Grid::new("regole").striped(true).num_columns(4).show(ui, |ui| {
        riga(ui, "polarizzabilità α (a0³)", r.f("h_alpha"), r.f("h_alpha_err"), Some(4.5), 3);
        riga(ui, "<1s|z|2p> (a0)", r.f("h_z12"), r.f("h_z12_err"), Some(0.7449), 4);
        riga(ui, "f(1s->2p)", r.f("h_f12"), None, Some(0.4162), 4);
        riga(ui, "A(2p->1s) (1e8/s)", r.f("h_A21").map(|x| x / 1e8), r.f("h_A21_err").map(|x| x / 1e8), Some(6.2649), 3);
        riga(ui, "σ Rayleigh 500 nm (1e-32 m²)", r.f("h_rayleigh500").map(|x| x / 1e-32), None, Some(9.289), 3);
        riga(ui, "σ risonante Lyα (1e-15 m²)", r.f("h_sigma_res").map(|x| x / 1e-15), None, Some(7.056), 3);
        riga(ui, "max fotoni/s diffusi (1e8)", r.f("h_rate_max").map(|x| x / 1e8), None, Some(3.132), 3);
    });
}

fn spiegazione(ui: &mut egui::Ui, lab: &Lab) {
    let ne = lab.scelta.elettroni();
    ui.label(egui::RichText::new("Che cosa stai guardando").strong());
    ui.label(egui::RichText::new(
        "Ogni punto è una PALLINA-ONDA dell'elettrone. A ogni tick fa solo due cose: un passo a caso \
         (il 'tremolio' quantistico, D = ħ/2mₑ) e poi si duplica o sparisce secondo l'energia nel punto \
         dove si trova (attrazione dei nuclei, repulsione degli altri elettroni). Qui le palline sono anche \
         guidate dall'orbitale che il livello sotto (l'atomo di idrogeno) ha già misurato: la guida cambia \
         solo il rumore, non la risposta.").small());
    ui.label(egui::RichText::new(
        "La nuvola è la DENSITÀ elettronica |ψ|² (quello che misura la diffrazione dei raggi X), non tanti \
         elettroni: l'elettrone è uno solo per colore. Senza guida la nuvola sarebbe ψ stessa.").small().color(GOLD));
    if ne == 2 {
        ui.label(egui::RichText::new(
            "Due elettroni: una pallina è UNA coppia di punti (azzurro ↑ e arancio ↓) che si muovono insieme \
             nello spazio a 6 dimensioni. Per questo non c'è un surfer per elettrone: il surfer è una \
             configurazione intera (le due stelle bianche unite da un filo).").small());
    }
    ui.label(egui::RichText::new(
        "Surfer (pallina di Bohm): per uno stato stazionario reale la fase è costante, quindi la sua \
         velocità ∇S è zero: resta FERMO nel punto dove è nato, scelto a caso dalla nuvola.").small());
    if matches!(lab.scelta, Scelta::H2p | Scelta::H3d | Scelta::H2s | Scelta::H3p) {
        ui.label(egui::RichText::new(
            "Stato eccitato: le superfici viola sono i NODI. Una pallina che ci passa sopra muore. I piani \
             vengono dalla simmetria; le sfere (2s, 3p) le hanno trovate le palline: il nodo è dove le due \
             tasche, dentro e fuori, hanno la stessa energia.").small().color(Color32::from_rgb(190, 160, 255)));
    }
    if lab.vibra {
        ui.label(egui::RichText::new(
            "Vibrazione: i nuclei oscillano con la frequenza calcolata dalla curvatura della curva E(R) \
             (ampiezza del punto zero), rallentata di ~1e14 volte; le palline seguono i nuclei.").small());
    }
}

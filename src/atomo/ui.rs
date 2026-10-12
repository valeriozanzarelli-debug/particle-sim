//! Interface: a clean screen (the atom, a time bar, a scale bar, ☰), everything else in one panel hidden at
//! start (Tab or ☰). Italian, simple words.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::{Color32, RichText};

use crate::camera::{Bersaglio, Cam};
use crate::carbonio::Stato;
use crate::{costanti as k, Lab, Modo};

const CYAN: Color32 = Color32::from_rgb(90, 220, 255);
const GOLD: Color32 = Color32::from_rgb(255, 205, 80);
const GREY: Color32 = Color32::from_rgb(150, 155, 170);
const GREEN: Color32 = Color32::from_rgb(140, 230, 140);
const RED: Color32 = Color32::from_rgb(255, 120, 100);
const BG: Color32 = Color32::from_rgba_premultiplied(8, 10, 18, 190);

/// Reference values (cited, never used by the simulation).
pub const E_ESATTA_C: f64 = -37.8450; // Chakravorty et al., Phys. Rev. A 47, 3649 (1993), non-relativistic
pub const E_ESATTA_CPIU: f64 = -37.4313; // same source
pub const IP_NIST_EV: f64 = 11.2603; // NIST ASD
pub const ECC_NIST_EV: f64 = 7.482; // NIST ASD, 2s²2p3s ³P° − 2s²2p² ³P, J-averaged
pub const PROTONE_RMS_FM: f64 = 0.84; // charge radius of the proton (CODATA 2018: 0.8414 fm)

pub fn fmt_t(s: f64) -> String {
    let a = s.abs();
    let (v, u) = if a < 1e-15 { (s * 1e18, "as") } else if a < 1e-12 { (s * 1e15, "fs") } else if a < 1e-9 { (s * 1e12, "ps") }
                 else if a < 1e-6 { (s * 1e9, "ns") } else { (s * 1e6, "µs") };
    format!("{:.3} {}", v, u)
}

pub fn fmt_len(m: f64) -> String {
    let a = m.abs();
    if a >= 1e-10 { format!("{} Å", nice(m * 1e10)) } else if a >= 1e-12 { format!("{} pm", nice(m * 1e12)) } else { format!("{} fm", nice(m * 1e15)) }
}

fn nice(v: f64) -> String { if v >= 10.0 { format!("{:.0}", v) } else if v >= 1.0 { format!("{:.1}", v) } else { format!("{:.2}", v) } }

/// A round length (1, 2, 5 × 10ⁿ metres) close to `m`.
fn tondo(m: f64) -> f64 {
    let e = m.log10().floor();
    let b = 10f64.powf(e);
    let f = m / b;
    let g = if f < 1.5 { 1.0 } else if f < 3.5 { 2.0 } else if f < 7.5 { 5.0 } else { 10.0 };
    g * b
}

fn livello(r_bohr: f64) -> &'static str {
    let m = r_bohr * k::bohr_m();
    if m > 2e-11 { "atomo" } else if m > 2e-13 { "tra gli elettroni" } else if m > 3e-15 { "vicino al nucleo" } else { "dentro il nucleo" }
}

pub fn draw_hud(mut ctx: EguiContexts, mut lab: ResMut<Lab>, mut cam: ResMut<Cam>, keys: Res<ButtonInput<KeyCode>>, auto: Res<crate::Auto>) {
    let lab = &mut *lab;
    let ctx = ctx.ctx_mut();
    if !ctx.wants_keyboard_input() && !auto.attivo {
        if keys.just_pressed(KeyCode::Tab) { lab.panel_open = !lab.panel_open; }
        if keys.just_pressed(KeyCode::Space) { lab.modo = if lab.modo == Modo::Avanti { Modo::Pausa } else { Modo::Avanti }; }
        if keys.just_pressed(KeyCode::KeyF) { lab.spara = true; }
        if keys.just_pressed(KeyCode::Digit1) { cam.vai(Bersaglio::Atomo, 9.0); }
        if keys.just_pressed(KeyCode::Digit2) { vai_nucleo(lab, &mut cam); }
        if keys.just_pressed(KeyCode::Digit3) { vai_elettrone(lab, &mut cam); }
    }
    // ── top left: menu and the three zoom levels ──
    egui::Area::new(egui::Id::new("menu")).anchor(egui::Align2::LEFT_TOP, [10.0, 10.0]).show(ctx, |ui| {
        ui.horizontal(|ui| {
            let b = egui::Button::new(RichText::new("☰").size(16.0)).fill(Color32::from_black_alpha(90)).selected(lab.panel_open);
            if ui.add(b).on_hover_text("pannello (Tab) · pausa (Spazio) · spara un fotone (F) · 1 atomo · 2 nucleo · 3 elettrone").clicked() {
                lab.panel_open = !lab.panel_open;
            }
            let bt = cam.bersaglio;
            if ui.add(egui::Button::new("1 atomo").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Atomo)).clicked() { cam.vai(Bersaglio::Atomo, 9.0); }
            if ui.add(egui::Button::new("2 nucleo").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Nucleo)).clicked() { vai_nucleo(lab, &mut cam); }
            if ui.add(egui::Button::new("3 elettrone").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Elettrone)).clicked() { vai_elettrone(lab, &mut cam); }
        });
    });
    // ── top centre: the state and the last event ──
    let p = lab.storia.frame(lab.vista).clone();
    egui::Area::new(egui::Id::new("stato")).anchor(egui::Align2::CENTER_TOP, [0.0, 10.0]).show(ctx, |ui| {
        egui::Frame::none().fill(BG).rounding(6.0).inner_margin(egui::Margin::symmetric(10.0, 4.0)).show(ui, |ui| {
            let col = match p.stato { Stato::Fondamentale => CYAN, Stato::Eccitato => GOLD, Stato::Ione => RED };
            ui.label(RichText::new(p.stato.nome()).monospace().color(col));
            if let Some(e) = p.registro.last() { ui.label(RichText::new(e).small().color(GREY)); }
        });
    });
    // ── bottom left: the scale bar ──
    let bpp = cam.bohr_per_px();
    let target_m = 140.0 * bpp * k::bohr_m();
    let len_m = tondo(target_m);
    let px = (len_m / k::bohr_m() / bpp) as f32;
    egui::Area::new(egui::Id::new("scala")).anchor(egui::Align2::LEFT_BOTTOM, [16.0, -70.0]).show(ctx, |ui| {
        let (rect, _) = ui.allocate_exact_size(egui::vec2(px.max(10.0) + 4.0, 26.0), egui::Sense::hover());
        let pnt = ui.painter();
        let y = rect.bottom() - 4.0;
        let x0 = rect.left() + 2.0;
        let stroke = egui::Stroke::new(2.0, Color32::from_rgb(220, 225, 235));
        pnt.line_segment([egui::pos2(x0, y), egui::pos2(x0 + px, y)], stroke);
        pnt.line_segment([egui::pos2(x0, y - 5.0), egui::pos2(x0, y + 1.0)], stroke);
        pnt.line_segment([egui::pos2(x0 + px, y - 5.0), egui::pos2(x0 + px, y + 1.0)], stroke);
        pnt.text(egui::pos2(x0, y - 8.0), egui::Align2::LEFT_BOTTOM, format!("{}   ·  {}", fmt_len(len_m), livello(cam.r())),
                 egui::FontId::monospace(13.0), Color32::from_rgb(220, 225, 235));
    });
    barra_tempo(ctx, lab, &p);
    if lab.panel_open { pannello(ctx, lab, &mut cam); }
}

fn vai_nucleo(lab: &mut Lab, cam: &mut Cam) {
    lab.mostra_protone = true;
    cam.vai(Bersaglio::Nucleo, 3.0e-15 / k::bohr_m());
}

fn vai_elettrone(lab: &mut Lab, cam: &mut Cam) {
    if !lab.mostra_elettrone || lab.dirac.is_none() {
        lab.mostra_elettrone = true;
        lab.scegli_elettrone();
        lab.reset_dirac();
    }
    cam.vai(Bersaglio::Elettrone, 30.0 * k::alpha_fs());
}

fn barra_tempo(ctx: &egui::Context, lab: &mut Lab, p: &crate::mondo::Piccolo) {
    egui::Area::new(egui::Id::new("tempo")).anchor(egui::Align2::CENTER_BOTTOM, [0.0, -10.0]).show(ctx, |ui| {
        egui::Frame::none().fill(BG).rounding(8.0).inner_margin(egui::Margin::symmetric(12.0, 6.0)).show(ui, |ui| {
            ui.horizontal(|ui| {
                let back = ui.add(egui::Button::new("⏮").selected(lab.modo == Modo::Indietro)).on_hover_text("indietro nel tempo");
                if back.clicked() { lab.modo = if lab.modo == Modo::Indietro { Modo::Pausa } else { Modo::Indietro }; }
                let play = if lab.modo == Modo::Avanti { "⏸" } else { "▶" };
                if ui.button(play).on_hover_text("avanti / pausa (Spazio)").clicked() {
                    lab.modo = if lab.modo == Modo::Avanti { Modo::Pausa } else { Modo::Avanti };
                }
                if ui.button("⏭").on_hover_text("velocità ×10").clicked() { lab.velocita = (lab.velocita * 10.0).min(1e-9); }
                let mut lv = lab.velocita.log10();
                ui.add(egui::Slider::new(&mut lv, -20.0..=-9.0).show_value(false)).on_hover_text("velocità del tempo");
                lab.velocita = 10f64.powf(lv);
                ui.label(RichText::new(format!("1 s = {}", fmt_t(lab.velocita))).monospace().color(GREY));
                ui.separator();
                let (a, b) = (lab.storia.inizio, lab.storia.ultimo());
                let mut v = lab.vista;
                let r = ui.add(egui::Slider::new(&mut v, a..=b.max(a + 1)).show_value(false)).on_hover_text("la storia conservata");
                if r.changed() { lab.modo = Modo::Pausa; lab.vista = v.min(b); }
                ui.label(RichText::new(format!("t = {}", fmt_t(p.t_s))).monospace().color(CYAN));
                ui.separator();
                ui.label("T");
                ui.add(egui::Slider::new(&mut lab.t_kelvin, 0.0..=3000.0).suffix(" K").max_decimals(0));
            });
        });
    });
}

fn pannello(ctx: &egui::Context, lab: &mut Lab, cam: &mut Cam) {
    let mut open = true;
    egui::Window::new("Carbonio dalle palline").open(&mut open).default_pos([12.0, 46.0]).default_width(470.0)
        .default_height(820.0).resizable(true).vscroll(true).show(ctx, |ui| {
        let lv = lab.lv.clone();
        ui.label(RichText::new("Un atomo di carbonio fatto solo dalle nostre regole.").color(GOLD));
        ui.label("Ingressi: ħ, e, mₑ, ε₀, c, la carica Z = 6 e la massa del nucleo. Nessun orbitale da tabella.");
        ui.separator();
        // ── photons ──
        ui.label(RichText::new("FOTONE").strong().color(GOLD));
        ui.horizontal(|ui| {
            ui.label("λ");
            ui.add(egui::Slider::new(&mut lab.lambda_nm, 40.0..=2000.0).logarithmic(true).suffix(" nm").max_decimals(2));
        });
        ui.horizontal(|ui| {
            if ui.button(format!("2p→3s ({:.1} nm, calcolata)", lv.lambda_ecc_nm())).clicked() { lab.lambda_nm = lv.lambda_ecc_nm(); }
            if ui.button(format!("ionizza (< {:.1} nm)", lv.lambda_ip_nm())).clicked() { lab.lambda_nm = (lv.lambda_ip_nm() * 0.8 * 100.0).round() / 100.0; }
            if ui.button("visibile 500 nm").clicked() { lab.lambda_nm = 500.0; }
        });
        ui.horizontal(|ui| {
            ui.label("da");
            egui::ComboBox::from_id_salt("dir").selected_text(["sinistra", "davanti", "destra", "dietro", "sotto", "sopra"][lab.dir_sparo % 6]).show_ui(ui, |ui| {
                for (i, n) in ["sinistra", "davanti", "destra", "dietro", "sotto", "sopra"].iter().enumerate() { ui.selectable_value(&mut lab.dir_sparo, i, *n); }
            });
            if ui.add(egui::Button::new(RichText::new("  Spara (F)  ").strong()).fill(Color32::from_rgb(40, 60, 110))).clicked() { lab.spara = true; }
        });
        let e_ev = crate::mondo::hartree_da_nm(lab.lambda_nm) * k::hartree_ev();
        ui.label(RichText::new(format!("energia del fotone {:.3} eV · il fotone è disegnato ~1000 volte più corto del vero", e_ev)).small().color(GREY));
        ui.label(RichText::new("Assorbito solo dove le palline hanno calcolato un livello: 2p→3s (entro l'errore del calcolo) \
                                o sopra la soglia di ionizzazione. Altrimenti passa.").small().color(GREY));
        ui.separator();
        // ── what is shown ──
        ui.label(RichText::new("VISTA").strong().color(GOLD));
        ui.horizontal(|ui| {
            ui.checkbox(&mut lab.mostra_nuvola, "nuvola");
            ui.checkbox(&mut lab.mostra_protone, "mostra protone vero");
            let mut e = lab.mostra_elettrone;
            if ui.checkbox(&mut e, "mostra elettrone vero").changed() {
                lab.mostra_elettrone = e;
                if e { lab.scegli_elettrone(); lab.reset_dirac(); }
            }
        });
        ui.horizontal(|ui| {
            if ui.button("vai al nucleo (2)").clicked() { vai_nucleo(lab, cam); }
            if ui.button("vai all'elettrone (3)").clicked() { vai_elettrone(lab, cam); }
            if ui.button("atomo (1)").clicked() { cam.vai(Bersaglio::Atomo, 9.0); }
        });
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut lab.grandezza, 0.3..=3.0).text("palline"));
            ui.add(egui::Slider::new(&mut lab.luce, 0.2..=4.0).text("luce"));
            ui.checkbox(&mut lab.bagliore, "bagliore");
        });
        let conf: Vec<String> = lab.confini.iter().map(|c| format!("{:.2} Å", c * k::bohr_ang())).collect();
        ui.label(RichText::new(format!("colori dei gusci: dai minimi della distribuzione radiale delle palline ({}). \
            arancio = guscio interno, azzurro = esterno, viola = 3s", if conf.is_empty() { "—".into() } else { conf.join(", ") })).small().color(GREY));
        ui.separator();
        // ── numbers ──
        ui.label(RichText::new("NUMERI (dalle nostre corse, assets/atomo/risultati.txt)").strong().color(GOLD));
        let e_c = lab.r.e_c();
        let pass = |ok: bool| if ok { RichText::new("PASS").color(GREEN) } else { RichText::new("FAIL").color(RED) };
        egui::Grid::new("num").striped(true).show(ui, |ui| {
            ui.label("energia C (DMC)");
            ui.label(format!("{:.4} ± {:.4} Eₕ", e_c, lab.r.o("dmc_c_err", f64::NAN)));
            ui.label(format!("esatta {:.4} ({:+.2} %)", E_ESATTA_C, 100.0 * (e_c - E_ESATTA_C) / E_ESATTA_C.abs()));
            ui.label(pass(((e_c - E_ESATTA_C) / E_ESATTA_C).abs() < 0.01));
            ui.end_row();
            let ip = lv.ip * k::hartree_ev();
            ui.label("ionizzazione");
            ui.label(format!("{:.2} ± {:.2} eV", ip, lv.ip_err * k::hartree_ev()));
            ui.label(format!("misura {:.2} eV ({:+.1} %)", IP_NIST_EV, 100.0 * (ip - IP_NIST_EV) / IP_NIST_EV));
            ui.label(pass(((ip - IP_NIST_EV) / IP_NIST_EV).abs() < 0.03));
            ui.end_row();
            let ex = lv.de_ecc * k::hartree_ev();
            ui.label("2p → 3s");
            ui.label(format!("{:.2} ± {:.2} eV", ex, lv.de_ecc_err * k::hartree_ev()));
            ui.label(format!("misura {:.2} eV ({:+.1} %)", ECC_NIST_EV, 100.0 * (ex - ECC_NIST_EV) / ECC_NIST_EV));
            ui.label(pass(((ex - ECC_NIST_EV) / ECC_NIST_EV).abs() < 0.03));
            ui.end_row();
            ui.label("energia ora");
            ui.label(format!("{:.3} Eₕ", lab.e_media));
            ui.label(RichText::new("stima viva, rumorosa").color(GREY));
            ui.end_row();
            let ni = lab.nodo_info;
            ui.label("protone (nodo B=1)");
            ui.label(format!("B = {:.4}, E/12π² = {:.4}", ni[1], ni[2]));
            ui.label(format!("raggio {:.2} fm, massa {:.0} MeV", ni[3], ni[4]));
            ui.label("");
            ui.end_row();
            ui.label("elettrone zig-zag");
            ui.label(format!("tremolio ±{:.3} ħ/mc", lab.zbw_ampiezza));
            ui.label(format!("= ±{:.0} fm", lab.zbw_ampiezza * 386.159));
            ui.label("");
            ui.end_row();
        });
        ui.label(RichText::new(format!("fotoni: {} assorbiti · {} ionizzazioni · {} passati · {} emessi",
            lab.storia.frame(lab.vista).assorbiti, lab.storia.frame(lab.vista).ionizzati,
            lab.storia.frame(lab.vista).passati, lab.storia.frame(lab.vista).emessi)).small());
        ui.label(RichText::new(format!("{} palline disegnate · {:.0} fps · storia {} frame, {:.0} MB",
            lab.palline_disegnate, lab.fps, lab.storia.frames.len(), lab.storia.memoria_mb())).small().color(GREY));
        ui.separator();
        ui.label(RichText::new("DEBITI (messi a mano, da togliere)").strong().color(GOLD));
        for d in [
            "Nodi fissi: l'antisimmetria dei 6 elettroni viene dai nodi della guida (un determinante per spin). Errore piccolo ma non zero.",
            "Il nucleo del carbonio ha B = 12: qui si vede UN protone (nodo B = 1). Il nucleo intero non è calcolato.",
            "I due numeri del modello di Skyrme (F_π, e) sono presi da Adkins–Nappi–Witten (1983), non calcolati.",
            "L'elettrone zig-zag è libero: il pacchetto non sente il nucleo e la sua larghezza è scelta.",
            "Quanto è probabile l'assorbimento (forza della riga) e quanto dura lo stato eccitato non sono calcolati: l'attesa è compressa (~3 s).",
            "Diffusione elastica (Rayleigh) non calcolata: fuori risonanza il fotone passa.",
            "Temperatura: il nucleo si muove come un atomo di gas a 1 atm (urti come sfere dure). La foschia è solo visiva.",
            "Il tremolio delle palline della nuvola è la regola (tempo immaginario), non un moto vero: la nuvola vera di uno stato stazionario è ferma.",
        ] { ui.label(RichText::new(format!("• {d}")).small()); }
    });
    if !open { lab.panel_open = false; }
}

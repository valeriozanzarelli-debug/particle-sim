//! Interface: a clean screen (the atom, a time bar, a scale bar, ☰), everything else in one panel hidden at
//! start (Tab or ☰). Italian, simple words.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::{Color32, RichText};

use crate::camera::{Bersaglio, Cam};
use crate::{costanti as k, zz, Lab, Modo, Scena};

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

fn livello(m: f64) -> &'static str {
    if m > 3e-8 { "luce" } else if m > 2e-11 { "atomo" } else if m > 2e-13 { "tra gli elettroni" } else if m > 3e-15 { "vicino al nucleo" } else { "dentro il nucleo" }
}

pub fn draw_hud(mut ctx: EguiContexts, mut lab: ResMut<Lab>, mut cam: ResMut<Cam>, keys: Res<ButtonInput<KeyCode>>, auto: Res<crate::Auto>) {
    let lab = &mut *lab;
    let ctx = ctx.ctx_mut();
    if !ctx.wants_keyboard_input() && !auto.attivo {
        if keys.just_pressed(KeyCode::Tab) { lab.panel_open = !lab.panel_open; }
        if keys.just_pressed(KeyCode::Space) { lab.modo = if lab.modo == Modo::Avanti { Modo::Pausa } else { Modo::Avanti }; }
        if keys.just_pressed(KeyCode::KeyF) { lab.spara = true; }
        if keys.just_pressed(KeyCode::Digit1) { vai(lab, &mut cam, Bersaglio::Atomo); }
        if keys.just_pressed(KeyCode::Digit2) { vai(lab, &mut cam, Bersaglio::Nucleo); }
        if keys.just_pressed(KeyCode::Digit3) { vai(lab, &mut cam, Bersaglio::Elettrone); }
    }
    // ── top left: menu, the scenes, the zoom levels ──
    egui::Area::new(egui::Id::new("menu")).anchor(egui::Align2::LEFT_TOP, [10.0, 10.0]).show(ctx, |ui| {
        ui.horizontal(|ui| {
            let b = egui::Button::new(RichText::new("☰").size(16.0)).fill(Color32::from_black_alpha(90)).selected(lab.panel_open);
            if ui.add(b).on_hover_text("pannello (Tab) · pausa (Spazio) · spara un fotone (F) · 1 vista intera · 2 nucleo · 3 elettrone").clicked() {
                lab.panel_open = !lab.panel_open;
            }
            for (s, n) in [(Scena::Fenditura, "doppia fenditura"), (Scena::Atomo, "atomo zig-zag"), (Scena::Carbonio, "carbonio (DMC)")] {
                if ui.add(egui::Button::new(n).fill(Color32::from_black_alpha(90)).selected(lab.scena == s)).clicked() {
                    lab.cambia_scena(s);
                    let r = lab.distanza_tipica(Bersaglio::Atomo);
                    cam.vai(Bersaglio::Atomo, r);
                    cam.ln_r = cam.ln_r_target;
                }
            }
            ui.separator();
            let bt = cam.bersaglio;
            if lab.scena != Scena::Fenditura {
                if ui.add(egui::Button::new("1 intero").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Atomo)).clicked() { vai(lab, &mut cam, Bersaglio::Atomo); }
                if ui.add(egui::Button::new("2 nucleo").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Nucleo)).clicked() { vai(lab, &mut cam, Bersaglio::Nucleo); }
                if ui.add(egui::Button::new("3 elettrone").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Elettrone)).clicked() { vai(lab, &mut cam, Bersaglio::Elettrone); }
            }
        });
    });
    // ── top centre: what is happening ──
    let riga = stato_scena(lab);
    egui::Area::new(egui::Id::new("stato")).anchor(egui::Align2::CENTER_TOP, [0.0, 44.0]).show(ctx, |ui| {
        egui::Frame::none().fill(BG).rounding(6.0).inner_margin(egui::Margin::symmetric(10.0, 4.0)).show(ui, |ui| {
            ui.label(RichText::new(riga.0).monospace().color(CYAN));
            if !riga.1.is_empty() { ui.label(RichText::new(riga.1).small().color(GREY)); }
        });
    });
    // ── bottom left: the scale bar ──
    let mpu = lab.metri_per_unita();
    let bpp = cam.bohr_per_px();
    let len_m = tondo(140.0 * bpp * mpu);
    let px = (len_m / mpu / bpp) as f32;
    egui::Area::new(egui::Id::new("scala")).anchor(egui::Align2::LEFT_BOTTOM, [16.0, -70.0]).show(ctx, |ui| {
        let (rect, _) = ui.allocate_exact_size(egui::vec2(px.max(10.0) + 4.0, 26.0), egui::Sense::hover());
        let pnt = ui.painter();
        let y = rect.bottom() - 4.0;
        let x0 = rect.left() + 2.0;
        let stroke = egui::Stroke::new(2.0, Color32::from_rgb(220, 225, 235));
        pnt.line_segment([egui::pos2(x0, y), egui::pos2(x0 + px, y)], stroke);
        pnt.line_segment([egui::pos2(x0, y - 5.0), egui::pos2(x0, y + 1.0)], stroke);
        pnt.line_segment([egui::pos2(x0 + px, y - 5.0), egui::pos2(x0 + px, y + 1.0)], stroke);
        pnt.text(egui::pos2(x0, y - 8.0), egui::Align2::LEFT_BOTTOM, format!("{}   ·  {}", fmt_len(len_m), livello(cam.r() * mpu)),
                 egui::FontId::monospace(13.0), Color32::from_rgb(220, 225, 235));
    });
    barra_tempo(ctx, lab);
    if lab.panel_open {
        match lab.scena {
            Scena::Carbonio => pannello(ctx, lab, &mut cam),
            Scena::Fenditura => pannello_fenditura(ctx, lab),
            Scena::Atomo => pannello_atomo(ctx, lab, &mut cam),
        }
    }
}

fn vai(lab: &mut Lab, cam: &mut Cam, b: Bersaglio) {
    match b {
        Bersaglio::Nucleo => lab.mostra_protone = true,
        Bersaglio::Elettrone if lab.scena == Scena::Carbonio && (!lab.mostra_elettrone || lab.dirac.is_none()) => {
            lab.mostra_elettrone = true;
            lab.scegli_elettrone();
            lab.reset_dirac();
        }
        _ => {}
    }
    let r = lab.distanza_tipica(b);
    cam.vai(b, r);
}

fn stato_scena(lab: &Lab) -> (String, String) {
    match lab.scena {
        Scena::Carbonio => {
            let p = lab.storia.frame(lab.vista);
            (format!("CARBONIO (regola DMC, nuvola media) · {}", p.stato.nome()), p.registro.last().cloned().unwrap_or_default())
        }
        Scena::Fenditura => {
            let f = lab.fen.as_ref().unwrap();
            let fr = f.frange.map(|(m, e)| format!("prima frangia: {:.1} celle (geometria esatta {:.1}, {:+.0} %)", m, e, 100.0 * (m - e) / e)).unwrap_or_default();
            (format!("DOPPIA FENDITURA · luce = palline a spin 1 della regola zig-zag · tick {} · fotone n. {}", f.f.ticks, f.lanci), fr)
        }
        Scena::Atomo => {
            let o = lab.zz.as_ref().map(|z| z.lock().unwrap().out.clone());
            match o {
                Some(o) if o.pronto => (
                    format!("Nd⁵⁹⁺ (Z = 60, un elettrone) · regola zig-zag dal vivo · tick {} · nello stato 1s: {:.0} %", o.ticks, 100.0 * o.p1s),
                    format!("dipolo ⟨z⟩ = {:+.3} ħ/mc · campo del fotone {:+.4} · {:.0} ms per tick", o.dipolo, o.campo, o.ms_tick)),
                _ => ("Nd⁵⁹⁺: le palline stanno partendo...".into(), String::new()),
            }
        }
    }
}

fn barra_tempo(ctx: &egui::Context, lab: &mut Lab) {
    egui::Area::new(egui::Id::new("tempo")).anchor(egui::Align2::CENTER_BOTTOM, [0.0, -10.0]).show(ctx, |ui| {
        egui::Frame::none().fill(BG).rounding(8.0).inner_margin(egui::Margin::symmetric(12.0, 6.0)).show(ui, |ui| {
            ui.horizontal(|ui| {
                let back = ui.add(egui::Button::new("⏮").selected(lab.modo == Modo::Indietro)).on_hover_text("indietro nel tempo");
                if back.clicked() { lab.modo = if lab.modo == Modo::Indietro { Modo::Pausa } else { Modo::Indietro }; }
                let play = if lab.modo == Modo::Avanti { "⏸" } else { "▶" };
                if ui.button(play).on_hover_text("avanti / pausa (Spazio)").clicked() {
                    lab.modo = if lab.modo == Modo::Avanti { Modo::Pausa } else { Modo::Avanti };
                }
                match lab.scena {
                    Scena::Carbonio => {
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
                        let t = lab.storia.frame(lab.vista).t_s;
                        ui.label(RichText::new(format!("t = {}", fmt_t(t))).monospace().color(CYAN));
                    }
                    Scena::Fenditura => {
                        if ui.button("⏭").on_hover_text("4 tick per frame").clicked() { lab.passi_frame = if lab.passi_frame < 4.0 { 4.0 } else { 1.0 }; }
                        ui.add(egui::Slider::new(&mut lab.passi_frame, 0.1..=4.0).logarithmic(true).show_value(false)).on_hover_text("tick per frame");
                        let spt = lab.secondi_per_tick();
                        ui.label(RichText::new(format!("1 s = {}", fmt_t(spt * 60.0 * lab.passi_frame))).monospace().color(GREY));
                        let t = lab.fen.as_ref().map(|f| f.f.ticks).unwrap_or(0) as f64 * spt;
                        ui.label(RichText::new(format!("t = {}", fmt_t(t))).monospace().color(CYAN));
                    }
                    Scena::Atomo => {
                        if ui.button("⏭").on_hover_text("più veloce").clicked() { lab.ritmo = if lab.ritmo < 60.0 { 60.0 } else { 25.0 }; }
                        ui.add(egui::Slider::new(&mut lab.ritmo, 1.0..=60.0).logarithmic(true).show_value(false)).on_hover_text("tick al secondo (al massimo quanto regge il PC)");
                        let spt = lab.secondi_per_tick();
                        let ticks = lab.zz.as_ref().map(|z| z.lock().unwrap().out.ticks).unwrap_or(0);
                        ui.label(RichText::new(format!("1 s ≤ {}", fmt_t(spt * lab.ritmo))).monospace().color(GREY));
                        ui.label(RichText::new(format!("t = {}", fmt_t(ticks as f64 * spt))).monospace().color(CYAN));
                    }
                }
                ui.separator();
                ui.label("T");
                ui.add(egui::Slider::new(&mut lab.t_kelvin, 0.0..=3000.0).suffix(" K").max_decimals(0));
            });
        });
    });
}

fn pannello_fenditura(ctx: &egui::Context, lab: &mut Lab) {
    let mut open = true;
    egui::Window::new("Doppia fenditura (regola zig-zag)").open(&mut open).default_pos([12.0, 80.0]).default_width(440.0).show(ctx, |ui| {
        ui.label(RichText::new("La prima prova del laboratorio nuovo: il fotone fatto di palline a spin 1 (F = E + icB) che \
            avanzano alla velocità della luce, la stessa regola dell'elettrone senza massa.").color(GOLD));
        ui.label("Un tick: ogni pallina divide il suo stato secondo lo spin lungo x, y, z; la parte +1 va avanti di una cella, \
                  la parte −1 indietro, la parte 0 resta. Il muro assorbe; ai bordi una spugna assorbe ciò che arriva.");
        ui.separator();
        let f = lab.fen.as_ref().unwrap();
        ui.label(format!("λ = {} celle, fenditure a {:.0} celle, larghe {:.0}; schermo a {} celle dal muro",
            f.f.lambda, f.f.fend[1] - f.f.fend[0], f.f.larghezza, f.f.x_schermo - f.f.x_muro));
        if let Some((m, e)) = f.frange {
            let dev = (m - e) / e;
            ui.label(RichText::new(format!("prima frangia luminosa: misurata {:.1} celle dall'asse, geometria esatta {:.1} ({:+.1} %)", m, e, 100.0 * dev))
                .color(if dev.abs() < 0.05 { GREEN } else { RED }));
            ui.label(RichText::new("Il reticolo non è isotropo: in diagonale la luce delle palline è più lenta, le frange stanno più strette. \
                Misurato: −25 % con 8 celle per λ, −13 % con 12, −11 % con 16 (si avvicina, piano).").small().color(GREY));
        } else {
            ui.label(RichText::new("le frange si misurano quando il fotone ha passato lo schermo").small().color(GREY));
        }
        ui.label(RichText::new("La luce del laboratorio del fotone (500 nm) è disegnata con 12 celle per λ; le fenditure sono lunghe \
            (il problema è piano): ne disegniamo un pezzo.").small().color(GREY));
        ui.separator();
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut lab.grandezza, 0.3..=3.0).text("palline"));
            ui.add(egui::Slider::new(&mut lab.luce, 0.2..=4.0).text("luce"));
            ui.checkbox(&mut lab.bagliore, "bagliore");
        });
        ui.label(RichText::new(format!("{} palline · {:.0} fps", lab.palline_disegnate, lab.fps)).small().color(GREY));
    });
    if !open { lab.panel_open = false; }
}

fn pannello_atomo(ctx: &egui::Context, lab: &mut Lab, cam: &mut Cam) {
    let mut open = true;
    let o = lab.zz.as_ref().map(|z| z.lock().unwrap().out.clone());
    egui::Window::new("Atomo zig-zag (dal vivo)").open(&mut open).default_pos([12.0, 80.0]).default_width(470.0).vscroll(true).show(ctx, |ui| {
        ui.label(RichText::new("Un elettrone legato a un nucleo dalla sola fase di Coulomb, con la stessa regola della luce \
            più la massa: ogni tick le palline avanzano a c e girano (L ↔ R). Gira dal vivo, e all'indietro esattamente.").color(GOLD));
        ui.label(RichText::new("Perché neodimio (Z = 60) e non idrogeno: il raggio dell'idrogeno è 137 ħ/mc e le celle devono essere \
            ≲ 0,5 ħ/mc: servirebbero ~10¹⁰ celle. Il problema dipende solo da Zα: con Z = 60 l'atomo sta in 48³ celle.").small().color(GREY));
        ui.separator();
        ui.label(RichText::new("FOTONE").strong().color(GOLD));
        if let Some(o) = &o {
            let res = o.omega_res * zz::MC2_KEV;
            if lab.e_kev <= 0.0 { lab.e_kev = res; }
            ui.horizontal(|ui| {
                ui.label("energia");
                ui.add(egui::Slider::new(&mut lab.e_kev, 10.0..=150.0).logarithmic(true).suffix(" keV").max_decimals(1));
            });
            ui.horizontal(|ui| {
                if ui.button(format!("1s → 2p ({:.1} keV, dalla regola)", res)).clicked() { lab.e_kev = res; }
                if ui.button("fuori risonanza (×1,3)").clicked() { lab.e_kev = 1.3 * res; }
                if ui.add(egui::Button::new(RichText::new("  Spara (F)  ").strong()).fill(Color32::from_rgb(40, 60, 110))).clicked() { lab.spara = true; }
            });
            let lam_pm = 1239.841_984 / (lab.e_kev * 1e3) * 1e3;
            ui.label(RichText::new(format!("λ = {:.1} pm = {:.0} ħ/mc: 50 volte l'atomo. Il pacchetto è disegnato 6 volte più corto; \
                all'atomo il suo campo entra come fase e E·r (forma di dipolo, giusta perché λ ≫ atomo).", lam_pm, lam_pm * 1e-12 / zz::LAMBDA_C_M)).small().color(GREY));
            ui.label(RichText::new("Se l'energia è quella giusta, la regola stessa sposta l'elettrone verso il 2p: la nuvola resta in una \
                SOVRAPPOSIZIONE 1s + 2p che oscilla su e giù (il dipolo), e quel dipolo manda fuori la luce (gusci oro/azzurro).").small());
            ui.separator();
            ui.label(RichText::new("NUMERI").strong().color(GOLD));
            egui::Grid::new("zz").striped(true).show(ui, |ui| {
                ui.label("legame 1s (regola radiale, a = 0,5)");
                ui.label(format!("{:.4} mc² = {:.1} keV", o.legame_1s, o.legame_1s * zz::MC2_KEV));
                ui.end_row();
                ui.label("Dirac esatto, Z = 60");
                let za = zz::z_alpha();
                ui.label(format!("{:.4} mc² ({:+.0} %: celle grosse)", 1.0 - (1.0 - za * za).sqrt(), 100.0 * (o.legame_1s / (1.0 - (1.0 - za * za).sqrt()) - 1.0)));
                ui.end_row();
                ui.label("nello stato 1s");
                ui.label(format!("{:.1} %", 100.0 * o.p1s));
                ui.end_row();
                ui.label("dipolo ⟨z⟩");
                ui.label(format!("{:+.3} ħ/mc", o.dipolo));
                ui.end_row();
                ui.label("costo");
                ui.label(format!("{:.0} ms per tick · {} palline · {:.0} fps", o.ms_tick, lab.palline_disegnate, lab.fps));
                ui.end_row();
            });
            ui.label(RichText::new("Lo stato di partenza è l'1s della forma radiale della regola messo sul reticolo cubico: su celle così \
                grosse è solo ~60 % l'1s della regola 3-D, il resto si vede come tremolio in più.").small().color(GREY));
        }
        ui.separator();
        ui.label(RichText::new("NUCLEO (nodo B = 1)").strong().color(GOLD));
        ui.checkbox(&mut lab.mostra_protone, "mostra protone vero");
        let ni = lab.nodo_info;
        let n2 = lab.nodo_info2;
        ui.label(format!("B = {:.4}, E/12π² = {:.4}, raggio {:.2} fm, massa {:.0} MeV", ni[1], ni[2], ni[3], ni[4]));
        ui.label(format!("gira (isorotazione, spin ½): Λ → M_Δ − M_N = {:.0} MeV (misura 293); respira: ħω = {:.0} MeV (stima di Derrick)", n2[1], n2[4]));
        ui.label(RichText::new(format!("orologio del nucleo rallentato: 1 s di schermo = {:.1e} s", lab.rallentamento_nucleo())).small().color(GREY));
        if ui.button("vai al nucleo (2)").clicked() { vai(lab, cam, Bersaglio::Nucleo); }
        ui.separator();
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut lab.grandezza, 0.3..=3.0).text("palline"));
            ui.add(egui::Slider::new(&mut lab.luce, 0.2..=4.0).text("luce"));
            ui.checkbox(&mut lab.bagliore, "bagliore");
        });
        if ui.button("ricomincia dall'1s").clicked() {
            if let Some(z) = &lab.zz { z.lock().unwrap().ricomincia = true; }
        }
        ui.label(RichText::new("DEBITI: il campo del fotone è un impulso forte (come un laser), non un fotone solo; l'emissione \
            spontanea (vuoto quantistico) non c'è: la sovrapposizione non decade da sola. Il nucleo di Nd ha B = 144: si vede un protone. \
            La temperatura qui non fa nulla di visibile (in 10⁻¹⁹ s un nucleo a 300 K si sposta di 10⁻⁵ fm).").small());
    });
    if !open { lab.panel_open = false; }
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
            if ui.button("vai al nucleo (2)").clicked() { vai(lab, cam, Bersaglio::Nucleo); }
            if ui.button("vai all'elettrone (3)").clicked() { vai(lab, cam, Bersaglio::Elettrone); }
            if ui.button("atomo (1)").clicked() { vai(lab, cam, Bersaglio::Atomo); }
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

//! Interface: a clean screen (the scene, a time bar, a scale bar, ☰), everything else in one panel hidden at
//! start (Tab or ☰). Italian, simple words.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::{Color32, RichText};

use crate::camera::{Bersaglio, Cam};
use crate::fotone3d::Polar;
use crate::fotscena::{self, DIREZIONI};
use crate::{zz, Lab, Modo, Scena};

const CYAN: Color32 = Color32::from_rgb(90, 220, 255);
const GOLD: Color32 = Color32::from_rgb(255, 205, 80);
const GREY: Color32 = Color32::from_rgb(150, 155, 170);
const GREEN: Color32 = Color32::from_rgb(140, 230, 140);
const RED: Color32 = Color32::from_rgb(255, 120, 100);
const BG: Color32 = Color32::from_rgba_premultiplied(8, 10, 18, 190);

/// hc in eV·nm (CODATA) and c.
pub const HC_EV_NM: f64 = 1239.841_984;
pub const C_M_S: f64 = 299_792_458.0;

/// The name of the band of a wavelength (m).
pub fn banda(l: f64) -> &'static str {
    if l > 1.0 { "onde radio" } else if l > 1e-3 { "microonde" } else if l > 7.5e-7 { "infrarosso" } else if l > 3.8e-7 { "luce visibile" }
    else if l > 1e-8 { "ultravioletto" } else if l > 1e-11 { "raggi X" } else { "raggi gamma" }
}

pub fn fmt_freq(f: f64) -> String {
    let (v, u) = if f >= 1e18 { (f / 1e18, "EHz") } else if f >= 1e15 { (f / 1e15, "PHz") } else if f >= 1e12 { (f / 1e12, "THz") }
                 else if f >= 1e9 { (f / 1e9, "GHz") } else if f >= 1e6 { (f / 1e6, "MHz") } else { (f / 1e3, "kHz") };
    format!("{} {}", nice(v), u)
}

pub fn fmt_energia(ev: f64) -> String {
    if ev >= 1e6 { format!("{} MeV", nice(ev / 1e6)) } else if ev >= 1e3 { format!("{} keV", nice(ev / 1e3)) } else if ev >= 1.0 { format!("{} eV", nice(ev)) }
    else if ev >= 1e-3 { format!("{} meV", nice(ev * 1e3)) } else { format!("{} µeV", nice(ev * 1e6)) }
}

pub fn fmt_t(s: f64) -> String {
    let a = s.abs();
    let (v, u) = if a < 1e-15 { (s * 1e18, "as") } else if a < 1e-12 { (s * 1e15, "fs") } else if a < 1e-9 { (s * 1e12, "ps") }
                 else if a < 1e-6 { (s * 1e9, "ns") } else { (s * 1e6, "µs") };
    format!("{:.3} {}", v, u)
}

pub fn fmt_len(m: f64) -> String {
    let a = m.abs();
    if a >= 1e3 { format!("{} km", nice(m * 1e-3)) } else if a >= 1.0 { format!("{} m", nice(m)) } else if a >= 1e-2 { format!("{} cm", nice(m * 1e2)) }
    else if a >= 1e-3 { format!("{} mm", nice(m * 1e3)) } else if a >= 1e-6 { format!("{} µm", nice(m * 1e6)) } else if a >= 1e-9 { format!("{} nm", nice(m * 1e9)) }
    else if a >= 1e-12 { format!("{} pm", nice(m * 1e12)) } else { format!("{} fm", nice(m * 1e15)) }
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
    if m > 3e-8 { "luce" } else if m > 2e-12 { "atomo" } else if m > 3e-15 { "vicino al nucleo" } else { "dentro il nucleo" }
}

pub fn draw_hud(mut ctx: EguiContexts, mut lab: ResMut<Lab>, mut cam: ResMut<Cam>, keys: Res<ButtonInput<KeyCode>>, auto: Res<crate::Auto>) {
    let lab = &mut *lab;
    let ctx = ctx.ctx_mut();
    if !ctx.wants_keyboard_input() && !auto.attivo {
        if keys.just_pressed(KeyCode::Tab) { lab.panel_open = !lab.panel_open; }
        if keys.just_pressed(KeyCode::Space) { lab.modo = if lab.modo == Modo::Avanti { Modo::Pausa } else { Modo::Avanti }; }
        if keys.just_pressed(KeyCode::KeyF) { lab.spara = true; }
        if keys.just_pressed(KeyCode::KeyM) { lab.misura = true; }
        if keys.just_pressed(KeyCode::Digit1) { vai(lab, &mut cam, Bersaglio::Atomo); }
        if lab.scena != Scena::Fotone {
            if keys.just_pressed(KeyCode::Digit2) { vai(lab, &mut cam, Bersaglio::Nucleo); }
            if keys.just_pressed(KeyCode::Digit3) { vai(lab, &mut cam, Bersaglio::Elettrone); }
        }
    }
    // ── top left: menu, the scenes, the zoom levels ──
    egui::Area::new(egui::Id::new("menu")).anchor(egui::Align2::LEFT_TOP, [10.0, 10.0]).show(ctx, |ui| {
        ui.horizontal(|ui| {
            let b = egui::Button::new(RichText::new("☰").size(16.0)).fill(Color32::from_black_alpha(90)).selected(lab.panel_open);
            if ui.add(b).on_hover_text("pannello (Tab) · pausa (Spazio) · lancia un fotone (F) · 1 vista intera · 2 nucleo · 3 elettrone").clicked() {
                lab.panel_open = !lab.panel_open;
            }
            for (s, n) in [(Scena::Fotone, "fotone libero"), (Scena::Atomo, "atomo zig-zag"), (Scena::Orbitali, "orbitali")] {
                if ui.add(egui::Button::new(n).fill(Color32::from_black_alpha(90)).selected(lab.scena == s)).clicked() && lab.scena != s {
                    lab.cambia_scena(s);
                    if s == Scena::Fotone { cam.yaw = 0.42; cam.pitch = 0.28; }
                    let r = lab.distanza_tipica(Bersaglio::Atomo);
                    cam.vai(Bersaglio::Atomo, r);
                    cam.ln_r = cam.ln_r_target;
                }
            }
            ui.separator();
            let bt = cam.bersaglio;
            if lab.scena != Scena::Fotone {
                if ui.add(egui::Button::new("1 intero").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Atomo)).clicked() { vai(lab, &mut cam, Bersaglio::Atomo); }
                if ui.add(egui::Button::new("2 nucleo").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Nucleo)).clicked() { vai(lab, &mut cam, Bersaglio::Nucleo); }
                if ui.add(egui::Button::new("3 elettrone").fill(Color32::from_black_alpha(90)).selected(bt == Bersaglio::Elettrone)).clicked() { vai(lab, &mut cam, Bersaglio::Elettrone); }
            } else {
                ui.checkbox(&mut lab.segui, RichText::new("segui il fotone").color(Color32::from_rgb(220, 225, 235)));
            }
        });
    });
    // ── top centre: what is happening ──
    let mut riga = stato_scena(lab);
    if let Some((n, Some(p))) = clic_scena(lab) {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() * lab.metri_per_unita();
        riga.1 = format!("{} · MISURA: clic n. {} a ({:.1}, {:.1}, {:.1}) unità, {} dal centro{}", riga.1, n, p[0], p[1], p[2], fmt_len(r),
            if lab.opz.collassa { " (collasso)" } else { " (solo registrata)" });
    }
    egui::Area::new(egui::Id::new("stato")).anchor(egui::Align2::CENTER_TOP, [0.0, 44.0]).show(ctx, |ui| {
        egui::Frame::none().fill(BG).rounding(6.0).inner_margin(egui::Margin::symmetric(10.0, 4.0)).show(ui, |ui| {
            ui.label(RichText::new(riga.0).monospace().color(CYAN));
            if !riga.1.is_empty() { ui.label(RichText::new(riga.1).monospace().color(GREY)); }
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
    // the orbital view: when a new state set changes the size, the camera follows
    if lab.scena == Scena::Orbitali {
        if let Some(r) = lab.orb.as_ref().and_then(|o| { let g = o.lock().unwrap(); if g.out.pronto { Some(g.out.r) } else { None } }) {
            if (r - lab.orb_r).abs() > 1e-6 { lab.orb_r = r; if cam.bersaglio == Bersaglio::Atomo { cam.vai(Bersaglio::Atomo, 2.4 * r); } }
        }
    }
    barra_tempo(ctx, lab);
    if lab.panel_open {
        match lab.scena {
            Scena::Fotone => pannello_fotone(ctx, lab),
            Scena::Atomo => pannello_atomo(ctx, lab, &mut cam),
            Scena::Orbitali => pannello_orbitali(ctx, lab),
        }
    }
}

fn vai(lab: &mut Lab, cam: &mut Cam, b: Bersaglio) {
    if b == Bersaglio::Nucleo { lab.mostra_protone = true; }
    let r = lab.distanza_tipica(b);
    cam.vai(b, r);
}

/// Clicks of the detector in the scene shown: (how many, the last one).
fn clic_scena(lab: &Lab) -> Option<(usize, Option<[f64; 3]>)> {
    match lab.scena {
        Scena::Fotone => lab.fot.as_ref().map(|z| { let o = &z.lock().unwrap().out; (o.n_clic, o.ultimo_clic) }),
        Scena::Atomo => lab.zz.as_ref().map(|z| { let o = &z.lock().unwrap().out; (o.n_clic, o.ultimo_clic) }),
        Scena::Orbitali => lab.orb.as_ref().map(|z| { let o = &z.lock().unwrap().out; (o.n_clic, o.ultimo_clic) }),
    }
}

/// The MEASUREMENT block of every panel.
fn pannello_misura(ui: &mut egui::Ui, lab: &mut Lab) {
    ui.separator();
    ui.label(RichText::new("MISURA E SURFISTI").strong().color(GOLD));
    ui.horizontal(|ui| {
        if ui.add(egui::Button::new(RichText::new("  MISURA (M)  ").strong()).fill(Color32::from_rgb(40, 100, 60))).clicked() { lab.misura = true; }
        if ui.button("cancella i clic").clicked() { lab.azzera_clic = true; }
        ui.checkbox(&mut lab.opz.mostra_clic, "mostra i clic passati");
    });
    let etichetta = match lab.scena {
        Scena::Fotone => "dopo la misura: il fotone è assorbito (ne parte un altro)",
        Scena::Atomo => "dopo la misura: collassa (un pacchetto stretto nel punto del clic)",
        Scena::Orbitali => "dopo la misura: collassa (non disponibile qui: usa la scena atomo)",
    };
    ui.add_enabled(lab.scena != Scena::Orbitali, egui::Checkbox::new(&mut lab.opz.collassa, etichetta));
    if !lab.opz.collassa || lab.scena == Scena::Orbitali {
        ui.label(RichText::new("Senza collasso il clic è solo registrato: una vista non distruttiva, non una misura vera (una misura \
            vera cambia l'onda).").small().color(GREY));
    }
    ui.checkbox(&mut lab.opz.surfisti, "surfisti (traiettorie di Bohm, interpretazione)");
    ui.label(RichText::new("Una misura di posizione dà UN clic, in un punto scelto con probabilità |ψ|² in quell'istante (per il fotone: \
        l'energia |F|²). Tanti clic, uno dopo l'altro, rifanno la forma dell'onda. I surfisti sono palline portate dalla corrente \
        dell'onda (de Broglie–Bohm): un'interpretazione, non servono alla regola; le palline dell'onda sono l'elettrone (o la luce).").small().color(GREY));
}

fn stato_scena(lab: &Lab) -> (String, String) {
    match lab.scena {
        Scena::Fotone => {
            let o = lab.fot.as_ref().map(|z| z.lock().unwrap().out.clone());
            match o {
                Some(o) if o.pronto => {
                    let lam = lab.lambda_m;
                    let v = o.v.map(|v| format!("{:.4} c", v)).unwrap_or_else(|| "…".into());
                    let h = o.elicita.map(|h| format!("{:+.3}", h)).unwrap_or_else(|| "…".into());
                    (format!("UN FOTONE di palline a spin 1 · {} · λ = {} · ν = {} · energia {} · velocità misurata {} · elicità {}",
                             banda(lam), fmt_len(lam), fmt_freq(C_M_S / lam), fmt_energia(HC_EV_NM * 1e-9 / lam), v, h),
                     format!("{} · {:.0} celle per λ · larghezza del fascio {} · tick {} · fotone n. {}",
                             o.par.polar.nome(), o.par.cpl(), fmt_len(o.w * lab.cella_m()), o.ticks, o.lanci))
                }
                _ => ("Il fotone sta partendo...".into(), String::new()),
            }
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
        Scena::Orbitali => {
            let o = lab.orb.as_ref().map(|z| z.lock().unwrap().out.clone());
            match o {
                Some(o) if o.pronto => {
                    let nomi: Vec<String> = o.livelli.iter().map(|l| l.0.nome()).collect();
                    let leg: Vec<String> = o.livelli.iter().map(|l| format!("{} {:.2} keV", l.0.nome(), (1.0 - l.1) * zz::MC2_KEV)).collect();
                    (format!("Nd⁵⁹⁺ · ORBITALI dalla forma radiale della regola zig-zag · {}", nomi.join(" + ")),
                     format!("legami {} · {} · t = {}", leg.join(", "),
                        if o.livelli.len() > 1 { format!("battimento ogni {}", fmt_t(o.periodo * zz::T_UNIT_S)) } else { format!("la fase gira ogni {}", fmt_t(o.periodo * zz::T_UNIT_S)) },
                        fmt_t(o.t * zz::T_UNIT_S)))
                }
                _ => ("Orbitali: la regola radiale sta ascoltando i livelli...".into(), String::new()),
            }
        }
    }
}

fn barra_tempo(ctx: &egui::Context, lab: &mut Lab) {
    egui::Area::new(egui::Id::new("tempo")).anchor(egui::Align2::CENTER_BOTTOM, [0.0, -10.0]).show(ctx, |ui| {
        egui::Frame::none().fill(BG).rounding(8.0).inner_margin(egui::Margin::symmetric(12.0, 6.0)).show(ui, |ui| {
            ui.horizontal(|ui| {
                let back = ui.add(egui::Button::new("⏮").selected(lab.modo == Modo::Indietro)).on_hover_text("indietro nel tempo (esatto)");
                if back.clicked() { lab.modo = if lab.modo == Modo::Indietro { Modo::Pausa } else { Modo::Indietro }; }
                let play = if lab.modo == Modo::Avanti { "⏸" } else { "▶" };
                if ui.button(play).on_hover_text("avanti / pausa (Spazio)").clicked() {
                    lab.modo = if lab.modo == Modo::Avanti { Modo::Pausa } else { Modo::Avanti };
                }
                let spt = lab.secondi_per_tick();
                match lab.scena {
                    Scena::Fotone => {
                        if ui.button("⏭").on_hover_text("più veloce").clicked() { lab.ritmo_fot = if lab.ritmo_fot < 60.0 { 60.0 } else { 20.0 }; }
                        ui.add(egui::Slider::new(&mut lab.ritmo_fot, 1.0..=60.0).logarithmic(true).show_value(false)).on_hover_text("tick al secondo (al massimo quanto regge il PC)");
                        let ticks = lab.fot.as_ref().map(|z| z.lock().unwrap().out.ticks).unwrap_or(0);
                        ui.label(RichText::new(format!("1 s ≤ {}", fmt_t(spt * lab.ritmo_fot))).monospace().color(GREY));
                        ui.label(RichText::new(format!("t = {}", fmt_t(ticks as f64 * spt))).monospace().color(CYAN));
                        if ui.add(egui::Button::new(RichText::new(" Lancia (F) ").strong()).fill(Color32::from_rgb(40, 60, 110))).clicked() { lab.spara = true; }
                    }
                    Scena::Orbitali => {
                        ui.add(egui::Slider::new(&mut lab.ritmo_orb, 1.0..=60.0).logarithmic(true).show_value(false)).on_hover_text("passi al secondo");
                        let t = lab.orb.as_ref().map(|z| z.lock().unwrap().out.t).unwrap_or(0.0);
                        ui.label(RichText::new(format!("t = {}", fmt_t(t * zz::T_UNIT_S))).monospace().color(CYAN));
                    }
                    Scena::Atomo => {
                        if ui.button("⏭").on_hover_text("più veloce").clicked() { lab.ritmo = if lab.ritmo < 60.0 { 60.0 } else { 25.0 }; }
                        ui.add(egui::Slider::new(&mut lab.ritmo, 1.0..=60.0).logarithmic(true).show_value(false)).on_hover_text("tick al secondo (al massimo quanto regge il PC)");
                        let ticks = lab.zz.as_ref().map(|z| z.lock().unwrap().out.ticks).unwrap_or(0);
                        ui.label(RichText::new(format!("1 s ≤ {}", fmt_t(spt * lab.ritmo))).monospace().color(GREY));
                        ui.label(RichText::new(format!("t = {}", fmt_t(ticks as f64 * spt))).monospace().color(CYAN));
                    }
                }
                ui.separator();
                if ui.add(egui::Button::new(RichText::new(" MISURA (M) ").strong()).fill(Color32::from_rgb(40, 100, 60))).on_hover_text("un clic del rivelatore, con probabilità |ψ|²").clicked() { lab.misura = true; }
            });
        });
    });
}

fn pannello_fotone(ctx: &egui::Context, lab: &mut Lab) {
    let mut open = true;
    let o = lab.fot.as_ref().map(|z| z.lock().unwrap().out.clone());
    egui::Window::new("Un fotone libero (regola zig-zag)").open(&mut open).default_pos([12.0, 80.0]).default_width(470.0).vscroll(true).show(ctx, |ui| {
        ui.label(RichText::new("Un solo fotone fatto di palline a spin 1 (F = E + icB), in una scatola vuota in 3-D. \
            Ogni tick ogni pallina divide il suo stato secondo lo spin lungo x, y, z: la parte +1 avanza di una cella, \
            la parte −1 torna indietro, la parte 0 resta. Nient'altro: è la stessa regola dell'elettrone senza massa.").color(GOLD));
        ui.label(RichText::new("Colore = la direzione del campo elettrico E nel piano di traverso al volo (una ruota di colori). \
            Luminosità = E² (le creste). Polarizzazione lineare: due colori che si alternano cresta dopo cresta. \
            Circolare: una vite colorata che gira mentre vola (a destra o a sinistra secondo l'elicità).").small());
        ui.separator();
        ui.label(RichText::new("IL FOTONE DA LANCIARE").strong().color(GOLD));
        ui.horizontal(|ui| {
            ui.label("λ vera");
            let mut lg = lab.lambda_m.log10();
            ui.add(egui::Slider::new(&mut lg, -12.0..=3.0).show_value(false)).on_hover_text("dalle onde radio (km) ai raggi gamma (pm)");
            lab.lambda_m = 10f64.powf(lg);
            ui.label(RichText::new(format!("{} · {} · {} · {}", fmt_len(lab.lambda_m), banda(lab.lambda_m), fmt_freq(C_M_S / lab.lambda_m),
                fmt_energia(HC_EV_NM * 1e-9 / lab.lambda_m))).color(CYAN));
        });
        ui.horizontal(|ui| {
            for (t, l) in [("radio 1 m", 1.0), ("microonde 1 cm", 1e-2), ("visibile 500 nm", 500e-9), ("raggi X 0,1 nm", 1e-10), ("gamma 1 pm", 1e-12)] {
                if ui.small_button(t).clicked() { lab.lambda_m = l; }
            }
        });
        ui.label(RichText::new("La regola senza massa non ha una lunghezza sua: un fotone radio e un fotone gamma sono LO STESSO \
            disegno di palline, cambia solo quanto è grande una cella (la scala in basso a sinistra) e quanto dura un tick. \
            Cambiare λ qui non rilancia nulla: cambia solo la scala. L'unica cosa che cambia davvero con λ è come la luce \
            parla con la materia (qui non c'è materia).").small().color(GREY));
        let p = &mut lab.fot_par;
        ui.horizontal(|ui| {
            ui.label("qualità");
            egui::ComboBox::from_id_salt("qual").selected_text(fotscena::QUALITA[p.qual].2).show_ui(ui, |ui| {
                for (i, q) in fotscena::QUALITA.iter().enumerate() { ui.selectable_value(&mut p.qual, i, q.2); }
            });
        });
        ui.horizontal(|ui| {
            ui.label("larghezza al fuoco w0");
            ui.add(egui::Slider::new(&mut p.w0_lambda, 0.7..=3.0).suffix(" λ").max_decimals(2));
        });
        ui.horizontal(|ui| {
            ui.label("polarizzazione");
            egui::ComboBox::from_id_salt("pol").selected_text(p.polar.nome()).show_ui(ui, |ui| {
                for q in [Polar::Orizzontale, Polar::Verticale, Polar::CircolarePiu, Polar::CircolareMeno] { ui.selectable_value(&mut p.polar, q, q.nome()); }
            });
        });
        ui.horizontal(|ui| {
            ui.label("direzione");
            egui::ComboBox::from_id_salt("dir").selected_text(DIREZIONI[p.dir % DIREZIONI.len()].0).show_ui(ui, |ui| {
                for (i, (n, _)) in DIREZIONI.iter().enumerate() { ui.selectable_value(&mut p.dir, i, *n); }
            });
        });
        if ui.add(egui::Button::new(RichText::new("  Lancia (F)  ").strong()).fill(Color32::from_rgb(40, 60, 110))).clicked() { lab.spara = true; }
        ui.label(RichText::new(format!("La scatola: {n}³ celle = {l:.0} λ per lato, periodica e vuota. Il fotone parte lontano dalle \
            facce, si stringe fino al fuoco al centro e poi si allarga (diffrazione); quando ha attraversato ne parte un altro uguale. \
            Qualità e forma valgono dal prossimo lancio.", n = p.n(), l = p.n() as f64 / p.cpl())).small().color(GREY));
        ui.separator();
        ui.label(RichText::new("MISURATO DAL VIVO (dal campo, non disegnato)").strong().color(GOLD));
        if let Some(o) = &o {
            let lam = o.par.cpl();
            let w0 = o.par.w0_lambda * lam;
            let th_att = lam / (std::f64::consts::PI * w0);
            egui::Grid::new("fot").striped(true).show(ui, |ui| {
                ui.label("velocità del centro dell'energia");
                match o.v {
                    Some(v) => {
                        ui.label(format!("{:.4} c", v));
                        ui.label(RichText::new(format!("onde continue, stesso impulso: {:.4} c", o.v_continuo)).color(GREY));
                        ui.label(if ((v - o.v_continuo) / o.v_continuo).abs() < 0.005 { RichText::new("entro 0,5 %").color(GREEN) } else { RichText::new("fuori 0,5 %").color(RED) });
                    }
                    None => { ui.label("…"); ui.label(""); ui.label(""); }
                }
                ui.end_row();
                ui.label("elicità");
                match o.elicita {
                    Some(h) => {
                        ui.label(format!("{:+.3}", h));
                        ui.label(RichText::new(format!("attesa {:+.0}", o.par.polar.elicita())).color(GREY));
                        ui.label(if (h - o.par.polar.elicita()).abs() < 0.02 { RichText::new("ok").color(GREEN) } else { RichText::new("no").color(RED) });
                    }
                    None => { ui.label("…"); ui.label(""); ui.label(""); }
                }
                ui.end_row();
                ui.label("apertura del fascio θ");
                match o.theta {
                    Some(t) => {
                        ui.label(format!("{:.3} rad", t));
                        let ka = std::f64::consts::TAU / lam;
                        ui.label(RichText::new(format!("λ/(π w0) = {:.3} · la regola: {:.3}", th_att, th_att * ka / ka.tan())).color(GREY))
                            .on_hover_text("di lato la regola trasporta la luce con v⊥ = q⊥/tan(ka) invece di q⊥/(ka): l'apertura è più piccola di ka/tan(ka)");
                        ui.label(if ((t - th_att) / th_att).abs() < 0.05 { RichText::new("entro 5 %").color(GREEN) } else { RichText::new("fuori 5 %").color(RED) });
                    }
                    None => { ui.label("…"); ui.label(""); ui.label(""); }
                }
                ui.end_row();
                ui.label("larghezza ora");
                ui.label(format!("{:.2} λ = {}", o.w / lam, fmt_len(o.w * lab.cella_m())));
                ui.label(RichText::new(format!("al fuoco w0 = {:.2} λ", w0 / lam)).color(GREY));
                ui.label("");
                ui.end_row();
            });
            ui.label(RichText::new("La velocità è quella del centro dell'energia |F|²: un impulso stretto contiene onde un po' oblique, \
                quindi lungo l'asse va appena meno di c (anche nelle onde continue: il confronto giusto è con la colonna grigia). \
                Lungo un asse del reticolo la regola porta la luce a c esatto; in diagonale il reticolo la rallenta (debito).").small().color(GREY));
            ui.label(RichText::new(format!("costo: {:.0} ms per tick · {} palline · {:.0} fps", o.ms_tick, lab.palline_disegnate, lab.fps)).small().color(GREY));
        }
        ui.separator();
        ui.label(RichText::new("STRATI (tutti dallo stesso campo F = E + icB delle palline)").strong().color(GOLD));
        let v = &mut lab.vista_fot;
        ui.horizontal(|ui| {
            ui.checkbox(&mut v.traccianti, "palline che corrono col campo");
            ui.checkbox(&mut v.fronti, RichText::new("fronti d'onda").color(GOLD));
            ui.checkbox(&mut v.elica, "elica di E");
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut v.e_vettori, "tratti E");
            ui.checkbox(&mut v.b_vettori, RichText::new("tratti B").color(Color32::from_rgb(110, 160, 255)));
            ui.checkbox(&mut v.poynting, RichText::new("Poynting S = E×B").color(GOLD));
            ui.checkbox(&mut v.spin, RichText::new("spin").color(GREEN));
            ui.checkbox(&mut v.sezione, "sezione");
        });
        ui.label(RichText::new("Palline: nascono dove il fotone ha energia (più energia, più palline) e a ogni tick si muovono con \
            la velocità dell'energia del campo, v = 2 E×cB/(E² + c²B²) (c in un'onda piana); colore = direzione di E, luce = E². \
            Fronti d'onda: fogli sottili dove E lungo la polarizzazione ha una cresta (oro) o una valle (azzurro). \
            Elica: la punta di E lungo l'asse del fascio; luce circolare = una vite che gira mentre vola, lineare = un'onda piatta. \
            Tratti: E (colori) e B (blu) dall'asse. Poynting: dove va l'energia. Spin: E(t−1)×E(t), verde lungo il volo. \
            Il disegno non tocca la fisica: tutto è letto dal campo delle palline della regola.").small().color(GREY));
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut lab.grandezza, 0.3..=3.0).text("palline"));
            ui.add(egui::Slider::new(&mut lab.luce, 0.2..=4.0).text("luce"));
            ui.checkbox(&mut lab.bagliore, "bagliore");
        });
        pannello_misura(ui, lab);
        ui.label(RichText::new("DEBITI: l'ampiezza è quella di un'onda classica normalizzata (la regola è lineare: un fotone e un \
            impulso forte hanno la stessa forma; il \"clic\" di un solo fotone non c'è ancora). Il reticolo cubico rallenta un po' \
            l'apertura di lato del fascio (fattore ka/tan ka: −21 % con 8 celle per λ, −6 % con 16, −2 % con 24) e rallenta la luce \
            obliqua (0,88 c in diagonale con 8 celle per λ). Lungo un asse va a c esatto.").small());
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
        pannello_misura(ui, lab);
        ui.label(RichText::new("DEBITI: il campo del fotone è un impulso forte (come un laser), non un fotone solo; l'emissione \
            spontanea (vuoto quantistico) non c'è: la sovrapposizione non decade da sola. Il nucleo di Nd ha B = 144: si vede un protone.").small());
    });
    if !open { lab.panel_open = false; }
}

fn editor_stato(ui: &mut egui::Ui, id: &str, s: &mut crate::orbitali::Stato) {
    ui.horizontal(|ui| {
        ui.label("n");
        for n in 1..=4u32 { if ui.selectable_label(s.n == n, n.to_string()).clicked() { s.n = n; } }
        ui.separator();
        for (l, nome) in ["s", "p", "d", "f"].iter().enumerate() {
            let l = l as u32;
            if l < s.n && ui.selectable_label(s.l == l, *nome).clicked() { s.l = l; }
        }
        ui.separator();
        egui::ComboBox::from_id_salt(format!("{id}j")).selected_text(format!("j = {}/2", s.j2)).show_ui(ui, |ui| {
            ui.selectable_value(&mut s.j2, 2 * s.l + 1, format!("{}/2", 2 * s.l + 1));
            if s.l > 0 { ui.selectable_value(&mut s.j2, 2 * s.l - 1, format!("{}/2", 2 * s.l - 1)); }
        });
        egui::ComboBox::from_id_salt(format!("{id}m")).selected_text(format!("m = {:+}/2", s.m2)).show_ui(ui, |ui| {
            let j2 = s.j2 as i32;
            let mut m = -j2;
            while m <= j2 { ui.selectable_value(&mut s.m2, m, format!("{:+}/2", m)); m += 2; }
        });
    });
    *s = s.aggiusta();
}

fn pannello_orbitali(ctx: &egui::Context, lab: &mut Lab) {
    use crate::orbitali::{Scelta, Stato};
    let mut open = true;
    let o = lab.orb.as_ref().map(|z| z.lock().unwrap().out.clone());
    egui::Window::new("Orbitali di un elettrone (regola zig-zag)").open(&mut open).default_pos([12.0, 80.0]).default_width(480.0).vscroll(true).show(ctx, |ui| {
        ui.label(RichText::new("Un solo elettrone attorno al nucleo di neodimio (Z = 60): gli stati sono quelli che dà la regola zig-zag \
            nella sua forma radiale (stessa cella, stessa fase di Coulomb), messi in 3-D con le funzioni d'angolo di Dirac. \
            Una sovrapposizione evolve con le frequenze della regola: la nuvola batte.").color(GOLD));
        ui.label(RichText::new("Da lontano: gusci morbidi a tre livelli di densità (colore = fase dell'onda: gira nel tempo e attorno a z \
            come e^{imφ}). Da vicino: palline-surfista portate dalla corrente dell'onda (v = j/ρ ≤ c), sei con la scia.").small().color(GREY));
        ui.separator();
        ui.label(RichText::new("SCELTE RAPIDE").strong().color(GOLD));
        ui.horizontal_wrapped(|ui| {
            for t in ["1s", "2s", "2p", "3d", "4f", "1s+2p", "2p+3d", "3d+4f", "2s+2p"] {
                if ui.button(t).clicked() {
                    if let Some(s) = crate::orbitali::scelta_da(t) { lab.orb_scelta = s.clone(); if let Some(z) = &lab.orb { z.lock().unwrap().scegli = Some(s); } }
                }
            }
        });
        ui.separator();
        ui.label(RichText::new("STATO SU MISURA").strong().color(GOLD));
        let sc = &mut lab.orb_scelta;
        if sc.stati.is_empty() { sc.stati.push((Stato::new(1, 0, 1), 1.0)); }
        editor_stato(ui, "a", &mut sc.stati[0].0);
        let mut due = sc.stati.len() > 1;
        ui.checkbox(&mut due, "sovrapposizione con un secondo stato");
        if due && sc.stati.len() < 2 { sc.stati.push((Stato::new(2, 1, 1), 1.0)); }
        if !due { sc.stati.truncate(1); }
        if sc.stati.len() > 1 {
            editor_stato(ui, "b", &mut sc.stati[1].0);
            let mut w = sc.stati[1].1 * sc.stati[1].1 / (sc.stati[0].1 * sc.stati[0].1 + sc.stati[1].1 * sc.stati[1].1);
            ui.add(egui::Slider::new(&mut w, 0.05..=0.95).text("peso del secondo"));
            sc.stati[0].1 = (1.0 - w).sqrt();
            sc.stati[1].1 = w.sqrt();
        }
        if ui.add(egui::Button::new(RichText::new("  Mostra  ").strong()).fill(Color32::from_rgb(40, 60, 110))).clicked() {
            let s: Scelta = lab.orb_scelta.clone();
            if let Some(z) = &lab.orb { z.lock().unwrap().scegli = Some(s); }
        }
        ui.separator();
        if let Some(o) = &o {
            ui.label(RichText::new("LIVELLI (dalla regola radiale, cella 0,5 ħ/mc)").strong().color(GOLD));
            egui::Grid::new("liv").striped(true).show(ui, |ui| {
                ui.label("stato"); ui.label("legame, regola"); ui.label("Dirac esatto"); ui.end_row();
                for (s, w, d) in &o.livelli {
                    ui.label(s.nome());
                    ui.label(format!("{:.2} keV", (1.0 - w) * zz::MC2_KEV));
                    ui.label(format!("{:.2} keV ({:+.0} %)", (1.0 - d) * zz::MC2_KEV, 100.0 * ((1.0 - w) / (1.0 - d) - 1.0)));
                    ui.end_row();
                }
            });
            ui.label(RichText::new(format!("raggio mostrato {:.0} ħ/mc = {} · {:.0} ms per passo", o.r, fmt_len(o.r * zz::LAMBDA_C_M), o.ms)).small().color(GREY));
        }
        pannello_misura(ui, lab);
        ui.label(RichText::new("DEBITI: le celle grosse (0,5 ħ/mc) legano troppo (1s +40 % rispetto a Dirac; gli stati alti meno). \
            Gli orbitali n = 3, 4 non stanno nella scatola dal vivo (48³): la forma viene dalla regola radiale, esatta per un \
            elettrone in un campo centrale; qui niente fotone e niente tremolio di Zitterbewegung (stati di sola energia positiva). \
            Gli atomi con più elettroni aspettano la regola a molti elettroni degli agglomerati: non li fingiamo con orbitali \
            di tabella.").small());
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut lab.grandezza, 0.3..=3.0).text("palline"));
            ui.add(egui::Slider::new(&mut lab.luce, 0.2..=4.0).text("luce"));
            ui.checkbox(&mut lab.mostra_protone, "protone al centro");
        });
    });
    if !open { lab.panel_open = false; }
}

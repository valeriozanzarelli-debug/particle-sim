//! Screen of the collision: a thin time bar at the bottom (phase, wavelength, probability per photon,
//! slow-motion factors) and a hidden panel (Tab or ☰) with the computed numbers vs NIST and the counters.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::Color32;

use crate::config::*;
use crate::costanti as k;
use crate::fisica::Fase;
use crate::Lab;

const CYAN: Color32 = Color32::from_rgb(90, 220, 255);
const GOLD: Color32 = Color32::from_rgb(255, 205, 80);
const GREY: Color32 = Color32::from_rgb(150, 155, 170);
const ORANGE: Color32 = Color32::from_rgb(255, 160, 70);

/// "1.23·10⁻⁴" style.
fn sci(x: f64) -> String {
    if x == 0.0 { return "0".into(); }
    let e = x.abs().log10().floor() as i32;
    if (-2..3).contains(&e) { return format!("{x:.3}"); }
    let m = x / 10f64.powi(e);
    const SUP: [char; 10] = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
    let mut s = String::new();
    if e < 0 { s.push('⁻'); }
    for c in e.abs().to_string().chars() { s.push(SUP[c.to_digit(10).unwrap() as usize]); }
    format!("{m:.2}·10{s}")
}

pub fn draw_hud(mut ctx: EguiContexts, mut lab: ResMut<Lab>, keys: Res<ButtonInput<KeyCode>>) {
    let lab = &mut *lab;
    let ctx = ctx.ctx_mut();
    if keys.just_pressed(KeyCode::Tab) && !ctx.wants_keyboard_input() { lab.panel_open = !lab.panel_open; }
    if keys.just_pressed(KeyCode::Space) && !ctx.wants_keyboard_input() { lab.sc.pausa = !lab.sc.pausa; }
    egui::Area::new(egui::Id::new("menu_button"))
        .anchor(egui::Align2::LEFT_TOP, [10.0, 10.0])
        .show(ctx, |ui| {
            let b = egui::Button::new(egui::RichText::new("☰").size(16.0)).fill(Color32::from_black_alpha(90)).selected(lab.panel_open);
            if ui.add(b).on_hover_text("pannello (Tab)").clicked() { lab.panel_open = !lab.panel_open; }
        });
    barra(ctx, lab);
    if lab.panel_open {
        let mut open = true;
        egui::Window::new("Fotone contro atomo di elio")
            .open(&mut open)
            .default_width(430.0)
            .default_pos([12.0, 48.0])
            .show(ctx, |ui| { egui::ScrollArea::vertical().max_height(760.0).show(ui, |ui| panel(ui, lab)); });
        if !open { lab.panel_open = false; }
    }
}

/// The time bar: what is happening, how much it is slowed down, the wavelength and P per photon.
fn barra(ctx: &egui::Context, lab: &mut Lab) {
    let sc = &mut lab.sc;
    egui::Area::new(egui::Id::new("barra"))
        .anchor(egui::Align2::CENTER_BOTTOM, [0.0, -10.0])
        .show(ctx, |ui| {
            egui::Frame::none().fill(Color32::from_rgba_unmultiplied(10, 12, 20, 185)).rounding(6.0)
                .inner_margin(egui::Margin::symmetric(12.0, 6.0)).show(ui, |ui| {
                    ui.set_width(900.0);
                    let (testo, frac) = match sc.fase {
                        Fase::Attesa => ("atomo nello stato fondamentale 1s² ¹S".to_string(), sc.t_fase / PAUSA),
                        Fase::Volo => {
                            let f = sc.fin.as_ref().unwrap();
                            let t_tot = (Z_END - Z_START) / C_SCR;
                            let ta = sc.t - f.t0;
                            let s = if sc.b2 > 0.0 && sc.b2 < 1.0 {
                                format!("ASSORBIMENTO: |b|² = {:.2} — la nuvola oscilla alla frequenza della riga (rallentata {}×)",
                                        sc.b2, sci(sc.rallentamento_osc()))
                            } else if f.preso && !f.risonante && f.fatto {
                                "diffusione di Rayleigh: il fotone cambia direzione, l'atomo resta in 1s²".to_string()
                            } else {
                                "fotone in arrivo (luce rallentata, fotone non in scala)".to_string()
                            };
                            (s, (ta / t_tot).min(1.0))
                        }
                        Fase::Eccitato => (format!("atomo ECCITATO 1s2p ¹P — emetterà dopo {:.3} ns (estratto dalla vita media {:.3} ns; attesa rallentata {}×)",
                                                   sc.attesa_fisica * 1e9, sc.luce.vita_s() * 1e9, sci(S_DECAY)), sc.t_fase / sc.attesa.max(1e-9)),
                        Fase::Emissione => (format!("EMISSIONE spontanea: |b|² = {:.2}, nuovo fotone in direzione da sin²θ (oscillazione rallentata {}×)",
                                                    sc.b2, sci(sc.rallentamento_osc())), (sc.b2 - 1.0).abs()),
                    };
                    ui.label(egui::RichText::new(testo).color(CYAN));
                    ui.add(egui::ProgressBar::new(frac.clamp(0.0, 1.0) as f32).desired_height(5.0));
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new(format!("λ = {:.4} nm   ·   P(per fotone) = {}   ·   inviati {}  assorbiti {}  passati {}  emessi {}",
                            sc.lambda_nm, sci(sc.prob()), sc.cont.inviati, sc.cont.assorbiti, sc.cont.passati, sc.cont.emessi))
                            .monospace().small().color(GOLD));
                    });
                    ui.horizontal(|ui| {
                        if ui.small_button(if sc.pausa { "▶ riprendi" } else { "⏸ pausa" }).clicked() { sc.pausa = !sc.pausa; }
                        if ui.small_button("risonanza (λ calcolata)").clicked() { sc.lambda_nm = sc.luce.lambda_nm; }
                        if ui.small_button("58,4334 nm (NIST)").clicked() { sc.lambda_nm = 58.4334; }
                        if ui.small_button("500 nm").clicked() { sc.lambda_nm = 500.0; }
                        let mut l = sc.lambda_nm;
                        // only a real user action changes λ (the slider rounds what it shows)
                        let r = ui.add(egui::Slider::new(&mut l, 30.0..=800.0).logarithmic(true).text("λ nm").max_decimals(4));
                        if r.changed() && (r.dragged() || r.has_focus()) { sc.lambda_nm = l; }
                        // detuning from the computed line, in units of the natural width Γ = A
                        let w0 = sc.luce.omega0_s();
                        let mut det = (crate::elio::LuceElio::omega_s(sc.lambda_nm) - w0) / sc.luce.a_s;
                        if det.abs() <= 30.0 || det.is_nan() {
                            let r = ui.add(egui::Slider::new(&mut det, -30.0..=30.0).text("Δω/Γ"));
                            if r.changed() && (r.dragged() || r.has_focus()) {
                                let w = w0 + det * sc.luce.a_s;
                                sc.lambda_nm = 2.0 * std::f64::consts::PI * k::C_L / w * 1e9;
                            }
                        } else {
                            ui.label(egui::RichText::new(format!("Δω/Γ = {}", sci(det))).small().color(GREY));
                        }
                    });
                });
        });
}

fn riga(ui: &mut egui::Ui, nome: &str, calc: String, rif: &str, diff: String) {
    ui.label(nome);
    ui.label(egui::RichText::new(calc).color(CYAN));
    ui.label(egui::RichText::new(rif).color(GOLD));
    ui.label(egui::RichText::new(diff).color(GREY));
    ui.end_row();
}

fn panel(ui: &mut egui::Ui, lab: &mut Lab) {
    ui.horizontal(|ui| {
        ui.checkbox(&mut lab.mostra_e[0], egui::RichText::new("elettrone spin su").color(CYAN));
        ui.checkbox(&mut lab.mostra_e[1], egui::RichText::new("elettrone spin giù").color(ORANGE));
        ui.checkbox(&mut lab.mostra_surfer, "surfer");
    });
    let sc = &lab.sc;
    let r = &lab.res;
    let l = &sc.luce;
    let f = |key: &str| r.f(key).unwrap_or(f64::NAN);
    let ev = k::hartree_ev();
    ui.separator();
    ui.label(egui::RichText::new("Calcolato dalle palline  vs  riferimento").strong());
    egui::Grid::new("tab").striped(true).num_columns(4).show(ui, |ui| {
        ui.label("");
        ui.label(egui::RichText::new("palline").color(CYAN));
        ui.label(egui::RichText::new("riferimento").color(GOLD));
        ui.label(egui::RichText::new("differenza").color(GREY));
        ui.end_row();
        riga(ui, "E(1s² ¹S), Eₕ", format!("{:.5} ± {:.5}", l.e_g, f("he1p_Eg_comb_err")), "−2,90372 (esatta)",
             format!("{:+.5}", l.e_g + 2.903_724));
        riga(ui, "E(1s2p ¹P), Eₕ", format!("{:.5} ± {:.5}", l.e_e, f("he1p_E_err")), "−2,12384 (esatta)",
             format!("{:+.5} (nodo fisso)", l.e_e + 2.123_843));
        riga(ui, "ΔE della riga", format!("{:.4} ± {:.4} eV", l.de * ev, f("he1p_de_err") * ev), "21,2180 eV (NIST)",
             format!("{:+.1} meV", (l.de * ev - 21.2180) * 1e3));
        riga(ui, "λ della riga", format!("{:.3} ± {:.3} nm", l.lambda_nm, f("he1p_lambda_err")), "58,4334 nm (NIST)",
             format!("{:+.3} %", (l.lambda_nm / 58.4334 - 1.0) * 100.0));
        riga(ui, "⟨1¹S|z₁+z₂|2¹P⟩", format!("{:.4} ± {:.4} a₀", f("he1p_d"), f("he1p_d_err")), "0,4208 a₀ (da f)",
             format!("{:+.1} %", (f("he1p_d") / 0.4208 - 1.0) * 100.0));
        riga(ui, "forza d'oscillatore f", format!("{:.4} ± {:.4}", l.f, f("he1p_f_err")), "0,2762 (NIST)",
             format!("{:+.1} %", (l.f / 0.2762 - 1.0) * 100.0));
        riga(ui, "A (2¹P→1¹S)", format!("{} s⁻¹", sci(l.a_s)), "1,799·10⁹ s⁻¹ (NIST)", format!("{:+.1} %", (l.a_s / 1.799e9 - 1.0) * 100.0));
        riga(ui, "vita media", format!("{:.3} ns", l.vita_s() * 1e9), "0,555 ns", format!("{:+.1} %", (l.vita_s() * 1e9 / 0.555 - 1.0) * 100.0));
        riga(ui, "polarizzabilità α", format!("{:.3} ± {:.3} a₀³", l.alpha0, f("he_alpha_err")), "1,3832 a₀³ (esatta)",
             format!("{:+.1} %", (l.alpha0 / 1.3832 - 1.0) * 100.0));
    });
    ui.separator();
    let sres = 3.0 * (l.lambda_nm * 1e-9).powi(2) / (2.0 * std::f64::consts::PI);
    ui.label(egui::RichText::new("Probabilità per fotone").strong());
    ui.label(format!("σ(λ = {:.4} nm) = {} m²   →   P = σ / (π w0²/2), w0 = λ:  {}", sc.lambda_nm, sci(l.sigma(sc.lambda_nm)), sci(sc.prob())));
    ui.label(egui::RichText::new(format!(
        "al centro della riga σ = 3λ²/2π = {} m², P = 3/π² = 0,304. A 500 nm: P = {} (diffusione di Rayleigh con la α calcolata: \
         il fotone passa dritto quasi sempre). Alla λ del NIST (58,4334 nm) il nostro atomo assorbe con P = {}: la nostra riga è \
         spostata di {:.0} larghezze naturali Γ — lo scarto di {:.1} meV è 10⁴ volte ħΓ.",
        sci(sres), sci(l.prob(500.0)), sci(l.prob(58.4334)),
        (crate::elio::LuceElio::omega_s(58.4334) - l.omega0_s()).abs() / l.a_s, (l.de * ev - 21.2180).abs() * 1e3)).small().color(GREY));
    ui.separator();
    ui.label(egui::RichText::new("Contatori").strong());
    ui.label(format!("fotoni inviati {}   ·   assorbiti {}   ·   passati dritti {}   ·   emessi {}   ·   diffusi (Rayleigh) {}",
        sc.cont.inviati, sc.cont.assorbiti, sc.cont.passati, sc.cont.emessi, sc.cont.diffusi));
    ui.separator();
    ui.label(egui::RichText::new("Tempi e scale (niente è in tempo reale)").strong());
    ui.label(egui::RichText::new(format!(
        "• oscillazione della nuvola: ΔE/ħ = {} rad/s, mostrata a {:.1} rad/s: rallentata {}×. Il fotone è disegnato con λ = {:.0} unità \
         (vera: {:.1} nm = {:.0} a₀, cioè {:.0} volte l'atomo) e la sua onda passa sull'atomo con la stessa frequenza della nuvola.\n\
         • attesa nello stato eccitato rallentata {}× (vita {:.3} ns → {:.1} s in media).\n\
         • un vero assorbimento dura ~vita media ≈ {} oscillazioni; qui ne vedi ~3. Il pacchetto vero è lungo ~c/Γ ≈ {:.0} cm.",
        sci(l.omega0_s()), sc.omega_scr(), sci(sc.rallentamento_osc()), crate::fisica::lambda_draw(sc.lambda_nm),
        sc.lambda_nm, sc.lambda_nm * 1e-9 / k::bohr_m(), sc.lambda_nm * 1e-9 / k::bohr_m() / 5.0,
        sci(S_DECAY), l.vita_s() * 1e9, l.vita_s() * S_DECAY,
        sci(l.vita_s() * l.omega0_s() / std::f64::consts::TAU), k::C_L / l.a_s * 100.0)).small().color(GREY));
    ui.separator();
    ui.label(egui::RichText::new("Che cosa stai guardando (onestamente)").strong());
    ui.label(egui::RichText::new(
        "• Le palline azzurre (spin su) e arancio (spin giù) sono le due nuvole calcolate dal vivo con la regola delle palline: \
         1s² ¹S e 1s2p ¹P (nodo fisso dai nostri orbitali 1s e 2p, k_b dal minimo della nostra energia variazionale).\n\
         • Durante assorbimento ed emissione la nuvola è |a ψ_g + b e^(−iΔE t) ψ_e|²: palline prese dalle due nuvole e tenute con \
         il peso dell'interferenza, che usa le nostre guide per il rapporto ψ_e/ψ_g. La freccia gialla sul nucleo è il          dipolo istantaneo ⟨z₁+z₂⟩ della nuvola disegnata (×4): va su e giù alla frequenza della riga. È la funzione d'onda CONDIZIONATA dell'atomo \
         (dato dove sta il fotone). Con un solo fotone, la media su tutte le possibilità (matrice densità ridotta) NON oscilla: \
         l'oscillazione del dipolo si vede solo così, o con luce classica.\n\
         • |b|² sale mentre il pacchetto passa (dato che è stato assorbito: probabilità P) — la rampa è disegnata, non calcolata \
         risolvendo l'equazione di Schrödinger con il campo.\n\
         • Il surfer (due stelle bianche: UNA configurazione dei due elettroni) è fermo negli stati stazionari (onda reale: \
         velocità di Bohm zero) e si muove con la velocità di Bohm della sovrapposizione disegnata (calcolata con le nostre guide). \
         Va avanti e indietro alla frequenza della riga; il suo spostamento netto non è esatto, perché la rampa di |b|² è compressa.\n\
         • Il fotone emesso è un'onda di dipolo (densità ∝ sin²θ attorno alla polarizzazione); la sua pallina-fotone parte in una \
         direzione estratta da sin²θ.").small().color(GREY));
}

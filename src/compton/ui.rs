//! Panel of the Compton viewer (hidden by default: Tab or ☰). Italian, for a 20-year-old.

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::{Color32, RichText};

use crate::config::*;
use crate::fisica::*;
use crate::scena::{statistica, N_BIN_MU, N_BIN_T};
use crate::{Lab, Vista};

const CYAN: Color32 = Color32::from_rgb(90, 220, 255);
const GOLD: Color32 = Color32::from_rgb(255, 205, 80);
const GREY: Color32 = Color32::from_rgb(150, 155, 170);
const PINK: Color32 = Color32::from_rgb(255, 120, 230);

const PRESET: [(&str, f64); 6] = [
    ("luce visibile 2,48 eV (500 nm)", 2.48e-3),
    ("Mo Kα 17,48 keV (Compton 1923)", 17.48),
    ("100 keV", 100.0),
    ("γ Cs-137 662 keV", 662.0),
    ("γ Co-60 1,17 MeV", 1173.2),
    ("γ Co-60 1,33 MeV", 1332.5),
];

fn fmt_e(e_kev: f64) -> String {
    if e_kev < 1.0 { format!("{:.3} eV", e_kev * 1e3) } else if e_kev < 1000.0 { format!("{e_kev:.2} keV") } else { format!("{:.4} MeV", e_kev / 1e3) }
}

fn fmt_len(pm: f64) -> String {
    if pm >= 1e5 { format!("{:.1} nm", pm / 1e3) } else { format!("{pm:.4} pm") }
}

pub fn draw_hud(mut ctx: EguiContexts, mut lab: ResMut<Lab>, keys: Res<ButtonInput<KeyCode>>) {
    let lab = &mut *lab;
    let ctx = ctx.ctx_mut();
    if !ctx.wants_keyboard_input() {
        if keys.just_pressed(KeyCode::Tab) { lab.panel_open = !lab.panel_open; }
        if keys.just_pressed(KeyCode::Space) { lab.paused = !lab.paused; }
        if keys.just_pressed(KeyCode::KeyN) { lab.nuovo_colpo(); }
    }
    egui::Area::new(egui::Id::new("menu_button"))
        .anchor(egui::Align2::LEFT_TOP, [10.0, 10.0])
        .show(ctx, |ui| {
            let b = egui::Button::new(RichText::new("☰").size(16.0)).fill(Color32::from_black_alpha(90)).selected(lab.panel_open);
            if ui.add(b).on_hover_text("pannello (Tab) · pausa (Spazio) · nuovo colpo (N)").clicked() {
                lab.panel_open = !lab.panel_open;
            }
        });
    egui::Area::new(egui::Id::new("status"))
        .anchor(egui::Align2::CENTER_BOTTOM, [0.0, -10.0])
        .show(ctx, |ui| {
            egui::Frame::none().fill(Color32::from_rgba_unmultiplied(10, 12, 20, 170)).rounding(6.0)
                .inner_margin(egui::Margin::symmetric(10.0, 5.0)).show(ui, |ui| {
                    ui.label(RichText::new(stato(lab)).monospace().color(CYAN));
                });
        });
    if lab.panel_open {
        let mut open = true;
        egui::Window::new("Compton: un fotone contro un elettrone libero")
            .open(&mut open)
            .default_width(430.0)
            .default_height(820.0)
            .default_pos([12.0, 48.0])
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().max_height(780.0).show(ui, |ui| panel(ui, lab));
            });
        if !open { lab.panel_open = false; }
    }
}

fn stato(lab: &Lab) -> String {
    let c = &lab.colpo;
    let head = format!("{} · λ {}", fmt_e(lab.e_kev), fmt_len(lab.kn.lam_pm));
    let t = lab.clock;
    let what = if t < c.t_c {
        "prima dell'urto: il fotone arriva, l'elettrone è fermo".to_string()
    } else {
        match c.evento() {
            Some(e) => format!("fotone θ {:.1}°, λ' {} · elettrone {:.1}°, T {}",
                               e.theta.to_degrees(), fmt_len(e.lam_out), e.theta_e.to_degrees(), fmt_e(e.t_e)),
            None => "il fotone è passato dritto, l'elettrone resta fermo".to_string(),
        }
    };
    format!("{head} · colpo {} · {what}{}", c.id, if lab.paused { " · PAUSA" } else { "" })
}

fn panel(ui: &mut egui::Ui, lab: &mut Lab) {
    // ── energy ──
    ui.label(RichText::new("Energia del fotone").strong());
    ui.horizontal_wrapped(|ui| {
        for (nome, e) in PRESET {
            if ui.selectable_label((lab.e_kev - e).abs() < 1e-9 * e, nome).clicked() { lab.cambia_energia(e); }
        }
    });
    let mut lg = lab.e_kev.log10();
    if ui.add(egui::Slider::new(&mut lg, -3.0..=3.5).text("log₁₀ E/keV").show_value(false)).changed() {
        lab.cambia_energia(10f64.powf(lg));
    }
    ui.horizontal(|ui| {
        ui.label("polarizzazione:");
        let mut p = lab.pol;
        ui.radio_value(&mut p, Pol::NonPolarizzata, "non polarizzata");
        ui.radio_value(&mut p, Pol::LineareX, "lineare (freccia)");
        if p != lab.pol { lab.pol = p; lab.stat = None; lab.nuovo_colpo(); }
    });
    ui.checkbox(&mut lab.solo_diffusi, "mostra solo i colpi in cui il fotone devia (condizionato)")
        .on_hover_text("Se spento, i surfisti nascono con la vera |ψ|²: il fotone quasi sempre passa dritto (p ≈ 10⁻¹¹).");
    ui.horizontal(|ui| {
        if ui.button(if lab.paused { "▶ riprendi" } else { "⏸ pausa" }).clicked() { lab.paused = !lab.paused; }
        if ui.button("nuovo colpo (N)").clicked() { lab.nuovo_colpo(); }
        ui.checkbox(&mut lab.auto, "automatico");
    });
    ui.add(egui::Slider::new(&mut lab.speed, 0.05..=3.0).text("velocità").logarithmic(true));
    ui.horizontal(|ui| {
        ui.checkbox(&mut lab.mostra_surfisti, "surfisti");
        ui.checkbox(&mut lab.mostra_coppie, "coppie correlate");
    });
    ui.horizontal(|ui| {
        ui.label("vista:");
        ui.radio_value(&mut lab.vista, Vista::Sezione, "sezione nel piano dell'urto");
        ui.radio_value(&mut lab.vista, Vista::Meta, "3-D a metà");
        ui.radio_value(&mut lab.vista, Vista::Intera, "3-D intera");
    });
    ui.add(egui::Slider::new(&mut lab.gain_onda, 0.2..=4.0).text("luminosità onda diffusa (solo visiva)").logarithmic(true));
    ui.add(egui::Slider::new(&mut lab.ball_scale, 0.4..=3.0).text("dimensione palline"));
    ui.separator();
    evento(ui, lab);
    ui.separator();
    grafici(ui, lab);
    ui.separator();
    numeri(ui, lab);
    ui.separator();
    storia(ui);
    ui.separator();
    spiegazione(ui);
}

fn riga(ui: &mut egui::Ui, a: &str, b: String) {
    ui.label(a);
    ui.label(RichText::new(b).color(CYAN));
    ui.end_row();
}

fn numeri(ui: &mut egui::Ui, lab: &Lab) {
    let kn = &lab.kn;
    let lam_m = kn.lam_pm * 1e-12;
    let slow = C_M_S / (C_SCENE as f64 * lam_m);
    let p = lab.pk.p_scatter(kn);
    let t_show = (R_END - Z_START) as f64 / C_SCENE as f64 / slow;
    egui::Grid::new("num").striped(true).num_columns(2).show(ui, |ui| {
        riga(ui, "λ entrante", fmt_len(kn.lam_pm));
        riga(ui, "λ' a 90° / a 180°", format!("{} / {}", fmt_len(lam_out(kn.lam_pm, PI_2)), fmt_len(lam_out(kn.lam_pm, std::f64::consts::PI))));
        riga(ui, "spostamento relativo Δλ/λ max", format!("{:.3e}", 2.0 * LAMBDA_C_PM / kn.lam_pm));
        riga(ui, "k = E/mec²", format!("{:.4e}", kn.k));
        riga(ui, "σ Klein–Nishina", format!("{:.4e} m²  = {:.4} σ_T", kn.sigma_m2(), kn.sigma_m2() / sigma_thomson()));
        riga(ui, "raggio del disco √(σ/π)", format!("{:.3e} m = {:.2e} λ", (kn.sigma_m2() / std::f64::consts::PI).sqrt(),
                                                     (kn.sigma_m2() / std::f64::consts::PI).sqrt() / lam_m));
        riga(ui, "probabilità di deviare, per fotone", format!("{p:.3e}"));
        riga(ui, "bordo Compton (T max dell'elettrone)", fmt_e(compton_edge(lab.e_kev)));
        riga(ui, "rallentamento", format!("{slow:.2e} volte (c = {C_SCENE} λ al secondo)"));
        riga(ui, "allargamento del pacchetto-elettrone", format!("τ = {:.2e} s; durata della scena {:.2e} s", tau_spread(S_E_M), t_show));
        riga(ui, "INPUT non dalle palline", format!("Klein–Nishina (QED); largh. fotone {S_GAMMA_LAMBDA} λ, elettrone {} nm", S_E_M * 1e9));
    });
    ui.label(RichText::new(format!("Disegno non in scala: l'elettrone è disegnato largo {S_E_DRAW} λ; il guscio diffuso e la nuvola \
        dell'elettrone che rincula sono disegnati circa {:.0e} volte più luminosi del loro peso vero ({p:.1e}).", 1.0 / p)).small().color(GREY));
}

const PI_2: f64 = std::f64::consts::FRAC_PI_2;

fn evento(ui: &mut egui::Ui, lab: &Lab) {
    let c = &lab.colpo;
    ui.label(RichText::new(format!("Colpo {}  (deviati {} su {} mostrati)", c.id, lab.n_diffusi, lab.n_colpi)).strong());
    let b = c.b_m[0].hypot(c.b_m[1]);
    ui.label(format!("distanza trasversa dei due surfisti quando si incontrano: b = {b:.3e} m"));
    match c.evento() {
        None => { ui.label(RichText::new("b fuori dal disco σ: il fotone passa dritto, l'elettrone resta fermo.").color(GREY)); }
        Some(e) => {
            let p_in = lab.e_kev;
            let r = [e.pg[0] + e.pe[0], e.pg[1] + e.pe[1], e.pg[2] + e.pe[2] - p_in];
            egui::Grid::new("ev").striped(true).num_columns(2).show(ui, |ui| {
                riga(ui, "fotone: θ, φ", format!("{:.2}°, {:.1}°", e.theta.to_degrees(), e.phi.to_degrees()));
                riga(ui, "λ' (formula di Compton)", format!("{}  (Δλ = {})", fmt_len(e.lam_out), fmt_len(e.lam_out - lab.kn.lam_pm)));
                riga(ui, "E' del fotone", fmt_e(e.e_out));
                riga(ui, "elettrone: angolo, T, v/c", format!("{:.2}°, {}, {:.4}", e.theta_e.to_degrees(), fmt_e(e.t_e), e.beta_e));
                riga(ui, "tan φe = cot(θ/2)/(1+k)", format!("{:.2}°", electron_angle_formula(lab.e_kev, e.theta).to_degrees()));
                riga(ui, "|p'γ + pe − p| / p", format!("{:.1e}", (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt() / p_in));
            });
        }
    }
}

fn grafici(ui: &mut egui::Ui, lab: &mut Lab) {
    ui.horizontal(|ui| {
        if ui.button("lancia 20 000 colpi (regola dei surfisti)").clicked() {
            lab.stat = Some(statistica(lab.e_kev, lab.pol, 20_000, 2026));
        }
        ui.label(RichText::new("istogrammi contro Klein–Nishina").small().color(GREY));
    });
    let Some(s) = &lab.stat else { return };
    ui.label(format!("cos θ del fotone deviato (40 classi): χ² = {:.1} su {} gradi di libertà", s.chi2_mu, N_BIN_MU - 1));
    plot(ui, &s.h_mu, &s.kn_mu, None, "θ = 180°", "θ = 0°");
    let edge = compton_edge(s.e_kev);
    ui.label(format!("energia cinetica dell'elettrone (50 classi): χ² = {:.1};  bordo Compton {}  (massimo visto {})",
                     s.chi2_t, fmt_e(edge), fmt_e(s.t_max_seen)));
    plot(ui, &s.h_t, &s.kn_t, Some(edge / (1.05 * edge)), "0", &fmt_e(1.05 * edge));
    ui.label(RichText::new(format!("conservazione dell'impulso evento per evento: max |Δp|/p = {:.1e}", s.max_resid)).small());
    let _ = N_BIN_T;
}

/// Bars = surfer events, line = Klein–Nishina expectation; optional vertical marker (fraction of width).
fn plot(ui: &mut egui::Ui, h: &[f64], e: &[f64], mark: Option<f64>, left: &str, right: &str) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width().min(410.0), 110.0), egui::Sense::hover());
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 3.0, Color32::from_rgb(14, 17, 28));
    let ymax = h.iter().chain(e.iter()).cloned().fold(1.0, f64::max) * 1.08;
    let n = h.len();
    let w = rect.width() / n as f32;
    let y = |v: f64| rect.bottom() - (v / ymax) as f32 * (rect.height() - 14.0);
    for (i, v) in h.iter().enumerate() {
        let x0 = rect.left() + i as f32 * w;
        p.rect_filled(egui::Rect::from_min_max(egui::pos2(x0 + 0.5, y(*v)), egui::pos2(x0 + w - 0.5, rect.bottom())), 0.0,
                      Color32::from_rgb(70, 120, 170));
    }
    let pts: Vec<egui::Pos2> = e.iter().enumerate().map(|(i, v)| egui::pos2(rect.left() + (i as f32 + 0.5) * w, y(*v))).collect();
    p.add(egui::Shape::line(pts, egui::Stroke::new(1.6, GOLD)));
    if let Some(m) = mark {
        let x = rect.left() + m as f32 * rect.width();
        p.line_segment([egui::pos2(x, rect.top()), egui::pos2(x, rect.bottom())], egui::Stroke::new(1.2, PINK));
    }
    let f = egui::FontId::proportional(10.0);
    p.text(rect.left_top() + egui::vec2(3.0, 2.0), egui::Align2::LEFT_TOP, left, f.clone(), GREY);
    p.text(rect.right_top() + egui::vec2(-3.0, 2.0), egui::Align2::RIGHT_TOP, right, f, GREY);
}

fn storia(ui: &mut egui::Ui) {
    ui.label(RichText::new("Compton 1923 (Phys. Rev. 21, 483): raggi X Mo Kα su grafite").strong());
    egui::Grid::new("storia").striped(true).num_columns(4).show(ui, |ui| {
        for h in ["angolo", "Δλ misurato", "Δλ nostro", "diff."] { ui.label(RichText::new(h).weak()); }
        ui.end_row();
        for &(ang, lam, dl) in COMPTON_1923.iter() {
            let ours = lam_out(lam, ang.to_radians()) - lam;
            ui.label(format!("{ang:.0}°"));
            ui.label(format!("{dl:.1} pm"));
            ui.label(RichText::new(format!("{ours:.3} pm")).color(CYAN));
            ui.label(RichText::new(format!("{:+.0} %", 100.0 * (dl - ours) / ours)).color(GOLD));
            ui.end_row();
        }
    });
    ui.label(RichText::new("Il valore misurato è quello citato dall'articolo (0,708 → 0,730 Å); la precisione dello \
        spettrometro del 1923 era di circa 0,002 Å. Le misure a 45° e 135° (Phys. Rev. 22, 409) non sono riportate qui \
        perché non abbiamo potuto ricontrollarne i numeri.").small().color(GREY));
}

fn spiegazione(ui: &mut egui::Ui) {
    ui.label(RichText::new("Che cosa stai guardando").strong());
    for t in [
        "PRIMA. A sinistra arriva l'onda del fotone: tante palline-onda che volano alla velocità della luce, più luminose \
         sulle creste (distanti λ). Al centro c'è l'onda dell'elettrone fermo (palline azzurre). Dentro ognuna delle due \
         onde c'è un surfista: la pallina bianca (fotone) e quella azzurra (elettrone). Il surfista dell'elettrone è fermo.",
        "L'URTO. L'onda in uscita non è più \"un'onda per il fotone e un'onda per l'elettrone\": è una somma di coppie. \
         Per ogni direzione θ del fotone c'è un pezzo d'onda del fotone che va lì E un pezzo d'onda dell'elettrone che \
         rincula nella direzione coniugata. Energia e impulso conservati fissano tutto: λ' − λ = (h/mec)(1 − cos θ) \
         con h/mec = 2,42631 pm, e tan φe = cot(θ/2)/(1 + E/mec²). Questa parte è derivata, esatta.",
        "QUANTA onda va in ogni direzione NON viene dalle nostre palline: la prendiamo dalla QED (formula di \
         Klein–Nishina, 1929), compresa la dipendenza dalla polarizzazione. È un ingresso, non un risultato.",
        "Solo un pezzettino minuscolo dell'onda viene deviato (σ_T = 6,652·10⁻²⁹ m² a bassa energia, meno ad alta): \
         quasi tutta l'onda del fotone passa dritta, e la vedi andare avanti. Il guscio colorato e la nuvola azzurra \
         dell'elettrone che rincula sono disegnati enormemente più luminosi del vero, se no non li vedresti.",
        "I SURFISTI scelgono UNA coppia. Quando il surfista-fotone arriva all'elettrone, conta la sua distanza \
         trasversa b dal surfista-elettrone: se b cade nel disco di area σ (raggio √(σ/π), pochi femtometri) il fotone \
         entra nella parte deviata; più è piccolo b, più va all'indietro (mantiene il suo \"posto in fila\" nel flusso, \
         come nella regola base al bordo); l'angolo attorno all'asse di b sceglie l'azimut. Altrimenti passa dritto e \
         l'elettrone resta fermo. Siccome i surfisti nascono con la densità |ψ|², gli angoli escono con la statistica \
         di Klein–Nishina (verificato: χ² con 20 000 colpi). La vera mappa di Bohm dell'onda completa a 6 dimensioni \
         sarebbe un'altra mappa con la stessa statistica: le singole traiettorie non sono garantite, le frequenze sì.",
        "DOPO. Fotone ed elettrone sono entangled: le loro direzioni sono legate (linea rosa). Ogni coppia somma \
         esattamente l'impulso del fotone iniziale, evento per evento (non solo in media). Le linee rosa tenui \
         uniscono altre coppie possibili dell'onda: l'onda le contiene tutte, i surfisti ne seguono una.",
        "COLORI e BANDE. Il colore del guscio dice di quanto si è allungata l'onda (giallo chiaro = niente, magenta = \
         2 h/mec all'indietro). Le creste distano λ'(θ): all'indietro le bande sono più larghe. Con la luce visibile \
         l'allungamento è 10⁻⁶ λ (limite di Thomson): non si vede. Con i raggi gamma è più grande di λ stessa.",
    ] {
        ui.label(RichText::new(t).small());
        ui.add_space(3.0);
    }
}

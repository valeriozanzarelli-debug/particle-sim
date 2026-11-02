//! GRAFICA: the drawing settings and the continuous level of detail (Valerio, 13/10: "farlo in 3-D", "allontanandoti lo
//! vedi sempre più come lo stai rappresentando ora, con dettagli diversi", "surfisti chiari con una mini scia").
//! Drawing only: nothing here touches a rule.
//!
//! LEVEL OF DETAIL, one interface for every level. A level of the world (the electron's surfers inside an atom today;
//! later atoms inside a molecule, molecules inside a droplet, quarks inside a nucleus) is a PROVIDER of glows: it says
//! how many glows it has and over what radius (`Gruppo`). `lod` turns that into what the screen can hold: the fraction
//! kept (a stable subset, the same glows stay while zooming), how much bigger each kept glow gets, and a brightness
//! factor, so that the light per pixel stays the same at every distance (no white ball) while the detail drops.
//! It also returns the SCREEN-SPACE ERROR of the level (`Lod::errore_px`, the size of the smallest detail it still
//! resolves, in pixels): a parent level can replace a child by its own coarser glows when that error is below a pixel,
//! and a child can ask to be refined where the error is large or a physics event happens (`Lod::raffina`).

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use egui::{Color32, RichText};

use crate::fotscena::QUALITA;
use crate::Lab;

/// Light glows drawn at most per light box, by preset (bassa … ultra): the drawing's resolution, never the dynamics.
pub const BUDGET_LUCE: [usize; 5] = [30_000, 45_000, 60_000, 90_000, 150_000];

#[derive(Resource, Clone, Copy, Debug)]
pub struct Grafica {
    /// the settings window is open (G or the ⚙ button)
    pub aperto: bool,
    /// preset: index into fotscena::QUALITA (resolution of the light rule) and BUDGET_LUCE
    pub preset: usize,
    /// length of the surfers' glow tail (screen seconds)
    pub scia_s: f32,
    /// strength of the level of detail (0 = every glow always, 1 default, 3 strong)
    pub lod: f32,
    /// brightness of the surfer glows
    pub bagliori: f32,
    /// slow rotation of the camera when nobody touches the mouse (R)
    pub auto_ruota: bool,
    /// depth cue: nearer glows brighter, farther dimmer (0 off … 1)
    pub profondita: f32,
    /// the light around the atom is cut open towards the camera (you see the shells from inside, in 3-D)
    pub taglio: bool,
    /// depth (scene units, ħ/mc in the atom) over which the cut-open light fades behind the atom
    pub taglio_prof: f32,
    /// the old flat slab of the scattered view (2-D), instead of the 3-D volume
    pub fetta: bool,
    /// the old electron cloud magnified ×6 in the scattered view (off: the atom is drawn at its true size)
    pub zoom_atomo: bool,
}

impl Default for Grafica {
    fn default() -> Self {
        Grafica { aperto: false, preset: 0, scia_s: 0.3, lod: 1.0, bagliori: 1.0, auto_ruota: true, profondita: 0.6,
                  taglio: true, taglio_prof: 90.0, fetta: false, zoom_atomo: false }
    }
}

impl Grafica {
    /// From the command line: `--scia S`, `--lod X`, `--fetta`, `--zoom-atomo`, `--senza-taglio`, `--senza-rotazione`,
    /// `--qualita bassa|media|medio-alta|alta|ultra`.
    pub fn da_argomenti(args: &[String]) -> Self {
        let arg = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
        let ha = |k: &str| args.iter().any(|a| a == k);
        let mut g = Grafica::default();
        if let Some(v) = arg("--scia").and_then(|v| v.parse::<f32>().ok()) { g.scia_s = v.clamp(0.0, 2.0); }
        if let Some(v) = arg("--lod").and_then(|v| v.parse::<f32>().ok()) { g.lod = v.clamp(0.0, 5.0); }
        if let Some(q) = arg("--qualita").and_then(|v| crate::fotscena::preset(&v)) { g.preset = q; }
        g.fetta = ha("--fetta");
        g.zoom_atomo = ha("--zoom-atomo");
        g.taglio = !ha("--senza-taglio");
        g.auto_ruota = !ha("--senza-rotazione");
        g
    }
    pub fn budget_luce(&self) -> usize { BUDGET_LUCE[self.preset.min(4)] }
}

/// What the drawing tells a scene worker every frame (physics coordinates).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VistaWorker {
    /// unit vector from the scene's centre towards the eye
    pub occhio: [f64; 3],
    /// camera distance and scene units per pixel at the focus
    pub r_cam: f64,
    pub bpp: f64,
    /// fraction of the surfers kept by the level of detail and how much their glows grow (for the tails)
    pub lod_f: f32,
    pub scia_s: f64,
    pub taglio: bool,
    pub taglio_prof: f64,
    pub fetta: bool,
    pub zoom_atomo: bool,
    pub budget: usize,
}

impl Default for VistaWorker {
    fn default() -> Self {
        VistaWorker { occhio: [0.0, -1.0, 0.0], r_cam: 11.0, bpp: 0.01, lod_f: 1.0, scia_s: 0.3, taglio: true, taglio_prof: 90.0,
                      fetta: false, zoom_atomo: false, budget: BUDGET_LUCE[0] }
    }
}

impl VistaWorker {
    /// Changed enough that a paused worker should draw again.
    pub fn diversa(&self, o: &VistaWorker) -> bool {
        let d = self.occhio[0] * o.occhio[0] + self.occhio[1] * o.occhio[1] + self.occhio[2] * o.occhio[2];
        d < 0.9995 || (self.r_cam / o.r_cam).ln().abs() > 0.05 || (self.lod_f / o.lod_f.max(1e-9)).ln().abs() > 0.05
            || self.scia_s != o.scia_s || self.taglio != o.taglio || self.taglio_prof != o.taglio_prof || self.fetta != o.fetta
            || self.zoom_atomo != o.zoom_atomo || self.budget != o.budget
    }
}

/// A level of the world seen as a provider of glows.
#[derive(Clone, Copy, Debug)]
pub struct Gruppo {
    /// how many glows the level has at full detail
    pub n: usize,
    /// radius (scene units) over which they are spread
    pub raggio: f64,
    /// length (scene units) of each glow's tail, 0 without tails
    pub scia: f64,
}

/// The level of detail of a group at the present view.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lod {
    /// fraction of the glows kept (a stable subset: chiave(i) < f)
    pub f: f32,
    /// brightness factor of each kept glow
    pub alfa: f32,
    /// size factor of each kept glow (fewer, bigger)
    pub cresci: f32,
    /// mean number of glows over a pixel of the group (≤ the target)
    pub sovrapposizione: f64,
    /// screen-space error: size (px) of the smallest detail still drawn = mean spacing of the kept glows
    pub errore_px: f64,
}

impl Lod {
    /// A child level asks to be refined (drawn by its own finer level) when its detail is coarser than this.
    pub fn raffina(&self, soglia_px: f64) -> bool { self.errore_px > soglia_px }
}

/// Glows over a pixel allowed at strength 1, on average over the group's disc (its centre holds more).
pub const SOVRAPPOSIZIONE: f64 = 16.0;
/// Never fewer than this many glows for a group (a far atom stays a small knot of glows, not one dot).
pub const N_MIN: f64 = 150.0;
/// Faintest a kept glow gets when the group is too small on the screen for N_MIN glows.
pub const ALFA_MIN: f64 = 0.6;
/// Largest growth of a kept glow.
pub const CRESCI_MAX: f64 = 2.5;

/// Level of detail of a group whose glows are `s_px` pixels in radius, at `bpp` scene units per pixel.
/// Mean overlap at full detail m₁ = 1.5 n (s² + 0.15 s·L) / R² (glows of radius s with a tapered tail L, on a disc of
/// radius R, 1.5 for the depth of a ball). Over the target K = 16/strength: keep f, grow each glow by g = f^(−1/6) (≤ 2.5),
/// so that m₁ f g² = K; below N_MIN glows, keep N_MIN and dim them instead (not below ALFA_MIN).
pub fn lod(g: Gruppo, bpp: f64, s_px: f64, forza: f64) -> Lod {
    let n = g.n.max(1) as f64;
    let r_px = (g.raggio / bpp.max(1e-30)).max(1.0);
    let l_px = g.scia / bpp.max(1e-30);
    let m1 = 1.5 * n * (s_px * s_px + 0.15 * s_px * l_px) / (r_px * r_px);
    let spaziatura = |f: f64| r_px * (4.0 / (3.0 * n * f)).cbrt();
    if forza <= 0.0 || m1 <= SOVRAPPOSIZIONE / forza.max(1e-9) {
        return Lod { f: 1.0, alfa: 1.0, cresci: 1.0, sovrapposizione: m1, errore_px: spaziatura(1.0) };
    }
    let k = SOVRAPPOSIZIONE / forza;
    let mut f = (k / m1).powf(1.5);
    let mut gr = f.powf(-1.0 / 6.0);
    if gr > CRESCI_MAX { gr = CRESCI_MAX; f = k / (m1 * gr * gr); }
    let mut alfa = 1.0;
    let fmin = (N_MIN / n).min(1.0);
    // a far atom stays a visible knot: its glows are dimmed at most to ALFA_MIN (a few pixels may saturate)
    if f < fmin { alfa = (k / (m1 * fmin * gr * gr)).max(ALFA_MIN); f = fmin; }
    Lod { f: f as f32, alfa: alfa as f32, cresci: gr as f32, sovrapposizione: m1 * f * gr * gr * alfa, errore_px: spaziatura(f) }
}

/// The stable key of glow i in [0, 1): the same integer hash (PCG) as `chiave` in istanze.wgsl.
#[inline]
pub fn chiave(i: u32) -> f32 {
    let s = i.wrapping_mul(747_796_405).wrapping_add(2_891_336_453);
    let w = ((s >> ((s >> 28) + 4)) ^ s).wrapping_mul(277_803_737);
    ((w >> 22) ^ w) as f32 / 4_294_967_296.0
}

/// Weight of glow i under the kept fraction f: 1 well inside, fading to 0 at the edge (no popping while zooming).
#[inline]
pub fn peso(i: u32, f: f32) -> f32 {
    if f >= 1.0 { return 1.0; }
    ((f - chiave(i)) / (0.15 * f + 1e-6)).clamp(0.0, 1.0)
}

// ── THE SCENE AS A GRAPH OF LEVELS, AND THE FRAME BUDGET ──────────────────────────────────────────────────────────────
//
// Each node is a level of the world that can provide glows. Its COARSE state is stored (the atom: ψ on its lattice);
// its FINE detail (the surfers, their tails, the sparkles of the light) is generated from it on demand, as much as the
// screen asks. Walking the graph from the root: a node whose screen-space error is below a pixel is drawn by its own
// coarse glows and its children are not visited; otherwise its children are refined. Today the graph has one node, the
// atom (its surfers are its glows); a molecule or a droplet will be parents that aggregate atoms the same way.

/// A node of the scene graph: a level of the world, its glows at full detail, its sub-levels.
#[derive(Clone, Debug)]
pub struct Nodo {
    pub nome: &'static str,
    /// centre (scene units) and the group of glows of this level
    pub centro: [f64; 3],
    pub gruppo: Gruppo,
    pub figli: Vec<Nodo>,
}

/// What the walk decided for one node.
#[derive(Clone, Copy, Debug)]
pub struct Decisione {
    pub nome: &'static str,
    pub lod: Lod,
    /// the node asks its children to be drawn (its detail on the screen is coarser than `soglia_px`)
    pub raffina: bool,
}

/// Walk the graph: every visited node gets its level of detail; children are visited only where the parent asks.
pub fn percorri(nodo: &Nodo, bpp: f64, s_px: f64, forza: f64, soglia_px: f64, out: &mut Vec<Decisione>) {
    let l = lod(nodo.gruppo, bpp, s_px, forza);
    let raffina = l.raffina(soglia_px) && !nodo.figli.is_empty();
    out.push(Decisione { nome: nodo.nome, lod: l, raffina });
    if raffina { for f in &nodo.figli { percorri(f, bpp, s_px, forza, soglia_px, out); } }
}

/// The frame budget: when frames take longer than the target, the level of detail is made stronger (fewer glows) before
/// the frame rate drops; when there is room it relaxes back. Never in unattended pictures (they must be reproducible).
#[derive(Resource, Clone, Copy, Debug)]
pub struct Budget {
    /// target time per frame (ms)
    pub obiettivo_ms: f64,
    /// smoothed measured time per frame (ms)
    pub misurato_ms: f64,
    /// extra strength of the level of detail (1 = none)
    pub fattore: f64,
}

impl Default for Budget { fn default() -> Self { Budget { obiettivo_ms: 33.0, misurato_ms: 16.0, fattore: 1.0 } } }

impl Budget {
    /// One frame of `dt_ms`: update the measure and the factor (1 … 8).
    pub fn aggiorna(&mut self, dt_ms: f64) {
        self.misurato_ms += 0.05 * (dt_ms.min(500.0) - self.misurato_ms);
        if self.misurato_ms > 1.1 * self.obiettivo_ms { self.fattore *= 1.02; }
        else if self.misurato_ms < 0.8 * self.obiettivo_ms { self.fattore /= 1.01; }
        self.fattore = self.fattore.clamp(1.0, 8.0);
    }
}

pub fn budget_sys(time: Res<Time>, auto: Res<crate::Auto>, mut b: ResMut<Budget>) {
    if auto.attivo { b.fattore = 1.0; return; }
    b.aggiorna(time.delta_secs_f64() * 1e3);
}

/// Keys: G opens the settings, R toggles the slow rotation. Then the settings window.
pub fn pannello(mut ctx: EguiContexts, mut g: ResMut<Grafica>, mut lab: ResMut<Lab>, keys: Res<ButtonInput<KeyCode>>, auto: Res<crate::Auto>,
                mut budget: ResMut<Budget>) {
    let ctx = ctx.ctx_mut();
    if !ctx.wants_keyboard_input() && !auto.attivo {
        if keys.just_pressed(KeyCode::KeyG) { g.aperto = !g.aperto; }
        if keys.just_pressed(KeyCode::KeyR) { g.auto_ruota = !g.auto_ruota; }
    }
    if auto.attivo && !lab.panel_open { return; }
    egui::Area::new(egui::Id::new("grafica_bottone")).anchor(egui::Align2::RIGHT_TOP, [-10.0, 10.0]).show(ctx, |ui| {
        let b = egui::Button::new(RichText::new("⚙ grafica (G)").size(14.0)).fill(Color32::from_black_alpha(90)).selected(g.aperto);
        if ui.add(b).on_hover_text("impostazioni grafiche · R rotazione lenta").clicked() { g.aperto = !g.aperto; }
    });
    if !g.aperto { return; }
    let mut open = true;
    let gold = Color32::from_rgb(255, 205, 80);
    let grey = Color32::from_rgb(150, 155, 170);
    egui::Window::new("Impostazioni grafica").open(&mut open).default_pos([1180.0, 50.0]).default_width(380.0).show(ctx, |ui| {
        ui.label(RichText::new("QUALITÀ").strong().color(gold));
        for (i, q) in QUALITA.iter().enumerate() {
            if ui.selectable_label(g.preset == i, q.2).clicked() { g.preset = i; lab.fot_par.qual = i; }
        }
        ui.label(RichText::new(format!("Cambia solo la risoluzione: celle per λ e lato della scatola del fotone (vale al prossimo \
            lancio), e quanti bagliori di luce si disegnano al massimo ({}). La regola e le sue costanti non cambiano. \
            Il reticolo dell'atomo (48³, celle 0,5 ħ/mc) è fissato dalla fisica, non dal preset.", g.budget_luce())).small().color(grey));
        ui.separator();
        ui.label(RichText::new("DISEGNO").strong().color(gold));
        ui.add(egui::Slider::new(&mut g.scia_s, 0.0..=1.0).text("scia dei surfisti (s)").max_decimals(2));
        ui.add(egui::Slider::new(&mut g.lod, 0.0..=3.0).text("livello di dettaglio (0 = tutto)").max_decimals(1));
        ui.add(egui::Slider::new(&mut g.bagliori, 0.3..=3.0).text("luminosità dei bagliori").max_decimals(1));
        ui.add(egui::Slider::new(&mut g.profondita, 0.0..=1.0).text("profondità (vicino più chiaro)").max_decimals(2));
        ui.checkbox(&mut g.auto_ruota, "rotazione lenta quando non tocchi niente (R)");
        ui.checkbox(&mut g.taglio, "luce attorno all'atomo aperta verso di te (3-D)");
        if g.taglio { ui.add(egui::Slider::new(&mut g.taglio_prof, 20.0..=300.0).text("spessore visibile dietro (ħ/mc)").max_decimals(0)); }
        ui.checkbox(&mut g.fetta, "vecchia fetta piatta della vista diffusa (2-D)");
        ui.checkbox(&mut g.zoom_atomo, "vecchia nuvola ingrandita ×6 al centro");
        ui.add(egui::Slider::new(&mut budget.obiettivo_ms, 12.0..=100.0).text("tempo per fotogramma voluto (ms)").max_decimals(0));
        ui.label(RichText::new(format!("ora {:.0} ms per fotogramma; dettaglio ridotto ×{:.1} per restare fluido", budget.misurato_ms, budget.fattore)).small().color(grey));
        ui.label(RichText::new("Da lontano si vedono gli stessi bagliori, meno e un po' più grandi: la luce per pixel resta la \
            stessa, così niente palla bianca. Trascina per girare, rotella per lo zoom.").small().color(grey));
    });
    if !open { g.aperto = false; }
}

#[cfg(test)]
mod prove {
    use super::*;

    /// Near: everything kept; far: fewer, bigger glows, overlap never above the target, at least N_MIN glows,
    /// and the error grows smoothly with distance.
    #[test]
    fn lod_continuo() {
        let g = Gruppo { n: 20_000, raggio: 3.0, scia: 0.0 };
        let px = |r: f64| 2.0 * r * (0.4f64).tan() / 900.0;
        let vicino = lod(g, px(11.0), 3.7, 1.0);
        assert_eq!(vicino.f, 1.0);
        let mut prima = vicino;
        for k in 1..60 {
            let r = 11.0 * 1.08f64.powi(k);
            let l = lod(g, px(r), 3.7, 1.0);
            assert!(l.sovrapposizione <= SOVRAPPOSIZIONE * 1.0001 || (l.alfa as f64 - ALFA_MIN).abs() < 1e-6, "r {r}: {l:?}");
            assert!(l.f <= prima.f + 1e-6 && l.cresci >= prima.cresci - 1e-6, "r {r}: {l:?} dopo {prima:?}");
            assert!(l.f as f64 * 20_000.0 >= N_MIN - 1e-2);
            prima = l;
        }
        assert!(prima.f < 0.05 && prima.cresci > 1.5);
        assert_eq!(lod(g, px(400.0), 3.7, 0.0).f, 1.0);
    }

    #[test]
    fn grafo_e_budget() {
        let atomo = Nodo { nome: "atomo", centro: [0.0; 3], gruppo: Gruppo { n: 20_000, raggio: 3.0, scia: 0.0 }, figli: vec![] };
        let goccia = Nodo { nome: "goccia", centro: [0.0; 3], gruppo: Gruppo { n: 5_000, raggio: 300.0, scia: 0.0 }, figli: vec![atomo] };
        let mut v = vec![];
        percorri(&goccia, 1000.0, 3.7, 1.0, 1.0, &mut v);
        assert_eq!(v.len(), 1, "da lontano la goccia basta: {v:?}");
        v.clear();
        percorri(&goccia, 0.01, 3.7, 1.0, 1.0, &mut v);
        assert_eq!(v.len(), 2, "da vicino si scende all'atomo: {v:?}");
        let mut b = Budget::default();
        for _ in 0..400 { b.aggiorna(60.0); }
        assert!(b.fattore > 4.0);
        for _ in 0..2000 { b.aggiorna(10.0); }
        assert_eq!(b.fattore, 1.0);
    }

    #[test]
    fn chiave_uniforme() {
        let n = 100_000;
        let sotto = (0..n).filter(|&i| chiave(i) < 0.25).count() as f64 / n as f64;
        assert!((sotto - 0.25).abs() < 0.01, "{sotto}");
        assert_eq!(peso(5, 1.0), 1.0);
    }
}

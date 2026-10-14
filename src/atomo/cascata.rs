//! HYDROGEN with several orbitals: a photon is absorbed (1s → 2p, or 1s → 3p) and the atom comes back down in a
//! CASCADE (2p → 1s; 3p → 1s or 3p → 2s, and the 2s stays: it has no one-photon way down). docs/ATOMO.md §10.4.
//!
//! What is computed (at start, ~1 s, from the zig-zag rule only — inputs α and the electron mass):
//!   • the LEVELS 1s, 2s, 2p, 3s, 3p, 3d: the radial form of the same rule (src/zigzag/walk.rs `Radial`, equal
//!     cells of 0.5 ħ/mc, Z = 1, the Coulomb phase −Zα/r), stationary states by "listening with a phase per link";
//!   • the LINES: their differences (Lyman-α, Lyman-β, Balmer-α ...), as energies and wavelengths;
//!   • the RATES: Einstein A = (4/3) α ω³ |⟨f|r|i⟩|² of each one-photon (dipole) line, with the radial integral
//!     taken over the rule's own radial functions (large component g; the small f is ~α/2 of g and left out:
//!     declared) and the exact angular factor l_max/(2l+1); hence the lifetimes and the branching ratios.
//!
//! What is drawn, NOT simulated in time (declared on screen):
//!   • the electron as glows = surfers carried by the current of the wave (v = j/ρ), plus a small Metropolis
//!     correction so they keep |ψ|² while the drawn state changes faster than the real one;
//!   • absorption: the photon is aimed and on the line, it is absorbed with certainty (a director's choice: in
//!     nature a single photon is absorbed with a probability set by the cross-section); during its passage the
//!     state goes from 1s to 2p following the Rabi formula of a π pulse;
//!   • emission: WHEN the jump happens is drawn at random with the computed lifetime (exponential), WHICH line with
//!     the computed branching ratios; the jump itself lasts a fixed screen time; the outgoing light is a glow shell
//!     with the dipole pattern sin²θ, its rings at the line's wavelength shrunk to the size of the atom;
//!   • two clocks: the wave's beat (ω₂ − ω₁ ~ 10¹⁶ rad/s) is shown slowed so that it is seen, the waiting times
//!     (nanoseconds) on the slow-motion clock chosen in the panel.

use std::sync::Arc;

use rayon::prelude::*;

use crate::cplx::C;
use crate::orbitali::{omega_km, Radiale};
use crate::render::{phase_color, schermo, Ball};
use crate::rng::Rng;
use crate::walk::{surfer_velocity, Radial};
use crate::zz::{dirac_a_chirale, ALPHA, MC2_KEV, T_UNIT_S};

/// bohr in λ_C (1/α).
pub const BOHR: f64 = 1.0 / ALPHA;
/// cell of the radial rule (λ_C) and the box (bohr): the 3s holds 99.99 % within 36 a0
pub const CELLA: f64 = 1.0 / 32.0;
pub const R_SCATOLA_BOHR: f64 = 36.0;
pub const HC_EV_NM: f64 = 1239.841_984;

/// The levels used: (name, n, l). j = l + ½ (the fine-structure doublets, ~10⁻⁵ eV apart, are not split: declared).
pub const LIVELLI: [(&str, u32, u32); 6] = [("1s", 1, 0), ("2s", 2, 0), ("2p", 2, 1), ("3s", 3, 0), ("3p", 3, 1), ("3d", 3, 2)];

#[derive(Clone)]
pub struct Livello { pub nome: &'static str, pub n: u32, pub l: u32, pub kappa: i32, pub rad: Radiale, pub dirac: f64 }

impl Livello {
    pub fn omega(&self) -> f64 { self.rad.omega }
    /// Binding energy (eV).
    pub fn legame_ev(&self) -> f64 { (1.0 - self.rad.omega) * MC2_KEV * 1e3 }
}

/// One one-photon line i → f with its rate.
#[derive(Clone, Debug)]
pub struct Riga { pub da: usize, pub a: usize, pub ev: f64, pub nm: f64, pub a_s: f64, pub nome: String }

/// All the hydrogen numbers of the scene.
#[derive(Clone)]
pub struct Spettro { pub liv: Vec<Livello>, pub righe: Vec<Riga> }

/// Exact Dirac energy (mc²) of level n, κ for a given Zα (comparison only, never used by the scene).
pub fn dirac_esatto_za(n: u32, kappa: i32, za: f64) -> f64 {
    let k = kappa.abs() as f64;
    let d = k - (k * k - za * za).sqrt();
    1.0 / (1.0 + (za / (n as f64 - d)).powi(2)).sqrt()
}

/// The n-th state of partial wave κ of the radial rule at Z = 1 (same recipe as `orbitali::radiale`, bindings of
/// hydrogen: from 3·10⁻⁵ to 2·10⁻⁶ mc², log-spaced).
pub fn radiale_h(kappa: i32, n: u32, l: u32) -> Radiale {
    // 1) find the state on cells of 0.5 ħ/mc (cheap): listen at many bindings, keep the distinct states
    let w = Radial::new(kappa, 0.5, R_SCATOLA_BOHR * BOHR, ALPHA);
    let mut trovati: Vec<f64> = Vec::new();
    for k in 0..28 {
        let b = 3.6e-5 * (2.4e-6f64 / 3.6e-5).powf(k as f64 / 27.0);
        let (om, _, res) = w.ascolta(1.0 - b);
        if om >= 1.0 || om <= 0.0 || !(res < 1e-8) { continue; }
        let bind = 1.0 - om;
        if trovati.iter().all(|t| ((1.0 - t) - bind).abs() > 0.005 * bind) { trovati.push(om); }
    }
    trovati.sort_by(|a, b| a.total_cmp(b));
    let mut omega = trovati.get((n - l - 1) as usize).copied().unwrap_or(1.0);
    // 2) follow it down to cells of CELLA (θ = mc·a/ħ → small: on cells of 0.5 the walk's kinetic mass is
    //    tan θ/θ = 1.093, +9 % on every binding; at 1/32 it is 1.0003), each time listening from the last value
    let mut a = 0.5;
    loop {
        a *= 0.5;
        let w = Radial::new(kappa, a, R_SCATOLA_BOHR * BOHR, ALPHA);
        let (om, gf, _) = w.ascolta(omega);
        omega = om;
        if a <= CELLA * 1.0001 { return Radiale { omega, gf, r: w.r.clone(), h: w.h.clone() }; }
    }
}

/// ⟨f|r|i⟩ radial (λ_C), large components.
pub fn radiale_r(a: &Radiale, b: &Radiale) -> f64 {
    a.gf.iter().zip(&b.gf).zip(a.r.iter().zip(&a.h)).map(|((x, y), (r, h))| x[0] * y[0] * r * h).sum()
}

impl Spettro {
    pub fn calcola() -> Self {
        let liv: Vec<Livello> = LIVELLI.par_iter().map(|&(nome, n, l)| {
            let kappa = -(l as i32 + 1);
            Livello { nome, n, l, kappa, rad: radiale_h(kappa, n, l), dirac: dirac_esatto_za(n, kappa, ALPHA) }
        }).collect();
        let mut righe = vec![];
        for i in 0..liv.len() {
            for f in 0..liv.len() {
                let (a, b) = (&liv[i], &liv[f]);
                // same n: microwave lines (2p→2s ~10⁻³ s⁻¹), whose energy here is mostly the rule's s-level shift: left out
                if b.omega() >= a.omega() || (a.l as i32 - b.l as i32).abs() != 1 || a.n == b.n { continue; }
                let w = a.omega() - b.omega();
                let r = radiale_r(&a.rad, &b.rad);
                let ang = a.l.max(b.l) as f64 / (2 * a.l + 1) as f64;
                let a_unit = 4.0 / 3.0 * ALPHA * w.powi(3) * r * r * ang;
                let ev = w * MC2_KEV * 1e3;
                righe.push(Riga { da: i, a: f, ev, nm: HC_EV_NM / ev, a_s: a_unit / T_UNIT_S, nome: nome_riga(a.n, b.n) });
            }
        }
        Spettro { liv, righe }
    }
    /// Total one-photon decay rate of level i (1/s); 0 for 1s and 2s.
    pub fn gamma(&self, i: usize) -> f64 { self.righe.iter().filter(|r| r.da == i).map(|r| r.a_s).sum() }
    pub fn riga(&self, da: usize, a: usize) -> Option<&Riga> { self.righe.iter().find(|r| r.da == da && r.a == a) }
}

fn nome_riga(n1: u32, n2: u32) -> String {
    let serie = match n2 { 1 => "Lyman", 2 => "Balmer", 3 => "Paschen", _ => "?" };
    let g = ["α", "β", "γ", "δ"][((n1 - n2 - 1) as usize).min(3)];
    format!("{serie}-{g}")
}

/// Real colour of a wavelength (nm) where the eye sees it; ultraviolet drawn violet-white (declared false colour).
pub fn colore_nm(nm: f64) -> [f32; 3] {
    let l = nm as f32;
    let (r, g, b) = if l < 380.0 { (0.62, 0.45, 1.0) } else if l < 440.0 { ((440.0 - l) / 60.0, 0.0, 1.0) } else if l < 490.0 { (0.0, (l - 440.0) / 50.0, 1.0) }
        else if l < 510.0 { (0.0, 1.0, (510.0 - l) / 20.0) } else if l < 580.0 { ((l - 510.0) / 70.0, 1.0, 0.0) } else if l < 645.0 { (1.0, (645.0 - l) / 65.0, 0.0) }
        else if l < 780.0 { (1.0, 0.0, 0.0) } else { (0.7, 0.1, 0.05) };
    [r.max(0.12), g.max(0.08), b.max(0.08)]
}

// ───────────────────────────── the scene ─────────────────────────────

/// What the panel asks the hydrogen scene.
#[derive(Clone, Debug, Default)]
pub struct Richiesta {
    /// fire a photon on this line (index of the upper level reached from the current one)
    pub lancia: Option<usize>,
    /// ns of the waiting times per screen second
    pub ns_al_secondo: Option<f64>,
    pub ricomincia: bool,
}

/// A photon on its way in (along +x, polarised along z) or a glow shell on its way out.
#[derive(Clone, Debug)]
struct Luce { t0: f64, nm: f64, uscente: bool, durata: f64 }

#[derive(Clone, Copy, Debug, PartialEq)]
enum Fase {
    /// in a level, waiting (decay drawn at `salto`, screen seconds, if it decays)
    Fermo { liv: usize, salto: Option<f64>, verso: Option<usize> },
    /// the incoming photon is passing: from `da` to `a`, started at t0
    Assorbe { da: usize, a: usize, t0: f64 },
    /// the jump down from `da` to `a`, started at t0
    Emette { da: usize, a: usize, t0: f64 },
}

/// Screen seconds: the photon's passage and the jump down.
pub const DURATA_ASSORBE: f64 = 2.5;
pub const DURATA_SALTO: f64 = 2.0;
/// Shown beat: one period of the wave's beat in this many screen seconds (the real one is ~0.4 fs).
pub const BATTITO_SCHERMO: f64 = 1.2;
pub const N_SURF: usize = 4500;
/// speed of the light shells on screen (bohr per screen second) and their life
const V_LUCE: f64 = 9.0;
const VITA_LUCE: f64 = 7.0;

pub struct Mondo {
    pub sp: Arc<Spettro>,
    fase: Fase,
    t: f64,
    /// the wave's own clock (λ_C/c): the phases e^{−i(ω−1)τ}
    tau: f64,
    surf: Vec<[f64; 3]>,
    luci: Vec<Luce>,
    pub ns_al_secondo: f64,
    rng: Rng,
    /// lines emitted so far (for the screen)
    pub emesse: Vec<String>,
    /// the light was fired on a line the current level cannot absorb
    pub avviso: String,
    /// screen time of the last photon fired
    pub lancio: Option<f64>,
}

/// Amplitudes of the drawn state: (level, complex coefficient without the e^{−iωτ}).
type Stato = Vec<(usize, f64)>;

impl Mondo {
    pub fn new() -> Self { Self::con(Arc::new(Spettro::calcola())) }
    pub fn con(sp: Arc<Spettro>) -> Self {
        let mut rng = Rng::new(2026);
        let surf = (0..N_SURF).map(|_| { let r = esponenziale_1s(&mut rng); r }).collect();
        Mondo { sp, fase: Fase::Fermo { liv: 0, salto: None, verso: None }, t: 0.0, tau: 0.0, surf, luci: vec![], ns_al_secondo: 1.0, rng, emesse: vec![], avviso: String::new(), lancio: None }
    }

    pub fn livello(&self) -> usize { match self.fase { Fase::Fermo { liv, .. } => liv, Fase::Assorbe { a, .. } => a, Fase::Emette { a, .. } => a } }

    /// The drawn superposition now.
    fn stato(&self) -> Stato {
        let s = |x: f64| x.clamp(0.0, 1.0);
        match self.fase {
            Fase::Fermo { liv, .. } => vec![(liv, 1.0)],
            Fase::Assorbe { da, a, t0 } => {
                // Rabi π pulse: θ grows with the pulse area, smooth (sin² envelope) over the passage
                let u = s((self.t - t0) / DURATA_ASSORBE);
                let th = std::f64::consts::PI * (u - (std::f64::consts::TAU * u).sin() / std::f64::consts::TAU);
                vec![(da, (0.5 * th).cos()), (a, (0.5 * th).sin())]
            }
            Fase::Emette { da, a, t0 } => {
                let u = s((self.t - t0) / DURATA_SALTO);
                let th = std::f64::consts::PI * (u - (std::f64::consts::TAU * u).sin() / std::f64::consts::TAU);
                vec![(da, (0.5 * th).cos()), (a, (0.5 * th).sin())]
            }
        }
    }

    /// ψ (Dirac, 4 components) of the drawn state at p (λ_C), with the wave clock τ.
    fn psi(&self, st: &Stato, p: [f64; 3]) -> [C; 4] {
        let mut o = [C::ZERO; 4];
        for &(i, c) in st {
            if c.abs() < 1e-6 { continue; }
            let lv = &self.sp.liv[i];
            let v = spinore_m(&lv.rad, lv.kappa, 1, p);
            let f = C::expi(-(lv.omega() - 1.0) * self.tau).scale(c);
            for k in 0..4 { o[k] += v[k] * f; }
        }
        o
    }

    pub fn richiesta(&mut self, r: Richiesta) {
        if let Some(v) = r.ns_al_secondo { self.ns_al_secondo = v.max(1e-3); }
        if r.ricomincia { *self = Mondo { ns_al_secondo: self.ns_al_secondo, ..Mondo::con(self.sp.clone()) }; return; }
        if let Some(alto) = r.lancia {
            match self.fase {
                Fase::Fermo { liv, .. } if self.sp.riga(alto, liv).is_some() => {
                    let nm = self.sp.riga(alto, liv).unwrap().nm;
                    self.luci.push(Luce { t0: self.t, nm, uscente: false, durata: DURATA_ASSORBE });
                    self.fase = Fase::Assorbe { da: liv, a: alto, t0: self.t + ARRIVO };
                    self.lancio = Some(self.t);
                    self.avviso.clear();
                }
                Fase::Fermo { liv, .. } => {
                    self.avviso = format!("dal {} non c'è una riga a un fotone verso il {}: il fotone passa senza essere assorbito", self.sp.liv[liv].nome, self.sp.liv[alto].nome);
                    let nm = HC_EV_NM / ((self.sp.liv[alto].omega() - self.sp.liv[liv].omega()).abs() * MC2_KEV * 1e3).max(1e-9);
                    self.luci.push(Luce { t0: self.t, nm, uscente: false, durata: 0.0 });
                }
                _ => { self.avviso = "aspetta: l'atomo sta già cambiando livello".into(); }
            }
        }
    }

    /// One frame of `dt` screen seconds.
    pub fn passo(&mut self, dt: f64) {
        self.t += dt;
        // the wave's clock: the largest beat among the drawn levels in BATTITO_SCHERMO seconds
        let w21 = self.sp.liv[2].omega() - self.sp.liv[0].omega();
        self.tau += dt * std::f64::consts::TAU / w21 / BATTITO_SCHERMO;
        // state machine
        match self.fase {
            Fase::Assorbe { a, t0, .. } if self.t >= t0 + DURATA_ASSORBE => { self.fase = self.attesa(a); }
            Fase::Emette { a, t0, .. } if self.t >= t0 + DURATA_SALTO => { self.fase = self.attesa(a); }
            Fase::Fermo { liv, salto: Some(ts), verso: Some(f) } if self.t >= ts => {
                let r = self.sp.riga(liv, f).cloned();
                if let Some(r) = r {
                    self.luci.push(Luce { t0: self.t, nm: r.nm, uscente: true, durata: DURATA_SALTO });
                    self.emesse.push(format!("{} {}→{} {:.1} nm", r.nome, self.sp.liv[liv].nome, self.sp.liv[f].nome, r.nm));
                    if self.emesse.len() > 6 { self.emesse.remove(0); }
                }
                self.fase = Fase::Emette { da: liv, a: f, t0: self.t };
            }
            _ => {}
        }
        self.luci.retain(|l| self.t - l.t0 < VITA_LUCE);
        // the glows: surfers of the drawn state + Metropolis on |ψ|²
        let st = self.stato();
        let dtau = dt * std::f64::consts::TAU / w21 / BATTITO_SCHERMO;
        let seme = self.rng.next_u64();
        let n_max = st.iter().filter(|x| x.1.abs() > 1e-3).map(|x| self.sp.liv[x.0].n).max().unwrap_or(1) as f64;
        let passo_m = 0.3 * n_max * n_max * BOHR;
        let me = &*self;
        let nuovi: Vec<[f64; 3]> = me.surf.par_iter().enumerate().map(|(k, p)| {
            let mut rng = Rng::new(seme ^ (k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let mut p = *p;
            // sub-steps: at most 0.04 of the cylinder radius per step (the 1s current circles z at α sinθ: the inner glows
            // turn fast, and a coarse step would push them outwards or inwards)
            let v0 = surfer_velocity(&dirac_a_chirale(me.psi(&st, p)));
            let rc = (p[0] * p[0] + p[1] * p[1]).sqrt().max(0.03 * BOHR);
            let vv = (v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]).sqrt();
            let sub = ((vv * dtau / (0.04 * rc)).ceil() as usize).clamp(1, 60);
            let h = dtau / sub as f64;
            for _ in 0..sub {
                let v1 = surfer_velocity(&dirac_a_chirale(me.psi(&st, p)));
                let q = [p[0] + 0.5 * h * v1[0], p[1] + 0.5 * h * v1[1], p[2] + 0.5 * h * v1[2]];
                let v2 = surfer_velocity(&dirac_a_chirale(me.psi(&st, q)));
                for c in 0..3 { p[c] += h * v2[c]; }
            }
            // Metropolis: keeps the glows on |ψ|² while the drawn state changes faster than the real one
            // (a step that does not depend on where the glow is, so the proposal is symmetric)
            let passo = passo_m;
            let rho = |x: [f64; 3]| me.psi(&st, x).iter().map(|c| c.norm2()).sum::<f64>();
            let mut r0 = rho(p);
            for _ in 0..2 {
                let q = [p[0] + passo * rng.gauss(), p[1] + passo * rng.gauss(), p[2] + passo * rng.gauss()];
                let r1 = rho(q);
                if r1 >= r0 * rng.uni() { p = q; r0 = r1; }
            }
            p
        }).collect();
        self.surf = nuovi;
    }

    /// Wait in level `liv`: draw the jump time from its computed lifetime and the line from its branching ratios.
    fn attesa(&mut self, liv: usize) -> Fase {
        let g = self.sp.gamma(liv);
        if g <= 0.0 { return Fase::Fermo { liv, salto: None, verso: None }; }
        let tau_ns = 1e9 / g;
        let attesa = -(1.0 - self.rng.uni()).ln() * tau_ns / self.ns_al_secondo;
        let mut u = self.rng.uni() * g;
        let mut verso = None;
        for r in self.sp.righe.iter().filter(|r| r.da == liv) { u -= r.a_s; if u <= 0.0 { verso = Some(r.a); break; } }
        let verso = verso.or_else(|| self.sp.righe.iter().filter(|r| r.da == liv).last().map(|r| r.a));
        Fase::Fermo { liv, salto: Some(self.t + attesa), verso }
    }

    /// Screen seconds since the last photon was fired.
    pub fn dal_lancio(&self) -> Option<f64> { self.lancio.map(|t| self.t - t) }

    /// Radius of the view (bohr).
    pub fn raggio_vista(&self) -> f64 { 18.0 }

    /// Glows (electron, nucleus) and the light (incoming packet, outgoing shells).
    pub fn disegna(&self) -> (Vec<Ball>, Vec<Ball>) {
        let st = self.stato();
        let mut p: Vec<Ball> = self.surf.par_iter().map(|x| {
            let s = self.psi(&st, *x);
            let f = if s[0].norm2() >= s[1].norm2() { s[0].arg() } else { s[1].arg() };
            (schermo([x[0] / BOHR, x[1] / BOHR, x[2] / BOHR]).to_array(), 1.0, phase_color(f, 0.42))
        }).collect();
        p.push((schermo([0.0; 3]).to_array(), 2.2, [1.0, 0.6, 0.25, 0.9]));
        let mut onde = vec![];
        for (k, l) in self.luci.iter().enumerate() {
            let c = colore_nm(l.nm);
            // rings of the wave at the line's wavelength shrunk ×(λ / 2.2 a0) (drawing choice, declared)
            let lam = 2.2;
            let mut rng = Rng::new(k as u64 * 7919 + (l.t0 * 1000.0) as u64);
            let eta = self.t - l.t0;
            if !l.uscente {
                // incoming packet along +x, polarised along z: a slab of glows, fronts every λ_drawn, fading as it is absorbed
                let x0 = -V_LUCE * (ARRIVO - eta);
                let assorbito = if l.durata > 0.0 { ((eta - ARRIVO) / l.durata).clamp(0.0, 1.0) } else { 0.0 };
                let w = 1.0 - assorbito;
                if w <= 0.01 || x0 > 40.0 { continue; }
                for _ in 0..2600 {
                    let x = x0 + 3.0 * lam * rng.gauss() * 0.5;
                    let (y, z) = (4.0 * rng.gauss(), 4.0 * rng.gauss());
                    let fase = std::f64::consts::TAU * (x - x0) / lam;
                    let a = (fase.cos().powi(2) * (-(x - x0).powi(2) / (2.0 * (1.5 * lam).powi(2))).exp()) as f32 * 0.55 * w as f32;
                    if a < 0.02 { continue; }
                    onde.push((schermo([x, y, z]).to_array(), 0.9, [c[0], c[1], c[2], a]));
                }
            } else {
                // outgoing shell: radius V·t, dipole pattern sin²θ about z, brightness ∝ |c_up c_down| of the jump
                let r0 = V_LUCE * eta;
                let amp = (std::f64::consts::PI * (eta / l.durata).clamp(0.0, 1.0)).sin().max(if eta > l.durata { 0.0 } else { 0.05 });
                for _ in 0..3200 {
                    let ct = 2.0 * rng.uni() - 1.0;
                    let sin2 = 1.0 - ct * ct;
                    if rng.uni() > sin2 { continue; }
                    let ph = std::f64::consts::TAU * rng.uni();
                    let st = sin2.sqrt();
                    // the emission lasted DURATA_SALTO: the shell is V·durata thick, with rings every λ_drawn
                    let d = V_LUCE * l.durata.min(eta) * rng.uni();
                    let r = r0 - d;
                    if r < 0.3 { continue; }
                    let inviluppo = (std::f64::consts::PI * d / (V_LUCE * l.durata)).sin();
                    let anello = (std::f64::consts::TAU * d / lam).cos().powi(2);
                    let a = (0.6 * inviluppo * anello) as f32 * (1.0 - (eta / VITA_LUCE) as f32).max(0.0);
                    let _ = amp;
                    if a < 0.02 { continue; }
                    onde.push((schermo([r * st * ph.cos(), r * st * ph.sin(), r * ct]).to_array(), 0.9, [c[0], c[1], c[2], a]));
                }
            }
        }
        (p, onde)
    }

    /// The two lines at the top of the screen.
    pub fn righe(&self) -> (String, String) {
        let sp = &self.sp;
        let liv = |i: usize| sp.liv[i].nome;
        let a = match self.fase {
            Fase::Fermo { liv: l, salto: None, .. } if l == 0 => "IDROGENO nello stato 1s · premi F (o un bottone del pannello) per lanciare un fotone".to_string(),
            Fase::Fermo { liv: l, salto: None, .. } => format!("IDROGENO nel {}: METASTABILE — nessuna riga a un fotone verso il basso; la via a due fotoni (~0,12 s in natura) non è calcolata (debito)", liv(l)),
            Fase::Fermo { liv: l, salto: Some(ts), verso } => format!("IDROGENO nel {} · vita media {:.2} ns (calcolata) · salto tra {:.1} s di schermo verso il {}", liv(l), 1e9 / sp.gamma(l), (ts - self.t).max(0.0), verso.map(liv).unwrap_or("?")),
            Fase::Assorbe { da, a, t0 } => if self.t < t0 { format!("Il fotone arriva: {} → {} ({:.1} nm)", liv(da), liv(a), sp.riga(a, da).map(|r| r.nm).unwrap_or(0.0)) }
                else { format!("ASSORBE {} → {} · {} · l'onda batte alla frequenza della riga (mostrata rallentata)", liv(da), liv(a), sp.riga(a, da).map(|r| format!("{} {:.2} nm, {:.3} eV", r.nome, r.nm, r.ev)).unwrap_or_default()) },
            Fase::Emette { da, a, .. } => format!("RIEMETTE {} → {} · {} (guscio di luce, figura sin²θ attorno al dipolo)", liv(da), liv(a), sp.riga(da, a).map(|r| format!("{} {:.2} nm, {:.3} eV", r.nome, r.nm, r.ev)).unwrap_or_default()),
        };
        let b = format!("rallentatore 1 s = {} ns · battito dell'onda rallentato ×{:.0e} · emesse: {}{}", self.ns_al_secondo,
            BATTITO_SCHERMO * (sp.liv[2].omega() - sp.liv[0].omega()) / T_UNIT_S / std::f64::consts::TAU,
            if self.emesse.is_empty() { "—".to_string() } else { self.emesse.join(", ") },
            if self.avviso.is_empty() { String::new() } else { format!(" · {}", self.avviso) });
        (a, b)
    }
}

/// Screen seconds from the launch to the photon reaching the atom.
const ARRIVO: f64 = 2.0;

impl Default for Mondo { fn default() -> Self { Self::new() } }

/// A point drawn from the 1s density (start of the glows): r with r² e^{−2r} (bohr), as λ_C.
fn esponenziale_1s(rng: &mut Rng) -> [f64; 3] {
    // Gamma(3, 1/2): sum of three exponentials of mean 1/2
    let r = -0.5 * ((1.0 - rng.uni()) * (1.0 - rng.uni()) * (1.0 - rng.uni())).ln();
    let ct = 2.0 * rng.uni() - 1.0;
    let ph = std::f64::consts::TAU * rng.uni();
    let st = (1.0 - ct * ct).sqrt();
    [r * st * ph.cos() * BOHR, r * st * ph.sin() * BOHR, r * ct * BOHR]
}

/// The 3-D Dirac spinor of a radial state (as `orbitali::spinore`, any mesh).
fn spinore_m(r: &Radiale, kappa: i32, m2: i32, p: [f64; 3]) -> [C; 4] {
    let rr = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(1e-9);
    let ns = r.gf.len();
    if rr >= r.r[ns - 1] { return [C::ZERO; 4]; }
    // equal cells: the index directly
    let i = (((rr - r.r[0]) / r.h[0]).floor().max(0.0) as usize).min(ns - 2);
    let t = ((rr - r.r[i]) / (r.r[i + 1] - r.r[i])).clamp(0.0, 1.0);
    let g = (r.gf[i][0] * (1.0 - t) + r.gf[i + 1][0] * t) / rr;
    let f = (r.gf[i][1] * (1.0 - t) + r.gf[i + 1][1] * t) / rr;
    let ct = p[2] / rr;
    let phi = p[1].atan2(p[0]);
    let up = omega_km(kappa, m2, ct, phi);
    let lo = omega_km(-kappa, m2, ct, phi);
    [up[0].scale(g), up[1].scale(g), C::I * lo[0].scale(f), C::I * lo[1].scale(f)]
}

#[cfg(test)]
mod test {
    use super::*;

    /// C1: the levels of the rule at Z = 1 against Dirac; C2: A(2p → 1s) against 6.2649·10⁸ s⁻¹ (exam value, not an
    /// input); C3: branching 3p → 2s against 0.118; Lyman-α and Balmer-α wavelengths printed.
    #[test]
    fn cascata_livelli_e_ritmi() {
        let t = std::time::Instant::now();
        let sp = Spettro::calcola();
        for l in &sp.liv {
            let d = (1.0 - l.dirac) * MC2_KEV * 1e3;
            println!("{}: legame {:.6} eV, Dirac {:.6} eV ({:+.4} %)", l.nome, l.legame_ev(), d, 100.0 * (l.legame_ev() / d - 1.0));
            assert!((l.legame_ev() / d - 1.0).abs() < 1e-3, "{} fuori dallo 0,1 %", l.nome);
        }
        for r in &sp.righe { println!("{} {}→{}: {:.4} eV {:.2} nm A = {:.4e} s⁻¹", r.nome, sp.liv[r.da].nome, sp.liv[r.a].nome, r.ev, r.nm, r.a_s); }
        let a21 = sp.riga(2, 0).unwrap().a_s;
        assert!((a21 / 6.2649e8 - 1.0).abs() < 0.01, "A(2p→1s) = {a21:e}");
        let g3p = sp.gamma(4);
        let b = sp.riga(4, 1).unwrap().a_s / g3p;
        println!("vita 2p {:.3} ns, 3p {:.3} ns, ramo 3p→2s {:.4}; {:.1} s", 1e9 / sp.gamma(2), 1e9 / g3p, b, t.elapsed().as_secs_f64());
        assert!((b / 0.118 - 1.0).abs() < 0.02);
        assert_eq!(sp.gamma(1), 0.0);
        assert!((sp.riga(2, 0).unwrap().nm - 121.50).abs() < 0.2);
    }

    /// The scene: absorb Lyman-α, wait, come back down with a Lyman-α shell.
    #[test]
    fn cascata_scena() {
        let mut m = Mondo::new();
        m.ns_al_secondo = 10.0;
        m.richiesta(Richiesta { lancia: Some(2), ..Default::default() });
        for _ in 0..(60.0 * 12.0) as usize { m.passo(1.0 / 60.0); if m.emesse.len() == 1 { break; } }
        assert_eq!(m.emesse.len(), 1, "{:?}", m.righe());
        assert!(m.emesse[0].starts_with("Lyman-α"));
        for _ in 0..60 { m.passo(1.0 / 60.0); }
        let (p, o) = m.disegna();
        assert!(p.len() > N_SURF && !o.is_empty());
        // the glows sample the drawn state: ⟨r⟩ of 1s is 1.5 a0
        for _ in 0..240 { m.passo(1.0 / 60.0); }
        assert_eq!(m.livello(), 0);
        let r: f64 = m.surf.iter().map(|x| (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt() / BOHR).sum::<f64>() / N_SURF as f64;
        assert!((r / 1.5 - 1.0).abs() < 0.08, "<r> = {r}");
    }
}

impl crate::tavola::Pannello {
    /// The upper level of the line the photon is tuned to (default: 2p, Lyman-α).
    pub fn riga_lancio(&self) -> usize { if self.riga == 0 { 2 } else { self.riga } }
}

/// The hydrogen block of the periodic-table window: photon lines, slow motion, the computed spectrum.
pub fn pannello(ui: &mut bevy_egui::egui::Ui, lab: &mut crate::Lab) {
    use bevy_egui::egui::{self, Color32, RichText};
    let oro = Color32::from_rgb(255, 205, 80);
    let tenue = Color32::from_rgb(150, 155, 170);
    let o = lab.tav.out();
    let Some(sp) = o.as_ref().and_then(|o| o.spettro.clone()) else {
        ui.label(RichText::new("idrogeno: sto calcolando i livelli con la regola…").color(tenue));
        return;
    };
    let ora = o.map(|o| o.livello).unwrap_or(0);
    ui.label(RichText::new("IDROGENO: ASSORBE E RIEMETTE").strong().color(oro));
    ui.label(RichText::new("Scegli la riga del fotone (sono le righe CALCOLATE dalla regola) e premi lancia (F). L'atomo assorbe solo una riga \
        che parte dal livello in cui si trova.").small().color(tenue));
    ui.horizontal_wrapped(|ui| {
        for r in sp.righe.iter().filter(|r| r.ev > 1.0) {
            let (alto, basso) = (r.da, r.a);
            let testo = format!("{} {}→{} {:.1} nm", r.nome, sp.liv[basso].nome, sp.liv[alto].nome, r.nm);
            let c = colore_nm(r.nm);
            let col = Color32::from_rgb((c[0] * 200.0) as u8, (c[1] * 200.0) as u8, (c[2] * 200.0) as u8);
            let ok = basso == ora;
            let b = egui::Button::new(RichText::new(testo).color(Color32::WHITE)).fill(if ok { col } else { Color32::from_rgb(45, 45, 55) }).selected(lab.tav.riga_lancio() == alto);
            if ui.add(b).on_hover_text(if ok { "lancia questo fotone" } else { "l'atomo non è nel livello di partenza di questa riga: il fotone passerebbe senza essere assorbito" }).clicked() {
                lab.tav.riga = alto;
                lab.spara = true;
            }
        }
    });
    if lab.tav.ns_al_secondo <= 0.0 { lab.tav.ns_al_secondo = 1.0; }
    ui.horizontal(|ui| {
        let r = ui.add(egui::Slider::new(&mut lab.tav.ns_al_secondo, 0.05..=50.0).logarithmic(true).text("ns per secondo (rallentatore)"));
        if r.changed() { let v = lab.tav.ns_al_secondo; lab.tav.richiedi(|q| q.ns_al_secondo = Some(v)); }
        if ui.button("ricomincia (1s)").clicked() { lab.tav.richiedi(|q| q.ricomincia = true); }
    });
    egui::CollapsingHeader::new(RichText::new("livelli e righe calcolati (regola zig-zag radiale, Z = 1)").color(oro)).default_open(false).show(ui, |ui| {
        egui::Grid::new("livelli_h").striped(true).show(ui, |ui| {
            ui.label("livello"); ui.label("legame, regola"); ui.label("Dirac (nucleo fermo)"); ui.label("vita"); ui.end_row();
            for (i, l) in sp.liv.iter().enumerate() {
                let d = (1.0 - l.dirac) * MC2_KEV * 1e3;
                ui.label(l.nome);
                ui.label(format!("{:.5} eV", l.legame_ev()));
                ui.label(format!("{:.5} eV ({:+.3} %)", d, 100.0 * (l.legame_ev() / d - 1.0)));
                let g = sp.gamma(i);
                ui.label(if i == 0 { "stabile".into() } else if g > 0.0 { format!("{:.3} ns", 1e9 / g) } else { "metastabile (2 fotoni: non calcolato)".into() });
                ui.end_row();
            }
        });
        egui::Grid::new("righe_h").striped(true).show(ui, |ui| {
            ui.label("riga"); ui.label("λ"); ui.label("A (1/s)"); ui.label("ramo"); ui.end_row();
            for r in &sp.righe {
                ui.label(format!("{} {}→{}", r.nome, sp.liv[r.da].nome, sp.liv[r.a].nome));
                ui.label(format!("{:.2} nm", r.nm));
                ui.label(format!("{:.4e}", r.a_s));
                ui.label(format!("{:.1} %", 100.0 * r.a_s / sp.gamma(r.da)));
                ui.end_row();
            }
        });
        ui.label(RichText::new("Nucleo fermo (massa infinita): con il rinculo del protone (la sua massa è un ingresso ammesso, non ancora messo) \
            tutte le λ crescono dello 0,054 % (Lyman-α 121,57 nm misurata). Doppietti di struttura fine non separati (j = l + ½).").small().color(tenue));
    });
    ui.label(RichText::new("CALCOLATO: livelli, righe, vite medie, rami (A = 4/3 α ω³ |⟨f|r|i⟩|² sulle funzioni radiali della regola). \
        DISEGNATO: il fotone puntato è assorbito con certezza (regia), il passaggio 1s→2p segue un impulso π (Rabi), il momento del salto è \
        pescato a caso con la vita calcolata, il salto dura 2 s di schermo, il guscio di luce ha anelli alla λ della riga rimpicciolita fino \
        alla taglia dell'atomo (figura sin²θ del dipolo vera), il battito dell'onda è rallentato. I bagliori sono surfisti della corrente \
        dell'onda (v = j/ρ) con una correzione Metropolis su |ψ|².").small().color(tenue));
}

//! LIGHT AND ATOM COUPLED, both with the zig-zag rule (semi-classical Maxwell–Dirac on the lattice).
//!
//! • The LIGHT around the atom is a box of spin-1 walkers (fotone3d.rs):
//!   128³ cells of b = 4 λ_C (512 λ_C, ~6.5 wavelengths of the 1s → 2p light), absorbing sponge at the faces.
//!   One light tick = b/c = 8 atom ticks.
//! • The ELECTRON is the Dirac walk of zz.rs (48³ cells of 0.5 λ_C at the centre of the light box).
//! • Light → electron: the electric field of the light at the atom (a Gaussian average over σ = 2 light cells)
//!   enters the electron's phase as e E·r (dipole form: the light wavelength is ~55 λ_C, the atom ~5).
//! • Electron → light: the electron's current I = −e ⟨α⟩ (from the spinor, averaged over the 8 atom ticks of a
//!   light tick) is a source of the light: ∂F/∂t = −ic∇×F − J, with J spread over the same Gaussian.
//!   So the light the atom sends out, and the atom's relaxation, come out of the coupled rule: nothing scripted.
//! • Semi-classical: only an atom in a SUPERPOSITION (1s + 2p, an oscillating dipole) radiates. A pure excited
//!   state would not decay (true spontaneous emission needs the quantized light field): a debt.
//! • The coupling is made K times stronger (e → K e, both ways) so that the decay is seen in about a minute
//!   instead of hours; the decay rate scales as K². The checks compare with formulas that use the same K.
//! • Timing of the coupling (docs/EMISSIONE.md §3.2, docs/ATOMO.md §3.7): `Accoppiamento::RITARDATO` (the original
//!   one) reads the field BEFORE the light tick and deposits AFTER it, non-transverse: the electron's own near field
//!   comes back half a tick late, a fake friction. `Accoppiamento::CENTRATO` reads the centred field ½(U F_n + F_{n+1})
//!   (predictor of the current + corrector) and deposits only the TRANSVERSE part of the source, P_⊥(k) = 1 − n nᵀ with
//!   n(k) the axis of the light rule's two-tick Bloch rotation (its exact longitudinal mode), computed once by FFT.

#![allow(dead_code)]

use rayon::prelude::*;

use crate::fotone3d::{Fotone3D, Lancio, Polar};

/// Light box: cells per side, cell size (λ_C), sponge width (cells), smearing σ (light cells), coupling boost.
pub const N_LUCE: usize = 128;
pub const B_LUCE: f64 = 4.0;
pub const SPUGNA: usize = 12;
pub const SIGMA_S: f64 = 2.0;
pub const K_ACC: f64 = 10.0;
/// Atom ticks per light tick (b / a).
pub const SUB: usize = 8;
/// Backward probe (cells from the centre, along −x).
pub const SONDA: usize = 40;

/// How the atom and the light exchange field and current within one light tick.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Accoppiamento {
    /// the atom sees ½(U F_n + F_{n+1}) (predictor of the current, then a corrector phase) instead of F_n
    pub centrato: bool,
    /// the deposit is projected with P_⊥(k) = 1 − n(k) n(k)ᵀ (Coulomb gauge of the lattice light)
    pub trasverso: bool,
}
impl Accoppiamento {
    /// the original coupling (until 2026-10-14): field before the light tick, deposit after, not transverse
    pub const RITARDATO: Accoppiamento = Accoppiamento { centrato: false, trasverso: false };
    /// the patch of docs/EMISSIONE.md §7
    pub const CENTRATO: Accoppiamento = Accoppiamento { centrato: true, trasverso: true };
    pub fn nome(self) -> &'static str {
        match (self.centrato, self.trasverso) {
            (false, false) => "ritardato",
            (true, false) => "centrato",
            (false, true) => "trasverso",
            (true, true) => "centrato+trasverso",
        }
    }
}
/// The coupling used by `Esperimento::new` / `nuovo` (the Laboratorio and the export). CENTRATO since 2026-10-14
/// (lab/out/emissione_lab_results.json: D1 1.004, E1 0.968, F1, K1 1.006 PASS; L1 FAIL; rule of
/// lab/out/emissione_lab_amend1.json, changed after the results and declared).
pub const ACCOPPIAMENTO: Accoppiamento = Accoppiamento::CENTRATO;

/// The 1s of `Esperimento::new` / `nuovo` kept in the electron's valley k = 0 (whole cube only; docs/ATOMO.md §3.7).
pub const VALLE_ZERO: bool = true;

pub fn alpha() -> f64 { crate::zz::ALPHA }
/// The coupling charge e = K √(4πα) (Heaviside–Lorentz, ħ = c = m = 1).
pub fn e_eff() -> f64 { K_ACC * (4.0 * std::f64::consts::PI * alpha()).sqrt() }

pub struct LuceAtomo {
    pub f: Fotone3D,
    maschera: Vec<f32>,
    /// cells of the Gaussian average around the atom and their weights (Σ = 1)
    pesi: Vec<(usize, f64)>,
    /// time series at the backward probe: (time λ_C/c, E_z)
    pub sonda: Vec<(f64, f64)>,
    /// energy taken by the sponge so far, and (t, light energy in the box + absorbed) every light tick
    pub assorbita: f64,
    pub serie_e: Vec<(f64, f64)>,
    pub t: f64,
    /// transverse deposit patterns (None: the plain Gaussian deposit)
    pub trasv: Option<std::sync::Arc<Trasverso>>,
    /// EXACT SCATTERED FIELD (docs/ATOMO.md §3.8): when Some, `f` holds only the FREE light (the pulse, no source) and
    /// `sc` only the light sourced by the atom's current. The light rule, the sponge and the deposit are all linear,
    /// so f + sc is the single-box field; the atom feels f + sc. `sc` IS the light the atom scatters and emits.
    pub sc: Option<Fotone3D>,
}

impl LuceAtomo {
    pub fn new() -> Self {
        let n = N_LUCE;
        let f = Fotone3D::new(n);
        let mut maschera = vec![1.0f32; n * n * n];
        maschera.par_iter_mut().enumerate().for_each(|(s, m)| {
            let c = [s / (n * n), (s / n) % n, s % n];
            let mut d2: f64 = 0.0;
            // a smooth cubic profile: one crossing of the sponge attenuates by ~e^{−7.5} (the light must not come
            // back around the periodic box: the atom's own light is ~10⁻³ of the pulse)
            for q in c { let d = (SPUGNA as f64 - q.min(n - 1 - q) as f64).max(0.0) / SPUGNA as f64; d2 += d * d * d; }
            *m = (-2.5 * d2).exp() as f32;
        });
        let mut pesi = vec![];
        let r = (3.0 * SIGMA_S).ceil() as i64;
        let h = (n / 2) as i64;
        for i in -r..=r { for j in -r..=r { for k in -r..=r {
            let d2 = (i * i + j * j + k * k) as f64;
            let w = (-d2 / (2.0 * SIGMA_S * SIGMA_S)).exp();
            if w < 1e-4 { continue; }
            pesi.push(((((h + i) as usize) * n + (h + j) as usize) * n + (h + k) as usize, w));
        }}}
        let tot: f64 = pesi.iter().map(|p| p.1).sum();
        for p in pesi.iter_mut() { p.1 /= tot; }
        LuceAtomo { f, maschera, pesi, sonda: vec![], assorbita: 0.0, serie_e: vec![], t: 0.0, trasv: None, sc: None }
    }

    /// The same box with the transverse deposit (patterns computed once per process, then shared).
    pub fn new_trasversa() -> Self {
        let mut l = Self::new();
        l.trasv = Some(Trasverso::condiviso(&l.pesi));
        l
    }

    /// The Gaussian weights (cell, weight) of the average around the atom.
    pub fn pesi_pub(&self) -> Vec<(usize, f64)> { self.pesi.clone() }

    /// Field at the atom made by the deposit of the source I of one light tick: −M I / b² (M = ⟨w|pattern⟩).
    pub fn campo_proprio(&self, corrente: [f64; 3]) -> [f64; 3] {
        let b2 = B_LUCE * B_LUCE;
        match &self.trasv {
            Some(t) => [0, 1, 2].map(|c| -(0..3).map(|q| t.m[c][q] * corrente[q]).sum::<f64>() / b2),
            None => {
                let w2: f64 = self.pesi.iter().map(|p| p.1 * p.1).sum();
                [0, 1, 2].map(|c| -w2 * corrente[c] / b2)
            }
        }
    }

    /// Split the box into free light (`f`) + scattered light (`sc`, starts empty). Call before any source is deposited.
    pub fn con_diffusa(&mut self) { if self.sc.is_none() { let mut g = Fotone3D::new(N_LUCE); g.ticks = self.f.ticks; self.sc = Some(g); } }

    /// The total field at cell s (f + sc when split).
    #[inline]
    pub fn tot(&self, s: usize) -> [crate::fotone3d::Cf; 3] {
        match &self.sc {
            None => self.f.f[s],
            Some(g) => { let (a, b) = (self.f.f[s], g.f[s]); [0, 1, 2].map(|c| crate::fotone3d::Cf { re: a[c].re + b[c].re, im: a[c].im + b[c].im }) }
        }
    }

    /// Free propagation of the light by one tick (no source, no sponge).
    pub fn propaga(&mut self) { self.f.tick(); if let Some(g) = self.sc.as_mut() { g.tick(); } }

    /// The source of one light tick: F −= I w / b² (Gaussian) or F −= Σ_c I_c T_c / b² (transverse).
    pub fn deposita(&mut self, corrente: [f64; 3]) {
        let b2 = B_LUCE * B_LUCE;
        // the source goes into the scattered box when the light is split
        let campo = match self.sc.as_mut() { Some(g) => &mut g.f, None => &mut self.f.f };
        match &self.trasv {
            None => {
                let b2 = b2 as f32;
                for &(s, w) in &self.pesi {
                    for c in 0..3 { campo[s][c].re -= (corrente[c] * w) as f32 / b2; }
                }
            }
            Some(t) => {
                let k = [0, 1, 2].map(|c| (corrente[c] / b2) as f32);
                campo.par_iter_mut().zip(t.t.par_iter()).for_each(|(p, tt)| {
                    for q in 0..3 { p[q].re -= k[0] * tt[0][q] + k[1] * tt[1][q] + k[2] * tt[2][q]; }
                });
            }
        }
    }

    /// Sponge, energy bookkeeping, probe: the end of a light tick.
    pub fn chiudi_tick(&mut self) {
        let m = &self.maschera;
        let (tolta, resta): (f64, f64) = match self.sc.as_mut() {
            None => self.f.f.par_iter_mut().zip(m.par_iter()).map(|(p, &w)| {
                let e = p.iter().map(|c| c.norm2() as f64).sum::<f64>();
                if w < 1.0 { for c in p.iter_mut() { c.re *= w; c.im *= w; } }
                let w2 = (w as f64) * (w as f64);
                ((1.0 - w2) * e, w2 * e)
            }).reduce(|| (0.0, 0.0), |x, y| (x.0 + y.0, x.1 + y.1)),
            // split: the sponge multiplies both boxes (linear); the energy is that of the total field f + sc
            Some(g) => self.f.f.par_iter_mut().zip(g.f.par_iter_mut()).zip(m.par_iter()).map(|((p, q), &w)| {
                let e = (0..3).map(|c| { let (re, im) = ((p[c].re + q[c].re) as f64, (p[c].im + q[c].im) as f64); re * re + im * im }).sum::<f64>();
                if w < 1.0 { for c in p.iter_mut().chain(q.iter_mut()) { c.re *= w; c.im *= w; } }
                let w2 = (w as f64) * (w as f64);
                ((1.0 - w2) * e, w2 * e)
            }).reduce(|| (0.0, 0.0), |x, y| (x.0 + y.0, x.1 + y.1)),
        };
        let v = 0.5 * B_LUCE.powi(3);
        self.assorbita += v * tolta;
        self.t += B_LUCE;
        self.serie_e.push((self.t, v * resta + self.assorbita));
        let n = N_LUCE;
        let s = (((n / 2 - SONDA) * n + n / 2) * n) + n / 2;
        self.sonda.push((self.t, self.tot(s)[2].re as f64));
    }

    /// Fourier factor of the Gaussian average at wave number k (λ_C⁻¹): the radiated power is reduced by W².
    pub fn fattore_forma(k: f64) -> f64 { (-0.5 * (k * SIGMA_S * B_LUCE).powi(2)).exp() }

    /// Electric field of the light at the atom (Gaussian average).
    pub fn campo_atomo(&self) -> [f64; 3] {
        let mut e = [0.0; 3];
        for &(s, w) in &self.pesi { let p = self.tot(s); for c in 0..3 { e[c] += w * p[c].re as f64; } }
        e
    }

    /// Electric field of the FREE light only at the atom (the pulse as it would be without the atom; = campo_atomo
    /// when the box is not split).
    pub fn campo_libero_atomo(&self) -> [f64; 3] {
        let mut e = [0.0; 3];
        for &(s, w) in &self.pesi { for c in 0..3 { e[c] += w * self.f.f[s][c].re as f64; } }
        e
    }

    /// One light tick with the source current I (charge × velocity, summed over the atom, averaged over the tick):
    /// propagate, deposit, sponge (the original order).
    pub fn passo(&mut self, corrente: [f64; 3]) {
        self.propaga();
        self.deposita(corrente);
        self.chiudi_tick();
    }

    /// A photon pulse of frequency ω (mc²/ħ) arriving along +x, polarised along z, with e·E ≈ `e0` at its peak.
    pub fn impulso(&mut self, omega: f64, e0: f64) {
        let n = N_LUCE;
        let lambda = std::f64::consts::TAU / omega / B_LUCE; // cells
        let sigma = 1.0 * lambda;
        let l = Lancio { lambda, w0: 1.2 * lambda, sigma, dir: [1.0, 0.0, 0.0], polar: Polar::Verticale,
                         t_fuoco: (n as f64 / 2.0 - 2.7 * sigma - 2.0).max(4.0) };
        let mut g = Fotone3D::new(n);
        g.lancia(&l);
        // linear light: |F| = √2 |E| at a crest
        let fmax = g.f.par_iter().map(|p| (p[0].norm2() + p[1].norm2() + p[2].norm2()).sqrt()).reduce(|| 0.0, f32::max).max(1e-30);
        let k = (std::f64::consts::SQRT_2 * e0 / e_eff()) as f32 / fmax;
        self.f.f.par_iter_mut().zip(g.f.par_iter()).for_each(|(a, b)| for c in 0..3 { a[c].re += k * b[c].re; a[c].im += k * b[c].im; });
    }

    /// Energy of the light in the box, ½ Σ |F|² b³ (units mc²/(e²)·…: only ratios are used).
    pub fn energia(&self) -> f64 {
        match &self.sc {
            None => 0.5 * B_LUCE.powi(3) * self.f.f.par_iter().map(|p| p.iter().map(|c| c.norm2() as f64).sum::<f64>()).sum::<f64>(),
            Some(_) => 0.5 * B_LUCE.powi(3) * (0..self.f.f.len()).into_par_iter().map(|s| self.tot(s).iter().map(|c| c.norm2() as f64).sum::<f64>()).sum::<f64>(),
        }
    }

    /// Energy of the scattered light alone (0 when not split).
    pub fn energia_diffusa(&self) -> f64 {
        self.sc.as_ref().map(|g| 0.5 * B_LUCE.powi(3) * g.f.par_iter().map(|p| p.iter().map(|c| c.norm2() as f64).sum::<f64>()).sum::<f64>()).unwrap_or(0.0)
    }
}

/// max |E| (real part of F) of a light box, over all cells or only farther than 8 light cells from the atom.
pub fn e_max(f: &[[crate::fotone3d::Cf; 3]], lontano: bool) -> f64 {
    let n = N_LUCE;
    let h = (n / 2) as i64;
    f.par_iter().enumerate().map(|(s, p)| {
        if lontano {
            let d = [(s / (n * n)) as i64 - h, ((s / n) % n) as i64 - h, (s % n) as i64 - h];
            if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= 64 { return 0.0; }
        }
        ((p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re) as f64).sqrt()
    }).reduce(|| 0.0, f64::max)
}

/// Frequency of a signal from its zero crossings after time `da` (linear interpolation): π per half period.
pub fn frequenza(serie: &[(f64, f64)], da: f64) -> Option<f64> {
    let v: Vec<&(f64, f64)> = serie.iter().filter(|p| p.0 >= da).collect();
    let mut zeri = vec![];
    for w in v.windows(2) {
        let (a, b) = (w[0], w[1]);
        if (a.1 <= 0.0 && b.1 > 0.0) || (a.1 >= 0.0 && b.1 < 0.0) {
            zeri.push(a.0 + (b.0 - a.0) * a.1 / (a.1 - b.1));
        }
    }
    if zeri.len() < 4 { return None; }
    Some(std::f64::consts::PI * (zeri.len() - 1) as f64 / (zeri[zeri.len() - 1] - zeri[0]))
}

/// Semi-classical decay: ln(p/(1−p)) falls with slope −Γ (dp/dt = −Γ p (1 − p)); fit after time `da`.
pub fn decadimento(serie: &[(f64, f64)], da: f64) -> Option<f64> {
    let p: Vec<(f64, f64)> = serie.iter().filter(|q| q.0 >= da && q.1 > 1e-4 && q.1 < 1.0 - 1e-4).map(|q| (q.0, (q.1 / (1.0 - q.1)).ln())).collect();
    if p.len() < 10 { return None; }
    let (a, b) = crate::fotscena::retta(&p);
    let _ = a;
    Some(-b)
}

/// Larmor prediction of Γ (same K, same Gaussian average): Γ = e² ω³ |z₁₂|² W² / (3π), with |z₁₂|² from the
/// measured dipole amplitude D and population p: D = 2 √(p(1−p)) |z₁₂|.
pub fn gamma_larmor(omega: f64, d_ampiezza: f64, p: f64) -> f64 {
    let z2 = d_ampiezza * d_ampiezza / (4.0 * p * (1.0 - p)).max(1e-12);
    let w2 = LuceAtomo::fattore_forma(omega).powi(2);
    e_eff().powi(2) * omega.powi(3) * z2 * w2 / (3.0 * std::f64::consts::PI)
}

/// The coupled experiment: one atom, its light, the measured series. Used by the scene and by the tests.
pub struct Esperimento {
    pub at: crate::zz::Atomo,
    pub luce: LuceAtomo,
    /// the 1s of the 3-D rule (listened) and its frequency; the 2p½ (m = ½) listened, and its frequency
    pub phi1s: Vec<[crate::cplx::C; 4]>,
    pub omega1s: f64,
    pub phi2p: Vec<[crate::cplx::C; 4]>,
    pub omega2p: f64,
    /// (t, excited population 1 − |⟨1s|ψ⟩|²) every light tick, and (t, ⟨z⟩) every atom tick
    pub serie_p: Vec<(f64, f64)>,
    pub serie_d: Vec<(f64, f64)>,
    /// (t, ⟨α_z⟩) every atom tick: the current of the walk (it should be d⟨z⟩/dt)
    pub serie_j: Vec<(f64, f64)>,
    /// time the last pulse was launched (λ_C/c)
    pub t_impulso: Option<f64>,
    /// ⟨r⟩ at the start of the current light tick
    r0: Option<[f64; 3]>,
    sub: usize,
    /// how field and current are exchanged (`ACCOPPIAMENTO` by default)
    pub acc: Accoppiamento,
    /// sources of the past light ticks (for the predictor), the predicted one and the field of the current tick
    jstoria: Vec<[f64; 3]>,
    jpred: [f64; 3],
    campo_tick: [f64; 3],
    /// coupling boost actually used by `passo` (K_ACC by default; 1 = the true charge)
    pub k_acc: f64,
    /// the light box split into free + scattered (docs/ATOMO.md §3.8); kept by `ricomincia`
    pub diffusa: bool,
}

impl Esperimento {
    /// Prepare: the start atom, then its clean 1s found by filter diagonalisation over 2 × `ascolto` ticks
    /// (docs/BATTIMENTO_1S.md), the 2p½ listened, the light box empty.
    pub fn new(ascolto: usize) -> Self { Self::nuovo(ascolto, None) }

    /// The same on one colour of the cube (`Some`, docs/SOTTORETICOLO.md: 4 times less work, the colours never
    /// exchange amplitude) or on the whole cube (`None`).
    pub fn nuovo(ascolto: usize, col: crate::zz::Colore) -> Self { Self::nuovo_prep(ascolto, col, crate::zz::NUCLEO, VALLE_ZERO && col.is_none()) }

    /// The same with the nucleus `nu` (docs/NUCLEO_1S.md §13–14).
    pub fn nuovo_nucleo(ascolto: usize, col: crate::zz::Colore, nu: crate::zz::Nucleo) -> Self { Self::nuovo_prep(ascolto, col, nu, false) }

    /// The same, and with `valle_zero` the 1s (smooth, learned nucleus) kept only in the electron's valley k = 0
    /// (docs/SURFISTI_ESATTI.md §13: then the surfers turn around the nucleus). Whole cube only.
    pub fn nuovo_prep(ascolto: usize, col: crate::zz::Colore, nu: crate::zz::Nucleo, valle_zero: bool) -> Self {
        let mut at = crate::zz::Atomo::new_in(col);
        at.con_nucleo(nu);
        // the 1s of the 3-D rule: one EXACT level of the 1s multiplet (the 4 valleys of the cubic walk mixed by the
        // nucleus), the one least touched by the first cells. Before (ascolta 0.885–0.915, short window, every
        // tick) the "1s" was a mixture of three levels and p_exc beat with period ~476 even without light.
        let (phi, om) = if nu == crate::zz::Nucleo::Nudo { let (s1, _) = at.uno_s_esatto(2 * ascolto); (s1.v, s1.omega) } else { let (v, w, _) = at.uno_s_liscio(2 * ascolto); (v, w) };
        let phi = if valle_zero { assert!(col.is_none()); crate::zz::proietta_valle_zero(&phi) } else { phi };
        // the 2p½: the radial 2p½ put on the lattice is already (almost) one level; listened briefly
        let mut a2 = crate::zz::Atomo::con_stato_in(1, 2, col);
        a2.con_nucleo(nu);
        let (phi2, om2) = a2.ascolta(0.96, 0.99, ascolto / 3);
        at.imposta(phi.clone());
        let mut x = Esperimento { at, luce: LuceAtomo::new(), phi1s: phi, omega1s: om, phi2p: phi2, omega2p: om2, serie_p: vec![], serie_d: vec![], serie_j: vec![],
                      t_impulso: None, r0: None, sub: 0, acc: ACCOPPIAMENTO, jstoria: vec![], jpred: [0.0; 3], campo_tick: [0.0; 3], k_acc: K_ACC, diffusa: false };
        x.ricomincia();
        x
    }

    /// From states already found (whole cube): the 1s and a 2p with their frequencies; no listening.
    pub fn da_stati(phi1s: Vec<[crate::cplx::C; 4]>, omega1s: f64, phi2p: Vec<[crate::cplx::C; 4]>, omega2p: f64, nu: crate::zz::Nucleo) -> Self {
        let mut at = crate::zz::Atomo::new_in(None);
        at.con_nucleo(nu);
        at.imposta(phi1s.clone());
        let mut x = Esperimento { at, luce: LuceAtomo::new(), phi1s, omega1s, phi2p, omega2p, serie_p: vec![], serie_d: vec![], serie_j: vec![],
                      t_impulso: None, r0: None, sub: 0, acc: ACCOPPIAMENTO, jstoria: vec![], jpred: [0.0; 3], campo_tick: [0.0; 3], k_acc: K_ACC, diffusa: false };
        x.ricomincia();
        x
    }

    /// Split the light into free + scattered (exact, the rule is linear) and start again from the 1s.
    pub fn con_diffusa(mut self) -> Self {
        self.diffusa = true;
        self.ricomincia();
        self
    }

    /// Change the coupling and start again from the 1s (empty light box).
    pub fn con_accoppiamento(mut self, acc: Accoppiamento) -> Self {
        self.acc = acc;
        self.ricomincia();
        self
    }

    /// The 1s → 2p level difference of the 3-D rule (listened).
    pub fn omega_21(&self) -> f64 { self.omega2p - self.omega1s }

    /// Population of the 2p½ (m = ½) of the 3-D rule.
    pub fn p2p(&self) -> f64 { self.at.proiezione(&self.phi2p) }

    /// Start again from the listened 1s with an empty light box.
    pub fn ricomincia(&mut self) {
        self.at.imposta(self.phi1s.clone());
        self.luce = if self.acc.trasverso { LuceAtomo::new_trasversa() } else { LuceAtomo::new() };
        if self.diffusa { self.luce.con_diffusa(); }
        self.jstoria.clear();
        self.jpred = [0.0; 3];
        self.campo_tick = [0.0; 3];
        self.serie_p.clear();
        self.serie_d.clear();
        self.serie_j.clear();
        self.t_impulso = None;
        self.r0 = None;
        self.sub = 0;
    }

    pub fn spara(&mut self, omega: f64, e0: f64) {
        self.luce.impulso(omega, e0);
        self.t_impulso = Some(self.at.t());
    }

    /// Start from an arbitrary atom state (same lattice), empty light box, clean history.
    pub fn parti_da(&mut self, psi: Vec<[crate::cplx::C; 4]>) {
        self.ricomincia();
        self.at.imposta(psi);
    }

    /// Predicted source of the coming light tick (linear prediction, as lab/emiss_run.py `sc_predict`).
    fn predici(&self) -> [f64; 3] {
        let h = &self.jstoria;
        let d = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        match h.len() {
            0 => [0.0; 3],
            1 => h[0],
            _ => {
                let (y, z) = (h[h.len() - 2], h[h.len() - 1]);
                if d(y, y) <= 0.0 { return z; }
                if h.len() == 2 { return [0, 1, 2].map(|c| 2.0 * z[c] - y[c]); }
                let x = h[h.len() - 3];
                let c = (d([x[0] + z[0], x[1] + z[1], x[2] + z[2]], y) / (2.0 * d(y, y))).clamp(-1.0, 1.0);
                [0, 1, 2].map(|q| 2.0 * c * z[q] - y[q])
            }
        }
    }

    /// One atom tick (and a light tick every SUB atom ticks).
    pub fn passo(&mut self) {
        let e = e_eff() * self.k_acc / K_ACC;
        if self.sub == 0 {
            self.r0 = Some(self.at.posizione());
            if self.acc.centrato {
                // F' = U F_n first; the atom sees ½(F' + F_{n+1}) = field(F') + ½ (self field of the predicted source)
                self.luce.propaga();
                self.jpred = self.predici();
                let ea = self.luce.campo_atomo();
                let sp = self.luce.campo_proprio(self.jpred);
                self.campo_tick = [0, 1, 2].map(|c| ea[c] + 0.5 * sp[c]);
            } else {
                self.campo_tick = self.luce.campo_atomo();
            }
        }
        let ea = self.campo_tick;
        self.at.campo_est = [e * ea[0], e * ea[1], e * ea[2]];
        self.at.avanti();
        self.sub += 1;
        let r = self.at.posizione();
        self.serie_d.push((self.at.t(), r[2]));
        if self.sub == SUB {
            // the current of the light tick: the charge that moved, −e Δ⟨r⟩/Δt (the walk's conserved flux of
            // balls between cells; the instantaneous R†σR − L†σL of one time slice is not it on this lattice)
            let r0 = self.r0.take().unwrap();
            let dt = SUB as f64 * crate::zz::A;
            let i = [-e * (r[0] - r0[0]) / dt, -e * (r[1] - r0[1]) / dt, -e * (r[2] - r0[2]) / dt];
            self.serie_j.push((self.at.t(), (r[2] - r0[2]) / dt));
            if self.acc.centrato {
                // corrector: the atom felt ½ selffield(I_pred) but the true source was I: a phase e ½ ΔE·r over the tick
                let sp = self.luce.campo_proprio([0, 1, 2].map(|c| i[c] - self.jpred[c]));
                let v = [0, 1, 2].map(|c| e * 0.5 * sp[c] * dt);
                kick(&mut self.at, v);
                self.luce.deposita(i);
            } else {
                self.luce.propaga();
                self.luce.deposita(i);
            }
            self.luce.chiudi_tick();
            self.jstoria.push(i);
            if self.jstoria.len() > 3 { self.jstoria.remove(0); }
            self.sub = 0;
            let p = 1.0 - self.at.proiezione(&self.phi1s);
            self.serie_p.push((self.at.t(), p));
        }
    }

    /// Amplitude of ⟨z⟩ oscillation (√2 × RMS around its mean) over the last `finestra` time units.
    pub fn ampiezza_dipolo(&self, finestra: f64) -> f64 {
        let t1 = self.at.t();
        let v: Vec<f64> = self.serie_d.iter().filter(|p| p.0 >= t1 - finestra).map(|p| p.1).collect();
        if v.len() < 4 { return 0.0; }
        let m = v.iter().sum::<f64>() / v.len() as f64;
        (2.0 * v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
    }
}

/// ψ ← e^{−i v·r} ψ (a uniform-field kick, v = e E Δt).
pub fn kick(at: &mut crate::zz::Atomo, v: [f64; 3]) {
    use crate::zz::{sito, x_of};
    let col = at.colore;
    at.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
        let c = sito(s, col);
        let ph = crate::cplx::C::expi(-(v[0] * x_of(c[0]) + v[1] * x_of(c[1]) + v[2] * x_of(c[2])));
        for q in p.iter_mut() { *q = *q * ph; }
    });
}

/// Amplitude of the component at frequency ω of a series after time `da` (Hann-weighted projection).
pub fn ampiezza_a(serie: &[(f64, f64)], omega: f64, da: f64, a: f64) -> f64 {
    let v: Vec<&(f64, f64)> = serie.iter().filter(|p| p.0 >= da && p.0 <= a).collect();
    if v.len() < 8 { return 0.0; }
    let (t0, t1) = (v[0].0, v[v.len() - 1].0);
    let (mut re, mut im, mut w) = (0.0, 0.0, 0.0);
    for p in &v {
        let h = (std::f64::consts::PI * (p.0 - t0) / (t1 - t0)).sin().powi(2);
        re += h * p.1 * (omega * p.0).cos();
        im += h * p.1 * (omega * p.0).sin();
        w += h;
    }
    2.0 * (re * re + im * im).sqrt() / w
}

/// Larmor power of an oscillating dipole −e·D cos ωt seen through the Gaussian average: e² ω⁴ D² W² / (12π).
pub fn potenza_larmor(omega: f64, d: f64) -> f64 {
    e_eff().powi(2) * omega.powi(4) * d * d * LuceAtomo::fattore_forma(omega).powi(2) / (12.0 * std::f64::consts::PI)
}

/// TRANSVERSE DEPOSIT (docs/EMISSIONE.md §3.2). The two-tick light walk at wave vector k is a real rotation R(k)
/// (one translation along j is exp(−k_j G_j), (G_j)_ab = ε_jab, a rotation by k_j about axis j; ticks alternate the
/// axis order x,y,z / z,y,x). Its axis n(k) is the lattice's longitudinal (static, Coulomb) mode exactly, so the light
/// part of a source is P_⊥(k) = 1 − n nᵀ applied to it. The Gaussian deposit w e_c is projected once by FFT:
/// pattern T_c = IFFT[P_⊥(k) e_c W(k)] (real and even, it fills the box: the longitudinal tail goes as 1/r³).
/// M_ab = Σ_s w_s T_b(s)_a is the instantaneous self field at the atom per unit source.
pub struct Trasverso {
    /// per light cell: t[s][c][q] = component q of the pattern of the source along c
    pub t: Vec<[[f32; 3]; 3]>,
    pub m: [[f64; 3]; 3],
    /// largest |imaginary part| / largest |real part| of the patterns (should be round-off)
    pub immaginaria: f64,
}

/// Two-tick rotation of the light rule at wave vector k (cells⁻¹), starting from an even tick.
pub fn rotazione_luce(k: [f64; 3]) -> [[f64; 3]; 3] {
    fn rot(j: usize, a: f64) -> [[f64; 3]; 3] {
        // exp(−a G_j) = 1 − sin a G_j + (1 − cos a) G_j²
        let mut g = [[0.0; 3]; 3];
        let (p, q) = ((j + 1) % 3, (j + 2) % 3);
        g[p][q] = 1.0;
        g[q][p] = -1.0;
        let mut g2 = [[0.0; 3]; 3];
        for i in 0..3 { for l in 0..3 { for m in 0..3 { g2[i][l] += g[i][m] * g[m][l]; } } }
        let (s, c) = a.sin_cos();
        let mut r = [[0.0; 3]; 3];
        for i in 0..3 { for l in 0..3 { r[i][l] = if i == l { 1.0 } else { 0.0 } - s * g[i][l] + (1.0 - c) * g2[i][l]; } }
        r
    }
    fn mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let mut r = [[0.0; 3]; 3];
        for i in 0..3 { for l in 0..3 { for m in 0..3 { r[i][l] += a[i][m] * b[m][l]; } } }
        r
    }
    let mut r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for t in 0..2u64 { for j in crate::walk::order(t) { r = mul(rot(j, k[j]), r); } }
    r
}

/// P_⊥(k) = 1 − n nᵀ of the light rule (2/3·1 where the rotation is the identity, e.g. k = 0).
pub fn proiettore_trasverso(k: [f64; 3]) -> [[f64; 3]; 3] {
    let r = rotazione_luce(k);
    let mut ax = [r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]];
    let mut nr = (ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]).sqrt();
    if nr < 1e-9 {
        let tr = r[0][0] + r[1][1] + r[2][2];
        if tr > 2.9 { return [[2.0 / 3.0, 0.0, 0.0], [0.0, 2.0 / 3.0, 0.0], [0.0, 0.0, 2.0 / 3.0]]; }
        // rotation by π: R = 2 n nᵀ − 1, take the largest column of (R + 1)/2
        let c = (0..3).max_by(|&a, &b| r[a][a].total_cmp(&r[b][b])).unwrap();
        ax = [0.5 * (r[0][c] + if c == 0 { 1.0 } else { 0.0 }), 0.5 * (r[1][c] + if c == 1 { 1.0 } else { 0.0 }), 0.5 * (r[2][c] + if c == 2 { 1.0 } else { 0.0 })];
        nr = (ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]).sqrt();
    }
    let n = [ax[0] / nr, ax[1] / nr, ax[2] / nr];
    let mut p = [[0.0; 3]; 3];
    for i in 0..3 { for l in 0..3 { p[i][l] = if i == l { 1.0 } else { 0.0 } - n[i] * n[l]; } }
    p
}

/// In-place radix-2 FFT of one line (sign −1 forward, +1 inverse, no normalisation).
fn fft_riga(v: &mut [crate::cplx::C], segno: f64) {
    use crate::cplx::C;
    let n = v.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 { j ^= bit; bit >>= 1; }
        j |= bit;
        if i < j { v.swap(i, j); }
    }
    let mut len = 2;
    while len <= n {
        let w = C::expi(segno * std::f64::consts::TAU / len as f64);
        for i in (0..n).step_by(len) {
            let mut wk = C::ONE;
            for k in 0..len / 2 {
                let (a, b) = (v[i + k], v[i + k + len / 2] * wk);
                v[i + k] = a + b;
                v[i + k + len / 2] = a - b;
                wk = wk * w;
            }
        }
        len <<= 1;
    }
}

/// 3-D FFT of an n³ array (C order), sign as `fft_riga`.
fn fft3(a: &mut [crate::cplx::C], n: usize, segno: f64) {
    use crate::cplx::C;
    a.par_chunks_mut(n).for_each(|r| fft_riga(r, segno));
    for ax in 0..2 {
        // lines along x (ax 0, stride n²) or y (ax 1, stride n)
        let stride = if ax == 0 { n * n } else { n };
        let starts: Vec<usize> = (0..n * n).map(|q| if ax == 0 { q } else { (q / n) * n * n + q % n }).collect();
        let lines: Vec<Vec<C>> = starts.par_iter().map(|&s0| {
            let mut l: Vec<C> = (0..n).map(|i| a[s0 + i * stride]).collect();
            fft_riga(&mut l, segno);
            l
        }).collect();
        for (s0, l) in starts.iter().zip(lines) { for i in 0..n { a[s0 + i * stride] = l[i]; } }
    }
}

impl Trasverso {
    pub fn calcola(pesi: &[(usize, f64)]) -> Trasverso {
        use crate::cplx::C;
        let n = N_LUCE;
        let nn = n * n * n;
        let mut w = vec![C::ZERO; nn];
        for &(s, p) in pesi { w[s] = C::real(p); }
        fft3(&mut w, n, -1.0);
        let kk = |i: usize| std::f64::consts::TAU * (if i < n / 2 { i as f64 } else { i as f64 - n as f64 });
        let p: Vec<[[f64; 3]; 3]> = (0..nn).into_par_iter().map(|s| proiettore_trasverso([kk(s / (n * n)) / n as f64 * 1.0, kk((s / n) % n) / n as f64, kk(s % n) / n as f64])).collect();
        let mut t = vec![[[0.0f32; 3]; 3]; nn];
        let mut imm: f64 = 0.0;
        let mut re_max: f64 = 0.0;
        for c in 0..3 { for q in 0..3 {
            let mut a: Vec<C> = w.par_iter().zip(p.par_iter()).map(|(wk, pk)| wk.scale(pk[q][c])).collect();
            fft3(&mut a, n, 1.0);
            let k = 1.0 / nn as f64;
            for s in 0..nn {
                t[s][c][q] = (a[s].re * k) as f32;
                imm = imm.max((a[s].im * k).abs());
                re_max = re_max.max((a[s].re * k).abs());
            }
        }}
        let mut m = [[0.0; 3]; 3];
        for a in 0..3 { for b in 0..3 { m[a][b] = pesi.iter().map(|&(s, ww)| ww * t[s][b][a] as f64).sum(); } }
        Trasverso { t, m, immaginaria: imm / re_max.max(1e-300) }
    }

    /// Computed once per process (the Gaussian weights depend only on the constants of this file).
    pub fn condiviso(pesi: &[(usize, f64)]) -> std::sync::Arc<Trasverso> {
        static T: std::sync::OnceLock<std::sync::Arc<Trasverso>> = std::sync::OnceLock::new();
        T.get_or_init(|| std::sync::Arc::new(Trasverso::calcola(pesi))).clone()
    }
}

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
        LuceAtomo { f, maschera, pesi, sonda: vec![], assorbita: 0.0, serie_e: vec![], t: 0.0 }
    }

    /// Fourier factor of the Gaussian average at wave number k (λ_C⁻¹): the radiated power is reduced by W².
    pub fn fattore_forma(k: f64) -> f64 { (-0.5 * (k * SIGMA_S * B_LUCE).powi(2)).exp() }

    /// Electric field of the light at the atom (Gaussian average).
    pub fn campo_atomo(&self) -> [f64; 3] {
        let mut e = [0.0; 3];
        for &(s, w) in &self.pesi { for c in 0..3 { e[c] += w * self.f.f[s][c].re as f64; } }
        e
    }

    /// One light tick with the source current I (charge × velocity, summed over the atom, averaged over the tick).
    pub fn passo(&mut self, corrente: [f64; 3]) {
        self.f.tick();
        let b2 = (B_LUCE * B_LUCE) as f32;
        for &(s, w) in &self.pesi {
            for c in 0..3 { self.f.f[s][c].re -= (corrente[c] * w) as f32 / b2; }
        }
        let m = &self.maschera;
        let (tolta, resta): (f64, f64) = self.f.f.par_iter_mut().zip(m.par_iter()).map(|(p, &w)| {
            let e = p.iter().map(|c| c.norm2() as f64).sum::<f64>();
            if w < 1.0 { for c in p.iter_mut() { c.re *= w; c.im *= w; } }
            let w2 = (w as f64) * (w as f64);
            ((1.0 - w2) * e, w2 * e)
        }).reduce(|| (0.0, 0.0), |x, y| (x.0 + y.0, x.1 + y.1));
        let v = 0.5 * B_LUCE.powi(3);
        self.assorbita += v * tolta;
        self.t += B_LUCE;
        self.serie_e.push((self.t, v * resta + self.assorbita));
        let n = N_LUCE;
        let s = (((n / 2 - SONDA) * n + n / 2) * n) + n / 2;
        self.sonda.push((self.t, self.f.f[s][2].re as f64));
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
        0.5 * B_LUCE.powi(3) * self.f.f.par_iter().map(|p| p.iter().map(|c| c.norm2() as f64).sum::<f64>()).sum::<f64>()
    }
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
}

impl Esperimento {
    /// Prepare: the start atom, then its clean 1s found by filter diagonalisation over 2 × `ascolto` ticks
    /// (docs/BATTIMENTO_1S.md), the 2p½ listened, the light box empty.
    pub fn new(ascolto: usize) -> Self { Self::nuovo(ascolto, None) }

    /// The same on one colour of the cube (`Some`, docs/SOTTORETICOLO.md: 4 times less work, the colours never
    /// exchange amplitude) or on the whole cube (`None`).
    pub fn nuovo(ascolto: usize, col: crate::zz::Colore) -> Self {
        let mut at = crate::zz::Atomo::new_in(col);
        // the 1s of the 3-D rule: one EXACT level of the 1s multiplet (the 4 valleys of the cubic walk mixed by the
        // nucleus), the one least touched by the first cells. Before (ascolta 0.885–0.915, short window, every
        // tick) the "1s" was a mixture of three levels and p_exc beat with period ~476 even without light.
        let (s1, _) = at.uno_s_esatto(2 * ascolto);
        let (phi, om) = (s1.v, s1.omega);
        // the 2p½: the radial 2p½ put on the lattice is already (almost) one level; listened briefly
        let a2 = crate::zz::Atomo::con_stato_in(1, 2, col);
        let (phi2, om2) = a2.ascolta(0.96, 0.99, ascolto / 3);
        at.imposta(phi.clone());
        Esperimento { at, luce: LuceAtomo::new(), phi1s: phi, omega1s: om, phi2p: phi2, omega2p: om2, serie_p: vec![], serie_d: vec![], serie_j: vec![],
                      t_impulso: None, r0: None, sub: 0 }
    }

    /// The 1s → 2p level difference of the 3-D rule (listened).
    pub fn omega_21(&self) -> f64 { self.omega2p - self.omega1s }

    /// Population of the 2p½ (m = ½) of the 3-D rule.
    pub fn p2p(&self) -> f64 { self.at.proiezione(&self.phi2p) }

    /// Start again from the listened 1s with an empty light box.
    pub fn ricomincia(&mut self) {
        self.at.imposta(self.phi1s.clone());
        self.luce = LuceAtomo::new();
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

    /// One atom tick (and a light tick every SUB atom ticks).
    pub fn passo(&mut self) {
        let e = e_eff();
        let ea = self.luce.campo_atomo();
        self.at.campo_est = [e * ea[0], e * ea[1], e * ea[2]];
        if self.r0.is_none() { self.r0 = Some(self.at.posizione()); }
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
            self.luce.passo(i);
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

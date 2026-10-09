//! Physics of one photon–electron collision (Compton scattering on a FREE electron at rest).
//!
//! What is derived and what is an INPUT:
//! • Kinematics (energy and momentum conservation, relativistic electron): derived, exact.
//!   λ' − λ = (h/mₑc)(1 − cos θ),  E' = E / (1 + k(1 − cos θ)),  k = E/mₑc².
//! • How much of the wave scatters and into which direction: the Klein–Nishina cross-section
//!   of QED (1929). This is an INPUT rule from QED, NOT derived from our balls.
//! • The surfers: the photon surfer is born ∝ |ψ|² inside the incoming packet; it ends in the
//!   scattered part iff its transverse distance b from the electron surfer lies inside the disk
//!   of area σ_tot around it (rank in the incoming flux below the scattered fraction), and then
//!   keeps its rank in the flux: the rank of b picks θ (smallest b → most backward), the azimuth
//!   of b picks φ. This map preserves |ψ|² (equivariance), so the surfers reproduce the
//!   Klein–Nishina statistics; the true Bohmian map of the full 6-D wave would be another
//!   measure-preserving map (same statistics, different single trajectories).
//!
//! Units: energies keV, momenta keV/c, lengths pm (unless stated), cross-sections in r_e² or m².

use std::f64::consts::PI;

/// mₑc² (CODATA 2018), keV.
pub const ME_C2_KEV: f64 = 510.998_950_00;
/// hc (CODATA 2018, exact h and c), keV·pm.
pub const HC_KEV_PM: f64 = 1239.841_984_332;
/// Compton wavelength of the electron h/mₑc (CODATA 2018), pm.
pub const LAMBDA_C_PM: f64 = 2.426_310_238_67;
/// Classical electron radius (CODATA 2018), m.
pub const R_E_M: f64 = 2.817_940_326_2e-15;
/// ħ (exact), J·s, and mₑ (CODATA 2018), kg: only for the spreading time of the electron packet.
pub const HBAR: f64 = 1.054_571_817e-34;
pub const ME_KG: f64 = 9.109_383_7015e-31;
pub const C_M_S: f64 = 299_792_458.0;

/// Thomson cross-section 8π/3 r_e², m².
pub fn sigma_thomson() -> f64 { 8.0 * PI / 3.0 * R_E_M * R_E_M }

pub fn lambda_pm(e_kev: f64) -> f64 { HC_KEV_PM / e_kev }

/// Scattered photon energy (keV) at polar angle θ.
pub fn e_out(e: f64, theta: f64) -> f64 { e / (1.0 + e / ME_C2_KEV * (1.0 - theta.cos())) }

/// Compton's formula: scattered wavelength (pm) for incoming wavelength `lam` (pm).
pub fn lam_out(lam: f64, theta: f64) -> f64 { lam + LAMBDA_C_PM * (1.0 - theta.cos()) }

/// Electron recoil angle from the formula tan φ = cot(θ/2)/(1 + E/mc²) (radians, 0 … π/2).
pub fn electron_angle_formula(e: f64, theta: f64) -> f64 {
    let half = 0.5 * theta;
    (half.cos()).atan2(half.sin() * (1.0 + e / ME_C2_KEV))
}

/// Compton edge: maximum electron kinetic energy (θ = 180°), keV.
pub fn compton_edge(e: f64) -> f64 { let k = e / ME_C2_KEV; e * 2.0 * k / (1.0 + 2.0 * k) }

/// Unit vector of direction (θ from +z, φ from +x).
pub fn dir(theta: f64, phi: f64) -> [f64; 3] {
    let (st, ct) = theta.sin_cos();
    let (sp, cp) = phi.sin_cos();
    [st * cp, st * sp, ct]
}

/// One scattering event: the correlated pair (photon at θ, φ; electron at the conjugate angle).
#[derive(Clone, Copy, Debug)]
pub struct Evento {
    pub theta: f64,
    pub phi: f64,
    /// scattered photon energy (keV) and wavelength (pm)
    pub e_out: f64,
    pub lam_out: f64,
    /// momenta (keV/c) of the scattered photon and of the electron
    pub pg: [f64; 3],
    pub pe: [f64; 3],
    /// electron kinetic energy (keV), polar angle (rad), speed β = v/c
    pub t_e: f64,
    pub theta_e: f64,
    pub beta_e: f64,
}

/// Energy–momentum conservation for an incoming photon of energy `e` along +z hitting an electron
/// at rest, photon leaving at (θ, φ). The electron takes exactly the missing momentum.
pub fn cinematica(e: f64, theta: f64, phi: f64) -> Evento {
    let e2 = e_out(e, theta);
    let n = dir(theta, phi);
    let pg = [e2 * n[0], e2 * n[1], e2 * n[2]];
    let pe = [-pg[0], -pg[1], e - pg[2]];
    let p2 = pe[0] * pe[0] + pe[1] * pe[1] + pe[2] * pe[2];
    let p = p2.sqrt();
    let etot = (p2 + ME_C2_KEV * ME_C2_KEV).sqrt();
    // kinetic energy from the momentum, without cancellation: p²/(E + m)
    let t_e = p2 / (etot + ME_C2_KEV);
    let theta_e = (pe[0].hypot(pe[1])).atan2(pe[2]);
    Evento { theta, phi, e_out: e2, lam_out: lam_out(lambda_pm(e), theta), pg, pe, t_e, theta_e, beta_e: p / etot }
}

/// Polarization of the incoming photon.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pol {
    /// linear, electric field along x
    LineareX,
    /// unpolarized (average over the two linear polarizations)
    NonPolarizzata,
}

/// Klein–Nishina dσ/dΩ in units of r_e². φ is measured from the polarization direction (x);
/// outgoing polarizations are summed (the detector does not look at them).
pub fn kn_domega(k: f64, theta: f64, phi: f64, pol: Pol) -> f64 {
    let p = 1.0 / (1.0 + k * (1.0 - theta.cos()));
    let s2 = theta.sin().powi(2);
    match pol {
        Pol::NonPolarizzata => 0.5 * p * p * (p + 1.0 / p - s2),
        Pol::LineareX => 0.5 * p * p * (p + 1.0 / p - 2.0 * s2 * phi.cos().powi(2)),
    }
}

/// dσ/dμ (μ = cos θ), integrated over φ, in r_e².
pub fn kn_dmu(k: f64, mu: f64) -> f64 {
    let p = 1.0 / (1.0 + k * (1.0 - mu));
    PI * p * p * (p + 1.0 / p - (1.0 - mu * mu))
}

/// Total Klein–Nishina cross-section in r_e² (closed form; for k < 0.005 the series, which avoids
/// the cancellation of the closed form: σ/σ_T = 1 − 2k + 26k²/5 − 133k³/10 + 1144k⁴/35 …).
pub fn kn_total_formula(k: f64) -> f64 {
    let st = 8.0 * PI / 3.0;
    if k < 0.005 {
        return st * (1.0 - 2.0 * k + 5.2 * k * k - 13.3 * k.powi(3) + 1144.0 / 35.0 * k.powi(4));
    }
    let l = (1.0 + 2.0 * k).ln();
    2.0 * PI * ((1.0 + k) / (k * k) * (2.0 * (1.0 + k) / (1.0 + 2.0 * k) - l / k) + l / (2.0 * k)
        - (1.0 + 3.0 * k) / (1.0 + 2.0 * k).powi(2))
}

// 8-point Gauss–Legendre on [-1, 1]
const GL_X: [f64; 8] = [
    -0.960_289_856_497_536_2, -0.796_666_477_413_626_7, -0.525_532_409_916_329_0, -0.183_434_642_495_649_8,
    0.183_434_642_495_649_8, 0.525_532_409_916_329_0, 0.796_666_477_413_626_7, 0.960_289_856_497_536_2,
];
const GL_W: [f64; 8] = [
    0.101_228_536_290_376_3, 0.222_381_034_453_374_5, 0.313_706_645_877_887_3, 0.362_683_783_378_362_0,
    0.362_683_783_378_362_0, 0.313_706_645_877_887_3, 0.222_381_034_453_374_5, 0.101_228_536_290_376_3,
];

fn gl(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let h = 0.5 * (b - a);
    let m = 0.5 * (a + b);
    (0..8).map(|i| GL_W[i] * f(m + h * GL_X[i])).sum::<f64>() * h
}

const N_TAB: usize = 2048;

/// The Klein–Nishina table of one photon energy: cumulative cross-section from θ = 180°
/// (μ = −1) forward, used by the surfer rule to turn a rank into an angle.
#[derive(Clone, Debug)]
pub struct Kn {
    pub e_kev: f64,
    pub k: f64,
    pub lam_pm: f64,
    /// cumulative ∫_{-1}^{μ_i} dσ/dμ dμ (r_e²) at μ_i = −1 + 2i/N
    cum: Vec<f64>,
    /// total cross-section, r_e² (numerical integral of the table)
    pub sigma_re2: f64,
}

impl Kn {
    pub fn new(e_kev: f64) -> Self {
        let k = e_kev / ME_C2_KEV;
        let mut cum = Vec::with_capacity(N_TAB + 1);
        cum.push(0.0);
        let mut acc = 0.0;
        for i in 0..N_TAB {
            let a = -1.0 + 2.0 * i as f64 / N_TAB as f64;
            let b = -1.0 + 2.0 * (i + 1) as f64 / N_TAB as f64;
            acc += gl(|m| kn_dmu(k, m), a, b);
            cum.push(acc);
        }
        Kn { e_kev, k, lam_pm: lambda_pm(e_kev), cum, sigma_re2: acc }
    }

    /// Total cross-section, m².
    pub fn sigma_m2(&self) -> f64 { self.sigma_re2 * R_E_M * R_E_M }

    /// Fraction of the cross-section with cos θ < μ.
    pub fn cdf_mu(&self, mu: f64) -> f64 {
        let mu = mu.clamp(-1.0, 1.0);
        let x = (mu + 1.0) * 0.5 * N_TAB as f64;
        let i = (x.floor() as usize).min(N_TAB - 1);
        let a = -1.0 + 2.0 * i as f64 / N_TAB as f64;
        (self.cum[i] + gl(|m| kn_dmu(self.k, m), a, mu)) / self.sigma_re2
    }

    /// μ = cos θ with cdf_mu(μ) = v (v in [0,1]): inverse of the table (bisection + Newton).
    pub fn mu_of_rank(&self, v: f64) -> f64 {
        let y = v.clamp(0.0, 1.0) * self.sigma_re2;
        let i = self.cum.partition_point(|&c| c <= y).saturating_sub(1).min(N_TAB - 1);
        let (mut lo, mut hi) = (-1.0 + 2.0 * i as f64 / N_TAB as f64, -1.0 + 2.0 * (i + 1) as f64 / N_TAB as f64);
        let target = y - self.cum[i];
        let a = lo;
        let mut mu = 0.5 * (lo + hi);
        for _ in 0..60 {
            let g = gl(|m| kn_dmu(self.k, m), a, mu) - target;
            if g.abs() < 1e-17 * self.sigma_re2 { break; }
            if g > 0.0 { hi = mu } else { lo = mu }
            let d = kn_dmu(self.k, mu);
            let next = mu - g / d;
            mu = if next >= lo && next <= hi { next } else { 0.5 * (lo + hi) };
            if hi - lo < 1e-15 { break; }
        }
        mu
    }

    /// Azimuth φ ∈ [0, 2π) whose conditional cumulative (at polar angle θ) equals w.
    /// Density ∝ A − B cos²φ with A = P + 1/P, B = 2 sin²θ (linear polarization along x).
    pub fn phi_of_rank(&self, theta: f64, w: f64, pol: Pol) -> f64 {
        let w = w.clamp(0.0, 1.0);
        if pol == Pol::NonPolarizzata { return 2.0 * PI * w; }
        let p = 1.0 / (1.0 + self.k * (1.0 - theta.cos()));
        let a = p + 1.0 / p;
        let b = 2.0 * theta.sin().powi(2);
        let norm = 2.0 * PI * (a - 0.5 * b);
        let g = |f: f64| (a * f - b * (0.5 * f + 0.25 * (2.0 * f).sin())) / norm;
        let (mut lo, mut hi) = (0.0, 2.0 * PI);
        let mut f = 2.0 * PI * w;
        for _ in 0..80 {
            let r = g(f) - w;
            if r.abs() < 1e-16 { break; }
            if r > 0.0 { hi = f } else { lo = f }
            let d = (a - b * f.cos().powi(2)) / norm;
            let next = if d > 1e-300 { f - r / d } else { 0.5 * (lo + hi) };
            f = if next >= lo && next <= hi { next } else { 0.5 * (lo + hi) };
            if hi - lo < 1e-15 { break; }
        }
        f
    }
}

/// The incoming joint wave, as far as the surfers need it: transverse rms widths of the photon
/// packet and of the electron packet (both INPUTS, stated in the panel), in metres.
#[derive(Clone, Copy, Debug)]
pub struct Pacchetti {
    pub s_gamma_m: f64,
    pub s_e_m: f64,
}

impl Pacchetti {
    /// rms width per axis of the relative transverse coordinate b = x_γ − x_e (m).
    pub fn s_rel(&self) -> f64 { self.s_gamma_m.hypot(self.s_e_m) }
    /// Probability that ONE photon of this packet scatters on the electron (first order):
    /// σ × (probability per unit area of b at b = 0) = σ / (2π s²).
    pub fn p_scatter(&self, kn: &Kn) -> f64 { kn.sigma_m2() / (2.0 * PI * self.s_rel().powi(2)) }
}

/// What the two surfers do when the photon surfer reaches the electron surfer.
#[derive(Clone, Copy, Debug)]
pub enum Esito {
    /// photon passes straight, electron stays
    Passa,
    /// correlated pair
    Diffuso(Evento),
}

/// THE SURFER RULE. `b` = transverse position (m) of the photon surfer relative to the electron
/// surfer when they meet. Rank in the incoming flux: u = F(b) = 1 − exp(−b²/2s²) (fraction of the
/// relative probability inside radius b). If u < p (scattered fraction) the pair scatters:
/// rank v = u/p → θ (v = 0 → θ = 180°), azimuth of b → rank w → φ.
pub fn regola(kn: &Kn, pk: &Pacchetti, b: [f64; 2], pol: Pol) -> Esito {
    let s = pk.s_rel();
    let u = -(-(b[0] * b[0] + b[1] * b[1]) / (2.0 * s * s)).exp_m1();
    let p = pk.p_scatter(kn);
    if u >= p { return Esito::Passa; }
    let v = u / p;
    let mu = kn.mu_of_rank(v);
    let theta = mu.clamp(-1.0, 1.0).acos();
    let alpha = b[1].atan2(b[0]).rem_euclid(2.0 * PI);
    let phi = kn.phi_of_rank(theta, alpha / (2.0 * PI), pol);
    Esito::Diffuso(cinematica(kn.e_kev, theta, phi))
}

/// Small deterministic generator (SplitMix64), so every run with the same seed is identical.
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self { Rng(seed ^ 0x9E37_79B9_7F4A_7C15) }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// uniform in [0, 1)
    pub fn uni(&mut self) -> f64 { (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64) }
    pub fn gauss(&mut self) -> f64 {
        let u1 = 1.0 - self.uni();
        let u2 = self.uni();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

/// Birth of the relative transverse coordinate b (m), ∝ |ψ|² (2-D Gaussian of rms s per axis).
/// `solo_diffusi`: draw it conditioned on the scattering disk (u < p) — the same density,
/// restricted, used to show/collect scattering events without waiting 1/p shots.
pub fn nascita_b(rng: &mut Rng, pk: &Pacchetti, p: f64, solo_diffusi: bool) -> [f64; 2] {
    let s = pk.s_rel();
    let u = if solo_diffusi { p * rng.uni() } else { rng.uni() };
    // invert u = 1 − exp(−b²/2s²) without losing digits when u is tiny
    let b = s * (-2.0 * (-u).ln_1p()).sqrt();
    let a = 2.0 * PI * rng.uni();
    [b * a.cos(), b * a.sin()]
}

/// Electron packet spreading time τ = 2 mₑ s²/ħ (s).
pub fn tau_spread(s_e_m: f64) -> f64 { 2.0 * ME_KG * s_e_m * s_e_m / HBAR }

/// Compton 1923 (Phys. Rev. 21, 483): Mo Kα scattered by graphite at 90°.
/// Values as quoted from the paper (λ primary 0.708 Å, modified line 0.730 Å, i.e. shift 0.022 Å,
/// against the predicted 0.0242 Å) — see docs/COMPTON.md for the caveat on the source.
pub const COMPTON_1923: [(f64, f64, f64); 1] = [(90.0, 70.8, 2.2)];

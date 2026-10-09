//! THE RULE OF THE ELECTRON BALLS (diffusion Monte Carlo written as a rule per ball).
//!
//! A ball is ONE point in the configuration space of all the electrons of the system
//! (1 electron: a point in 3-D; H2 or He: a pair of points, 6 numbers; …).
//! Each tick of duration τ (atomic units: ħ = m_e = e = 4πε0 = 1) every ball, on its own:
//!   1. JITTER — a random Gaussian step of variance 2Dτ per coordinate, D = ħ/2m_e (= ½ a.u.);
//!      with a guide (optional, see below) it also SLIDES by τ·∇ln ψ_g towards where the guide is larger;
//!   2. LIFE — it duplicates or disappears with probability set by its local energy:
//!      weight w = exp(−τ (E_loc − E_ref)), copies = ⌊w + u⌋, u uniform in [0,1).
//!      Without a guide E_loc is just the Coulomb potential energy V at its position
//!      (nuclei and other electrons). E_ref is one global number that keeps the population stable.
//!   3. NODE RULE (excited states / fermions): if it crosses a node surface it dies; if it is near one
//!      it dies with the probability that its random path touched it, exp(−2 d d′ / τ)
//!      (d, d′ = distances from the node before and after the step).
//! After equilibration the cloud of balls IS the ground-state wave function of its region (without a guide
//! the ball density is ψ; with a guide ψ_g it is ψ_g·ψ, i.e. ≈ |ψ|² when the guide is good) and
//! E_ref is the energy. The guide never changes the answer (except through its nodes): it only reduces
//! the noise. Our guides are built from OUR lower level: the 1s orbital fitted from the hydrogen cloud,
//! scaled by Z (exact Coulomb scaling), plus the electron–electron cusp exp(r/2/(1+r)) that the Coulomb
//! law imposes when two electrons touch.

use rayon::prelude::*;

// ───────────────────────── random numbers (xoshiro256++, deterministic per chunk) ─────────────────────────
#[derive(Clone)]
pub struct Rng([u64; 4]);

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut s = [0u64; 4];
        for v in s.iter_mut() {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            *v = x ^ (x >> 31);
        }
        Rng(s)
    }
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.0;
        let r = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        r
    }
    /// Uniform in [0, 1).
    pub fn uni(&mut self) -> f64 { (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64) }
    /// Standard normal (Box–Muller).
    pub fn gauss(&mut self) -> f64 {
        let u1 = 1.0 - self.uni();
        let u2 = self.uni();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

// ───────────────────────── the system the balls live in ─────────────────────────
#[derive(Clone, Copy, Debug)]
pub struct Nucleo {
    pub pos: [f64; 3],
    pub z: f64,
}

/// Optional guide (importance sampling). Built only from our own lower-level orbital.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Guida {
    Nessuna,
    /// every electron in φ(r) = Σ_A exp(−a Z_A |r − R_A|) (our 1s, scaled by Z), × e–e cusp if `jastrow`
    Orbitali { a: f64, jastrow: bool },
    /// one electron, nucleus at the origin: ψ_g = A(x)·Π_k (r − r_k)·exp(−a r), with the angular factor
    /// A = 1, z or x·y taken from the mirror symmetry of the state and r_k its radial nodes; `a` is the
    /// decay measured on the cloud of a first, unguided generation of balls.
    Eccitato { ang: Angolare, nodi: [f64; 2], nn: usize, a: f64 },
    /// two electrons (helium 1s2p ¹P, M = 0): ψ_g = [φ_a(r1)φ_b(r2) + φ_b(r1)φ_a(r2)] × e–e cusp, with our
    /// hydrogen-like orbitals φ_a = e^{−k_a r} (1s) and φ_b = z e^{−k_b r} (2p_z). Spatially symmetric
    /// (singlet), odd under z → −z of both electrons. Its node is the fixed-node rule of this state.
    S1P { ka: f64, kb: f64 },
    /// every electron in our 1s with an effective charge ζ (screening), corrected near the nucleus so that
    /// the electron–nucleus cusp is exact: φ = exp(−ζ r − (Z − ζ) r/(1 + r)), × e–e cusp. One nucleus at the
    /// origin. (ζ from the minimum of our variational energy.)
    Schermato { zeta: f64 },
}

/// Angular factor of a one-electron guide (from mirror symmetry).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Angolare { S, Z, XY }


/// Node surfaces (one-electron excited states; coordinates relative to the origin).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Vincolo {
    /// z > 0 (odd under the z mirror: p_z-like)
    PianoZ,
    /// x > 0
    PianoX,
    /// y > 0
    PianoY,
    /// r < R (inner pocket of a radial node)
    SferaDentro(f64),
    /// r > R (outer pocket of a radial node)
    SferaFuori(f64),
}

impl Vincolo {
    /// Signed distance from the node surface, positive on the allowed side.
    pub fn dist(&self, p: &[f64]) -> f64 {
        match *self {
            Vincolo::PianoZ => p[2],
            Vincolo::PianoX => p[0],
            Vincolo::PianoY => p[1],
            Vincolo::SferaDentro(r) => r - norm(p),
            Vincolo::SferaFuori(r) => norm(p) - r,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Sistema {
    pub nuclei: Vec<Nucleo>,
    /// +1 / −1 per electron
    pub spin: Vec<i8>,
    pub guida: Guida,
    pub vincoli: Vec<Vincolo>,
    /// uniform electric field along z (a.u.): adds F·z per electron, −Z F z per nucleus
    pub campo_z: f64,
}

#[inline]
fn norm(p: &[f64]) -> f64 { (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() }

impl Sistema {
    pub fn ne(&self) -> usize { self.spin.len() }

    /// Coulomb potential energy (hartree) of one configuration x = [x1,y1,z1, x2,…].
    pub fn potenziale(&self, x: &[f64]) -> f64 {
        let ne = self.ne();
        let mut v = 0.0;
        for i in 0..ne {
            let r = &x[3 * i..3 * i + 3];
            for n in &self.nuclei {
                let d = [r[0] - n.pos[0], r[1] - n.pos[1], r[2] - n.pos[2]];
                v -= n.z / norm(&d).max(1e-12);
            }
            for j in 0..i {
                let s = &x[3 * j..3 * j + 3];
                let d = [r[0] - s[0], r[1] - s[1], r[2] - s[2]];
                v += 1.0 / norm(&d).max(1e-12);
            }
            v += self.campo_z * r[2];
        }
        for (a, na) in self.nuclei.iter().enumerate() {
            v -= self.campo_z * na.z * na.pos[2];
            for nb in &self.nuclei[..a] {
                let d = [na.pos[0] - nb.pos[0], na.pos[1] - nb.pos[1], na.pos[2] - nb.pos[2]];
                v += na.z * nb.z / norm(&d);
            }
        }
        v
    }

    /// Local energy E_loc = (Hψ_g)/ψ_g and drift ∇ln ψ_g (written into `drift`). Without guide: V, 0.
    pub fn valuta(&self, x: &[f64], drift: &mut [f64]) -> f64 {
        let v = self.potenziale(x);
        let (a, jas) = match self.guida {
            Guida::Nessuna => {
                drift.iter_mut().for_each(|d| *d = 0.0);
                return v;
            }
            Guida::Orbitali { a, jastrow } => (a, jastrow),
            Guida::S1P { ka, kb } => {
                let lap = self.s1p_parte_orbitale(x, ka, kb, drift);
                return self.con_jastrow(x, drift, lap, v);
            }
            Guida::Schermato { zeta } => {
                let c = self.nuclei[0].z - zeta;
                let mut lap = 0.0;
                for i in 0..self.ne() {
                    let p = &x[3 * i..3 * i + 3];
                    let r = norm(p).max(1e-12);
                    let d1 = -zeta - c / (1.0 + r).powi(2);
                    let d2 = 2.0 * c / (1.0 + r).powi(3);
                    for k in 0..3 { drift[3 * i + k] = d1 * p[k] / r; }
                    lap += d2 + 2.0 * d1 / r;
                }
                return self.con_jastrow(x, drift, lap, v);
            }
            Guida::Eccitato { ang, nodi, nn, a } => {
                let p = &x[0..3];
                let r = norm(p).max(1e-12);
                let mut lap = -2.0 * a / r; // ∇² ln ψ
                let mut g = [-a * p[0] / r, -a * p[1] / r, -a * p[2] / r];
                for k in nodi.iter().take(nn) {
                    let q = r - k;
                    for c in 0..3 { g[c] += p[c] / r / q; }
                    lap += -1.0 / (q * q) + 2.0 / (r * q);
                }
                match ang {
                    Angolare::S => {}
                    Angolare::Z => { g[2] += 1.0 / p[2]; lap -= 1.0 / (p[2] * p[2]); }
                    Angolare::XY => {
                        g[0] += 1.0 / p[0];
                        g[1] += 1.0 / p[1];
                        lap -= 1.0 / (p[0] * p[0]) + 1.0 / (p[1] * p[1]);
                    }
                }
                drift[0..3].copy_from_slice(&g);
                return -0.5 * (lap + g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) + v;
            }
        };
        let ne = self.ne();
        let mut lap_sum = 0.0; // Σ_i ∇_i² ln ψ
        for i in 0..ne {
            let r = &x[3 * i..3 * i + 3];
            // orbital: φ = Σ_A exp(−k_A r_A), shifted for stability
            let mut kr = [0.0f64; 8];
            let mut dv = [[0.0f64; 3]; 8];
            let mut rr = [0.0f64; 8];
            let nn = self.nuclei.len().min(8);
            let mut m = f64::INFINITY;
            for (k, n) in self.nuclei.iter().take(nn).enumerate() {
                let d = [r[0] - n.pos[0], r[1] - n.pos[1], r[2] - n.pos[2]];
                let dd = norm(&d).max(1e-12);
                rr[k] = dd;
                dv[k] = [d[0] / dd, d[1] / dd, d[2] / dd];
                kr[k] = a * n.z;
                m = m.min(kr[k] * dd);
            }
            let (mut phi, mut g, mut lap) = (0.0, [0.0f64; 3], 0.0);
            for k in 0..nn {
                let e = (-(kr[k] * rr[k] - m)).exp();
                phi += e;
                for c in 0..3 { g[c] -= kr[k] * e * dv[k][c]; }
                lap += (kr[k] * kr[k] - 2.0 * kr[k] / rr[k]) * e;
            }
            let gl = [g[0] / phi, g[1] / phi, g[2] / phi];
            lap_sum += lap / phi - (gl[0] * gl[0] + gl[1] * gl[1] + gl[2] * gl[2]);
            drift[3 * i..3 * i + 3].copy_from_slice(&gl);
        }
        if jas { return self.con_jastrow(x, drift, lap_sum, v); }
        let g2: f64 = drift.iter().map(|d| d * d).sum();
        -0.5 * (lap_sum + g2) + v
    }

    /// Adds the electron–electron cusp factor exp(Σ c r_ij/(1+r_ij)) to drift and Laplacian, returns E_loc.
    fn con_jastrow(&self, x: &[f64], drift: &mut [f64], mut lap_sum: f64, v: f64) -> f64 {
        let ne = self.ne();
        for i in 0..ne {
            for j in 0..i {
                let c = if self.spin[i] == self.spin[j] { 0.25 } else { 0.5 };
                let d = [x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]];
                let r = norm(&d).max(1e-12);
                let up = c / (1.0 + r).powi(2); // u'
                let upp = -2.0 * c / (1.0 + r).powi(3); // u''
                for k in 0..3 {
                    let gk = up * d[k] / r;
                    drift[3 * i + k] += gk;
                    drift[3 * j + k] -= gk;
                }
                lap_sum += 2.0 * (upp + 2.0 * up / r);
            }
        }
        let g2: f64 = drift.iter().map(|d| d * d).sum();
        -0.5 * (lap_sum + g2) + v
    }

    /// The two terms of the 1s2p guide, A = φ_a(r1)φ_b(r2) = z2·e_A and B = φ_b(r1)φ_a(r2) = z1·e_B, with
    /// the exponentials e_A, e_B multiplied by the same factor e^{m} (m = smallest exponent) so that
    /// nothing underflows. Returns (e_A, e_B, m, r1, r2).
    pub fn s1p_termini(x: &[f64], ka: f64, kb: f64) -> (f64, f64, f64, f64, f64) {
        let r1 = norm(&x[0..3]).max(1e-12);
        let r2 = norm(&x[3..6]).max(1e-12);
        let (ea, eb) = (ka * r1 + kb * r2, kb * r1 + ka * r2);
        let m = ea.min(eb);
        ((-(ea - m)).exp(), (-(eb - m)).exp(), m, r1, r2)
    }

    /// Orbital part of the 1s2p guide: writes ∇ ln ψ_orb into `drift`, returns Σ_i ∇_i² ln ψ_orb.
    fn s1p_parte_orbitale(&self, x: &[f64], ka: f64, kb: f64, drift: &mut [f64]) -> f64 {
        let (ex_a, ex_b, _, r1, r2) = Self::s1p_termini(x, ka, kb);
        let (a, b) = (ex_a * x[5], ex_b * x[2]);
        let psi = a + b;
        let u1 = [x[0] / r1, x[1] / r1, x[2] / r1];
        let u2 = [x[3] / r2, x[4] / r2, x[5] / r2];
        let zhat = [0.0, 0.0, 1.0];
        let mut lap = 0.0;
        // electron 1: A has φ_a(r1), B has φ_b(r1)
        let mut g = [0.0; 3];
        for c in 0..3 { g[c] = -ka * u1[c] * a + ex_b * (zhat[c] - kb * x[2] * u1[c]); }
        let l1 = (ka * ka - 2.0 * ka / r1) * a + (kb * kb - 4.0 * kb / r1) * b;
        for c in 0..3 { drift[c] = g[c] / psi; }
        lap += l1 / psi - (drift[0] * drift[0] + drift[1] * drift[1] + drift[2] * drift[2]);
        // electron 2: A has φ_b(r2), B has φ_a(r2)
        for c in 0..3 { g[c] = ex_a * (zhat[c] - kb * x[5] * u2[c]) - ka * u2[c] * b; }
        let l2 = (kb * kb - 4.0 * kb / r2) * a + (ka * ka - 2.0 * ka / r2) * b;
        for c in 0..3 { drift[3 + c] = g[c] / psi; }
        lap += l2 / psi - (drift[3] * drift[3] + drift[4] * drift[4] + drift[5] * drift[5]);
        lap
    }

    /// Sign of the guide (only `Eccitato` guides change sign: across their nodes).
    pub fn segno(&self, x: &[f64]) -> f64 {
        match self.guida {
            Guida::Eccitato { ang, nodi, nn, .. } => {
                let p = &x[0..3];
                let r = norm(p);
                let mut s = match ang { Angolare::S => 1.0, Angolare::Z => p[2], Angolare::XY => p[0] * p[1] };
                for k in nodi.iter().take(nn) { s *= r - k; }
                s.signum()
            }
            Guida::S1P { ka, kb } => {
                let (ea, eb, ..) = Self::s1p_termini(x, ka, kb);
                (ea * x[5] + eb * x[2]).signum()
            }
            _ => 1.0,
        }
    }

    /// ln|ψ_g| (only needed by the accept/reject rule).
    pub fn ln_guida(&self, x: &[f64]) -> f64 {
        match self.guida {
            Guida::Nessuna => 0.0,
            Guida::Eccitato { ang, nodi, nn, a } => {
                let p = &x[0..3];
                let r = norm(p);
                let mut l = -a * r;
                for k in nodi.iter().take(nn) { l += (r - k).abs().max(1e-300).ln(); }
                l + match ang { Angolare::S => 0.0, Angolare::Z => p[2].abs().max(1e-300).ln(),
                                Angolare::XY => (p[0] * p[1]).abs().max(1e-300).ln() }
            }
            Guida::Schermato { zeta } => {
                let c = self.nuclei[0].z - zeta;
                let mut l = 0.0;
                for i in 0..self.ne() {
                    let r = norm(&x[3 * i..3 * i + 3]);
                    l += -zeta * r - c * r / (1.0 + r);
                }
                let r = norm(&[x[0] - x[3], x[1] - x[4], x[2] - x[5]]);
                l + 0.5 * r / (1.0 + r)
            }
            Guida::S1P { ka, kb } => {
                let (ea, eb, m, ..) = Self::s1p_termini(x, ka, kb);
                let r = norm(&[x[0] - x[3], x[1] - x[4], x[2] - x[5]]);
                (ea * x[5] + eb * x[2]).abs().max(1e-300).ln() - m + 0.5 * r / (1.0 + r)
            }
            Guida::Orbitali { a, jastrow } => {
                let ne = self.ne();
                let mut l = 0.0;
                for i in 0..ne {
                    let r = &x[3 * i..3 * i + 3];
                    let mut m = f64::INFINITY;
                    let mut ks = [0.0f64; 8];
                    for (k, n) in self.nuclei.iter().take(8).enumerate() {
                        ks[k] = a * n.z * norm(&[r[0] - n.pos[0], r[1] - n.pos[1], r[2] - n.pos[2]]);
                        m = m.min(ks[k]);
                    }
                    let nn = self.nuclei.len().min(8);
                    l += -m + ks[..nn].iter().map(|k| (-(k - m)).exp()).sum::<f64>().ln();
                }
                if jastrow {
                    for i in 0..ne {
                        for j in 0..i {
                            let c = if self.spin[i] == self.spin[j] { 0.25 } else { 0.5 };
                            let r = norm(&[x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]]);
                            l += c * r / (1.0 + r);
                        }
                    }
                }
                l
            }
        }
    }

    /// Smallest signed distance from the node rules (∞ if none).
    pub fn dentro(&self, x: &[f64]) -> f64 {
        self.vincoli.iter().map(|c| c.dist(&x[0..3])).fold(f64::INFINITY, f64::min)
    }
}

// ───────────────────────── the population of balls ─────────────────────────
pub struct Dmc {
    pub sys: Sistema,
    pub tau: f64,
    /// target population
    pub n0: usize,
    /// ball positions, 3·ne numbers per ball
    pub x: Vec<f64>,
    el: Vec<f64>,
    dr: Vec<f64>,
    pub e_ref: f64,
    e_bar: f64,
    pub passi: u64,
    seed: u64,
    /// last estimates: mixed ⟨E_loc⟩ and growth (from the population change)
    pub e_mix: f64,
    pub e_growth: f64,
    /// Optional per-ball accept/reject of the move (Reynolds et al. 1982): with a guide, the ball keeps
    /// its step with probability min(1, ψ_g²(new)·G(old←new) / ψ_g²(old)·G(new←old)), else it stays.
    /// It removes most of the finite-tick error; the branching then uses τ·(acceptance).
    pub metropolis: bool,
    acc: f64,
}

const CHUNK: usize = 256;

/// Drift of one ball, capped per electron (Umrigar–Nightingale–Runge 1993): v·(√(1+2v²τ)−1)/(v²τ).
fn drift_capped(d: &[f64], tau: f64, out: &mut [f64]) {
    for e in 0..d.len() / 3 {
        let v = &d[3 * e..3 * e + 3];
        let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let f = if v2 * tau > 1e-12 { (-1.0 + (1.0 + 2.0 * v2 * tau).sqrt()) / (v2 * tau) } else { 1.0 };
        for k in 0..3 { out[3 * e + k] = f * v[k]; }
    }
}
/// Relaxation time (a.u.) of the population control.
const T_POP: f64 = 1.0;

impl Dmc {
    pub fn new(sys: Sistema, tau: f64, n0: usize, seed: u64) -> Self {
        let ne = sys.ne();
        let mut rng = Rng::new(seed ^ 0xA5A5_5A5A);
        let mut x = Vec::with_capacity(n0 * 3 * ne);
        let mut k = 0usize;
        while x.len() < n0 * 3 * ne {
            let mut w = vec![0.0; 3 * ne];
            for i in 0..ne {
                let n = sys.nuclei[(k + i) % sys.nuclei.len()];
                let s = 1.0 / n.z.max(0.5);
                for c in 0..3 { w[3 * i + c] = n.pos[c] + s * rng.gauss(); }
            }
            // for node rules: start well inside the allowed pocket (scale by a random radius)
            if !sys.vincoli.is_empty() {
                for c in 0..3 { w[c] *= 1.0 + 6.0 * rng.uni(); }
                if sys.dentro(&w) <= 0.0 { continue; }
            }
            x.extend_from_slice(&w);
            k += 1;
        }
        let mut d = Dmc { sys, tau, n0, x, el: vec![], dr: vec![], e_ref: 0.0, e_bar: 0.0, passi: 0, seed,
                          e_mix: 0.0, e_growth: 0.0, metropolis: false, acc: 1.0 };
        d.ricalcola();
        let n = d.n();
        d.e_ref = d.el.iter().sum::<f64>() / n as f64;
        d.e_bar = d.e_ref;
        d
    }

    pub fn n(&self) -> usize { self.x.len() / (3 * self.sys.ne()) }

    /// Recompute the cached local energies and drifts (after the nuclei or the field changed).
    pub fn ricalcola(&mut self) {
        let m = 3 * self.sys.ne();
        let n = self.n();
        self.el = vec![0.0; n];
        self.dr = vec![0.0; n * m];
        let sys = &self.sys;
        self.el.par_iter_mut().zip(self.dr.par_chunks_mut(m)).zip(self.x.par_chunks(m))
            .for_each(|((e, d), x)| *e = sys.valuta(x, d));
        // balls that are now outside the allowed region die
        if !self.sys.vincoli.is_empty() {
            let keep: Vec<usize> = (0..n).filter(|&i| self.sys.dentro(&self.x[i * m..(i + 1) * m]) > 0.0).collect();
            if keep.len() < n {
                self.x = keep.iter().flat_map(|&i| self.x[i * m..(i + 1) * m].to_vec()).collect();
                self.dr = keep.iter().flat_map(|&i| self.dr[i * m..(i + 1) * m].to_vec()).collect();
                self.el = keep.iter().map(|&i| self.el[i]).collect();
            }
        }
    }

    /// One tick for every ball. Returns (mixed estimate, growth estimate) of the energy.
    pub fn passo(&mut self) -> (f64, f64) {
        let m = 3 * self.sys.ne();
        let n = self.n();
        let (tau, e_ref, sys) = (self.tau, self.e_ref, &self.sys);
        let sq = tau.sqrt();
        let clip = 10.0 / sq;
        let step = self.passi;
        let seed = self.seed;
        let guided = sys.guida != Guida::Nessuna;
        let has_nodes = !sys.vincoli.is_empty();
        let metro = self.metropolis && guided;
        let nodal_guide = matches!(sys.guida, Guida::Eccitato { .. } | Guida::S1P { .. });
        let tau_eff = if metro { tau * self.acc } else { tau };
        struct Out { x: Vec<f64>, dr: Vec<f64>, el: Vec<f64>, sw: f64, swe: f64, acc: f64, tried: f64 }
        let outs: Vec<Out> = (0..n.div_ceil(CHUNK)).into_par_iter().map(|c| {
            let mut rng = Rng::new(seed.wrapping_mul(0x2545_F491_4F6C_DD1D) ^ step.wrapping_mul(0x9E37_79B9) ^ ((c as u64) << 40));
            let mut o = Out { x: vec![], dr: vec![], el: vec![], sw: 0.0, swe: 0.0, acc: 0.0, tried: 0.0 };
            let mut y = vec![0.0; m];
            let mut dy = vec![0.0; m];
            let mut vx = vec![0.0; m];
            let mut vy = vec![0.0; m];
            for i in c * CHUNK..((c + 1) * CHUNK).min(n) {
                let xi = &self.x[i * m..(i + 1) * m];
                let di = &self.dr[i * m..(i + 1) * m];
                // 1. jitter (+ slide along the guide)
                if guided { drift_capped(di, tau, &mut vx); } else { vx.iter_mut().for_each(|v| *v = 0.0); }
                for j in 0..m { y[j] = xi[j] + tau * vx[j] + sq * rng.gauss(); }
                // 3. node rule
                if has_nodes {
                    if sys.dentro(&y) <= 0.0 { continue; }
                    let mut surv = 1.0;
                    for cst in &sys.vincoli {
                        let (a, b) = (cst.dist(&xi[0..3]), cst.dist(&y[0..3]));
                        surv *= 1.0 - (-2.0 * a * b / tau).exp();
                    }
                    if rng.uni() > surv { continue; }
                }
                // fixed-node rule for guides with nodes: a step across a node of the guide is not allowed
                let crosses = nodal_guide && sys.segno(&y) != sys.segno(xi);
                if crosses && !metro { continue; }
                let mut e1 = sys.valuta(&y, &mut dy);
                let e0 = self.el[i];
                if metro {
                    // optional accept/reject of the move, from the guide only
                    drift_capped(&dy, tau, &mut vy);
                    let (mut gf, mut gb) = (0.0, 0.0);
                    for j in 0..m {
                        gf += (y[j] - xi[j] - tau * vx[j]).powi(2);
                        gb += (xi[j] - y[j] - tau * vy[j]).powi(2);
                    }
                    let lr = 2.0 * (sys.ln_guida(&y) - sys.ln_guida(xi)) + (gf - gb) / (2.0 * tau);
                    let p = if crosses { 0.0 } else { lr.min(0.0).exp() };
                    o.tried += 1.0;
                    o.acc += p;
                    if rng.uni() >= p {
                        y.copy_from_slice(xi);
                        dy.copy_from_slice(di);
                        e1 = e0;
                    }
                }
                // 2. life: duplicate or disappear
                let em = (0.5 * (e0 + e1) - e_ref).clamp(-clip, clip);
                let w = (-tau_eff * em).exp();
                o.sw += w;
                o.swe += w * e1;
                let copies = ((w + rng.uni()) as usize).min(8);
                for _ in 0..copies {
                    o.x.extend_from_slice(&y);
                    o.dr.extend_from_slice(&dy);
                    o.el.push(e1);
                }
            }
            o
        }).collect();
        let (acc, tried): (f64, f64) = outs.iter().fold((0.0, 0.0), |s, o| (s.0 + o.acc, s.1 + o.tried));
        if tried > 0.0 { self.acc = acc / tried; }
        let (mut sw, mut swe) = (0.0, 0.0);
        let mut x = Vec::with_capacity(self.x.len() + m * 64);
        let mut dr = Vec::with_capacity(self.dr.len() + m * 64);
        let mut el = Vec::with_capacity(n + 64);
        for o in outs {
            sw += o.sw;
            swe += o.swe;
            x.extend(o.x);
            dr.extend(o.dr);
            el.extend(o.el);
        }
        self.x = x;
        self.dr = dr;
        self.el = el;
        self.passi += 1;
        let e_mix = if sw > 0.0 { swe / sw } else { self.e_ref };
        let e_growth = self.e_ref - (sw.max(1e-300) / n as f64).ln() / tau_eff;
        self.e_mix = e_mix;
        self.e_growth = e_growth;
        // population control: E_ref follows the energy and pulls the population back to n0
        let est = if guided { e_mix } else { e_growth };
        let k = (tau_eff / T_POP).min(1.0);
        self.e_bar += k * (est - self.e_bar);
        let nn = self.n().max(1) as f64;
        self.e_ref = self.e_bar - (nn / self.n0 as f64).ln() / T_POP;
        if self.n() == 0 {
            panic!("all balls died (tau too large for the node rule?)");
        }
        (e_mix, e_growth)
    }

    /// Positions of electron `e` of every ball (for drawing and histograms).
    pub fn elettrone(&self, e: usize) -> impl Iterator<Item = [f64; 3]> + '_ {
        let m = 3 * self.sys.ne();
        self.x.chunks(m).map(move |w| [w[3 * e], w[3 * e + 1], w[3 * e + 2]])
    }
}

// ───────────────────────── statistics ─────────────────────────
/// Mean and statistical error by blocking: the largest error among 32, 16 and 8 blocks of consecutive
/// ticks (conservative: the slow clouds of excited states stay correlated for tens of a.u. of time).
pub fn media_errore(v: &[f64]) -> (f64, f64) {
    let n = v.len();
    let m = v.iter().sum::<f64>() / n.max(1) as f64;
    if n < 64 { return (m, f64::NAN); }
    let mut err: f64 = 0.0;
    for nb in [32usize, 16, 8] {
        let bl = n / nb;
        let means: Vec<f64> = (0..nb).map(|b| v[b * bl..(b + 1) * bl].iter().sum::<f64>() / bl as f64).collect();
        let mm = means.iter().sum::<f64>() / nb as f64;
        let var = means.iter().map(|x| (x - mm) * (x - mm)).sum::<f64>() / (nb - 1) as f64;
        err = err.max((var / nb as f64).sqrt());
    }
    (m, err)
}

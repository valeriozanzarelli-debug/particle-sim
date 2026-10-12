//! THE CARBON ATOM FROM THE ELECTRON-BALL RULE (fixed-node diffusion Monte Carlo written per ball).
//!
//! A ball is ONE point in the 18-dimensional space of the 6 electrons (15 for the ion C⁺). Each tick τ
//! (atomic units: ħ = mₑ = e = 4πε₀ = 1) every ball, on its own, does exactly what the hydrogen, H₂ and
//! helium balls do in src/molecole/dmc.rs:
//!   1. JITTER (variance τ per coordinate, D = ħ/2mₑ) + SLIDE τ·∇ln|ψ_g| along the guide;
//!   2. ACCEPT/REJECT the move with the guide only (Reynolds et al. 1982), and NEVER across a node of
//!      the guide (fixed-node rule: this is how fermions with 3 or more electrons are handled);
//!   3. LIFE: copies = ⌊exp(−τ(E_loc − E_ref)) + u⌋.
//!
//! THE GUIDE is built from our own hydrogen-like orbitals (nothing tabulated):
//!   1s = exp(−ζ₁r − (Z−ζ₁) r/(1+r))      (the Coulomb cusp ∂lnφ/∂r = −Z on the nucleus is exact)
//!   2s = (1 − (Z−ζ_s) r) e^{−ζ_s r}       (hydrogen-like, cusp exact)
//!   2p = (x, y, z) e^{−ζ_p r}             (hydrogen-like)
//!   3s = (1 − (Z−ζ₃) r + d₃ r²) e^{−ζ₃ r}  (hydrogen-like, cusp exact; only for the excited state)
//! antisymmetrised with one Slater determinant per spin, times the electron–electron cusp factor
//! exp(Σ a r/(1+b r)) (a = ½ antiparallel, ¼ parallel: both are the Coulomb law at contact).
//! The exponents ζ are NOT taken from tables: they are the minimum of our own variational energy
//! (correlated-sampling VMC, `ottimizza`).
//!
//! Spin. The ³P ground state of carbon (1s²2s²2p²) is written as its M_S = 1 component: 4 up
//! (1s 2s 2p_x 2p_y) + 2 down (1s 2s). The Hamiltonian has no spin, so every M_S component has the same
//! energy, and M_S = 1 is a SINGLE determinant pair (M_S = 0 would need two).

use rayon::prelude::*;

pub use rng::Rng;
mod rng {
    /// xoshiro256++ (same generator as src/molecole/dmc.rs).
    #[derive(Clone, Debug)]
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
        pub fn uni(&mut self) -> f64 { (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64) }
        pub fn gauss(&mut self) -> f64 {
            let u1 = 1.0 - self.uni();
            let u2 = self.uni();
            (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
        }
    }
}

/// Orbital labels (index into `Guida::orb`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orb { S1, S2, Px, Py, Pz, S3 }

/// The three states the viewer can show.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stato {
    /// C 1s²2s²2p² ³P (ground state)
    Fondamentale,
    /// C 1s²2s²2p3s ³P° (the first excited state reached by an allowed photon: 2p → 3s)
    Eccitato,
    /// C⁺ 1s²2s²2p ²P° (the ion, after photoionization)
    Ione,
}

impl Stato {
    pub fn up(self) -> &'static [Orb] {
        match self {
            Stato::Fondamentale => &[Orb::S1, Orb::S2, Orb::Px, Orb::Py],
            Stato::Eccitato => &[Orb::S1, Orb::S2, Orb::Px, Orb::S3],
            Stato::Ione => &[Orb::S1, Orb::S2, Orb::Px],
        }
    }
    pub fn dn(self) -> &'static [Orb] { &[Orb::S1, Orb::S2] }
    pub fn ne(self) -> usize { self.up().len() + self.dn().len() }
    pub fn nome(self) -> &'static str {
        match self {
            Stato::Fondamentale => "C  1s²2s²2p²  ³P",
            Stato::Eccitato => "C*  1s²2s²2p3s  ³P°",
            Stato::Ione => "C⁺  1s²2s²2p  ²P°",
        }
    }
}

/// Variational parameters of the guide (exponents in 1/bohr).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Param {
    pub z1: f64,
    pub zs: f64,
    pub zp: f64,
    pub z3: f64,
    pub d3: f64,
    /// Jastrow range parameter b (fixed, not optimised)
    pub b: f64,
}

impl Param {
    /// Starting point: the bare hydrogen-like exponents Z/n (no screening at all).
    pub fn idrogenoidi(z: f64) -> Self {
        Param { z1: z, zs: z / 2.0, zp: z / 2.0, z3: z / 3.0, d3: 2.0 * z * z / 27.0, b: 1.0 }
    }
    pub fn get(&self, k: usize) -> f64 { [self.z1, self.zs, self.zp, self.z3, self.d3][k] }
    pub fn set(&mut self, k: usize, v: f64) {
        match k { 0 => self.z1 = v, 1 => self.zs = v, 2 => self.zp = v, 3 => self.z3 = v, _ => self.d3 = v }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Guida {
    pub z: f64,
    pub stato: Stato,
    pub p: Param,
}

#[inline]
fn norm3(p: &[f64]) -> f64 { (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() }

/// Radial s function P(r) e^{−ζr}, P = 1 − c r + d r²: (value, d/dr, d²/dr²).
#[inline]
fn s_poly(r: f64, z: f64, c: f64, d: f64) -> (f64, f64, f64) {
    let e = (-z * r).exp();
    let p = 1.0 - c * r + d * r * r;
    let p1 = -c + 2.0 * d * r;
    let p2 = 2.0 * d;
    (p * e, (p1 - z * p) * e, (p2 - 2.0 * z * p1 + z * z * p) * e)
}

impl Guida {
    /// Orbital k at point q: (value, gradient, laplacian).
    #[inline]
    pub fn orb(&self, k: Orb, q: &[f64]) -> (f64, [f64; 3], f64) {
        let r = norm3(q).max(1e-12);
        let u = [q[0] / r, q[1] / r, q[2] / r];
        let zz = self.z;
        let p = &self.p;
        let radial = |f: f64, f1: f64, f2: f64| (f, [f1 * u[0], f1 * u[1], f1 * u[2]], f2 + 2.0 * f1 / r);
        match k {
            Orb::S1 => {
                let g1 = p.z1 + (zz - p.z1) / ((1.0 + r) * (1.0 + r));
                let g2 = -2.0 * (zz - p.z1) / (1.0 + r).powi(3);
                let f = (-(p.z1 * r + (zz - p.z1) * r / (1.0 + r))).exp();
                radial(f, -g1 * f, (g1 * g1 - g2) * f)
            }
            Orb::S2 => {
                let (f, f1, f2) = s_poly(r, p.zs, zz - p.zs, 0.0);
                radial(f, f1, f2)
            }
            Orb::S3 => {
                let (f, f1, f2) = s_poly(r, p.z3, zz - p.z3, p.d3);
                radial(f, f1, f2)
            }
            Orb::Px | Orb::Py | Orb::Pz => {
                let a = match k { Orb::Px => 0, Orb::Py => 1, _ => 2 };
                let e = (-p.zp * r).exp();
                let x = q[a];
                let mut g = [-p.zp * x * u[0] * e, -p.zp * x * u[1] * e, -p.zp * x * u[2] * e];
                g[a] += e;
                (x * e, g, x * e * (p.zp * p.zp - 4.0 * p.zp / r))
            }
        }
    }

    pub fn ne(&self) -> usize { self.stato.ne() }

    /// Local energy E_loc = (Hψ_g)/ψ_g, drift ∇ln|ψ_g| (into `drift`), sign of ψ_g and ln|ψ_g|.
    /// Configuration x = [up electrons…, down electrons…], 3 numbers each.
    pub fn valuta(&self, x: &[f64], drift: &mut [f64]) -> (f64, f64, f64) {
        let ne = self.ne();
        let nu = self.stato.up().len();
        let mut lap = [0.0f64; 6];
        let mut gd = [[0.0f64; 3]; 6];
        let mut sign = 1.0;
        let mut lnpsi = 0.0;
        for (start, orbs) in [(0usize, self.stato.up()), (nu, self.stato.dn())] {
            let m = orbs.len();
            let mut a = [[0.0f64; 5]; 5];
            let mut ga = [[[0.0f64; 3]; 5]; 5];
            let mut la = [[0.0f64; 5]; 5];
            for i in 0..m {
                let q = &x[3 * (start + i)..3 * (start + i) + 3];
                for (k, &o) in orbs.iter().enumerate() {
                    let (v, g, l) = self.orb(o, q);
                    a[i][k] = v;
                    ga[i][k] = g;
                    la[i][k] = l;
                }
            }
            let Some((inv, det)) = inverti(&a, m) else {
                drift.iter_mut().for_each(|d| *d = 0.0);
                return (f64::NAN, 0.0, f64::NEG_INFINITY);
            };
            sign *= det.signum();
            lnpsi += det.abs().ln();
            for i in 0..m {
                let e = start + i;
                for k in 0..m {
                    let w = inv[k][i];
                    for c in 0..3 { gd[e][c] += ga[i][k][c] * w; }
                    lap[e] += la[i][k] * w;
                }
            }
        }
        // Jastrow (electron–electron cusp) and potential energy
        let mut gj = [[0.0f64; 3]; 6];
        let mut lj = [0.0f64; 6];
        let mut v = 0.0;
        let b = self.p.b;
        for i in 0..ne {
            let qi = &x[3 * i..3 * i + 3];
            v -= self.z / norm3(qi).max(1e-12);
            for j in i + 1..ne {
                let qj = &x[3 * j..3 * j + 3];
                let d = [qi[0] - qj[0], qi[1] - qj[1], qi[2] - qj[2]];
                let r = norm3(&d).max(1e-12);
                v += 1.0 / r;
                let par = (i < nu) == (j < nu);
                let aa = if par { 0.25 } else { 0.5 };
                let den = 1.0 + b * r;
                lnpsi += aa * r / den;
                let u1 = aa / (den * den);
                let u2 = -2.0 * aa * b / (den * den * den);
                for c in 0..3 {
                    gj[i][c] += u1 * d[c] / r;
                    gj[j][c] -= u1 * d[c] / r;
                }
                let l = u2 + 2.0 * u1 / r;
                lj[i] += l;
                lj[j] += l;
            }
        }
        let mut kin = 0.0;
        for e in 0..ne {
            let g = [gd[e][0] + gj[e][0], gd[e][1] + gj[e][1], gd[e][2] + gj[e][2]];
            let cross = gd[e][0] * gj[e][0] + gd[e][1] * gj[e][1] + gd[e][2] * gj[e][2];
            let gj2 = gj[e][0] * gj[e][0] + gj[e][1] * gj[e][1] + gj[e][2] * gj[e][2];
            kin += -0.5 * (lap[e] + 2.0 * cross + lj[e] + gj2);
            drift[3 * e..3 * e + 3].copy_from_slice(&g);
        }
        (kin + v, sign, lnpsi)
    }

    /// ln|ψ_g| and sign only (for reweighting).
    pub fn ln_psi(&self, x: &[f64]) -> (f64, f64) {
        let mut d = [0.0; 18];
        let (_, s, l) = self.valuta(x, &mut d[..3 * self.ne()]);
        (l, s)
    }
}

/// Inverse and determinant of the top-left m×m block (Gauss–Jordan with partial pivoting).
fn inverti(a: &[[f64; 5]; 5], m: usize) -> Option<([[f64; 5]; 5], f64)> {
    let mut w = *a;
    let mut inv = [[0.0f64; 5]; 5];
    for i in 0..m { inv[i][i] = 1.0; }
    let mut det = 1.0;
    for c in 0..m {
        let mut piv = c;
        for r in c + 1..m { if w[r][c].abs() > w[piv][c].abs() { piv = r; } }
        if w[piv][c].abs() < 1e-300 { return None; }
        if piv != c {
            w.swap(piv, c);
            inv.swap(piv, c);
            det = -det;
        }
        let p = w[c][c];
        det *= p;
        for k in 0..m { w[c][k] /= p; inv[c][k] /= p; }
        for r in 0..m {
            if r != c {
                let f = w[r][c];
                if f != 0.0 {
                    for k in 0..m { w[r][k] -= f * w[c][k]; inv[r][k] -= f * inv[c][k]; }
                }
            }
        }
    }
    Some((inv, det))
}

// ───────────────────────── the population of balls ─────────────────────────
#[derive(Clone)]
pub struct Dmc {
    pub g: Guida,
    pub tau: f64,
    pub n0: usize,
    pub x: Vec<f64>,
    el: Vec<f64>,
    dr: Vec<f64>,
    sg: Vec<f64>,
    lp: Vec<f64>,
    pub e_ref: f64,
    e_bar: f64,
    pub passi: u64,
    pub seed: u64,
    pub e_mix: f64,
    pub acc: f64,
}

const CHUNK: usize = 128;
/// Relaxation time (a.u.) of the population control.
const T_POP: f64 = 1.0;

fn drift_capped(d: &[f64], tau: f64, out: &mut [f64]) {
    for e in 0..d.len() / 3 {
        let v = &d[3 * e..3 * e + 3];
        let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let f = if v2 * tau > 1e-12 { (-1.0 + (1.0 + 2.0 * v2 * tau).sqrt()) / (v2 * tau) } else { 1.0 };
        for k in 0..3 { out[3 * e + k] = f * v[k]; }
    }
}

/// Starting configurations: electrons near the nucleus at the radius of their hydrogen-like shell.
fn partenza(g: &Guida, rng: &mut Rng) -> Vec<f64> {
    let mut w = Vec::with_capacity(3 * g.ne());
    for o in g.stato.up().iter().chain(g.stato.dn()) {
        let s = match o { Orb::S1 => 1.0 / g.z, Orb::S3 => 3.0, _ => 1.0 };
        for _ in 0..3 { w.push(s * rng.gauss()); }
    }
    w
}

impl Dmc {
    pub fn new(g: Guida, tau: f64, n0: usize, seed: u64) -> Self {
        let ne = g.ne();
        let mut rng = Rng::new(seed ^ 0xC0FF_EE00);
        let mut x = Vec::with_capacity(n0 * 3 * ne);
        while x.len() < n0 * 3 * ne {
            let w = partenza(&g, &mut rng);
            let (l, _) = g.ln_psi(&w);
            if l.is_finite() { x.extend_from_slice(&w); }
        }
        let mut d = Dmc { g, tau, n0, x, el: vec![], dr: vec![], sg: vec![], lp: vec![], e_ref: 0.0, e_bar: 0.0,
                          passi: 0, seed, e_mix: 0.0, acc: 1.0 };
        d.ricalcola();
        // a short VMC (no branching) puts the balls on |ψ_g|² before the rule starts
        d.vmc_passi(200);
        let n = d.n();
        d.e_ref = d.el.iter().sum::<f64>() / n as f64;
        d.e_bar = d.e_ref;
        d.e_mix = d.e_ref;
        d
    }

    pub fn n(&self) -> usize { self.x.len() / (3 * self.g.ne()) }

    pub fn ricalcola(&mut self) {
        let m = 3 * self.g.ne();
        let n = self.n();
        self.el = vec![0.0; n];
        self.sg = vec![0.0; n];
        self.lp = vec![0.0; n];
        self.dr = vec![0.0; n * m];
        let g = &self.g;
        self.el.par_iter_mut().zip(self.sg.par_iter_mut()).zip(self.lp.par_iter_mut())
            .zip(self.dr.par_chunks_mut(m)).zip(self.x.par_chunks(m))
            .for_each(|((((e, s), l), d), x)| { let r = g.valuta(x, d); *e = r.0; *s = r.1; *l = r.2; });
    }

    /// Change the state (ground / excited / ion): the electron count may change, so the balls restart
    /// from the new guide (a fresh population with the same seed stream).
    pub fn cambia(&mut self, g: Guida) {
        let seed = self.seed.wrapping_add(self.passi);
        *self = Dmc::new(g, self.tau, self.n0, seed);
    }

    /// Switch the guide keeping the balls where they are (`x` already has the right number of electrons):
    /// the cloud then RELAXES into the new state under the same rule (no jump, no restart).
    pub fn da_configurazioni(&mut self, g: Guida, x: Vec<f64>) {
        self.g = g;
        self.x = x;
        self.ricalcola();
        let m = 3 * g.ne();
        let keep: Vec<usize> = (0..self.n()).filter(|&i| self.el[i].is_finite() && self.lp[i].is_finite()).collect();
        if keep.len() < self.n() {
            self.x = keep.iter().flat_map(|&i| self.x[i * m..(i + 1) * m].to_vec()).collect();
            self.ricalcola();
        }
        let n = self.n().max(1);
        // robust start of E_ref: the median local energy (a few balls sit near the new nodes)
        let mut e: Vec<f64> = self.el.clone();
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.e_ref = e[n / 2];
        self.e_bar = self.e_ref;
        self.e_mix = self.e_ref;
    }

    /// Plain Metropolis-Langevin steps on |ψ_g|² (no branching): used to equilibrate and for VMC.
    pub fn vmc_passi(&mut self, k: usize) {
        for _ in 0..k { self.muovi(false); }
    }

    /// One tick for every ball. Returns the mixed energy estimate of this tick.
    pub fn passo(&mut self) -> f64 { self.muovi(true) }

    fn muovi(&mut self, branch: bool) -> f64 {
        let m = 3 * self.g.ne();
        let n = self.n();
        let (tau, e_ref, g) = (self.tau, self.e_ref, &self.g);
        let sq = tau.sqrt();
        let clip = 0.2 / tau;
        let step = self.passi;
        let seed = self.seed;
        let tau_eff = tau * self.acc;
        struct Out { x: Vec<f64>, dr: Vec<f64>, el: Vec<f64>, sg: Vec<f64>, lp: Vec<f64>, sw: f64, swe: f64, acc: f64 }
        let (x0, dr0, el0, sg0, lp0) = (&self.x, &self.dr, &self.el, &self.sg, &self.lp);
        let outs: Vec<Out> = (0..n.div_ceil(CHUNK)).into_par_iter().map(|c| {
            let mut rng = Rng::new(seed.wrapping_mul(0x2545_F491_4F6C_DD1D) ^ step.wrapping_mul(0x9E37_79B9) ^ ((c as u64) << 40)
                                   ^ if branch { 0 } else { 0x5555_0000_0000_0000 });
            let cap = (CHUNK * 3 / 2).min(n);
            let mut o = Out { x: Vec::with_capacity(cap * m), dr: Vec::with_capacity(cap * m), el: Vec::with_capacity(cap),
                              sg: Vec::with_capacity(cap), lp: Vec::with_capacity(cap), sw: 0.0, swe: 0.0, acc: 0.0 };
            let mut y = [0.0f64; 18];
            let mut dy = [0.0f64; 18];
            let mut vx = [0.0f64; 18];
            let mut vy = [0.0f64; 18];
            for i in c * CHUNK..((c + 1) * CHUNK).min(n) {
                let xi = &x0[i * m..(i + 1) * m];
                let di = &dr0[i * m..(i + 1) * m];
                drift_capped(di, tau, &mut vx[..m]);
                for j in 0..m { y[j] = xi[j] + tau * vx[j] + sq * rng.gauss(); }
                let (mut e1, s1, l1) = g.valuta(&y[..m], &mut dy[..m]);
                let (e0, s0, l0) = (el0[i], sg0[i], lp0[i]);
                let mut s_new = s1;
                let mut l_new = l1;
                // accept / reject from the guide only; never across a node (fixed-node rule)
                let p = if !e1.is_finite() || s1 != s0 { 0.0 } else {
                    drift_capped(&dy[..m], tau, &mut vy[..m]);
                    let (mut gf, mut gb) = (0.0, 0.0);
                    for j in 0..m {
                        gf += (y[j] - xi[j] - tau * vx[j]).powi(2);
                        gb += (xi[j] - y[j] - tau * vy[j]).powi(2);
                    }
                    (2.0 * (l1 - l0) + (gf - gb) / (2.0 * tau)).min(0.0).exp()
                };
                o.acc += p;
                if rng.uni() >= p {
                    y[..m].copy_from_slice(xi);
                    dy[..m].copy_from_slice(di);
                    e1 = e0;
                    s_new = s0;
                    l_new = l0;
                }
                let copies = if branch {
                    let em = (0.5 * (e0 + e1) - e_ref).clamp(-clip, clip);
                    let w = (-tau_eff * em).exp();
                    o.sw += w;
                    o.swe += w * e1;
                    ((w + rng.uni()) as usize).min(4)
                } else {
                    o.sw += 1.0;
                    o.swe += e1;
                    1
                };
                for _ in 0..copies {
                    o.x.extend_from_slice(&y[..m]);
                    o.dr.extend_from_slice(&dy[..m]);
                    o.el.push(e1);
                    o.sg.push(s_new);
                    o.lp.push(l_new);
                }
            }
            o
        }).collect();
        let acc: f64 = outs.iter().map(|o| o.acc).sum::<f64>() / n.max(1) as f64;
        self.acc = 0.9 * self.acc + 0.1 * acc;
        let (mut sw, mut swe) = (0.0, 0.0);
        let mut x = Vec::with_capacity(self.x.len() + m * 64);
        let mut dr = Vec::with_capacity(self.dr.len() + m * 64);
        let mut el = Vec::with_capacity(n + 64);
        let mut sg = Vec::with_capacity(n + 64);
        let mut lp = Vec::with_capacity(n + 64);
        for o in outs {
            sw += o.sw;
            swe += o.swe;
            x.extend(o.x);
            dr.extend(o.dr);
            el.extend(o.el);
            sg.extend(o.sg);
            lp.extend(o.lp);
        }
        self.x = x;
        self.dr = dr;
        self.el = el;
        self.sg = sg;
        self.lp = lp;
        self.passi += 1;
        let e_mix = if sw > 0.0 { swe / sw } else { self.e_ref };
        self.e_mix = e_mix;
        if branch {
            let k = (tau_eff / T_POP).min(1.0);
            self.e_bar += k * (e_mix - self.e_bar);
            let nn = self.n().max(1) as f64;
            self.e_ref = self.e_bar - (nn / self.n0 as f64).ln() / T_POP;
        }
        assert!(self.n() > 0, "all balls died");
        e_mix
    }

    /// Positions of electron `e` of every ball.
    pub fn elettrone(&self, e: usize) -> impl Iterator<Item = [f64; 3]> + '_ {
        let m = 3 * self.g.ne();
        self.x.chunks(m).map(move |w| [w[3 * e], w[3 * e + 1], w[3 * e + 2]])
    }
}

/// Mean and error by blocking (largest error among 32, 16, 8 blocks), as in src/molecole/dmc.rs.
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

// ───────────────────────── variational exponents (our own VMC) ─────────────────────────
/// Which parameters each state optimises.
fn liberi(s: Stato) -> &'static [usize] {
    match s {
        Stato::Fondamentale | Stato::Ione => &[0, 1, 2],
        Stato::Eccitato => &[3, 4],
    }
}

/// Correlated-sampling VMC energy of parameters `p` on configurations drawn from `p0`.
fn energia_pesata(z: f64, s: Stato, p: &Param, conf: &[f64], l0: &[f64]) -> (f64, f64) {
    let g = Guida { z, stato: s, p: *p };
    let m = 3 * s.ne();
    let r: Vec<(f64, f64)> = conf.par_chunks(m).zip(l0.par_iter()).map(|(x, &l)| {
        let mut d = [0.0; 18];
        let (e, _, lp) = g.valuta(x, &mut d[..m]);
        if !e.is_finite() { return (0.0, 0.0); }
        let w = (2.0 * (lp - l)).exp();
        (w, w * e)
    }).collect();
    let sw: f64 = r.iter().map(|v| v.0).sum();
    let sw2: f64 = r.iter().map(|v| v.0 * v.0).sum();
    let swe: f64 = r.iter().map(|v| v.1).sum();
    (swe / sw, sw * sw / sw2.max(1e-300))
}

/// Minimise our variational energy over the free exponents of state `s`, starting from `p`.
/// Returns the optimised parameters and the final VMC energy (± error).
pub fn ottimizza(z: f64, s: Stato, mut p: Param, n_conf: usize, giri: usize, seed: u64, log: bool) -> (Param, f64, f64) {
    for giro in 0..giri {
        // sample |ψ_p|²
        let mut d = Dmc::new(Guida { z, stato: s, p }, 0.01, 2000, seed + giro as u64);
        d.vmc_passi(300);
        let mut conf = Vec::new();
        let mut l0 = Vec::new();
        while conf.len() < n_conf * 3 * s.ne() {
            d.vmc_passi(10);
            conf.extend_from_slice(&d.x);
            l0.extend_from_slice(&d.lp);
        }
        let n = l0.len();
        let (e_start, _) = energia_pesata(z, s, &p, &conf, &l0);
        // coordinate descent with shrinking steps, staying where the reweighting is trustworthy
        let mut best = e_start;
        for &step in &[0.08, 0.04, 0.02, 0.01] {
            for _sweep in 0..3 {
                for &k in liberi(s) {
                    for dir in [1.0, -1.0] {
                        loop {
                            let mut q = p;
                            let v = p.get(k);
                            q.set(k, v * (1.0 + dir * step));
                            let (e, neff) = energia_pesata(z, s, &q, &conf, &l0);
                            if e < best && neff > 0.5 * n as f64 { best = e; p = q; } else { break; }
                        }
                    }
                }
            }
        }
        if log { eprintln!("  VMC {:?} giro {}: E {:.4} -> {:.4}  {:?}", s, giro, e_start, best, p); }
    }
    // final plain VMC energy at the optimum
    let mut d = Dmc::new(Guida { z, stato: s, p }, 0.01, 2000, seed + 999);
    d.vmc_passi(300);
    let mut es = Vec::new();
    for _ in 0..600 { es.push(d.muovi(false)); }
    let (e, err) = media_errore(&es);
    (p, e, err)
}

/// Run fixed-node DMC: `eq` ticks to equilibrate, then `k` ticks measured. Returns (E, error, balls).
pub fn misura_dmc(g: Guida, tau: f64, n0: usize, eq: usize, k: usize, seed: u64) -> (f64, f64, Dmc) {
    let mut d = Dmc::new(g, tau, n0, seed);
    for _ in 0..eq { d.passo(); }
    let mut es = Vec::with_capacity(k);
    for _ in 0..k { es.push(d.passo()); }
    let (e, err) = media_errore(&es);
    (e, err, d)
}

/// Radial histogram of all the electrons of all the balls: r² ρ(r) on bins of width `dr` up to `rmax`.
pub fn istogramma_radiale(d: &Dmc, dr: f64, rmax: f64) -> Vec<f64> {
    let nb = (rmax / dr) as usize;
    let mut h = vec![0.0; nb];
    for e in 0..d.g.ne() {
        for p in d.elettrone(e) {
            let r = norm3(&p);
            let i = (r / dr) as usize;
            if i < nb { h[i] += 1.0; }
        }
    }
    h
}

/// All the shell boundaries found by the balls: the minima of the (smoothed) radial distribution between
/// consecutive clear maxima. Only radial distances enter.
pub fn confini(h: &[f64], dr: f64) -> Vec<f64> {
    let n = h.len();
    if n < 20 { return vec![]; }
    let s: Vec<f64> = (0..n).map(|i| {
        let a = i.saturating_sub(3);
        let b = (i + 4).min(n);
        h[a..b].iter().sum::<f64>() / (b - a) as f64
    }).collect();
    let top = s.iter().cloned().fold(0.0, f64::max);
    if top <= 0.0 { return vec![]; }
    let w = 10usize;
    let maxima: Vec<usize> = (1..n - 1).filter(|&i| {
        let a = i.saturating_sub(w);
        let b = (i + w + 1).min(n);
        s[i] > 0.02 * top && s[a..b].iter().all(|&v| v <= s[i])
    }).collect();
    let mut out = vec![];
    for p in maxima.windows(2) {
        let (a, b) = (p[0], p[1]);
        let imin = (a..=b).min_by(|&i, &j| s[i].partial_cmp(&s[j]).unwrap()).unwrap();
        // a real valley: clearly below both neighbouring maxima
        if s[imin] < 0.8 * s[a].min(s[b]) { out.push((imin as f64 + 0.5) * dr); }
    }
    out
}

/// Shell boundary found by the balls: the minimum of the radial distribution between its first two
/// maxima (K shell / L shell). Nothing is labelled: only the radial distance of each ball enters.
pub fn confine_guscio(h: &[f64], dr: f64) -> Option<f64> {
    // light smoothing
    let s: Vec<f64> = (0..h.len()).map(|i| {
        let a = i.saturating_sub(2);
        let b = (i + 3).min(h.len());
        h[a..b].iter().sum::<f64>() / (b - a) as f64
    }).collect();
    let imax1 = (1..s.len() - 1).find(|&i| s[i] >= s[i - 1] && s[i] > s[i + 1] && s[i] > 0.0)?;
    let mut imin = imax1;
    for i in imax1..s.len() - 1 {
        if s[i] < s[imin] { imin = i; }
        if s[i] > 1.5 * s[imin] && i > imin + 2 { return Some((imin as f64 + 0.5) * dr); }
    }
    None
}

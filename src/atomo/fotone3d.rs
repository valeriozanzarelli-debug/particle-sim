//! ONE FREE PHOTON IN 3-D WITH THE ZIG-ZAG RULE (spin-1 balls, F = E + icB, 3 complex numbers per cell).
//!
//! The same rule as `Light3D` in src/zigzag/walk.rs, written for speed (f32, one fused pass per axis, two
//! buffers): one tick = translate along x, y, z in turn (order x,y,z on even ticks, z,y,x on odd ticks); along
//! axis j the part of F on the +1 eigenvector of the spin matrix S_j moves one cell forward, the −1 part one
//! cell back, the 0 part stays. No mass, no flip: this is light. Every step is unitary, so the tick is undone
//! exactly (`indietro`). The box is a periodic cube; the photon is launched far from the faces and relaunched
//! before it reaches them.
//!
//! The launch: a Gaussian pulse (waist w0 at the focus, length σ along the flight) built mode by mode from
//! Maxwell plane waves: for every wave vector q, E ⊥ q and cB = q̂ × E, so the pulse has no static part and
//! no wave going backwards. The polarization ε (horizontal, vertical, circular ±) is the same for all modes
//! (projected ⊥ q). F = ½[(ε + i q̂×ε) e^{i(q·x−ωt)} + (ε* + i q̂×ε*) e^{−i(q·x−ωt)}]: for ε = (e₁ + i e₂)/√2
//! the second part vanishes (only positive frequencies: helicity +1).

use rayon::prelude::*;

use crate::cplx::C;
use crate::walk::order;

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct Cf { pub re: f32, pub im: f32 }

impl Cf {
    pub const ZERO: Cf = Cf { re: 0.0, im: 0.0 };
    #[inline] pub const fn new(re: f32, im: f32) -> Self { Cf { re, im } }
    #[inline] pub fn norm2(self) -> f32 { self.re * self.re + self.im * self.im }
}

const S2: f32 = std::f32::consts::FRAC_1_SQRT_2;

/// Polarization of the launched photon.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Polar { Orizzontale, Verticale, CircolarePiu, CircolareMeno }

impl Polar {
    pub fn nome(self) -> &'static str {
        match self {
            Polar::Orizzontale => "orizzontale",
            Polar::Verticale => "verticale",
            Polar::CircolarePiu => "circolare + (elicità +1)",
            Polar::CircolareMeno => "circolare − (elicità −1)",
        }
    }
    /// Expected helicity.
    pub fn elicita(self) -> f64 { match self { Polar::CircolarePiu => 1.0, Polar::CircolareMeno => -1.0, _ => 0.0 } }
}

/// What is launched (lengths in cells, times in ticks; c = 1 cell per tick).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lancio {
    pub lambda: f64,
    /// waist (1/e² radius of the intensity) at the focus
    pub w0: f64,
    /// length of the pulse: amplitude envelope exp(−s²/2σ²) along the flight
    pub sigma: f64,
    /// direction of flight (unit vector)
    pub dir: [f64; 3],
    pub polar: Polar,
    /// ticks from the launch to the focus (the focus is the centre of the box)
    pub t_fuoco: f64,
}

pub fn norma(v: [f64; 3]) -> [f64; 3] { let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt(); [v[0] / n, v[1] / n, v[2] / n] }
pub fn croce(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]] }
pub fn punto(a: [f64; 3], b: [f64; 3]) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }

/// Transverse basis of a direction: e₁ horizontal (⊥ to the vertical z and to d), e₂ = d × e₁; e₁ × e₂ = d.
pub fn base_trasversa(d: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let e1 = if d[2].abs() > 0.999 { [1.0, 0.0, 0.0] } else { norma(croce([0.0, 0.0, 1.0], d)) };
    (e1, croce(d, e1))
}

pub struct Fotone3D {
    pub n: usize,
    pub f: Vec<[Cf; 3]>,
    buf: Vec<[Cf; 3]>,
    pub ticks: u64,
}

/// Measured from the field: energy, centre of energy and its spread (cells, relative to the box centre).
#[derive(Clone, Copy, Debug, Default)]
pub struct Misura {
    #[allow(dead_code)]
    pub energia: f64,
    pub centro: [f64; 3],
    /// covariance of the energy distribution
    pub cov: [[f64; 3]; 3],
}

impl Misura {
    /// ⟨(δ·u)²⟩ along a unit vector.
    pub fn varianza(&self, u: [f64; 3]) -> f64 {
        let mut s = 0.0;
        for a in 0..3 { for b in 0..3 { s += u[a] * self.cov[a][b] * u[b]; } }
        s
    }
}

impl Fotone3D {
    pub fn new(n: usize) -> Self { Fotone3D { n, f: vec![[Cf::ZERO; 3]; n * n * n], buf: vec![[Cf::ZERO; 3]; n * n * n], ticks: 0 } }

    /// Cell coordinate relative to the box centre.
    #[inline] pub fn x(&self, i: usize) -> f64 { i as f64 - (self.n / 2) as f64 }

    /// One translation along `ax`; `verso` = +1 the rule, −1 its exact inverse.
    ///
    /// With a = F_{j+1}, b = F_{j+2} (j = ax): the +1 part (a − ib)/√2 comes from the cell behind (m), the −1
    /// part (a + ib)/√2 from the cell ahead (p), the 0 part F_j stays. Written out, the new components are
    /// a' = ½(a_m + a_p) + ½ i (b_p − b_m),  b' = ½(b_m + b_p) + ½ i (a_m − a_p),  F_j' = F_j.
    fn shift(&mut self, ax: usize, verso: i32) {
        let n = self.n;
        let (j1, j2) = ((ax + 1) % 3, (ax + 2) % 3);
        let src = &self.f;
        #[inline(always)]
        fn nuovo(q: &[Cf; 3], m: &[Cf; 3], p: &[Cf; 3], ax: usize, j1: usize, j2: usize) -> [Cf; 3] {
            let (am, bm, ap, bp) = (m[j1], m[j2], p[j1], p[j2]);
            let mut o = [Cf::ZERO; 3];
            o[ax] = q[ax];
            o[j1] = Cf::new(0.5 * (am.re + ap.re) - 0.5 * (bp.im - bm.im), 0.5 * (am.im + ap.im) + 0.5 * (bp.re - bm.re));
            o[j2] = Cf::new(0.5 * (bm.re + bp.re) - 0.5 * (am.im - ap.im), 0.5 * (bm.im + bp.im) + 0.5 * (am.re - ap.re));
            o
        }
        self.buf.par_chunks_mut(n).enumerate().for_each(|(riga, out)| {
            let (i, j) = (riga / n, riga % n);
            let rs = &src[riga * n..(riga + 1) * n];
            if ax == 2 {
                for k in 0..n {
                    let (km, kp) = (if k == 0 { n - 1 } else { k - 1 }, if k + 1 == n { 0 } else { k + 1 });
                    let (m, p) = if verso > 0 { (&rs[km], &rs[kp]) } else { (&rs[kp], &rs[km]) };
                    out[k] = nuovo(&rs[k], m, p, ax, j1, j2);
                }
            } else {
                let (rm, rp) = if ax == 0 {
                    (((i + n - 1) % n) * n + j, ((i + 1) % n) * n + j)
                } else {
                    (i * n + (j + n - 1) % n, i * n + (j + 1) % n)
                };
                let (rm, rp) = (&src[rm * n..(rm + 1) * n], &src[rp * n..(rp + 1) * n]);
                let (rm, rp) = if verso > 0 { (rm, rp) } else { (rp, rm) };
                for k in 0..n { out[k] = nuovo(&rs[k], &rm[k], &rp[k], ax, j1, j2); }
            }
        });
        std::mem::swap(&mut self.f, &mut self.buf);
    }

    pub fn tick(&mut self) {
        for ax in order(self.ticks) { self.shift(ax, 1); }
        self.ticks += 1;
    }

    /// The exact inverse of the last tick.
    pub fn indietro(&mut self) {
        if self.ticks == 0 { return; }
        self.ticks -= 1;
        let o = order(self.ticks);
        for ax in o.iter().rev() { self.shift(*ax, -1); }
    }

    /// Build the pulse of `l` (focus at the box centre, `t_fuoco` ticks after the launch); total Σ|F|² = 1.
    pub fn lancia(&mut self, l: &Lancio) {
        let n = self.n;
        let nf = n as f64;
        let tau = std::f64::consts::TAU;
        let d = norma(l.dir);
        let (e1, e2) = base_trasversa(d);
        // ε = a e₁ + b e₂ (complex)
        let (a, b) = match l.polar {
            Polar::Orizzontale => (C::real(1.0), C::ZERO),
            Polar::Verticale => (C::ZERO, C::real(1.0)),
            Polar::CircolarePiu => (C::real(S2 as f64), C::new(0.0, S2 as f64)),
            Polar::CircolareMeno => (C::real(S2 as f64), C::new(0.0, -S2 as f64)),
        };
        let eps: [C; 3] = [0, 1, 2].map(|c| a.scale(e1[c]) + b.scale(e2[c]));
        let k0 = tau / l.lambda;
        // the modes kept: A(q) = exp(−q⊥² w0²/4 − (q∥ − k0)² σ²/2), cut where it is below ~e⁻¹⁸
        let r = (8.5 / l.w0).max(6.0 / l.sigma).min(std::f64::consts::PI);
        let mr = (r * nf / tau).ceil() as i64;
        let m0: [i64; 3] = [0, 1, 2].map(|c| (k0 * d[c] * nf / tau).round() as i64);
        let nm = (2 * mr + 1) as usize;
        let qv = |c: usize, im: usize| tau * (m0[c] - mr + im as i64) as f64 / nf;
        // coefficients: 6 complex per mode (positive part, and the conjugate of the negative part)
        let mut coef = vec![[C::ZERO; 6]; nm * nm * nm];
        coef.par_iter_mut().enumerate().for_each(|(idx, cf)| {
            let (ia, ib, ic) = (idx / (nm * nm), (idx / nm) % nm, idx % nm);
            let q = [qv(0, ia), qv(1, ib), qv(2, ic)];
            let qn = punto(q, q).sqrt();
            if qn < 1e-9 { return; }
            let qh = [q[0] / qn, q[1] / qn, q[2] / qn];
            let qpar = punto(q, d);
            let qperp2 = (qn * qn - qpar * qpar).max(0.0);
            let amp = (-qperp2 * l.w0 * l.w0 / 4.0 - (qpar - k0).powi(2) * l.sigma * l.sigma / 2.0).exp();
            if amp < 1e-8 { return; }
            // ε ⊥ q
            let qe = eps[0].scale(qh[0]) + eps[1].scale(qh[1]) + eps[2].scale(qh[2]);
            let ep: [C; 3] = [0, 1, 2].map(|c| eps[c] - qe.scale(qh[c]));
            // q̂ × ε⊥
            let cr: [C; 3] = [0, 1, 2].map(|c| {
                let (c1, c2) = ((c + 1) % 3, (c + 2) % 3);
                ep[c2].scale(qh[c1]) - ep[c1].scale(qh[c2])
            });
            // positive part (ε + i q̂×ε) e^{+i|q| t_f}; negative part, conjugated: (ε − i q̂×ε) e^{+i|q| t_f}
            let ph = C::expi(qn * l.t_fuoco).scale(0.5 * amp);
            for c in 0..3 {
                cf[c] = (ep[c] + C::I * cr[c]) * ph;
                cf[3 + c] = (ep[c] - C::I * cr[c]) * ph;
            }
        });
        // separable inverse DFT on the kept modes: z, then y, then x
        let ex = |c: usize| -> Vec<C> {
            let mut v = vec![C::ZERO; nm * n];
            for im in 0..nm { for i in 0..n { v[im * n + i] = C::expi(qv(c, im) * self.x(i)); } }
            v
        };
        let (ezx, eyx, exx) = (ex(2), ex(1), ex(0));
        let mut t1 = vec![[C::ZERO; 6]; nm * nm * n];
        t1.par_chunks_mut(n).enumerate().for_each(|(ab, out)| {
            for ic in 0..nm {
                let cf = &coef[ab * nm + ic];
                if cf.iter().all(|z| *z == C::ZERO) { continue; }
                for (k, o) in out.iter_mut().enumerate() {
                    let e = ezx[ic * n + k];
                    for q in 0..6 { o[q] += cf[q] * e; }
                }
            }
        });
        let mut t2 = vec![[C::ZERO; 6]; nm * n * n];
        t2.par_chunks_mut(n * n).enumerate().for_each(|(ia, out)| {
            for ib in 0..nm {
                let row = &t1[(ia * nm + ib) * n..(ia * nm + ib + 1) * n];
                for j in 0..n {
                    let e = eyx[ib * n + j];
                    for k in 0..n { for q in 0..6 { out[j * n + k][q] += row[k][q] * e; } }
                }
            }
        });
        self.f.par_chunks_mut(n * n).enumerate().for_each(|(i, out)| {
            let mut acc = vec![[C::ZERO; 6]; n * n];
            for ia in 0..nm {
                let e = exx[ia * n + i];
                let pl = &t2[ia * n * n..(ia + 1) * n * n];
                for (o, p) in acc.iter_mut().zip(pl) { for q in 0..6 { o[q] += p[q] * e; } }
            }
            for (o, v) in out.iter_mut().zip(&acc) {
                for c in 0..3 {
                    let z = v[c] + v[3 + c].conj();
                    o[c] = Cf::new(z.re as f32, z.im as f32);
                }
            }
        });
        // normalise: Σ|F|² = 1
        let e: f64 = self.f.par_iter().map(|p| p.iter().map(|c| c.norm2() as f64).sum::<f64>()).sum();
        let k = (1.0 / e.max(1e-300)).sqrt() as f32;
        self.f.par_iter_mut().for_each(|p| for c in p.iter_mut() { c.re *= k; c.im *= k; });
        self.ticks = 0;
    }

    /// Energy Σ|F|², its centre and covariance; displacements taken as the nearest image of `rif`.
    pub fn misura(&self, rif: [f64; 3]) -> Misura {
        let n = self.n;
        let nf = n as f64;
        let wrap = |v: f64| v - nf * (v / nf).round();
        let s: [f64; 10] = self.f.par_chunks(n * n).enumerate().map(|(i, pl)| {
            let mut a = [0.0f64; 10];
            let dx = wrap(self.x(i) - rif[0]);
            for j in 0..n {
                let dy = wrap(self.x(j) - rif[1]);
                for k in 0..n {
                    let p = &pl[j * n + k];
                    let w = (p[0].norm2() + p[1].norm2() + p[2].norm2()) as f64;
                    if w == 0.0 { continue; }
                    let dz = wrap(self.x(k) - rif[2]);
                    a[0] += w;
                    a[1] += w * dx; a[2] += w * dy; a[3] += w * dz;
                    a[4] += w * dx * dx; a[5] += w * dy * dy; a[6] += w * dz * dz;
                    a[7] += w * dx * dy; a[8] += w * dx * dz; a[9] += w * dy * dz;
                }
            }
            a
        }).reduce(|| [0.0; 10], |mut a, b| { for q in 0..10 { a[q] += b[q]; } a });
        let w = s[0].max(1e-300);
        let m = [s[1] / w, s[2] / w, s[3] / w];
        let c2 = |q: usize, a: usize, b: usize| s[q] / w - m[a] * m[b];
        let cov = [[c2(4, 0, 0), c2(7, 0, 1), c2(8, 0, 2)], [c2(7, 0, 1), c2(5, 1, 1), c2(9, 1, 2)], [c2(8, 0, 2), c2(9, 1, 2), c2(6, 2, 2)]];
        Misura { energia: s[0], centro: [rif[0] + m[0], rif[1] + m[1], rif[2] + m[2]], cov }
    }

    /// The electric field E = Re F per cell.
    pub fn campo_e(&self) -> Vec<[f32; 3]> { self.f.par_iter().map(|p| [p[0].re, p[1].re, p[2].re]).collect() }

    /// Helicity from how E turns in one tick: Σ (E_prima × E_ora)·d / (sin ω · Σ|E⊥|²), with ω = 2π/λ.
    /// +1: E turns counter-clockwise about the flight (spin along d), −1 clockwise, 0 linear.
    pub fn elicita(&self, prima: &[[f32; 3]], d: [f64; 3], lambda: f64) -> f64 {
        let (e1, e2) = base_trasversa(d);
        let (e1, e2) = ([e1[0] as f32, e1[1] as f32, e1[2] as f32], [e2[0] as f32, e2[1] as f32, e2[2] as f32]);
        let pr = |e: &[f32; 3], u: &[f32; 3]| e[0] * u[0] + e[1] * u[1] + e[2] * u[2];
        let (x, w) = self.f.par_iter().zip(prima.par_iter()).map(|(p, a)| {
            let b = [p[0].re, p[1].re, p[2].re];
            let (a1, a2, b1, b2) = (pr(a, &e1), pr(a, &e2), pr(&b, &e1), pr(&b, &e2));
            ((a1 * b2 - a2 * b1) as f64, (0.5 * (a1 * a1 + a2 * a2 + b1 * b1 + b2 * b2)) as f64)
        }).reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
        x / ((std::f64::consts::TAU / lambda).sin() * w.max(1e-300))
    }

    /// Exact continuum prediction of the speed of the energy centre for the launched spectrum:
    /// c ⟨q̂·d⟩ weighted with |A(q)|² (a focused pulse is a little slower than c along its axis).
    pub fn velocita_continuo(l: &Lancio, n: usize) -> f64 {
        let tau = std::f64::consts::TAU;
        let nf = n as f64;
        let d = norma(l.dir);
        let k0 = tau / l.lambda;
        let mr = 40i64;
        let (mut s, mut w) = (0.0, 0.0);
        for a in -mr..=mr { for b in -mr..=mr { for c in -mr..=mr {
            let q = [k0 * d[0] + tau * a as f64 / nf, k0 * d[1] + tau * b as f64 / nf, k0 * d[2] + tau * c as f64 / nf];
            let qn = punto(q, q).sqrt();
            let qpar = punto(q, d);
            let amp2 = (2.0 * (-(qn * qn - qpar * qpar).max(0.0) * l.w0 * l.w0 / 4.0 - (qpar - k0).powi(2) * l.sigma * l.sigma / 2.0)).exp();
            // (the mode q = 0 has no direction and no weight: skipped, it gave 0/0)
            if qn < 1e-12 { continue; }
            s += amp2 * qpar / qn;
            w += amp2;
        }}}
        s / w
    }
}

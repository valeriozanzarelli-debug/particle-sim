//! THE ZIG-ZAG ATOM: one electron bound to a nucleus by the zig-zag rule (the same rule as src/zigzag,
//! docs/ZIGZAG.md), run LIVE in 3-D, driven by the field of a photon, exactly reversible in time.
//!
//! One tick (length a/c, units ħ = c = mₑ = 1, lengths in λ_C = ħ/mc):
//!   1. TRANSLATE along x, y, z in turn (order alternated): the right-handed part moves along its spin, the
//!      left-handed part against it, one cell (= speed c);
//!   2. ROTATE: the mass flips L ↔ R by θ = a (the electron mass), and the phase turns by V·a, where
//!      V = −Zα/r (Coulomb phase of a point nucleus) + e E(t)·r (the photon's electric field at the atom,
//!      dipole form: correct because the photon's wavelength is ~70 times the atom).
//! Every operation is unitary, so the rule runs BACKWARDS exactly (`indietro`): going back in time is the
//! same rule undone, not a recording.
//!
//! Why not hydrogen: hydrogen's Bohr radius is 137 λ_C and the cells must be ≲ 0.4 λ_C: (5·137/0.4)³ ≈ 5·10⁹
//! cells. The Dirac–Coulomb problem depends only on Zα, so we take the hydrogen-like ion with Z = 60
//! (neodymium, Nd⁵⁹⁺, Zα = 0.4378, Bohr radius 2.28 λ_C ≈ 880 fm): the biggest Z·α that this PC can hold
//! with 2p inside the box (docs/ZIGZAG.md §5.6 uses the same ion).

use rayon::prelude::*;

use crate::rng::Rng;
use crate::cplx::C;
use crate::walk::{order, surfer_velocity, HALF, Radial};

pub const Z_ION: f64 = 60.0;
pub const ALPHA: f64 = 7.297_352_5693e-3;
pub fn z_alpha() -> f64 { Z_ION * ALPHA }
/// cells per side, cell size (λ_C)
pub const N: usize = 48;
pub const A: f64 = 0.5;
/// metres per λ_C and seconds per (λ_C / c)
pub const LAMBDA_C_M: f64 = 3.861_592_6796e-13;
pub const T_UNIT_S: f64 = 1.288_088_668_19e-21;
pub const MC2_KEV: f64 = 510.998_950;

#[inline]
fn wrap(i: usize, d: isize, n: usize) -> usize { ((i as isize + d).rem_euclid(n as isize)) as usize }
#[inline]
fn nb(s: usize, j: usize, d: isize, n: usize) -> usize {
    let (i, k, l) = (s / (n * n), (s / n) % n, s % n);
    match j {
        0 => (wrap(i, d, n) * n + k) * n + l,
        1 => (i * n + wrap(k, d, n)) * n + l,
        _ => (i * n + k) * n + wrap(l, d, n),
    }
}

/// The photon's drive: a pulse E(t) = E₀ e^{−((t−t₀)/σ)²} cos(ω (t − t₀)) along z (units mc²/(e λ_C)).
#[derive(Clone, Copy, Debug)]
pub struct Impulso { pub e0: f64, pub omega: f64, pub t0: f64, pub sigma: f64 }

impl Impulso {
    pub fn campo(&self, t: f64) -> f64 {
        let x = (t - self.t0) / self.sigma;
        if x.abs() > 4.0 { return 0.0; }
        self.e0 * (-x * x).exp() * (self.omega * (t - self.t0)).cos()
    }
}

#[derive(Clone)]
pub struct Atomo {
    /// (L₀, L₁, R₀, R₁) per site
    pub psi: Vec<[C; 4]>,
    tmp: Vec<[C; 4]>,
    /// Coulomb potential at each site and z of each site
    pc: Vec<C>,
    zs: Vec<f64>,
    pub ticks: u64,
    pub impulsi: Vec<Impulso>,
    /// field of the live light at the atom, already times the charge: e·E (units mc²/(e λ_C)·e), any direction
    pub campo_est: [f64; 3],
    /// the 1s state the atom started in (for the population of the ground state)
    pub psi0: Vec<[C; 4]>,
}

pub fn x_of(i: usize) -> f64 { (i as f64 - N as f64 / 2.0 + 0.5) * A }

impl Atomo {
    /// The atom in its ground state: the 1s (m = +½) wave of the radial form of the same rule (same a,
    /// same Zα), put on the cubic lattice.
    pub fn new() -> Self { Self::con_stato(-1, 1) }

    /// The atom started in level (κ, n) of the radial rule (κ = −1: s½; κ = +1: p½), m = +½.
    pub fn con_stato(kappa: i32, nn: u32) -> Self {
        let n = N;
        let za = z_alpha();
        let lv = livello(kappa, nn);
        let mut psi = vec![[C::ZERO; 4]; n * n * n];
        let mut vc = vec![0.0; n * n * n];
        let mut zs = vec![0.0; n * n * n];
        psi.par_iter_mut().zip(vc.par_iter_mut()).zip(zs.par_iter_mut()).enumerate().for_each(|(s, ((p, v), z))| {
            let x = [x_of(s / (n * n)), x_of((s / n) % n), x_of(s % n)];
            let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            *v = -za / r;
            *z = x[2];
            let d = spinore(&lv, kappa, x);
            *p = dirac_a_chirale(d);
        });
        let norm: f64 = psi.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum();
        let k = 1.0 / norm.sqrt();
        psi.par_iter_mut().for_each(|p| for c in p.iter_mut() { *c = c.scale(k); });
        let psi0 = psi.clone();
        let pc = vc.iter().map(|v| C::expi(-v * A)).collect();
        Atomo { tmp: vec![[C::ZERO; 4]; n * n * n], psi, pc, zs, ticks: 0, impulsi: vec![], campo_est: [0.0; 3], psi0 }
    }

    pub fn t(&self) -> f64 { self.ticks as f64 * A }
    pub fn campo(&self, t: f64) -> f64 { self.impulsi.iter().map(|i| i.campo(t)).sum() }

    fn shift(&mut self, j: usize, verso: isize) {
        let n = N;
        let [u, d] = HALF[j];
        self.tmp.par_iter_mut().zip(self.psi.par_iter()).for_each(|(t, p)| {
            t[0] = u[0].conj() * p[0] + u[1].conj() * p[1];
            t[1] = d[0].conj() * p[0] + d[1].conj() * p[1];
            t[2] = u[0].conj() * p[2] + u[1].conj() * p[3];
            t[3] = d[0].conj() * p[2] + d[1].conj() * p[3];
        });
        let tmp = &self.tmp;
        self.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
            // forward (verso = +1): L u₊ moves −1, L u₋ +1, R u₊ +1, R u₋ −1; backward: the opposite
            let lp = tmp[nb(s, j, verso, n)][0];
            let lm = tmp[nb(s, j, -verso, n)][1];
            let rp = tmp[nb(s, j, -verso, n)][2];
            let rm = tmp[nb(s, j, verso, n)][3];
            p[0] = u[0] * lp + d[0] * lm;
            p[1] = u[1] * lp + d[1] * lm;
            p[2] = u[0] * rp + d[0] * rm;
            p[3] = u[1] * rp + d[1] * rm;
        });
    }

    fn ruota(&mut self, t: f64, verso: f64) {
        let (ct, st) = (A.cos(), A.sin());
        let e = self.campo(t);
        // phases: Coulomb part per site (fixed), field part per z-slice (N numbers per tick)
        // V = e E·r: the prescribed pulse (along z) plus the live light at the atom (any direction)
        let ce = self.campo_est;
        let pz: Vec<C> = (0..N).map(|k| C::expi(-verso * (e + ce[2]) * x_of(k) * A)).collect();
        let px: Vec<C> = (0..N).map(|k| C::expi(-verso * ce[0] * x_of(k) * A)).collect();
        let py: Vec<C> = (0..N).map(|k| C::expi(-verso * ce[1] * x_of(k) * A)).collect();
        let trasverso = ce[0] != 0.0 || ce[1] != 0.0;
        let pc = &self.pc;
        let mi = C::new(0.0, -verso * st);
        self.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
            let c0 = if verso > 0.0 { pc[s] } else { pc[s].conj() };
            let ph = if trasverso { c0 * pz[s % N] * px[s / (N * N)] * py[(s / N) % N] } else { c0 * pz[s % N] };
            for q in 0..2 {
                let (l, r) = (p[q], p[2 + q]);
                if verso > 0.0 {
                    p[q] = ph * (l * ct + mi * r);
                    p[2 + q] = ph * (r * ct + mi * l);
                } else {
                    // undo: first the phase, then the inverse flip
                    let (l, r) = (ph * l, ph * r);
                    p[q] = l * ct + mi * r;
                    p[2 + q] = r * ct + mi * l;
                }
            }
        });
    }

    /// One tick forward.
    pub fn avanti(&mut self) {
        let t = self.t();
        for j in order(self.ticks) { self.shift(j, 1); }
        self.ruota(t, 1.0);
        self.ticks += 1;
    }

    /// One tick backward: the last tick undone exactly.
    pub fn indietro(&mut self) {
        if self.ticks == 0 { return; }
        self.ticks -= 1;
        let t = self.t();
        self.ruota(t, -1.0);
        let o = order(self.ticks);
        for &j in o.iter().rev() { self.shift(j, -1); }
    }

    /// Spinor at a point (λ_C), trilinear.
    pub fn campiona(&self, p: [f64; 3]) -> [C; 4] {
        let n = N;
        let f = |x: f64| (x / A + n as f64 / 2.0 - 0.5).rem_euclid(n as f64);
        let (fx, fy, fz) = (f(p[0]), f(p[1]), f(p[2]));
        let (i0, j0, k0) = (fx.floor() as usize % n, fy.floor() as usize % n, fz.floor() as usize % n);
        let (tx, ty, tz) = (fx - fx.floor(), fy - fy.floor(), fz - fz.floor());
        let mut out = [C::ZERO; 4];
        for (di, wx) in [(0usize, 1.0 - tx), (1, tx)] {
            for (dj, wy) in [(0usize, 1.0 - ty), (1, ty)] {
                for (dk, wz) in [(0usize, 1.0 - tz), (1, tz)] {
                    let s = (((i0 + di) % n) * n + (j0 + dj) % n) * n + (k0 + dk) % n;
                    let w = wx * wy * wz;
                    for q in 0..4 { out[q] += self.psi[s][q] * w; }
                }
            }
        }
        out
    }

    /// Position measurement with collapse: the wave becomes a narrow packet (width 2 cells) at `p`, with the
    /// spinor the wave had there; the rule carries on from it (it spreads again).
    pub fn collassa(&mut self, p: [f64; 3]) {
        let n = N;
        let mut chi = self.campiona(p);
        let m: f64 = chi.iter().map(|c| c.norm2()).sum::<f64>().sqrt();
        if m < 1e-300 { chi = [C::ZERO, C::ZERO, C::ONE, C::ZERO]; } else { for c in chi.iter_mut() { *c = c.scale(1.0 / m); } }
        let sig = 2.0 * A;
        let lato = n as f64 * A;
        let wrap = |d: f64| d - lato * (d / lato).round();
        self.psi.par_iter_mut().enumerate().for_each(|(s, q)| {
            let x = [x_of(s / (n * n)), x_of((s / n) % n), x_of(s % n)];
            let d2: f64 = (0..3).map(|c| wrap(x[c] - p[c]).powi(2)).sum();
            let g = (-d2 / (2.0 * sig * sig)).exp();
            for c in 0..4 { q[c] = chi[c].scale(g); }
        });
        let norm: f64 = self.psi.par_iter().map(|q| q.iter().map(|c| c.norm2()).sum::<f64>()).sum();
        let k = 1.0 / norm.max(1e-300).sqrt();
        self.psi.par_iter_mut().for_each(|q| for c in q.iter_mut() { *c = c.scale(k); });
    }

    /// ⟨α⟩ = Σ (R†σR − L†σL): the probability current of the walk summed over the atom (units of c);
    /// the electric current is −e⟨α⟩.
    pub fn corrente(&self) -> [f64; 3] {
        self.psi.par_iter().map(|p| surfer_velocity_num(p)).reduce(|| [0.0; 3], |a, b| [a[0] + b[0], a[1] + b[1], a[2] + b[2]])
    }

    /// ⟨ψ|φ⟩ projection probability |⟨φ|ψ⟩|² onto a state φ of the lattice.
    pub fn proiezione(&self, phi: &[[C; 4]]) -> f64 {
        let s = self.psi.par_iter().zip(phi.par_iter()).map(|(a, b)| {
            let mut acc = C::ZERO;
            for q in 0..4 { acc += b[q].conj() * a[q]; }
            acc
        }).reduce(|| C::ZERO, |a, b| a + b);
        s.norm2()
    }

    /// A level of THIS 3-D lattice rule, by listening. The start state (no light) is evolved for `ticks` ticks;
    /// meanwhile ⟨ψ₀|ψ(t)⟩ is recorded and Σ ψ(t) e^{iωt} (Hann window) is kept at frequencies every 0.0025 in
    /// [lo, hi]. The spectrum of ⟨ψ₀|ψ(t)⟩ shows the levels the start state contains (for the radial 1s put on the
    /// cubic lattice: 0.901 (60 %), 0.8405 (26 %), 0.882 (13 %)); the strongest peak in [lo, hi] is taken and its
    /// filtered wave returned with its frequency (parabolic peak). Long windows separate close levels.
    pub fn ascolta(&self, lo: f64, hi: f64, ticks: usize) -> (Vec<[C; 4]>, f64) {
        let mut a = self.clone();
        a.impulsi.clear();
        a.campo_est = [0.0; 3];
        let p0 = a.psi.clone();
        let nc = ((hi - lo) / 0.0025).round() as usize + 1;
        let cand: Vec<f64> = (0..nc).map(|k| lo + 0.0025 * k as f64).collect();
        let mut acc: Vec<Vec<[C; 4]>> = vec![vec![[C::ZERO; 4]; a.psi.len()]; cand.len()];
        let mut serie = Vec::with_capacity(ticks);
        let tt = ticks as f64 * A;
        for t in 0..ticks {
            let w = (std::f64::consts::PI * (t as f64 + 0.5) / ticks as f64).sin().powi(2);
            let time = a.t();
            let c = a.psi.par_iter().zip(p0.par_iter()).map(|(x, y)| { let mut s = C::ZERO; for q in 0..4 { s += y[q].conj() * x[q]; } s })
                .reduce(|| C::ZERO, |u, v| u + v);
            serie.push((time, c));
            let fasi: Vec<C> = cand.iter().map(|om| C::expi(om * time).scale(w)).collect();
            acc.par_iter_mut().zip(fasi.par_iter()).for_each(|(v, ph)| for (s, p) in v.iter_mut().zip(a.psi.iter()) { for q in 0..4 { s[q] += p[q] * *ph; } });
            a.avanti();
        }
        // spectrum of the overlap on a fine grid, strongest peak in [lo, hi]
        let spettro = |om: f64| -> f64 {
            let mut s = C::ZERO;
            for (t, v) in &serie { let h = (std::f64::consts::PI * t / tt).sin().powi(2); s += *v * C::expi(om * t).scale(h); }
            s.norm2()
        };
        let fine: Vec<(f64, f64)> = (0..=((hi - lo) / 0.0002) as usize).map(|k| { let om = lo + 0.0002 * k as f64; (om, spettro(om)) }).collect();
        let kb = (1..fine.len() - 1).max_by(|&i, &j| fine[i].1.total_cmp(&fine[j].1)).unwrap_or(0);
        let mut om = fine[kb].0;
        if kb > 0 && kb + 1 < fine.len() {
            let (y0, y1, y2) = (fine[kb - 1].1, fine[kb].1, fine[kb + 1].1);
            let den = y0 - 2.0 * y1 + y2;
            if den.abs() > 0.0 { om += 0.0002 * 0.5 * (y0 - y2) / den; }
        }
        let kc = (0..cand.len()).min_by(|&i, &j| (cand[i] - om).abs().total_cmp(&(cand[j] - om).abs())).unwrap();
        let mut v = std::mem::take(&mut acc[kc]);
        let nn: f64 = v.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum();
        let k = 1.0 / nn.max(1e-300).sqrt();
        v.par_iter_mut().for_each(|p| for c in p.iter_mut() { *c = c.scale(k); });
        (v, om)
    }

    /// Replace the wave (e.g. by the listened 1s); the start state follows.
    pub fn imposta(&mut self, psi: Vec<[C; 4]>) {
        self.psi = psi;
        self.psi0 = self.psi.clone();
        self.ticks = 0;
    }

    /// ⟨r⟩ (λ_C): the position of the electron, all three components.
    pub fn posizione(&self) -> [f64; 3] {
        let n = N;
        self.psi.par_iter().enumerate().map(|(s, p)| {
            let w: f64 = p.iter().map(|c| c.norm2()).sum();
            [w * x_of(s / (n * n)), w * x_of((s / n) % n), w * x_of(s % n)]
        }).reduce(|| [0.0; 3], |a, b| [a[0] + b[0], a[1] + b[1], a[2] + b[2]])
    }

    /// ⟨z⟩ (λ_C): the dipole of the atom (it radiates when it oscillates).
    pub fn dipolo(&self) -> f64 {
        self.psi.par_iter().zip(self.zs.par_iter()).map(|(p, z)| p.iter().map(|c| c.norm2()).sum::<f64>() * z).sum()
    }

    /// Probability of still being in the starting 1s state, |⟨1s|ψ⟩|².
    pub fn p_fondamentale(&self) -> f64 {
        let s = self.psi.par_iter().zip(self.psi0.par_iter()).map(|(a, b)| {
            let mut acc = C::ZERO;
            for q in 0..4 { acc += b[q].conj() * a[q]; }
            acc
        }).reduce(|| C::ZERO, |a, b| a + b);
        s.norm2()
    }

    /// Balls born with |ψ|² (random points, jittered inside their cell: never the lattice itself).
    pub fn nascite(&self, k: usize, rng: &mut Rng) -> Vec<[f64; 3]> {
        let w: Vec<f64> = self.psi.iter().map(|p| p.iter().map(|c| c.norm2()).sum()).collect();
        let mut cdf = Vec::with_capacity(w.len());
        let mut s = 0.0;
        for x in &w { s += x; cdf.push(s); }
        let n = N;
        (0..k).map(|_| {
            let u = rng.uni() * s;
            let i = cdf.partition_point(|&c| c < u).min(w.len() - 1);
            [x_of(i / (n * n)) + A * (rng.uni() - 0.5), x_of((i / n) % n) + A * (rng.uni() - 0.5), x_of(i % n) + A * (rng.uni() - 0.5)]
        }).collect()
    }

    /// Move the surfers one tick along v = j/ρ of the walk (midpoint rule; `verso` = −1 backward).
    pub fn muovi_surfisti(&self, s: &mut [[f64; 3]], verso: f64) {
        s.par_iter_mut().for_each(|p| {
            let v1 = surfer_velocity(&self.campiona(*p));
            let q = [p[0] + 0.5 * verso * v1[0] * A, p[1] + 0.5 * verso * v1[1] * A, p[2] + 0.5 * verso * v1[2] * A];
            let v2 = surfer_velocity(&self.campiona(q));
            for k in 0..3 { p[k] += verso * v2[k] * A; }
        });
    }

    #[allow(dead_code)]
    /// Phase of the dominant component at a point (for the colour).
    pub fn fase(&self, p: [f64; 3]) -> f32 {
        let s = self.campiona(p);
        let mut best = 0;
        for q in 1..4 { if s[q].norm2() > s[best].norm2() { best = q; } }
        s[best].arg() as f32
    }
}

/// Numerator of the surfer velocity: R†σR − L†σL (the probability current of one site, units of c).
pub fn surfer_velocity_num(p: &[C; 4]) -> [f64; 3] {
    let (l0, l1, r0, r1) = (p[0], p[1], p[2], p[3]);
    let sx = |a: C, b: C| 2.0 * (a.conj() * b).re;
    let sy = |a: C, b: C| 2.0 * (a.conj() * b).im;
    let sz = |a: C, b: C| a.norm2() - b.norm2();
    [sx(r0, r1) - sx(l0, l1), sy(r0, r1) - sy(l0, l1), sz(r0, r1) - sz(l0, l1)]
}

/// Dirac representation (φ₀, φ₁, χ₀, χ₁) → chiral (L₀, L₁, R₀, R₁): L = (φ − χ)/√2, R = (φ + χ)/√2.
pub fn dirac_a_chirale(d: [C; 4]) -> [C; 4] {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    [(d[0] - d[2]).scale(s), (d[1] - d[3]).scale(s), (d[0] + d[2]).scale(s), (d[1] + d[3]).scale(s)]
}

/// A level of partial wave κ (n-th ring) from the radial form of the same rule, on the same lattice.
pub struct Livello { pub omega: f64, pub gf: Vec<[f64; 2]> }

pub fn livello(kappa: i32, n: u32) -> Livello {
    let za = z_alpha();
    let rmax = 40.0;
    let w = Radial::new(kappa, A, rmax, za);
    // guess from the Bohr scale (only chooses which ring to listen to)
    let guess = 1.0 - 0.5 * za * za / (n * n) as f64 * 1.2;
    let (om, gf) = w.eigen(guess, 8);
    Livello { omega: om, gf }
}

/// Dirac spinor (Dirac representation) at p (λ_C) for radial (g, f), κ = −1 (s½) or +1 (p½), m = +½.
pub fn spinore(lv: &Livello, kappa: i32, p: [f64; 3]) -> [C; 4] {
    let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(1e-9);
    let x = r / A - 0.5;
    let i = (x.floor().max(0.0) as usize).min(lv.gf.len() - 2);
    let t = (x - i as f64).clamp(0.0, 1.0);
    let g = (lv.gf[i][0] * (1.0 - t) + lv.gf[i + 1][0] * t) / r;
    let f = (lv.gf[i][1] * (1.0 - t) + lv.gf[i + 1][1] * t) / r;
    let y00 = 0.5 / std::f64::consts::PI.sqrt();
    // spin-angle functions Ω for j = ½, m = +½: Ω_{κ=−1} = Y₀₀ (1, 0); Ω_{κ=+1} = −(σ·r̂) Ω_{κ=−1}
    let (cz, cx, cy) = (p[2] / r, p[0] / r, p[1] / r);
    let om_s = [C::real(y00), C::ZERO];
    let om_p = [C::real(-y00 * cz), C::new(-y00 * cx, -y00 * cy)];
    let (up, lo) = if kappa < 0 { (om_s, om_p) } else { (om_p, om_s) };
    [up[0].scale(g), up[1].scale(g), C::I * lo[0].scale(f), C::I * lo[1].scale(f)]
}

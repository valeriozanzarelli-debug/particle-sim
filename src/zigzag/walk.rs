//! The ZIG-ZAG rule (the same rule as `lab/zz_core.py`, see docs/ZIGZAG.md).
//!
//! Every wave ball moves at c and rotates its internal state. One tick (length a/c):
//! 1. TRANSLATE along x, y, z in turn (order x,y,z on even ticks, z,y,x on odd ticks): the ball's state
//!    is split on the eigenvectors of its spin matrix S_j; the +1 part moves one cell along +j, the −1
//!    part one cell along −j, the 0 part (light only) stays.
//!    spin ½ right-handed R: S_j = σ_j; left-handed L: S_j = −σ_j; spin 1 (light, F = E + icB):
//!    (S_k)_ij = −i ε_kij.
//! 2. ROTATE: L' = e^{−iVdt}(cos θ L − i sin θ R), R' = e^{−iVdt}(cos θ R − i sin θ L) with
//!    θ = m c a / ħ (the mass flip) and V the Coulomb phase of a point nucleus.
//! Units: ħ = c = mₑ = 1, lattice spacing a, tick = a.

use rayon::prelude::*;

use crate::cplx::C;

const S2: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// Eigenvectors (σ_j = +1, σ_j = −1) of the Pauli matrix σ_j: `HALF[j][0]` = u₊, `HALF[j][1]` = u₋.
pub const HALF: [[[C; 2]; 2]; 3] = [
    [[C::real(S2), C::real(S2)], [C::real(S2), C::real(-S2)]],
    [[C::real(S2), C::new(0.0, S2)], [C::real(S2), C::new(0.0, -S2)]],
    [[C::ONE, C::ZERO], [C::ZERO, C::ONE]],
];

/// Eigenvectors (+1, −1) of the spin-1 matrix S_j ((S_k)_ij = −i ε_kij): u± = (e_{j+1} ± i e_{j+2})/√2.
pub fn one_basis(j: usize) -> [[C; 3]; 2] {
    let (j1, j2) = ((j + 1) % 3, (j + 2) % 3);
    let mut up = [C::ZERO; 3];
    let mut dn = [C::ZERO; 3];
    up[j1] = C::real(S2);
    up[j2] = C::new(0.0, S2);
    dn[j1] = C::real(S2);
    dn[j2] = C::new(0.0, -S2);
    [up, dn]
}

/// Axis order of tick `n` (alternated: two ticks form a palindrome → second-order rule).
pub fn order(n: u64) -> [usize; 3] { if n % 2 == 0 { [0, 1, 2] } else { [2, 1, 0] } }

#[inline]
fn wrap(i: usize, d: isize, n: usize) -> usize { ((i as isize + d).rem_euclid(n as isize)) as usize }

/// Neighbour index of site `s` along axis `j` by `d` cells (periodic cube of side n).
#[inline]
fn nb(s: usize, j: usize, d: isize, n: usize) -> usize {
    let (i, k, l) = (s / (n * n), (s / n) % n, s % n);
    match j {
        0 => (wrap(i, d, n) * n + k) * n + l,
        1 => (i * n + wrap(k, d, n)) * n + l,
        _ => (i * n + k) * n + wrap(l, d, n),
    }
}

/// Electron (or any massive spin-½ particle) on a periodic cube: ψ = (L₀, L₁, R₀, R₁) per site.
pub struct Dirac3D {
    pub n: usize,
    pub a: f64,
    pub theta: f64,
    pub psi: Vec<[C; 4]>,
    tmp: Vec<[C; 4]>,
    phase: Option<Vec<C>>,
    pub ticks: u64,
}

impl Dirac3D {
    pub fn new(n: usize, a: f64, theta: f64) -> Self {
        Self { n, a, theta, psi: vec![[C::ZERO; 4]; n * n * n], tmp: vec![[C::ZERO; 4]; n * n * n], phase: None, ticks: 0 }
    }

    /// Static potential V (units m c²) at every site → phase e^{−iVa} per tick.
    pub fn set_potential(&mut self, v: &[f64]) {
        self.phase = Some(v.iter().map(|&x| C::expi(-x * self.a)).collect());
    }

    /// Coordinate of lattice index i (cube centred on the origin, the origin between sites).
    pub fn x(&self, i: usize) -> f64 { (i as f64 - self.n as f64 / 2.0 + 0.5) * self.a }

    fn shift(&mut self, j: usize) {
        let n = self.n;
        let [u, d] = HALF[j];
        self.tmp.par_iter_mut().zip(self.psi.par_iter()).for_each(|(t, p)| {
            // projections on u₊ and u₋ for the L block (0,1) and the R block (2,3)
            t[0] = u[0].conj() * p[0] + u[1].conj() * p[1];
            t[1] = d[0].conj() * p[0] + d[1].conj() * p[1];
            t[2] = u[0].conj() * p[2] + u[1].conj() * p[3];
            t[3] = d[0].conj() * p[2] + d[1].conj() * p[3];
        });
        let tmp = &self.tmp;
        self.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
            // L (h = −1): u₊ part moves −1 → read from s+1; u₋ part moves +1 → read from s−1
            let lp = tmp[nb(s, j, 1, n)][0];
            let lm = tmp[nb(s, j, -1, n)][1];
            // R (h = +1): u₊ part moves +1 → read from s−1; u₋ part moves −1 → read from s+1
            let rp = tmp[nb(s, j, -1, n)][2];
            let rm = tmp[nb(s, j, 1, n)][3];
            p[0] = u[0] * lp + d[0] * lm;
            p[1] = u[1] * lp + d[1] * lm;
            p[2] = u[0] * rp + d[0] * rm;
            p[3] = u[1] * rp + d[1] * rm;
        });
    }

    pub fn tick(&mut self) {
        for j in order(self.ticks) { self.shift(j); }
        let (ct, st) = (self.theta.cos(), self.theta.sin());
        let mi = C::new(0.0, -st);
        let ph = self.phase.as_deref();
        self.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
            let e = ph.map_or(C::ONE, |v| v[s]);
            for q in 0..2 {
                let (l, r) = (p[q], p[2 + q]);
                p[q] = e * (l * ct + mi * r);
                p[2 + q] = e * (r * ct + mi * l);
            }
        });
        self.ticks += 1;
    }

    pub fn norm2(&self) -> f64 { self.psi.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum() }

    /// Spinor at a point (trilinear interpolation of the components), for the surfers.
    pub fn sample(&self, p: [f64; 3]) -> [C; 4] {
        let n = self.n;
        let f = |x: f64| (x / self.a + n as f64 / 2.0 - 0.5).rem_euclid(n as f64);
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
}

/// Surfer velocity v = (R†σR − L†σL)/(|L|²+|R|²) (units of c). |v| ≤ 1 for any spinor.
pub fn surfer_velocity(p: &[C; 4]) -> [f64; 3] {
    let (l0, l1, r0, r1) = (p[0], p[1], p[2], p[3]);
    let rho = l0.norm2() + l1.norm2() + r0.norm2() + r1.norm2();
    if rho < 1e-300 { return [0.0; 3]; }
    let sx = |a: C, b: C| 2.0 * (a.conj() * b).re;
    let sy = |a: C, b: C| 2.0 * (a.conj() * b).im;
    let sz = |a: C, b: C| a.norm2() - b.norm2();
    [(sx(r0, r1) - sx(l0, l1)) / rho, (sy(r0, r1) - sy(l0, l1)) / rho, (sz(r0, r1) - sz(l0, l1)) / rho]
}

/// Light on a periodic cube: F = E + icB (3 components) per site, no flip.
pub struct Light3D {
    pub n: usize,
    pub f: Vec<[C; 3]>,
    tmp: Vec<[C; 2]>,
    pub ticks: u64,
}

impl Light3D {
    pub fn new(n: usize) -> Self { Self { n, f: vec![[C::ZERO; 3]; n * n * n], tmp: vec![[C::ZERO; 2]; n * n * n], ticks: 0 } }

    fn shift(&mut self, j: usize) {
        let n = self.n;
        let [u, d] = one_basis(j);
        self.tmp.par_iter_mut().zip(self.f.par_iter()).for_each(|(t, p)| {
            t[0] = u[0].conj() * p[0] + u[1].conj() * p[1] + u[2].conj() * p[2];
            t[1] = d[0].conj() * p[0] + d[1].conj() * p[1] + d[2].conj() * p[2];
        });
        let tmp = &self.tmp;
        self.f.par_iter_mut().enumerate().for_each(|(s, p)| {
            // F' = F + u₊ (c₊ shifted − c₊) + u₋ (c₋ shifted − c₋): the 0 part stays
            let dp = tmp[nb(s, j, -1, n)][0] - tmp[s][0];
            let dm = tmp[nb(s, j, 1, n)][1] - tmp[s][1];
            for q in 0..3 { p[q] += u[q] * dp + d[q] * dm; }
        });
    }

    pub fn tick(&mut self) {
        for j in order(self.ticks) { self.shift(j); }
        self.ticks += 1;
    }

    pub fn norm2(&self) -> f64 { self.f.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum() }
}

/// 1-D walk (Feynman checkerboard): R moves +1, L moves −1, then the flip. Periodic.
pub struct Walk1D {
    pub theta: f64,
    pub r: Vec<C>,
    pub l: Vec<C>,
}

impl Walk1D {
    pub fn tick(&mut self) {
        let n = self.r.len();
        self.r.rotate_right(1);
        self.l.rotate_left(1);
        let (ct, st) = (self.theta.cos(), self.theta.sin());
        let mi = C::new(0.0, -st);
        for i in 0..n {
            let (r, l) = (self.r[i], self.l[i]);
            self.r[i] = r * ct + mi * l;
            self.l[i] = l * ct + mi * r;
        }
    }
}

/// Two-tick Bloch matrix of the electron walk at dimensionless momentum k a (4×4, rows/cols L₀ L₁ R₀ R₁).
pub fn bloch_dirac(ka: [f64; 3], theta: f64) -> [[C; 4]; 4] {
    let mut w = ident4();
    for tick in 0..2u64 {
        let mut t = ident4();
        for j in order(tick) {
            let [u, d] = HALF[j];
            // one axis: exp(−i h k_j σ_j) = P diag(e^{−i h k}, e^{+i h k}) P†, h = −1 (L), +1 (R)
            let mut m = [[C::ZERO; 4]; 4];
            for (blk, h) in [(0usize, -1.0), (2usize, 1.0)] {
                let (ep, em) = (C::expi(-h * ka[j]), C::expi(h * ka[j]));
                for a in 0..2 {
                    for b in 0..2 {
                        m[blk + a][blk + b] = u[a] * ep * u[b].conj() + d[a] * em * d[b].conj();
                    }
                }
            }
            t = mul4(&m, &t);
        }
        let (ct, st) = (theta.cos(), theta.sin());
        let mut cm = [[C::ZERO; 4]; 4];
        for q in 0..2 {
            cm[q][q] = C::real(ct);
            cm[2 + q][2 + q] = C::real(ct);
            cm[q][2 + q] = C::new(0.0, -st);
            cm[2 + q][q] = C::new(0.0, -st);
        }
        w = mul4(&mul4(&cm, &t), &w);
    }
    w
}

/// Positive energy (units m c²) of the electron walk at momentum p (units m c) along x, from the
/// two-tick Bloch matrix: along an axis both spin states are degenerate, so Re tr W₂ = 4 cos 2ω.
pub fn walk_energy(p: f64, theta: f64) -> f64 {
    let w = bloch_dirac([p * theta, 0.0, 0.0], theta);
    let tr = (0..4).fold(C::ZERO, |s, i| s + w[i][i]);
    (tr.re / 4.0).clamp(-1.0, 1.0).acos() / 2.0 / theta
}

fn ident4() -> [[C; 4]; 4] {
    let mut m = [[C::ZERO; 4]; 4];
    for (i, row) in m.iter_mut().enumerate() { row[i] = C::ONE; }
    m
}

pub fn mul4(a: &[[C; 4]; 4], b: &[[C; 4]; 4]) -> [[C; 4]; 4] {
    let mut m = [[C::ZERO; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut s = C::ZERO;
            for k in 0..4 { s += a[i][k] * b[k][j]; }
            m[i][j] = s;
        }
    }
    m
}

// ─────────────────────────────── radial walk (hydrogen) ───────────────────────────────

/// One Dirac partial wave κ in the Coulomb field of a point nucleus (charge Z):
/// sites r_i = (i + ½) a; per site (χ₊ moving outwards, χ₋ moving inwards); at r = 0 the inward walker
/// turns into the outward one with a sign (g(0) = 0), at r_max the outward one turns back.
/// Rotation per tick: exp(−i a (m σ_z + κ/r σ_x − Zα/r)) acting on (g, f).
pub struct Radial {
    pub a: f64,
    pub n: usize,
    /// rotation in the (χ₊, χ₋) basis per site
    pub m: Vec<[[C; 2]; 2]>,
    /// half rotation (for the state at the symmetric point of the tick)
    pub mh: Vec<[[C; 2]; 2]>,
}

/// (g, f) = U χ with U = [u₊ u₋], u± = (1, ±i)/√2 (eigenvectors of σ_y).
const UCHI: [[C; 2]; 2] = [[C::real(S2), C::real(S2)], [C::new(0.0, S2), C::new(0.0, -S2)]];

fn rot_chi(phi: f64, nx: f64, nz: f64, v: f64) -> [[C; 2]; 2] {
    // e^{−iv} (cos φ − i sin φ (nx σx + nz σz)) in (g,f), then U† M U
    let (s, c) = phi.sin_cos();
    let e = C::expi(-v);
    let m = [[C::new(c, -s * nz), C::new(0.0, -s * nx)], [C::new(0.0, -s * nx), C::new(c, s * nz)]];
    let mut out = [[C::ZERO; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            let mut acc = C::ZERO;
            for k in 0..2 {
                for l in 0..2 { acc += UCHI[k][i].conj() * m[k][l] * UCHI[l][j]; }
            }
            out[i][j] = e * acc;
        }
    }
    out
}

impl Radial {
    pub fn new(kappa: i32, a: f64, rmax: f64, z_alpha: f64) -> Self {
        let n = (rmax / a).ceil() as usize;
        let mut m = Vec::with_capacity(n);
        let mut mh = Vec::with_capacity(n);
        for i in 0..n {
            let r = (i as f64 + 0.5) * a;
            let (hx, hz) = (kappa as f64 / r, 1.0);
            let nrm = (hx * hx + hz * hz).sqrt();
            let v = -z_alpha / r;
            m.push(rot_chi(a * nrm, hx / nrm, hz / nrm, v * a));
            mh.push(rot_chi(0.5 * a * nrm, hx / nrm, hz / nrm, 0.5 * v * a));
        }
        Self { a, n, m, mh }
    }

    /// One tick on a state stored as x[2i + s].
    pub fn apply(&self, x: &[C], y: &mut [C]) {
        let n = self.n;
        for i in 0..n {
            let p = if i >= 1 { x[2 * (i - 1)] } else { -x[1] };
            let q = if i + 1 < n { x[2 * (i + 1) + 1] } else { -x[2 * (n - 1)] };
            let m = &self.m[i];
            y[2 * i] = m[0][0] * p + m[0][1] * q;
            y[2 * i + 1] = m[1][0] * p + m[1][1] * q;
        }
    }

    /// Banded matrix W − λ I (kl = ku = 3) in the interleaved order 2i + s.
    fn band(&self, lam: C) -> Band {
        let n2 = 2 * self.n;
        let mut b = Band::new(n2, 3, 3);
        for i in 0..self.n {
            let m = &self.m[i];
            // source of P_i and Q_i
            let (sp, fp) = if i >= 1 { (2 * (i - 1), 1.0) } else { (1, -1.0) };
            let (sq, fq) = if i + 1 < self.n { (2 * (i + 1) + 1, 1.0) } else { (2 * (self.n - 1), -1.0) };
            for s in 0..2 {
                let row = 2 * i + s;
                b.add(row, sp, m[s][0] * fp);
                b.add(row, sq, m[s][1] * fq);
                b.add(row, row, -lam);
            }
        }
        b
    }

    /// Infinite listening at frequency ω (units m c²): inverse iteration on W at λ = e^{−iωa}.
    /// Returns (ω of the state, χ at the symmetric point of the tick, as (g, f) per site).
    pub fn eigen(&self, omega_guess: f64, iters: usize) -> (f64, Vec<[f64; 2]>) {
        let lam0 = C::expi(-omega_guess * self.a);
        let mut lu = self.band(lam0);
        lu.factor();
        let n2 = 2 * self.n;
        let mut x: Vec<C> = (0..n2).map(|i| C::real(1.0 + 0.001 * (i % 7) as f64)).collect();
        normalize(&mut x);
        let mut y = vec![C::ZERO; n2];
        let mut lam = lam0;
        for _ in 0..iters {
            lu.solve(&mut x);
            normalize(&mut x);
            self.apply(&x, &mut y);
            lam = dot(&x, &y);
        }
        let omega = -lam.arg() / self.a;
        // symmetric point: u = Mh† χ, then (g, f) = U u, global phase making g real
        let mut gf = vec![[C::ZERO; 2]; self.n];
        for i in 0..self.n {
            let (c0, c1) = (x[2 * i], x[2 * i + 1]);
            let h = &self.mh[i];
            let u0 = h[0][0].conj() * c0 + h[1][0].conj() * c1;
            let u1 = h[0][1].conj() * c0 + h[1][1].conj() * c1;
            gf[i] = [(u0 + u1) * S2, C::I * (u0 - u1) * S2];
        }
        let kmax = (0..self.n).max_by(|&i, &j| gf[i][0].norm2().total_cmp(&gf[j][0].norm2())).unwrap_or(0);
        let ph = C::expi(-gf[kmax][0].arg());
        let norm: f64 = gf.iter().map(|v| v[0].norm2() + v[1].norm2()).sum::<f64>() * self.a;
        let s = 1.0 / norm.sqrt();
        let out = gf.iter().map(|v| [(v[0] * ph).re * s, (v[1] * ph).re * s]).collect();
        (omega, out)
    }
}

/// ⟨x|y⟩ with compensated (Kahan) summation: the level is a tiny phase on top of the rest-energy phase,
/// so the round-off of a plain sum over ~10⁶ terms (~10⁻¹⁴) would show up at 10⁻⁸ eV.
pub fn dot(x: &[C], y: &[C]) -> C {
    let (mut s, mut c) = (C::ZERO, C::ZERO);
    for (a, b) in x.iter().zip(y) {
        let v = a.conj() * *b - c;
        let t = s + v;
        c = (t - s) - v;
        s = t;
    }
    s
}

fn normalize(x: &mut [C]) {
    let s = x.iter().map(|c| c.norm2()).sum::<f64>().sqrt();
    for c in x.iter_mut() { *c = c.scale(1.0 / s); }
}

/// Complex band matrix with partial pivoting LU (LAPACK gbtrf-style, row storage).
struct Band {
    n: usize,
    kl: usize,
    /// upper bandwidth after pivoting = kl + ku
    ku2: usize,
    w: usize,
    a: Vec<C>,
    piv: Vec<usize>,
}

impl Band {
    fn new(n: usize, kl: usize, ku: usize) -> Self {
        let ku2 = kl + ku;
        let w = kl + ku2 + 1;
        Self { n, kl, ku2, w, a: vec![C::ZERO; n * w], piv: vec![0; n] }
    }
    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize { i * self.w + (j + self.kl - i) }
    fn add(&mut self, i: usize, j: usize, v: C) { let k = self.idx(i, j); self.a[k] += v; }
    #[inline]
    fn get(&self, i: usize, j: usize) -> C {
        if j + self.kl < i || j > i + self.ku2 { C::ZERO } else { self.a[self.idx(i, j)] }
    }
    #[inline]
    fn set(&mut self, i: usize, j: usize, v: C) { let k = self.idx(i, j); self.a[k] = v; }

    fn factor(&mut self) {
        let n = self.n;
        for k in 0..n {
            let last = (k + self.kl).min(n - 1);
            let mut p = k;
            let mut best = self.get(k, k).norm2();
            for i in k + 1..=last {
                let v = self.get(i, k).norm2();
                if v > best { best = v; p = i; }
            }
            self.piv[k] = p;
            let jmax = (k + self.ku2).min(n - 1);
            if p != k {
                for j in k..=jmax {
                    let (x, y) = (self.get(k, j), self.get(p, j));
                    self.set(k, j, y);
                    self.set(p, j, x);
                }
            }
            let d = self.get(k, k);
            let d = if d.norm2() == 0.0 { C::real(1e-300) } else { d };
            for i in k + 1..=last {
                let f = self.get(i, k) / d;
                self.set(i, k, f);
                if f.norm2() == 0.0 { continue; }
                for j in k + 1..=jmax {
                    let v = self.get(i, j) - f * self.get(k, j);
                    self.set(i, j, v);
                }
            }
        }
    }

    fn solve(&self, b: &mut [C]) {
        let n = self.n;
        for k in 0..n {
            let p = self.piv[k];
            if p != k { b.swap(k, p); }
            let last = (k + self.kl).min(n - 1);
            for i in k + 1..=last {
                let f = self.get(i, k);
                let bk = b[k];
                b[i] -= f * bk;
            }
        }
        for k in (0..n).rev() {
            let jmax = (k + self.ku2).min(n - 1);
            let mut s = b[k];
            for j in k + 1..=jmax { s -= self.get(k, j) * b[j]; }
            b[k] = s / self.get(k, k);
        }
    }
}

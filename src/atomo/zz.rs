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

/// ONE SUBLATTICE ("colore", docs/SOTTORETICOLO.md, docs/NUCLEO_1S.md §12.1). In one tick every amplitude moves one
/// cell along x, one along y and one along z, so (i+j) mod 2 and (j+k) mod 2 never change: the cube is 4 lattices
/// that never exchange anything (not even through the nucleus, the field or the mass flip, which are on-site).
/// `Some([cx, cy])` keeps only the sites with (i+j)%2 = cx and (j+k)%2 = cy: N³/4 sites, a BCC lattice of side 2a.
/// Storage between ticks: index s = (i·H + j/2)·H + k/2 (H = N/2); the parities of j and k follow from i and the
/// label. Inside a tick (after the x or y shift) the amplitude sits on another colour: the label flips
/// (x: cx; y: cx and cy; z: cy) and after x, y, z (or z, y, x) it is back.
pub type Colore = Option<[u8; 2]>;
pub const H: usize = N / 2;
/// The colour of the Laboratorio (default after the exam of docs/SOTTORETICOLO.md): i ≡ j ≡ k (mod 2). It sees the
/// nucleus along the (1,1,1) diagonal: its nearest sites are (23,23,23) and (24,24,24).
pub const COLORE_LAB: [u8; 2] = [0, 0];
/// Whether the Laboratorio runs one colour or the whole cube by default. OFF: the exam of docs/SOTTORETICOLO.md
/// found that one colour alone is ANISOTROPIC (a z field also moves it along x and y by 5 % and 14 % of ⟨z⟩; the 4
/// colours cancel each other), so the coupled Laboratorio changes by 0.4 % in the dipole (E2 FAIL). The one-colour
/// rule itself is exact (E1). Overrides: RETICOLO (set by `--reticolo colore|pieno` or `--stati colore`) or the
/// environment variable ATOMO_RETICOLO=colore|pieno.
pub const SOTTORETICOLO: bool = false;
/// 0: default; 1: one colour; 2: whole cube.
pub static RETICOLO: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);
pub fn colore_lab() -> Colore {
    let env = std::env::var("ATOMO_RETICOLO").unwrap_or_default();
    let uno = match (RETICOLO.load(std::sync::atomic::Ordering::Relaxed), env.as_str()) {
        (1, _) => true,
        (2, _) => false,
        (_, "colore") => true,
        (_, "pieno") => false,
        _ => SOTTORETICOLO,
    };
    if uno { Some(COLORE_LAB) } else { None }
}
pub fn n_siti(col: Colore) -> usize { if col.is_some() { N * N * N / 4 } else { N * N * N } }
#[inline]
pub fn colore_di(c: [usize; 3]) -> [u8; 2] { [((c[0] + c[1]) & 1) as u8, ((c[1] + c[2]) & 1) as u8] }
#[inline]
fn decodifica(s: usize, lab: [u8; 2]) -> [usize; 3] {
    let i = s / (H * H);
    let j = 2 * ((s / H) % H) + ((i + lab[0] as usize) & 1);
    let k = 2 * (s % H) + ((j + lab[1] as usize) & 1);
    [i, j, k]
}
#[inline]
fn codifica(c: [usize; 3]) -> usize { (c[0] * H + c[1] / 2) * H + c[2] / 2 }
/// Cell indices (i, j, k) of storage site `s` between ticks.
#[inline]
pub fn sito(s: usize, col: Colore) -> [usize; 3] {
    match col { Some(l) => decodifica(s, l), None => [s / (N * N), (s / N) % N, s % N] }
}
/// Storage index of cell (i, j, k), if it belongs to the layout.
#[inline]
pub fn indice(c: [usize; 3], col: Colore) -> Option<usize> {
    match col {
        Some(l) => if colore_di(c) == l { Some(codifica(c)) } else { None },
        None => Some((c[0] * N + c[1]) * N + c[2]),
    }
}
#[inline]
fn gira(l: [u8; 2], j: usize) -> [u8; 2] { match j { 0 => [l[0] ^ 1, l[1]], 1 => [l[0] ^ 1, l[1] ^ 1], _ => [l[0], l[1] ^ 1] } }
/// A whole-cube wave restricted to one colour (compact storage), not renormalised.
#[allow(dead_code)]
pub fn estrai_colore(full: &[[C; 4]], col: [u8; 2]) -> Vec<[C; 4]> {
    (0..N * N * N / 4).map(|s| full[indice(decodifica(s, col), None).unwrap()]).collect()
}
/// A one-colour wave put back on the whole cube (zero on the other colours).
#[allow(dead_code)]
pub fn espandi_colore(v: &[[C; 4]], col: [u8; 2]) -> Vec<[C; 4]> {
    let mut out = vec![[C::ZERO; 4]; N * N * N];
    for (s, p) in v.iter().enumerate() { out[indice(decodifica(s, col), None).unwrap()] = *p; }
    out
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
    /// None: the whole cube; Some(colour): one sublattice only (docs/SOTTORETICOLO.md)
    pub colore: Colore,
    /// the colour the amplitude sits on right now (changes inside a tick)
    lab: [u8; 2],
}

/// The nucleus seen by the walk. `Nudo`: −Zα/r on every site. `Imparato` (docs/NUCLEO_1S.md §13–14, used by
/// docs/EMISSIONE.md §11): the potential of the 8 sites at r0 = √3/2·a times s0 and of the 24 sites at r1 = √11/2·a
/// times s1 (2 numbers learned for Z = 60 with an internal criterion: one 1s line, on the rule's blind pair U); at
/// other a they follow the line s0(a) = 0.5546 + 0.1011 (0.5 − a), s1(a) = 1.0035 − 0.0051 (0.5 − a) (phase 4).
/// It closes the 4 valleys of the cubic walk (1s split 0.0196 → 0.0002 mc²). Debt: 2 learned numbers, valid for Z = 60.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Nucleo { Nudo, Imparato }
impl Nucleo {
    pub fn nome(self) -> &'static str { match self { Nucleo::Nudo => "nudo", Nucleo::Imparato => "imparato (s0 0.554593, s1 1.003549)" } }
    /// (s0, s1) at the cell size A of this lattice
    pub fn fattori(self) -> (f64, f64) {
        match self {
            Nucleo::Nudo => (1.0, 1.0),
            // at a = 0.5 the phase-3 numbers themselves (the phase-4 line passes through them)
            Nucleo::Imparato => if (A - 0.5).abs() < 1e-12 { (0.554593, 1.003549) } else { (0.5546 + 0.1011 * (0.5 - A), 1.0035 - 0.0051 * (0.5 - A)) },
        }
    }
}
/// The nucleus of every new `Atomo` (Laboratorio, export).
pub const NUCLEO: Nucleo = Nucleo::Imparato;

/// Coulomb potential at a site, with the nucleus `nu`.
pub fn potenziale(x: [f64; 3], nu: Nucleo) -> f64 {
    let r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    let r = r2.sqrt();
    let (s0, s1) = nu.fattori();
    // r² in units of (a/2)²: 3 for the nearest 8 sites, 11 for the next 24
    let q = (r2 / (0.25 * A * A)).round() as i64;
    let f = if (r2 / (0.25 * A * A) - q as f64).abs() < 1e-6 { match q { 3 => s0, 11 => s1, _ => 1.0 } } else { 1.0 };
    -z_alpha() / r * f
}

pub fn x_of(i: usize) -> f64 { (i as f64 - N as f64 / 2.0 + 0.5) * A }

impl Atomo {
    /// The atom in its ground state: the 1s (m = +½) wave of the radial form of the same rule (same a,
    /// same Zα), put on the cubic lattice.
    #[allow(dead_code)]
    pub fn new() -> Self { Self::con_stato(-1, 1) }
    /// The same on one colour (`None`: the whole cube).
    pub fn new_in(col: Colore) -> Self { Self::con_stato_in(-1, 1, col) }

    /// The atom started in level (κ, n) of the radial rule (κ = −1: s½; κ = +1: p½), m = +½.
    #[allow(dead_code)]
    pub fn con_stato(kappa: i32, nn: u32) -> Self { Self::con_stato_in(kappa, nn, None) }

    /// The same, on one colour (the radial wave restricted to its sites and renormalised) or on the whole cube.
    pub fn con_stato_in(kappa: i32, nn: u32, col: Colore) -> Self {
        let ns = n_siti(col);
        let za = z_alpha();
        let lv = livello(kappa, nn);
        let mut psi = vec![[C::ZERO; 4]; ns];
        let mut vc = vec![0.0; ns];
        let mut zs = vec![0.0; ns];
        psi.par_iter_mut().zip(vc.par_iter_mut()).zip(zs.par_iter_mut()).enumerate().for_each(|(s, ((p, v), z))| {
            let c = sito(s, col);
            let x = [x_of(c[0]), x_of(c[1]), x_of(c[2])];
            let _ = za;
            *v = potenziale(x, NUCLEO);
            *z = x[2];
            let d = spinore(&lv, kappa, x);
            *p = dirac_a_chirale(d);
        });
        let norm: f64 = psi.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum();
        let k = 1.0 / norm.sqrt();
        psi.par_iter_mut().for_each(|p| for c in p.iter_mut() { *c = c.scale(k); });
        let psi0 = psi.clone();
        let pc = vc.iter().map(|v| C::expi(-v * A)).collect();
        Atomo { tmp: vec![[C::ZERO; 4]; ns], psi, pc, zs, ticks: 0, impulsi: vec![], campo_est: [0.0; 3], psi0, colore: col, lab: col.unwrap_or([0, 0]) }
    }

    /// Change the nucleus (the Coulomb phases of every site); the wave is untouched.
    pub fn con_nucleo(&mut self, nu: Nucleo) {
        let col = self.colore;
        self.pc = (0..self.psi.len()).map(|s| {
            let c = sito(s, col);
            C::expi(-potenziale([x_of(c[0]), x_of(c[1]), x_of(c[2])], nu) * A)
        }).collect();
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
        let mossa = |p: &mut [C; 4], lp: C, lm: C, rp: C, rm: C| {
            p[0] = u[0] * lp + d[0] * lm;
            p[1] = u[1] * lp + d[1] * lm;
            p[2] = u[0] * rp + d[0] * rm;
            p[3] = u[1] * rp + d[1] * rm;
        };
        match self.colore {
            None => self.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
                // forward (verso = +1): L u₊ moves −1, L u₋ +1, R u₊ +1, R u₋ −1; backward: the opposite
                let lp = tmp[nb(s, j, verso, n)][0];
                let lm = tmp[nb(s, j, -verso, n)][1];
                let rp = tmp[nb(s, j, -verso, n)][2];
                let rm = tmp[nb(s, j, verso, n)][3];
                mossa(p, lp, lm, rp, rm);
            }),
            Some(_) => {
                // one colour: the target sites carry the flipped label, the sources (one cell away along j) the old one
                let nuova = gira(self.lab, j);
                self.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
                    let c = decodifica(s, nuova);
                    let src = |dd: isize| { let mut q = c; q[j] = wrap(q[j], dd, n); codifica(q) };
                    let (a, b) = (src(verso), src(-verso));
                    mossa(p, tmp[a][0], tmp[b][1], tmp[b][2], tmp[a][3]);
                });
                self.lab = nuova;
            }
        }
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
        let col = self.colore;
        debug_assert!(col.map(|c| c == self.lab).unwrap_or(true));
        self.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
            let c0 = if verso > 0.0 { pc[s] } else { pc[s].conj() };
            let q = sito(s, col);
            let ph = if trasverso { c0 * pz[q[2]] * px[q[0]] * py[q[1]] } else { c0 * pz[q[2]] };
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
        let mut wtot = 0.0;
        let mut usati = 0;
        let mut media = [C::ZERO; 4];
        for (di, wx) in [(0usize, 1.0 - tx), (1, tx)] {
            for (dj, wy) in [(0usize, 1.0 - ty), (1, ty)] {
                for (dk, wz) in [(0usize, 1.0 - tz), (1, tz)] {
                    // one colour: only 2 of the 8 corners (opposite ones) are sites; their weights are renormalised
                    let Some(s) = indice([(i0 + di) % n, (j0 + dj) % n, (k0 + dk) % n], self.colore) else { continue };
                    let w = wx * wy * wz;
                    wtot += w;
                    usati += 1;
                    for q in 0..4 { out[q] += self.psi[s][q] * w; media[q] += self.psi[s][q]; }
                }
            }
        }
        if self.colore.is_some() {
            if wtot > 1e-9 { for c in out.iter_mut() { *c = c.scale(1.0 / wtot); } } else { for (o, m) in out.iter_mut().zip(media) { *o = m.scale(1.0 / usati as f64); } }
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
        let col = self.colore;
        self.psi.par_iter_mut().enumerate().for_each(|(s, q)| {
            let i = sito(s, col);
            let x = [x_of(i[0]), x_of(i[1]), x_of(i[2])];
            let d2: f64 = (0..3).map(|c| wrap(x[c] - p[c]).powi(2)).sum();
            let g = (-d2 / (2.0 * sig * sig)).exp();
            for c in 0..4 { q[c] = chi[c].scale(g); }
        });
        let norm: f64 = self.psi.par_iter().map(|q| q.iter().map(|c| c.norm2()).sum::<f64>()).sum();
        let k = 1.0 / norm.max(1e-300).sqrt();
        self.psi.par_iter_mut().for_each(|q| for c in q.iter_mut() { *c = c.scale(k); });
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
        assert_eq!(psi.len(), self.psi.len(), "wave and atom on different lattices (whole cube / one colour)");
        self.psi = psi;
        self.psi0 = self.psi.clone();
        self.ticks = 0;
    }

    /// ⟨r⟩ (λ_C): the position of the electron, all three components.
    pub fn posizione(&self) -> [f64; 3] {
        let col = self.colore;
        self.psi.par_iter().enumerate().map(|(s, p)| {
            let w: f64 = p.iter().map(|c| c.norm2()).sum();
            let c = sito(s, col);
            [w * x_of(c[0]), w * x_of(c[1]), w * x_of(c[2])]
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
        let col = self.colore;
        (0..k).map(|_| {
            let u = rng.uni() * s;
            let i = cdf.partition_point(|&c| c < u).min(w.len() - 1);
            match col {
                None => [x_of(i / (n * n)) + A * (rng.uni() - 0.5), x_of((i / n) % n) + A * (rng.uni() - 0.5), x_of(i % n) + A * (rng.uni() - 0.5)],
                Some(_) => {
                    // one colour: a site stands for 4 cells, its Wigner–Seitz cell of the BCC lattice
                    // (|δx|, |δy|, |δz| ≤ a and |δx| + |δy| + |δz| ≤ 3a/2), sampled uniformly by rejection
                    let c = sito(i, col);
                    loop {
                        let d = [A * (2.0 * rng.uni() - 1.0), A * (2.0 * rng.uni() - 1.0), A * (2.0 * rng.uni() - 1.0)];
                        if d[0].abs() + d[1].abs() + d[2].abs() <= 1.5 * A { break [x_of(c[0]) + d[0], x_of(c[1]) + d[1], x_of(c[2]) + d[2]]; }
                    }
                }
            }
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

/// Jacobi eigen-decomposition of a real symmetric matrix (row-major n×n): (eigenvalues, eigenvectors as columns).
fn jacobi(mut a: Vec<f64>, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut v = vec![0.0; n * n];
    for i in 0..n { v[i * n + i] = 1.0; }
    for _sweep in 0..100 {
        let off: f64 = (0..n).flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j))).map(|(i, j)| a[i * n + j].powi(2)).sum();
        if off < 1e-30 { break; }
        for p in 0..n { for q in p + 1..n {
            let apq = a[p * n + q];
            if apq.abs() < 1e-300 { continue; }
            let th = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
            let t = th.signum() / (th.abs() + (th * th + 1.0).sqrt());
            let t = if th == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for k in 0..n {
                let (akp, akq) = (a[k * n + p], a[k * n + q]);
                a[k * n + p] = c * akp - s * akq;
                a[k * n + q] = s * akp + c * akq;
            }
            for k in 0..n {
                let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                a[p * n + k] = c * apk - s * aqk;
                a[q * n + k] = s * apk + c * aqk;
            }
            for k in 0..n {
                let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                v[k * n + p] = c * vkp - s * vkq;
                v[k * n + q] = s * vkp + c * vkq;
            }
        }}
    }
    ((0..n).map(|i| a[i * n + i]).collect(), v)
}

/// One period of the rule (two ticks: the axis order alternates), from an even tick.
fn u2(psi: &[[C; 4]], base: &Atomo) -> Vec<[C; 4]> {
    let mut a = base.clone();
    a.impulsi.clear();
    a.campo_est = [0.0; 3];
    a.psi = psi.to_vec();
    a.ticks = 0;
    a.avanti();
    a.avanti();
    a.psi
}

fn norma(v: &[[C; 4]]) -> f64 { sovrapp(v, v).re.sqrt() }

fn sovrapp(a: &[[C; 4]], b: &[[C; 4]]) -> C {
    a.par_iter().zip(b.par_iter()).map(|(x, y)| { let mut s = C::ZERO; for q in 0..4 { s += x[q].conj() * y[q]; } s })
        .reduce(|| C::ZERO, |u, v| u + v)
}

/// A level found by `ritz`: frequency, residual ‖U₂v − e^{−iω}v‖ (0 for an exact level of the rule), the wave.
pub struct Livello3D { pub omega: f64, pub resid: f64, pub v: Vec<[C; 4]> }

/// FILTER DIAGONALISATION of the 3-D rule (the "listening" done properly): the start state is evolved without
/// light for `ticks` ticks; on EVEN ticks only (the rule repeats every two ticks, its axis order alternates)
/// Σ h(t) e^{iωt} ψ(t) is kept at frequencies lo, lo+dw, …, hi (Hann window h). These filtered waves span the
/// levels near [lo, hi]; the two-tick rule U₂ restricted to them is diagonalised (Rayleigh–Ritz). Each Ritz
/// vector comes with its residual: the levels are separated by their exact frequencies, not by the window width.
fn ritz(start: &Atomo, lo: f64, hi: f64, dw: f64, ticks: usize) -> Vec<Livello3D> {
    let mut a = start.clone();
    a.impulsi.clear();
    a.campo_est = [0.0; 3];
    a.ticks = 0;
    let nc = ((hi - lo) / dw).round() as usize + 1;
    let cand: Vec<f64> = (0..nc).map(|k| lo + dw * k as f64).collect();
    let len = a.psi.len();
    let mut acc: Vec<Vec<[C; 4]>> = vec![vec![[C::ZERO; 4]; len]; nc];
    let ticks = ticks + ticks % 2;
    for t in 0..ticks {
        if t % 2 == 0 {
            let w = (std::f64::consts::PI * (t as f64 + 1.0) / (ticks as f64 + 2.0)).sin().powi(2);
            let time = a.t();
            let fasi: Vec<C> = cand.iter().map(|om| C::expi(om * time).scale(w)).collect();
            let p = &a.psi;
            acc.par_iter_mut().zip(fasi.par_iter()).for_each(|(v, ph)| for (s, x) in v.iter_mut().zip(p.iter()) { for q in 0..4 { s[q] += x[q] * *ph; } });
        }
        a.avanti();
    }
    // orthonormal basis (modified Gram–Schmidt, dropping dependent directions)
    let mut qs: Vec<Vec<[C; 4]>> = vec![];
    let n0 = acc.iter().map(|v| norma(v)).fold(0.0, f64::max);
    for mut v in acc {
        for _ in 0..2 { for q in &qs {
            let c = sovrapp(q, &v);
            v.par_iter_mut().zip(q.par_iter()).for_each(|(x, y)| for k in 0..4 { x[k] -= c * y[k]; });
        }}
        let nv = norma(&v);
        if nv < 1e-7 * n0 { continue; }
        v.par_iter_mut().for_each(|x| for c in x.iter_mut() { *c = c.scale(1.0 / nv); });
        qs.push(v);
    }
    let m = qs.len();
    let uq: Vec<Vec<[C; 4]>> = qs.iter().map(|q| u2(q, start)).collect();
    let mut mm = vec![C::ZERO; m * m];
    for i in 0..m { for j in 0..m { mm[i * m + j] = sovrapp(&qs[i], &uq[j]); } }
    // M is (nearly) unitary hence normal: the joint eigenvectors of its Hermitian parts (M+M†)/2 and (M−M†)/2i
    // are its eigenvectors; diagonalise H = Re-part + γ Im-part (γ breaks the cos/sin degeneracies) as a real
    // symmetric 2m×2m matrix [[Re, −Im], [Im, Re]].
    let g = 0.37;
    let h = |i: usize, j: usize| -> C {
        let (a, b) = (mm[i * m + j], mm[j * m + i].conj());
        let herm = (a + b).scale(0.5);
        let anti = (a - b) * C::new(0.0, -0.5);
        herm + anti.scale(g)
    };
    let n2 = 2 * m;
    let mut r = vec![0.0; n2 * n2];
    for i in 0..m { for j in 0..m {
        let z = h(i, j);
        r[i * n2 + j] = z.re; r[(i + m) * n2 + j + m] = z.re;
        r[i * n2 + j + m] = -z.im; r[(i + m) * n2 + j] = z.im;
    }}
    let (_, vecs) = jacobi(r, n2);
    let mut out: Vec<Livello3D> = vec![];
    for k in 0..n2 {
        let c: Vec<C> = (0..m).map(|i| C::new(vecs[i * n2 + k], vecs[(i + m) * n2 + k])).collect();
        let nn: f64 = c.iter().map(|z| z.norm2()).sum::<f64>().sqrt();
        if nn < 1e-9 { continue; }
        let mut v = vec![[C::ZERO; 4]; len];
        for (ci, q) in c.iter().zip(qs.iter()) {
            let ci = ci.scale(1.0 / nn);
            v.par_iter_mut().zip(q.par_iter()).for_each(|(x, y)| for kk in 0..4 { x[kk] += ci * y[kk]; });
        }
        let nv = norma(&v);
        v.par_iter_mut().for_each(|x| for c in x.iter_mut() { *c = c.scale(1.0 / nv); });
        // the real doubling gives each eigenvector twice (c and i·c): keep one
        if out.iter().any(|o| sovrapp(&o.v, &v).norm2() > 0.5) { continue; }
        let uv = u2(&v, start);
        let lam = sovrapp(&v, &uv);
        let omega = -lam.arg();
        let res: f64 = uv.par_iter().zip(v.par_iter()).map(|(x, y)| (0..4).map(|q| (x[q] - lam * y[q]).norm2()).sum::<f64>()).sum::<f64>().sqrt();
        out.push(Livello3D { omega, resid: res, v });
    }
    out.sort_by(|a, b| a.resid.total_cmp(&b.resid));
    out
}


/// Weight of a wave inside r < `r` (λ_C) from the nucleus.
pub fn peso_interno(v: &[[C; 4]], r: f64, col: Colore) -> f64 {
    let (mut w1, mut tot) = (0.0, 0.0);
    for (s, p) in v.iter().enumerate() {
        let c = sito(s, col);
        let x = [x_of(c[0]), x_of(c[1]), x_of(c[2])];
        let w: f64 = p.iter().map(|c| c.norm2()).sum();
        tot += w;
        if x[0] * x[0] + x[1] * x[1] + x[2] * x[2] < r * r { w1 += w; }
    }
    w1 / tot
}

impl Atomo {
    /// The exact levels of the 3-D rule in [lo, hi] reached from this state (filter diagonalisation, see `ritz`).
    pub fn livelli_esatti(&self, lo: f64, hi: f64, dw: f64, ticks: usize) -> Vec<Livello3D> { ritz(self, lo, hi, dw, ticks) }

    /// THE 1s OF THE 3-D RULE, clean (docs/BATTIMENTO_1S.md, exam B1–B4 PASS). On this cubic lattice the 1s is a
    /// multiplet: the two-tick walk has 4 equivalent valleys (k = 0 and (π/2, π/2, 0) with its permutations) and the
    /// point nucleus mixes them, so the smooth radial 1s is spread over 4 exact levels (0.8404, 0.8863, 0.8993,
    /// 0.8996). From this state the exact levels in [0.82, 0.94] are found (residual < 1e-6, weight > 1 % here) and
    /// the one least touched by the first cells (least weight inside r < 1 λ_C) is returned, with the others.
    pub fn uno_s_esatto(&self, ticks: usize) -> (Livello3D, Vec<Livello3D>) {
        let lv = self.livelli_esatti(0.82, 0.94, 0.005, ticks);
        let mut sel: Vec<Livello3D> = lv.into_iter().filter(|l| l.resid < 1e-6 && sovrapp(&l.v, &self.psi).norm2() > 0.01).collect();
        sel.sort_by(|a, b| a.omega.total_cmp(&b.omega));
        let k = (0..sel.len()).min_by(|&i, &j| peso_interno(&sel[i].v, 1.0, self.colore).total_cmp(&peso_interno(&sel[j].v, 1.0, self.colore))).expect("no exact level near the 1s");
        let s = sel.remove(k);
        (s, sel)
    }
}

/// ⟨a| r |b⟩ on the lattice (λ_C): the dipole matrix element, three complex components.
pub fn dipolo_tra(a: &[[C; 4]], b: &[[C; 4]], col: Colore) -> [C; 3] {
    a.par_iter().zip(b.par_iter()).enumerate().map(|(s, (x, y))| {
        let mut o = C::ZERO;
        for q in 0..4 { o += x[q].conj() * y[q]; }
        let c = sito(s, col);
        [o.scale(x_of(c[0])), o.scale(x_of(c[1])), o.scale(x_of(c[2]))]
    }).reduce(|| [C::ZERO; 3], |u, v| [u[0] + v[0], u[1] + v[1], u[2] + v[2]])
}

/// An exact level reached from a lower one by the dipole: the level and ⟨lower| r |level⟩.
pub struct LivelloDipolo { pub lv: Livello3D, pub d: [C; 3] }
impl LivelloDipolo {
    pub fn d2(&self) -> f64 { self.d.iter().map(|c| c.norm2()).sum() }
}

impl Atomo {
    /// THE 2p LEVELS THE 1s TALKS TO, exact (docs/EMISSIONE.md §7, docs/ATOMO.md §3.7): the start state is
    /// r_axis · φ (the 1s pushed by a field along `asse`, which contains exactly the levels the light along that axis
    /// can reach), then the same filter diagonalisation as `uno_s_esatto` in [lo, hi]. Levels with residual < `tol` and
    /// weight > 1e-3 in the start are returned with their dipole ⟨φ| r |level⟩, the strongest |d_asse|² first.
    pub fn livelli_di_dipolo(&self, phi: &[[C; 4]], asse: usize, lo: f64, hi: f64, dw: f64, ticks: usize, tol: f64) -> Vec<LivelloDipolo> {
        let col = self.colore;
        let mut start = self.clone();
        let mut v: Vec<[C; 4]> = phi.par_iter().enumerate().map(|(s, p)| {
            let x = x_of(sito(s, col)[asse]);
            [p[0].scale(x), p[1].scale(x), p[2].scale(x), p[3].scale(x)]
        }).collect();
        let nv = norma(&v);
        v.par_iter_mut().for_each(|p| for c in p.iter_mut() { *c = c.scale(1.0 / nv); });
        start.imposta(v.clone());
        let lv = start.livelli_esatti(lo, hi, dw, ticks);
        let mut out: Vec<LivelloDipolo> = lv.into_iter().filter(|l| l.resid < tol && sovrapp(&l.v, &v).norm2() > 1e-3)
            .map(|l| { let d = dipolo_tra(phi, &l.v, col); LivelloDipolo { lv: l, d } }).collect();
        out.sort_by(|a, b| b.d[asse].norm2().total_cmp(&a.d[asse].norm2()));
        out
    }
}

impl Atomo {
    /// THE SMOOTH 1s WITH THE LEARNED NUCLEUS (docs/NUCLEO_1S.md §13–14): the 4 valley levels are pushed within ~2·10⁻⁴
    /// of each other, too close for a filter of `ticks` to tell apart; each single Ritz vector is still a valley mixture
    /// (checkerboard). The start state (the smooth radial 1s) is projected on ALL the exact levels within 0.0015 of the
    /// strongest one in [0.82, 0.94] (the criterion of docs/EMISSIONE.md §11.1): smooth, but it beats slowly with the
    /// residual split. Returns (state, weighted frequency, weight of the start state in the group).
    pub fn uno_s_liscio(&self, ticks: usize) -> (Vec<[C; 4]>, f64, f64) {
        let lv = self.livelli_esatti(0.82, 0.94, 0.005, ticks);
        let pesi: Vec<(f64, C)> = lv.iter().map(|l| (l.omega, sovrapp(&l.v, &self.psi))).collect();
        let k = (0..lv.len()).max_by(|&i, &j| pesi[i].1.norm2().total_cmp(&pesi[j].1.norm2())).expect("no level near the 1s");
        let w0 = lv[k].omega;
        let mut v = vec![[C::ZERO; 4]; self.psi.len()];
        let (mut wsum, mut om) = (0.0, 0.0);
        for (l, (w, c)) in lv.iter().zip(pesi.iter()) {
            if (w - w0).abs() > 0.0015 { continue; }
            wsum += c.norm2();
            om += c.norm2() * w;
            v.par_iter_mut().zip(l.v.par_iter()).for_each(|(x, y)| for q in 0..4 { x[q] += *c * y[q]; });
        }
        let nv = norma(&v);
        v.par_iter_mut().for_each(|p| for c in p.iter_mut() { *c = c.scale(1.0 / nv); });
        (v, om / wsum, wsum)
    }
}

impl Atomo {
    /// One tick forward that also returns the EXACT probability flux of the walk across every link (x, x + e_j) during
    /// the tick (whole cube only; docs/SURFISTI_ESATTI.md: the amplitude moving +1 at x minus the one moving −1 at x+e_j).
    pub fn avanti_flusso(&mut self) -> [Vec<f64>; 3] {
        assert!(self.colore.is_none(), "flux only on the whole cube");
        let n = N;
        let t = self.t();
        let mut f = [vec![0.0; self.psi.len()], vec![0.0; self.psi.len()], vec![0.0; self.psi.len()]];
        for j in order(self.ticks) {
            self.shift(j, 1);
            let tmp = &self.tmp;
            f[j].par_iter_mut().enumerate().for_each(|(s, v)| {
                let a = tmp[s][1].norm2() + tmp[s][2].norm2();
                let b = { let q = &tmp[nb(s, j, 1, n)]; q[0].norm2() + q[3].norm2() };
                *v = a - b;
            });
        }
        self.ruota(t, 1.0);
        self.ticks += 1;
        f
    }
}

/// Azimuthal flow of the 1s as in docs/SURFISTI_ESATTI.md (`vphi_from_flux`): link fluxes over two ticks centred on the
/// cells, J_φ = (−y F_x + x F_y)/r_c, mask 1 < r < 6, |z| < r/2: Σ J_φ r/r_c / Σ ρ / (Zα). Dirac's Bohm–Dirac flow of the
/// 1s gives 1 (v_φ = Zα c sinθ). The atom is advanced by two ticks.
pub fn vphi_rapporto(at: &mut Atomo) -> f64 {
    let n = N;
    let rho: Vec<f64> = at.psi.iter().map(|p| p.iter().map(|c| c.norm2()).sum()).collect();
    let f1 = at.avanti_flusso();
    let f2 = at.avanti_flusso();
    let (mut num, mut den) = (0.0, 0.0);
    for s in 0..rho.len() {
        let c = sito(s, None);
        let x = [x_of(c[0]), x_of(c[1]), x_of(c[2])];
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        let rc = (x[0] * x[0] + x[1] * x[1]).sqrt();
        if !(r > 1.0 && r < 6.0 && x[2].abs() < 0.5 * r) { continue; }
        let fc = |j: usize| 0.25 * (f1[j][s] + f2[j][s] + f1[j][nb(s, j, -1, n)] + f2[j][nb(s, j, -1, n)]);
        let jp = (-x[1] * fc(0) + x[0] * fc(1)) / rc;
        num += jp * r / rc;
        den += rho[s];
    }
    num / den / z_alpha()
}

/// Keep only the electron's own valley k = 0 (the Fourier cells nearer to k ∈ {0, π}³ than to the other three valleys
/// (π/2, π/2, 0)+{0, π}³ and permutations, as `valley_split` in lab/surf_f2_valli.py), renormalised. Whole cube.
pub fn proietta_valle_zero(psi: &[[C; 4]]) -> Vec<[C; 4]> {
    let n = N;
    let tw = std::f64::consts::TAU / n as f64;
    let dft = |a: &mut Vec<[C; 4]>, segno: f64| {
        let m: Vec<C> = (0..n * n).map(|q| C::expi(segno * tw * ((q / n) * (q % n)) as f64)).collect();
        for ax in 0..3 {
            let stride = [n * n, n, 1][ax];
            let starts: Vec<usize> = (0..n * n * n).filter(|s| (s / stride) % n == 0).collect();
            let lines: Vec<Vec<[C; 4]>> = starts.par_iter().map(|&s0| {
                (0..n).map(|k| { let mut o = [C::ZERO; 4]; for i in 0..n { let w = m[k * n + i]; let v = a[s0 + i * stride]; for q in 0..4 { o[q] += w * v[q]; } } o }).collect()
            }).collect();
            for (s0, l) in starts.iter().zip(lines) { for i in 0..n { a[s0 + i * stride] = l[i]; } }
        }
    };
    let mut a = psi.to_vec();
    dft(&mut a, -1.0);
    let h = std::f64::consts::FRAC_PI_2;
    let valli = [[0.0, 0.0, 0.0], [h, h, 0.0], [h, 0.0, h], [0.0, h, h]];
    let kv = |i: usize| { let k = tw * i as f64; if k > std::f64::consts::PI { k - std::f64::consts::TAU } else { k } };
    a.par_iter_mut().enumerate().for_each(|(s, p)| {
        let k = [kv(s / (n * n)), kv((s / n) % n), kv(s % n)];
        let dist = |q: [f64; 3]| {
            let mut best = f64::MAX;
            for m in 0..8 {
                let mut d = 0.0;
                for c in 0..3 {
                    let qq = q[c] + std::f64::consts::PI * ((m >> c) & 1) as f64;
                    let x = (k[c] - qq + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI;
                    d += x * x;
                }
                best = best.min(d);
            }
            best
        };
        let d0 = dist(valli[0]);
        if valli[1..].iter().any(|v| dist(*v) < d0) { *p = [C::ZERO; 4]; }
    });
    dft(&mut a, 1.0);
    let nv = norma(&a);
    a.par_iter_mut().for_each(|p| for c in p.iter_mut() { *c = c.scale(1.0 / nv); });
    a
}

/// Weight of a whole-cube wave in valley k = 0.
pub fn peso_valle_zero(psi: &[[C; 4]]) -> f64 {
    let p = proietta_valle_zero(psi);
    let nn: f64 = psi.iter().map(|q| q.iter().map(|c| c.norm2()).sum::<f64>()).sum();
    sovrapp(&p, psi).norm2() / nn
}

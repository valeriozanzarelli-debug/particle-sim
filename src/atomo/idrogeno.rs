//! IDROGENO VERO A SCALA DI BOHR (docs/IDROGENO_BOHR.md) — the hydrogen atom (Z = 1) live in 3-D with the rule
//! learned by the agglomerator (assets/atomo/agglomerato_zz.json) far from the proton and the zig-zag rule itself
//! (4 components, point Coulomb phase every tick) in a nested fine block around the proton.
//!
//! Units: ħ = c = μ = 1, μ = mₑ m_p/(mₑ + m_p) the reduced mass (lengths ħ/μc, times ħ/μc², energies μc²).
//! The Dirac–Coulomb problem with the reduced mass is the same problem in these units, so the learned levels
//! (cells in λ_C) are used as they are; eV = (energy) × μc² (the recoil enters only here: declared debt).
//!
//! Pieces (all derived from the rule or from the agglomerator's own descriptor convention; knobs are listed in
//! the doc):
//!   * `Livello`: a learned level (one Pauli component, spin part 1e-10 dropped and reported) on a cubic grid,
//!     one step = e^{−iVτ/2} · IFFT · e^{iφ₀} e^{−iG(k)} · FFT · e^{−iVτ/2}, V = cell average of −α/r (debt #3);
//!   * `Fine`: the zig-zag walk on a patch of B_f³ fine cells per coarse cell around the proton;
//!   * the glue (`prolunga`, `restringi`): fine cell values = conservative reconstruction of the coarse block
//!     means (the restriction of the prolonged field gives back the coarse values exactly), small component from
//!     the walk's own positive-energy eigenvectors; coarse value = block mean of the large component.
//!   * `Atomo`: coarse level + fine block, one coarse step = 32 fine ticks (Berger–Oliger with the halo as the
//!     light cone of one coarse step, as in docs/RETICOLO_UNIVERSO.md).
#![allow(dead_code)]

use rayon::prelude::*;

use crate::cplx::C;

pub const ALPHA: f64 = 7.297_352_5693e-3;
pub const MC2_EV: f64 = 510_998.950;
pub const MP_ME: f64 = 1836.152_673_43;
/// metres per λ_C (electron) and seconds per ħ/(mₑc²)
pub const LAMBDA_C_M: f64 = 3.861_592_6796e-13;
pub const T_UNIT_S: f64 = 1.288_088_668_19e-21;

/// μ/mₑ
pub fn mu() -> f64 { MP_ME / (1.0 + MP_ME) }
/// energy unit μc² in eV
pub fn ev() -> f64 { MC2_EV * mu() }
/// Bohr radius in units ħ/μc
pub fn a0() -> f64 { 1.0 / ALPHA }

/// Sommerfeld–Dirac energy minus the rest energy (units μc²) — CHECK ONLY, never an input of the rule.
pub fn dirac(n: u32, kappa: i32, za: f64) -> f64 {
    let k = kappa.abs() as f64;
    let d = k - (k * k - za * za).sqrt();
    1.0 / (1.0 + (za / (n as f64 - d)).powi(2)).sqrt() - 1.0
}
/// Dirac × reduced mass references in eV: 1s½, 2s½ (= 2p½), centroid of 2p ((2·2p½ + 4·2p3/2)/6)
pub fn riferimenti_ev() -> (f64, f64, f64) {
    let za = ALPHA;
    let e1 = dirac(1, -1, za) * ev();
    let e2s = dirac(2, -1, za) * ev();
    let e2p = (2.0 * dirac(2, 1, za) + 4.0 * dirac(2, -2, za)) / 6.0 * ev();
    (e1, e2s, e2p)
}

// =====================================================================================================
//  a minimal JSON reader (only what the asset files need)
// =====================================================================================================
#[derive(Clone, Debug)]
pub enum J { Null, B(bool), N(f64), S(String), A(Vec<J>), O(Vec<(String, J)>) }

impl J {
    pub fn get(&self, k: &str) -> &J {
        if let J::O(v) = self { for (kk, x) in v { if kk == k { return x; } } }
        panic!("json: key {k} missing")
    }
    pub fn at(&self, i: usize) -> &J { if let J::A(v) = self { &v[i] } else { panic!("json: not an array") } }
    pub fn num(&self) -> f64 { if let J::N(x) = self { *x } else { panic!("json: not a number") } }
    pub fn arr(&self) -> &Vec<J> { if let J::A(v) = self { v } else { panic!("json: not an array") } }
    pub fn nums(&self) -> Vec<f64> { self.arr().iter().map(|x| x.num()).collect() }
}

pub fn json_parse(s: &str) -> J {
    let b = s.as_bytes();
    let mut i = 0usize;
    let v = jv(b, &mut i);
    v
}
fn ws(b: &[u8], i: &mut usize) { while *i < b.len() && (b[*i] as char).is_whitespace() { *i += 1; } }
fn jv(b: &[u8], i: &mut usize) -> J {
    ws(b, i);
    match b[*i] {
        b'{' => {
            *i += 1;
            let mut v = vec![];
            loop {
                ws(b, i);
                if b[*i] == b'}' { *i += 1; break; }
                let k = if let J::S(s) = jv(b, i) { s } else { panic!("json key") };
                ws(b, i);
                *i += 1; // ':'
                let x = jv(b, i);
                v.push((k, x));
                ws(b, i);
                if b[*i] == b',' { *i += 1; }
            }
            J::O(v)
        }
        b'[' => {
            *i += 1;
            let mut v = vec![];
            loop {
                ws(b, i);
                if b[*i] == b']' { *i += 1; break; }
                v.push(jv(b, i));
                ws(b, i);
                if b[*i] == b',' { *i += 1; }
            }
            J::A(v)
        }
        b'"' => {
            *i += 1;
            let mut s = String::new();
            while b[*i] != b'"' {
                if b[*i] == b'\\' { *i += 1; }
                s.push(b[*i] as char);
                *i += 1;
            }
            *i += 1;
            J::S(s)
        }
        b't' => { *i += 4; J::B(true) }
        b'f' => { *i += 5; J::B(false) }
        b'n' => { *i += 4; J::Null }
        b'N' => { *i += 3; J::N(f64::NAN) }
        b'I' => { *i += 8; J::N(f64::INFINITY) }
        _ => {
            let st = *i;
            if b[*i] == b'-' && b.get(*i + 1) == Some(&b'I') { *i += 9; return J::N(f64::NEG_INFINITY); }
            while *i < b.len() && (b[*i] == b'-' || b[*i] == b'+' || b[*i] == b'.' || b[*i] == b'e' || b[*i] == b'E'
                || b[*i].is_ascii_digit()) { *i += 1; }
            J::N(std::str::from_utf8(&b[st..*i]).unwrap().parse().unwrap())
        }
    }
}

// =====================================================================================================
//  FFT: Stockham autosort, radices 4, 2, 3, batched over W lanes; 3-D transforms on n³ (x slowest)
// =====================================================================================================
#[derive(Clone)]
pub struct Fft { pub n: usize, radici: Vec<usize> }

const LANES: usize = 8;

#[derive(Clone, Copy)]
struct SendPtr(*mut C);
unsafe impl Send for SendPtr {}
unsafe impl Sync for SendPtr {}

impl Fft {
    pub fn new(n: usize) -> Self {
        let mut m = n;
        let mut r = vec![];
        while m % 4 == 0 { r.push(4); m /= 4; }
        while m % 2 == 0 { r.push(2); m /= 2; }
        while m % 3 == 0 { r.push(3); m /= 3; }
        assert!(m == 1, "fft: n = {n} must be 2^a 3^b");
        Fft { n, radici: r }
    }

    /// in-place transform of W interleaved lines: element (i, lane) at buf[i*W + lane]; sign −1 forward
    fn linee(&self, buf: &mut [C], tmp: &mut [C], w: usize, inversa: bool) {
        let n = self.n;
        let sgn = if inversa { 1.0 } else { -1.0 };
        let mut ns = 1usize;
        let mut src_is_buf = true;
        for &r in &self.radici {
            let (src, dst): (&[C], &mut [C]) = if src_is_buf { (&*buf, &mut *tmp) } else { (&*tmp, &mut *buf) };
            let m = n / r;
            for j in 0..m {
                let k = j % ns;
                let ang = sgn * 2.0 * std::f64::consts::PI * k as f64 / (ns * r) as f64;
                let base = (j / ns) * ns * r + k;
                let w1 = C::expi(ang);
                let w2 = w1 * w1;
                let w3 = w2 * w1;
                for l in 0..w {
                    let x0 = src[j * w + l];
                    let x1 = src[(j + m) * w + l] * w1;
                    if r == 2 {
                        dst[base * w + l] = x0 + x1;
                        dst[(base + ns) * w + l] = x0 - x1;
                    } else if r == 4 {
                        let x2 = src[(j + 2 * m) * w + l] * w2;
                        let x3 = src[(j + 3 * m) * w + l] * w3;
                        let a = x0 + x2;
                        let b = x0 - x2;
                        let c = x1 + x3;
                        let d = x1 - x3;
                        // ω = e^{sgn i π/2}: forward −i
                        let id = C::new(-sgn * d.im, sgn * d.re); // sgn·i·d
                        dst[base * w + l] = a + c;
                        dst[(base + ns) * w + l] = b + id;
                        dst[(base + 2 * ns) * w + l] = a - c;
                        dst[(base + 3 * ns) * w + l] = b - id;
                    } else {
                        let x2 = src[(j + 2 * m) * w + l] * w2;
                        let s = x1 + x2;
                        let d = x1 - x2;
                        let h = (3f64).sqrt() / 2.0;
                        let t = x0 - s.scale(0.5);
                        let id = C::new(-sgn * h * d.im, sgn * h * d.re);
                        dst[base * w + l] = x0 + s;
                        dst[(base + ns) * w + l] = t + id;
                        dst[(base + 2 * ns) * w + l] = t - id;
                    }
                }
            }
            ns *= r;
            src_is_buf = !src_is_buf;
        }
        if !src_is_buf { buf.copy_from_slice(&tmp[..buf.len()]); }
    }

    /// 3-D transform of a[(x n + y) n + z]; the inverse is normalized (1/n³)
    pub fn tre(&self, a: &mut [C], inversa: bool) {
        let n = self.n;
        assert_eq!(a.len(), n * n * n);
        let p = SendPtr(a.as_mut_ptr());
        let w = LANES.min(n);
        // axis z: for each x plane, chunks of w consecutive y-lines (gather/scatter)
        (0..n).into_par_iter().for_each(|x| {
            let p = p;
            let mut buf = vec![C::ZERO; n * w];
            let mut tmp = vec![C::ZERO; n * w];
            let mut y0 = 0;
            while y0 < n {
                let ww = w.min(n - y0);
                for l in 0..ww { for z in 0..n { unsafe { buf[z * ww + l] = *p.0.add((x * n + y0 + l) * n + z); } } }
                self.linee(&mut buf[..n * ww], &mut tmp[..n * ww], ww, inversa);
                for l in 0..ww { for z in 0..n { unsafe { *p.0.add((x * n + y0 + l) * n + z) = buf[z * ww + l]; } } }
                y0 += ww;
            }
        });
        // axis y: for each x plane, chunks of w consecutive z (contiguous)
        (0..n).into_par_iter().for_each(|x| {
            let p = p;
            let mut buf = vec![C::ZERO; n * w];
            let mut tmp = vec![C::ZERO; n * w];
            let mut z0 = 0;
            while z0 < n {
                let ww = w.min(n - z0);
                for y in 0..n { for l in 0..ww { unsafe { buf[y * ww + l] = *p.0.add((x * n + y) * n + z0 + l); } } }
                self.linee(&mut buf[..n * ww], &mut tmp[..n * ww], ww, inversa);
                for y in 0..n { for l in 0..ww { unsafe { *p.0.add((x * n + y) * n + z0 + l) = buf[y * ww + l]; } } }
                z0 += ww;
            }
        });
        // axis x: for each y, chunks of w consecutive z
        (0..n).into_par_iter().for_each(|y| {
            let p = p;
            let mut buf = vec![C::ZERO; n * w];
            let mut tmp = vec![C::ZERO; n * w];
            let mut z0 = 0;
            while z0 < n {
                let ww = w.min(n - z0);
                for x in 0..n { for l in 0..ww { unsafe { buf[x * ww + l] = *p.0.add((x * n + y) * n + z0 + l); } } }
                self.linee(&mut buf[..n * ww], &mut tmp[..n * ww], ww, inversa);
                for x in 0..n { for l in 0..ww { unsafe { *p.0.add((x * n + y) * n + z0 + l) = buf[x * ww + l]; } } }
                z0 += ww;
            }
        });
        if inversa {
            let k = 1.0 / (n * n * n) as f64;
            a.par_iter_mut().for_each(|c| *c = c.scale(k));
        }
    }
}

// =====================================================================================================
//  small vector helpers (parallel)
// =====================================================================================================
pub fn dot(a: &[C], b: &[C]) -> C {
    a.par_chunks(4096).zip(b.par_chunks(4096)).map(|(x, y)| {
        let mut s = C::ZERO;
        for i in 0..x.len() { s += x[i].conj() * y[i]; }
        s
    }).reduce(|| C::ZERO, |p, q| p + q)
}
pub fn norm2(a: &[C]) -> f64 { a.par_iter().map(|c| c.norm2()).sum() }
/// a ← a + k b
pub fn axpy(a: &mut [C], k: C, b: &[C]) { a.par_iter_mut().zip(b.par_iter()).for_each(|(x, y)| *x += k * *y); }
pub fn scal(a: &mut [C], k: C) { a.par_iter_mut().for_each(|x| *x = *x * k); }

// =====================================================================================================
//  a learned level of the agglomerator on a cubic periodic grid
// =====================================================================================================
/// The learned level as read from the asset (spin-averaged scalar stencil g_j = (G₀₀ + G₁₁)/2).
#[derive(Clone)]
pub struct Regola {
    pub b_fine: usize,
    pub h: f64,
    pub phi0: f64,
    pub offs: Vec<[i32; 3]>,
    pub g: Vec<C>,
    /// max |G₀₁| + max |G₀₀ − G₁₁| over the stencil (the spin part that is dropped)
    pub spin: f64,
    pub theta_fine: f64,
}

pub fn carica_regola(path: &str, catena: &str, b_fine: usize) -> Regola {
    let s = std::fs::read_to_string(path).expect("agglomerato_zz.json");
    let j = json_parse(&s);
    let ch = j.get(catena);
    let th = ch.get("theta_fine").num();
    for lv in ch.get("levels").arr() {
        if lv.get("B_fine_cells").num() as usize != b_fine { continue; }
        let offs: Vec<[i32; 3]> = lv.get("offsets").arr().iter().map(|o| {
            let v = o.nums();
            [v[0] as i32, v[1] as i32, v[2] as i32]
        }).collect();
        let re = lv.get("generator_G").get("re");
        let im = lv.get("generator_G").get("im");
        let comp = |a: usize, b: usize| -> Vec<C> {
            let r = re.at(a).at(b).nums();
            let i = im.at(a).at(b).nums();
            r.iter().zip(i.iter()).map(|(x, y)| C::new(*x, *y)).collect()
        };
        let (g00, g11, g01) = (comp(0, 0), comp(1, 1), comp(0, 1));
        let g: Vec<C> = g00.iter().zip(g11.iter()).map(|(a, b)| (*a + *b).scale(0.5)).collect();
        let s01 = g01.iter().map(|c| c.abs()).fold(0.0, f64::max);
        let sd = g00.iter().zip(g11.iter()).map(|(a, b)| (*a - *b).abs()).fold(0.0, f64::max);
        return Regola { b_fine, h: lv.get("cell_lambdaC").num(), phi0: lv.get("rest_phase_phi0").num(), offs, g,
                        spin: s01 + sd, theta_fine: th };
    }
    panic!("level B = {b_fine} not in {catena}")
}

/// Cell average of 1/r over the cube of side h centred at (x, y, z) (cells near the nucleus sub-sampled).
/// the field of the nucleus seen by the electron, V(r) (units μc²): point Coulomb −Zα/r, or (diagnostic only,
/// env IDRO_POT=gauss) a smooth Gaussian well of depth 2 Ha and width a₀/2, to test the glue without a cusp
pub fn pot(r: f64, z: f64) -> f64 {
    if std::env::var("IDRO_POT").map(|v| v == "gauss").unwrap_or(false) {
        let s = 0.5 * a0();
        -z * 8.0 * ALPHA * ALPHA * (-(r * r) / (2.0 * s * s)).exp()
    } else { -z * ALPHA / r }
}

pub fn media_pot(x: [f64; 3], h: f64, nsub: usize, z: f64) -> f64 {
    let mut s = 0.0;
    for a in 0..nsub { for b in 0..nsub { for c in 0..nsub {
        let u = [(a as f64 + 0.5) / nsub as f64 - 0.5, (b as f64 + 0.5) / nsub as f64 - 0.5, (c as f64 + 0.5) / nsub as f64 - 0.5];
        let p = [x[0] + u[0] * h, x[1] + u[1] * h, x[2] + u[2] * h];
        s += pot((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt(), z);
    }}}
    s / (nsub * nsub * nsub) as f64
}

pub fn media_inv_r(x: [f64; 3], h: f64, nsub: usize) -> f64 {
    let mut s = 0.0;
    for a in 0..nsub {
        for b in 0..nsub {
            for c in 0..nsub {
                let u = [(a as f64 + 0.5) / nsub as f64 - 0.5, (b as f64 + 0.5) / nsub as f64 - 0.5,
                         (c as f64 + 0.5) / nsub as f64 - 0.5];
                let p = [x[0] + u[0] * h, x[1] + u[1] * h, x[2] + u[2] * h];
                s += 1.0 / (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            }
        }
    }
    s / (nsub * nsub * nsub) as f64
}

/// One learned level on an n³ periodic box, the nucleus at the centre of cell c0 = n/2.
pub struct Livello {
    pub n: usize,
    pub h: f64,
    pub tau: f64,
    pub phi0: f64,
    pub c0: usize,
    pub fft: Fft,
    /// G(k) (real: the stencil is Hermitian) and the step multiplier e^{iφ₀} e^{−iG(k)}
    pub gk: Vec<f64>,
    pub molt: Vec<C>,
    /// potential V (μc²) per cell and the half kick e^{−iVτ/2}
    pub v: Vec<f64>,
    pub calcio: Vec<C>,
    pub z: f64,
    /// 0: nucleus at the centre of cell c0; 0.5: nucleus at the corner between cells c0−1 and c0
    pub off: f64,
    /// parity-odd part of the learned stencil: max |Im g_j| (learning noise; kept, reported)
    pub im_g: f64,
    pub spin: f64,
    /// the stencil (kept for the velocity operator ∂G/∂k)
    pub offs: Vec<[i32; 3]>,
    pub g: Vec<C>,
}

impl Livello {
    pub fn new(reg: &Regola, n: usize, z: f64) -> Self { Self::con_nucleo(reg, n, z, 0.0) }

    pub fn con_nucleo(reg: &Regola, n: usize, z: f64, off: f64) -> Self {
        let fft = Fft::new(n);
        let mut a = vec![C::ZERO; n * n * n];
        let w = |o: i32| -> usize { o.rem_euclid(n as i32) as usize };
        for (o, g) in reg.offs.iter().zip(reg.g.iter()) {
            a[(w(o[0]) * n + w(o[1])) * n + w(o[2])] += *g;
        }
        // G(k) = Σ_j g_j e^{+i k·j}: unnormalized inverse transform
        fft.tre(&mut a, true);
        let nn = (n * n * n) as f64;
        let gk: Vec<f64> = a.par_iter().map(|c| c.re * nn).collect();
        let molt: Vec<C> = gk.par_iter().map(|g| C::expi(reg.phi0 - g)).collect();
        let im_g = reg.g.iter().map(|c| c.im.abs()).fold(0.0, f64::max);
        let mut lv = Livello { n, h: reg.h, tau: reg.h, phi0: reg.phi0, c0: n / 2, fft, gk, molt, v: vec![0.0; n * n * n],
                               calcio: vec![C::ONE; n * n * n], z, off, im_g, spin: reg.spin,
                               offs: reg.offs.clone(), g: reg.g.clone() };
        lv.metti_potenziale(z);
        lv
    }

    pub fn x(&self, i: usize) -> f64 { (i as f64 - self.c0 as f64 + self.off) * self.h }
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize { (i * self.n + j) * self.n + k }

    /// V = −Zα ⟨1/r⟩_cell (the nucleus at the centre of cell c0); z = 0: free
    pub fn metti_potenziale(&mut self, z: f64) {
        let n = self.n;
        let h = self.h;
        let c0 = self.c0 as i64;
        let off = self.off;
        self.z = z;
        let tau = self.tau;
        self.v.par_iter_mut().zip(self.calcio.par_iter_mut()).enumerate().for_each(|(s, (v, k))| {
            let (i, j, l) = (s / (n * n), (s / n) % n, s % n);
            let d = [(i as i64 - c0), (j as i64 - c0), (l as i64 - c0)];
            let x = [(d[0] as f64 + off) * h, (d[1] as f64 + off) * h, (d[2] as f64 + off) * h];
            let near = d.iter().map(|q| (*q as f64 + off).abs().ceil() as i64).max().unwrap();
            let punto = std::env::var("IDRO_PUNTO").is_ok() && (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) > 0.0;
            let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            *v = if near <= 4 && !punto { media_pot(x, h, 24, z) } else { pot(r, z) };
            *k = C::expi(-0.5 * *v * tau);
        });
        if std::env::var("IDRO_CONTATTO").is_ok() {
            // diagnostic only: second-order covariance of the cell mean, Hermitian part −(h²/24)⟨∇²V⟩_cell (contact)
            let dv = -std::f64::consts::PI * z * ALPHA / (6.0 * h);
            let cells: Vec<[i64; 3]> = if off == 0.0 { vec![[0, 0, 0]] } else {
                let mut v = vec![]; for a in [-1i64, 0] { for b in [-1i64, 0] { for c in [-1i64, 0] { v.push([a, b, c]); } } } v };
            let w = 1.0 / cells.len() as f64;
            for d in cells {
                let s = (((c0 + d[0]) as usize * n) + (c0 + d[1]) as usize) * n + (c0 + d[2]) as usize;
                self.v[s] += dv * w;
                self.calcio[s] = C::expi(-0.5 * self.v[s] * tau);
            }
        }
    }

    /// one step of the learned rule with the Coulomb half kicks
    pub fn passo(&self, f: &mut [C]) {
        f.par_iter_mut().zip(self.calcio.par_iter()).for_each(|(a, k)| *a = *a * *k);
        self.fft.tre(f, false);
        f.par_iter_mut().zip(self.molt.par_iter()).for_each(|(a, m)| *a = *a * *m);
        self.fft.tre(f, true);
        f.par_iter_mut().zip(self.calcio.par_iter()).for_each(|(a, k)| *a = *a * *k);
    }

    /// the rule's own velocity operator v_a(k) = ∂G/∂k_a (units c; k in 1/cell, G in rad per step, h = τ)
    pub fn velocita_k(&self) -> [Vec<f64>; 3] {
        let n = self.n;
        let w = |o: i32| -> usize { o.rem_euclid(n as i32) as usize };
        let mk = |a: usize| -> Vec<f64> {
            let mut arr = vec![C::ZERO; n * n * n];
            for (o, g) in self.offs.iter().zip(self.g.iter()) {
                arr[(w(o[0]) * n + w(o[1])) * n + w(o[2])] += *g * C::new(0.0, o[a] as f64);
            }
            self.fft.tre(&mut arr, true);
            let nn = (n * n * n) as f64;
            arr.iter().map(|c| c.re * nn).collect()
        };
        [mk(0), mk(1), mk(2)]
    }

    /// probability current of the rule j_a = Re(f* v_a f) per cell (units c), with the velocity operator above
    pub fn corrente(&self, f: &[C], vk: &[Vec<f64>; 3]) -> [Vec<f64>; 3] {
        let mut fk = f.to_vec();
        self.fft.tre(&mut fk, false);
        let one = |a: usize| -> Vec<f64> {
            let mut t: Vec<C> = fk.par_iter().zip(vk[a].par_iter()).map(|(x, v)| x.scale(*v)).collect();
            self.fft.tre(&mut t, true);
            t.par_iter().zip(f.par_iter()).map(|(vt, x)| (x.conj() * *vt).re).collect()
        };
        [one(0), one(1), one(2)]
    }

    /// n steps of the FREE learned rule at once (exact: the multiplier to the n-th power) — used for time blocking
    pub fn molt_n(&self, m: f64) -> Vec<C> {
        self.gk.par_iter().map(|g| C::expi(m * (self.phi0 - g))).collect()
    }

    /// H f = IFFT(G/τ FFT f) + V f (Hermitian generator of the level with the potential)
    pub fn ham(&self, f: &[C], out: &mut [C]) {
        out.copy_from_slice(f);
        self.fft.tre(out, false);
        let tau = self.tau;
        out.par_iter_mut().zip(self.gk.par_iter()).for_each(|(a, g)| *a = a.scale(*g / tau));
        self.fft.tre(out, true);
        out.par_iter_mut().zip(f.par_iter()).zip(self.v.par_iter()).for_each(|((o, x), v)| *o += x.scale(*v));
    }

    /// energy read from the RULE: E = −arg(⟨f|W f⟩ e^{−iφ₀})/τ, and the residual |Wf − λ f|/|f|
    pub fn energia_regola(&self, f: &[C]) -> (f64, f64) {
        let mut w = f.to_vec();
        self.passo(&mut w);
        let nf = norm2(f);
        let lam = dot(f, &w).scale(1.0 / nf);
        let mut r = w;
        axpy(&mut r, C::ZERO - lam, f);
        let e = -(lam * C::expi(-self.phi0)).arg() / self.tau;
        (e, (norm2(&r) / nf).sqrt())
    }

    /// preconditioner 1/(G/τ + σ) in Fourier space
    pub fn precond(&self, f: &mut [C], sigma: f64) {
        self.fft.tre(f, false);
        let tau = self.tau;
        f.par_iter_mut().zip(self.gk.par_iter()).for_each(|(a, g)| *a = a.scale(1.0 / (g / tau + sigma)));
        self.fft.tre(f, true);
    }

    /// initial guesses (hydrogenic shapes, only to start the solver)
    pub fn guess(&self, kind: &str) -> Vec<C> {
        let n = self.n;
        let a0 = a0();
        let mut f: Vec<C> = (0..n * n * n).into_par_iter().map(|s| {
            let (i, j, l) = (s / (n * n), (s / n) % n, s % n);
            let (x, y, z) = (self.x(i), self.x(j), self.x(l));
            let r = (x * x + y * y + z * z).sqrt() / a0;
            let v = match kind { "1s" => (-r).exp(), "2s" => (1.0 - r / 2.0) * (-r / 2.0).exp(), "2p" => z / a0 * (-r / 2.0).exp(),
                                 _ => panic!() };
            C::real(v)
        }).collect();
        let k = 1.0 / norm2(&f).sqrt();
        scal(&mut f, C::real(k));
        f
    }
}

/// parity ⟨ψ|Pψ⟩/⟨ψ|ψ⟩ with P the inversion through the nucleus (off = 0.5: i → n−1−i; off = 0: i → n−i)
pub fn parita(lv: &Livello, f: &[C]) -> f64 {
    let n = lv.n;
    let m = |i: usize| if lv.off > 0.25 { n - 1 - i } else { (n - i) % n };
    let s: C = (0..n * n * n).into_par_iter().map(|s| {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        f[s].conj() * f[lv.idx(m(i), m(j), m(k))]
    }).reduce(|| C::ZERO, |a, b| a + b);
    s.re / norm2(f)
}

/// save / load a state as complex f32 little endian (scratch only, never committed)
pub fn salva_stato(path: &str, f: &[C]) {
    let mut b: Vec<u8> = Vec::with_capacity(f.len() * 8);
    for c in f { b.extend_from_slice(&(c.re as f32).to_le_bytes()); b.extend_from_slice(&(c.im as f32).to_le_bytes()); }
    std::fs::write(path, b).unwrap();
}
pub fn carica_stato(path: &str) -> Vec<C> {
    let b = std::fs::read(path).unwrap();
    b.chunks(8).map(|q| C::new(f32::from_le_bytes([q[0], q[1], q[2], q[3]]) as f64, f32::from_le_bytes([q[4], q[5], q[6], q[7]]) as f64)).collect()
}

// =====================================================================================================
//  small Hermitian eigenproblem (complex Jacobi)
// =====================================================================================================
pub fn herm_eig(a: &mut Vec<Vec<C>>) -> (Vec<f64>, Vec<Vec<C>>) {
    let n = a.len();
    let mut v = vec![vec![C::ZERO; n]; n];
    for i in 0..n { v[i][i] = C::ONE; }
    for _sweep in 0..60 {
        let mut off = 0.0;
        for p in 0..n { for q in p + 1..n { off += a[p][q].norm2(); } }
        if off < 1e-30 { break; }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q].abs() < 1e-300 { continue; }
                // phase: make a[p][q] real positive (column q × e^{−iφ}, row q × e^{iφ})
                let ph = C::expi(-a[p][q].arg());
                for k in 0..n { a[k][q] = a[k][q] * ph; }
                for k in 0..n { a[q][k] = a[q][k] * ph.conj(); }
                for k in 0..n { v[k][q] = v[k][q] * ph; }
                let apq = a[p][q].re;
                let th = (a[q][q].re - a[p][p].re) / (2.0 * apq);
                let t = th.signum() / (th.abs() + (th * th + 1.0).sqrt());
                let t = if th == 0.0 { 1.0 } else { t };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = akp.scale(c) - akq.scale(s);
                    a[k][q] = akp.scale(s) + akq.scale(c);
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = apk.scale(c) - aqk.scale(s);
                    a[q][k] = apk.scale(s) + aqk.scale(c);
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = vkp.scale(c) - vkq.scale(s);
                    v[k][q] = vkp.scale(s) + vkq.scale(c);
                }
            }
        }
    }
    let ev: Vec<f64> = (0..n).map(|i| a[i][i].re).collect();
    (ev, v)
}

// =====================================================================================================
//  lowest eigenvector of H in a deflated subspace: single-vector LOBPCG with the kinetic preconditioner
// =====================================================================================================
pub struct Autostato { pub e_gen: f64, pub e_regola: f64, pub resid_regola: f64, pub resid: f64, pub iter: usize, pub v: Vec<C> }

fn proietta_fuori(f: &mut [C], y: &[&Vec<C>]) {
    for q in y { let c = dot(q, f); axpy(f, C::ZERO - c, q); }
}

pub fn lobpcg(lv: &Livello, x0: Vec<C>, defl: &[&Vec<C>], sigma: f64, tol: f64, maxit: usize, log: bool) -> Autostato {
    let n3 = x0.len();
    let mut x = x0;
    proietta_fuori(&mut x, defl);
    let k = 1.0 / norm2(&x).sqrt();
    scal(&mut x, C::real(k));
    let mut hx = vec![C::ZERO; n3];
    lv.ham(&x, &mut hx);
    let mut p: Option<(Vec<C>, Vec<C>)> = None;
    let mut lam = dot(&x, &hx).re;
    let mut it = 0;
    let mut rn = 1.0;
    let (mut best, mut best_it) = (f64::MAX, 0usize);
    while it < maxit {
        it += 1;
        // residual and preconditioned direction
        let mut w = hx.clone();
        axpy(&mut w, C::real(-lam), &x);
        rn = norm2(&w).sqrt();
        if log && (it % 10 == 1) { eprintln!("    lobpcg it {it}: λ = {:.12e}  |r| = {:.3e}", lam, rn / lam.abs()); }
        if rn / lam.abs() < tol { break; }
        // floating-point floor: stop when the residual has not halved in 12 iterations
        if rn < best { best = rn; best_it = it; }
        if it > best_it + 12 && rn > 0.5 * best { break; }
        lv.precond(&mut w, sigma);
        proietta_fuori(&mut w, defl);
        let mut hw = vec![C::ZERO; n3];
        lv.ham(&w, &mut hw);
        // basis (x, w, p) orthonormalized with their H images
        let mut base: Vec<(Vec<C>, Vec<C>)> = vec![(x.clone(), hx.clone()), (w, hw)];
        if let Some(pp) = p.take() { base.push(pp); }
        let mut ort: Vec<(Vec<C>, Vec<C>)> = vec![];
        for (mut v, mut hv) in base.into_iter() {
            for _rep in 0..2 {
                for (q, hq) in ort.iter() {
                    let c = dot(q, &v);
                    axpy(&mut v, C::ZERO - c, q);
                    axpy(&mut hv, C::ZERO - c, hq);
                }
            }
            let nv = norm2(&v).sqrt();
            if nv < 1e-10 { continue; }
            scal(&mut v, C::real(1.0 / nv));
            scal(&mut hv, C::real(1.0 / nv));
            ort.push((v, hv));
        }
        let m = ort.len();
        let mut a = vec![vec![C::ZERO; m]; m];
        for i in 0..m { for j in 0..m { a[i][j] = dot(&ort[i].0, &ort[j].1); } }
        for i in 0..m { for j in i..m { let s = (a[i][j] + a[j][i].conj()).scale(0.5); a[i][j] = s; a[j][i] = s.conj(); } }
        let (ev, vv) = herm_eig(&mut a);
        let mut imin = 0;
        for i in 0..m { if ev[i] < ev[imin] { imin = i; } }
        let y: Vec<C> = (0..m).map(|i| vv[i][imin]).collect();
        let mut xn = vec![C::ZERO; n3];
        let mut hxn = vec![C::ZERO; n3];
        let mut pn = vec![C::ZERO; n3];
        let mut hpn = vec![C::ZERO; n3];
        for i in 0..m {
            axpy(&mut xn, y[i], &ort[i].0);
            axpy(&mut hxn, y[i], &ort[i].1);
            if i > 0 { axpy(&mut pn, y[i], &ort[i].0); axpy(&mut hpn, y[i], &ort[i].1); }
        }
        let nx = norm2(&xn).sqrt();
        scal(&mut xn, C::real(1.0 / nx));
        scal(&mut hxn, C::real(1.0 / nx));
        x = xn;
        hx = hxn;
        p = Some((pn, hpn));
        lam = ev[imin];
    }
    // fix the global phase: the largest component real positive
    let (mut imax, mut amax) = (0, 0.0);
    for (i, c) in x.iter().enumerate() { if c.norm2() > amax { amax = c.norm2(); imax = i; } }
    let ph = C::expi(-x[imax].arg());
    scal(&mut x, ph);
    let (er, rr) = lv.energia_regola(&x);
    Autostato { e_gen: lam, e_regola: er, resid_regola: rr, resid: rn / lam.abs(), iter: it, v: x }
}

// =====================================================================================================
//  4×4 complex helpers for the walk's Bloch matrices (positive-energy relation χ = Σ_j G_j k_j φ)
// =====================================================================================================
type M4 = [[C; 4]; 4];
fn m4_id() -> M4 { let mut m = [[C::ZERO; 4]; 4]; for i in 0..4 { m[i][i] = C::ONE; } m }
fn m4_mul(a: &M4, b: &M4) -> M4 {
    let mut c = [[C::ZERO; 4]; 4];
    for i in 0..4 { for j in 0..4 { let mut s = C::ZERO; for k in 0..4 { s += a[i][k] * b[k][j]; } c[i][j] = s; } }
    c
}
fn m4_dag(a: &M4) -> M4 { let mut c = [[C::ZERO; 4]; 4]; for i in 0..4 { for j in 0..4 { c[i][j] = a[j][i].conj(); } } c }
fn m4_inv(a: &M4) -> M4 {
    let mut m = *a;
    let mut inv = m4_id();
    for col in 0..4 {
        let mut piv = col;
        for r in col + 1..4 { if m[r][col].abs() > m[piv][col].abs() { piv = r; } }
        m.swap(col, piv);
        inv.swap(col, piv);
        let d = C::ONE / m[col][col];
        for j in 0..4 { m[col][j] = m[col][j] * d; inv[col][j] = inv[col][j] * d; }
        for r in 0..4 {
            if r == col { continue; }
            let f = m[r][col];
            for j in 0..4 { m[r][j] = m[r][j] - f * m[col][j]; inv[r][j] = inv[r][j] - f * inv[col][j]; }
        }
    }
    inv
}
fn m2_inv(a: [[C; 2]; 2]) -> [[C; 2]; 2] {
    let d = C::ONE / (a[0][0] * a[1][1] - a[0][1] * a[1][0]);
    [[a[1][1] * d, C::ZERO - a[0][1] * d], [C::ZERO - a[1][0] * d, a[0][0] * d]]
}

/// one-axis translation of a Weyl spinor in k-space: e^{−i h k σ_j} (h = −1 for L, +1 for R)
fn trasl_weyl(j: usize, h: f64, k: f64) -> [[C; 2]; 2] {
    let (c, s) = ((h * k).cos(), (h * k).sin());
    // e^{−i h k σ_j} = cos(hk) I − i sin(hk) σ_j
    let mi = C::new(0.0, -s);
    match j {
        0 => [[C::real(c), mi], [mi, C::real(c)]],
        1 => [[C::real(c), C::real(-s)], [C::real(s), C::real(c)]],
        _ => [[C::new(c, -s), C::ZERO], [C::ZERO, C::new(c, s)]],
    }
}

/// two-tick Bloch matrix W(k) = M₁ M₀ of the 3-D walk (state (L₀, L₁, R₀, R₁)), k = k a per axis
pub fn bloch_due_tick(k: [f64; 3], th: f64) -> M4 {
    let mut w = m4_id();
    for nt in 0..2u64 {
        let ord = if nt == 0 { [0, 1, 2] } else { [2, 1, 0] };
        let mut t = m4_id();
        for &j in &ord {
            let tl = trasl_weyl(j, -1.0, k[j]);
            let tr = trasl_weyl(j, 1.0, k[j]);
            let mut m = [[C::ZERO; 4]; 4];
            for a in 0..2 { for b in 0..2 { m[a][b] = tl[a][b]; m[2 + a][2 + b] = tr[a][b]; } }
            t = m4_mul(&m, &t);
        }
        let (c, s) = (th.cos(), th.sin());
        let mut r = [[C::ZERO; 4]; 4];
        for a in 0..2 {
            r[a][a] = C::real(c);
            r[2 + a][2 + a] = C::real(c);
            r[a][2 + a] = C::new(0.0, -s);
            r[2 + a][a] = C::new(0.0, -s);
        }
        w = m4_mul(&m4_mul(&r, &t), &w);
    }
    w
}

/// χ = T(k) φ on the positive band of the walk (projector from the matrix sign of (W − W†)/2i)
pub fn relazione_positiva(k: [f64; 3], th: f64) -> [[C; 2]; 2] {
    let w = bloch_due_tick(k, th);
    let wd = m4_dag(&w);
    let mut a = [[C::ZERO; 4]; 4];
    for i in 0..4 { for j in 0..4 { a[i][j] = (w[i][j] - wd[i][j]) * C::new(0.0, -0.5); } }
    // matrix sign by Newton: S ← (S + S⁻¹)/2
    let mut nrm = 0.0;
    for i in 0..4 { for j in 0..4 { nrm += a[i][j].norm2(); } }
    let mut s = a;
    for i in 0..4 { for j in 0..4 { s[i][j] = s[i][j].scale(1.0 / nrm.sqrt()); } }
    for _ in 0..100 {
        let si = m4_inv(&s);
        let mut d = 0.0;
        for i in 0..4 { for j in 0..4 { let v = (s[i][j] + si[i][j]).scale(0.5); d += (v - s[i][j]).norm2(); s[i][j] = v; } }
        if d < 1e-30 { break; }
    }
    // positive band: Im λ < 0 → sign −1: P₊ = (I − S)/2
    let mut p = [[C::ZERO; 4]; 4];
    for i in 0..4 { for j in 0..4 { p[i][j] = ((if i == j { C::ONE } else { C::ZERO }) - s[i][j]).scale(0.5); } }
    // basis: P₊ (e_a + e_{2+a})/√2
    let r2 = 0.5f64.sqrt();
    let mut u = [[C::ZERO; 2]; 4];
    for a in 0..2 { for i in 0..4 { u[i][a] = (p[i][a] + p[i][2 + a]).scale(r2); } }
    let mut phi = [[C::ZERO; 2]; 2];
    let mut chi = [[C::ZERO; 2]; 2];
    for i in 0..2 { for a in 0..2 { phi[i][a] = (u[i][a] + u[2 + i][a]).scale(r2); chi[i][a] = (u[2 + i][a] - u[i][a]).scale(r2); } }
    let pi = m2_inv(phi);
    let mut t = [[C::ZERO; 2]; 2];
    for i in 0..2 { for j in 0..2 { for q in 0..2 { t[i][j] += chi[i][q] * pi[q][j]; } } }
    t
}

/// G_j = ∂T/∂k_j at k = 0 (symmetric difference)
pub fn g_positiva(th: f64) -> [[[C; 2]; 2]; 3] {
    let dk = 1e-4;
    let mut g = [[[C::ZERO; 2]; 2]; 3];
    for j in 0..3 {
        let mut kp = [0.0; 3];
        kp[j] = dk;
        let mut km = [0.0; 3];
        km[j] = -dk;
        let (tp, tm) = (relazione_positiva(kp, th), relazione_positiva(km, th));
        for a in 0..2 { for b in 0..2 { g[j][a][b] = (tp[a][b] - tm[a][b]).scale(0.5 / dk); } }
    }
    g
}

/// Conservative reconstruction weights: the average over sub-cell s (of b) of the quartic whose cell averages on
/// the 5 coarse cells q = −2..2 are the coarse means (primitive interpolated through 6 edges). Σ_s w_q(s)/b = δ_q0.
/// order of the conservative reconstruction (odd number of coarse means; IDRO_ORD, default 5)
pub fn ordine_ricostruzione() -> usize { std::env::var("IDRO_ORD").ok().and_then(|v| v.parse().ok()).unwrap_or(5) }
/// conservative reconstruction of any odd order w: the primitive through the w+1 cell edges is interpolated by a
/// polynomial of degree w, and its differences on the b fine sub-cells give the fine means (R∘P = I exactly)
pub fn pesi_ricostruzione_ord(b: usize, w: usize) -> Vec<Vec<f64>> {
    let hw = (w / 2) as f64;
    let nodes: Vec<f64> = (0..=w).map(|e| e as f64 - hw).collect();
    let ell = |e: usize, x: f64| -> f64 {
        let mut v = 1.0;
        for (q, &xq) in nodes.iter().enumerate() { if q != e { v *= (x - xq) / (nodes[e] - xq); } }
        v
    };
    (0..b).map(|s| {
        let (x0, x1) = (s as f64 / b as f64, (s + 1) as f64 / b as f64);
        (0..w).map(|q| { let mut acc = 0.0; for e in 0..=w { if e > q { acc += ell(e, x1) - ell(e, x0); } } acc * b as f64 }).collect()
    }).collect()
}
#[allow(dead_code)]
pub fn pesi_ricostruzione(b: usize) -> Vec<[f64; 5]> {
    let nodes = [-2.0f64, -1.0, 0.0, 1.0, 2.0, 3.0];
    let ell = |e: usize, x: f64| -> f64 {
        let mut v = 1.0;
        for (q, &xq) in nodes.iter().enumerate() { if q != e { v *= (x - xq) / (nodes[e] - xq); } }
        v
    };
    (0..b).map(|s| {
        let (x0, x1) = (s as f64 / b as f64, (s + 1) as f64 / b as f64);
        let mut w = [0.0; 5];
        for q in 0..5 {
            // q index 0..4 ↔ cell offset −2..2 ; F_e depends on m_q iff edge e > cell q (edge e at node index e)
            let mut acc = 0.0;
            for e in 0..6 { if e > q { acc += ell(e, x1) - ell(e, x0); } }
            w[q] = acc * b as f64;
        }
        w
    }).collect()
}

// =====================================================================================================
//  the fine block: the zig-zag walk itself on a patch around the proton
// =====================================================================================================
/// Rotation angle of the walk in the fine block. A walk with angle θ on cells a (units ħ/mc) has, at small k,
/// 1/m_eff = a·cot θ; the bare choice θ = a gives θ cot θ: at a = h/32 = 0.268 the mass is 2.4 % too heavy. The angle is
/// DERIVED so that the block's kinetic coefficient equals the coarse rule's (read from G at the smallest k):
/// a·cot(ang) = κ_rule, and the rest phase is restored with the uniform factor e^{−i(a − ang)} per tick.
/// IDRO_MASSA=nuda → ang = a (bare walk, the old dev runs); IDRO_MASSA=fisica → κ = 1 (exact physical mass).
pub fn angolo_fine(lv: &Livello, a: f64) -> f64 {
    let modo = std::env::var("IDRO_MASSA").unwrap_or_default();
    if modo == "nuda" { return a; }
    let kappa = if modo == "fisica" { 1.0 } else { kappa_regola(lv) };
    (a / kappa).atan()
}
/// kinetic coefficient of the coarse rule at k → 0: G(k)/(k²h/2) at the two smallest k of the box (k per cell
/// 2π/N and 4π/N, mean of the three axes), extrapolated with κ(k) = κ₀ − c k² (so it does not depend on the box)
pub fn kappa_regola(lv: &Livello) -> f64 {
    let n = lv.n;
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let kap = |m: usize| {
        let g = lv.gk[lv.idx(m, 0, 0)] + lv.gk[lv.idx(0, m, 0)] + lv.gk[lv.idx(0, 0, m)];
        (g / 3.0) / ((m as f64 * k).powi(2) / (2.0 * lv.h))
    };
    (4.0 * kap(1) - kap(2)) / 3.0
}

pub struct Fine {
    /// fine cells per coarse cell, interior and halo (coarse cells per axis), patch side (fine cells)
    pub b: usize,
    pub ni: usize,
    pub nh: usize,
    pub p: usize,
    /// θ = a_f (mass 1)
    pub th: f64,
    /// rotation angle of the walk in the block (derived: its kinetic mass matches the coarse rule's, see `angolo_fine`)
    pub ang: f64,
    /// first coarse cell of the patch, relative to the nucleus cell c0 (can be negative); coarse grid side
    pub cs: i64,
    pub nc: usize,
    pub c0: usize,
    /// (L₀, L₁, R₀, R₁) per site
    pub psi: [Vec<C>; 4],
    tmp: [Vec<C>; 4],
    /// e^{−iV a_f} (point Coulomb at the site) and e^{−iVτ/2} (for the derived coarse kick)
    pub fase: Vec<C>,
    pub mezzo: Vec<C>,
    pub gj: [[[C; 2]; 2]; 3],
    /// conservative reconstruction weights (2·hw+1 coarse means → b fine values per cell)
    pub w5: Vec<Vec<f64>>,
    pub hw: usize,
    pub tick: u64,
    pub z: f64,
    /// diagnostics: spin-down large component / up in the interior at the last restriction
    pub giu_su: f64,
}

#[inline]
fn shift_val(ax: usize, h: isize, am: C, bm: C, ap: C, bp: C) -> (C, C) {
    // am, bm = values at i − h (content moving +h arrives), ap, bp = values at i + h (content moving −h arrives)
    let _ = h;
    match ax {
        2 => (am, bp),
        0 => (((am + bm) + (ap - bp)).scale(0.5), ((am + bm) - (ap - bp)).scale(0.5)),
        _ => {
            let i = C::I;
            ((am - i * bm + ap + i * bp).scale(0.5), (i * am + bm - i * ap + bp).scale(0.5))
        }
    }
}

impl Fine {
    /// patch of ni interior coarse cells (odd if the nucleus is at a coarse cell centre, even at a corner) plus
    /// nh halo cells on every side; b fine cells per coarse cell; lv gives the coarse grid
    pub fn new(lv: &Livello, ni: usize, nh: usize, b: usize, z: f64) -> Self {
        let t = ni + 2 * nh;
        assert!(if lv.off == 0.0 { t % 2 == 1 } else { t % 2 == 0 }, "patch must be centred on the nucleus");
        let p = t * b;
        let th = lv.h / b as f64;
        let cs = -((t / 2) as i64);
        let np3 = p * p * p;
        let af = th;
        let ang = angolo_fine(lv, th);
        let resto = C::expi(-(af - ang));
        let half = p as f64 / 2.0;
        let tau = lv.tau;
        let (fase, mezzo): (Vec<C>, Vec<C>) = (0..np3).into_par_iter().map(|s| {
            let (i, j, k) = (s / (p * p), (s / p) % p, s % p);
            let x = [(i as f64 + 0.5 - half) * af, (j as f64 + 0.5 - half) * af, (k as f64 + 0.5 - half) * af];
            let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            let v = pot(r, z);
            (C::expi(-v * af) * resto, C::expi(-0.5 * v * tau))
        }).unzip();
        let z4 = || [vec![C::ZERO; np3], vec![C::ZERO; np3], vec![C::ZERO; np3], vec![C::ZERO; np3]];
        Fine { b, ni, nh, p, th, ang, cs, nc: lv.n, c0: lv.c0, psi: z4(), tmp: z4(), fase, mezzo, gj: g_positiva(ang),
               w5: pesi_ricostruzione_ord(b, ordine_ricostruzione()), hw: ordine_ricostruzione() / 2, tick: 0, z, giu_su: 0.0 }
    }

    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize { (i * self.p + j) * self.p + k }
    /// fine site centre (units ħ/μc), nucleus at the origin
    pub fn x(&self, i: usize) -> f64 { (i as f64 + 0.5 - self.p as f64 / 2.0) * self.th }
    fn interno(&self, i: usize) -> bool { i >= self.nh * self.b && i < (self.nh + self.ni) * self.b }
    /// coarse grid index of patch coarse cell q (periodic)
    fn cg(&self, q: i64) -> usize { (self.c0 as i64 + self.cs + q).rem_euclid(self.nc as i64) as usize }

    /// PROLONGATION: refill the patch (halo only, or everything) from the coarse field
    pub fn prolunga(&mut self, coarse: &[C], tutto: bool) {
        let (p, b, nc) = (self.p, self.b, self.nc);
        let t = p / b;
        let hw = self.hw;
        let wn = t + 2 * hw;
        let bb = (b as f64).powf(1.5);
        // coarse means on the window (patch coarse cells −hw .. t+hw−1)
        let mut m = vec![C::ZERO; wn * wn * wn];
        for a in 0..wn { for c in 0..wn { for d in 0..wn {
            let (ga, gc, gd) = (self.cg(a as i64 - hw as i64), self.cg(c as i64 - hw as i64), self.cg(d as i64 - hw as i64));
            m[(a * wn + c) * wn + d] = coarse[(ga * nc + gc) * nc + gd].scale(1.0 / bb);
        }}}
        let w5 = &self.w5;
        // stage 1: z
        let mut a1 = vec![C::ZERO; wn * wn * p];
        a1.par_chunks_mut(p).enumerate().for_each(|(ac, row)| {
            for fz in 0..p {
                let (cz, sz) = (fz / b + hw, fz % b);
                let mut s = C::ZERO;
                for q in 0..2 * hw + 1 { s += m[ac * wn + cz + q - hw].scale(w5[sz][q]); }
                row[fz] = s;
            }
        });
        // stage 2: y
        let mut a2 = vec![C::ZERO; wn * p * p];
        a2.par_chunks_mut(p * p).enumerate().for_each(|(a, pl)| {
            for fy in 0..p {
                let (cy, sy) = (fy / b + hw, fy % b);
                for fz in 0..p {
                    let mut s = C::ZERO;
                    for q in 0..2 * hw + 1 { s += a1[(a * wn + cy + q - hw) * p + fz].scale(w5[sy][q]); }
                    pl[fy * p + fz] = s;
                }
            }
        });
        // stage 3: x → φ (spin up)
        let mut phi = vec![C::ZERO; p * p * p];
        phi.par_chunks_mut(p * p).enumerate().for_each(|(fx, pl)| {
            let (cx, sx) = (fx / b + hw, fx % b);
            for s in 0..p * p {
                let mut v = C::ZERO;
                for q in 0..2 * hw + 1 { v += a2[(cx + q - hw) * p * p + s].scale(w5[sx][q]); }
                pl[s] = v;
            }
        });
        // small component χ = Σ_j G_j (−i ∂_j φ), φ = (f, 0); then (L, R) = ((φ − χ), (φ + χ))/√2
        let gj = self.gj;
        let r2 = 0.5f64.sqrt();
        let (ni, nh) = (self.ni, self.nh);
        let lo = nh * b;
        let hi = (nh + ni) * b;
        let ptrs: [SendPtr; 4] = [SendPtr(self.psi[0].as_mut_ptr()), SendPtr(self.psi[1].as_mut_ptr()),
                                  SendPtr(self.psi[2].as_mut_ptr()), SendPtr(self.psi[3].as_mut_ptr())];
        (0..p).into_par_iter().for_each(|i| {
            let ptrs = ptrs;
            for j in 0..p {
                for k in 0..p {
                    let inside = i >= lo && i < hi && j >= lo && j < hi && k >= lo && k < hi;
                    if inside && !tutto { continue; }
                    let at = |a: usize, c: usize, d: usize| phi[(a * p + c) * p + d];
                    let der = |ax: usize| -> C {
                        let (ii, jj, kk) = (i as isize, j as isize, k as isize);
                        let get = |di: isize| -> C {
                            let (a, c, d) = match ax { 0 => (ii + di, jj, kk), 1 => (ii, jj + di, kk), _ => (ii, jj, kk + di) };
                            let cl = |v: isize| v.clamp(0, p as isize - 1) as usize;
                            at(cl(a), cl(c), cl(d))
                        };
                        (get(1) - get(-1)).scale(0.5)
                    };
                    let f = at(i, j, k);
                    let mut chi = [C::ZERO; 2];
                    for ax in 0..3 {
                        let d = der(ax) * C::new(0.0, -1.0);
                        chi[0] += gj[ax][0][0] * d;
                        chi[1] += gj[ax][1][0] * d;
                    }
                    let s = (i * p + j) * p + k;
                    unsafe {
                        *ptrs[0].0.add(s) = (f - chi[0]).scale(r2);
                        *ptrs[1].0.add(s) = (C::ZERO - chi[1]).scale(r2);
                        *ptrs[2].0.add(s) = (f + chi[0]).scale(r2);
                        *ptrs[3].0.add(s) = chi[1].scale(r2);
                    }
                }
            }
        });
    }

    /// RESTRICTION of the interior: coarse value = block mean of the large component φ↑ = (L₀ + R₀)/√2 × b^{3/2};
    /// with `calcio`, the block mean of e^{−iVτ/2} φ↑ (the coarse step's half kick read from the fine field)
    pub fn restringi(&mut self, calcio: bool) -> Vec<(usize, C)> {
        let (p, b) = (self.p, self.b);
        let r2 = 0.5f64.sqrt();
        let bb = (b as f64).powf(1.5);
        let mut out = vec![];
        let mut su = 0.0;
        let mut giu = 0.0;
        for qa in 0..self.ni { for qb in 0..self.ni { for qc in 0..self.ni {
            let (sa, sb, sc) = ((self.nh + qa) * b, (self.nh + qb) * b, (self.nh + qc) * b);
            let (s, d): (C, C) = (sa..sa + b).into_par_iter().map(|i| {
                let mut s = C::ZERO;
                let mut d = C::ZERO;
                for j in sb..sb + b { for k in sc..sc + b {
                    let q = (i * p + j) * p + k;
                    let mut u = (self.psi[0][q] + self.psi[2][q]).scale(r2);
                    let mut w = (self.psi[1][q] + self.psi[3][q]).scale(r2);
                    if calcio { u = u * self.mezzo[q]; w = w * self.mezzo[q]; }
                    s += u;
                    d += w;
                }}
                (s, d)
            }).reduce(|| (C::ZERO, C::ZERO), |x, y| (x.0 + y.0, x.1 + y.1));
            su += s.norm2();
            giu += d.norm2();
            let g = (self.cg(self.nh as i64 + qa as i64), self.cg(self.nh as i64 + qb as i64), self.cg(self.nh as i64 + qc as i64));
            out.push(((g.0 * self.nc + g.1) * self.nc + g.2, s.scale(1.0 / bb)));
        }}}
        if !calcio { self.giu_su = (giu / su.max(1e-300)).sqrt(); }
        out
    }

    /// one tick of the zig-zag walk; t = 1..=b inside the coarse step: only the light cone of the interior
    /// ([t, p − t) per axis after this tick) is computed
    pub fn tick(&mut self, t: usize) {
        let p = self.p;
        let ord = if self.tick % 2 == 0 { [0usize, 1, 2] } else { [2, 1, 0] };
        let mut fatto = [false; 3];
        for &ax in &ord {
            fatto[ax] = true;
            let mut lo = [0usize; 3];
            let mut hi = [0usize; 3];
            for a in 0..3 {
                if fatto[a] { lo[a] = t; hi[a] = p - t; } else { lo[a] = t - 1; hi[a] = p - t + 1; }
            }
            let src = &self.psi;
            let dst: [SendPtr; 4] = [SendPtr(self.tmp[0].as_mut_ptr()), SendPtr(self.tmp[1].as_mut_ptr()),
                                     SendPtr(self.tmp[2].as_mut_ptr()), SendPtr(self.tmp[3].as_mut_ptr())];
            let stride = match ax { 0 => p * p, 1 => p, _ => 1 };
            (lo[0]..hi[0]).into_par_iter().for_each(|i| {
                let dst = dst;
                for j in lo[1]..hi[1] {
                    for k in lo[2]..hi[2] {
                        let s = (i * p + j) * p + k;
                        // L moves with h = −1: content moving +h comes from s − h·stride = s + stride
                        let (lm, lp) = (s + stride, s - stride);
                        let (a, bq) = shift_val(ax, -1, src[0][lm], src[1][lm], src[0][lp], src[1][lp]);
                        let (rm, rp) = (s - stride, s + stride);
                        let (c, d) = shift_val(ax, 1, src[2][rm], src[3][rm], src[2][rp], src[3][rp]);
                        unsafe {
                            *dst[0].0.add(s) = a;
                            *dst[1].0.add(s) = bq;
                            *dst[2].0.add(s) = c;
                            *dst[3].0.add(s) = d;
                        }
                    }
                }
            });
            std::mem::swap(&mut self.psi, &mut self.tmp);
        }
        // rotation (mass flip) and the Coulomb phase
        let (ct, st) = (self.ang.cos(), self.ang.sin());
        let ptr: [SendPtr; 4] = [SendPtr(self.psi[0].as_mut_ptr()), SendPtr(self.psi[1].as_mut_ptr()),
                                 SendPtr(self.psi[2].as_mut_ptr()), SendPtr(self.psi[3].as_mut_ptr())];
        let fase = &self.fase;
        let mis = C::new(0.0, -st);
        (t..p - t).into_par_iter().for_each(|i| {
            let ptr = ptr;
            for j in t..p - t {
                for k in t..p - t {
                    let s = (i * p + j) * p + k;
                    unsafe {
                        let (l0, l1, r0, r1) = (*ptr[0].0.add(s), *ptr[1].0.add(s), *ptr[2].0.add(s), *ptr[3].0.add(s));
                        let f = fase[s];
                        *ptr[0].0.add(s) = (l0.scale(ct) + mis * r0) * f;
                        *ptr[1].0.add(s) = (l1.scale(ct) + mis * r1) * f;
                        *ptr[2].0.add(s) = (r0.scale(ct) + mis * l0) * f;
                        *ptr[3].0.add(s) = (r1.scale(ct) + mis * l1) * f;
                    }
                }
            }
        });
        self.tick += 1;
    }

    /// coarse grid indices of the interior cells
    pub fn celle_coperte(&self) -> Vec<usize> {
        let mut v = vec![];
        for a in 0..self.ni { for b in 0..self.ni { for c in 0..self.ni {
            let g = (self.cg((self.nh + a) as i64), self.cg((self.nh + b) as i64), self.cg((self.nh + c) as i64));
            v.push((g.0 * self.nc + g.1) * self.nc + g.2);
        }}}
        v
    }

    /// b ticks = one coarse step
    pub fn passo(&mut self) { for t in 1..=self.b { self.tick(t); } }

    /// multiply the whole patch (all four components) by a real factor (reflux)
    pub fn scala_interno(&mut self, k: f64) {
        for c in 0..4 { self.psi[c].par_iter_mut().for_each(|x| *x = x.scale(k)); }
    }

    /// probability in the interior (fine), and the scalar density ψ†βψ = 2 Re(L†R) summed over the interior
    pub fn interno_norma_beta(&self) -> (f64, f64) {
        let p = self.p;
        let (lo, hi) = (self.nh * self.b, (self.nh + self.ni) * self.b);
        (lo..hi).into_par_iter().map(|i| {
            let (mut n, mut bt) = (0.0, 0.0);
            for j in lo..hi { for k in lo..hi {
                let s = (i * p + j) * p + k;
                let (l0, l1, r0, r1) = (self.psi[0][s], self.psi[1][s], self.psi[2][s], self.psi[3][s]);
                n += l0.norm2() + l1.norm2() + r0.norm2() + r1.norm2();
                bt += 2.0 * ((l0.conj() * r0).re + (l1.conj() * r1).re);
            }}
            (n, bt)
        }).reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1))
    }
}

// =====================================================================================================
//  the atom: coarse learned level + (optional) fine zig-zag block around the proton
// =====================================================================================================
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Calcio { DalFine, Grosso }

pub struct Atomo {
    pub lv: Livello,
    pub psi: Vec<C>,
    pub fine: Option<Fine>,
    pub passi: u64,
    pub calcio: Calcio,
    /// second-order covariance of the cell mean, mean(Vψ) − V̄ψ̄ ≈ (h²/12) ∇V·∇ψ, on the cells not covered by the
    /// fine block (derived from the descriptor = block mean; off by default)
    pub covarianza: bool,
    /// charge balance at the seam (IDRO_REFLUX=1) and the last correction factor − 1 (diagnostic)
    pub reflux: bool,
    pub ultimo_reflux: f64,
}

/// (h²/12) Σ_j ∂_jV ∂_jf with central differences, f ← f − i (τ/2) (that), skipping masked cells
pub fn calcio_covarianza(lv: &Livello, f: &mut [C], coperte: &[usize]) {
    let n = lv.n;
    let src = f.to_vec();
    let w = |i: usize, d: isize| -> usize { ((i as isize + d).rem_euclid(n as isize)) as usize };
    let fac = C::new(0.0, -0.5 * lv.tau / 48.0);
    f.par_iter_mut().enumerate().for_each(|(s, out)| {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        let mut acc = C::ZERO;
        for ax in 0..3 {
            let (sp, sm) = match ax {
                0 => ((w(i, 1) * n + j) * n + k, (w(i, -1) * n + j) * n + k),
                1 => ((i * n + w(j, 1)) * n + k, (i * n + w(j, -1)) * n + k),
                _ => ((i * n + j) * n + w(k, 1), (i * n + j) * n + w(k, -1)),
            };
            acc += (src[sp] - src[sm]).scale(lv.v[sp] - lv.v[sm]);
        }
        *out += fac * acc;
    });
    for &s in coperte { f[s] = src[s]; }
}

impl Atomo {
    pub fn new(lv: Livello, psi: Vec<C>, fine: Option<Fine>, calcio: Calcio) -> Self {
        let reflux = std::env::var("IDRO_REFLUX").map(|v| v != "0").unwrap_or(false);
        let mut a = Atomo { lv, psi, fine, passi: 0, calcio, covarianza: false, reflux, ultimo_reflux: 0.0 };
        if let Some(f) = a.fine.as_mut() {
            f.prolunga(&a.psi, true);
            for (s, v) in f.restringi(false) { a.psi[s] = v; }
        }
        a
    }

    /// one coarse step (τ = h): halo ← coarse (time n); coarse ← rule (covered cells' kick from the fine field);
    /// fine ← b ticks; covered coarse cells ← fine (time n + 1)
    pub fn passo(&mut self) { self.passo_osserva(&mut |_f: &Fine, _t: usize| {}); }

    /// the same step, calling `oss(fine, t)` after every fine tick t = 1..b (for the zoom export)
    pub fn passo_osserva(&mut self, oss: &mut dyn FnMut(&Fine, usize)) {
        let mut g = self.psi.clone();
        g.par_iter_mut().zip(self.lv.calcio.par_iter()).for_each(|(a, k)| *a = *a * *k);
        let coperte: Vec<usize> = match self.fine.as_ref() { Some(f) => f.celle_coperte(), None => vec![] };
        if self.covarianza { calcio_covarianza(&self.lv, &mut g, &coperte); }
        if let Some(f) = self.fine.as_mut() {
            f.prolunga(&self.psi, false);
            if self.calcio == Calcio::DalFine { for (s, v) in f.restringi(true) { g[s] = v; } }
        }
        self.lv.fft.tre(&mut g, false);
        g.par_iter_mut().zip(self.lv.molt.par_iter()).for_each(|(a, m)| *a = *a * *m);
        self.lv.fft.tre(&mut g, true);
        g.par_iter_mut().zip(self.lv.calcio.par_iter()).for_each(|(a, k)| *a = *a * *k);
        if self.covarianza { calcio_covarianza(&self.lv, &mut g, &coperte); }
        // reflux: the charge that crosses the interface is the one the coarse rule moves. The covered cells get, in total,
        // exactly the charge Q(before) − Q(uncovered, after); the fine block is rescaled by the same factor, so the
        // fine region keeps its own shape (cusp, jitter) but cannot create or lose charge at the seam
        let q_prima = if self.reflux { norm2(&self.psi) } else { 0.0 };
        if let Some(f) = self.fine.as_mut() {
            for t in 1..=f.b { f.tick(t); oss(f, t); }
            let r = f.restringi(false);
            if self.reflux {
                let q_cov_g: f64 = r.iter().map(|(s, _)| g[*s].norm2()).sum();
                let q_unc = norm2(&g) - q_cov_g;
                let q_r: f64 = r.iter().map(|(_, v)| v.norm2()).sum();
                let fatt = ((q_prima - q_unc).max(0.0) / q_r).sqrt();
                self.ultimo_reflux = fatt - 1.0;
                for (s, v) in r { g[s] = v.scale(fatt); }
                f.scala_interno(fatt);
            } else {
                for (s, v) in r { g[s] = v; }
            }
        }
        self.psi = g;
        self.passi += 1;
    }

    pub fn t(&self) -> f64 { self.passi as f64 * self.lv.tau }
}


// =====================================================================================================
//  export helpers (float16 little endian, as the other exports in ripresa/dati)
// =====================================================================================================
pub fn f16_bits(x: f32) -> u16 {
    let b = x.to_bits();
    let sign = ((b >> 16) & 0x8000) as u16;
    let exp = ((b >> 23) & 0xff) as i32;
    let man = b & 0x7f_ffff;
    if exp == 255 { return sign | 0x7c00 | if man != 0 { 0x200 } else { 0 }; }
    let e = exp - 127 + 15;
    if e >= 31 { return sign | 0x7c00; }
    if e <= 0 {
        if e < -10 { return sign; }
        let m = (man | 0x80_0000) >> (1 - e);
        let r = (m >> 13) + ((m >> 12) & 1);
        return sign | r as u16;
    }
    let mut v = ((e as u32) << 10) | (man >> 13);
    if (man >> 12) & 1 == 1 { v += 1; }
    sign | v as u16
}

pub fn scrivi_f16(w: &mut impl std::io::Write, data: &[f64], scala: f64) {
    let mut buf = Vec::with_capacity(data.len() * 2);
    for &x in data { buf.extend_from_slice(&f16_bits((x * scala) as f32).to_le_bytes()); }
    w.write_all(&buf).unwrap();
}

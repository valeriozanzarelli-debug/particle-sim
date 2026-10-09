//! The BASE rule v2 in 3-D — the default rule of the program for the
//! single-photon experiments (free flight, double slit, polarizers, Airy disk,
//! Arago spot).
//!
//! Port of `lab/base_rule.py` (classes `Mask3`, `Obj3`, `Scene3D`, `Photons`
//! with dim = 3) and of the 3-D field sums of `lab/base_kernel/src/lib.rs`.
//! Spec in Italian: docs/REGOLA_BASE.md (§1 the rule, §2 how many balls,
//! §3 the objects). Frozen parameters (v2): η = 2, ε = 10⁻³, δ = 0.02 λ,
//! Dormand–Prince 5(4) with 10⁻⁴ λ per step and steps ≤ 2 λ, edge zone λ/2.
//!
//! A photon = many WAVE BALLS + ONE PHOTON BALL.
//!   • Every object (planes ⊥ z) re-emits the wave that reaches it as a fan:
//!     origins on a regular grid ds0 = λ/(2η) over the part where
//!     |B| ≥ ε·max|B| (B = T·E_t, the outgoing boundary value), directions on
//!     the Nyquist grid d(sinθ) = λ/(ηW) over the whole forward half-space;
//!     evanescent levels ℓ = 1..ℓmax (ds0/2^ℓ, heard only for R < R0·2^(2−ℓ))
//!     are made on demand by tiles and kept (they are the same for every photon).
//!   • The photon ball is born on the source with probability |B|², hears the
//!     fan of the last object crossed, sums the exact vector Rayleigh–Sommerfeld
//!     wavelets (E_t = −2 Σ B_t ∂_z g, E_z = 2 Σ B_t·∇_t g, g = e^{ikR}/4πR),
//!     computes j = Σ_a Im(E_a* ∇E_a) and moves at c along j. At an absorbing
//!     object: local coin min(1, Φ_out/Φ_in); within λ/2 of an edge it keeps
//!     its rank in the flux across the edge (v2). Screen → click.
//!
//! What is new here with respect to the Python reference (same arithmetic,
//! only faster): the boundary values B of an object are read from the
//! previous fan on a regular lattice; when the planes are ≥ 2·R0 apart only
//! level 0 is heard there and the sum over the origins (all on a lattice) is a
//! discrete convolution, done with an FFT (`lattice_tangential`; tested
//! against the direct sum in `base3d_validation.rs`).
//!
//! The number of wave balls the rule uses is huge (10⁷–10¹² per photon, see
//! `Fan::ndir` and `Bench3::ball_counts`): the guidance only needs the origins
//! (the directions decide who hears whom, and in 3-D the photon hears the fan
//! of the last plane crossed). What is DRAWN is a random subset (`render3d.rs`).

use std::collections::{HashMap, VecDeque};
use std::f64::consts::PI;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, RwLock};

use bevy::math::Vec3;
use rayon::prelude::*;

use crate::pilot::{hash2, uniform, Z};

pub type J2 = [Z; 2];
/// 2×2 complex (Jones) matrix, row-major.
pub type M2 = [[Z; 2]; 2];

// ─── The rule (frozen, v2) ────────────────────────────────────────────────────

/// Parameters of the rule, in units of λ (docs/REGOLA_BASE.md §2.3).
#[derive(Clone, Copy, Debug)]
pub struct RuleV2 {
    pub eta: f64,
    pub eps: f64,
    pub delta: f64,
    pub atol: f64,
    pub hmax: f64,
    pub hmin: f64,
    pub edge_zone: f64,
}

pub const RULE: RuleV2 = RuleV2 { eta: 2.0, eps: 1e-3, delta: 0.02, atol: 1e-4, hmax: 2.0, hmin: 1e-6, edge_zone: 0.5 };

/// Largest number of fine-level origins kept in the tile caches (all fans).
const TILE_CACHE_CAP: usize = 250_000;
/// Origins per partial sum of level 0 (fixed: the result does not depend on threads).
const CHUNK: usize = 512;
/// From this many level-0 origins on, the partial sums run in parallel.
const PAR_MIN: usize = 4096;
/// Points kept in the photon ball's trail.
pub const TRAIL_POINTS: usize = 900;
pub const MAX_STEPS: usize = 20_000;

// ─── Small complex helpers ────────────────────────────────────────────────────

/// sin and cos of x (|x| < 2⁵⁰): Cody–Waite reduction by π/2 and the fdlibm
/// minimax kernels on [−π/4, π/4]. Error ≲ 2·10⁻¹⁶·max(1, |x|·10⁻¹⁶·…) — the
/// same as the library's to a few ulp (tested in base3d_validation.rs), about
/// 4× faster (it is the inner loop of every wave-ball sum).
#[inline(always)]
pub fn sincos(x: f64) -> (f64, f64) {
    const TWO_OVER_PI: f64 = 0.636_619_772_367_581_343_08;
    const P1: f64 = 1.570_796_326_734_125_614_17e0;
    const P2: f64 = 6.077_100_506_303_965_976_60e-11;
    const P3: f64 = 2.022_266_248_795_950_631_54e-21;
    const MAGIC: f64 = 6_755_399_441_055_744.0; // 1.5·2⁵²: round to nearest integer
    let q = (x * TWO_OVER_PI + MAGIC) - MAGIC;
    let r = ((x - q * P1) - q * P2) - q * P3;
    let z = r * r;
    let s = r + r * z * (-1.666_666_666_666_663_243_48e-1 + z * (8.333_333_333_322_489_461_24e-3
        + z * (-1.984_126_982_985_794_931_34e-4 + z * (2.755_731_370_707_006_767_89e-6
        + z * (-2.505_076_025_340_686_341_95e-8 + z * 1.589_690_995_211_550_102_21e-10)))));
    let c = 1.0 - 0.5 * z + z * z * (4.166_666_666_666_660_190_37e-2 + z * (-1.388_888_888_887_410_957_49e-3
        + z * (2.480_158_728_947_672_941_78e-5 + z * (-2.755_731_435_139_066_330_35e-7
        + z * (2.087_572_321_298_174_827_90e-9 + z * -1.135_964_755_778_819_482_65e-11)))));
    match (q as i64) & 3 {
        0 => (s, c),
        1 => (c, -s),
        2 => (-s, -c),
        _ => (-c, s),
    }
}

#[inline(always)]
fn expi_fast(t: f64) -> Z {
    let (s, c) = sincos(t);
    Z::new(c, s)
}

#[inline]
fn cdiv(a: Z, b: Z) -> Z {
    let d = b.norm2();
    Z::new((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d)
}
#[inline]
fn csqrt(a: Z) -> Z {
    let r = a.abs();
    let re = (0.5 * (r + a.re)).max(0.0).sqrt();
    let im = (0.5 * (r - a.re)).max(0.0).sqrt();
    Z::new(re, if a.im < 0.0 { -im } else { im })
}
#[inline]
fn cexp(a: Z) -> Z { Z::expi(a.im).s(a.re.exp()) }
#[inline]
fn j2n(a: J2) -> f64 { a[0].norm2() + a[1].norm2() }
#[inline]
fn mat_apply(m: &M2, v: J2) -> J2 {
    [m[0][0] * v[0] + m[0][1] * v[1], m[1][0] * v[0] + m[1][1] * v[1]]
}
/// Projector a a† (linear polarizer with real axis a).
pub fn projector(a: J2) -> M2 {
    [[a[0] * a[0].conj(), a[0] * a[1].conj()], [a[1] * a[0].conj(), a[1] * a[1].conj()]]
}
pub fn identity() -> M2 { [[Z::new(1.0, 0.0), Z::ZERO], [Z::ZERO, Z::new(1.0, 0.0)]] }

// ─── Levels ───────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Lev {
    pub k: f64,
    pub r0: f64,
    pub ds0: f64,
    pub lmax: usize,
}

impl Lev {
    pub fn new(lam: f64, r: &RuleV2) -> Self {
        let r0 = (1.0 / r.eps).ln() * lam / (2.0 * PI);
        let ds0 = lam / (2.0 * r.eta);
        let lmax = ((r0 / (r.delta * lam)).log2().ceil() as usize) + 1;
        Lev { k: 2.0 * PI / lam, r0, ds0, lmax }
    }
    /// weight of level l for an origin at distance R (smoothstep in log2 R; Σ_l = 1)
    #[inline]
    pub fn w(&self, l: usize, r: f64) -> f64 {
        let u = (self.r0 / r).log2();
        if u <= -1.0 {
            return if l == 0 { 1.0 } else { 0.0 };
        }
        let c = u.ceil();
        let f = u - (c - 1.0);
        let s = f * f * (3.0 - 2.0 * f);
        let lm = self.lmax as i64;
        let a = (c as i64).min(lm);
        let b = (c as i64 + 1).min(lm);
        let li = l as i64;
        let mut w = 0.0;
        if a == li { w += 1.0 - s; }
        if b == li { w += s; }
        w
    }
    #[inline]
    pub fn rmax(&self, l: usize) -> f64 { if l == 0 { f64::INFINITY } else { self.r0 * 2f64.powi(2 - l as i32) } }
    #[inline]
    pub fn ds(&self, l: usize) -> f64 { self.ds0 / 2f64.powi(l as i32) }
}

// ─── Masks and objects ────────────────────────────────────────────────────────

/// Where an object transmits (planes ⊥ z): whole plane, inside / outside a
/// circle (cx, cy, r), union of rectangles (x0, x1, y0, y1).
#[derive(Clone, Debug)]
pub enum Mask {
    Full,
    CircIn(f64, f64, f64),
    CircOut(f64, f64, f64),
    Rects(Vec<[f64; 4]>),
}

impl Mask {
    #[inline]
    pub fn inside(&self, x: f64, y: f64) -> bool {
        match self {
            Mask::Full => true,
            Mask::CircIn(cx, cy, r) => (x - cx).powi(2) + (y - cy).powi(2) < r * r,
            Mask::CircOut(cx, cy, r) => (x - cx).powi(2) + (y - cy).powi(2) >= r * r,
            Mask::Rects(v) => v.iter().any(|q| x >= q[0] && x < q[1] && y >= q[2] && y < q[3]),
        }
    }

    /// (area inside, centroid) of the cell [x0,x1]×[y0,y1]
    fn cell(&self, x0: f64, x1: f64, y0: f64, y1: f64) -> (f64, f64, f64) {
        let full = ((x1 - x0) * (y1 - y0), 0.5 * (x0 + x1), 0.5 * (y0 + y1));
        match self {
            Mask::Full => full,
            Mask::Rects(v) => {
                let (mut a, mut mx, mut my) = (0.0, 0.0, 0.0);
                for q in v {
                    let l = x0.max(q[0]);
                    let r = x1.min(q[1]);
                    let b = y0.max(q[2]);
                    let t = y1.min(q[3]);
                    if r > l && t > b {
                        let ar = (r - l) * (t - b);
                        a += ar;
                        mx += ar * 0.5 * (l + r);
                        my += ar * 0.5 * (b + t);
                    }
                }
                if a > 0.0 { (a, mx / a, my / a) } else { (0.0, 0.0, 0.0) }
            }
            Mask::CircIn(cx, cy, r) | Mask::CircOut(cx, cy, r) => {
                let inside = matches!(self, Mask::CircIn(..));
                let dmax = ((x0 - cx).abs().max((x1 - cx).abs())).hypot((y0 - cy).abs().max((y1 - cy).abs()));
                let dmin = (cx.clamp(x0, x1) - cx).hypot(cy.clamp(y0, y1) - cy);
                if dmax <= *r { return if inside { full } else { (0.0, 0.0, 0.0) }; }
                if dmin >= *r { return if inside { (0.0, 0.0, 0.0) } else { full }; }
                // boundary cell: 8×8 sub-samples
                let n = 8;
                let (mut c, mut mx, mut my) = (0usize, 0.0, 0.0);
                for i in 0..n {
                    let x = x0 + (i as f64 + 0.5) * (x1 - x0) / n as f64;
                    for j in 0..n {
                        let y = y0 + (j as f64 + 0.5) * (y1 - y0) / n as f64;
                        if self.inside(x, y) {
                            c += 1;
                            mx += x;
                            my += y;
                        }
                    }
                }
                if c == 0 { (0.0, 0.0, 0.0) } else { (full.0 * c as f64 / (n * n) as f64, mx / c as f64, my / c as f64) }
            }
        }
    }

    pub fn hull(&self) -> Option<[f64; 4]> {
        match self {
            Mask::Full | Mask::CircOut(..) => None,
            Mask::CircIn(cx, cy, r) => Some([cx - r, cx + r, cy - r, cy + r]),
            Mask::Rects(v) => Some([
                v.iter().map(|q| q[0]).fold(f64::INFINITY, f64::min),
                v.iter().map(|q| q[1]).fold(f64::NEG_INFINITY, f64::max),
                v.iter().map(|q| q[2]).fold(f64::INFINITY, f64::min),
                v.iter().map(|q| q[3]).fold(f64::NEG_INFINITY, f64::max),
            ]),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Gaussian beam (or plane wave) on a plane, possibly through an aperture mask.
    Source,
    /// Absorbing screen with openings (slits, circular aperture).
    Aperture,
    /// Absorbing linear polarizer (index into `Lab::polarizers`).
    Polarizer(usize),
    /// Detection screen (click).
    Screen,
    /// Passive plane: records the crossing, touches nothing (tests).
    Plane,
}

#[derive(Clone, Debug)]
pub struct SrcBeam {
    /// Waist radius (∞ = plane wave), optional different waist along y.
    pub w0: f64,
    pub w0y: Option<f64>,
    /// Distance of the waist downstream of the source plane.
    pub zw: f64,
    pub jones: J2,
}

#[derive(Clone, Debug)]
pub struct Obj {
    pub name: String,
    pub kind: Kind,
    pub z: f64,
    pub mask: Mask,
    /// Transmission Jones matrix (None: does not transmit).
    pub t: Option<M2>,
    /// Optional transmission per rectangle of a `Rects` mask (which-path polarizers).
    pub rect_t: Vec<Option<M2>>,
    pub absorbing: bool,
    pub src: Option<SrcBeam>,
}

impl Obj {
    pub fn new(name: &str, kind: Kind, z: f64, mask: Mask) -> Self {
        let t = matches!(kind, Kind::Source | Kind::Aperture | Kind::Polarizer(_)).then(identity);
        let absorbing = matches!(kind, Kind::Source | Kind::Aperture | Kind::Polarizer(_));
        Obj { name: name.into(), kind, z, mask, t, rect_t: Vec::new(), absorbing, src: None }
    }

    /// Gaussian beam (paraxial, waist zw downstream) × Jones (Python `Obj3.src_B`).
    pub fn src_b(&self, lam: f64, x: f64, y: f64) -> J2 {
        let s = self.src.as_ref().expect("source");
        let zrel = -s.zw;
        let u = match s.w0y {
            Some(wy) => gauss_beam(lam, x * x, s.w0, zrel, false) * gauss_beam(lam, y * y, wy, zrel, false),
            None => gauss_beam(lam, x * x + y * y, s.w0, zrel, true),
        };
        [u * s.jones[0], u * s.jones[1]]
    }

    /// B = T · E_t (T of the rectangle the point is in, or the nearest one).
    fn transmit(&self, x: f64, y: f64, e: J2) -> J2 {
        let mut t = self.t.as_ref();
        if let (Mask::Rects(v), false) = (&self.mask, self.rect_t.is_empty()) {
            let d = |q: &[f64; 4]| (q[0] - x).max(x - q[1]).max(0.0).hypot((q[2] - y).max(y - q[3]).max(0.0));
            let best = (0..v.len()).min_by(|&a, &b| d(&v[a]).partial_cmp(&d(&v[b])).unwrap()).unwrap();
            if let Some(Some(m)) = self.rect_t.get(best) { t = Some(m); }
        }
        match t { Some(m) => mat_apply(m, e), None => [Z::ZERO; 2] }
    }
}

/// Paraxial Gaussian beam on a plane at `zrel` from its waist; r2 = transverse
/// distance²; 3-D (`dim3`) or cylindrical (1 transverse axis). w0 = ∞: plane wave.
pub fn gauss_beam(lam: f64, r2: f64, w0: f64, zrel: f64, dim3: bool) -> Z {
    if !w0.is_finite() {
        return Z::new(1.0, 0.0);
    }
    let k = 2.0 * PI / lam;
    let zr = PI * w0 * w0 / lam;
    let q0 = Z::new(0.0, -zr);
    let q = Z::new(zrel, -zr);
    let ratio = cdiv(q0, q);
    let amp = if dim3 { ratio } else { csqrt(ratio) };
    // exp(i k r² / (2q))
    let e = cdiv(Z::new(0.0, k * r2 / 2.0), q);
    amp * cexp(e)
}

// ─── Fans ─────────────────────────────────────────────────────────────────────

/// One wave-ball origin: position on the plane and amplitude a = B · cell area.
#[derive(Clone, Copy, Debug)]
pub struct Org {
    pub p: [f64; 2],
    pub a: J2,
}

#[inline]
fn cr_weights(f: f64) -> [f64; 4] {
    let f2 = f * f;
    let f3 = f2 * f;
    [(-f3 + 2.0 * f2 - f) * 0.5, (3.0 * f3 - 5.0 * f2 + 2.0) * 0.5, (-3.0 * f3 + 4.0 * f2 + f) * 0.5, (f3 - f2) * 0.5]
}

const TILE: i64 = 16;
type TileKey = (usize, i64, i64);

/// Fine-level origins, made on demand by tiles of TILE×TILE cells and kept.
#[derive(Default)]
struct TileCache {
    map: RwLock<HashMap<TileKey, Arc<Vec<Org>>>>,
}

/// All the wave balls emitted by one object (one plane).
pub struct Fan {
    pub obj: usize,
    pub z0: f64,
    /// Box of the origins (x0, x1, y0, y1).
    pub g: [f64; 4],
    pub mask: Mask,
    sx0: f64,
    sy0: f64,
    sds: f64,
    nx: usize,
    ny: usize,
    /// B sampled at λ/16 (Catmull–Rom interpolation for the fine levels).
    samp: Vec<J2>,
    bthr: f64,
    /// Level-0 origins.
    pub l0: Vec<Org>,
    /// The same, as arrays (padded to a multiple of 4) for the far-field sum.
    soa: Soa,
    cache: TileCache,
    /// Extent W (x, y) and directions per origin (Nyquist, whole half-space).
    pub w: [f64; 2],
    pub ndir: usize,
    /// Light-clock when the on-axis front reaches the object, and the plane
    /// where these balls end (next object or the screen): for drawing.
    pub t_emit: f64,
    pub z_kill: f64,
}

impl Fan {
    fn interp(&self, x: f64, y: f64) -> J2 {
        let u = (x - self.sx0) / self.sds;
        let v = (y - self.sy0) / self.sds;
        let iu = u.floor();
        let iv = v.floor();
        let wu = cr_weights(u - iu);
        let wv = cr_weights(v - iv);
        let iu = iu as i64;
        let iv = iv as i64;
        let mut b = [Z::ZERO; 2];
        for n in 0..4 {
            let jy = iv - 1 + n as i64;
            if jy < 0 || jy >= self.ny as i64 { continue; }
            for m in 0..4 {
                let jx = iu - 1 + m as i64;
                if jx < 0 || jx >= self.nx as i64 { continue; }
                let w = wu[m] * wv[n];
                let s = self.samp[jy as usize * self.nx + jx as usize];
                b[0] += s[0].s(w);
                b[1] += s[1].s(w);
            }
        }
        b
    }

    /// origins of the cells i0..i1 × j0..j1 (cell (i, j) = [g0 + i ds, …] × [g2 + j ds, …], clipped to g)
    fn origins_cells(&self, ds: f64, i0: i64, i1: i64, j0: i64, j1: i64, out: &mut Vec<Org>) {
        for j in j0.max(0)..j1 {
            let y0 = self.g[2] + j as f64 * ds;
            if y0 >= self.g[3] { break; }
            let y1 = (y0 + ds).min(self.g[3]);
            for i in i0.max(0)..i1 {
                let x0 = self.g[0] + i as f64 * ds;
                if x0 >= self.g[1] { break; }
                let x1 = (x0 + ds).min(self.g[1]);
                let (area, cx, cy) = self.mask.cell(x0, x1, y0, y1);
                if area <= 0.0 { continue; }
                let bb = self.interp(cx, cy);
                if j2n(bb).sqrt() < self.bthr { continue; }
                out.push(Org { p: [cx, cy], a: [bb[0].s(area), bb[1].s(area)] });
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn tiles(&self, cached: &AtomicUsize, l: usize, ds: f64, xa: f64, xb: f64, ya: f64, yb: f64,
             out: &mut Vec<Arc<Vec<Org>>>) {
        let xlo = xa.max(self.g[0]);
        let xhi = xb.min(self.g[1]);
        let ylo = ya.max(self.g[2]);
        let yhi = yb.min(self.g[3]);
        if xhi <= xlo || yhi <= ylo { return; }
        let t0 = (((xlo - self.g[0]) / ds).floor() as i64).max(0) / TILE;
        let t1 = (((xhi - self.g[0]) / ds).ceil() as i64 - 1).max(0) / TILE;
        let u0 = (((ylo - self.g[2]) / ds).floor() as i64).max(0) / TILE;
        let u1 = (((yhi - self.g[2]) / ds).ceil() as i64 - 1).max(0) / TILE;
        for ty in u0..=u1 {
            for tx in t0..=t1 {
                let key = (l, tx, ty);
                if let Some(v) = self.cache.map.read().unwrap().get(&key) {
                    out.push(v.clone());
                    continue;
                }
                let mut v = Vec::new();
                self.origins_cells(ds, tx * TILE, (tx + 1) * TILE, ty * TILE, (ty + 1) * TILE, &mut v);
                let n = v.len();
                let arc = Arc::new(v);
                let mut m = self.cache.map.write().unwrap();
                if cached.load(Ordering::Relaxed) > TILE_CACHE_CAP {
                    // bounded memory: the tiles are rebuilt identically when needed again
                    let freed: usize = m.values().map(|t| t.len()).sum();
                    m.clear();
                    cached.fetch_sub(freed.min(cached.load(Ordering::Relaxed)), Ordering::Relaxed);
                }
                if m.insert(key, arc.clone()).is_none() {
                    cached.fetch_add(n, Ordering::Relaxed);
                }
                out.push(arc);
            }
        }
    }
}

/// E (3) and ∂_b E_a ([a*3 + b]) from the tangential boundary values (Bx, By):
///   E_x = −2 Σ Bx ∂z g,  E_y = −2 Σ By ∂z g,  E_z = 2 Σ (Bx ∂x g + By ∂y g),  g = e^{ikR}/(4πR)
#[inline(always)]
fn add3(k: f64, w: f64, r: [f64; 3], a: J2, e: &mut [Z; 3], de: &mut [Z; 9]) {
    let rr = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    let inv = 1.0 / rr;
    let ex = expi_fast(k * rr).s(w * inv * (0.25 / PI));
    let h = Z::new(-inv, k) * ex; // g'(R)
    let h2 = Z::new(-k * k + 2.0 * inv * inv, -2.0 * k * inv) * ex; // g''(R)
    let hr = h.s(inv);
    let q = (h2 - hr).s(inv * inv);
    let (bx, by) = (a[0], a[1]);
    let bxh = bx * hr;
    let byh = by * hr;
    let bxq = bx * q;
    let byq = by * q;
    e[0] += bxh.s(-2.0 * r[2]);
    e[1] += byh.s(-2.0 * r[2]);
    e[2] += (bxh.s(r[0]) + byh.s(r[1])).s(2.0);
    let sq = bxq.s(r[0]) + byq.s(r[1]);
    for b in 0..3 {
        let mut dx = bxq.s(r[b] * r[2]);
        let mut dy = byq.s(r[b] * r[2]);
        if b == 2 { dx += bxh; dy += byh; }
        de[b] += dx.s(-2.0);
        de[3 + b] += dy.s(-2.0);
        let mut dz = sq.s(r[b]);
        if b == 0 { dz += bxh; }
        if b == 1 { dz += byh; }
        de[6 + b] += dz.s(2.0);
    }
}

/// Tangential E only (level 0, weight 1): the kernel of the lattice convolution.
#[inline(always)]
fn kz(k: f64, rx: f64, ry: f64, dz: f64) -> Z {
    let rr = (rx * rx + ry * ry + dz * dz).sqrt();
    let inv = 1.0 / rr;
    let ex = expi_fast(k * rr).s(inv * (0.25 / PI));
    let h = Z::new(-inv, k) * ex;
    h.s(inv * -2.0 * dz)
}

type Fld = ([Z; 3], [Z; 9]);

/// Level-0 origins as arrays, for the far-field sum (4 lanes; AVX2 + FMA when
/// the CPU has them). Same arithmetic as `add3` with weight 1, per origin;
/// the lanes are added in a fixed order (deterministic on a given machine).
#[derive(Default)]
struct Soa {
    n: usize,
    px: Vec<f64>,
    py: Vec<f64>,
    a: [Vec<f64>; 4],
}

impl Soa {
    fn new(l0: &[Org]) -> Self {
        let n = l0.len().div_ceil(4) * 4;
        let mut s = Soa { n, px: vec![1e6; n], py: vec![1e6; n], a: [vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]] };
        for (i, o) in l0.iter().enumerate() {
            s.px[i] = o.p[0];
            s.py[i] = o.p[1];
            s.a[0][i] = o.a[0].re;
            s.a[1][i] = o.a[0].im;
            s.a[2][i] = o.a[1].re;
            s.a[3][i] = o.a[1].im;
        }
        s
    }

    fn sum_far(&self, a: usize, b: usize, p: [f64; 3], dz: f64, k: f64) -> Fld {
        #[cfg(target_arch = "x86_64")]
        {
            if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma") {
                // SAFETY: the CPU supports the enabled features (checked above)
                return unsafe { self.sum_far_avx2(a, b, p, dz, k) };
            }
        }
        self.sum_far_lanes(a, b, p, dz, k)
    }

    /// The 4-lane sum with AVX2/FMA intrinsics (4 origins per instruction).
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2,fma")]
    unsafe fn sum_far_avx2(&self, a: usize, b: usize, p: [f64; 3], dz: f64, k: f64) -> Fld {
        use std::arch::x86_64::*;
        macro_rules! set { ($v:expr) => { _mm256_set1_pd($v) }; }
        let (px0, py0, dzv, kv) = (set!(p[0]), set!(p[1]), set!(dz), set!(k));
        let dz2 = set!(dz * dz);
        let one = set!(1.0);
        let (two_over_pi, magic) = (set!(0.636_619_772_367_581_343_08), set!(6_755_399_441_055_744.0));
        let (p1, p2, p3) = (set!(1.570_796_326_734_125_614_17e0), set!(6.077_100_506_303_965_976_60e-11), set!(2.022_266_248_795_950_631_54e-21));
        let sc = [-1.666_666_666_666_663_243_48e-1, 8.333_333_333_322_489_461_24e-3, -1.984_126_982_985_794_931_34e-4,
                  2.755_731_370_707_006_767_89e-6, -2.505_076_025_340_686_341_95e-8, 1.589_690_995_211_550_102_21e-10];
        let cc = [4.166_666_666_666_660_190_37e-2, -1.388_888_888_887_410_957_49e-3, 2.480_158_728_947_672_941_78e-5,
                  -2.755_731_435_139_066_330_35e-7, 2.087_572_321_298_174_827_90e-9, -1.135_964_755_778_819_482_65e-11];
        let c4pi = set!(0.25 / PI);
        let (m2, p2v, half) = (set!(-2.0), set!(2.0), set!(0.5));
        let int1 = _mm256_set1_epi64x(1);
        let int2 = _mm256_set1_epi64x(2);
        let zero = _mm256_setzero_pd();
        let mut acc = [zero; 24];
        let mut i = a;
        while i + 4 <= b {
            let rx = _mm256_sub_pd(px0, _mm256_loadu_pd(self.px.as_ptr().add(i)));
            let ry = _mm256_sub_pd(py0, _mm256_loadu_pd(self.py.as_ptr().add(i)));
            let rr2 = _mm256_fmadd_pd(rx, rx, _mm256_fmadd_pd(ry, ry, dz2));
            let rr = _mm256_sqrt_pd(rr2);
            let inv = _mm256_div_pd(one, rr);
            // e^{ikR}: Cody–Waite reduction by π/2, fdlibm kernels, quadrant by bits
            let x = _mm256_mul_pd(kv, rr);
            let qm = _mm256_fmadd_pd(x, two_over_pi, magic);
            let quad = _mm256_castpd_si256(qm);
            let q = _mm256_sub_pd(qm, magic);
            let r = _mm256_fnmadd_pd(q, p3, _mm256_fnmadd_pd(q, p2, _mm256_fnmadd_pd(q, p1, x)));
            let z = _mm256_mul_pd(r, r);
            let mut ps = set!(sc[5]);
            for c in sc[..5].iter().rev() { ps = _mm256_fmadd_pd(ps, z, set!(*c)); }
            let sn = _mm256_fmadd_pd(_mm256_mul_pd(r, z), ps, r);
            let mut pc = set!(cc[5]);
            for c in cc[..5].iter().rev() { pc = _mm256_fmadd_pd(pc, z, set!(*c)); }
            let cs = _mm256_fmadd_pd(_mm256_mul_pd(z, z), pc, _mm256_fnmadd_pd(half, z, one));
            let swap = _mm256_castsi256_pd(_mm256_cmpeq_epi64(_mm256_and_si256(quad, int1), int1));
            let s0 = _mm256_blendv_pd(sn, cs, swap);
            let c0 = _mm256_blendv_pd(cs, sn, swap);
            let sgn_s = _mm256_castsi256_pd(_mm256_slli_epi64(_mm256_and_si256(quad, int2), 62));
            let sgn_c = _mm256_castsi256_pd(_mm256_slli_epi64(_mm256_and_si256(_mm256_add_epi64(quad, int1), int2), 62));
            let m = _mm256_mul_pd(inv, c4pi);
            let exr = _mm256_mul_pd(_mm256_xor_pd(c0, sgn_c), m);
            let exi = _mm256_mul_pd(_mm256_xor_pd(s0, sgn_s), m);
            // h = (−1/R + ik)·ex ; h2 = (−k² + 2/R², −2k/R)·ex ; hr = h/R ; q = (h2 − hr)/R²
            let ninv = _mm256_sub_pd(zero, inv);
            let hre = _mm256_fnmadd_pd(kv, exi, _mm256_mul_pd(ninv, exr));
            let him = _mm256_fmadd_pd(kv, exr, _mm256_mul_pd(ninv, exi));
            let inv2 = _mm256_mul_pd(inv, inv);
            let ar = _mm256_fmadd_pd(p2v, inv2, set!(-k * k));
            let ai = _mm256_mul_pd(set!(-2.0 * k), inv);
            let h2r = _mm256_fnmadd_pd(ai, exi, _mm256_mul_pd(ar, exr));
            let h2i = _mm256_fmadd_pd(ai, exr, _mm256_mul_pd(ar, exi));
            let hrr = _mm256_mul_pd(hre, inv);
            let hri = _mm256_mul_pd(him, inv);
            let qr = _mm256_mul_pd(_mm256_sub_pd(h2r, hrr), inv2);
            let qi = _mm256_mul_pd(_mm256_sub_pd(h2i, hri), inv2);
            let bxr = _mm256_loadu_pd(self.a[0].as_ptr().add(i));
            let bxi = _mm256_loadu_pd(self.a[1].as_ptr().add(i));
            let byr = _mm256_loadu_pd(self.a[2].as_ptr().add(i));
            let byi = _mm256_loadu_pd(self.a[3].as_ptr().add(i));
            // (no closures here: they would not inherit the target features)
            macro_rules! cm {
                ($ar:expr, $ai:expr, $br:expr, $bi:expr) => {
                    (_mm256_fnmadd_pd($ai, $bi, _mm256_mul_pd($ar, $br)), _mm256_fmadd_pd($ai, $br, _mm256_mul_pd($ar, $bi)))
                };
            }
            macro_rules! add {
                ($q:expr, $v:expr) => { acc[$q] = _mm256_add_pd(acc[$q], $v) };
            }
            let (bxhr, bxhi) = cm!(bxr, bxi, hrr, hri);
            let (byhr, byhi) = cm!(byr, byi, hrr, hri);
            let (bxqr, bxqi) = cm!(bxr, bxi, qr, qi);
            let (byqr, byqi) = cm!(byr, byi, qr, qi);
            let m2dz = _mm256_mul_pd(m2, dzv);
            add!(0, _mm256_mul_pd(bxhr, m2dz));
            add!(1, _mm256_mul_pd(bxhi, m2dz));
            add!(2, _mm256_mul_pd(byhr, m2dz));
            add!(3, _mm256_mul_pd(byhi, m2dz));
            add!(4, _mm256_mul_pd(_mm256_fmadd_pd(bxhr, rx, _mm256_mul_pd(byhr, ry)), p2v));
            add!(5, _mm256_mul_pd(_mm256_fmadd_pd(bxhi, rx, _mm256_mul_pd(byhi, ry)), p2v));
            let sqr = _mm256_fmadd_pd(bxqr, rx, _mm256_mul_pd(byqr, ry));
            let sqi = _mm256_fmadd_pd(bxqi, rx, _mm256_mul_pd(byqi, ry));
            let r3 = [rx, ry, dzv];
            for bb in 0..3 {
                let rbdz = _mm256_mul_pd(r3[bb], dzv);
                let (mut dxr, mut dxi) = (_mm256_mul_pd(bxqr, rbdz), _mm256_mul_pd(bxqi, rbdz));
                let (mut dyr, mut dyi) = (_mm256_mul_pd(byqr, rbdz), _mm256_mul_pd(byqi, rbdz));
                if bb == 2 {
                    dxr = _mm256_add_pd(dxr, bxhr);
                    dxi = _mm256_add_pd(dxi, bxhi);
                    dyr = _mm256_add_pd(dyr, byhr);
                    dyi = _mm256_add_pd(dyi, byhi);
                }
                let (mut dzr, mut dzi) = (_mm256_mul_pd(sqr, r3[bb]), _mm256_mul_pd(sqi, r3[bb]));
                if bb == 0 { dzr = _mm256_add_pd(dzr, bxhr); dzi = _mm256_add_pd(dzi, bxhi); }
                if bb == 1 { dzr = _mm256_add_pd(dzr, byhr); dzi = _mm256_add_pd(dzi, byhi); }
                add!(6 + 2 * bb, _mm256_mul_pd(dxr, m2));
                add!(7 + 2 * bb, _mm256_mul_pd(dxi, m2));
                add!(12 + 2 * bb, _mm256_mul_pd(dyr, m2));
                add!(13 + 2 * bb, _mm256_mul_pd(dyi, m2));
                add!(18 + 2 * bb, _mm256_mul_pd(dzr, p2v));
                add!(19 + 2 * bb, _mm256_mul_pd(dzi, p2v));
            }
            i += 4;
        }
        let mut lanes = [[0.0f64; 4]; 24];
        for (q, l) in lanes.iter_mut().enumerate() { _mm256_storeu_pd(l.as_mut_ptr(), acc[q]); }
        let tot = |q: usize| lanes[q][0] + lanes[q][1] + lanes[q][2] + lanes[q][3];
        let mut f: Fld = ([Z::ZERO; 3], [Z::ZERO; 9]);
        for c in 0..3 { f.0[c] = Z::new(tot(2 * c), tot(2 * c + 1)); }
        for c in 0..9 { f.1[c] = Z::new(tot(6 + 2 * c), tot(7 + 2 * c)); }
        f
    }

    #[inline(always)]
    fn sum_far_lanes(&self, a: usize, b: usize, p: [f64; 3], dz: f64, k: f64) -> Fld {
        const L: usize = 4;
        const TWO_OVER_PI: f64 = 0.636_619_772_367_581_343_08;
        const P1: f64 = 1.570_796_326_734_125_614_17e0;
        const P2: f64 = 6.077_100_506_303_965_976_60e-11;
        const P3: f64 = 2.022_266_248_795_950_631_54e-21;
        const MAGIC: f64 = 6_755_399_441_055_744.0;
        let c4pi = 0.25 / PI;
        // accumulators: E (re, im) x3, dE (re, im) x9, per lane
        let mut acc = [[0.0f64; L]; 24];
        let mut i = a;
        while i + L <= b {
            let px = &self.px[i..i + L];
            let py = &self.py[i..i + L];
            let (axr, axi, ayr, ayi) = (&self.a[0][i..i + L], &self.a[1][i..i + L], &self.a[2][i..i + L], &self.a[3][i..i + L]);
            let mut v = [[0.0f64; L]; 24];
            for l in 0..L {
                let rx = p[0] - px[l];
                let ry = p[1] - py[l];
                let rr = (rx * rx + ry * ry + dz * dz).sqrt();
                let inv = 1.0 / rr;
                // e^{ikR} (same reduction and kernels as `sincos`, branch-free)
                let x = k * rr;
                let qm = x * TWO_OVER_PI + MAGIC;
                let quad = qm.to_bits() & 3;
                let q = qm - MAGIC;
                let r = ((x - q * P1) - q * P2) - q * P3;
                let z = r * r;
                let sn = r + r * z * (-1.666_666_666_666_663_243_48e-1 + z * (8.333_333_333_322_489_461_24e-3
                    + z * (-1.984_126_982_985_794_931_34e-4 + z * (2.755_731_370_707_006_767_89e-6
                    + z * (-2.505_076_025_340_686_341_95e-8 + z * 1.589_690_995_211_550_102_21e-10)))));
                let cs = 1.0 - 0.5 * z + z * z * (4.166_666_666_666_660_190_37e-2 + z * (-1.388_888_888_887_410_957_49e-3
                    + z * (2.480_158_728_947_672_941_78e-5 + z * (-2.755_731_435_139_066_330_35e-7
                    + z * (2.087_572_321_298_174_827_90e-9 + z * -1.135_964_755_778_819_482_65e-11)))));
                let swap = quad & 1 == 1;
                let (s0, c0) = if swap { (cs, sn) } else { (sn, cs) };
                let ssg = if quad & 2 == 2 { -1.0 } else { 1.0 };
                let csg = if (quad + 1) & 2 == 2 { -1.0 } else { 1.0 };
                let m = inv * c4pi;
                let (exr, exi) = (c0 * csg * m, s0 * ssg * m);
                // h = (-1/R + ik) ex ; h2 = (-k^2 + 2/R^2, -2k/R) ex
                let (hr_, hi_) = (-inv * exr - k * exi, -inv * exi + k * exr);
                let (ar, ai) = (-k * k + 2.0 * inv * inv, -2.0 * k * inv);
                let (h2r, h2i) = (ar * exr - ai * exi, ar * exi + ai * exr);
                let (hrr, hri) = (hr_ * inv, hi_ * inv);
                let inv2 = inv * inv;
                let (qr, qi) = ((h2r - hrr) * inv2, (h2i - hri) * inv2);
                let (bxr, bxi, byr, byi) = (axr[l], axi[l], ayr[l], ayi[l]);
                let (bxhr, bxhi) = (bxr * hrr - bxi * hri, bxr * hri + bxi * hrr);
                let (byhr, byhi) = (byr * hrr - byi * hri, byr * hri + byi * hrr);
                let (bxqr, bxqi) = (bxr * qr - bxi * qi, bxr * qi + bxi * qr);
                let (byqr, byqi) = (byr * qr - byi * qi, byr * qi + byi * qr);
                let r3 = [rx, ry, dz];
                v[0][l] = bxhr * (-2.0 * dz);
                v[1][l] = bxhi * (-2.0 * dz);
                v[2][l] = byhr * (-2.0 * dz);
                v[3][l] = byhi * (-2.0 * dz);
                v[4][l] = (bxhr * rx + byhr * ry) * 2.0;
                v[5][l] = (bxhi * rx + byhi * ry) * 2.0;
                let (sqr, sqi) = (bxqr * rx + byqr * ry, bxqi * rx + byqi * ry);
                for bb in 0..3 {
                    let rb = r3[bb];
                    let (mut dxr, mut dxi) = (bxqr * (rb * dz), bxqi * (rb * dz));
                    let (mut dyr, mut dyi) = (byqr * (rb * dz), byqi * (rb * dz));
                    if bb == 2 { dxr += bxhr; dxi += bxhi; dyr += byhr; dyi += byhi; }
                    let (mut dzr, mut dzi) = (sqr * rb, sqi * rb);
                    if bb == 0 { dzr += bxhr; dzi += bxhi; }
                    if bb == 1 { dzr += byhr; dzi += byhi; }
                    v[6 + 2 * bb][l] = dxr * -2.0;
                    v[7 + 2 * bb][l] = dxi * -2.0;
                    v[12 + 2 * bb][l] = dyr * -2.0;
                    v[13 + 2 * bb][l] = dyi * -2.0;
                    v[18 + 2 * bb][l] = dzr * 2.0;
                    v[19 + 2 * bb][l] = dzi * 2.0;
                }
            }
            for (q, vq) in v.iter().enumerate() {
                for l in 0..L { acc[q][l] += vq[l]; }
            }
            i += L;
        }
        let tot = |q: usize| acc[q][0] + acc[q][1] + acc[q][2] + acc[q][3];
        let mut f: Fld = ([Z::ZERO; 3], [Z::ZERO; 9]);
        for c in 0..3 { f.0[c] = Z::new(tot(2 * c), tot(2 * c + 1)); }
        for c in 0..9 { f.1[c] = Z::new(tot(6 + 2 * c), tot(7 + 2 * c)); }
        f
    }
}

#[inline]
fn fld_add(a: &mut Fld, b: &Fld) {
    for i in 0..3 { a.0[i] += b.0[i]; }
    for i in 0..9 { a.1[i] += b.1[i]; }
}

// ─── FFT (radix 2) for the lattice convolutions ───────────────────────────────

fn fft(buf: &mut [Z], tw: &[Z]) {
    let n = buf.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 { j ^= bit; bit >>= 1; }
        j |= bit;
        if i < j { buf.swap(i, j); }
    }
    let mut len = 2;
    while len <= n {
        let stride = n / len;
        for s in (0..n).step_by(len) {
            for q in 0..len / 2 {
                let w = tw[q * stride];
                let u = buf[s + q];
                let v = buf[s + q + len / 2] * w;
                buf[s + q] = u + v;
                buf[s + q + len / 2] = u - v;
            }
        }
        len <<= 1;
    }
}

fn twiddles(n: usize, inverse: bool) -> Vec<Z> {
    let sg = if inverse { 1.0 } else { -1.0 };
    (0..n / 2).map(|j| Z::expi(sg * 2.0 * PI * j as f64 / n as f64)).collect()
}

fn transpose(a: &[Z], rows: usize, cols: usize) -> Vec<Z> {
    let mut out = vec![Z::ZERO; a.len()];
    out.par_chunks_mut(rows).enumerate().for_each(|(c, col)| {
        for r in 0..rows { col[r] = a[r * cols + c]; }
    });
    out
}

/// 2-D FFT of a rows×cols array (row-major, both powers of two).
fn fft2(a: Vec<Z>, rows: usize, cols: usize, inverse: bool) -> Vec<Z> {
    let mut a = a;
    let twc = twiddles(cols, inverse);
    a.par_chunks_mut(cols).for_each(|r| fft(r, &twc));
    let mut t = transpose(&a, rows, cols);
    drop(a);
    let twr = twiddles(rows, inverse);
    t.par_chunks_mut(rows).for_each(|r| fft(r, &twr));
    transpose(&t, cols, rows)
}

// ─── The bench: objects and fans of one laboratory ────────────────────────────

/// Birth table of the photon ball on the source: |B|² on the λ/16 cells.
struct Birth {
    x0: f64,
    y0: f64,
    nx: usize,
    cdf: Vec<f64>,
}

/// Ball count of one fan (docs/REGOLA_BASE.md §2): origins (level 0),
/// directions per origin, balls of level 0 and of all levels.
#[derive(Clone, Debug)]
pub struct BallCount {
    pub obj: String,
    pub n0: usize,
    pub ndir: usize,
    pub balls0: f64,
    pub balls_all: f64,
}

/// The wave of one laboratory: identical for every photon, built once.
pub struct Bench3 {
    pub key: String,
    pub lam: f64,
    pub lev: Lev,
    pub rule: RuleV2,
    pub objs: Vec<Obj>,
    pub fans: Vec<Fan>,
    pub fan_of: Vec<Option<usize>>,
    pub src: usize,
    pub screen_z: f64,
    /// Photons leaving this box are lost: [x0, x1, y0, y1, z0, z1].
    pub bbox: [f64; 6],
    /// Fraction of the source beam's |B|² on the opaque part of the source mask.
    pub wall_frac: f64,
    pub build_secs: f64,
    birth: Birth,
    cached: AtomicUsize,
    flux: std::sync::Mutex<Option<((u64, usize), Vec<f32>)>>,
}

fn arange(a: f64, b: f64, step: f64) -> Vec<f64> {
    let n = ((b - a) / step).ceil().max(0.0) as usize;
    (0..n).map(|i| a + i as f64 * step).collect()
}

impl Bench3 {
    /// Build the fans of the objects, in the order of the light.
    /// `objs[src]` must be the source; objects are planes ⊥ z sorted by z.
    pub fn build(key: String, lam: f64, objs: Vec<Obj>, bbox: [f64; 6]) -> Self {
        let t0 = std::time::Instant::now();
        let rule = RULE;
        let lev = Lev::new(lam, &rule);
        let src = objs.iter().position(|o| o.kind == Kind::Source).expect("a source");
        let screen_z = objs.iter().find(|o| o.kind == Kind::Screen).map_or(f64::INFINITY, |o| o.z);
        let mut b = Bench3 {
            key, lam, lev, rule, objs: Vec::new(), fans: Vec::new(), fan_of: Vec::new(), src, screen_z, bbox,
            wall_frac: f64::NAN, build_secs: 0.0, birth: Birth { x0: 0.0, y0: 0.0, nx: 0, cdf: Vec::new() },
            cached: AtomicUsize::new(0), flux: std::sync::Mutex::new(None),
        };
        for o in objs { b.add(o); }
        // drawing: emission times and kill planes
        let zsrc = b.objs[b.src].z;
        let nf = b.fans.len();
        for m in 0..nf {
            let z_kill = if m + 1 < nf { b.fans[m + 1].z0 } else { screen_z };
            b.fans[m].t_emit = b.fans[m].z0 - zsrc;
            b.fans[m].z_kill = z_kill;
        }
        b.make_birth();
        b.build_secs = t0.elapsed().as_secs_f64();
        b
    }

    fn add(&mut self, o: Obj) {
        let i = self.objs.len();
        self.objs.push(o.clone());
        self.fan_of.push(None);
        if matches!(o.kind, Kind::Screen | Kind::Plane) { return; }
        let lam = self.lam;
        let ds0 = self.lev.ds0;
        let sds = ds0 / 4.0;
        let eps = self.rule.eps;
        let prev = if o.kind == Kind::Source { None } else { self.fans.len().checked_sub(1) };
        if o.kind != Kind::Source && prev.is_none() { return; }
        let mut bx = match prev {
            None => {
                let s = o.src.as_ref().expect("source beam");
                let ext = |w0: f64| {
                    if !w0.is_finite() { return 50.0 * lam; }
                    let zr = PI * w0 * w0 / lam;
                    4.0 * w0 * (1.0 + (s.zw / zr).powi(2)).sqrt()
                };
                let (ex, ey) = (ext(s.w0), ext(s.w0y.unwrap_or(s.w0)));
                [-ex, ex, -ey, ey]
            }
            Some(p) => {
                let pf = &self.fans[p];
                let m = 0.5 * (o.z - pf.z0) + 2.0 * lam;
                [pf.g[0] - m, pf.g[1] + m, pf.g[2] - m, pf.g[3] + m]
            }
        };
        let hull = o.mask.hull();
        let clip = |b: [f64; 4]| match hull {
            Some(h) => [b[0].max(h[0]), b[1].min(h[1]), b[2].max(h[2]), b[3].min(h[3])],
            None => b,
        };
        bx = clip(bx);
        // B on a lattice (x0, y0, step, nx, ny)
        let bvals = |me: &Self, x0: f64, y0: f64, step: f64, nx: usize, ny: usize| -> Vec<J2> {
            match prev {
                None => (0..nx * ny).into_par_iter()
                    .map(|q| o.src_b(lam, x0 + (q % nx) as f64 * step, y0 + (q / nx) as f64 * step)).collect(),
                Some(p) => {
                    let e = me.lattice_tangential(p, o.z, x0, y0, step, nx, ny);
                    e.into_par_iter().enumerate()
                        .map(|(q, v)| o.transmit(x0 + (q % nx) as f64 * step, y0 + (q / nx) as f64 * step, v)).collect()
                }
            }
        };
        // where |B| >= ε·max (grow the box while the light reaches its border)
        let (mut xs, mut ys, mut on) = (Vec::new(), Vec::new(), Vec::new());
        for _ in 0..8 {
            xs = arange(bx[0] + ds0 / 2.0, bx[1], ds0);
            ys = arange(bx[2] + ds0 / 2.0, bx[3], ds0);
            if xs.is_empty() || ys.is_empty() { return; }
            let (nx, ny) = (xs.len(), ys.len());
            let bv = bvals(self, xs[0], ys[0], ds0, nx, ny);
            let mag: Vec<f64> = bv.iter().enumerate()
                .map(|(q, v)| if o.mask.inside(xs[q % nx], ys[q / nx]) { j2n(*v).sqrt() } else { 0.0 }).collect();
            let mx = mag.iter().cloned().fold(0.0, f64::max);
            if mx <= 1e-300 { return; }
            on = mag.iter().map(|&m| m >= eps * mx).collect::<Vec<bool>>();
            let edge = (0..nx).any(|c| on[c] || on[(ny - 1) * nx + c]) || (0..ny).any(|r| on[r * nx] || on[r * nx + nx - 1]);
            let at_hull = hull.is_some_and(|h| (0..4).all(|q| (bx[q] - h[q]).abs() < 1e-9));
            if !edge || at_hull { break; }
            let (cx, cy) = (0.5 * (bx[0] + bx[1]), 0.5 * (bx[2] + bx[3]));
            let (hx, hy) = (0.75 * (bx[1] - bx[0]), 0.75 * (bx[3] - bx[2]));
            bx = clip([cx - hx, cx + hx, cy - hy, cy + hy]);
        }
        let nxs = xs.len();
        let (mut ix0, mut ix1, mut iy0, mut iy1) = (usize::MAX, 0, usize::MAX, 0);
        for (q, &f) in on.iter().enumerate() {
            if f {
                let (c, r) = (q % nxs, q / nxs);
                ix0 = ix0.min(c); ix1 = ix1.max(c); iy0 = iy0.min(r); iy1 = iy1.max(r);
            }
        }
        let g = [bx[0].max(xs[ix0] - ds0), bx[1].min(xs[ix1] + ds0), bx[2].max(ys[iy0] - ds0), bx[3].min(ys[iy1] + ds0)];
        let sx0 = g[0] - 3.0 * sds;
        let sy0 = g[2] - 3.0 * sds;
        let nx = ((g[1] - g[0]) / sds).ceil() as usize + 7;
        let ny = ((g[3] - g[2]) / sds).ceil() as usize + 7;
        let samp = bvals(self, sx0, sy0, sds, nx, ny);
        let mut bmax = 0.0f64;
        for r in 0..ny {
            for c in 0..nx {
                if o.mask.inside(sx0 + c as f64 * sds, sy0 + r as f64 * sds) {
                    bmax = bmax.max(j2n(samp[r * nx + c]).sqrt());
                }
            }
        }
        let mut fan = Fan {
            obj: i, z0: o.z, g, mask: o.mask.clone(), sx0, sy0, sds, nx, ny, samp, bthr: eps * bmax,
            l0: Vec::new(), soa: Soa::default(), cache: TileCache::default(), w: [0.0; 2], ndir: 0, t_emit: 0.0, z_kill: 0.0,
        };
        let mut l0 = Vec::new();
        let ni = ((g[1] - g[0]) / ds0).ceil() as i64;
        let nj = ((g[3] - g[2]) / ds0).ceil() as i64;
        fan.origins_cells(ds0, 0, ni, 0, nj, &mut l0);
        if !l0.is_empty() {
            let (mut a, mut bb, mut c, mut d) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
            for q in &l0 { a = a.min(q.p[0]); bb = bb.max(q.p[0]); c = c.min(q.p[1]); d = d.max(q.p[1]); }
            fan.w = [bb - a + ds0, d - c + ds0];
            let eta = self.rule.eta;
            fan.ndir = (PI / 4.0 * (2.0 * eta * fan.w[0] / lam) * (2.0 * eta * fan.w[1] / lam)).ceil() as usize;
        }
        fan.soa = Soa::new(&l0);
        fan.l0 = l0;
        self.fan_of[i] = Some(self.fans.len());
        self.fans.push(fan);
    }

    /// E_x, E_y of fan `fi` on the lattice (x0 + i·step, y0 + j·step, z).
    /// Planes ≥ 2·R0 apart: only level 0 with weight 1 → discrete convolution
    /// of the origins (a lattice of spacing ds0 = m·step) with the kernel,
    /// done by FFT; origins off the lattice (cells cut by an edge) and close
    /// planes are summed directly.
    pub fn lattice_tangential(&self, fi: usize, z: f64, x0: f64, y0: f64, step: f64, nx: usize, ny: usize) -> Vec<J2> {
        let fan = &self.fans[fi];
        let dz = z - fan.z0;
        let k = self.lev.k;
        let ds0 = self.lev.ds0;
        let direct = |pts: &[Org], out: &mut Vec<J2>| {
            out.par_iter_mut().enumerate().for_each(|(q, v)| {
                let (x, y) = (x0 + (q % nx) as f64 * step, y0 + (q / nx) as f64 * step);
                for o in pts {
                    let kk = kz(k, x - o.p[0], y - o.p[1], dz);
                    v[0] += kk * o.a[0];
                    v[1] += kk * o.a[1];
                }
            });
        };
        let mut out = vec![[Z::ZERO; 2]; nx * ny];
        let m = (ds0 / step).round();
        if dz < 2.0 * self.lev.r0 + 1e-9 || m < 1.0 || (m * step - ds0).abs() > 1e-9 * ds0 {
            out.par_iter_mut().enumerate().for_each(|(q, v)| {
                let (e, _) = self.field([x0 + (q % nx) as f64 * step, y0 + (q / nx) as f64 * step, z], fi);
                *v = [e[0], e[1]];
            });
            return out;
        }
        let m = m as usize;
        let ox0 = fan.g[0] + 0.5 * ds0;
        let oy0 = fan.g[2] + 0.5 * ds0;
        let mut on: Vec<(usize, usize, J2)> = Vec::new();
        let mut off: Vec<Org> = Vec::new();
        for o in &fan.l0 {
            let fi_ = (o.p[0] - ox0) / ds0;
            let fj_ = (o.p[1] - oy0) / ds0;
            if (fi_ - fi_.round()).abs() < 1e-7 && (fj_ - fj_.round()).abs() < 1e-7 && fi_ > -0.5 && fj_ > -0.5 {
                on.push((fi_.round() as usize, fj_.round() as usize, o.a));
            } else {
                off.push(*o);
            }
        }
        if !on.is_empty() {
            let nox = on.iter().map(|q| q.0).max().unwrap() + 1;
            let noy = on.iter().map(|q| q.1).max().unwrap() + 1;
            let (tax, tay) = (m * (nox - 1) + 1, m * (noy - 1) + 1);
            let (lx, ly) = (tax - 1 + nx, tay - 1 + ny);
            let (cx, cy) = (lx.next_power_of_two(), ly.next_power_of_two());
            let mut kern = vec![Z::ZERO; cx * cy];
            let (dx0, dy0) = (x0 - ox0, y0 - oy0);
            kern.par_chunks_mut(cx).enumerate().for_each(|(uy, row)| {
                if uy >= ly { return; }
                let ry = dy0 + (uy as f64 - (tay - 1) as f64) * step;
                for (ux, v) in row.iter_mut().enumerate().take(lx) {
                    let rx = dx0 + (ux as f64 - (tax - 1) as f64) * step;
                    *v = kz(k, rx, ry, dz);
                }
            });
            let kf = fft2(kern, cy, cx, false);
            let norm = 1.0 / (cx * cy) as f64;
            for comp in 0..2 {
                let mut a = vec![Z::ZERO; cx * cy];
                for &(i, j, amp) in &on { a[(m * j) * cx + m * i] = amp[comp]; }
                let af = fft2(a, cy, cx, false);
                let prod: Vec<Z> = af.par_iter().zip(kf.par_iter()).map(|(p, q)| *p * *q).collect();
                drop(af);
                let r = fft2(prod, cy, cx, true);
                out.par_chunks_mut(nx).enumerate().for_each(|(qy, row)| {
                    let src = &r[(tay - 1 + qy) * cx + tax - 1..];
                    for (qx, v) in row.iter_mut().enumerate() { v[comp] = src[qx].s(norm); }
                });
            }
        }
        if !off.is_empty() { direct(&off, &mut out); }
        out
    }

    /// Birth table on the source fan: |B|² on the λ/16 cells touching the opening.
    fn make_birth(&mut self) {
        let Some(fi) = self.fan_of[self.src] else { return };
        let o = &self.objs[self.src];
        let lam = self.lam;
        let sds = self.lev.ds0 / 4.0;
        let g = self.fans[fi].g;
        let xs = arange(g[0] + sds / 2.0, g[1], sds);
        let ys = arange(g[2] + sds / 2.0, g[3], sds);
        let nx = xs.len();
        let p: Vec<f64> = (0..nx * ys.len()).into_par_iter().map(|q| {
            let (x, y) = (xs[q % nx], ys[q / nx]);
            let mut touch = false;
            for ddx in [-0.5, 0.0, 0.5] {
                for ddy in [-0.5, 0.0, 0.5] {
                    touch |= o.mask.inside(x + ddx * sds, y + ddy * sds);
                }
            }
            if touch { j2n(o.src_b(lam, x, y)) } else { 0.0 }
        }).collect();
        let mut acc = 0.0;
        let mut cdf: Vec<f64> = p.iter().map(|v| { acc += v; acc }).collect();
        let tot = acc.max(1e-300);
        cdf.iter_mut().for_each(|v| *v /= tot);
        self.birth = Birth { x0: xs[0], y0: ys[0], nx, cdf };
        // fraction of the beam on the opaque part (over the whole beam)
        let s = o.src.as_ref().unwrap();
        if s.w0.is_finite() {
            let ext = |w0: f64| { let zr = PI * w0 * w0 / lam; 4.0 * w0 * (1.0 + (s.zw / zr).powi(2)).sqrt() };
            let (bx, by) = (ext(s.w0), ext(s.w0y.unwrap_or(s.w0)));
            let (hx, hy) = (sds.min(bx / 200.0), sds.min(by / 200.0));
            let xx = arange(-bx + hx / 2.0, bx, hx);
            let yy = arange(-by + hy / 2.0, by, hy);
            let (tot, ins) = yy.par_iter().map(|&y| {
                let (mut t, mut i) = (0.0, 0.0);
                for &x in &xx {
                    let v = j2n(o.src_b(lam, x, y));
                    t += v;
                    if o.mask.inside(x, y) { i += v; }
                }
                (t, i)
            }).reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
            self.wall_frac = 1.0 - ins / tot;
        }
    }

    /// E and ∂_b E_a of fan `fi` at p (what a photon ball hearing that fan reads).
    pub fn field(&self, p: [f64; 3], fi: usize) -> Fld {
        let mut acc: Fld = ([Z::ZERO; 3], [Z::ZERO; 9]);
        let Some(fan) = self.fans.get(fi) else { return acc };
        let lv = self.lev;
        let dz = p[2] - fan.z0;
        if dz <= 1e-12 { return acc; }
        let two_r0 = 2.0 * lv.r0;
        if dz >= two_r0 {
            // every level-0 origin is at R >= 2 R0: weight 1, no fine level (rmax(1) = 2 R0)
            let n = fan.soa.n;
            let nch = n.div_ceil(CHUNK);
            let part = |c: usize| fan.soa.sum_far(c * CHUNK, ((c + 1) * CHUNK).min(n), p, dz, lv.k);
            if fan.l0.len() >= PAR_MIN {
                let parts: Vec<Fld> = (0..nch).into_par_iter().map(part).collect();
                for q in &parts { fld_add(&mut acc, q); }
            } else {
                for c in 0..nch { fld_add(&mut acc, &part(c)); }
            }
            return acc;
        }
        let part = |c: &[Org]| -> Fld {
            let mut f: Fld = ([Z::ZERO; 3], [Z::ZERO; 9]);
            for o in c {
                let r = [p[0] - o.p[0], p[1] - o.p[1], dz];
                let rr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                let w = if rr2 >= two_r0 * two_r0 { 1.0 } else { lv.w(0, rr2.sqrt()) };
                if w == 0.0 { continue; }
                add3(lv.k, w, r, o.a, &mut f.0, &mut f.1);
            }
            f
        };
        if fan.l0.len() >= PAR_MIN {
            let parts: Vec<Fld> = fan.l0.par_chunks(CHUNK).map(part).collect();
            for q in &parts { fld_add(&mut acc, q); }
        } else {
            for c in fan.l0.chunks(CHUNK) { fld_add(&mut acc, &part(c)); }
        }
        let mut scratch: Vec<Arc<Vec<Org>>> = Vec::new();
        for l in 1..=lv.lmax {
            let rm = lv.rmax(l);
            if dz >= rm { continue; }
            let half = (rm * rm - dz * dz).sqrt();
            scratch.clear();
            fan.tiles(&self.cached, l, lv.ds(l), p[0] - half, p[0] + half, p[1] - half, p[1] + half, &mut scratch);
            let rm2 = rm * rm;
            let mut f: Fld = ([Z::ZERO; 3], [Z::ZERO; 9]);
            for t in scratch.iter() {
                for o in t.iter() {
                    let r = [p[0] - o.p[0], p[1] - o.p[1], dz];
                    let rr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    if rr2 >= rm2 { continue; }
                    let w = lv.w(l, rr2.sqrt());
                    if w == 0.0 { continue; }
                    add3(lv.k, w, r, o.a, &mut f.0, &mut f.1);
                }
            }
            fld_add(&mut acc, &f);
        }
        acc
    }

    /// Current j = Σ_a Im(E_a* ∇E_a) and |E|² of fan `fi` at p (None: no wave).
    pub fn current(&self, p: [f64; 3], fi: Option<usize>) -> ([f64; 3], f64) {
        let Some(fi) = fi else { return ([0.0; 3], 0.0) };
        let (e, de) = self.field(p, fi);
        let mut j = [0.0; 3];
        for (b, jb) in j.iter_mut().enumerate() {
            *jb = (0..3).map(|a| (e[a].conj() * de[a * 3 + b]).im).sum();
        }
        (j, e.iter().map(|v| v.norm2()).sum())
    }

    /// Wave balls the rule uses, per fan.
    pub fn ball_counts(&self) -> Vec<BallCount> {
        let lmax = self.lev.lmax as i32;
        self.fans.iter().filter(|f| !f.l0.is_empty()).map(|f| {
            let n0 = f.l0.len();
            let all: f64 = (0..=lmax).map(|l| n0 as f64 * 4f64.powi(l)).sum();
            BallCount {
                obj: self.objs[f.obj].name.clone(), n0, ndir: f.ndir,
                balls0: n0 as f64 * f.ndir as f64, balls_all: all * f.ndir as f64,
            }
        }).collect()
    }

    /// Origins kept in the tile caches (fine levels), all fans.
    pub fn cached_origins(&self) -> usize { self.cached.load(Ordering::Relaxed) }

    /// `screen_flux_x` over |y| ≤ 30 λ, computed once per (half, bins).
    pub fn screen_flux_x_cached(&self, half: f64, bins: usize) -> Vec<f32> {
        let key = (half.to_bits(), bins);
        let mut g = self.flux.lock().unwrap();
        if let Some((k, v)) = g.as_ref() {
            if *k == key { return v.clone(); }
        }
        let v = self.screen_flux_x(half, bins, 30.0 * self.lam, 60);
        *g = Some((key, v.clone()));
        v
    }

    /// x-marginal of the flux j_z on the screen (bins over ±half, y over ±yh):
    /// what the clicks must build up. From the last fan (the program's own wave).
    pub fn screen_flux_x(&self, half: f64, bins: usize, yh: f64, ny: usize) -> Vec<f32> {
        let Some(fi) = self.fans.len().checked_sub(1) else { return vec![0.0; bins] };
        let z = self.screen_z;
        let sub = 2;
        let n = bins * sub;
        let w = 2.0 * half / n as f64;
        let fine: Vec<f64> = (0..n).into_par_iter().map(|i| {
            let x = -half + (i as f64 + 0.5) * w;
            (0..ny).map(|q| {
                let y = -yh + (q as f64 + 0.5) * 2.0 * yh / ny as f64;
                self.current([x, y, z], Some(fi)).0[2].max(0.0)
            }).sum::<f64>()
        }).collect();
        let coarse: Vec<f64> = fine.chunks(sub).map(|c| c.iter().sum()).collect();
        let s: f64 = coarse.iter().sum();
        coarse.iter().map(|v| if s > 0.0 { (v / s) as f32 } else { 0.0 }).collect()
    }
}

// ─── The photon ball ──────────────────────────────────────────────────────────

const DP_A: [[f64; 6]; 7] = [
    [0.0; 6],
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0],
    [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0],
    [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0],
];
const DP_B5: [f64; 7] = [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0];
const DP_B4: [f64; 7] = [5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0];

/// How a photon ball ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    /// Clicked on the screen.
    Click,
    /// Hit the opaque part of object i (absorbed whole).
    Absorbed(usize),
    /// Removed by the local flux coin at object i.
    Coin(usize),
    /// Removed by the flux coin at the source.
    SourceCoin,
    /// Left the laboratory.
    Lost,
    /// Step limit.
    Stuck,
}

#[inline] fn add_s(a: [f64; 3], s: f64, b: [f64; 3]) -> [f64; 3] { [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]] }
#[inline] fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
#[inline] fn norm3(a: [f64; 3]) -> f64 { (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt() }

/// An edge zone crossed by the photon ball (v2): distance u from the edge
/// into the opening, edge point, unit vector into the opening, zone
/// half-width, curvature of a circular edge (area weight).
struct EdgeHit { u: f64, e: [f64; 3], dv: [f64; 3], w: f64, curv: f64 }

#[derive(Clone)]
pub struct Ball3 {
    pub bench: Arc<Bench3>,
    pub x: [f64; 3],
    pub d: [f64; 3],
    /// Light-clock of the photon ball (path length at c).
    pub t: f64,
    h: f64,
    k1: Option<[f64; 3]>,
    /// Fan it hears (the last object crossed).
    pub hear: Option<usize>,
    pub steps: usize,
    pub end: Option<End>,
    pub born: [f64; 3],
    /// Last local coin: (p, survived).
    pub last_coin: Option<(f64, bool)>,
    /// Recent path: (light-clock, position).
    pub trail: VecDeque<(f64, [f64; 3])>,
    /// Crossings of passive planes: (object index, position).
    pub rec: Vec<(usize, [f64; 3])>,
    key: u64,
    draws: u64,
}

impl Ball3 {
    fn rnd(&mut self) -> f64 {
        let u = uniform(self.key, self.draws);
        self.draws += 1;
        u
    }

    /// A photon ball born on the source with probability |B|², flown δ along
    /// the normal, then the source coin (Φ_out/Φ_in, edge zone v2 if near an edge).
    pub fn new(bench: Arc<Bench3>, seed: u64, id: u64) -> Self {
        let b = bench.clone();
        let src = &b.objs[b.src];
        let delta = b.rule.delta * b.lam;
        let sds = b.lev.ds0 / 4.0;
        let fan = b.fan_of[b.src];
        let mut me = Ball3 {
            bench, x: [0.0, 0.0, src.z], d: [0.0, 0.0, 1.0], t: delta, h: delta / 4.0, k1: None, hear: fan, steps: 0,
            end: None, born: [0.0; 3], last_coin: None, trail: VecDeque::new(), rec: Vec::new(),
            key: hash2(seed ^ 0xBA5E_3D00_0000_0001, id), draws: 0,
        };
        if b.birth.cdf.is_empty() {
            me.end = Some(End::SourceCoin);
            return me;
        }
        // rejection: a point of a cell touching the opening, kept if inside it
        let (mut x, mut y) = (0.0, 0.0);
        for _ in 0..100_000 {
            let u = me.rnd();
            let k = b.birth.cdf.partition_point(|&c| c < u).min(b.birth.cdf.len() - 1);
            x = b.birth.x0 + (k % b.birth.nx) as f64 * sds + (me.rnd() - 0.5) * sds;
            y = b.birth.y0 + (k / b.birth.nx) as f64 * sds + (me.rnd() - 0.5) * sds;
            if src.mask.inside(x, y) { break; }
        }
        let x0 = [x, y, src.z];
        me.born = x0;
        let lam = b.lam;
        let phi_in = b.lev.k * j2n(src.src_b(lam, x, y));
        let mut x1 = [x, y, src.z + delta];
        let (j, _) = b.current(x1, fan);
        let mut pr = (j[2] / phi_in.max(1e-300)).clamp(0.0, 1.0);
        if let Some(eh) = me.edge(src, x0) {
            let kk = b.lev.k;
            let (np, p) = me.transport(&eh, |q: [f64; 3]| kk * j2n(src.src_b(lam, q[0], q[1])), |q| b.current(q, fan).0);
            x1 = np;
            pr = p;
        }
        let ok = me.rnd() < pr;
        me.last_coin = Some((pr, ok));
        me.x = x1;
        me.push_trail();
        if !ok {
            me.end = Some(End::SourceCoin);
            me.t = 0.0;
        }
        me
    }

    pub fn pos(&self) -> Vec3 { Vec3::new(self.x[0] as f32, self.x[1] as f32, self.x[2] as f32) }

    /// Position at light-clock `t` (≤ the ball's own clock), from the trail.
    pub fn pos_at(&self, t: f64) -> Vec3 {
        let v = |p: [f64; 3]| Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32);
        let n = self.trail.len();
        if n == 0 { return self.pos(); }
        let i = self.trail.partition_point(|q| q.0 <= t);
        if i == 0 { return v(self.trail[0].1); }
        if i >= n { return v(self.trail[n - 1].1); }
        let (ta, pa) = self.trail[i - 1];
        let (tb, pb) = self.trail[i];
        let f = if tb > ta { (t - ta) / (tb - ta) } else { 1.0 };
        v(add_s(pa, f, sub3(pb, pa)))
    }

    fn push_trail(&mut self) {
        self.trail.push_back((self.t, self.x));
        if self.trail.len() > TRAIL_POINTS { self.trail.pop_front(); }
    }

    /// Velocity: unit vector along j (its own direction where no wave is heard).
    fn vel(&self, x: [f64; 3], own: [f64; 3]) -> [f64; 3] {
        let (j, rho) = self.bench.current(x, self.hear);
        let nj = norm3(j);
        if rho > 1e-300 && nj > 1e-12 * self.bench.lev.k * rho.max(1e-300) {
            [j[0] / nj, j[1] / nj, j[2] / nj]
        } else {
            own
        }
    }

    /// First object crossed by the chord a→b (sources and passive planes excluded).
    fn hits(&self, a: [f64; 3], b: [f64; 3]) -> Option<(usize, f64)> {
        let dz = b[2] - a[2];
        if dz <= 0.0 { return None; }
        let mut best: Option<(usize, f64)> = None;
        for (i, o) in self.bench.objs.iter().enumerate() {
            if matches!(o.kind, Kind::Source | Kind::Plane) { continue; }
            let tt = (o.z - a[2]) / dz;
            if (0.0..=1.0).contains(&tt) && best.map_or(true, |q| tt < q.1) { best = Some((i, tt)); }
        }
        best
    }

    /// One Dormand–Prince attempt (accepted or rejected), with the events of
    /// the objects crossed. Python `Photons.run`, one photon.
    pub fn attempt(&mut self) {
        let b = self.bench.clone();
        let r = b.rule;
        let lam = b.lam;
        let (atol, hmin, hmax, delta) = (r.atol * lam, r.hmin * lam, r.hmax * lam, r.delta * lam);
        // cap: do not fly more than (distance to the next object + δ) in one step
        let cap = match self.hits(self.x, add_s(self.x, 1e6, self.d)) {
            Some((_, f)) => f * 1e6 * norm3(self.d),
            None => f64::INFINITY,
        } + delta;
        let h = self.h.min(cap).min(hmax);
        let mut ks = [[0.0; 3]; 7];
        ks[0] = match self.k1 { Some(k) => k, None => self.vel(self.x, self.d) };
        for s in 1..7 {
            let mut y = self.x;
            for (q, kq) in ks.iter().enumerate().take(s) { y = add_s(y, h * DP_A[s][q], *kq); }
            ks[s] = self.vel(y, self.d);
        }
        let (mut y5, mut y4) = (self.x, self.x);
        for s in 0..7 {
            y5 = add_s(y5, h * DP_B5[s], ks[s]);
            y4 = add_s(y4, h * DP_B4[s], ks[s]);
        }
        let err = norm3(sub3(y5, y4));
        let acc = err <= atol || h <= hmin * 1.0001;
        let fac = (0.9 * (atol / err.max(1e-300)).powf(0.2)).clamp(0.2, 5.0);
        let hnew = (h * fac).clamp(hmin, hmax);
        let mut hit = None;
        if acc {
            let (xo, xn) = (self.x, y5);
            hit = self.hits(xo, xn);
            let fr = hit.map_or(f64::INFINITY, |q| q.1);
            let dz = xn[2] - xo[2];
            if dz > 0.0 {
                for (pl, o) in b.objs.iter().enumerate().filter(|(_, o)| o.kind == Kind::Plane) {
                    let fp = (o.z - xo[2]) / dz;
                    if (0.0..=1.0).contains(&fp) && fp <= fr { self.rec.push((pl, add_s(xo, fp, sub3(xn, xo)))); }
                }
            }
            if hit.is_none() {
                self.x = xn;
                self.d = ks[6];
                self.t += h;
                self.steps += 1;
                self.k1 = Some(ks[6]);
                self.push_trail();
            }
        }
        self.h = hnew;
        if let Some((oi, fr)) = hit {
            let (xo, xn) = (self.x, y5);
            let chord = sub3(xn, xo);
            let xc = add_s(xo, fr, chord);
            let n = norm3(chord).max(1e-300);
            let din = [chord[0] / n, chord[1] / n, chord[2] / n];
            let tc = self.t + fr * h;
            self.k1 = None;
            self.event(oi, xc, din, tc, delta);
        }
        if self.end.is_none() {
            let bb = b.bbox;
            let out = (0..3).any(|a| self.x[a] < bb[2 * a] || self.x[a] > bb[2 * a + 1]);
            if out { self.end = Some(End::Lost); } else if self.steps > MAX_STEPS { self.end = Some(End::Stuck); }
        }
    }

    fn event(&mut self, oi: usize, xc: [f64; 3], dn: [f64; 3], tc: f64, delta: f64) {
        let b = self.bench.clone();
        let o = &b.objs[oi];
        if o.kind == Kind::Screen {
            self.x = xc;
            self.t = tc;
            self.end = Some(End::Click);
            self.push_trail();
            return;
        }
        let tr = o.mask.inside(xc[0], xc[1]) && o.t.is_some();
        if !tr {
            self.x = xc;
            self.t = tc;
            self.end = Some(End::Absorbed(oi));
            self.push_trail();
            return;
        }
        let cosn = dn[2];
        let mut xd = add_s(xc, delta / cosn.abs().max(1e-3), dn);
        let newhear = b.fan_of[oi];
        if o.absorbing {
            let (jin, _) = b.current([xc[0], xc[1], xc[2] - 1e-7], self.hear);
            let (jout, _) = b.current(xd, newhear);
            let (fin, fout) = (jin[2], jout[2]);
            let mut p = if fin > 0.0 { (fout / fin.max(1e-300)).clamp(0.0, 1.0) } else if fout > 0.0 { 1.0 } else { 0.0 };
            if let Some(eh) = self.edge(o, xc) {
                let hear = self.hear;
                let oz = o.z;
                let (np, pp) = self.transport(&eh, |q: [f64; 3]| b.current([q[0], q[1], oz], hear).0[2], |q| b.current(q, newhear).0);
                xd = np;
                p = pp;
            }
            let ok = self.rnd() < p;
            self.last_coin = Some((p, ok));
            if !ok {
                self.x = xc;
                self.t = tc;
                self.end = Some(End::Coin(oi));
                self.push_trail();
                return;
            }
        }
        self.x = xd;
        self.d = dn;
        self.t = tc + delta;
        self.h = delta / 4.0;
        self.hear = newhear;
        self.push_trail();
    }

    /// Nearest edge of the transmitting part of `o` if xc is within the edge zone.
    fn edge(&self, o: &Obj, xc: [f64; 3]) -> Option<EdgeHit> {
        let ue = self.bench.rule.edge_zone * self.bench.lam;
        if ue <= 0.0 { return None; }
        let (x, y) = (xc[0], xc[1]);
        let mut best: Option<EdgeHit> = None;
        match &o.mask {
            Mask::Full => return None,
            Mask::CircIn(cx, cy, r) | Mask::CircOut(cx, cy, r) => {
                let (dx, dy) = (x - cx, y - cy);
                let rr = dx.hypot(dy).max(1e-12);
                let (ex, ey) = (dx / rr, dy / rr);
                let e = [cx + r * ex, cy + r * ey, xc[2]];
                best = Some(if matches!(o.mask, Mask::CircIn(..)) {
                    EdgeHit { u: r - rr, e, dv: [-ex, -ey, 0.0], w: ue.min(*r), curv: -1.0 / r }
                } else {
                    EdgeHit { u: rr - r, e, dv: [ex, ey, 0.0], w: ue, curv: 1.0 / r }
                });
            }
            Mask::Rects(v) => {
                let mut u = f64::INFINITY;
                for q in v {
                    let (x0, x1, y0, y1) = (q[0], q[1], q[2], q[3]);
                    if !(x >= x0 && x <= x1 && y >= y0 && y <= y1) { continue; }
                    let cands = [
                        (x - x0, [x0, y], [1.0, 0.0], x1 - x0),
                        (x1 - x, [x1, y], [-1.0, 0.0], x1 - x0),
                        (y - y0, [x, y0], [0.0, 1.0], y1 - y0),
                        (y1 - y, [x, y1], [0.0, -1.0], y1 - y0),
                    ];
                    for (d, ep, dv, wid) in cands {
                        if d < u {
                            u = d;
                            best = Some(EdgeHit { u: d, e: [ep[0], ep[1], xc[2]], dv: [dv[0], dv[1], 0.0], w: ue.min(0.5 * wid), curv: 0.0 });
                        }
                    }
                }
            }
        }
        best.filter(|h| h.u < h.w)
    }

    /// v2: monotone (rank-keeping) map across an edge zone. On the cut E + v·dv,
    /// v ∈ [−W, W]: F_in = cumulative incident flux (opening side, v ≥ 0),
    /// F_out = cumulative outgoing j_z at δ (both sides). Survival p =
    /// F_out(W)/F_in(W) (≤ 1); new position: F_out(v*) = p·F_in(u).
    fn transport(&self, eh: &EdgeHit, phi_in: impl Fn([f64; 3]) -> f64 + Sync,
                 field_out: impl Fn([f64; 3]) -> [f64; 3] + Sync) -> ([f64; 3], f64) {
        let delta = self.bench.rule.delta * self.bench.lam;
        let nc = 161;
        let g: Vec<f64> = (0..nc).map(|i| {
            let tt = -1.0 + 2.0 * i as f64 / (nc - 1) as f64;
            (6.0 * tt).sinh() / 6f64.sinh()
        }).collect();
        let v: Vec<f64> = g.iter().map(|q| eh.w * q).collect();
        let jz: Vec<f64> = v.par_iter().map(|&vv| {
            let p = add_s(add_s(eh.e, vv, eh.dv), delta, [0.0, 0.0, 1.0]);
            field_out(p)[2].max(0.0) * (1.0 + eh.curv * vv)
        }).collect();
        let cum = |f: &[f64], x: &[f64]| -> Vec<f64> {
            let mut out = vec![0.0; f.len()];
            for i in 1..f.len() { out[i] = out[i - 1] + 0.5 * (f[i] + f[i - 1]) * (x[i] - x[i - 1]); }
            out
        };
        let fo = cum(&jz, &v);
        let i0 = nc / 2;
        let vp = &v[i0..];
        let fin: Vec<f64> = vp.par_iter().map(|&vv| phi_in(add_s(eh.e, vv, eh.dv)).max(0.0) * (1.0 + eh.curv * vv)).collect();
        let fi = cum(&fin, vp);
        let p = (fo[nc - 1] / fi[fi.len() - 1].max(1e-300)).clamp(0.0, 1.0);
        let target = p * interp(eh.u, vp, &fi);
        let vs = interp(target, &fo, &v);
        (add_s(add_s(eh.e, vs, eh.dv), delta, [0.0, 0.0, 1.0]), p)
    }

    /// Run until light-clock `t` (or until the photon ends).
    pub fn advance_to(&mut self, t: f64) {
        while self.end.is_none() && self.t < t { self.attempt(); }
    }

    /// Run to the end.
    pub fn run(&mut self) -> End {
        while self.end.is_none() { self.attempt(); }
        self.end.unwrap()
    }
}

/// np.interp(x, xp, fp) for nondecreasing xp.
fn interp(x: f64, xp: &[f64], fp: &[f64]) -> f64 {
    let n = xp.len();
    if x <= xp[0] { return fp[0]; }
    if x >= xp[n - 1] { return fp[n - 1]; }
    let j = xp.partition_point(|&q| q <= x).saturating_sub(1).min(n - 2);
    let dx = xp[j + 1] - xp[j];
    if dx > 0.0 { fp[j] + (fp[j + 1] - fp[j]) * (x - xp[j]) / dx } else { fp[j] }
}

// ─── The laboratories of the program ──────────────────────────────────────────

impl Bench3 {
    /// Objects of a laboratory (planes ⊥ z, in the order of the light).
    pub fn objects(lab: &crate::physics::Lab) -> Vec<Obj> {
        use crate::physics::{Experiment, WhichPath};
        let z2 = |j: crate::complex::Jones| -> J2 { [Z::from_c(j[0]), Z::from_c(j[1])] };
        let jones = z2(lab.polarization.jones());
        let mut objs = Vec::new();
        let gauss_source = |z: f64, mask: Mask| {
            let mut o = Obj::new("source", Kind::Source, z, mask);
            o.src = Some(SrcBeam { w0: lab.waist as f64, w0y: None, zw: (lab.waist_z - z as f32) as f64, jones });
            o
        };
        match lab.experiment {
            Experiment::FreeBeam | Experiment::Bell => objs.push(gauss_source(lab.source_z as f64, Mask::Full)),
            Experiment::DoubleSlit => {
                objs.push(gauss_source(lab.source_z as f64, Mask::Full));
                let (c, h) = lab.slits();
                let hy = 0.5 * lab.slit_height as f64;
                let rects: Vec<[f64; 4]> = c.iter().map(|&cc| [(cc - h) as f64, (cc + h) as f64, -hy, hy]).collect();
                let mut wall = Obj::new("slits", Kind::Aperture, lab.wall_z as f64, Mask::Rects(rects));
                if lab.which_path != WhichPath::Off {
                    wall.rect_t = lab.slit_filters().iter().map(|f| f.map(|j| projector(z2(j)))).collect();
                }
                objs.push(wall);
            }
            Experiment::Polarizers => {
                objs.push(gauss_source(lab.source_z as f64, Mask::Full));
                let mut ps: Vec<(usize, f32, f32)> = lab.polarizers.iter().enumerate()
                    .filter(|(_, p)| p.enabled && p.z > lab.source_z && p.z < lab.screen_z)
                    .map(|(i, p)| (i, p.z, p.angle)).collect();
                ps.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
                for (i, z, a) in ps {
                    let mut o = Obj::new(&format!("polarizer {}", i + 1), Kind::Polarizer(i), z as f64, Mask::Full);
                    o.t = Some(projector(z2(crate::complex::linear(a))));
                    objs.push(o);
                }
            }
            Experiment::Airy => {
                // the illuminated aperture is the source (plane wave), as in lab/base_exp3d.py
                let mut o = Obj::new("aperture", Kind::Source, lab.wall_z as f64, Mask::CircIn(0.0, 0.0, 0.5 * lab.aperture_d as f64));
                o.src = Some(SrcBeam { w0: f64::INFINITY, w0y: None, zw: 0.0, jones });
                objs.push(o);
            }
            Experiment::Arago => {
                objs.push(gauss_source(lab.wall_z as f64, Mask::CircOut(0.0, 0.0, lab.disk_r as f64)));
                objs[0].name = "disk".into();
            }
        }
        objs.push(Obj::new("screen", Kind::Screen, lab.screen_z as f64, Mask::Full));
        objs
    }

    pub fn lab_box(lab: &crate::physics::Lab) -> [f64; 6] {
        let side = (150.0f64).max(4.0 * lab.screen_half as f64);
        let z0 = lab.source_z.min(lab.wall_z) as f64 - 1.0;
        [-side, side, -side, side, z0, lab.screen_z as f64 + 1.0]
    }

    pub fn from_lab(lab: &crate::physics::Lab) -> Self {
        Self::build(format!("{lab:?}"), lab.lambda as f64, Self::objects(lab), Self::lab_box(lab))
    }
}

// ─── Estimators shared by the panel and the tests (lab/base_exp2d.py) ─────────

/// Weighted least-squares polynomial fit (numpy.polyfit with weights w), degree 3.
/// Returns coefficients c0 + c1 x + c2 x² + c3 x³.
fn polyfit3(x: &[f64], y: &[f64], w: &[f64]) -> Option<[f64; 4]> {
    let mut a = [[0.0f64; 5]; 4];
    for i in 0..x.len() {
        let ww = w[i] * w[i];
        let p = [1.0, x[i], x[i] * x[i], x[i] * x[i] * x[i]];
        for r in 0..4 {
            for c in 0..4 { a[r][c] += ww * p[r] * p[c]; }
            a[r][4] += ww * p[r] * y[i];
        }
    }
    for c in 0..4 {
        let piv = (c..4).max_by(|&i, &j| a[i][c].abs().partial_cmp(&a[j][c].abs()).unwrap())?;
        a.swap(c, piv);
        if a[c][c].abs() < 1e-300 { return None; }
        for r in 0..4 {
            if r != c {
                let f = a[r][c] / a[c][c];
                for k in c..5 { a[r][k] -= f * a[c][k]; }
            }
        }
    }
    Some([a[0][4] / a[0][0], a[1][4] / a[1][1], a[2][4] / a[2][2], a[3][4] / a[3][3]])
}

/// Position of a local minimum of the click density (`min_position` of
/// lab/base_exp2d.py): weighted cubic fit of the histogram (30 bins in
/// guess ± half), divided by `norm(x)` (e.g. x for a radial density),
/// stationary point nearest to the guess.
pub fn min_position(x: &[f64], guess: f64, half: f64, norm: impl Fn(f64) -> f64) -> Option<f64> {
    let sel: Vec<f64> = x.iter().copied().filter(|v| (v - guess).abs() < half).collect();
    if sel.len() < 30 { return None; }
    let nb = 30;
    let mut h = vec![0.0f64; nb];
    let w = 2.0 * half / nb as f64;
    for v in &sel {
        let k = (((v - (guess - half)) / w).floor() as isize).clamp(0, nb as isize - 1) as usize;
        h[k] += 1.0;
    }
    let c: Vec<f64> = (0..nb).map(|i| -half + (i as f64 + 0.5) * w).collect();
    let nn: Vec<f64> = c.iter().map(|&ci| norm(ci + guess)).collect();
    let y: Vec<f64> = h.iter().zip(&nn).map(|(a, b)| a / b).collect();
    let wt: Vec<f64> = h.iter().zip(&nn).map(|(a, b)| b / a.max(1.0).sqrt()).collect();
    let p = polyfit3(&c, &y, &wt)?;
    // p' = c1 + 2 c2 x + 3 c3 x² ; p'' = 2 c2 + 6 c3 x
    let (qa, qb, qc) = (3.0 * p[3], 2.0 * p[2], p[1]);
    let mut roots = Vec::new();
    if qa.abs() < 1e-300 {
        if qb.abs() > 1e-300 { roots.push(-qc / qb); }
    } else {
        let disc = qb * qb - 4.0 * qa * qc;
        if disc >= 0.0 {
            let s = disc.sqrt();
            roots.push((-qb + s) / (2.0 * qa));
            roots.push((-qb - s) / (2.0 * qa));
        }
    }
    roots.into_iter()
        .filter(|&q| q.abs() < half && 2.0 * p[2] + 6.0 * p[3] * q > 0.0)
        .min_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap())
        .map(|q| guess + q)
}

/// Airy first zero radius L·tan(asin(1.2197 λ/D)) for the panel.
pub fn airy_first_zero(lam: f64, d: f64, l: f64) -> f64 { l * (1.2197 * lam / d).asin().tan() }

#[cfg(test)]
impl Bench3 {
    /// Level-0 sum of fan `fi` at p with the scalar kernel `add3` (reference
    /// for the 4-lane far-field sum).
    pub fn field_l0_scalar(&self, p: [f64; 3], fi: usize) -> ([Z; 3], [Z; 9]) {
        let fan = &self.fans[fi];
        let lv = self.lev;
        let dz = p[2] - fan.z0;
        let mut f: Fld = ([Z::ZERO; 3], [Z::ZERO; 9]);
        for o in &fan.l0 {
            let r = [p[0] - o.p[0], p[1] - o.p[1], dz];
            let rr = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
            add3(lv.k, lv.w(0, rr), r, o.a, &mut f.0, &mut f.1);
        }
        f
    }
}

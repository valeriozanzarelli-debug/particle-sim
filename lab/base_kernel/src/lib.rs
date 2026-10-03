//! base_kernel — the field sums of the BASE photon rule (docs/REGOLA_BASE.md, lab/base_rule.py).
//!
//! Nothing here decides anything: it only evaluates, for a photon ball at P, the sum of the
//! wavelets of the wave balls it hears (and their gradient). The rule (who is heard, the levels,
//! the blend, the kernels) is written in lab/base_rule.py, which also contains a numpy reference
//! of the same arithmetic (`field2_numpy`, `field3_numpy`) used to test this library.
//!
//! Lengths in units of λ is NOT assumed: k is passed in the scene.
//!
//! 2-D scene buffer (f64):
//!   [2, k, r0, ds0, lmax, nfans, nblock, tog_obj, tog_init]
//!   nblock × [obj, x1, z1, x2, z2]                      (every object's segment: rays crossing it are blocked)
//!   nfans  × [obj, cx, cz, tx, tz, nx, nz, glo, ghi, nint, (a, b)×nint,
//!             sp0, sds, nsamp, (re0, im0, re1, im1)×nsamp, bthr]
//! 3-D scene buffer (f64):
//!   [3, k, r0, ds0, lmax, nfans]
//!   nfans  × [obj, z0, gx0, gx1, gy0, gy1, mtype, nmp, mparams×nmp,
//!             sx0, sy0, sds, nx, ny, (re0, im0, re1, im1)×(nx·ny), bthr]
//!   mtype: 0 full plane, 1 inside circle [cx, cy, r], 2 outside circle [cx, cy, r],
//!          3 union of rectangles [n, (x0, x1, y0, y1)×n]
use std::f64::consts::{FRAC_PI_4, PI};
use std::ops::{Add, AddAssign, Mul, Sub};
use std::os::raw::c_void;

// ─── complex ────────────────────────────────────────────────────────────────────────
#[derive(Clone, Copy, Default, Debug)]
pub struct C {
    pub re: f64,
    pub im: f64,
}
impl C {
    #[inline] pub const fn new(re: f64, im: f64) -> C { C { re, im } }
    #[inline] pub fn expi(t: f64) -> C { let (s, c) = t.sin_cos(); C::new(c, s) }
    #[inline] pub fn s(self, k: f64) -> C { C::new(self.re * k, self.im * k) }
    #[inline] pub fn mul_i(self) -> C { C::new(-self.im, self.re) }
    #[inline] pub fn mul_neg_i(self) -> C { C::new(self.im, -self.re) }
    #[inline] pub fn norm2(self) -> f64 { self.re * self.re + self.im * self.im }
}
impl Add for C { type Output = C; #[inline] fn add(self, o: C) -> C { C::new(self.re + o.re, self.im + o.im) } }
impl Sub for C { type Output = C; #[inline] fn sub(self, o: C) -> C { C::new(self.re - o.re, self.im - o.im) } }
impl AddAssign for C { #[inline] fn add_assign(&mut self, o: C) { self.re += o.re; self.im += o.im; } }
impl Mul for C {
    type Output = C;
    #[inline] fn mul(self, o: C) -> C { C::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re) }
}

// ─── Hankel H0(1), H1(1) (same code as src/photon/pilot.rs, validated there to 2e-7) ───
const ASY0: [f64; 6] = [1.0, -1.0 / 8.0, 9.0 / 128.0, -225.0 / 3072.0, 11025.0 / 98304.0, -893025.0 / 3932160.0];
const ASY1: [f64; 6] = [1.0, 3.0 / 8.0, -15.0 / 128.0, 315.0 / 3072.0, -14175.0 / 98304.0, 1091475.0 / 3932160.0];
#[inline]
fn asy_series(a: &[f64; 6], v: f64) -> C {
    let mut s = C::new(a[5], 0.0);
    for k in (0..5).rev() { s = s.mul_i().s(v) + C::new(a[k], 0.0); }
    s
}
#[inline]
fn poly(t: f64, c: &[f64; 7]) -> f64 {
    let mut s = c[6];
    for k in (0..6).rev() { s = s * t + c[k]; }
    s
}
pub fn hankel01(u: f64) -> (C, C) {
    if u >= 12.0 {
        let v = 1.0 / u;
        let e = C::expi(u - FRAC_PI_4).s((2.0 / (PI * u)).sqrt());
        (e * asy_series(&ASY0, v), (e * asy_series(&ASY1, v)).mul_neg_i())
    } else if u >= 3.0 {
        let t = 3.0 / u;
        let r = 1.0 / u.sqrt();
        let f0 = poly(t, &[0.797_884_56, -0.000_000_77, -0.005_527_40, -0.000_095_12, 0.001_372_37, -0.000_728_05, 0.000_144_76]);
        let th0 = u + poly(t, &[-0.785_398_16, -0.041_663_97, -0.000_039_54, 0.002_625_73, -0.000_541_25, -0.000_293_33, 0.000_135_58]);
        let f1 = poly(t, &[0.797_884_56, 0.000_001_56, 0.016_596_67, 0.000_171_05, -0.002_495_11, 0.001_136_53, -0.000_200_33]);
        let th1 = u + poly(t, &[-2.356_194_49, 0.124_996_12, 0.000_056_50, -0.006_378_79, 0.000_743_48, 0.000_798_24, -0.000_291_66]);
        (C::expi(th0).s(r * f0), C::expi(th1).s(r * f1))
    } else {
        let t = (u / 3.0) * (u / 3.0);
        let lg = (2.0 / PI) * (0.5 * u).ln();
        let j0 = poly(t, &[1.0, -2.249_999_7, 1.265_620_8, -0.316_386_6, 0.044_447_9, -0.003_944_4, 0.000_210_0]);
        let y0 = lg * j0 + poly(t, &[0.367_466_91, 0.605_593_66, -0.743_503_84, 0.253_001_17, -0.042_612_14, 0.004_279_16, -0.000_248_46]);
        let j1 = u * poly(t, &[0.5, -0.562_499_85, 0.210_935_73, -0.039_542_89, 0.004_433_19, -0.000_317_61, 0.000_011_09]);
        let y1 = (lg * u * j1 + poly(t, &[-0.636_619_8, 0.221_209_1, 2.168_270_9, -1.316_482_7, 0.312_395_1, -0.040_097_6, 0.002_787_3])) / u;
        (C::new(j0, y0), C::new(j1, y1))
    }
}

// ─── levels: which wave balls a photon at distance R hears ──────────────────────────
#[derive(Clone, Copy)]
struct Lev {
    k: f64,
    r0: f64,
    ds0: f64,
    lmax: usize,
}
impl Lev {
    /// weight of level l for an origin at distance R (smoothstep blend in log2 R, see base_rule.py)
    #[inline]
    fn w(&self, l: usize, r: f64) -> f64 {
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
    /// level l is heard only for R < rmax(l)
    #[inline]
    fn rmax(&self, l: usize) -> f64 {
        if l == 0 { f64::INFINITY } else { self.r0 * 2f64.powi(2 - l as i32) }
    }
    #[inline]
    fn ds(&self, l: usize) -> f64 { self.ds0 / 2f64.powi(l as i32) }
}

// ─── Catmull–Rom interpolation of the smooth boundary value B (Jones, 2 complex) ────
#[inline]
fn cr_weights(f: f64) -> [f64; 4] {
    let f2 = f * f;
    let f3 = f2 * f;
    [(-f3 + 2.0 * f2 - f) * 0.5, (3.0 * f3 - 5.0 * f2 + 2.0) * 0.5, (-3.0 * f3 + 4.0 * f2 + f) * 0.5, (f3 - f2) * 0.5]
}

// ─── 2-D scene ──────────────────────────────────────────────────────────────────────
#[derive(Clone, Copy)]
struct Org2 {
    p: [f64; 2],
    a: [C; 2],
}
struct Fan2 {
    obj: i64,
    c: [f64; 2],
    t: [f64; 2],
    n: [f64; 2],
    glo: f64,
    ghi: f64,
    ints: Vec<[f64; 2]>,
    sp0: f64,
    sds: f64,
    samp: Vec<[C; 2]>,
    bthr: f64,
    l0: Vec<Org2>,
}
struct Block {
    obj: i64,
    p1: [f64; 2],
    p2: [f64; 2],
}
struct Scene2 {
    lev: Lev,
    fans: Vec<Fan2>,
    blocks: Vec<Block>,
    tog_obj: i64,
    tog_init: bool,
}

impl Fan2 {
    fn interp(&self, s: f64) -> [C; 2] {
        let u = (s - self.sp0) / self.sds;
        let j = u.floor();
        let w = cr_weights(u - j);
        let j = j as i64;
        let mut b = [C::default(); 2];
        for m in 0..4 {
            let idx = j - 1 + m as i64;
            if idx < 0 || idx >= self.samp.len() as i64 { continue; }
            let v = self.samp[idx as usize];
            b[0] += v[0].s(w[m]);
            b[1] += v[1].s(w[m]);
        }
        b
    }
    /// origins of level `ds` with tangent coordinate in [sa, sb] (cells clipped to [glo, ghi])
    fn origins(&self, ds: f64, sa: f64, sb: f64, out: &mut Vec<Org2>) {
        let lo = sa.max(self.glo);
        let hi = sb.min(self.ghi);
        if hi <= lo { return; }
        let k0 = ((lo - self.glo) / ds).floor() as i64;
        let k1 = ((hi - self.glo) / ds).ceil() as i64;
        for kk in k0.max(0)..k1 {
            let a = self.glo + kk as f64 * ds;
            let b = (a + ds).min(self.ghi);
            let mut len = 0.0;
            let mut mom = 0.0;
            for iv in &self.ints {
                let l = a.max(iv[0]);
                let r = b.min(iv[1]);
                if r > l {
                    len += r - l;
                    mom += (r - l) * 0.5 * (l + r);
                }
            }
            if len <= 0.0 { continue; }
            let s = mom / len;
            let bb = self.interp(s);
            if (bb[0].norm2() + bb[1].norm2()).sqrt() < self.bthr { continue; }
            out.push(Org2 {
                p: [self.c[0] + s * self.t[0], self.c[1] + s * self.t[1]],
                a: [bb[0].s(len), bb[1].s(len)],
            });
        }
    }
}

#[inline]
fn seg_cross(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> Option<f64> {
    // intersection of segment a→b with segment c→d; returns parameter t on a→b (0<t<1)
    let r = [b[0] - a[0], b[1] - a[1]];
    let s = [d[0] - c[0], d[1] - c[1]];
    let den = r[0] * s[1] - r[1] * s[0];
    if den.abs() < 1e-300 { return None; }
    let q = [c[0] - a[0], c[1] - a[1]];
    let t = (q[0] * s[1] - q[1] * s[0]) / den;
    let u = (q[0] * r[1] - q[1] * r[0]) / den;
    if t > 1e-12 && t < 1.0 - 1e-12 && u >= 0.0 && u <= 1.0 { Some(t) } else { None }
}

impl Scene2 {
    #[inline]
    fn present(&self, obj: i64, tau: f64, tog: f64) -> bool {
        if obj != self.tog_obj { return true; }
        if tau < tog { self.tog_init } else { !self.tog_init }
    }
    #[inline]
    fn blocked(&self, o: [f64; 2], p: [f64; 2], t: f64, tog: f64, own: i64) -> bool {
        for b in &self.blocks {
            if b.obj == own { continue; }
            if let Some(tt) = seg_cross(o, p, b.p1, b.p2) {
                let dist = (1.0 - tt) * ((p[0] - o[0]).hypot(p[1] - o[1]));
                if self.present(b.obj, t - dist, tog) { return true; }
            }
        }
        false
    }
    /// ψ, ∂xψ, ∂zψ (2 Jones components each) at p.
    fn field(&self, p: [f64; 2], t: f64, tog: f64, excl: i64, scratch: &mut Vec<Org2>) -> [C; 6] {
        let lv = self.lev;
        let k = lv.k;
        let mut out = [C::default(); 6];
        for fan in &self.fans {
            if fan.obj == excl { continue; }
            let d = [p[0] - fan.c[0], p[1] - fan.c[1]];
            let zn = d[0] * fan.n[0] + d[1] * fan.n[1];
            if zn <= 1e-12 { continue; }
            let sp = d[0] * fan.t[0] + d[1] * fan.t[1];
            for l in 0..=lv.lmax {
                let list: &[Org2] = if l == 0 {
                    &fan.l0
                } else {
                    let rm = lv.rmax(l);
                    if zn >= rm { continue; }
                    let half = (rm * rm - zn * zn).sqrt();
                    scratch.clear();
                    fan.origins(lv.ds(l), sp - half, sp + half, scratch);
                    &scratch[..]
                };
                for o in list {
                    let rx = p[0] - o.p[0];
                    let rz = p[1] - o.p[1];
                    let r = (rx * rx + rz * rz).sqrt();
                    let w = lv.w(l, r);
                    if w == 0.0 { continue; }
                    if !self.present(fan.obj, t - r, tog) { continue; }
                    if self.blocked(o.p, p, t, tog, fan.obj) { continue; }
                    let znn = rx * fan.n[0] + rz * fan.n[1];
                    let inv = 1.0 / r;
                    let (h0, h1) = hankel01(k * r);
                    let f = h1.s(inv);
                    let fp = (h0.s(k) - h1.s(2.0 * inv)).s(inv);
                    let c = 0.5 * k * w;
                    // G = (ik/2) zn F ; ∇G = (ik/2)[n F + zn Fp r/R]
                    let g = f.s(c * znn).mul_i();
                    let gx = (f.s(fan.n[0]) + fp.s(znn * rx * inv)).s(c).mul_i();
                    let gz = (f.s(fan.n[1]) + fp.s(znn * rz * inv)).s(c).mul_i();
                    for q in 0..2 {
                        out[q] += g * o.a[q];
                        out[2 + q] += gx * o.a[q];
                        out[4 + q] += gz * o.a[q];
                    }
                }
            }
        }
        out
    }
}

// ─── 3-D scene (object planes ⊥ z) ──────────────────────────────────────────────────
#[derive(Clone, Copy)]
struct Org3 {
    p: [f64; 2],
    a: [C; 2],
}
enum Mask3 {
    Full,
    CircIn(f64, f64, f64),
    CircOut(f64, f64, f64),
    Rects(Vec<[f64; 4]>),
}
impl Mask3 {
    #[inline]
    fn inside(&self, x: f64, y: f64) -> bool {
        match self {
            Mask3::Full => true,
            Mask3::CircIn(cx, cy, r) => (x - cx).powi(2) + (y - cy).powi(2) < r * r,
            Mask3::CircOut(cx, cy, r) => (x - cx).powi(2) + (y - cy).powi(2) >= r * r,
            Mask3::Rects(v) => v.iter().any(|q| x >= q[0] && x < q[1] && y >= q[2] && y < q[3]),
        }
    }
    /// (area inside, centroid) of the cell [x0,x1]×[y0,y1]
    fn cell(&self, x0: f64, x1: f64, y0: f64, y1: f64) -> (f64, f64, f64) {
        let full = ((x1 - x0) * (y1 - y0), 0.5 * (x0 + x1), 0.5 * (y0 + y1));
        let circ = |cx: f64, cy: f64, r: f64, inside: bool| -> Option<(f64, f64, f64)> {
            let dx = (cx - x0).max(0.0).max(x1 - cx).max(0.0);
            let dy = (cy - y0).max(0.0).max(y1 - cy).max(0.0);
            let dmax = ((x0 - cx).abs().max((x1 - cx).abs())).hypot((y0 - cy).abs().max((y1 - cy).abs()));
            let nx = cx.clamp(x0, x1) - cx;
            let ny = cy.clamp(y0, y1) - cy;
            let dmin = nx.hypot(ny);
            let _ = (dx, dy);
            if dmax <= r { return Some(if inside { full } else { (0.0, 0.0, 0.0) }); }
            if dmin >= r { return Some(if inside { (0.0, 0.0, 0.0) } else { full }); }
            None
        };
        match self {
            Mask3::Full => full,
            Mask3::Rects(v) => {
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
            Mask3::CircIn(cx, cy, r) | Mask3::CircOut(cx, cy, r) => {
                let inside = matches!(self, Mask3::CircIn(..));
                if let Some(v) = circ(*cx, *cy, *r, inside) { return v; }
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
                if c == 0 { (0.0, 0.0, 0.0) } else {
                    (full.0 * c as f64 / (n * n) as f64, mx / c as f64, my / c as f64)
                }
            }
        }
    }
}
struct Fan3 {
    #[allow(dead_code)]
    obj: i64,
    z0: f64,
    g: [f64; 4],
    mask: Mask3,
    sx0: f64,
    sy0: f64,
    sds: f64,
    nx: usize,
    ny: usize,
    samp: Vec<[C; 2]>,
    bthr: f64,
    l0: Vec<Org3>,
}
impl Fan3 {
    fn interp(&self, x: f64, y: f64) -> [C; 2] {
        let u = (x - self.sx0) / self.sds;
        let v = (y - self.sy0) / self.sds;
        let iu = u.floor();
        let iv = v.floor();
        let wu = cr_weights(u - iu);
        let wv = cr_weights(v - iv);
        let iu = iu as i64;
        let iv = iv as i64;
        let mut b = [C::default(); 2];
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
    fn origins(&self, ds: f64, xa: f64, xb: f64, ya: f64, yb: f64, out: &mut Vec<Org3>) {
        let xlo = xa.max(self.g[0]);
        let xhi = xb.min(self.g[1]);
        let ylo = ya.max(self.g[2]);
        let yhi = yb.min(self.g[3]);
        if xhi <= xlo || yhi <= ylo { return; }
        let i0 = ((xlo - self.g[0]) / ds).floor() as i64;
        let i1 = ((xhi - self.g[0]) / ds).ceil() as i64;
        let j0 = ((ylo - self.g[2]) / ds).floor() as i64;
        let j1 = ((yhi - self.g[2]) / ds).ceil() as i64;
        self.origins_cells(ds, i0, i1, j0, j1, out);
    }
    /// origins of the cells i0..i1 x j0..j1 (cell (i, j) = [g0 + i ds, ...] x [g2 + j ds, ...], clipped to the box g)
    fn origins_cells(&self, ds: f64, i0: i64, i1: i64, j0: i64, j1: i64, out: &mut Vec<Org3>) {
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
                if (bb[0].norm2() + bb[1].norm2()).sqrt() < self.bthr { continue; }
                out.push(Org3 { p: [cx, cy], a: [bb[0].s(area), bb[1].s(area)] });
            }
        }
    }
}

const TILE: i64 = 16;
type TileKey = (usize, i64, i64);
/// fine-level origins are made on demand, by tiles of TILE x TILE cells, and kept (they are the same for
/// every photon: the wave is the same)
struct TileCache {
    map: std::sync::RwLock<std::collections::HashMap<TileKey, std::sync::Arc<Vec<Org3>>>>,
    count: std::sync::atomic::AtomicUsize,
}
impl TileCache {
    fn new() -> Self {
        TileCache { map: std::sync::RwLock::new(std::collections::HashMap::new()), count: std::sync::atomic::AtomicUsize::new(0) }
    }
}
impl Fan3 {
    fn tiles(&self, cache: &TileCache, l: usize, ds: f64, xa: f64, xb: f64, ya: f64, yb: f64,
             out: &mut Vec<std::sync::Arc<Vec<Org3>>>) {
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
                if let Some(v) = cache.map.read().unwrap().get(&key) {
                    out.push(v.clone());
                    continue;
                }
                let mut v = Vec::new();
                self.origins_cells(ds, tx * TILE, (tx + 1) * TILE, ty * TILE, (ty + 1) * TILE, &mut v);
                let n = v.len();
                let arc = std::sync::Arc::new(v);
                let mut m = cache.map.write().unwrap();
                if cache.count.load(std::sync::atomic::Ordering::Relaxed) > 30_000_000 {
                    m.clear();
                    cache.count.store(0, std::sync::atomic::Ordering::Relaxed);
                }
                m.insert(key, arc.clone());
                cache.count.fetch_add(n, std::sync::atomic::Ordering::Relaxed);
                out.push(arc);
            }
        }
    }
}

struct Scene3 {
    lev: Lev,
    fans: Vec<Fan3>,
    cache: Vec<TileCache>,
}

/// E (3) and ∂_b E_a (a-major: [a*3 + b]) from the tangential boundary values (Bx, By):
///   E_x = −2 Σ Bx ∂z g,  E_y = −2 Σ By ∂z g,  E_z = 2 Σ (Bx ∂x g + By ∂y g),  g = e^{ikR}/(4πR)
#[inline(always)]
fn add3(k: f64, w: f64, r: [f64; 3], a: [C; 2], e: &mut [C; 3], de: &mut [C; 9]) {
    let rr = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    let inv = 1.0 / rr;
    let ex = C::expi(k * rr).s(w * inv * (0.25 / PI));
    let h = C::new(-inv, k) * ex; // g'(R)
    let h2 = C::new(-k * k + 2.0 * inv * inv, -2.0 * k * inv) * ex; // g''(R)
    let hr = h.s(inv);
    let q = (h2 - hr).s(inv * inv);
    let bx = a[0];
    let by = a[1];
    // gradient of g: g_a = r_a h/R ; Hessian H_ab = d_ab h/R + r_a r_b q
    let bxh = bx * hr;
    let byh = by * hr;
    let bxq = bx * q;
    let byq = by * q;
    // E_x = -2 Bx g_z, E_y = -2 By g_z, E_z = 2 (Bx g_x + By g_y)
    e[0] += bxh.s(-2.0 * r[2]);
    e[1] += byh.s(-2.0 * r[2]);
    e[2] += (bxh.s(r[0]) + byh.s(r[1])).s(2.0);
    let sq = bxq.s(r[0]) + byq.s(r[1]);
    for b in 0..3 {
        // d_b E_x = -2 Bx H_bz
        let mut dx = bxq.s(r[b] * r[2]);
        let mut dy = byq.s(r[b] * r[2]);
        if b == 2 { dx += bxh; dy += byh; }
        de[b] += dx.s(-2.0);
        de[3 + b] += dy.s(-2.0);
        // d_b E_z = 2 (Bx H_bx + By H_by)
        let mut dz = sq.s(r[b]);
        if b == 0 { dz += bxh; }
        if b == 1 { dz += byh; }
        de[6 + b] += dz.s(2.0);
    }
}

impl Scene3 {
    fn field(&self, p: [f64; 3], fi: i64, scratch: &mut Vec<std::sync::Arc<Vec<Org3>>>) -> ([C; 3], [C; 9]) {
        let mut e = [C::default(); 3];
        let mut de = [C::default(); 9];
        if fi < 0 || fi as usize >= self.fans.len() { return (e, de); }
        let fan = &self.fans[fi as usize];
        let lv = self.lev;
        let dz = p[2] - fan.z0;
        if dz <= 1e-12 { return (e, de); }
        let two_r0 = 2.0 * lv.r0;
        for o in &fan.l0 {
            let r = [p[0] - o.p[0], p[1] - o.p[1], dz];
            let rr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            let w = if rr2 >= two_r0 * two_r0 { 1.0 } else { lv.w(0, rr2.sqrt()) };
            if w == 0.0 { continue; }
            add3(lv.k, w, r, o.a, &mut e, &mut de);
        }
        for l in 1..=lv.lmax {
            let rm = lv.rmax(l);
            if dz >= rm { continue; }
            let half = (rm * rm - dz * dz).sqrt();
            scratch.clear();
            fan.tiles(&self.cache[fi as usize], l, lv.ds(l), p[0] - half, p[0] + half, p[1] - half, p[1] + half, scratch);
            let rm2 = rm * rm;
            for t in scratch.iter() {
                for o in t.iter() {
                    let r = [p[0] - o.p[0], p[1] - o.p[1], dz];
                    let rr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    if rr2 >= rm2 { continue; }
                    let w = lv.w(l, rr2.sqrt());
                    if w == 0.0 { continue; }
                    add3(lv.k, w, r, o.a, &mut e, &mut de);
                }
            }
        }
        (e, de)
    }
}

// ─── parsing ──────────────────────────────────────────────────────────────────────
struct Rd<'a> {
    b: &'a [f64],
    i: usize,
}
impl<'a> Rd<'a> {
    fn f(&mut self) -> f64 { let v = self.b[self.i]; self.i += 1; v }
    fn u(&mut self) -> usize { self.f() as usize }
    fn c2(&mut self) -> [C; 2] { let a = self.f(); let b = self.f(); let c = self.f(); let d = self.f(); [C::new(a, b), C::new(c, d)] }
}

enum Scene {
    D2(Scene2),
    D3(Scene3),
}

fn parse(buf: &[f64]) -> Scene {
    let mut r = Rd { b: buf, i: 0 };
    let dim = r.u();
    let k = r.f();
    let r0 = r.f();
    let ds0 = r.f();
    let lmax = r.u();
    let lev = Lev { k, r0, ds0, lmax };
    let nfans = r.u();
    if dim == 2 {
        let nblock = r.u();
        let tog_obj = r.f() as i64;
        let tog_init = r.f() != 0.0;
        let mut blocks = Vec::new();
        for _ in 0..nblock {
            let obj = r.f() as i64;
            blocks.push(Block { obj, p1: [r.f(), r.f()], p2: [r.f(), r.f()] });
        }
        let mut fans = Vec::new();
        for _ in 0..nfans {
            let obj = r.f() as i64;
            let c = [r.f(), r.f()];
            let t = [r.f(), r.f()];
            let n = [r.f(), r.f()];
            let glo = r.f();
            let ghi = r.f();
            let nint = r.u();
            let mut ints = Vec::new();
            for _ in 0..nint { ints.push([r.f(), r.f()]); }
            let sp0 = r.f();
            let sds = r.f();
            let ns = r.u();
            let samp: Vec<[C; 2]> = (0..ns).map(|_| r.c2()).collect();
            let bthr = r.f();
            let mut fan = Fan2 { obj, c, t, n, glo, ghi, ints, sp0, sds, samp, bthr, l0: Vec::new() };
            let mut l0 = Vec::new();
            fan.origins(ds0, f64::NEG_INFINITY, f64::INFINITY, &mut l0);
            fan.l0 = l0;
            fans.push(fan);
        }
        Scene::D2(Scene2 { lev, fans, blocks, tog_obj, tog_init })
    } else {
        let mut fans = Vec::new();
        for _ in 0..nfans {
            let obj = r.f() as i64;
            let z0 = r.f();
            let g = [r.f(), r.f(), r.f(), r.f()];
            let mtype = r.u();
            let nmp = r.u();
            let mp: Vec<f64> = (0..nmp).map(|_| r.f()).collect();
            let mask = match mtype {
                0 => Mask3::Full,
                1 => Mask3::CircIn(mp[0], mp[1], mp[2]),
                2 => Mask3::CircOut(mp[0], mp[1], mp[2]),
                _ => {
                    let n = mp[0] as usize;
                    Mask3::Rects((0..n).map(|i| [mp[1 + 4 * i], mp[2 + 4 * i], mp[3 + 4 * i], mp[4 + 4 * i]]).collect())
                }
            };
            let sx0 = r.f();
            let sy0 = r.f();
            let sds = r.f();
            let nx = r.u();
            let ny = r.u();
            let samp: Vec<[C; 2]> = (0..nx * ny).map(|_| r.c2()).collect();
            let bthr = r.f();
            let mut fan = Fan3 { obj, z0, g, mask, sx0, sy0, sds, nx, ny, samp, bthr, l0: Vec::new() };
            let mut l0 = Vec::new();
            fan.origins(ds0, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, &mut l0);
            fan.l0 = l0;
            fans.push(fan);
        }
        let cache = (0..fans.len()).map(|_| TileCache::new()).collect();
        Scene::D3(Scene3 { lev, fans, cache })
    }
}

// ─── C ABI ────────────────────────────────────────────────────────────────────────
#[no_mangle]
pub extern "C" fn scene_new(buf: *const f64, len: usize) -> *mut c_void {
    let b = unsafe { std::slice::from_raw_parts(buf, len) };
    Box::into_raw(Box::new(parse(b))) as *mut c_void
}

#[no_mangle]
pub extern "C" fn scene_free(p: *mut c_void) {
    if !p.is_null() {
        unsafe { drop(Box::from_raw(p as *mut Scene)); }
    }
}

/// Origins of `level` of fan `fi` inside the box [lo0,hi0]×[lo1,hi1] (2-D: tangent interval [lo0,hi0]).
/// Writes up to `cap` rows of [x, z_or_y, re0, im0, re1, im1]; returns the total count.
#[no_mangle]
pub extern "C" fn fan_origins(p: *const c_void, fi: usize, level: usize, lo0: f64, hi0: f64, lo1: f64, hi1: f64,
                              out: *mut f64, cap: usize) -> usize {
    let sc = unsafe { &*(p as *const Scene) };
    let o = unsafe { std::slice::from_raw_parts_mut(out, cap * 6) };
    let mut put = |i: usize, a: f64, b: f64, c: [C; 2]| {
        if i < cap {
            o[6 * i..6 * i + 6].copy_from_slice(&[a, b, c[0].re, c[0].im, c[1].re, c[1].im]);
        }
    };
    match sc {
        Scene::D2(s) => {
            let mut v = Vec::new();
            s.fans[fi].origins(s.lev.ds(level), lo0, hi0, &mut v);
            for (i, q) in v.iter().enumerate() { put(i, q.p[0], q.p[1], q.a); }
            v.len()
        }
        Scene::D3(s) => {
            let mut v = Vec::new();
            s.fans[fi].origins(s.lev.ds(level), lo0, hi0, lo1, hi1, &mut v);
            for (i, q) in v.iter().enumerate() { put(i, q.p[0], q.p[1], q.a); }
            v.len()
        }
    }
}

fn split<F: Fn(usize, usize) + Sync>(n: usize, nthreads: usize, f: F) {
    let nt = nthreads.max(1).min(n.max(1));
    if nt <= 1 {
        f(0, n);
        return;
    }
    let chunk = (n + nt - 1) / nt;
    std::thread::scope(|s| {
        for t in 0..nt {
            let a = t * chunk;
            let b = ((t + 1) * chunk).min(n);
            if a >= b { continue; }
            let f = &f;
            s.spawn(move || f(a, b));
        }
    });
}

struct SendPtr(*mut f64);
unsafe impl Sync for SendPtr {}
unsafe impl Send for SendPtr {}
impl SendPtr { #[inline] fn get(&self) -> *mut f64 { self.0 } }

/// 2-D: out[n×12] = ψ0, ψ1, ∂xψ0, ∂xψ1, ∂zψ0, ∂zψ1 (re, im)
#[no_mangle]
pub extern "C" fn field2(p: *const c_void, n: usize, pts: *const f64, tph: *const f64, tog: *const f64, excl: i64,
                         out: *mut f64, nthreads: usize) {
    let sc = match unsafe { &*(p as *const Scene) } { Scene::D2(s) => s, _ => return };
    let pts = unsafe { std::slice::from_raw_parts(pts, 2 * n) };
    let tph = unsafe { std::slice::from_raw_parts(tph, n) };
    let tog = unsafe { std::slice::from_raw_parts(tog, n) };
    let op = SendPtr(out);
    split(n, nthreads, |a, b| {
        let o = unsafe { std::slice::from_raw_parts_mut(op.get().add(12 * a), 12 * (b - a)) };
        let mut scratch = Vec::new();
        for i in a..b {
            let v = sc.field([pts[2 * i], pts[2 * i + 1]], tph[i], tog[i], excl, &mut scratch);
            for q in 0..6 {
                o[12 * (i - a) + 2 * q] = v[q].re;
                o[12 * (i - a) + 2 * q + 1] = v[q].im;
            }
        }
    });
}

/// 3-D: out[n×24] = E_x, E_y, E_z, then ∂_b E_a for a in x,y,z, b in x,y,z (re, im)
#[no_mangle]
pub extern "C" fn field3(p: *const c_void, n: usize, pts: *const f64, hear: *const f64, out: *mut f64, nthreads: usize) {
    let sc = match unsafe { &*(p as *const Scene) } { Scene::D3(s) => s, _ => return };
    let pts = unsafe { std::slice::from_raw_parts(pts, 3 * n) };
    let hear = unsafe { std::slice::from_raw_parts(hear, n) };
    let op = SendPtr(out);
    split(n, nthreads, |a, b| {
        let o = unsafe { std::slice::from_raw_parts_mut(op.get().add(24 * a), 24 * (b - a)) };
        let mut scratch = Vec::new();
        for i in a..b {
            let (e, de) = sc.field([pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]], hear[i] as i64, &mut scratch);
            for q in 0..3 {
                o[24 * (i - a) + 2 * q] = e[q].re;
                o[24 * (i - a) + 2 * q + 1] = e[q].im;
            }
            for q in 0..9 {
                o[24 * (i - a) + 6 + 2 * q] = de[q].re;
                o[24 * (i - a) + 6 + 2 * q + 1] = de[q].im;
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn blend_partition_of_unity() {
        let lv = Lev { k: 2.0 * PI, r0: 1.1, ds0: 0.25, lmax: 7 };
        for i in 0..2000 {
            let r = 10f64.powf(-3.0 + 5.0 * i as f64 / 2000.0);
            let s: f64 = (0..=7).map(|l| lv.w(l, r)).sum();
            assert!((s - 1.0).abs() < 1e-12, "r={r} sum={s}");
            for l in 1..=7 { if r >= lv.rmax(l) { assert_eq!(lv.w(l, r), 0.0); } }
        }
    }
}

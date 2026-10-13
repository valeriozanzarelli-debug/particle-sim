//! How the free photon is DRAWN (the physics is in fotone3d.rs; nothing here changes the field).
//!
//! The look follows the photon laboratory (src/photon/render3d.rs): WAVE BALLS with their own phase colour
//! (alpha blending, soft rim) streaming forward, plus a few bright SURFERS with trails.
//!
//! • ENERGY BALLS (docs/FOTONE_PALLINE.md): tens of thousands of soft balls, each a small piece of the energy of the
//!   field. They are born ONCE with probability ∝ u = E² + c²B² and then every tick move with the energy velocity
//!   v = 2 E × cB / (E² + c²B²) (|v| ≤ c), read from the lattice by trilinear interpolation. Energy is conserved,
//!   so the cloud stays ∝ u by itself (exam F1 of the study: correlation 0.998 after 4 λ): no re-birth, no ages.
//!   Each ball is a ROUND GLOW (Valerio, 13/10: no stretching, no tails, no special surfers) coloured by the
//!   direction of E around the flight axis (circular light: every colour once per λ; linear: two colours, crest and
//!   trough); its size and light grow with u (`Bagliore`, docs/BAGLIORI.md). At the edge of a circular beam the
//!   glows turn around the axis: that swirl is the spin of the photon (exam F2: ±0.999 ħ from the balls; BAGLIORI S1–S3).
//!   (Stretching along E and the 4-sample tail are still in the code, off: `Vista::allunga`, `SCIA_QUAL`.)
//! • WAVEFRONTS: thin glowing sheets where E along the polarization has a crest (gold) or a trough (cyan)
//!   along the flight; the peak position inside the cell comes from a parabola through three cells, so the sheets
//!   glide smoothly instead of jumping from cell to cell.
//! • E and B streaks (off by default): a few short smooth strokes from the axis along E (wheel colours) and B
//!   (blue).
//! • Poynting S = E × B, spin E(t−1) × E(t), cross-section (off by default).

use rayon::prelude::*;

use crate::fotone3d::{base_trasversa, punto, Cf, Fotone3D, Lancio, Polar};
use crate::render::{phase_color, schermo, Ball};
use crate::rng::{hash01, Rng};

/// default number of energy balls (the quality presets choose their own: fotscena::PALLINE_QUAL)
pub const N_TRACCIANTI: usize = 32_000;
pub const N_SURF: usize = 6;
/// past positions kept per energy ball (one every 2 ticks): the tail shows the real curved path (the swirl)
pub const SCIA_MAX: usize = 4;
pub const L_SCIA: usize = 48;

/// Which layers are drawn.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vista {
    pub traccianti: bool,
    pub fronti: bool,
    /// balls stretched along the local E (the direction of the field, read by eye)
    pub allunga: bool,
    pub e_vettori: bool,
    pub b_vettori: bool,
    pub poynting: bool,
    pub spin: bool,
    pub sezione: bool,
}

impl Default for Vista {
    fn default() -> Self {
        Vista { traccianti: true, fronti: false, allunga: false, e_vettori: false, b_vettori: false, poynting: false, spin: false, sezione: false }
    }
}

/// The glow of an energy ball (study "bagliori", docs/BAGLIORI.md): its SIZE and its LIGHT follow the local energy
/// density u = (E² + c²B²)/max, so the core of the beam reads brighter and fuller than its edge, and the light of one
/// ball stays low enough that the thousands seen through the beam head-on do not saturate into a white disc.
/// size = s0 + s1 √u, light = min(a0 + a1 u^p, amax), colour = phase colour × c (c < 1 keeps the colour out of the
/// tone mapper's white shoulder: at 1 the HDR tone mapping turns every hue pastel). Drawing knobs, not physics
/// (counted in the document). `ATOMO_BAGLIORE="s0,s1,a0,a1,p,amax,c"` overrides them (only to try values without
/// recompiling).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bagliore { pub s0: f32, pub s1: f32, pub a0: f32, pub a1: f32, pub p: f32, pub amax: f32, pub c: f32 }

pub const BAGLIORE: Bagliore = Bagliore { s0: 0.45, s1: 0.85, a0: 0.03, a1: 0.42, p: 1.0, amax: 0.5, c: 1.0 };

impl Bagliore {
    pub fn attuale() -> Bagliore {
        static B: std::sync::OnceLock<Bagliore> = std::sync::OnceLock::new();
        *B.get_or_init(|| {
            let v: Vec<f32> = std::env::var("ATOMO_BAGLIORE").ok().map(|s| s.split(',').filter_map(|x| x.trim().parse().ok()).collect()).unwrap_or_default();
            if v.len() == 7 { Bagliore { s0: v[0], s1: v[1], a0: v[2], a1: v[3], p: v[4], amax: v[5], c: v[6] } } else { BAGLIORE }
        })
    }
    /// (size factor, light) of a ball at relative energy density u ∈ [0, 1]
    #[inline]
    pub fn di(&self, u: f32) -> (f32, f32) {
        let u = u.clamp(0.0, 1.0);
        (self.s0 + self.s1 * u.sqrt(), (self.a0 + self.a1 * u.powf(self.p)).min(self.amax))
    }
}

const BLU: [f32; 3] = [0.40, 0.62, 1.0];
const ORO: [f32; 3] = [1.0, 0.80, 0.36];
const CIANO: [f32; 3] = [0.35, 0.85, 1.0];

#[inline]
fn d3(a: &[f32; 3], b: &[f32; 3]) -> f32 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
#[inline]
fn x3(a: &[f32; 3], b: &[f32; 3]) -> [f32; 3] { [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]] }
fn f32v(v: [f64; 3]) -> [f32; 3] { [v[0] as f32, v[1] as f32, v[2] as f32] }

/// F at a point (cells from the box centre), trilinear, periodic.
pub fn campo_in(f: &Fotone3D, p: [f32; 3]) -> [Cf; 3] {
    let n = f.n;
    let h = (n / 2) as f32;
    let mut i0 = [0usize; 3];
    let mut w = [0f32; 3];
    for c in 0..3 {
        let x = p[c] + h;
        let fl = x.floor();
        w[c] = x - fl;
        i0[c] = (fl as i64).rem_euclid(n as i64) as usize;
    }
    let mut out = [Cf::ZERO; 3];
    for dx in 0..2 {
        for dy in 0..2 {
            for dz in 0..2 {
                let wt = (if dx == 1 { w[0] } else { 1.0 - w[0] }) * (if dy == 1 { w[1] } else { 1.0 - w[1] }) * (if dz == 1 { w[2] } else { 1.0 - w[2] });
                if wt == 0.0 { continue; }
                let (i, j, k) = ((i0[0] + dx) % n, (i0[1] + dy) % n, (i0[2] + dz) % n);
                let q = &f.f[(i * n + j) * n + k];
                for c in 0..3 { out[c].re += wt * q[c].re; out[c].im += wt * q[c].im; }
            }
        }
    }
    out
}

#[allow(dead_code)]
/// E (an array of [f32; 3] per cell) at a point, trilinear, periodic.
pub fn e_in(ep: &[[f32; 3]], n: usize, p: [f32; 3]) -> [f32; 3] {
    let h = (n / 2) as f32;
    let mut i0 = [0usize; 3];
    let mut w = [0f32; 3];
    for c in 0..3 {
        let x = p[c] + h;
        let fl = x.floor();
        w[c] = x - fl;
        i0[c] = (fl as i64).rem_euclid(n as i64) as usize;
    }
    let mut out = [0f32; 3];
    for dx in 0..2 { for dy in 0..2 { for dz in 0..2 {
        let wt = (if dx == 1 { w[0] } else { 1.0 - w[0] }) * (if dy == 1 { w[1] } else { 1.0 - w[1] }) * (if dz == 1 { w[2] } else { 1.0 - w[2] });
        if wt == 0.0 { continue; }
        let q = &ep[(((i0[0] + dx) % n) * n + (i0[1] + dy) % n) * n + (i0[2] + dz) % n];
        for c in 0..3 { out[c] += wt * q[c]; }
    }}}
    out
}

/// Energy velocity of the field at a point: v = 2 E × B / (E² + B²) (|v| ≤ 1), and the energy density.
fn vel_energia(f: &Fotone3D, p: [f32; 3]) -> ([f32; 3], f32) {
    let q = campo_in(f, p);
    let e = [q[0].re, q[1].re, q[2].re];
    let b = [q[0].im, q[1].im, q[2].im];
    let u = d3(&e, &e) + d3(&b, &b);
    if u <= 0.0 { return ([0.0; 3], 0.0); }
    let s = x3(&e, &b);
    let mut v = [2.0 * s[0] / u, 2.0 * s[1] / u, 2.0 * s[2] / u];
    let m = d3(&v, &v).sqrt();
    if m > 1.0 { for c in v.iter_mut() { *c /= m; } }
    (v, u)
}

/// The tracer balls that ride the field.
pub struct Traccianti {
    pub pos: Vec<[f32; 3]>,
    pub vel: Vec<[f32; 3]>,
    /// past positions, newest first, one every 2 ticks
    pub scia: Vec<[[f32; 3]; SCIA_MAX]>,
    /// how many of them are drawn (quality preset)
    pub scia_n: usize,
    passi: u64,
    pub eta: Vec<u16>,
    pub vita: Vec<u16>,
    /// the surfers: their trails (last = where they are)
    pub surf: Vec<Vec<[f32; 3]>>,
    rng: Rng,
}

impl Traccianti {
    pub fn new(f: &Fotone3D, seme: u64, quante: usize) -> Self {
        let mut t = Traccianti { pos: vec![], vel: vec![], scia: vec![], scia_n: 2, passi: 0, eta: vec![], vita: vec![], surf: vec![], rng: Rng::new(seme) };
        let nuovi = semina(f, quante.max(N_SURF), &mut t.rng);
        for p in nuovi {
            t.pos.push(p);
            t.vel.push([0.0; 3]);
            t.scia.push([p; SCIA_MAX]);
            t.eta.push(u16::MAX / 2);
            t.vita.push(u16::MAX);
        }
        // the surfers start among the brightest balls
        for i in 0..N_SURF { t.surf.push(vec![t.pos[(i * 997) % t.pos.len()]]); }
        t
    }

    /// One tick: every ball moves with the energy velocity of the field (backwards in time: against it).
    /// Balls are never re-born while the photon flies (energy is conserved); only a ball that has lost the light
    /// completely (u < 1e-4 of the peak: the integration error of the far tail) is put back where the light is.
    pub fn avanza(&mut self, f: &Fotone3D, verso: f32, centro: [f64; 3]) {
        let imax = f.f.par_iter().map(|p| p[0].norm2() + p[1].norm2() + p[2].norm2()).reduce(|| 0.0, f32::max).max(1e-30);
        // the tail: every 2 ticks forward the history shifts; going back in time it is cleared
        self.passi += 1;
        if verso > 0.0 && self.passi % 2 == 0 {
            self.scia.par_iter_mut().zip(self.pos.par_iter()).for_each(|(s, p)| { s.rotate_right(1); s[0] = *p; });
        } else if verso < 0.0 {
            self.scia.par_iter_mut().zip(self.pos.par_iter()).for_each(|(s, p)| { *s = [*p; SCIA_MAX]; });
        }
        let morti: Vec<usize> = self.pos.par_iter_mut().zip(self.vel.par_iter_mut()).zip(self.eta.par_iter_mut()).zip(self.vita.par_iter())
            .enumerate().filter_map(|(i, (((p, v), eta), _vita))| {
                let (vv, u) = vel_energia(f, *p);
                if u < 1e-4 * imax { return Some(i); }
                for c in 0..3 { p[c] += verso * vv[c]; }
                *v = vv;
                *eta = eta.saturating_add(1);
                None
            }).collect();
        // the surfers: carried by the same flow; one that falls out of the light starts again near the centre
        for k in 0..self.surf.len() {
            let p = *self.surf[k].last().unwrap();
            let (v, u) = vel_energia(f, p);
            if u < 0.05 * imax {
                let c = [centro[0] as f32, centro[1] as f32, centro[2] as f32];
                let r = [self.rng.gauss() as f32 * 2.5, self.rng.gauss() as f32 * 2.5, self.rng.gauss() as f32 * 2.5];
                self.surf[k] = vec![[c[0] + r[0], c[1] + r[1], c[2] + r[2]]];
                continue;
            }
            let q = [p[0] + verso * v[0], p[1] + verso * v[1], p[2] + verso * v[2]];
            if verso > 0.0 { self.surf[k].push(q); if self.surf[k].len() > L_SCIA { self.surf[k].remove(0); } }
            else if self.surf[k].len() > 1 { self.surf[k].pop(); } else { self.surf[k] = vec![q]; }
        }
        if morti.is_empty() { return; }
        let nuovi = semina(f, morti.len(), &mut self.rng);
        for (i, p) in morti.into_iter().zip(nuovi) {
            self.pos[i] = p;
            self.vel[i] = [0.0; 3];
            self.scia[i] = [p; SCIA_MAX];
            self.eta[i] = 0;
            self.vita[i] = u16::MAX;
        }
    }
}

/// `k` points drawn with probability ∝ |F|² (a cell, then a random point inside it).
pub fn semina(f: &Fotone3D, k: usize, rng: &mut Rng) -> Vec<[f32; 3]> {
    let n = f.n;
    let imax = f.f.par_iter().map(|p| p[0].norm2() + p[1].norm2() + p[2].norm2()).reduce(|| 0.0, f32::max).max(1e-30);
    // every other cell along each axis (a block of 2³ cells), then a random point in the block
    let celle: Vec<(usize, f32)> = (0..n / 2).into_par_iter().flat_map_iter(|a| {
        let mut v = Vec::new();
        for b in 0..n / 2 { for c in 0..n / 2 {
            let s = ((2 * a) * n + 2 * b) * n + 2 * c;
            let p = &f.f[s];
            let w = p[0].norm2() + p[1].norm2() + p[2].norm2();
            if w > 0.01 * imax { v.push((s, w)); }
        }}
        v
    }).collect();
    if celle.is_empty() { return vec![[0.0; 3]; k]; }
    let mut cdf = Vec::with_capacity(celle.len());
    let mut acc = 0.0f64;
    for c in &celle { acc += c.1 as f64; cdf.push(acc); }
    let h = (n / 2) as f32;
    (0..k).map(|_| {
        let u = rng.uni() * acc;
        let i = cdf.partition_point(|&c| c < u).min(celle.len() - 1);
        let s = celle[i].0;
        let (a, b, c) = (s / (n * n), (s / n) % n, s % n);
        [a as f32 - h + 2.0 * rng.uni() as f32 - 1.0, b as f32 - h + 2.0 * rng.uni() as f32 - 1.0, c as f32 - h + 2.0 * rng.uni() as f32 - 1.0]
    }).collect()
}

/// A smooth stroke: overlapping balls every 0.25 cells from `p` to `p + u l`.
fn tratto(out: &mut Vec<Ball>, p: [f64; 3], u: [f64; 3], l: f64, col: [f32; 3], alfa: f32, size: f32) {
    let k = ((l / 0.25).ceil() as usize).max(2);
    for q in 0..=k {
        let t = l * q as f64 / k as f64;
        let a = alfa * (0.35 + 0.65 * q as f32 / k as f32);
        out.push((schermo([p[0] + u[0] * t, p[1] + u[1] * t, p[2] + u[2] * t]).to_array(), size, [col[0], col[1], col[2], a]));
    }
}

/// Everything drawn for the photon (scene units = cells, box centred on the origin).
/// Returns (glowing balls, energy balls drawn with alpha blending, the stretch of each energy ball).
pub fn palline(f: &Fotone3D, l: &Lancio, centro: [f64; 3], tr: &Traccianti, e_prima: Option<&[[f32; 3]]>, vista: Vista, surfisti: bool) -> (Vec<Ball>, Vec<Ball>, Vec<[f32; 3]>) {
    let n = f.n;
    let dir = l.dir;
    let (e1d, e2d) = base_trasversa(dir);
    let nhat = if l.polar == Polar::Verticale { e2d } else { e1d };
    let (e1, e2, nh, df) = (f32v(e1d), f32v(e2d), f32v(nhat), f32v(dir));
    let (imax, emax2) = f.f.par_iter()
        .map(|p| (p[0].norm2() + p[1].norm2() + p[2].norm2(), p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re))
        .reduce(|| (0.0, 0.0), |a, b| (a.0.max(b.0), a.1.max(b.1)));
    let imax = imax.max(1e-30);
    let emax = emax2.sqrt().max(1e-30);
    let lam = l.lambda;
    let h = (n / 2) as f32;
    let mut out: Vec<Ball> = Vec::with_capacity(40_000);
    let mut onde: Vec<Ball> = Vec::with_capacity(3 * tr.pos.len());
    let mut allunga: Vec<[f32; 3]> = Vec::with_capacity(3 * tr.pos.len());
    // ── tracers (wave balls) with their tails, and the surfers ──
    if vista.traccianti {
        // stretch: up to 0.35 λ along E at the peak of the field (scene units = cells)
        let l_max = (0.35 * lam) as f32;
        let ns = tr.scia_n.min(SCIA_MAX);
        let bg = Bagliore::attuale();
        let parti: Vec<Vec<(Ball, [f32; 3])>> = tr.pos.par_iter().zip(tr.scia.par_iter()).zip(tr.eta.par_iter()).map(|((p, sc), eta)| {
            let q = campo_in(f, *p);
            let e = [q[0].re, q[1].re, q[2].re];
            let b = [q[0].im, q[1].im, q[2].im];
            let u = (d3(&e, &e) + d3(&b, &b)) / imax;
            // colour = direction of E around the flight axis
            let ang = d3(&e, &e2).atan2(d3(&e, &e1)) as f64;
            // a re-born ball fades in over 6 ticks
            let fade = (*eta as f32 / 6.0).min(1.0);
            let (size, luce) = bg.di(u);
            let a = fade * luce;
            let c = phase_color(ang, a);
            // the palette is written in sRGB; the Lab renders in linear HDR, so it is converted (≈ x^2.2), otherwise
            // every hue comes out pastel and thousands of them average to a white disc
            let c = [c[0].powf(2.2) * bg.c, c[1].powf(2.2) * bg.c, c[2].powf(2.2) * bg.c, c[3]];
            let st = if vista.allunga {
                let k = l_max / emax;
                schermo([(e[0] * k) as f64, (e[1] * k) as f64, (e[2] * k) as f64]).to_array()
            } else { [0.0; 3] };
            let sp = |x: &[f32; 3]| schermo([x[0] as f64, x[1] as f64, x[2] as f64]).to_array();
            let mut v = Vec::with_capacity(1 + ns);
            // the tail first (drawn behind the head): the same ball where it was 2, 4, 6, 8 ticks ago
            for j in (0..ns).rev() {
                let w = 0.5 * (1.0 - j as f32 / ns as f32);
                v.push(((sp(&sc[j]), 0.85 - 0.08 * j as f32, [c[0], c[1], c[2], w * a]), [0.0; 3]));
            }
            v.push(((sp(p), size, c), st));
            v
        }).collect();
        for t in parti { for (ball, st) in t { onde.push(ball); allunga.push(st); } }
        // no separate surfers with trails (Valerio, 13/10): the energy balls already ride the flow
        for s in tr.surf.iter().filter(|_| surfisti && false) {
            let k = s.len();
            for (j, p) in s.iter().enumerate() {
                let x = (j + 1) as f32 / k as f32;
                let pos = schermo([p[0] as f64, p[1] as f64, p[2] as f64]).to_array();
                out.push((pos, 0.5 + 0.7 * x, [1.0, 0.95, 0.85, 0.12 + 0.6 * x * x]));
            }
            if let Some(p) = s.last() {
                let pos = schermo([p[0] as f64, p[1] as f64, p[2] as f64]).to_array();
                out.push((pos, 2.6, [1.0, 0.97, 0.9, 1.0]));
                out.push((pos, 6.0, [1.0, 0.9, 0.7, 0.18]));
            }
        }
    }
    // ── wavefront sheets: crests and troughs of E·n̂ along the flight ──
    if vista.fronti {
        let o = [dir[0].round() as i64, dir[1].round() as i64, dir[2].round() as i64];
        let olen = ((o[0] * o[0] + o[1] * o[1] + o[2] * o[2]) as f64).sqrt() as f32;
        let ni = n as i64;
        let vicino = |s: usize, sgn: i64| -> usize {
            let (a, b, c) = ((s / (n * n)) as i64, ((s / n) % n) as i64, (s % n) as i64);
            let w = |x: i64| x.rem_euclid(ni) as usize;
            (w(a + sgn * o[0]) * n + w(b + sgn * o[1])) * n + w(c + sgn * o[2])
        };
        let en = |s: usize| { let p = &f.f[s]; p[0].re * nh[0] + p[1].re * nh[1] + p[2].re * nh[2] };
        let soglia = 0.30 * emax;
        let parti: Vec<Vec<Ball>> = (0..n).into_par_iter().map(|i| {
            let mut v = Vec::new();
            for j in 0..n {
                for k in 0..n {
                    let s = (i * n + j) * n + k;
                    let b = en(s);
                    if b.abs() < soglia { continue; }
                    let (a, c) = (en(vicino(s, -1)), en(vicino(s, 1)));
                    let picco = if b > 0.0 { b >= a && b > c } else { b <= a && b < c };
                    if !picco { continue; }
                    let den = a - 2.0 * b + c;
                    let off = if den.abs() > 1e-30 { (0.5 * (a - c) / den).clamp(-0.5, 0.5) } else { 0.0 };
                    let col = if b > 0.0 { ORO } else { CIANO };
                    let alfa = 0.10 + 0.32 * (b.abs() / emax);
                    for r in 0..3u64 {
                        let q = s as u64 * 7 + r * 2;
                        let (u1, u2) = (hash01(q) - 0.5, hash01(q + 1) - 0.5);
                        let pos = [0, 1, 2].map(|c2| {
                            let base = [i, j, k][c2] as f32 - h;
                            base + off * o[c2] as f32 / olen.max(1.0) * olen + u1 * e1[c2] + u2 * e2[c2]
                        });
                        v.push((schermo([pos[0] as f64, pos[1] as f64, pos[2] as f64]).to_array(), 1.15, [col[0], col[1], col[2], alfa]));
                    }
                }
            }
            v
        }).collect();
        for p in parti { out.extend(p); }
    }
    // ── E and B streaks from the axis (off by default; diagnostics) ──
    if vista.e_vettori || vista.b_vettori {
        let lung = 0.9 * l.w0;
        let passo = 0.12;
        // E and B sampled one cell apart along the axis, then a Catmull–Rom curve between the samples
        // (smooth, no corners from the linear interpolation of the lattice)
        let campioni: Vec<([f64; 3], [f64; 3])> = (0..=n + 2).map(|i| {
            let s0 = i as f64 - 1.0 - (n as f64) / 2.0;
            let q = campo_in(f, f32v([0, 1, 2].map(|c| centro[c] + dir[c] * s0)));
            ([q[0].re as f64, q[1].re as f64, q[2].re as f64], [q[0].im as f64, q[1].im as f64, q[2].im as f64])
        }).collect();
        let cr = |p0: f64, p1: f64, p2: f64, p3: f64, t: f64| 0.5 * (2.0 * p1 + (p2 - p0) * t + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t * t + (3.0 * p1 - p0 - 3.0 * p2 + p3) * t * t * t);
        let mut s = -(n as f64) / 2.0;
        let mut k = 0usize;
        while s < n as f64 / 2.0 - 1.0 {
            let a = [0, 1, 2].map(|c| centro[c] + dir[c] * s);
            let x = s + (n as f64) / 2.0 + 1.0;
            let i = x.floor() as usize;
            let t = x - x.floor();
            let (c0, c1, c2, c3) = (&campioni[i - 1], &campioni[i], &campioni[i + 1], &campioni[(i + 2).min(n + 2)]);
            let e = [0, 1, 2].map(|c| cr(c0.0[c], c1.0[c], c2.0[c], c3.0[c], t));
            let b = [0, 1, 2].map(|c| cr(c0.1[c], c1.1[c], c2.1[c], c3.1[c], t));
            let m = punto(e, e).sqrt();
            let mb = punto(b, b).sqrt();
            let ang = punto(e, e2d).atan2(punto(e, e1d));
            // streaks every λ/4
            if k % ((lam / 4.0 / passo).round().max(1.0) as usize) == 0 {
                if vista.e_vettori && m > 0.08 * emax as f64 {
                    let pc = phase_color(ang, 1.0);
                    tratto(&mut out, a, [e[0] / m, e[1] / m, e[2] / m], lung * m / emax as f64, [0.7 * pc[0] + 0.3, 0.7 * pc[1] + 0.3, 0.7 * pc[2] + 0.3], 0.45, 0.7);
                }
                if vista.b_vettori && mb > 0.08 * emax as f64 {
                    tratto(&mut out, a, [b[0] / mb, b[1] / mb, b[2] / mb], lung * mb / emax as f64, [0.6 * BLU[0] + 0.4, 0.6 * BLU[1] + 0.4, 0.6 * BLU[2] + 0.4], 0.45, 0.7);
                }
            }
            s += passo;
            k += 1;
        }
    }
    // ── Poynting, spin and cross-section ──
    if vista.poynting || vista.spin {
        let g = ((lam / 2.0).round() as usize).max(3);
        let ng = n / g;
        let sin_w = (std::f64::consts::TAU / lam).sin() as f32;
        let parti: Vec<Vec<Ball>> = (0..ng).into_par_iter().map(|gi| {
            let mut v = Vec::new();
            for gj in 0..ng {
                for gk in 0..ng {
                    let (i, j, k) = (gi * g, gj * g, gk * g);
                    let s = (i * n + j) * n + k;
                    let p = &f.f[s];
                    let i_tot = p[0].norm2() + p[1].norm2() + p[2].norm2();
                    if i_tot < 0.2 * imax { continue; }
                    let pos = [i as f64 - h as f64, j as f64 - h as f64, k as f64 - h as f64];
                    let e = [p[0].re, p[1].re, p[2].re];
                    let b = [p[0].im, p[1].im, p[2].im];
                    if vista.poynting {
                        let sv = x3(&e, &b);
                        let m = d3(&sv, &sv).sqrt();
                        if m > 1e-30 {
                            let u = [(sv[0] / m) as f64, (sv[1] / m) as f64, (sv[2] / m) as f64];
                            tratto(&mut v, pos, u, 0.6 * lam * (m / (0.5 * imax)).min(1.5) as f64, ORO, 0.4, 0.6);
                        }
                    }
                    if vista.spin {
                        if let Some(ep) = e_prima {
                            let sd = x3(&ep[s], &e);
                            let mm = d3(&sd, &sd).sqrt();
                            let m = mm / sin_w.abs().max(1e-6);
                            if m > 0.05 * imax {
                                let u = [(sd[0] / mm) as f64, (sd[1] / mm) as f64, (sd[2] / mm) as f64];
                                let col = if d3(&sd, &df) >= 0.0 { [0.45, 1.0, 0.55] } else { [1.0, 0.40, 0.85] };
                                tratto(&mut v, [pos[0] + 0.5, pos[1] + 0.5, pos[2] + 0.5], u, 0.6 * lam * (m / (0.5 * imax)).min(1.5) as f64, col, 0.4, 0.6);
                            }
                        }
                    }
                }
            }
            v
        }).collect();
        for p in parti { out.extend(p); }
    }
    if vista.sezione {
        let cf = f32v(centro);
        let parti: Vec<Vec<Ball>> = (0..n).into_par_iter().map(|i| {
            let mut v = Vec::new();
            for j in 0..n {
                for k in 0..n {
                    let pos = [i as f32 - h, j as f32 - h, k as f32 - h];
                    let dp = [pos[0] - cf[0], pos[1] - cf[1], pos[2] - cf[2]];
                    if d3(&dp, &df).abs() >= 0.5 { continue; }
                    let p = &f.f[(i * n + j) * n + k];
                    let e = [p[0].re, p[1].re, p[2].re];
                    let e2v = d3(&e, &e);
                    if e2v < 0.01 * imax { continue; }
                    let ang = d3(&e, &e2).atan2(d3(&e, &e1)) as f64;
                    v.push((schermo([pos[0] as f64, pos[1] as f64, pos[2] as f64]).to_array(), 1.6, phase_color(ang, 0.05 + 0.6 * (e2v / imax).powf(0.6))));
                }
            }
            v
        }).collect();
        for p in parti { out.extend(p); }
    }
    // ── the box: dim balls along its 12 edges ──
    let lo = -h - 0.5;
    let hi = h - 0.5;
    let steps = 96;
    for e in 0..12 {
        for s in 0..=steps {
            let t = lo + (hi - lo) * s as f32 / steps as f32;
            let (u, v) = ([lo, hi][e & 1], [lo, hi][(e >> 1) & 1]);
            let p = match e / 4 { 0 => [t, u, v], 1 => [u, t, v], _ => [u, v, t] };
            out.push((schermo([p[0] as f64, p[1] as f64, p[2] as f64]).to_array(), 0.6, [0.45, 0.55, 0.75, 0.09]));
        }
    }
    (out, onde, allunga)
}

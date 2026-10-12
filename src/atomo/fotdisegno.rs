//! How the free photon is DRAWN (the physics is in fotone3d.rs; nothing here changes the field).
//!
//! The look follows the photon laboratory (src/photon/render3d.rs): WAVE BALLS with their own phase colour
//! (alpha blending, soft rim) streaming forward, plus a few bright SURFERS with trails.
//!
//! • TRACERS (the wave balls): tens of thousands of soft balls that ride the field. They are born where the photon has energy
//!   (probability ∝ |F|² = E² + c²B²) and every tick they move with the energy velocity of the field at their
//!   point, v = 2 E × cB / (E² + c²B²) (|v| ≤ c; = c in a plane wave), read from the lattice by trilinear
//!   interpolation. Each one carries a short fading tail, so the packet visibly streams forward. Colour = the phase
//!   of the wave at the ball, φ from E·n̂ now and one tick ago (E·n̂ = A cos φ): every colour once per λ, as in
//!   the photon laboratory; opacity ∝ √(|F|²/max).
//! • SURFERS: 6 bright white balls carried by the same energy flow, with long trails. Balls that fall where the photon has left, or reach the end of their life, are born again.
//! • WAVEFRONTS: thin glowing sheets where E along the polarization has a crest (gold) or a trough (cyan)
//!   along the flight; the peak position inside the cell comes from a parabola through three cells, so the sheets
//!   glide smoothly instead of jumping from cell to cell.
//! • HELIX: the tip of E drawn along the beam axis as a continuous glowing tube (overlapping balls every 0.2
//!   cells): circular light = a screw that turns as it flies, linear light = a flat wave.
//! • E and B streaks (off by default): a few short smooth strokes from the axis along E (wheel colours) and B
//!   (blue).
//! • Poynting S = E × B, spin E(t−1) × E(t), cross-section (off by default).

use rayon::prelude::*;

use crate::fotone3d::{base_trasversa, punto, Cf, Fotone3D, Lancio, Polar};
use crate::render::{phase_color, schermo, Ball};
use crate::rng::{hash01, Rng};

pub const N_TRACCIANTI: usize = 32_000;
pub const N_SURF: usize = 6;
pub const L_SCIA: usize = 48;

/// Which layers are drawn.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vista {
    pub traccianti: bool,
    pub fronti: bool,
    pub elica: bool,
    pub e_vettori: bool,
    pub b_vettori: bool,
    pub poynting: bool,
    pub spin: bool,
    pub sezione: bool,
}

impl Default for Vista {
    fn default() -> Self {
        Vista { traccianti: true, fronti: false, elica: true, e_vettori: false, b_vettori: false, poynting: false, spin: false, sezione: false }
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
    pub eta: Vec<u16>,
    pub vita: Vec<u16>,
    /// the surfers: their trails (last = where they are)
    pub surf: Vec<Vec<[f32; 3]>>,
    rng: Rng,
}

impl Traccianti {
    pub fn new(f: &Fotone3D, seme: u64) -> Self {
        let mut t = Traccianti { pos: vec![], vel: vec![], eta: vec![], vita: vec![], surf: vec![], rng: Rng::new(seme) };
        let nuovi = semina(f, N_TRACCIANTI, &mut t.rng);
        for p in nuovi {
            t.pos.push(p);
            t.vel.push([0.0; 3]);
            // random ages, so they do not all fade together
            let v = 30 + (t.rng.next_u64() % 50) as u16;
            t.eta.push((t.rng.next_u64() % v as u64) as u16);
            t.vita.push(v);
        }
        // the surfers start among the brightest balls
        for i in 0..N_SURF { t.surf.push(vec![t.pos[(i * 997) % t.pos.len()]]); }
        t
    }

    /// One tick: every ball moves with the energy velocity of the field (backwards in time: against it);
    /// balls that are dead (old, or where there is no more light) are born again where the light is.
    pub fn avanza(&mut self, f: &Fotone3D, verso: f32, centro: [f64; 3]) {
        let imax = f.f.par_iter().map(|p| p[0].norm2() + p[1].norm2() + p[2].norm2()).reduce(|| 0.0, f32::max).max(1e-30);
        let morti: Vec<usize> = self.pos.par_iter_mut().zip(self.vel.par_iter_mut()).zip(self.eta.par_iter_mut()).zip(self.vita.par_iter())
            .enumerate().filter_map(|(i, (((p, v), eta), vita))| {
                let (vv, u) = vel_energia(f, *p);
                if u < 0.002 * imax || *eta >= *vita { return Some(i); }
                for c in 0..3 { p[c] += verso * vv[c]; }
                *v = vv;
                *eta += 1;
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
            self.eta[i] = 0;
            self.vita[i] = 30 + (self.rng.next_u64() % 50) as u16;
        }
    }
}

/// `k` points drawn with probability ∝ |F|² (a cell, then a random point inside it).
fn semina(f: &Fotone3D, k: usize, rng: &mut Rng) -> Vec<[f32; 3]> {
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
/// Returns (glowing balls, wave balls drawn with alpha blending).
pub fn palline(f: &Fotone3D, l: &Lancio, centro: [f64; 3], tr: &Traccianti, e_prima: Option<&[[f32; 3]]>, vista: Vista) -> (Vec<Ball>, Vec<Ball>) {
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
    let mut out: Vec<Ball> = Vec::with_capacity(4 * N_TRACCIANTI + 40_000);
    let mut onde: Vec<Ball> = Vec::with_capacity(3 * N_TRACCIANTI);
    // ── tracers (wave balls) with their tails, and the surfers ──
    if vista.traccianti {
        let w = std::f64::consts::TAU as f32 / lam as f32;
        let (cw, sw) = (w.cos(), w.sin().max(1e-6));
        let parti: Vec<[Ball; 3]> = tr.pos.par_iter().zip(tr.vel.par_iter()).zip(tr.eta.par_iter().zip(tr.vita.par_iter())).map(|((p, v), (eta, vita))| {
            let q = campo_in(f, *p);
            let e = [q[0].re, q[1].re, q[2].re];
            let b = [q[0].im, q[1].im, q[2].im];
            let u = (d3(&e, &e) + d3(&b, &b)) / imax;
            // the phase of the wave at the ball: E·n̂ = A cos φ now, A cos(φ + ω) one tick ago
            let en = d3(&e, &nh);
            let fase = match e_prima {
                Some(ep) => { let ea = d3(&e_in(ep, n, *p), &nh); ((en * cw - ea) / sw).atan2(en) as f64 }
                None => d3(&e, &e2).atan2(d3(&e, &e1)) as f64,
            };
            let fade = ((*eta as f32 / 6.0).min(1.0)) * (((*vita as f32 - *eta as f32) / 8.0).clamp(0.0, 1.0));
            // two fronts per wavelength (φ = 0 and π, crest and trough): the packet is a stack of moving coloured crests
            let fronte = ((1.0 + (2.0 * fase).cos()) / 2.0).powi(14) as f32;
            let a = fade * (0.15 + 0.6 * u.sqrt()).min(0.75) * (0.02 + 0.98 * fronte);
            let c = phase_color(fase, a);
            let pp = |k: f32| schermo([(p[0] - k * v[0]) as f64, (p[1] - k * v[1]) as f64, (p[2] - k * v[2]) as f64]).to_array();
            [(pp(1.6), 0.6, [c[0], c[1], c[2], 0.25 * a]), (pp(0.8), 0.8, [c[0], c[1], c[2], 0.5 * a]), (pp(0.0), 1.0, c)]
        }).collect();
        for t in parti { onde.extend_from_slice(&t); }
        for s in &tr.surf {
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
    // ── the helix: the tip of E along the axis, a continuous tube ──
    if vista.elica || vista.e_vettori || vista.b_vettori {
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
            if vista.elica && m > 0.05 * emax as f64 {
                let tip = [0, 1, 2].map(|c| a[c] + e[c] / emax as f64 * lung);
                let mut col = phase_color(ang, 0.55 * (m / emax as f64).sqrt() as f32 + 0.1);
                for x in col.iter_mut().take(3) { *x = 0.85 * *x + 0.15; }
                out.push((schermo(tip).to_array(), 2.4, col));
            }
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
    (out, onde)
}

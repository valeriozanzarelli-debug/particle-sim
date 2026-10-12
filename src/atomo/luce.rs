//! LIGHT FROM THE ZIG-ZAG RULE: the photon as spin-1 balls (F = E + icB, 3 complex numbers per site),
//! the same rule as `Light3D` in src/zigzag/walk.rs on a box of any shape: one tick = translate along x, y, z
//! in turn (the +1 part of S_j moves one cell forward, the −1 part one cell back, the 0 part stays), no flip
//! (no mass). Here: the DOUBLE SLIT, the first validation scene of the zig-zag laboratory.
//!
//! Box nx × ny × nz, periodic in z with nz small and the packet uniform along z (the slits are long
//! along z: the problem is 2-D, and so is the textbook formula). A wall at x = x_w with two slits; the wall
//! absorbs (F = 0 on its cells), and a sponge on the x and y edges absorbs what reaches them (so nothing comes
//! back around the periodic box). The screen at x = x_s adds up |E|² tick after tick.

use rayon::prelude::*;

use crate::cplx::C;
use crate::walk::{one_basis, order};

#[derive(Clone)]
pub struct Fenditura {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub f: Vec<[C; 3]>,
    tmp: Vec<[C; 2]>,
    /// 1 = free, 0 = wall; sponge factors on the edges
    maschera: Vec<f32>,
    pub ticks: u64,
    /// wall x, slit centres (y) and width, screen x
    pub x_muro: usize,
    pub fend: [f64; 2],
    pub larghezza: f64,
    pub x_schermo: usize,
    pub lambda: f64,
    /// intensity accumulated on the screen, per y
    pub schermo: Vec<f64>,
}

pub const NX: usize = 176;
pub const NY: usize = 144;
pub const NZ: usize = 2;
pub const LAMBDA: f64 = 8.0;
pub const SPUGNA: usize = 14;

impl Fenditura {
    pub fn new(d: f64, lambda: f64) -> Self {
        let (nx, ny, nz) = (NX, NY, NZ);
        let x_muro = 48;
        let x_schermo = nx - SPUGNA - 6;
        let cy = ny as f64 / 2.0;
        let fend = [cy - d / 2.0, cy + d / 2.0];
        let larghezza = (lambda / 2.0).max(3.0);
        let mut maschera = vec![1.0f32; nx * ny * nz];
        for i in 0..nx {
            for j in 0..ny {
                // sponge: smooth absorption over SPUGNA cells at the x and y edges
                let dx = (SPUGNA as f64 - i.min(nx - 1 - i) as f64).max(0.0) / SPUGNA as f64;
                let dy = (SPUGNA as f64 - j.min(ny - 1 - j) as f64).max(0.0) / SPUGNA as f64;
                let mut m = (1.0 - 0.12 * (dx * dx + dy * dy)) as f32;
                // the wall: 2 cells thick, open only in the slits
                if i == x_muro || i == x_muro + 1 {
                    let y = j as f64 + 0.5;
                    let open = fend.iter().any(|c| (y - c).abs() < larghezza / 2.0);
                    if !open { m = 0.0; }
                }
                for k in 0..nz { maschera[(i * ny + j) * nz + k] = m; }
            }
        }
        let mut s = Fenditura { nx, ny, nz, f: vec![[C::ZERO; 3]; nx * ny * nz], tmp: vec![[C::ZERO; 2]; nx * ny * nz], maschera,
                                ticks: 0, x_muro, fend, larghezza, x_schermo, lambda, schermo: vec![0.0; ny] };
        s.lancia();
        s
    }

    /// A plane-wave packet moving along +x, polarised along z (E = ẑ cos, cB = −ŷ cos), wide in y, 2.5 λ long.
    pub fn lancia(&mut self) {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let k = std::f64::consts::TAU / self.lambda;
        let x0 = 24.0;
        let sx = 1.4 * self.lambda;
        let sy = 1.1 * (self.fend[1] - self.fend[0]);
        let cy = ny as f64 / 2.0;
        self.f.par_iter_mut().enumerate().for_each(|(s, p)| {
            let i = s / (ny * nz);
            let j = (s / nz) % ny;
            let (x, y) = (i as f64 - x0, j as f64 - cy);
            let env = (-x * x / (2.0 * sx * sx) - y * y / (2.0 * sy * sy)).exp();
            let c = (k * x).cos() * env;
            // E along z, cB = k̂ × E = x̂ × ẑ E = −ŷ E
            *p = [C::ZERO, C::new(0.0, -c), C::real(c)];
        });
        let _ = nx;
        self.ticks = 0;
    }

    /// What changes in time (for going back): the field, the clock and the screen.
    pub fn salva(&self) -> (Vec<[C; 3]>, u64, Vec<f64>) { (self.f.clone(), self.ticks, self.schermo.clone()) }
    pub fn ripristina(&mut self, s: &(Vec<[C; 3]>, u64, Vec<f64>)) {
        self.f.clone_from(&s.0);
        self.ticks = s.1;
        self.schermo.clone_from(&s.2);
    }
    /// Is this cell part of the wall?
    pub fn muro(&self, i: usize, j: usize) -> bool { self.maschera[(i * self.ny + j) * self.nz] == 0.0 }

    #[inline]
    fn idx(&self, i: usize, j: usize, k: usize) -> usize { (i * self.ny + j) * self.nz + k }

    fn shift(&mut self, ax: usize) {
        let [u, d] = one_basis(ax);
        self.tmp.par_iter_mut().zip(self.f.par_iter()).for_each(|(t, p)| {
            t[0] = u[0].conj() * p[0] + u[1].conj() * p[1] + u[2].conj() * p[2];
            t[1] = d[0].conj() * p[0] + d[1].conj() * p[1] + d[2].conj() * p[2];
        });
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let tmp = &self.tmp;
        self.f.par_iter_mut().enumerate().for_each(|(s, p)| {
            let (i, j, k) = (s / (ny * nz), (s / nz) % ny, s % nz);
            let (m1, p1) = match ax {
                0 => (((i + nx - 1) % nx * ny + j) * nz + k, ((i + 1) % nx * ny + j) * nz + k),
                1 => ((i * ny + (j + ny - 1) % ny) * nz + k, (i * ny + (j + 1) % ny) * nz + k),
                _ => ((i * ny + j) * nz + (k + nz - 1) % nz, (i * ny + j) * nz + (k + 1) % nz),
            };
            // F' = F + u₊ (c₊ from the cell behind − c₊ here) + u₋ (c₋ from the cell ahead − c₋ here)
            let dp = tmp[m1][0] - tmp[s][0];
            let dm = tmp[p1][1] - tmp[s][1];
            for q in 0..3 { p[q] += u[q] * dp + d[q] * dm; }
        });
    }

    pub fn tick(&mut self) {
        for ax in order(self.ticks) { self.shift(ax); }
        let m = &self.maschera;
        self.f.par_iter_mut().zip(m.par_iter()).for_each(|(p, &w)| if w < 1.0 { for c in p.iter_mut() { *c = c.scale(w as f64); } });
        // the screen collects |E|² (the energy that crosses it, tick after tick)
        let xs = self.x_schermo;
        for j in 0..self.ny {
            let mut e = 0.0;
            for k in 0..self.nz { let p = self.f[self.idx(xs, j, k)]; e += p.iter().map(|c| c.re * c.re).sum::<f64>(); }
            self.schermo[j] += e;
        }
        self.ticks += 1;
    }

    /// Electric energy (Re F)² per site (summed over z), for drawing.
    pub fn energia(&self) -> Vec<f32> {
        let nz = self.nz;
        (0..self.nx * self.ny).into_par_iter().map(|c| {
            (0..nz).map(|k| self.f[c * nz + k].iter().map(|v| (v.re * v.re) as f32).sum::<f32>()).sum()
        }).collect()
    }

    /// First-order bright fringe on the screen: measured distance from the axis (cells, mean of both sides)
    /// and the exact geometric one, where the two paths from the slits differ by one wavelength:
    /// √(L² + (y + d/2)²) − √(L² + (y − d/2)²) = λ.
    pub fn frange(&self) -> Option<(f64, f64)> {
        let cy = self.ny as f64 / 2.0;
        let l = (self.x_schermo - self.x_muro) as f64 - 1.0;
        let d = self.fend[1] - self.fend[0];
        let diff = |y: f64| (l * l + (y + d / 2.0).powi(2)).sqrt() - (l * l + (y - d / 2.0).powi(2)).sqrt() - self.lambda;
        let (mut lo, mut hi) = (0.0, cy);
        for _ in 0..60 { let m = 0.5 * (lo + hi); if diff(m) > 0.0 { hi = m; } else { lo = m; } }
        let y_exp = 0.5 * (lo + hi);
        let mut ys = vec![];
        for sgn in [1.0, -1.0] {
            let (a, b) = ((0.6 * y_exp) as usize, ((1.4 * y_exp).min(cy - SPUGNA as f64 - 1.0)) as usize);
            // peak with a parabola through the three highest cells
            let mut best = (0.0, a);
            for dy in a..b {
                let j = (cy + sgn * (dy as f64 + 0.5) - 0.5).round() as usize;
                if self.schermo[j] > best.0 { best = (self.schermo[j], dy); }
            }
            if best.0 <= 0.0 { return None; }
            let at = |dy: usize| self.schermo[(cy + sgn * (dy as f64 + 0.5) - 0.5).round() as usize];
            let (y0, y1, y2) = (at(best.1 - 1), at(best.1), at(best.1 + 1));
            let den = y0 - 2.0 * y1 + y2;
            let off = if den.abs() > 1e-300 { 0.5 * (y0 - y2) / den } else { 0.0 };
            ys.push(best.1 as f64 + 0.5 + off);
        }
        Some((0.5 * (ys[0] + ys[1]), y_exp))
    }
}

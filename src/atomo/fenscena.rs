//! The double slit with the zig-zag light (luce.rs), run in the frame (~5–10 ms per tick), with exact going
//! back: keyframes every `OGNI` ticks + replay (the wall and the sponge absorb, so the rule cannot simply be
//! undone as for the atom).

use crate::carbonio::Rng;
use crate::cplx::C;
use crate::luce::{Fenditura, SPUGNA};
use crate::render::{schermo, Ball};

pub const OGNI: u64 = 8;
pub const FINE: u64 = 300;
/// Drawn depth of the (infinitely long) slits, cells.
pub const PROF: f64 = 28.0;

pub struct ScenaFenditura {
    pub f: Fenditura,
    chiavi: Vec<(Vec<[C; 3]>, u64, Vec<f64>)>,
    pub lanci: u32,
    pub frange: Option<(f64, f64)>,
}

impl ScenaFenditura {
    pub fn new() -> Self {
        let f = Fenditura::new(36.0, 12.0);
        let k = f.salva();
        ScenaFenditura { f, chiavi: vec![k], lanci: 1, frange: None }
    }
    pub fn avanti(&mut self) {
        if self.f.ticks >= FINE {
            // the next photon: same packet again, the screen keeps counting
            let sch = self.f.schermo.clone();
            self.f.lancia();
            self.f.schermo = sch;
            self.chiavi.clear();
            self.lanci += 1;
        }
        self.f.tick();
        if self.f.ticks % OGNI == 0 { self.chiavi.push(self.f.salva()); }
        if self.f.ticks > 200 { self.frange = self.f.frange(); }
    }
    pub fn indietro(&mut self) {
        if self.f.ticks == 0 { return; }
        let target = self.f.ticks - 1;
        while self.chiavi.last().map(|k| k.1 > target).unwrap_or(false) { self.chiavi.pop(); }
        let Some(k) = self.chiavi.last().cloned() else { return };
        self.f.ripristina(&k);
        while self.f.ticks < target { self.f.tick(); }
    }
    /// Cell (i, j) → scene (box centred on the origin), z drawn over PROF cells.
    fn pos(&self, i: f64, j: f64, z: f64) -> [f64; 3] {
        // physics (x, depth, y): the slit plane faces the camera
        [i - self.f.nx as f64 / 2.0, z, j - self.f.ny as f64 / 2.0]
    }

    pub fn palline(&self) -> Vec<Ball> {
        let f = &self.f;
        let e = f.energia();
        let emax = e.iter().cloned().fold(0.0f32, f32::max).max(1e-20);
        let mut rng = Rng::new(f.ticks * 31 + 7);
        let mut out = Vec::with_capacity(60_000);
        let (nx, ny) = (f.nx, f.ny);
        // light: balls where the electric energy is high (the crests), several per cell spread along z
        for i in SPUGNA / 2..nx - SPUGNA / 2 {
            for j in SPUGNA / 2..ny - SPUGNA / 2 {
                let v = e[i * ny + j] / emax;
                if v < 0.06 { continue; }
                let k = (2.0 + 6.0 * v) as usize;
                for _ in 0..k {
                    let p = self.pos(i as f64 + rng.uni() - 0.5, j as f64 + rng.uni() - 0.5, PROF * (rng.uni() - 0.5));
                    out.push((schermo(p).to_array(), 1.0, [1.0, 0.82 + 0.15 * v, 0.40 + 0.5 * v, 0.10 + 0.35 * v.sqrt()]));
                }
            }
        }
        // the wall: dim blue-grey balls
        let mut rw = Rng::new(5);
        for j in 0..ny {
            if !f.muro(f.x_muro, j) || j < SPUGNA / 2 || j > ny - SPUGNA / 2 { continue; }
            for _ in 0..10 {
                let p = self.pos(f.x_muro as f64 + 2.0 * rw.uni() - 0.5, j as f64 + rw.uni() - 0.5, PROF * (rw.uni() - 0.5));
                out.push((schermo(p).to_array(), 1.3, [0.35, 0.45, 0.65, 0.12]));
            }
        }
        // the screen: clicks ∝ the energy collected so far
        let tot: f64 = f.schermo.iter().sum();
        if tot > 0.0 {
            let n_click = ((self.lanci as usize - 1) * 9000 + (f.ticks as usize * 30)).min(24_000);
            let mut cdf = Vec::with_capacity(ny);
            let mut s = 0.0;
            for v in &f.schermo { s += v; cdf.push(s); }
            let mut rc = Rng::new(11);
            for _ in 0..n_click {
                let u = rc.uni() * s;
                let j = cdf.partition_point(|&c| c < u).min(ny - 1);
                let p = self.pos(f.x_schermo as f64 + 0.6 * rc.uni(), j as f64 + rc.uni() - 0.5, PROF * (rc.uni() - 0.5));
                out.push((schermo(p).to_array(), 0.8, [0.55, 0.95, 1.0, 0.30]));
            }
        }
        out
    }
}

//! THE CORE RENDERER OF THE PROGRAM: GLOWS (Valerio, 13/10: the zig-zag atom's graphics are the base of everything:
//! orbitals, the photon, quarks and gluons, molecules, the droplet are drawn with the same glows, never a separate style).
//!
//! A scene builds a `Livello` (one level of the world: its glows, their tails, its group for the level of detail) and
//! hands it to the drawing (render_zz.rs → istanze.rs/istanze.wgsl). Nothing else: no lines, no meshes, no surfaces.
//!
//!   let mut l = Livello::new("atomo");
//!   l.bagliore(pos, 1.0, [r, g, b], intensita);            // a glow (pos in screen coordinates: render::schermo)
//!   l.scia(&passati, [r, g, b], intensita, peso);           // its short fading tail (newest position first)
//!   l.chiudi(raggio);                                       // the group for the level of detail (grafica::lod)
//!
//! Heads come first in `teste` so the GPU can thin them by their index with the stable key (grafica::chiave); tails are
//! thinned here by the caller with grafica::peso. A level knows its group (`Gruppo`), so a parent level (a molecule over
//! its atoms, a droplet over its molecules) can replace it by its own glows when its screen-space error is below a pixel
//! (grafica::percorri).

use crate::grafica::Gruppo;
use crate::render::Ball;

/// One level of the world drawn as glows.
#[derive(Clone, Debug, Default)]
pub struct Livello {
    pub nome: &'static str,
    /// the glows (additive): position, size factor, colour with alpha = intensity
    pub teste: Vec<Ball>,
    /// the tails: overlapping stretched glows and their stretch vectors (one step each: half-length = one step)
    pub scie: Vec<Ball>,
    pub scie_allunga: Vec<[f32; 3]>,
    /// the group of the level of detail (filled by `chiudi`)
    pub gruppo: Option<Gruppo>,
}

impl Livello {
    pub fn new(nome: &'static str) -> Self { Livello { nome, ..Default::default() } }

    /// A glow at `pos` (screen coordinates), size factor `raggio`, colour, intensity (0…1). Returns its index.
    pub fn bagliore(&mut self, pos: [f32; 3], raggio: f32, colore: [f32; 3], intensita: f32) -> usize {
        self.teste.push((pos, raggio, [colore[0], colore[1], colore[2], intensita]));
        self.teste.len() - 1
    }

    /// The short tail of a glow: `passati` are its last positions, newest first (the first is where the glow is now).
    /// One stretched glow per step, centred on the step and twice as long (it reaches the middle of the steps before and
/// after, plus its own radius), so neighbours overlap by half and the tail
    /// reads as one continuous light; along the tail (age t from 0 to 1) the brightness falls as (1 − t)^1.5 and the size
    /// as 1 − 0.6 t. `peso` multiplies the brightness (the level of detail's soft weight).
    pub fn scia(&mut self, passati: &[[f32; 3]], colore: [f32; 3], intensita: f32, peso: f32) {
        scia_in(&mut self.scie, &mut self.scie_allunga, passati, [colore[0], colore[1], colore[2], intensita], peso);
    }

    /// Close the level: its group for the level of detail (n glows over `raggio`, tails of mean length).
    pub fn chiudi(&mut self, raggio: f64) {
        let k = self.teste.len().max(1);
        let scia = if self.scie.is_empty() { 0.0 } else {
            let tot: f64 = self.scie_allunga.iter().map(|a| ((a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) as f64).sqrt()).sum();
            tot / k as f64
        };
        self.gruppo = Some(Gruppo { n: self.teste.len(), raggio, scia });
    }
}

/// The tail pieces (the same as `Livello::scia`, into given vectors).
pub fn scia_in(out: &mut Vec<Ball>, st: &mut Vec<[f32; 3]>, passati: &[[f32; 3]], c: [f32; 4], peso: f32) {
    let k = passati.len().saturating_sub(1);
    if k == 0 || peso <= 0.0 { return; }
    for j in 0..k {
        let (a, b) = (passati[j], passati[j + 1]);
        let t = (j as f32 + 0.5) / k as f32;
        out.push(([0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])], 1.0 - 0.6 * t,
                  [c[0], c[1], c[2], c[3] * peso * 0.7 * (1.0 - t).powf(1.5)]));
        st.push([a[0] - b[0], a[1] - b[1], a[2] - b[2]]);
    }
}

/// Past positions kept per glow.
pub const S_MAX: usize = 16;
/// The most tail pieces per frame.
pub const MAX_SCIA: f64 = 60_000.0;

/// The recent past of a group of glows, for their tails. `liscia` = 1 keeps the true path (glows riding a current);
/// below 1 the tail follows a smoothed path (a sampler's jagged walk would read as a wire: a declared drawing choice).
#[derive(Clone, Debug)]
pub struct Storie { buf: Vec<[[f32; 3]; S_MAX]>, testa: usize, pieni: usize, liscia: f32 }

impl Storie {
    pub fn new(pos: impl Iterator<Item = [f32; 3]>, liscia: f32) -> Self {
        Storie { buf: pos.map(|p| [p; S_MAX]).collect(), testa: 0, pieni: 1, liscia }
    }
    /// One step: `pos(i)` is where glow i is now.
    pub fn registra(&mut self, pos: impl Iterator<Item = [f32; 3]>) {
        let (t, u) = ((self.testa + 1) % S_MAX, self.testa);
        self.testa = t;
        self.pieni = (self.pieni + 1).min(S_MAX);
        let a = self.liscia;
        for (b, p) in self.buf.iter_mut().zip(pos) {
            let q = b[u];
            b[t] = [q[0] + a * (p[0] - q[0]), q[1] + a * (p[1] - q[1]), q[2] + a * (p[2] - q[2])];
        }
    }
    /// Glow i was moved by hand (reborn elsewhere): no tail across the jump.
    pub fn azzera(&mut self, i: usize, p: [f32; 3]) { self.buf[i] = [p; S_MAX]; }
    /// Forget the past of all the glows (the tails start again).
    pub fn ricomincia(&mut self) { self.pieni = 1; }
    /// The tails of the first glows of `teste` under the level of detail's kept fraction `lod_f`: `scia_s` seconds at `ritmo` steps per second.
    pub fn scie(&self, teste: &[Ball], lod_f: f32, scia_s: f64, ritmo: f64) -> (Vec<Ball>, Vec<[f32; 3]>) {
        let n = self.buf.len().min(teste.len());
        let k = ((scia_s * ritmo).round() as usize).min(S_MAX - 1).min(self.pieni.saturating_sub(1));
        if k == 0 || n == 0 { return (vec![], vec![]); }
        let f = lod_f.clamp(0.0, 1.0);
        let f = f * (MAX_SCIA / ((f as f64 * n as f64).max(1.0) * k as f64)).min(1.0) as f32;
        let cap = (f as f64 * n as f64 * k as f64) as usize + 16;
        let (mut balls, mut st) = (Vec::with_capacity(cap), Vec::with_capacity(cap));
        let mut passati = [[0.0f32; 3]; S_MAX];
        for i in 0..n {
            let w = crate::grafica::peso(i as u32, f);
            if w <= 0.0 { continue; }
            for (j, q) in passati.iter_mut().take(k + 1).enumerate() { *q = self.buf[i][(self.testa + S_MAX - j) % S_MAX]; }
            passati[0] = teste[i].0;
            scia_in(&mut balls, &mut st, &passati[..k + 1], teste[i].2, w);
        }
        (balls, st)
    }
}

#[cfg(test)]
mod prove {
    use super::*;

    /// The tail covers its path without gaps (each piece reaches the next) and fades and shrinks towards its end.
    #[test]
    fn scia_continua() {
        let mut l = Livello::new("prova");
        l.bagliore([0.0; 3], 1.0, [1.0; 3], 0.8);
        let passati: Vec<[f32; 3]> = (0..8).map(|i| [0.1 * i as f32, 0.0, 0.0]).collect();
        l.scia(&passati, [1.0; 3], 0.8, 1.0);
        assert_eq!(l.scie.len(), 7);
        for w in l.scie.windows(2) {
            assert!(w[1].2[3] < w[0].2[3] && w[1].1 < w[0].1);
        }
        for (b, a) in l.scie.iter().zip(&l.scie_allunga) {
            // half-length of the stretched glow = one step: each reaches the centres of its neighbours
            assert!((a[0].abs() - 0.1).abs() < 1e-6 && b.0[1] == 0.0);
        }
        l.chiudi(1.0);
        assert!((l.gruppo.unwrap().scia - 0.7).abs() < 1e-5);
    }
}

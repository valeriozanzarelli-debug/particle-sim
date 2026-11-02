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

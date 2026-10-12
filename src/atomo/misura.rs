//! MEASUREMENT (button "MISURA", key M), the same in every scene.
//!
//! A position measurement gives ONE click, at a point drawn with probability |ψ|² at that instant (Born rule;
//! for the photon: the energy density |F|²). The click flashes (a bright ball and a spreading shell of balls),
//! its position is logged on the HUD, and past clicks stay as a faint scatter: many measurements rebuild the
//! shape of the wave click after click.
//!
//! After the click the user chooses: COLLAPSE (the wave is replaced by a narrow packet at the click and the rule
//! carries on from there; for the photon the detector absorbs it), or just RECORD (a non-destructive view, not
//! a real measurement: the wave is left as it was).
//!
//! The surfers (balls carried by the current, with trails) are de Broglie–Bohm trajectories: an
//! interpretation, not something the rule needs. They are off unless asked for.

use crate::render::{schermo, Ball};

/// What the measurement and the surfers do (set in the panel, the same for every scene).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Opzioni {
    /// show the surfers (Bohm trajectories, an interpretation)
    pub surfisti: bool,
    /// keep the faint scatter of past clicks
    pub mostra_clic: bool,
    /// after a measurement: collapse the wave at the click (else only record)
    pub collassa: bool,
}

impl Default for Opzioni { fn default() -> Self { Opzioni { surfisti: false, mostra_clic: true, collassa: false } } }

/// Steps of the click animation.
pub const DURATA: u32 = 45;

/// The detector: all clicks so far (scene units) and the age of the last one (worker steps).
#[derive(Clone, Debug, Default)]
pub struct Rivelatore {
    pub clic: Vec<[f64; 3]>,
    pub eta: u32,
}

impl Rivelatore {
    pub fn registra(&mut self, p: [f64; 3]) {
        self.clic.push(p);
        self.eta = 0;
    }
    pub fn azzera(&mut self) { self.clic.clear(); self.eta = DURATA; }
    /// One step older; true while the flash is still running (the scene must be redrawn).
    pub fn invecchia(&mut self) -> bool {
        if self.eta < DURATA { self.eta += 1; true } else { false }
    }
    pub fn ultimo(&self) -> Option<[f64; 3]> { self.clic.last().copied() }

    /// Past clicks (faint green dots), the last one (bright) and its flash: a shell of balls that grows to
    /// radius `scala` and fades.
    pub fn disegna(&self, out: &mut Vec<Ball>, scala: f64, mostra_passati: bool) {
        let n = self.clic.len();
        if n == 0 { return; }
        if mostra_passati {
            for p in &self.clic[..n - 1] { out.push((schermo(*p).to_array(), 1.4, [0.40, 1.0, 0.45, 0.65])); }
        }
        let p = self.clic[n - 1];
        out.push((schermo(p).to_array(), 1.6, [0.75, 1.0, 0.75, 0.95]));
        if self.eta < DURATA {
            let x = self.eta as f64 / DURATA as f64;
            let fade = (1.0 - x) as f32;
            out.push((schermo(p).to_array(), 4.0 * (1.0 - 0.6 * x as f32), [0.9, 1.0, 0.9, fade]));
            let r = scala * x.powf(0.6);
            // a shell of balls on a Fibonacci sphere
            let k = 160;
            let oro = std::f64::consts::PI * (3.0 - 5f64.sqrt());
            for i in 0..k {
                let z = 1.0 - 2.0 * (i as f64 + 0.5) / k as f64;
                let rr = (1.0 - z * z).sqrt();
                let ph = oro * i as f64;
                let q = [p[0] + r * rr * ph.cos(), p[1] + r * rr * ph.sin(), p[2] + r * z];
                out.push((schermo(q).to_array(), 0.9, [0.55, 1.0, 0.65, 0.6 * fade]));
            }
        }
    }
}

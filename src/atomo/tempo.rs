//! Reversible time (the same design as src/photon/timeline.rs, written for the atom).
//!
//! The world is deterministic (mondo.rs): the state at frame n is the state at an earlier frame plus the
//! `Input`s applied since. The history keeps
//!   • per frame: the small part of the world (`Piccolo`: nucleus, photons, state, clocks…) and its Input;
//!   • per tick of the rule: the electron positions as f32 (what is drawn), ~100 KB;
//!   • every `OGNI` frames: a full keyframe (all the balls, f64), ~0.5 MB.
//! Going back shows the stored frames exactly as they were. Playing forward from a past frame rebuilds the
//! full state (keyframe + replay of the recorded Inputs, bit for bit) and the recorded future is dropped:
//! the simulation continues live from there (the same actions would give the same future again).
//! Memory is bounded: at most `MAX_FRAMES` frames (~100–150 MB).

use std::collections::VecDeque;

use crate::mondo::{Input, Livelli, Mondo, Piccolo};
use crate::risultati::Risultati;

pub const MAX_FRAMES: usize = 900;
pub const OGNI: usize = 12;
/// Snapshots drawn: one every this many ticks.
pub const PASSO_FOTO: u64 = 3;

/// Electron positions of one tick of the rule (relative to the nucleus), with the number of electrons.
pub struct Istantanea {
    pub tick: u64,
    pub ne: usize,
    pub pos: Vec<[f32; 3]>,
}

pub struct Storia {
    /// absolute index of the first stored frame
    pub inizio: u64,
    pub frames: VecDeque<(Piccolo, Input)>,
    pub chiavi: VecDeque<(u64, Mondo)>,
    pub foto: VecDeque<Istantanea>,
}

pub fn istantanea(m: &Mondo) -> Istantanea {
    let ne = m.dmc.g.ne();
    Istantanea { tick: m.p.tick, ne, pos: m.dmc.x.chunks(3).map(|c| [c[0] as f32, c[1] as f32, c[2] as f32]).collect() }
}

impl Storia {
    pub fn new(m: &Mondo) -> Self {
        let mut s = Storia { inizio: 0, frames: VecDeque::new(), chiavi: VecDeque::new(), foto: VecDeque::new() };
        s.frames.push_back((m.p.clone(), Input::default()));
        s.chiavi.push_back((0, m.clone()));
        s.foto.push_back(istantanea(m));
        s
    }
    pub fn ultimo(&self) -> u64 { self.inizio + self.frames.len() as u64 - 1 }
    pub fn frame(&self, n: u64) -> &Piccolo { &self.frames[(n.clamp(self.inizio, self.ultimo()) - self.inizio) as usize].0 }

    /// Record frame (state after `inp`).
    pub fn registra(&mut self, m: &Mondo, inp: Input) {
        self.frames.push_back((m.p.clone(), inp));
        let n = self.ultimo();
        if inp.ticks > 0 || self.foto.back().map(|f| f.tick) != Some(m.p.tick) { self.foto.push_back(istantanea(m)); }
        if n % OGNI as u64 == 0 { self.chiavi.push_back((n, m.clone())); }
        while self.frames.len() > MAX_FRAMES {
            self.frames.pop_front();
            self.inizio += 1;
        }
        // keep the keyframe that the oldest frame needs
        while self.chiavi.len() > 1 && self.chiavi[1].0 <= self.inizio { self.chiavi.pop_front(); }
        let t0 = self.frames.front().unwrap().0.tick;
        while self.foto.len() > 1 && self.foto[1].tick <= t0.saturating_sub(16) { self.foto.pop_front(); }
    }

    /// The full world at frame n (keyframe + exact replay), and the future after n dropped.
    pub fn riprendi_da(&mut self, n: u64, r: &Risultati, lv: &Livelli) -> Mondo {
        let n = n.clamp(self.inizio, self.ultimo());
        let (k0, mut m) = self.chiavi.iter().rev().find(|(k, _)| *k <= n).map(|(k, m)| (*k, m.clone()))
            .unwrap_or_else(|| (self.chiavi[0].0, self.chiavi[0].1.clone()));
        for f in k0 + 1..=n {
            let inp = self.frames[(f - self.inizio) as usize].1;
            m.passo(&inp, r, lv);
        }
        // drop the future
        let keep = (n - self.inizio + 1) as usize;
        self.frames.truncate(keep);
        while self.chiavi.back().map(|c| c.0 > n).unwrap_or(false) { self.chiavi.pop_back(); }
        let tick = m.p.tick;
        while self.foto.back().map(|f| f.tick > tick).unwrap_or(false) { self.foto.pop_back(); }
        m
    }

    /// The last `k` snapshots up to (and including) the tick of frame n.
    pub fn foto_fino_a(&self, tick: u64, k: usize) -> Vec<&Istantanea> {
        // every PASSO_FOTO-th tick: consecutive ticks of one ball are too close (they would draw little chains)
        let mut v: Vec<&Istantanea> = self.foto.iter().rev().filter(|f| f.tick <= tick && f.tick % PASSO_FOTO == 0).take(k).collect();
        v.reverse();
        v
    }

    pub fn memoria_mb(&self) -> f64 {
        let foto: usize = self.foto.iter().map(|f| f.pos.len() * 12).sum();
        let chiavi: usize = self.chiavi.iter().map(|(_, m)| m.dmc.x.len() * 8 * 2 + m.dmc.n() * 24).sum();
        (foto + chiavi) as f64 / 1e6
    }
}

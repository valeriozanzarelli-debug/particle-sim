//! Reversible time.
//!
//! The simulation is deterministic: every random draw comes from the seeded
//! `StdRng` inside `Sim` (or from a per-grain hash of the seed), and a step of
//! any size is exact. So the state at frame n is fully defined by the state at
//! an earlier frame plus the list of operations applied since.
//!
//! Memory design. A naive snapshot of every frame would cost
//! 60 000 grains × 64 B × 600 frames ≈ 2.3 GB. Instead the history keeps:
//!   • one tiny `Op` per frame (≤ 8 B: launch / step by d / run n photons…);
//!   • a full keyframe (`Sim` clone: grains, pair, statistics, RNG) every
//!     k frames, k chosen so that 600 frames of keyframes fit the budget
//!     `HISTORY_BUDGET_BYTES` (192 MB): k = ⌈600 · size(keyframe) / budget⌉,
//!     e.g. 60 000 grains ≈ 4.3 MB per keyframe → k = 14 → ~43 keyframes.
//!   • a keyframe after every operation that is expensive to replay
//!     (whole photons at high speed), which are small (no photon in flight).
//! Going to frame n = cloning the last keyframe ≤ n and replaying < k ops
//! (exactly, or with consecutive steps merged into one while scrubbing: see
//! `History::reconstruct`; a 60 000-grain step costs ~30 ms in a dev build).
//! The oldest keyframes are dropped when the history exceeds 600 frames or
//! the budget, so memory stays bounded whatever happens.
//!
//! Playing forward from a past frame re-runs the recorded operations from the
//! snapshot (identical result, deterministically). Changing anything in the
//! past (a parameter, emitting a photon) makes that frame the new present:
//! the future is discarded and the simulation continues live from there.

use std::collections::VecDeque;

use crate::config::*;
use crate::physics::{Experiment, Grain, Lab, Sim};

/// One frame's worth of change to the laboratory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Op {
    /// Emit a photon (or a pair).
    Launch,
    /// Light travels this far.
    Step(f32),
    /// Run this many whole photons to completion (high speed factors).
    Photons(u32),
    /// Run the photons in flight to completion.
    Finish,
}

pub fn apply(sim: &mut Sim, op: Op) {
    match op {
        Op::Launch => sim.launch(),
        Op::Step(d) => sim.step(d),
        Op::Photons(n) => (0..n).for_each(|_| sim.run_one()),
        Op::Finish => {
            let mut guard = 0;
            while sim.in_flight() && guard < 10_000 {
                sim.step(40.0);
                guard += 1;
            }
        }
    }
}

/// Approximate heap + inline size of a keyframe.
pub fn keyframe_bytes(sim: &Sim) -> usize {
    let grains: usize = sim.photons.iter().map(|p| p.grains.len()).sum();
    let bins = sim.lab.bins;
    std::mem::size_of::<Sim>()
        + grains * (std::mem::size_of::<Grain>() + 8)   // grains + screen hits (upper bound)
        + sim.photons.len() * bins * 24                  // per-photon screen accumulators
        + sim.stats.clicks.len() * 12
        + bins * 12
}

struct Key {
    frame: u64,
    sim:   Sim,
    bytes: usize,
}

/// Ring buffer of the recent past: ops for every frame, keyframes every k.
pub struct History {
    /// Absolute index of the oldest frame kept (= frame of `keys[0]`).
    base:     u64,
    /// `ops[i]` turns frame `base + i` into frame `base + i + 1`.
    ops:      VecDeque<Op>,
    keys:     VecDeque<Key>,
    /// Frame currently shown.
    cursor:   u64,
    bytes:    usize,
    interval: u64,
}

impl History {
    pub fn new(sim: &Sim) -> Self {
        let mut h = Self { base: 0, ops: VecDeque::new(), keys: VecDeque::new(), cursor: 0, bytes: 0, interval: 1 };
        h.push_key(sim);
        h
    }

    pub fn base(&self) -> u64 { self.base }
    pub fn head(&self) -> u64 { self.base + self.ops.len() as u64 }
    pub fn cursor(&self) -> u64 { self.cursor }
    pub fn at_head(&self) -> bool { self.cursor == self.head() }
    pub fn bytes(&self) -> usize { self.bytes }
    pub fn keyframes(&self) -> usize { self.keys.len() }
    pub fn interval(&self) -> u64 { self.interval }

    fn push_key(&mut self, sim: &Sim) {
        let bytes = keyframe_bytes(sim);
        let k = (HISTORY_FRAMES as f64 * bytes as f64 / HISTORY_BUDGET_BYTES as f64).ceil() as u64;
        self.interval = k.clamp(1, MAX_KEY_INTERVAL);
        self.bytes += bytes;
        self.keys.push_back(Key { frame: self.cursor, sim: sim.clone(), bytes });
    }

    /// Append the op that turned the present into `sim_after` (cursor at head).
    pub fn record(&mut self, op: Op, sim_after: &Sim) {
        debug_assert!(self.at_head());
        self.ops.push_back(op);
        self.cursor += 1;
        let last = self.keys.back().map_or(0, |k| k.frame);
        if matches!(op, Op::Photons(_) | Op::Finish) || self.cursor - last >= self.interval {
            self.push_key(sim_after);
        }
        // Drop the oldest keyframe while the rest still covers the window,
        // or while over budget.
        while self.keys.len() > 1
            && (self.head() - self.keys[1].frame >= HISTORY_FRAMES || self.bytes > HISTORY_BUDGET_BYTES)
            && self.keys[1].frame <= self.cursor
        {
            let old = self.keys.pop_front().unwrap();
            self.bytes -= old.bytes;
            let nb = self.keys[0].frame;
            self.ops.drain(..(nb - self.base) as usize);
            self.base = nb;
        }
    }

    /// The shown frame becomes the present: forget its future and keyframe
    /// it (its state may have just been edited).
    pub fn branch(&mut self, sim: &Sim) {
        self.ops.truncate((self.cursor - self.base) as usize);
        while self.keys.back().is_some_and(|k| k.frame >= self.cursor) {
            let k = self.keys.pop_back().unwrap();
            self.bytes -= k.bytes;
        }
        self.push_key(sim);
    }

    /// Show frame `target` (clamped) exactly, starting from `sim` if it is
    /// the exact state at the cursor: replays every recorded op.
    pub fn move_to(&mut self, sim: &mut Sim, target: u64) {
        let t = target.clamp(self.base, self.head());
        if t == self.cursor { return; }
        let key_between = self.keys.iter().any(|k| k.frame > self.cursor && k.frame <= t);
        if t > self.cursor && !key_between {
            for i in (self.cursor - self.base)..(t - self.base) {
                apply(sim, self.ops[i as usize]);
            }
        } else {
            *sim = self.reconstruct(t, true);
        }
        self.cursor = t;
    }

    /// Show frame `target` (clamped) as a preview: see `reconstruct`.
    pub fn preview(&mut self, sim: &mut Sim, target: u64) {
        let t = target.clamp(self.base, self.head());
        *sim = self.reconstruct(t, false);
        self.cursor = t;
    }

    /// State at frame `t`, rebuilt from the last keyframe ≤ t.
    /// `exact`: replay every recorded op — bit-identical to what was shown
    /// live, but up to k steps. Otherwise (preview, used while scrubbing)
    /// consecutive steps are merged into one: a step of any size is exact
    /// physics and meets the same events with the same random draws in the
    /// same order, so it is the same history up to float rounding, for the
    /// cost of about one step.
    fn reconstruct(&self, t: u64, exact: bool) -> Sim {
        let k = self.keys.iter().rev().find(|k| k.frame <= t).expect("history has a base keyframe");
        let mut sim = k.sim.clone();
        let mut pending = 0.0f32;
        for i in (k.frame - self.base)..(t - self.base) {
            match self.ops[i as usize] {
                Op::Step(d) if !exact => pending += d,
                op => {
                    if pending > 0.0 { sim.step(pending); pending = 0.0; }
                    apply(&mut sim, op);
                }
            }
        }
        if pending > 0.0 { sim.step(pending); }
        sim
    }
}

// ─── Time control ─────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mode { Rewind, Pause, Play, FastForward }

pub struct TimeCtl {
    pub mode:    Mode,
    /// Speed factor: ×1 = a photon crosses the bench in `CROSSING_SECONDS`.
    pub factor:  f32,
    pub history: History,
    /// False while the shown state is a preview (see `History::reconstruct`).
    pub exact:    bool,
    frame_accum:  f32,
    photon_accum: f32,
}

impl TimeCtl {
    pub fn new(sim: &Sim) -> Self {
        Self { mode: Mode::Play, factor: 1.0, history: History::new(sim), exact: true, frame_accum: 0.0, photon_accum: 0.0 }
    }

    /// Speed factor including ⏭.
    pub fn effective_factor(&self) -> f32 {
        if self.mode == Mode::FastForward { self.factor * FAST_FORWARD } else { self.factor }
    }

    /// Signed rate of time: −f rewinding, 0 paused, +f playing.
    pub fn signed_factor(&self) -> f32 {
        match self.mode {
            Mode::Rewind => -self.effective_factor(),
            Mode::Pause => 0.0,
            Mode::Play | Mode::FastForward => self.effective_factor(),
        }
    }

    /// History frames to cross this display frame when scrubbing.
    fn scrub_frames(&mut self) -> u64 {
        self.frame_accum += self.effective_factor().min(MAX_SCRUB_PER_FRAME);
        let n = self.frame_accum.floor();
        self.frame_accum -= n;
        n as u64
    }
}

/// Distance light covers per second at speed ×1: the bench in `CROSSING_SECONDS`.
pub fn base_speed(lab: &Lab) -> f32 {
    let bench = match lab.experiment {
        Experiment::Bell => lab.bell.dist_a.max(lab.bell.dist_b) + 15.0 - lab.source_z,
        _ => lab.screen_z - lab.source_z,
    } + 3.0 * lab.packet_len;
    bench / CROSSING_SECONDS
}

/// Light-clock: (distance in wavelengths, time in femtoseconds).
pub fn light_clock(sim: &Sim) -> (f32, f64) {
    let fs = sim.clock as f64 * NM_PER_UNIT as f64 * 1e-9 / C_M_PER_S * 1e15;
    (sim.clock / sim.lab.lambda, fs)
}

impl crate::LabState {
    /// Make the shown state exact (replaying from its keyframe if needed).
    fn materialize(&mut self) {
        if !self.time.exact {
            let h = &mut self.time.history;
            self.sim = h.reconstruct(h.cursor(), true);
            self.time.exact = true;
        }
    }

    /// Scrub to `target` (a preview while moving; cheap).
    pub fn scrub(&mut self, target: u64) {
        self.time.history.preview(&mut self.sim, target);
        self.time.exact = false;
    }

    /// Apply an op as the new present (discarding any future), and record it.
    pub fn do_op(&mut self, op: Op) {
        self.materialize();
        if !self.time.history.at_head() {
            self.time.history.branch(&self.sim);
        }
        apply(&mut self.sim, op);
        self.time.history.record(op, &self.sim);
    }

    /// The state shown was edited by hand: it becomes the present.
    pub fn edited(&mut self) {
        self.time.exact = true;
        self.time.history.branch(&self.sim);
    }

    /// Advance (or rewind) the laboratory by one display frame of `dt` s.
    pub fn advance(&mut self, dt: f32) {
        let dt = dt.min(MAX_DT);
        match self.time.mode {
            Mode::Pause => {}
            Mode::Rewind => {
                let n = self.time.scrub_frames();
                if n == 0 { return; }
                let h = &self.time.history;
                let target = h.cursor().saturating_sub(n).max(h.base());
                self.scrub(target);
                if target == self.time.history.base() { self.time.mode = Mode::Pause; }
            }
            Mode::Play | Mode::FastForward => {
                if self.time.history.at_head() {
                    self.materialize();
                    self.live(dt);
                    return;
                }
                // In the past: re-run the recorded future.
                let n = self.time.scrub_frames();
                if n == 0 { return; }
                let h = &self.time.history;
                let target = (h.cursor() + n).min(h.head());
                let at_head = target == h.head();
                if n == 1 {
                    self.materialize();
                    self.time.history.move_to(&mut self.sim, target);
                } else {
                    self.scrub(target);
                    if at_head { self.materialize(); }
                }
            }
        }
    }

    fn live(&mut self, dt: f32) {
        let f = self.time.effective_factor();
        if f >= TURBO_FACTOR {
            // Too fast to animate: whole photons per frame.
            if self.sim.in_flight() { self.do_op(Op::Finish); return; }
            if !self.auto { return; }
            self.time.photon_accum = (self.time.photon_accum + f / CROSSING_SECONDS * dt).min(f);
            let want = self.time.photon_accum.floor() as u32;
            if want == 0 { return; }
            let t0 = std::time::Instant::now();
            let mut n = 0;
            while n < want {
                self.sim.run_one();
                n += 1;
                if t0.elapsed().as_millis() > 40 { break; }
            }
            // Whatever did not fit this frame is not carried over forever.
            self.time.photon_accum = (self.time.photon_accum - n as f32).min(1.0);
            self.time.history.record(Op::Photons(n), &self.sim);
            return;
        }
        if !self.sim.in_flight() {
            if self.auto { self.do_op(Op::Launch); }
            return;
        }
        let dist = base_speed(&self.sim.lab) * f * dt;
        if dist > 0.0 { self.do_op(Op::Step(dist)); }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::Lab;

    fn snapshot_positions(sim: &Sim) -> Vec<[f32; 3]> {
        sim.photons.iter().flat_map(|p| p.grains.iter().map(|g| g.pos.to_array())).collect()
    }

    #[test]
    fn rewind_and_replay_are_exact() {
        let mut lab = Lab::preset(Experiment::DoubleSlit);
        lab.grains = 2_000;
        let mut sim = Sim::new(lab);
        let mut h = History::new(&sim);
        let mut trail = Vec::new();
        for i in 0..400 {
            let op = if !sim.in_flight() { Op::Launch } else { Op::Step(0.7 + (i % 5) as f32 * 0.3) };
            apply(&mut sim, op);
            h.record(op, &sim);
            trail.push((h.cursor(), snapshot_positions(&sim), sim.stats.emitted));
        }
        // Scrub backwards through every frame kept, then forward again.
        let mut shown = sim.clone();
        let base = h.base();
        for (frame, pos, emitted) in trail.iter().rev().filter(|t| t.0 >= base) {
            h.move_to(&mut shown, *frame);
            assert_eq!(&snapshot_positions(&shown), pos, "frame {frame}");
            assert_eq!(shown.stats.emitted, *emitted);
        }
        for (frame, pos, _) in trail.iter().filter(|t| t.0 >= base) {
            h.move_to(&mut shown, *frame);
            assert_eq!(&snapshot_positions(&shown), pos, "frame {frame}");
        }
        // Previews (merged steps) are the same history up to rounding.
        for (frame, pos, emitted) in trail.iter().rev().filter(|t| t.0 >= base) {
            h.preview(&mut shown, *frame);
            let p = snapshot_positions(&shown);
            assert_eq!(p.len(), pos.len(), "frame {frame}");
            assert_eq!(shown.stats.emitted, *emitted);
            let worst = p.iter().zip(pos).map(|(a, b)| (0..3).map(|i| (a[i] - b[i]).abs()).fold(0.0, f32::max)).fold(0.0, f32::max);
            assert!(worst < 1e-3, "frame {frame}: preview off by {worst}");
        }
        assert!(h.bytes() <= HISTORY_BUDGET_BYTES);
    }
}

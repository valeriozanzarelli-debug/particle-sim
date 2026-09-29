//! A photon made of grains.
//!
//! A photon is one wave packet, represented by N grains (default tens of
//! thousands). The grains are samples of the photon's wave: where they are
//! dense, the photon is likely to be found.
//!
//! The rule each grain follows — rotate and translate, and nothing else:
//!
//!   • Translate: move in a straight line at the speed of light c.
//!   • Rotate:    its phase arrow turns by 2π every wavelength λ travelled.
//!                The polarization vector E = Re[J·e^{i(φ − ωt)}] therefore
//!                traces a helix for circular light, a flat wave for linear.
//!   • At an aperture (Huygens): the grain becomes a secondary source and
//!                leaves in a new direction, keeping its phase.
//!   • Meeting other grains: grains never push each other (in vacuum light
//!                does not interact with light). They "see" each other only
//!                by comparing phases where they cross — interference:
//!                  I(x) = |Σ a_g|² − Σ|a_g|² = Σ_{g≠h} a_g·a_h*
//!                i.e. the sum over pairs of neighbours of the cosine of
//!                their phase difference.
//!
//! Quantum behaviour belongs to the photon as a whole: at a polarizer or a
//! wall the whole photon either passes or is absorbed (Born rule), and at the
//! screen it is detected at one single point, drawn from its interference
//! pattern. Many photons build the fringes up one dot at a time.
//!
//! Why this is physics and not decoration:
//!   – For free flight the grains are sampled from the beam's Wigner
//!     distribution; for a Gaussian beam it is positive and moves exactly
//!     like straight rays, so the grain density is exactly |ψ|² at all z.
//!   – Behind an aperture the grains are a Monte-Carlo estimate of the
//!     Huygens–Fresnel integral; the pair sum above is its unbiased estimator.

use bevy::math::Vec3;
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::*;

use crate::complex::{inner, jones_norm2, linear, Jones, C};
use crate::entangle::{self, Pair};

pub const TAU: f32 = std::f32::consts::TAU;
const DEG: f32 = std::f32::consts::PI / 180.0;

/// CHSH settings: Alice {0°, 45°}, Bob {22.5°, 67.5°}.
pub const CHSH_A: [f32; 2] = [0.0, 45.0 * DEG];
pub const CHSH_B: [f32; 2] = [22.5 * DEG, 67.5 * DEG];

// ─── Configuration ────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Experiment {
    /// A Gaussian beam flying freely: speed, spread, helicity.
    FreeBeam,
    /// Young's double slit, one photon at a time.
    DoubleSlit,
    /// Up to three absorbing polarizers: Malus's law.
    Polarizers,
    /// A polarization-entangled pair and a Bell (CHSH) test.
    Bell,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Polarization { H, V, D, A, CircPlus, CircMinus }

impl Polarization {
    pub const ALL: [Polarization; 6] = [Self::H, Self::V, Self::D, Self::A, Self::CircPlus, Self::CircMinus];

    pub fn jones(self) -> Jones {
        let s = std::f32::consts::FRAC_1_SQRT_2;
        match self {
            Self::H => [C::ONE, C::ZERO],
            Self::V => [C::ZERO, C::ONE],
            Self::D => [C::real(s), C::real(s)],
            Self::A => [C::real(s), C::real(-s)],
            // E turns from e1 toward e2: right-handed about k → helicity +1
            Self::CircPlus  => [C::real(s), C::new(0.0, s)],
            Self::CircMinus => [C::real(s), C::new(0.0, -s)],
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::H => "linear H (0°)",
            Self::V => "linear V (90°)",
            Self::D => "linear D (45°)",
            Self::A => "linear A (−45°)",
            Self::CircPlus  => "circular, helicity +1",
            Self::CircMinus => "circular, helicity −1",
        }
    }
}

/// Polarizers placed on the two slits to mark which path the photon took.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum WhichPath { Off, Parallel, Orthogonal }

#[derive(Clone, Copy, Debug)]
pub struct PolarizerCfg {
    pub z:       f32,
    pub angle:   f32,
    pub enabled: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct BellCfg {
    /// Alice's and Bob's polarizer angles (manual mode).
    pub angle_a: f32,
    pub angle_b: f32,
    /// Pick each setting at random from the CHSH set at measurement time.
    pub random_chsh: bool,
    /// Replace quantum entanglement with a local hidden-variable model.
    pub hidden_variables: bool,
    /// Distance of Alice's (+z) and Bob's (−z) polarizing beam splitter.
    pub dist_a: f32,
    pub dist_b: f32,
}

#[derive(Clone, Debug)]
pub struct Lab {
    pub experiment:   Experiment,
    /// Wavelength (simulation units).
    pub lambda:       f32,
    /// Grains per photon.
    pub grains:       usize,
    /// Beam waist radius w0 (1/e² of intensity) and where it sits.
    pub waist:        f32,
    pub waist_z:      f32,
    /// Longitudinal σ of the wave packet.
    pub packet_len:   f32,
    /// Packet centre at emission.
    pub source_z:     f32,
    pub polarization: Polarization,

    // Double slit
    pub wall_z:     f32,
    pub slit_sep:   f32,
    pub slit_width: f32,
    /// Half-angle of the Huygens re-emission cone (must cover the screen).
    pub theta_max:  f32,
    pub which_path: WhichPath,

    // Absorbing polarizers (Polarizers experiment)
    pub polarizers: [PolarizerCfg; 3],

    // Detection screen
    pub screen_z:    f32,
    pub screen_half: f32,
    pub bins:        usize,

    pub bell: BellCfg,
    pub seed: u64,
}

impl Lab {
    pub fn preset(experiment: Experiment) -> Self {
        let base = Lab {
            experiment,
            lambda: 1.0,
            grains: 30_000,
            waist: 3.0,
            waist_z: 0.0,
            packet_len: 6.0,
            source_z: -20.0,
            polarization: Polarization::H,
            wall_z: 0.0,
            slit_sep: 8.0,
            slit_width: 2.0,
            theta_max: 0.45,
            which_path: WhichPath::Off,
            polarizers: [
                PolarizerCfg { z: 20.0,  angle: 30.0 * DEG, enabled: true  },
                PolarizerCfg { z: 60.0,  angle: 45.0 * DEG, enabled: false },
                PolarizerCfg { z: 100.0, angle: 90.0 * DEG, enabled: false },
            ],
            screen_z: 160.0,
            screen_half: 60.0,
            bins: 240,
            bell: BellCfg {
                angle_a: 0.0,
                angle_b: 22.5 * DEG,
                random_chsh: true,
                hidden_variables: false,
                dist_a: 40.0,
                dist_b: 70.0,
            },
            seed: 1,
        };
        match experiment {
            Experiment::FreeBeam => Lab {
                grains: 40_000, source_z: 0.0, polarization: Polarization::CircPlus,
                screen_z: 200.0, ..base
            },
            Experiment::DoubleSlit => Lab {
                grains: 60_000, waist: 7.0, waist_z: 0.0, source_z: -40.0,
                polarization: Polarization::D, screen_z: 160.0, ..base
            },
            Experiment::Polarizers => Lab { grains: 20_000, ..base },
            Experiment::Bell => Lab {
                grains: 20_000, source_z: 0.0, packet_len: 4.0, ..base
            },
        }
    }

    #[inline] pub fn k(&self) -> f32 { TAU / self.lambda }
    /// Rayleigh range z_R = π w0² / λ.
    #[inline] pub fn rayleigh(&self) -> f32 { std::f32::consts::PI * self.waist * self.waist / self.lambda }
    /// Theoretical beam radius w(z) = w0 √(1 + (Δz/z_R)²).
    pub fn beam_radius(&self, z: f32) -> f32 {
        let dz = (z - self.waist_z) / self.rayleigh();
        self.waist * (1.0 + dz * dz).sqrt()
    }
    /// Slit centres (x) and half-width.
    pub fn slits(&self) -> ([f32; 2], f32) {
        ([-0.5 * self.slit_sep, 0.5 * self.slit_sep], 0.5 * self.slit_width)
    }
    /// Polarizer on each slit (which-path marker), as Jones projectors.
    pub fn slit_filters(&self) -> [Option<Jones>; 2] {
        match self.which_path {
            WhichPath::Off        => [None, None],
            WhichPath::Parallel   => [Some(linear(0.0)), Some(linear(0.0))],
            WhichPath::Orthogonal => [Some(linear(0.0)), Some(linear(90.0 * DEG))],
        }
    }
    pub fn bin_x(&self, b: usize) -> f32 {
        -self.screen_half + (b as f32 + 0.5) * 2.0 * self.screen_half / self.bins as f32
    }
    fn bin_of(&self, x: f32) -> Option<usize> {
        let f = (x + self.screen_half) / (2.0 * self.screen_half);
        (0.0..1.0).contains(&f).then(|| (f * self.bins as f32) as usize)
    }
}

// ─── Grains and photons ───────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Grain {
    pub pos:    Vec3,
    /// Unit direction of flight.
    pub dir:    Vec3,
    /// Point of the last event (birth or re-emission) and the phase there.
    pub origin: Vec3,
    pub phase0: f32,
    /// Monte-Carlo amplitude weight.
    pub amp:    f32,
    /// Polarization in the transverse basis of `dir`. Zero = undefined (entangled).
    pub jones:  Jones,
    pub alive:  bool,
    /// Passed through a slit as a Huygens secondary source.
    pub diffracted: bool,
    /// Next optical element this grain will meet (elements are met in order).
    pub next: u8,
}

impl Grain {
    /// Phase = phase at the last event + k · path since.
    #[inline]
    pub fn phase(&self, k: f32) -> f32 {
        self.phase0 + k * (self.pos - self.origin).length()
    }

    /// Instantaneous electric field direction·amplitude at light-clock `clock`.
    pub fn e_field(&self, k: f32, clock: f32) -> Vec3 {
        let (e1, e2) = transverse_basis(self.dir);
        let rot = C::expi(self.phase(k) - k * clock);
        e1 * (self.jones[0] * rot).re + e2 * (self.jones[1] * rot).re
    }
}

/// Transverse basis (e1, e2) of a beam direction u: e2 is ŷ made
/// perpendicular to u, e1 = e2 × u. For u = +ẑ: e1 = x̂, e2 = ŷ.
#[inline]
pub fn transverse_basis(u: Vec3) -> (Vec3, Vec3) {
    let e2 = (Vec3::Y - u * u.y).normalize();
    (e2.cross(u), e2)
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Element {
    Wall,
    /// Index into `lab.polarizers`, or Bell side (PBS: both outputs detected).
    Polarizer(usize),
    BellSplitter,
    Screen,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Fate { InFlight, AbsorbedWall, AbsorbedPolarizer(usize), Finished }

#[derive(Clone)]
pub struct Photon {
    pub id:     u64,
    /// 0 = travels +z (Alice), 1 = travels −z (Bob).
    pub side:   usize,
    pub grains: Vec<Grain>,
    elems:      Vec<(f32, Element)>,
    next_elem:  usize,
    pub fate:   Fate,
    /// False while the photon is half of an unmeasured entangled pair.
    pub pol_defined: bool,
    screen:     ScreenAccum,
}

#[derive(Default, Clone)]
struct ScreenAccum {
    sum:   Vec<Jones>,
    self_: Vec<f32>,
    count: Vec<u32>,
    /// (bin, y) of every landing, to place the click vertically.
    hits:  Vec<(u32, f32)>,
}

impl ScreenAccum {
    fn new(bins: usize) -> Self {
        Self { sum: vec![[C::ZERO; 2]; bins], self_: vec![0.0; bins], count: vec![0; bins], hits: Vec::new() }
    }
}

struct Deposit { bin: usize, a: Jones, w2: f32, y: f32 }

// ─── Results ──────────────────────────────────────────────────────────────────

#[derive(Default, Clone)]
pub struct BellStats {
    /// [a index][b index][++, +−, −+, −−] in random CHSH mode.
    pub chsh:   [[[u64; 4]; 2]; 2],
    /// [++, +−, −+, −−] at the manual angles.
    pub manual: [u64; 4],
}

impl BellStats {
    pub fn correlation(c: &[u64; 4]) -> Option<f32> {
        let n = c.iter().sum::<u64>();
        (n > 0).then(|| (c[0] as f32 + c[3] as f32 - c[1] as f32 - c[2] as f32) / n as f32)
    }
    /// S = E(a,b) − E(a,b′) + E(a′,b) + E(a′,b′)
    pub fn chsh_s(&self) -> Option<f32> {
        let e = |i: usize, j: usize| Self::correlation(&self.chsh[i][j]);
        Some(e(0, 0)? - e(0, 1)? + e(1, 0)? + e(1, 1)?)
    }
    /// P(Bob = + | Alice's setting i): must not depend on i (no signalling).
    pub fn bob_plus_given_a(&self, i: usize) -> Option<f32> {
        let (mut plus, mut n) = (0u64, 0u64);
        for j in 0..2 {
            let c = self.chsh[i][j];
            plus += c[0] + c[2];
            n += c.iter().sum::<u64>();
        }
        (n > 0).then(|| plus as f32 / n as f32)
    }
    pub fn pairs(&self) -> u64 {
        self.chsh.iter().flatten().flatten().sum::<u64>() + self.manual.iter().sum::<u64>()
    }
}

#[derive(Clone)]
pub struct Stats {
    pub emitted:        u64,
    pub detected:       u64,
    pub absorbed_wall:  u64,
    pub absorbed_pol:   [u64; 3],
    pub missed:         u64,
    /// Click histogram on the screen.
    pub hist:           Vec<u32>,
    /// Average over photons of each photon's own (normalized) intensity.
    pub mean_intensity: Vec<f64>,
    pub photons_in_mean: u64,
    /// Latest clicks, for drawing.
    pub clicks:         Vec<Vec3>,
    pub bell:           BellStats,
    /// Last Bell outcome pair and settings, for the HUD.
    pub last_bell:      Option<(f32, bool, f32, bool)>,
}

pub const MAX_CLICKS_KEPT: usize = 30_000;

impl Stats {
    fn new(bins: usize) -> Self {
        Self {
            emitted: 0, detected: 0, absorbed_wall: 0, absorbed_pol: [0; 3], missed: 0,
            hist: vec![0; bins], mean_intensity: vec![0.0; bins], photons_in_mean: 0,
            clicks: Vec::new(), bell: BellStats::default(), last_bell: None,
        }
    }
}

// ─── The simulation ───────────────────────────────────────────────────────────

pub struct Sim {
    pub lab:     Lab,
    pub photons: Vec<Photon>,
    pub pair:    Option<Pair>,
    /// Light-clock: distance light has travelled since the current launch.
    pub clock:   f32,
    pub stats:   Stats,
    rng:         StdRng,
    next_id:     u64,
}

/// A full copy of the laboratory state (grains, pair, statistics, RNG), used
/// by the reversible timeline (`timeline.rs`) for its keyframes. Not physics.
impl Clone for Sim {
    fn clone(&self) -> Self {
        Self {
            lab: self.lab.clone(),
            photons: self.photons.clone(),
            pair: self.pair.clone(),
            clock: self.clock,
            stats: self.stats.clone(),
            rng: self.rng.clone(),
            next_id: self.next_id,
        }
    }
}

impl Sim {
    pub fn new(lab: Lab) -> Self {
        let rng = StdRng::seed_from_u64(lab.seed);
        let stats = Stats::new(lab.bins);
        Self { lab, photons: Vec::new(), pair: None, clock: 0.0, stats, rng, next_id: 0 }
    }

    pub fn in_flight(&self) -> bool { !self.photons.is_empty() }

    pub fn reset_stats(&mut self) {
        self.stats = Stats::new(self.lab.bins);
    }

    /// Emit one photon (or one entangled pair).
    pub fn launch(&mut self) {
        self.clock = 0.0;
        let lab = self.lab.clone();
        match lab.experiment {
            Experiment::Bell => {
                let pair = if lab.bell.hidden_variables {
                    Pair::hidden(self.rng.gen_range(0.0..std::f32::consts::PI))
                } else {
                    Pair::phi_plus()
                };
                for side in 0..2 {
                    let jones = pair.local_jones();
                    let defined = pair.is_local();
                    let ph = self.emit(side, jones, defined);
                    self.photons.push(ph);
                }
                self.pair = Some(pair);
            }
            _ => {
                let ph = self.emit(0, lab.polarization.jones(), true);
                self.photons.push(ph);
            }
        }
        self.stats.emitted += 1;
    }

    fn emit(&mut self, side: usize, jones: Jones, pol_defined: bool) -> Photon {
        let lab = &self.lab;
        let sign = if side == 0 { 1.0 } else { -1.0 };
        let k = lab.k();
        let sig_x = 0.5 * lab.waist;                   // |ψ|² ∝ exp(−2r²/w0²)
        let sig_t = lab.lambda / (TAU * lab.waist);    // angular spread λ/(2π w0)
        let (waist_z, center_z) = (sign * lab.waist_z.abs(), lab.source_z);

        let mut grains = Vec::with_capacity(lab.grains);
        for _ in 0..lab.grains {
            let xw = normal(&mut self.rng) * sig_x;
            let yw = normal(&mut self.rng) * sig_x;
            let dir = Vec3::new(normal(&mut self.rng) * sig_t, normal(&mut self.rng) * sig_t, sign).normalize();
            let start_z = center_z + normal(&mut self.rng) * lab.packet_len;
            // Path length from the start point to the waist plane, along the ray
            let s = (waist_z - start_z) / dir.z;
            let pos = Vec3::new(xw, yw, waist_z) - dir * s;
            grains.push(Grain {
                pos, dir, origin: pos,
                phase0: -k * s,              // flat phase (0) on the waist plane
                amp: 1.0, jones, alive: true, diffracted: false, next: 0,
            });
        }

        // Optical elements along this photon's path, in travel order
        let mut elems: Vec<(f32, Element)> = Vec::new();
        match lab.experiment {
            Experiment::FreeBeam => {}
            Experiment::DoubleSlit => elems.push((lab.wall_z, Element::Wall)),
            Experiment::Polarizers => {
                for (i, p) in lab.polarizers.iter().enumerate() {
                    if p.enabled { elems.push((p.z, Element::Polarizer(i))); }
                }
            }
            Experiment::Bell => {
                let d = if side == 0 { lab.bell.dist_a } else { lab.bell.dist_b };
                elems.push((sign * d, Element::BellSplitter));
            }
        }
        let screen_z = match lab.experiment {
            Experiment::Bell => sign * (if side == 0 { lab.bell.dist_a } else { lab.bell.dist_b } + 15.0),
            _ => lab.screen_z,
        };
        elems.push((screen_z, Element::Screen));
        elems.sort_by(|a, b| (sign * a.0).partial_cmp(&(sign * b.0)).unwrap());

        let id = self.next_id;
        self.next_id += 1;
        Photon {
            id, side, grains, elems, next_elem: 0, fate: Fate::InFlight, pol_defined,
            screen: ScreenAccum::new(lab.bins),
        }
    }

    /// Advance every grain of every photon by `dist` (light travels `dist`).
    /// Any step size is exact: crossings are solved analytically.
    pub fn step(&mut self, dist: f32) {
        // 1 ─ photon-level events: the front of a photon reaches an element
        for p in 0..self.photons.len() {
            self.front_events(p, dist);
        }

        // 2 ─ grain-level motion (parallel)
        let lab = &self.lab;
        let ctx = StepCtx::new(lab);
        let seed = lab.seed;
        for ph in self.photons.iter_mut() {
            if ph.fate != Fate::InFlight { continue; }
            let elems = ph.elems.clone();
            let pid = ph.id;
            let deposits: Vec<Deposit> = ph.grains
                .par_iter_mut()
                .enumerate()
                .filter_map(|(gi, g)| {
                    if !g.alive { return None; }
                    advance(g, dist, &elems, &ctx, hash3(seed, pid, gi as u64))
                })
                .collect();
            for d in deposits {
                let s = &mut ph.screen;
                s.sum[d.bin][0] += d.a[0];
                s.sum[d.bin][1] += d.a[1];
                s.self_[d.bin] += d.w2;
                s.count[d.bin] += 1;
                s.hits.push((d.bin as u32, d.y));
            }
            if ph.grains.iter().all(|g| !g.alive) {
                ph.fate = Fate::Finished;
            }
        }
        self.clock += dist;

        // 3 ─ finalize photons that are done
        let mut i = 0;
        while i < self.photons.len() {
            if self.photons[i].fate != Fate::InFlight {
                let ph = self.photons.swap_remove(i);
                self.finalize(ph);
            } else {
                i += 1;
            }
        }
        if self.photons.is_empty() {
            self.pair = None;
        }
    }

    /// Launch a photon and run it to completion with large exact steps.
    pub fn run_one(&mut self) {
        self.launch();
        let big = 40.0;
        let mut guard = 0;
        while self.in_flight() && guard < 10_000 {
            self.step(big);
            guard += 1;
        }
    }

    fn front_events(&mut self, p: usize, dist: f32) {
        loop {
            let ph = &self.photons[p];
            if ph.fate != Fate::InFlight || ph.next_elem >= ph.elems.len() { return; }
            let (z, el) = ph.elems[ph.next_elem];
            let reaches = ph.grains.iter().any(|g| {
                g.alive && { let t = (z - g.pos.z) / g.dir.z; (0.0..=dist).contains(&t) }
            });
            if !reaches { return; }
            match el {
                Element::Wall => self.wall_decision(p),
                Element::Polarizer(i) => self.polarizer(p, i),
                Element::BellSplitter => self.bell_measure(p),
                Element::Screen => {}
            }
            self.photons[p].next_elem += 1;
        }
    }

    /// The photon meets the slit wall: it passes with the probability that
    /// its wave goes through the slits (and their polarizers), else it is
    /// absorbed whole.
    fn wall_decision(&mut self, p: usize) {
        let lab = &self.lab;
        let (centres, half) = lab.slits();
        let filters = lab.slit_filters();
        let ph = &self.photons[p];
        let (mut through, mut total) = (0.0f32, 0.0f32);
        for g in ph.grains.iter().filter(|g| g.alive) {
            let w = jones_norm2(g.jones);
            total += w;
            let t = (lab.wall_z - g.pos.z) / g.dir.z;
            let x = g.pos.x + g.dir.x * t;
            for s in 0..2 {
                if (x - centres[s]).abs() <= half {
                    through += match filters[s] {
                        Some(f) => inner(f, g.jones).norm2(),
                        None => w,
                    };
                }
            }
        }
        let p_pass = if total > 0.0 { through / total } else { 0.0 };
        if self.rng.gen::<f32>() >= p_pass {
            let ph = &mut self.photons[p];
            ph.grains.iter_mut().for_each(|g| g.alive = false);
            ph.fate = Fate::AbsorbedWall;
        }
    }

    /// Absorbing polarizer: pass with probability |⟨a|J⟩|², then J = a.
    fn polarizer(&mut self, p: usize, i: usize) {
        let a = linear(self.lab.polarizers[i].angle);
        let passed = self.measure_single(p, a);
        let ph = &mut self.photons[p];
        if passed {
            ph.grains.iter_mut().for_each(|g| g.jones = a);
        } else {
            ph.grains.iter_mut().for_each(|g| g.alive = false);
            ph.fate = Fate::AbsorbedPolarizer(i);
        }
    }

    fn measure_single(&mut self, p: usize, a: Jones) -> bool {
        let ph = &self.photons[p];
        let (mut num, mut den) = (0.0f32, 0.0f32);
        for g in ph.grains.iter().filter(|g| g.alive) {
            num += inner(a, g.jones).norm2();
            den += jones_norm2(g.jones);
        }
        let prob = if den > 0.0 { num / den } else { 0.0 };
        self.rng.gen::<f32>() < prob
    }

    /// Polarizing beam splitter of the Bell test: both outputs are detected.
    fn bell_measure(&mut self, p: usize) {
        let side = self.photons[p].side;
        let cfg = self.lab.bell;
        let (idx, angle) = if cfg.random_chsh {
            let i = self.rng.gen_range(0..2usize);
            (Some(i), if side == 0 { CHSH_A[i] } else { CHSH_B[i] })
        } else {
            (None, if side == 0 { cfg.angle_a } else { cfg.angle_b })
        };

        let Some(pair) = self.pair.as_mut() else { return };
        let r = self.rng.gen::<f32>();
        let (plus, other) = entangle::measure(pair, side, angle, r);

        // This photon now has a definite polarization
        let mine = if plus { linear(angle) } else { linear(angle + 90.0 * DEG) };
        {
            let ph = &mut self.photons[p];
            ph.pol_defined = true;
            ph.grains.iter_mut().for_each(|g| g.jones = mine);
        }
        // Collapse: the partner's polarization is fixed at this very instant
        if let Some(j) = other {
            if let Some(q) = self.photons.iter_mut().find(|q| q.side != side) {
                q.pol_defined = true;
                q.grains.iter_mut().for_each(|g| g.jones = j);
            }
        }

        let pair = self.pair.as_mut().unwrap();
        pair.record(side, idx, angle, plus);
        if let Some((ia, a, pa, ib, b, pb)) = pair.completed() {
            let o = match (pa, pb) { (true, true) => 0, (true, false) => 1, (false, true) => 2, (false, false) => 3 };
            match (ia, ib) {
                (Some(i), Some(j)) => self.stats.bell.chsh[i][j][o] += 1,
                _ => self.stats.bell.manual[o] += 1,
            }
            self.stats.last_bell = Some((a, pa, b, pb));
        }
    }

    fn finalize(&mut self, ph: Photon) {
        match ph.fate {
            Fate::AbsorbedWall => { self.stats.absorbed_wall += 1; return; }
            Fate::AbsorbedPolarizer(i) => { self.stats.absorbed_pol[i] += 1; return; }
            _ => {}
        }
        if self.lab.experiment == Experiment::Bell { return; }

        // This photon's own intensity on the screen
        let diffracted = self.lab.experiment == Experiment::DoubleSlit;
        let s = &ph.screen;
        let intensity: Vec<f32> = (0..self.lab.bins).map(|b| {
            if diffracted {
                // Interference = sum over pairs of different grains
                (s.sum[b][0].norm2() + s.sum[b][1].norm2() - s.self_[b]).max(0.0)
            } else {
                s.count[b] as f32
            }
        }).collect();
        let total: f32 = intensity.iter().sum();
        if total <= 0.0 {
            self.stats.missed += 1;
            return;
        }
        for (m, i) in self.stats.mean_intensity.iter_mut().zip(&intensity) {
            *m += (*i / total) as f64;
        }
        self.stats.photons_in_mean += 1;

        // Detection: the photon appears at ONE point, drawn from |ψ|²
        let mut r = self.rng.gen::<f32>() * total;
        let mut bin = self.lab.bins - 1;
        for (b, i) in intensity.iter().enumerate() {
            if r < *i { bin = b; break; }
            r -= *i;
        }
        let in_bin: Vec<f32> = s.hits.iter().filter(|h| h.0 as usize == bin).map(|h| h.1).collect();
        let y = if in_bin.is_empty() { 0.0 } else { in_bin[self.rng.gen_range(0..in_bin.len())] };
        let width = 2.0 * self.lab.screen_half / self.lab.bins as f32;
        let x = self.lab.bin_x(bin) + (self.rng.gen::<f32>() - 0.5) * width;

        self.stats.hist[bin] += 1;
        self.stats.detected += 1;
        if self.stats.clicks.len() >= MAX_CLICKS_KEPT { self.stats.clicks.remove(0); }
        self.stats.clicks.push(Vec3::new(x, y, self.lab.screen_z));
    }
}

// ─── Grain motion ─────────────────────────────────────────────────────────────

struct StepCtx {
    k: f32,
    lambda: f32,
    slit_centres: [f32; 2],
    slit_half: f32,
    filters: [Option<Jones>; 2],
    waist: f32,
    theta_max: f32,
    lab: Lab,
}

impl StepCtx {
    fn new(lab: &Lab) -> Self {
        let (slit_centres, slit_half) = lab.slits();
        Self {
            k: lab.k(), lambda: lab.lambda,
            slit_centres, slit_half, filters: lab.slit_filters(),
            waist: lab.waist, theta_max: lab.theta_max, lab: lab.clone(),
        }
    }
}

/// Translate (and rotate the phase) by `dist`, meeting the optical
/// elements strictly in order — no tolerance, so no grain can slip past one.
fn advance(g: &mut Grain, mut dist: f32, elems: &[(f32, Element)], ctx: &StepCtx, rnd: u64) -> Option<Deposit> {
    while g.alive {
        let Some(&(z, el)) = elems.get(g.next as usize) else {
            g.pos += g.dir * dist;
            return None;
        };
        let t = ((z - g.pos.z) / g.dir.z).max(0.0);
        if t > dist {
            g.pos += g.dir * dist;
            return None;
        }
        g.pos += g.dir * t;
        g.pos.z = z;
        dist -= t;
        g.next += 1;
        match el {
            Element::Wall => huygens(g, ctx, rnd),
            Element::Polarizer(_) | Element::BellSplitter => {}   // handled for the whole photon
            Element::Screen => {
                g.alive = false;
                return deposit(g, ctx);
            }
        }
    }
    None
}

/// At the wall: blocked, or re-emitted from the slit as a Huygens wavelet.
fn huygens(g: &mut Grain, ctx: &StepCtx, rnd: u64) {
    let x = g.pos.x;
    let Some(s) = (0..2).find(|&s| (x - ctx.slit_centres[s]).abs() <= ctx.slit_half) else {
        g.alive = false;
        return;
    };
    let phase_here = g.phase(ctx.k);
    if let Some(f) = ctx.filters[s] {
        let c = inner(f, g.jones);
        g.jones = [f[0] * c, f[1] * c];
    }
    // New direction in the x–z plane (the slits are long in y: 2-D problem)
    let u = (rnd >> 11) as f32 / (1u64 << 53) as f32;
    let theta = ctx.theta_max * (2.0 * u - 1.0);
    let sign = g.dir.z.signum();
    g.dir = Vec3::new(theta.sin(), g.dir.y, sign * theta.cos()).normalize();
    // Weight: 1/|U_inc| (grains sample |U|², the integral wants U) × obliquity
    let inc = (-(x * x) / (ctx.waist * ctx.waist)).exp();
    g.amp = 0.5 * (1.0 + theta.cos()) / inc;
    g.origin = g.pos;
    g.phase0 = phase_here - 0.25 * std::f32::consts::PI;   // 2-D wavelet phase
    g.diffracted = true;
}

fn deposit(g: &Grain, ctx: &StepCtx) -> Option<Deposit> {
    let bin = ctx.lab.bin_of(g.pos.x)?;
    if !g.diffracted {
        return Some(Deposit { bin, a: [C::ZERO; 2], w2: 0.0, y: g.pos.y });
    }
    // Huygens kernel e^{ikr}/√(λr), and the Jacobian dθ/dx = cos²θ / L
    let dx = g.pos.x - g.origin.x;
    let dz = g.pos.z - g.origin.z;
    let r = (dx * dx + dz * dz).sqrt();
    let cos_t = dz.abs() / r;
    let w = g.amp * dz.abs() / (cos_t * cos_t * (ctx.lambda * r).sqrt());
    let ph = C::expi(g.phase(ctx.k)) * w;
    let a = [g.jones[0] * ph, g.jones[1] * ph];
    Some(Deposit { bin, a, w2: a[0].norm2() + a[1].norm2(), y: g.pos.y })
}

// ─── Theory (for validation) ──────────────────────────────────────────────────

/// Exact 2-D Huygens–Fresnel integral over the two slits, normalized to sum 1.
pub fn double_slit_reference(lab: &Lab) -> Vec<f32> {
    let (centres, half) = lab.slits();
    let filters = lab.slit_filters();
    let j0 = lab.polarization.jones();
    let k = lab.k();
    let l = lab.screen_z - lab.wall_z;
    let m = 400;
    let mut out: Vec<f32> = (0..lab.bins).map(|b| {
        let x = lab.bin_x(b);
        let mut u = [C::ZERO; 2];
        for s in 0..2 {
            let j = match filters[s] {
                Some(f) => { let c = inner(f, j0); [f[0] * c, f[1] * c] }
                None => j0,
            };
            for i in 0..m {
                let xs = centres[s] - half + (i as f32 + 0.5) * 2.0 * half / m as f32;
                let inc = (-(xs * xs) / (lab.waist * lab.waist)).exp();
                let r = ((x - xs).powi(2) + l * l).sqrt();
                let kern = C::expi(k * r) * (inc * 0.5 * (1.0 + l / r) / (lab.lambda * r).sqrt());
                u[0] += j[0] * kern;
                u[1] += j[1] * kern;
            }
        }
        u[0].norm2() + u[1].norm2()
    }).collect();
    let s: f32 = out.iter().sum();
    out.iter_mut().for_each(|v| *v /= s);
    out
}

/// Probability that one photon crosses all enabled polarizers.
pub fn malus_theory(lab: &Lab) -> f32 {
    let mut j = lab.polarization.jones();
    let mut p = 1.0;
    for pol in lab.polarizers.iter().filter(|p| p.enabled) {
        let a = linear(pol.angle);
        p *= inner(a, j).norm2();
        j = a;
    }
    p
}

/// Fringe spacing of a profile (simulation units): the period P whose
/// cosine resonates most with the central part of the pattern. Robust to
/// the shot noise of a few photons, unlike peak picking.
pub fn fringe_spacing(lab: &Lab, profile: &[f32]) -> Option<f32> {
    let window = 0.8 * lab.screen_half;
    let mut central: Vec<(f32, f32)> = (0..lab.bins)
        .map(|b| (lab.bin_x(b), profile[b]))
        .filter(|(x, _)| x.abs() <= window)
        .collect();
    let mean = central.iter().map(|c| c.1).sum::<f32>() / central.len() as f32;
    central.iter_mut().for_each(|c| c.1 -= mean);
    let spacing = 2.0 * lab.screen_half / lab.bins as f32;
    let (mut best_p, mut best) = (0.0f32, 0.0f32);
    let mut p = 4.0 * spacing;
    while p < window {
        let (mut re, mut im) = (0.0f32, 0.0f32);
        for &(x, v) in &central {
            let a = TAU * x / p;
            re += v * a.cos();
            im += v * a.sin();
        }
        let power = re * re + im * im;
        if power > best { best = power; best_p = p; }
        p *= 1.002;
    }
    (best > 0.0).then_some(best_p)
}

/// Fringe visibility around the centre, (Imax − Imin)/(Imax + Imin),
/// measured as the relative strength of the component at the fringe
/// period over two fringes: I = I0 (1 + V cos(2πx/P + φ)) → V.
/// Unlike max/min it does not read shot noise as fringes.
#[cfg_attr(not(test), allow(dead_code))] // the HUD uses `visibility_at`
pub fn visibility(lab: &Lab, profile: &[f32]) -> f32 {
    let p = fringe_spacing(lab, &double_slit_reference(lab))
        .unwrap_or(lab.lambda * (lab.screen_z - lab.wall_z) / lab.slit_sep);
    visibility_at(lab, profile, p)
}

/// Same, with the fringe period `p` already known.
pub fn visibility_at(lab: &Lab, profile: &[f32], p: f32) -> f32 {
    let (mut re, mut im, mut dc) = (0.0f32, 0.0f32, 0.0f32);
    for b in (0..lab.bins).filter(|&b| lab.bin_x(b).abs() <= p) {
        let x = lab.bin_x(b);
        let a = TAU * x / p;
        re += profile[b] * a.cos();
        im += profile[b] * a.sin();
        dc += profile[b];
    }
    if dc <= 0.0 { return 0.0; }
    (2.0 * (re * re + im * im).sqrt() / dc).min(1.0)
}

/// Pearson correlation of two profiles.
pub fn correlation(a: &[f32], b: &[f32]) -> f32 {
    let n = a.len() as f32;
    let (ma, mb) = (a.iter().sum::<f32>() / n, b.iter().sum::<f32>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma).powi(2);
        sbb += (y - mb).powi(2);
    }
    sab / (saa * sbb).sqrt()
}

// ─── Randomness ───────────────────────────────────────────────────────────────

fn normal(rng: &mut StdRng) -> f32 {
    let u1: f32 = rng.gen_range(f32::EPSILON..1.0);
    let u2: f32 = rng.gen();
    (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
}

/// Deterministic per-grain random bits (reproducible under rayon).
fn hash3(a: u64, b: u64, c: u64) -> u64 {
    let mut x = a ^ b.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ c.wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
    x ^= x >> 30; x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 27; x = x.wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

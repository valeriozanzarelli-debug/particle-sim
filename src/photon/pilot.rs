//! The pilot rule ("two roles"): one LOCAL rule per ball. 2-D (x–z plane).
//!
//! Port of `lab/photon_rule.py` (spec: `docs/REGOLA_FOTONE.md` §3 and §8).
//!
//! A photon = many WAVE BALLS + ONE PHOTON BALL.
//!
//! Wave ball (always the same rule):
//!   • flies straight at c; phase = phase0 + k·|pos − origin|;
//!   • carries origin, direction, amplitude (complex Jones vector) fixed at birth;
//!   • at every object (source, slit wall, polarizer) it is re-emitted: the object
//!     becomes a FAN of new balls, origins on a regular grid (step ds), directions
//!     on a regular grid covering the whole forward half-plane (Θ = π/2 − 10⁻³);
//!     fan amplitude a_i = T · ψ_in(x_i) · ds. Opaque parts: the ball dies.
//!
//! Photon ball (one per photon):
//!   • born at the source at a random wave-ball position weighted by intensity;
//!   • every step it listens, for each origin of the last object crossed, to the
//!     two balls of that origin whose rays straddle it (hat weights) and sums the
//!     exact 2-D Rayleigh–Sommerfeld wavelet G = (ik/2)·Δz·H1(kR)/R and its
//!     gradient (polarization components separate);
//!   • current j = Σ_pol Im(ψ* ∇ψ); it moves at c along j (dx/dz = j_x/j_z, RK4,
//!     dz = min(2λ, max(0.01λ, 0.1·distance from the object)));
//!   • no neighbour → straight on;
//!   • objects: opaque → dies whole; absorbing (slit edges, polarizers) → local
//!     coin p = min(1, Φ_out/Φ_in) read at its own position δ behind the object;
//!     screen → the click is where it lands.
//!
//! No global pass anywhere: no `wall_decision`, no `measure_single`, no click
//! drawn from the whole pattern. Two random numbers per photon at most (birth,
//! one coin per absorbing object), both drawn from a per-photon hash.
//!
//! The wave is the same for every photon of a given laboratory, so the fans
//! (origins and amplitudes) are built once per laboratory (`Bench`) and shared.
//! The wave balls drawn on screen are made from those fans; the photon ball
//! "hears" a ball only if its ray is not blocked between origin and photon
//! (`Heard`), which in a headless run (no balls stored) is always true.

use std::collections::VecDeque;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};
use std::ops::{Add, AddAssign, Mul, Sub};
use std::sync::Arc;

use bevy::math::Vec3;
use rayon::prelude::*;

use crate::complex::{Jones, C};
use crate::physics::{Experiment, Grain, Lab};

// ─── Parameters (docs/REGOLA_FOTONE.md §3.3) ──────────────────────────────────

/// Half-aperture of every fan: the whole forward half-plane.
pub const THETA: f64 = FRAC_PI_2 - 1e-3;
/// Origin spacing, in λ: objects with sharp edges / smooth objects.
pub const DS_SHARP: f64 = 0.02;
pub const DS_SMOOTH: f64 = 0.1;
/// Flux-rule distance δ = 5·ds.
pub const DELTA_PER_DS: f64 = 5.0;
/// Directions per origin, bounds (31 = dθ 0.1, validated as variant H).
pub const NU_MIN: usize = 31;
pub const NU_MAX: usize = 5000;
/// Wave balls leaving |x| > this are dropped (drawing only: see `heard`).
pub const X_LIMIT: f32 = 100.0;
/// A ray is blocked for the photon ball if its ball died (against an object)
/// before reaching the photon's plane z — every object here is a plane ⊥ z.
/// Balls die AT an object, where the photon reads Φ_in: hence a tolerance.
pub const HEARD_TOL_LAMBDA: f64 = 0.01;
/// VISUAL ONLY: half-thickness in y of the wave-ball slab and of the clicks.
/// The rule lives in the x–z plane; this spread only makes the cloud visible.
pub const Y_VISUAL_BALLS: f32 = 1.5;
pub const Y_VISUAL_CLICKS: f32 = 12.0;
/// Points kept in the photon ball's trail, and their spacing (λ).
pub const TRAIL_POINTS: usize = 600;
const TRAIL_STEP: f64 = 0.25;
/// Origins per partial sum (fixed: the summation order never depends on threads).
const CHUNK: usize = 32;

// ─── Complex f64 ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Z {
    pub re: f64,
    pub im: f64,
}

impl Z {
    pub const ZERO: Z = Z { re: 0.0, im: 0.0 };
    #[inline] pub const fn new(re: f64, im: f64) -> Self { Self { re, im } }
    #[inline] pub fn expi(t: f64) -> Self { let (s, c) = t.sin_cos(); Self::new(c, s) }
    #[inline] pub fn conj(self) -> Self { Self::new(self.re, -self.im) }
    #[inline] pub fn norm2(self) -> f64 { self.re * self.re + self.im * self.im }
    #[inline] pub fn abs(self) -> f64 { self.norm2().sqrt() }
    #[inline] pub fn arg(self) -> f64 { self.im.atan2(self.re) }
    #[inline] pub fn s(self, k: f64) -> Self { Self::new(self.re * k, self.im * k) }
    /// i·z
    #[inline] pub fn mul_i(self) -> Self { Self::new(-self.im, self.re) }
    /// −i·z
    #[inline] pub fn mul_neg_i(self) -> Self { Self::new(self.im, -self.re) }
    #[inline] pub fn from_c(c: C) -> Self { Self::new(c.re as f64, c.im as f64) }
}

impl Add for Z { type Output = Z; #[inline] fn add(self, o: Z) -> Z { Z::new(self.re + o.re, self.im + o.im) } }
impl Sub for Z { type Output = Z; #[inline] fn sub(self, o: Z) -> Z { Z::new(self.re - o.re, self.im - o.im) } }
impl AddAssign for Z { #[inline] fn add_assign(&mut self, o: Z) { self.re += o.re; self.im += o.im; } }
impl Mul for Z {
    type Output = Z;
    #[inline] fn mul(self, o: Z) -> Z { Z::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re) }
}

/// Jones vector with f64 complex components.
pub type J2 = [Z; 2];

#[inline] fn j2(j: Jones) -> J2 { [Z::from_c(j[0]), Z::from_c(j[1])] }
/// ⟨a|b⟩
#[inline] fn inner2(a: J2, b: J2) -> Z { a[0].conj() * b[0] + a[1].conj() * b[1] }

// ─── Hankel functions H0⁽¹⁾, H1⁽¹⁾ ─────────────────────────────────────────────

/// Asymptotic-series coefficients a_k(ν) = Π_{m≤k}(4ν² − (2m−1)²) / (k! 8^k), ν = 0, 1.
/// Six terms (the Python lab uses four: rel. error ~7e-6 at u = 12; six: ~2e-7).
const ASY0: [f64; 6] = [1.0, -1.0 / 8.0, 9.0 / 128.0, -225.0 / 3072.0, 11025.0 / 98304.0, -893025.0 / 3932160.0];
const ASY1: [f64; 6] = [1.0, 3.0 / 8.0, -15.0 / 128.0, 315.0 / 3072.0, -14175.0 / 98304.0, 1091475.0 / 3932160.0];

/// Σ a_k (i v)^k (Horner).
#[inline]
fn asy_series(a: &[f64; 6], v: f64) -> Z {
    let mut s = Z::new(a[5], 0.0);
    for k in (0..5).rev() {
        s = s.mul_i().s(v) + Z::new(a[k], 0.0);
    }
    s
}

#[inline]
fn poly(t: f64, c: &[f64; 7]) -> f64 {
    let mut s = c[6];
    for k in (0..6).rev() { s = s * t + c[k]; }
    s
}

/// (H0⁽¹⁾(u), H1⁽¹⁾(u)) for u > 0.
///   u ≥ 12 : asymptotic series (6 terms);
///   3 ≤ u < 12 : Abramowitz–Stegun 9.4.3 / 9.4.6 (modulus–phase form);
///   u < 3  : Abramowitz–Stegun 9.4.1, 9.4.2 / 9.4.4, 9.4.5 (polynomials).
pub fn hankel01(u: f64) -> (Z, Z) {
    if u >= 12.0 {
        let v = 1.0 / u;
        let e = Z::expi(u - FRAC_PI_4).s((2.0 / (PI * u)).sqrt());
        // e^{i(u − 3π/4)} = −i · e^{i(u − π/4)}
        (e * asy_series(&ASY0, v), (e * asy_series(&ASY1, v)).mul_neg_i())
    } else if u >= 3.0 {
        let t = 3.0 / u;
        let r = 1.0 / u.sqrt();
        let f0 = poly(t, &[0.797_884_56, -0.000_000_77, -0.005_527_40, -0.000_095_12, 0.001_372_37, -0.000_728_05, 0.000_144_76]);
        let th0 = u + poly(t, &[-0.785_398_16, -0.041_663_97, -0.000_039_54, 0.002_625_73, -0.000_541_25, -0.000_293_33, 0.000_135_58]);
        let f1 = poly(t, &[0.797_884_56, 0.000_001_56, 0.016_596_67, 0.000_171_05, -0.002_495_11, 0.001_136_53, -0.000_200_33]);
        let th1 = u + poly(t, &[-2.356_194_49, 0.124_996_12, 0.000_056_50, -0.006_378_79, 0.000_743_48, 0.000_798_24, -0.000_291_66]);
        (Z::expi(th0).s(r * f0), Z::expi(th1).s(r * f1))
    } else {
        let t = (u / 3.0) * (u / 3.0);
        let lg = (2.0 / PI) * (0.5 * u).ln();
        let j0 = poly(t, &[1.0, -2.249_999_7, 1.265_620_8, -0.316_386_6, 0.044_447_9, -0.003_944_4, 0.000_210_0]);
        let y0 = lg * j0 + poly(t, &[0.367_466_91, 0.605_593_66, -0.743_503_84, 0.253_001_17, -0.042_612_14, 0.004_279_16, -0.000_248_46]);
        let j1 = u * poly(t, &[0.5, -0.562_499_85, 0.210_935_73, -0.039_542_89, 0.004_433_19, -0.000_317_61, 0.000_011_09]);
        let y1 = (lg * u * j1 + poly(t, &[-0.636_619_8, 0.221_209_1, 2.168_270_9, -1.316_482_7, 0.312_395_1, -0.040_097_6, 0.002_787_3])) / u;
        (Z::new(j0, y0), Z::new(j1, y1))
    }
}

/// Exact 2-D Rayleigh–Sommerfeld wavelet G = (ik/2)·Δz·H1(kR)/R and its
/// gradient, for the offset (dx, dz) from the origin.
#[inline]
pub fn kernel(k: f64, dx: f64, dz: f64) -> (Z, Z, Z) {
    let r = (dx * dx + dz * dz).sqrt();
    let inv = 1.0 / r;
    let (h0, h1) = hankel01(k * r);
    let f = h1.s(inv);
    let fp = (h0.s(k) - h1.s(2.0 * inv)).s(inv); // dF/dR
    let c = 0.5 * k;
    (
        f.s(c * dz).mul_i(),
        fp.s(c * dz * dx * inv).mul_i(),
        (f + fp.s(dz * dz * inv)).s(c).mul_i(),
    )
}

// ─── The local field ──────────────────────────────────────────────────────────

/// ψ and its gradient (both polarization components).
#[derive(Clone, Copy, Default, Debug)]
pub struct Field {
    pub u: J2,
    pub ux: J2,
    pub uz: J2,
}

impl Field {
    #[inline]
    fn add(&mut self, o: &Field) {
        for p in 0..2 { self.u[p] += o.u[p]; self.ux[p] += o.ux[p]; self.uz[p] += o.uz[p]; }
    }
    /// Probability current (per unit k), summed over polarizations.
    pub fn jx(&self) -> f64 { (0..2).map(|p| (self.u[p].conj() * self.ux[p]).im).sum() }
    pub fn jz(&self) -> f64 { (0..2).map(|p| (self.u[p].conj() * self.uz[p]).im).sum() }
    pub fn norm2(&self) -> f64 { self.u[0].norm2() + self.u[1].norm2() }
}

/// "Is ball (fan, origin, direction) heard by a photon ball on the plane z?"
/// (i.e. not dead against an object before z). None = every ray is open
/// (headless run: no balls stored).
pub type Heard<'a> = Option<&'a (dyn Fn(usize, usize, usize, f64) -> bool + Sync)>;

// ─── Fans and the bench ───────────────────────────────────────────────────────

#[derive(Clone, Debug)]
pub enum Obj {
    Source,
    Slits { centres: [f64; 2], half: f64, filters: [Option<J2>; 2] },
    Polarizer { idx: usize, axis: J2 },
}

/// All wave balls emitted by one object.
#[derive(Clone, Debug)]
pub struct Fan {
    pub obj: Obj,
    /// Plane of the object.
    pub ze: f64,
    /// Origins (regular grid) and their amplitudes a_i (Jones, includes ds).
    pub xs: Vec<f64>,
    pub amp: Vec<J2>,
    pub ds: f64,
    /// Directions θ_j = −Θ + (j + ½)·dθ, j < nu, dθ = 2Θ/nu.
    pub nu: usize,
    pub dth: f64,
    /// tan(Θ − dθ): inside this cone both straddling balls exist.
    t_in: f64,
    /// Flux-rule distance behind the object.
    pub delta: f64,
    /// Light-clock when the on-axis front reaches the object (fan emission).
    pub t_emit: f64,
    /// Where these balls end (next object, or the screen).
    pub z_kill: f64,
    /// Index of ball (0, 0) in the photon's ball list.
    pub offset: usize,
}

impl Fan {
    pub fn n_balls(&self) -> usize { self.xs.len() * self.nu }
}

/// The wave of one laboratory: identical for every photon, built once.
#[derive(Clone)]
pub struct Bench {
    pub key: String,
    pub experiment: Experiment,
    pub k: f64,
    pub lambda: f64,
    pub fans: Vec<Fan>,
    pub screen_z: f64,
    pub screen_half: f64,
    pub bins: usize,
    /// Source profile: Gaussian beam (w0, waist z, z_R), polarization.
    src: (f64, f64, f64),
    pol: J2,
    /// Cumulative |a_i|² of the source fan (photon ball's birth).
    src_cdf: Vec<f64>,
    /// Flux j_z of the wave on the screen (bins, normalized): what the
    /// clicks must build up. From the last fan, i.e. the program's own wave.
    pub screen_flux: Vec<f32>,
}

/// 2-D Gaussian beam at (x, z): √(w0/w) e^{−x²/w²} e^{ikx²/2R} e^{−iψ/2}.
fn gauss(k: f64, (w0, zw, zr): (f64, f64, f64), x: f64, z: f64) -> Z {
    let d = z - zw;
    let w = w0 * (1.0 + (d / zr).powi(2)).sqrt();
    let inv_r = d / (d * d + zr * zr);
    let amp = (w0 / w).sqrt() * (-(x * x) / (w * w)).exp();
    Z::expi(0.5 * k * x * x * inv_r - 0.5 * (d / zr).atan()).s(amp)
}

fn grid(lo: f64, hi: f64, ds: f64) -> (Vec<f64>, f64) {
    let n = (((hi - lo) / ds).round() as usize).max(1);
    let step = (hi - lo) / n as f64;
    ((0..n).map(|i| lo + (i as f64 + 0.5) * step).collect(), step)
}

impl Bench {
    pub fn key_of(lab: &Lab) -> String { format!("{lab:?}") }

    pub fn new(lab: &Lab) -> Self {
        let lam = lab.lambda as f64;
        let k = std::f64::consts::TAU / lam;
        let w0 = lab.waist as f64;
        let src = (w0, lab.waist_z as f64, PI * w0 * w0 / lam);
        let pol = j2(lab.polarization.jones());
        let screen_z = lab.screen_z as f64;
        let mut bench = Bench {
            key: Self::key_of(lab), experiment: lab.experiment, k, lambda: lam,
            fans: Vec::new(), screen_z, screen_half: lab.screen_half as f64, bins: lab.bins,
            src, pol, src_cdf: Vec::new(), screen_flux: Vec::new(),
        };

        // Source fan: the Gaussian beam on the source plane
        let zs = lab.source_z as f64;
        let ws = lab.beam_radius(lab.source_z) as f64;
        let ds = DS_SMOOTH * lam;
        let (xs, step) = grid(-3.5 * ws, 3.5 * ws, ds);
        let amp: Vec<J2> = xs.iter().map(|&x| {
            let g = gauss(k, src, x, zs).s(step);
            [g * pol[0], g * pol[1]]
        }).collect();
        bench.push_fan(Obj::Source, zs, xs, amp, step, lab);

        // The objects, in travel order
        let mut objs: Vec<(f64, Obj)> = Vec::new();
        match lab.experiment {
            Experiment::DoubleSlit => {
                let (c, h) = lab.slits();
                let f = lab.slit_filters();
                objs.push((lab.wall_z as f64, Obj::Slits {
                    centres: [c[0] as f64, c[1] as f64], half: h as f64,
                    filters: [f[0].map(j2), f[1].map(j2)],
                }));
            }
            Experiment::Polarizers => {
                for (idx, p) in lab.polarizers.iter().enumerate().filter(|(_, p)| p.enabled) {
                    objs.push((p.z as f64, Obj::Polarizer { idx, axis: j2(crate::complex::linear(p.angle)) }));
                }
            }
            Experiment::FreeBeam | Experiment::Bell | Experiment::Airy | Experiment::Arago => {}
        }
        objs.retain(|(z, _)| *z > zs && *z < screen_z);
        objs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        for (zo, obj) in objs {
            let m = bench.fans.len() - 1;
            let (xs, steps): (Vec<f64>, Vec<f64>) = match &obj {
                Obj::Slits { centres, half, .. } => {
                    let mut xs = Vec::new();
                    let mut st = Vec::new();
                    for c in centres {
                        let (x, s) = grid(c - half, c + half, DS_SHARP * lam);
                        st.extend(std::iter::repeat(s).take(x.len()));
                        xs.extend(x);
                    }
                    (xs, st)
                }
                Obj::Polarizer { .. } => {
                    let half = 4.0 * lab.beam_radius(zo as f32) as f64 + 2.0 * lam;
                    let (x, s) = grid(-half, half, DS_SMOOTH * lam);
                    let n = x.len();
                    (x, vec![s; n])
                }
                Obj::Source => unreachable!(),
            };
            // a_i = T · ψ_in(x_i) · ds, ψ_in read from the incoming fan at x_i
            let psi: Vec<J2> = xs.par_iter().map(|&x| bench.field(m, x, zo, None, false).u).collect();
            let amp: Vec<J2> = xs.iter().zip(&psi).zip(&steps).map(|((&x, &p), &s)| {
                let t = match &obj {
                    Obj::Slits { centres, filters, .. } => {
                        let side = if (x - centres[0]).abs() < (x - centres[1]).abs() { 0 } else { 1 };
                        match filters[side] {
                            Some(f) => { let c = inner2(f, p); [f[0] * c, f[1] * c] }
                            None => p,
                        }
                    }
                    Obj::Polarizer { axis, .. } => { let c = inner2(*axis, p); [axis[0] * c, axis[1] * c] }
                    Obj::Source => p,
                };
                [t[0].s(s), t[1].s(s)]
            }).collect();
            let ds = steps.first().copied().unwrap_or(DS_SMOOTH * lam);
            bench.push_fan(obj, zo, xs, amp, ds, lab);
        }

        // Kill planes, emission times, ball offsets
        let n = bench.fans.len();
        let mut offset = 0;
        for m in 0..n {
            let z_kill = if m + 1 < n { bench.fans[m + 1].ze } else { screen_z };
            let t_emit = if m == 0 { 0.0 } else { bench.fans[m - 1].t_emit + bench.fans[m].ze - bench.fans[m - 1].ze };
            let f = &mut bench.fans[m];
            f.z_kill = z_kill;
            f.t_emit = t_emit;
            f.offset = offset;
            offset += f.n_balls();
        }

        // Birth distribution of the photon ball: |a_i|² of the source fan
        let mut acc = 0.0;
        bench.src_cdf = bench.fans[0].amp.iter().map(|a| { acc += a[0].norm2() + a[1].norm2(); acc }).collect();
        bench.screen_flux = bench.compute_screen_flux(4);
        bench
    }

    fn push_fan(&mut self, obj: Obj, ze: f64, xs: Vec<f64>, amp: Vec<J2>, ds: f64, lab: &Lab) {
        let nu = ((lab.grains as f64 / xs.len() as f64).round() as usize).clamp(NU_MIN, NU_MAX);
        let dth = 2.0 * THETA / nu as f64;
        self.fans.push(Fan {
            obj, ze, xs, amp, ds, nu, dth, t_in: (THETA - dth).tan(),
            delta: DELTA_PER_DS * ds, t_emit: 0.0, z_kill: 0.0, offset: 0,
        });
    }

    pub fn total_balls(&self) -> usize { self.fans.iter().map(|f| f.n_balls()).sum() }

    /// Intensity of the source profile at x (for Φ_in at the source).
    fn source_intensity(&self, x: f64) -> f64 {
        let g = gauss(self.k, self.src, x, self.fans[0].ze);
        g.norm2() * (self.pol[0].norm2() + self.pol[1].norm2())
    }

    fn compute_screen_flux(&self, sub: usize) -> Vec<f32> {
        let m = self.fans.len() - 1;
        let n = self.bins * sub;
        let w = 2.0 * self.screen_half / n as f64;
        let fine: Vec<f64> = (0..n).into_par_iter().map(|i| {
            let x = -self.screen_half + (i as f64 + 0.5) * w;
            self.field(m, x, self.screen_z, None, false).jz().max(0.0)
        }).collect();
        let coarse: Vec<f64> = fine.chunks(sub).map(|c| c.iter().sum()).collect();
        let s: f64 = coarse.iter().sum();
        coarse.iter().map(|v| if s > 0.0 { (v / s) as f32 } else { 0.0 }).collect()
    }

    /// What the photon ball reads at (x, z) from the balls of fan `m`.
    /// Partial sums over fixed chunks of origins, added in order: the result
    /// does not depend on `par` or on the number of threads (exact replay).
    pub fn field(&self, m: usize, x: f64, z: f64, heard: Heard, par: bool) -> Field {
        let fan = &self.fans[m];
        let dz = z - fan.ze;
        if dz <= 0.0 { return Field::default(); }
        let n = fan.xs.len();
        let k = self.k;
        let part = |c: usize| -> Field {
            let mut acc = Field::default();
            for i in c * CHUNK..((c + 1) * CHUNK).min(n) {
                let dx = x - fan.xs[i];
                let w = if heard.is_none() && dx.abs() < fan.t_in * dz { 1.0 } else { hat(fan, m, i, dx, dz, heard) };
                if w == 0.0 { continue; }
                let (g, gx, gz) = kernel(k, dx, dz);
                let a = fan.amp[i];
                let (g, gx, gz) = (g.s(w), gx.s(w), gz.s(w));
                for p in 0..2 {
                    acc.u[p] += g * a[p];
                    acc.ux[p] += gx * a[p];
                    acc.uz[p] += gz * a[p];
                }
            }
            acc
        };
        let nch = n.div_ceil(CHUNK);
        let mut tot = Field::default();
        if par && nch > 1 {
            let parts: Vec<Field> = (0..nch).into_par_iter().map(part).collect();
            for p in &parts { tot.add(p); }
        } else {
            for c in 0..nch { tot.add(&part(c)); }
        }
        tot
    }

    /// dx/dz of the photon ball = j_x/j_z; `own` where nothing is heard.
    pub fn slope(&self, m: usize, x: f64, z: f64, own: f64, heard: Heard, par: bool) -> f64 {
        let f = self.field(m, x, z, heard, par);
        let a2 = f.norm2();
        let jz = f.jz();
        if a2 > 1e-300 && jz > 1e-12 * self.k * a2 { f.jx() / jz } else { own }
    }

    // ── wave balls (for drawing; the photon ball only needs `heard`) ──

    /// The balls of fan `m`, `age` after their emission.
    pub fn emit_balls(&self, m: usize, age: f32, seed: u64) -> Vec<Grain> {
        let fan = &self.fans[m];
        let nu = fan.nu;
        let ze = fan.ze as f32;
        let z_kill = fan.z_kill as f32;
        let units: Vec<(J2, f32, f32, bool)> = fan.amp.iter().map(|a| unit_jones(*a)).collect();
        (0..fan.n_balls()).into_par_iter().map(|b| {
            let (i, j) = (b / nu, b % nu);
            let th = -THETA + (j as f64 + 0.5) * fan.dth;
            let (jones, phase0, mag, live) = units[i];
            let y = Y_VISUAL_BALLS * (2.0 * uniform(hash2(seed ^ 0x5bd1_e995, (m as u64) << 40 | b as u64), 0) as f32 - 1.0);
            let origin = Vec3::new(fan.xs[i] as f32, y, ze);
            let mut g = Grain {
                pos: origin,
                dir: Vec3::new(th.sin() as f32, 0.0, th.cos() as f32),
                origin,
                phase0,
                amp: mag,
                jones: [C::new(jones[0].re as f32, jones[0].im as f32), C::new(jones[1].re as f32, jones[1].im as f32)],
                alive: live,
                diffracted: m > 0,
                next: 0,
            };
            fly(&mut g, age, z_kill);
            g
        }).collect()
    }
}

/// Hat window: the two balls of origin i whose rays straddle the photon.
#[inline]
fn hat(fan: &Fan, m: usize, i: usize, dx: f64, dz: f64, heard: Heard) -> f64 {
    let t = dx.atan2(dz);
    let q = (t + THETA) / fan.dth - 0.5;
    let jf = q.floor();
    let f = q - jf;
    let j = jf as i64;
    let z = fan.ze + dz;
    let ok = |jj: i64| jj >= 0 && (jj as usize) < fan.nu && heard.map_or(true, |h| h(m, i, jj as usize, z));
    match (ok(j), ok(j + 1)) {
        (true, true) => 1.0,
        (a, b) => (if a { 1.0 - f } else { 0.0 }) + (if b { f } else { 0.0 }),
    }
}

/// Unit Jones vector (common phase removed), display phase, |a|, alive.
fn unit_jones(a: J2) -> (J2, f32, f32, bool) {
    let mag = (a[0].norm2() + a[1].norm2()).sqrt();
    if mag == 0.0 { return ([Z::ZERO; 2], 0.0, 0.0, false); }
    let p = if a[0].norm2() >= a[1].norm2() { 0 } else { 1 };
    let ph = a[p].arg();
    let r = Z::expi(-ph).s(1.0 / mag);
    // far-field phase of (ik/2)·H1(kR): kR − π/4
    ([a[0] * r, a[1] * r], (ph - FRAC_PI_4) as f32, mag as f32, true)
}

/// A wave ball flies straight for `d`; it dies at its kill plane or when it
/// leaves the bench sideways (|x| > X_LIMIT).
pub fn fly(g: &mut Grain, d: f32, z_kill: f32) {
    if !g.alive || d <= 0.0 { return; }
    let tz = if g.dir.z > 0.0 { (z_kill - g.pos.z) / g.dir.z } else { f32::INFINITY };
    let tx = if g.dir.x != 0.0 { (X_LIMIT.copysign(g.dir.x) - g.pos.x) / g.dir.x } else { f32::INFINITY };
    let t = tz.min(tx);
    if t <= d {
        g.pos += g.dir * t.max(0.0);
        g.alive = false;
    } else {
        g.pos += g.dir * d;
    }
}

// ─── The photon ball ──────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Outcome {
    /// Landed on the screen at x (inside it).
    Click(f64),
    /// Landed outside the screen.
    Missed(f64),
    /// Absorbed by the wall (opaque part or the flux coin at a slit edge).
    Wall,
    /// Absorbed by polarizer `i` (flux coin).
    Polarizer(usize),
    /// Removed by the flux coin at the source (rare).
    Source,
}

#[derive(Clone)]
pub struct PilotPhoton {
    pub bench: Arc<Bench>,
    pub x: f64,
    pub z: f64,
    /// dx/dz of the last step.
    pub u: f64,
    /// Fan that guides it now (the last object crossed).
    pub fan: usize,
    pub started: bool,
    pub outcome: Option<Outcome>,
    pub trail: VecDeque<Vec3>,
    /// Fans whose wave balls have been emitted (live mode).
    pub emitted: usize,
    /// Last flux coin: (p, survived).
    pub last_coin: Option<(f64, bool)>,
    /// Path already travelled ahead of the light-clock.
    debt: f64,
    key: u64,
    coins: u64,
}

impl PilotPhoton {
    /// Born on the source plane, at a wave-ball position drawn with |a_i|².
    pub fn new(bench: Arc<Bench>, seed: u64, id: u64) -> Self {
        let key = hash2(seed, id);
        let f0 = &bench.fans[0];
        let total = *bench.src_cdf.last().unwrap_or(&1.0);
        let r = uniform(key, 0) * total;
        let i = bench.src_cdf.partition_point(|&c| c <= r).min(f0.xs.len() - 1);
        let x = f0.xs[i] + (uniform(key, 1) - 0.5) * f0.ds;
        let z = f0.ze;
        let mut trail = VecDeque::new();
        trail.push_back(Vec3::new(x as f32, 0.0, z as f32));
        Self {
            bench, x, z, u: 0.0, fan: 0, started: false, outcome: None, trail,
            emitted: 0, last_coin: None, debt: 0.0, key, coins: 0,
        }
    }

    pub fn pos(&self) -> Vec3 { Vec3::new(self.x as f32, 0.0, self.z as f32) }

    /// VISUAL ONLY: y of the click (the rule is 2-D).
    pub fn click_y(&self) -> f32 { Y_VISUAL_CLICKS * (2.0 * uniform(self.key, 999) as f32 - 1.0) }

    /// Move along the current for a path length `dist` (light-clock), or until
    /// z reaches `z_stop`. Returns true when the photon is done.
    pub fn advance(&mut self, dist: f64, heard: Heard, par: bool, z_stop: Option<f64>) -> bool {
        let b = self.bench.clone();
        let lam = b.lambda;
        let mut budget = dist - self.debt;
        self.debt = 0.0;
        if !self.started && budget > 0.0 {
            self.started = true;
            budget -= self.start(&b, heard, par);
        }
        while self.outcome.is_none() && budget > 1e-9 * lam {
            let fan = &b.fans[self.fan];
            let z_next = b.fans.get(self.fan + 1).map_or(b.screen_z, |f| f.ze);
            let mut z_end = z_next;
            if let Some(s) = z_stop {
                if self.z >= s { break; }
                z_end = z_end.min(s);
            }
            let mut dz = (2.0 * lam).min((0.1 * (self.z - fan.ze)).max(0.01 * lam));
            let mut hit = false;
            if dz >= z_end - self.z { dz = z_end - self.z; hit = true; }
            let stretch = (1.0 + self.u * self.u).sqrt();
            if dz * stretch > budget { dz = budget / stretch; hit = false; }
            self.rk4(&b, dz, heard, par);
            self.z = if hit { z_end } else { self.z + dz };
            budget -= dz * (1.0 + self.u * self.u).sqrt();
            self.push_trail(false);
            if hit && z_end == z_next {
                budget -= self.cross(&b, heard, par);
            }
        }
        if budget < 0.0 { self.debt = -budget; }
        self.outcome.is_some()
    }

    fn rk4(&mut self, b: &Bench, dz: f64, heard: Heard, par: bool) {
        let (x, z, m) = (self.x, self.z, self.fan);
        let k1 = b.slope(m, x, z, self.u, heard, par);
        let k2 = b.slope(m, x + 0.5 * dz * k1, z + 0.5 * dz, k1, heard, par);
        let k3 = b.slope(m, x + 0.5 * dz * k2, z + 0.5 * dz, k2, heard, par);
        let k4 = b.slope(m, x + dz * k3, z + dz, k3, heard, par);
        self.u = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        self.x += dz * self.u;
    }

    /// Local coin of the photon ball: survives with probability min(1, r).
    fn coin(&mut self, r: f64) -> bool {
        let p = if r.is_finite() { r.clamp(0.0, 1.0) } else { 0.0 };
        let ok = uniform(self.key, 2 + self.coins) < p;
        self.coins += 1;
        self.last_coin = Some((p, ok));
        ok
    }

    /// Source: straight for δ, then the flux coin Φ_out/Φ_in (Φ_in = k|U_src|²).
    fn start(&mut self, b: &Bench, heard: Heard, par: bool) -> f64 {
        let d = b.fans[0].delta;
        let phi_in = b.k * b.source_intensity(self.x);
        self.x += d * self.u;
        self.z = b.fans[0].ze + d;
        let phi_out = b.field(0, self.x, self.z, heard, par).jz();
        if !self.coin(phi_out / phi_in) { self.outcome = Some(Outcome::Source); }
        self.push_trail(true);
        d * (1.0 + self.u * self.u).sqrt()
    }

    /// The photon ball reaches the next object (or the screen). Returns the
    /// path used by the first step behind the object.
    fn cross(&mut self, b: &Bench, heard: Heard, par: bool) -> f64 {
        let old = self.fan;
        let m = old + 1;
        if m >= b.fans.len() {
            self.outcome = Some(if self.x.abs() <= b.screen_half { Outcome::Click(self.x) } else { Outcome::Missed(self.x) });
            return 0.0;
        }
        let fan = &b.fans[m];
        let (d, zo) = (fan.delta, fan.ze);
        // Φ_in: the incoming balls that reach the object (heard there), their
        // wavelets read at the photon's own position behind it, like Φ_out.
        let at_object = heard.map(|h| move |f: usize, i: usize, j: usize, _z: f64| h(f, i, j, zo));
        let heard_in: Heard = at_object.as_ref().map(|f| f as &(dyn Fn(usize, usize, usize, f64) -> bool + Sync));
        match &fan.obj {
            Obj::Slits { centres, half, .. } => {
                if !centres.iter().any(|c| (self.x - c).abs() <= *half) {
                    self.outcome = Some(Outcome::Wall);     // opaque: absorbed whole
                    return 0.0;
                }
                // straight first step, then the coin
                self.x += d * self.u;
                self.z = zo + d;
                self.fan = m;
                let phi_in = b.field(old, self.x, self.z, heard_in, par).jz();
                let phi_out = b.field(m, self.x, self.z, heard, par).jz();
                if !self.coin(phi_out / phi_in) { self.outcome = Some(Outcome::Wall); }
            }
            Obj::Polarizer { idx, .. } => {
                let idx = *idx;
                // first step along the OUTGOING current (implicit, 2 iterations)
                let x0 = self.x;
                let mut x1 = x0 + d * self.u;
                for _ in 0..2 { x1 = x0 + d * b.slope(m, x1, zo + d, self.u, heard, par); }
                self.u = (x1 - x0) / d;
                self.x = x1;
                self.z = zo + d;
                self.fan = m;
                let phi_in = b.field(old, self.x, self.z, heard_in, par).jz();
                let phi_out = b.field(m, self.x, self.z, heard, par).jz();
                if !self.coin(phi_out / phi_in) { self.outcome = Some(Outcome::Polarizer(idx)); }
            }
            Obj::Source => unreachable!(),
        }
        self.push_trail(true);
        d * (1.0 + self.u * self.u).sqrt()
    }

    fn push_trail(&mut self, force: bool) {
        let p = self.pos();
        let far = self.trail.back().map_or(true, |q| q.distance(p) as f64 >= TRAIL_STEP * self.bench.lambda);
        if force || far {
            self.trail.push_back(p);
            if self.trail.len() > TRAIL_POINTS { self.trail.pop_front(); }
        }
    }
}

#[cfg_attr(not(test), allow(dead_code))]
/// Run one photon to the end without wave balls (tests, high speed factors).
/// Identical to the live run: there every ray is open too (see `heard`).
pub fn run_headless(bench: &Arc<Bench>, seed: u64, id: u64) -> Outcome {
    let mut p = PilotPhoton::new(bench.clone(), seed, id);
    p.advance(f64::INFINITY, None, false, None);
    p.outcome.expect("a headless photon always ends")
}

// ─── Per-photon randomness ────────────────────────────────────────────────────

pub fn hash2(a: u64, b: u64) -> u64 {
    let mut x = a ^ b.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xD1B5_4A32_D192_ED03;
    x ^= x >> 30; x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 27; x = x.wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// Uniform in [0, 1), the n-th number of the stream `key`.
pub fn uniform(key: u64, n: u64) -> f64 {
    (hash2(key, n.wrapping_add(0x632B_E59B_D9B4_E019)) >> 11) as f64 / (1u64 << 53) as f64
}

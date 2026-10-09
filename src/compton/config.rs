//! Constants of the Compton viewer. 1 scene unit = 1 wavelength λ of the INCOMING photon.

// ─── Camera (names shared with the orbit camera module) ───────────────────────
pub const CAM_RADIUS_DEFAULT: f32 = 200.0;
pub const CAM_RADIUS_MIN: f32 = 5.0;
pub const CAM_RADIUS_MAX: f32 = 600.0;

// ─── The shot (scene units, screen seconds) ───────────────────────────────────
/// speed of light on screen: scene units (λ) per screen second
pub const C_SCENE: f32 = 10.0;
/// the incoming packet starts with its centre here (the electron packet is at the origin)
pub const Z_START: f32 = -55.0;
/// length of the photon packet (in λ: it contains this many wave crests)
pub const N_CYCLES: f32 = 6.0;
/// transverse rms width of the drawn photon packet (λ). Real width: `S_GAMMA_LAMBDA` λ.
pub const S_GAMMA_DRAW: f32 = 7.0;
/// rms width of the DRAWN electron packet (λ) — not to scale (its real width is `S_E_M`)
pub const S_E_DRAW: f32 = 1.6;
/// the shot ends when the scattered shell has this radius (λ)
pub const R_END: f32 = 80.0;
/// pause at the end of a shot before the next one (screen s)
pub const PAUSE_END: f32 = 2.5;

// ─── Inputs of the joint wave that are NOT derived from the balls ─────────────
/// transverse rms width of the photon packet, in units of its own wavelength
pub const S_GAMMA_LAMBDA: f64 = 10.0;
/// rms width of the electron packet at rest (m)
pub const S_E_M: f64 = 1.0e-9;

// ─── Balls drawn (before the cut-away view removes the near half) ─────────────
pub const N_IN: usize = 9000;
pub const N_SHELL: usize = 24000;
pub const N_E_REST: usize = 1500;
pub const N_E_RECOIL: usize = 6000;
/// how many correlated (photon direction, electron direction) pairs are joined by faint lines
pub const N_PAIR_LINES: usize = 14;

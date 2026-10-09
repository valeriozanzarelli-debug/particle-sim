//! Constants of the collision viewer. Screen units: the atom is drawn with 1 bohr = `ATOM_SCALE` units;
//! the photon is NOT to scale (its real wavelength, 58 nm = 1100 bohr, is 200 times the atom).

// ─── Camera (names shared with the orbit camera module) ───────────────────────
pub const CAM_RADIUS_DEFAULT: f32 = 92.0;
pub const CAM_RADIUS_MIN: f32 = 4.0;
pub const CAM_RADIUS_MAX: f32 = 300.0;

/// Screen units per bohr for the atom.
pub const ATOM_SCALE: f64 = 5.0;
/// Speed of light on screen (units per second).
pub const C_SCR: f64 = 24.0;
/// The incoming photon starts here (centre of its packet, screen z) and is forgotten past `Z_END`.
pub const Z_START: f64 = -80.0;
pub const Z_END: f64 = 95.0;
/// Drawn wavelength at the resonance (screen units); other wavelengths: × √(λ / 58.4334 nm).
pub const LAMBDA_DRAW_RES: f64 = 12.0;
/// The spontaneous-emission wait is shown slowed down by this factor (0.56 ns → ≈ 2.2 s).
pub const S_DECAY: f64 = 4.0e9;
/// Pause between two photons (s).
pub const PAUSA: f64 = 1.0;
/// Live balls of each cloud (ground 1¹S and excited 2¹P).
pub const N_POOL: usize = 2500;
/// Tick of the live clouds (a.u. of time) and ticks per frame.
pub const TAU_LIVE: f64 = 0.02;
pub const TICKS_PER_FRAME: usize = 1;
/// Wave balls drawn per photon.
pub const N_WAVE: usize = 5000;

//! Fundamental constants for the emergent particle simulation.
//! All parameters are tunable at runtime via the UI.

// ─── Simulation box ────────────────────────────────────────────────────────────
/// Cube half-side in simulation units.  Particles wrap (toroidal boundary).
pub const BOX_HALF: f32 = 40.0;
#[allow(dead_code)]
pub const BOX_SIZE: f32 = BOX_HALF * 2.0;

// ─── Default particle count ────────────────────────────────────────────────────
pub const N_DEFAULT: usize = 100_000;

// ─── PPS (Primordial Particle System) rule parameters ─────────────────────────
/// Base rotation angle per tick (radians). Positive = always turn a little.
pub const ALPHA: f32 = 0.15;
/// Additional rotation per neighbour (radians/neighbour). Drives self-organisation.
pub const BETA: f32 = 0.012;
/// Neighbour detection radius (simulation units).
pub const R_NEIGHBOR: f32 = 3.0;
/// Forward speed (sim units / tick).
pub const SPEED: f32 = 0.35;

// ─── Camera ────────────────────────────────────────────────────────────────────
pub const CAM_RADIUS_DEFAULT: f32 = 110.0;
pub const CAM_RADIUS_MIN: f32     = 1.0;
pub const CAM_RADIUS_MAX: f32     = 500.0;

// ─── Rendering ─────────────────────────────────────────────────────────────────
/// Neighbour count that maps to maximum colour saturation in the density palette.
pub const COLOR_SCALE_N: f32 = 28.0;

//! Constants of the molecule viewer (1 scene unit = 1 bohr = 0.529 Å).

// ─── Camera (names shared with the orbit camera module) ───────────────────────
pub const CAM_RADIUS_DEFAULT: f32 = 14.0;
pub const CAM_RADIUS_MIN: f32 = 1.0;
pub const CAM_RADIUS_MAX: f32 = 120.0;

// ─── Balls shown live ─────────────────────────────────────────────────────────
/// Target population of balls in the viewer (the measurements use 10 000 – 20 000).
pub const BALLS_LIVE: usize = 8000;
/// Time budget per frame for running the rule (seconds).
pub const STEP_BUDGET_S: f64 = 0.006;
/// Ticks of the rule per frame, at most.
pub const MAX_STEPS_PER_FRAME: usize = 20;
/// One vibration period of the nuclei lasts this many seconds on screen (slowed down ~10¹⁴ times).
pub const VIB_PERIOD_S: f64 = 2.5;

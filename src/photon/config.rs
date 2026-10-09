//! Constants of the photon laboratory.

// ─── Camera (names shared with the orbit camera module) ───────────────────────
pub const CAM_RADIUS_DEFAULT: f32 = 230.0;
pub const CAM_RADIUS_MIN: f32     = 1.0;
pub const CAM_RADIUS_MAX: f32     = 900.0;

// ─── Rendering ────────────────────────────────────────────────────────────────
/// Upper bound on E-field arrows drawn per frame.
pub const MAX_ARROWS: usize = 1500;
/// Arrow length per unit of |E|.
pub const ARROW_LEN: f32 = 2.5;
/// Grains drawn at most (all are simulated; see `render::update_meshes`).
pub const MAX_GRAINS_DRAWN: usize = 250_000;
/// Base 3-D rule: wave balls drawn per fan (random subset of the rule's balls).
pub const BALLS_PER_FAN: usize = 30_000;
/// Clicks drawn on the screen.
pub const MAX_CLICKS_DRAWN: usize = 20_000;

// ─── Physical scale shown in the HUD ──────────────────────────────────────────
/// One simulation wavelength unit corresponds to this many nanometres.
pub const NM_PER_UNIT: f32 = 500.0;

// ─── Time (see `timeline.rs`) ─────────────────────────────────────────────────
/// At speed ×1 a photon crosses the bench in about this many seconds.
pub const CROSSING_SECONDS: f32 = 4.0;
/// Range of the speed factor shown in the time bar.
pub const SPEED_MIN: f32 = 0.1;
pub const SPEED_MAX: f32 = 1000.0;
/// From this factor on, whole photons are run per frame without animation.
pub const TURBO_FACTOR: f32 = 50.0;
/// ⏭ multiplies the speed factor by this much.
pub const FAST_FORWARD: f32 = 10.0;
/// Longest frame step (s): a hiccup must not teleport the photon.
pub const MAX_DT: f32 = 1.0 / 20.0;
/// Frames of history kept for rewinding (at least this many, memory allowing).
pub const HISTORY_FRAMES: u64 = 600;
/// Memory budget of the keyframes of the history.
pub const HISTORY_BUDGET_BYTES: usize = 192 << 20;
/// Keyframes are never further apart than this (bounds the replay cost).
pub const MAX_KEY_INTERVAL: u64 = 30;
/// History frames crossed per display frame, at most, while scrubbing.
pub const MAX_SCRUB_PER_FRAME: f32 = 30.0;
/// Speed of light (m/s), for the light-clock readout.
pub const C_M_PER_S: f64 = 299_792_458.0;

// ─── Environment (see `haze.rs`) ──────────────────────────────────────────────
pub const T_DEFAULT: f32 = 293.0;
pub const T_MAX: f32 = 3000.0;
/// Haze points drawn at T_MAX; the count scales as T³ (photon density).
pub const HAZE_MAX_POINTS: usize = 6000;
/// Wien's displacement constant, µm·K.
pub const WIEN_UM_K: f32 = 2897.77;

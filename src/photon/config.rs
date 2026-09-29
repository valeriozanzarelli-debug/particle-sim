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
/// Clicks drawn on the screen.
pub const MAX_CLICKS_DRAWN: usize = 20_000;

// ─── Physical scale shown in the HUD ──────────────────────────────────────────
/// One simulation wavelength unit corresponds to this many nanometres.
pub const NM_PER_UNIT: f32 = 500.0;

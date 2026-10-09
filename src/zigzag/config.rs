//! Constants of the zig-zag viewer.

// ─── Camera (names shared with the orbit camera module) ───────────────────────
pub const CAM_RADIUS_DEFAULT: f32 = 120.0;
pub const CAM_RADIUS_MIN: f32 = 4.0;
pub const CAM_RADIUS_MAX: f32 = 600.0;

// ─── (a) light ────────────────────────────────────────────────────────────────
/// cells per side of the light box (1 cell = 1 scene unit)
pub const N_LIGHT: usize = 64;
/// wavelength of the light packet in cells
pub const LAMBDA_CELLS: f64 = 8.0;
/// width of the packet envelope across the motion (cells; along the motion 1.4×)
pub const SIGMA_LIGHT: f64 = 9.0;

// ─── (b) free electron ────────────────────────────────────────────────────────
pub const N_ELECTRON: usize = 80;
/// lattice spacing in λ_C = ħ/mc: the flip angle θ_m = a/λ_C (the electron mass, INPUT)
pub const A_ELECTRON: f64 = 0.25;
/// rms width of the packet (λ_C)
pub const SIGMA_ELECTRON: f64 = 2.5;
/// number of surfers
pub const N_SURFERS: usize = 14;
/// 1-D checkerboard picture: sites × ticks
pub const CHECK_W: usize = 300;
pub const CHECK_T: usize = 170;

// ─── (c) hydrogen ─────────────────────────────────────────────────────────────
/// scene units per Bohr radius
pub const A0_SCENE: f32 = 6.0;
/// balls of the cloud
pub const N_CLOUD: usize = 26_000;
pub const N_H_SURFERS: usize = 20;
/// physical seconds per screen second (default) for the atom: 2·10⁻¹⁷ s
pub const T_VIS_ATOM: f64 = 2.0e-17;

// ─── balls drawn by the live walks ────────────────────────────────────────────
pub const MAX_BALLS: usize = 40_000;

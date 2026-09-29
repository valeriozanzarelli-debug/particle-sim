//! Global state of the primordial particle simulation.

use bevy::prelude::*;
use rand::Rng;
use crate::config::*;

// ─── Particle world resource ───────────────────────────────────────────────────

#[derive(Resource)]
pub struct ParticleWorld {
    /// 3D positions, toroidal box [−BOX_HALF, +BOX_HALF]³.
    pub positions: Vec<Vec3>,
    /// Unit heading vectors (the direction each particle is moving).
    pub directions: Vec<Vec3>,
    /// Local neighbour count from the last step — used for colour.
    pub neighbor_counts: Vec<u8>,
    /// Number of particles.
    pub n: usize,

    // ── Tunable rule parameters ──────────────────────────────────────────────
    /// Base rotation per tick (radians).
    pub alpha: f32,
    /// Neighbour influence on rotation (radians / neighbour).
    pub beta: f32,
    /// Neighbour detection radius.
    pub r_neighbor: f32,
    /// Forward speed (sim units / tick).
    pub speed: f32,
    /// Half-side of the toroidal simulation cube.
    pub box_half: f32,

    // ── Timing ───────────────────────────────────────────────────────────────
    /// Total ticks simulated so far.
    pub ticks: u64,
    /// Number of simulation steps executed per rendered frame.
    pub steps_per_frame: u32,

    // ── Statistics ───────────────────────────────────────────────────────────
    pub mean_neighbors: f32,
    pub max_cluster: u8,

    /// Incremented each time `reset()` is called. The GPU layer watches this
    /// to know when it must re-upload CPU positions/directions.
    pub gpu_generation: u64,
}

impl ParticleWorld {
    pub fn new(n: usize) -> Self {
        let (positions, directions) = init_grid(n, BOX_HALF);
        let neighbor_counts = vec![0u8; n];
        Self {
            positions,
            directions,
            neighbor_counts,
            n,
            alpha:          ALPHA,
            beta:           BETA,
            r_neighbor:     R_NEIGHBOR,
            speed:          SPEED,
            box_half:       BOX_HALF,
            ticks:          0,
            steps_per_frame: 1,
            mean_neighbors: 0.0,
            max_cluster:    0,
            gpu_generation: 0,
        }
    }

    /// Reinitialise with a new particle count (cubic grid + random headings).
    pub fn reset(&mut self, n: usize) {
        let (positions, directions) = init_grid(n, self.box_half);
        self.positions       = positions;
        self.directions      = directions;
        self.neighbor_counts = vec![0u8; n];
        self.n               = n;
        self.ticks           = 0;
        self.mean_neighbors  = 0.0;
        self.max_cluster     = 0;
        self.gpu_generation += 1;
    }
}

impl Default for ParticleWorld {
    fn default() -> Self {
        Self::new(N_DEFAULT)
    }
}

// ─── Initialisation ────────────────────────────────────────────────────────────

/// Place N particles on a cubic grid with spacing ≈ 1 sim-unit inside the box.
/// Each particle gets a uniformly random unit-sphere heading.
fn init_grid(n: usize, box_half: f32) -> (Vec<Vec3>, Vec<Vec3>) {
    let mut rng = rand::thread_rng();
    let side = (n as f32).cbrt().ceil() as usize;
    let spacing = (2.0 * box_half * 0.85) / side as f32; // fill 85% of box
    let offset = -(side as f32) * spacing / 2.0;

    let mut positions  = Vec::with_capacity(n);
    let mut directions = Vec::with_capacity(n);

    'outer: for gz in 0..side {
        for gy in 0..side {
            for gx in 0..side {
                if positions.len() >= n { break 'outer; }
                // Slight jitter so particles don't start perfectly aligned
                let jitter = spacing * 0.05;
                positions.push(Vec3::new(
                    offset + gx as f32 * spacing + rng.gen_range(-jitter..jitter),
                    offset + gy as f32 * spacing + rng.gen_range(-jitter..jitter),
                    offset + gz as f32 * spacing + rng.gen_range(-jitter..jitter),
                ));
                directions.push(random_unit(&mut rng));
            }
        }
    }

    (positions, directions)
}

/// Uniformly random unit vector on S².
pub fn random_unit(rng: &mut impl Rng) -> Vec3 {
    loop {
        let v = Vec3::new(
            rng.gen_range(-1.0f32..1.0),
            rng.gen_range(-1.0f32..1.0),
            rng.gen_range(-1.0f32..1.0),
        );
        let l2 = v.length_squared();
        if l2 > 0.001 && l2 <= 1.0 {
            return v / l2.sqrt();
        }
    }
}

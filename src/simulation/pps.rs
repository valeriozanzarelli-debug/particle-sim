//! 3D Primordial Particle System (PPS) update rule.
//!
//! Rule per particle i:
//!   1. Build local frame (forward = dir_i, right, up).
//!   2. Count neighbours within R, split into left/right/up/down hemispheres.
//!   3. Turn yaw   = α + β · (N_r+N_l) · sign(N_r−N_l)
//!      Turn pitch = α + β · (N_u+N_d) · sign(N_u−N_d)
//!   4. Rotate heading by those angles (Rodrigues).
//!   5. Advance position, wrap toroidally.
//!
//! Reference: Schmickl et al., "A primordial Particle System" (2016),
//! extended here to 3D.

use bevy::prelude::*;
use rayon::prelude::*;

use crate::simulation::world::ParticleWorld;
use crate::simulation::spatial::{SpatialGrid, wrap_pos, min_image};

// ─── Main entry point ─────────────────────────────────────────────────────────

/// Run one (or more) PPS ticks. Called every Bevy FixedUpdate.
pub fn step_simulation(mut world: ResMut<ParticleWorld>) {
    let steps = world.steps_per_frame;
    for _ in 0..steps {
        tick(&mut world);
    }
}

fn tick(world: &mut ParticleWorld) {
    let n        = world.n;
    let alpha    = world.alpha;
    let beta     = world.beta;
    let r        = world.r_neighbor;
    let speed    = world.speed;
    let half     = world.box_half;

    let grid = SpatialGrid::build(&world.positions, r, half);

    // Immutable snapshots (needed for correct parallel reads)
    let pos_snap: &[Vec3] = &world.positions.clone();
    let dir_snap: &[Vec3] = &world.directions.clone();

    // Parallel update ─────────────────────────────────────────────────────────
    let results: Vec<(Vec3, Vec3, u8)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let pos = pos_snap[i];
            let dir = dir_snap[i];

            // Local orthonormal frame
            let forward = dir;
            let right = if forward.cross(Vec3::Y).length_squared() > 1e-6 {
                forward.cross(Vec3::Y).normalize()
            } else {
                forward.cross(Vec3::Z).normalize()
            };
            let up = right.cross(forward).normalize();

            // Neighbour scan
            let candidates = grid.candidates(pos, 0);
            let mut n_r = 0i32; let mut n_l = 0i32;
            let mut n_u = 0i32; let mut n_d = 0i32;
            let mut total = 0u32;

            for &j in candidates.iter() {
                let j = j as usize;
                if j == i { continue; }
                let delta = min_image(pos_snap[j] - pos, half);
                let d2 = delta.length_squared();
                if d2 < r * r && d2 > 1e-8 {
                    total += 1;
                    if delta.dot(right) >= 0.0 { n_r += 1; } else { n_l += 1; }
                    if delta.dot(up)    >= 0.0 { n_u += 1; } else { n_d += 1; }
                }
            }

            // Turning angles (PPS rule)
            let n_lr = (n_r + n_l) as f32;
            let n_ud = (n_u + n_d) as f32;
            let sign_lr = (n_r - n_l).signum() as f32;
            let sign_ud = (n_u - n_d).signum() as f32;

            let d_yaw   = alpha + beta * n_lr * sign_lr;
            let d_pitch = alpha + beta * n_ud * sign_ud;

            // Rodrigues rotation: first yaw around 'up', then pitch around 'right'
            let new_dir = rodrigues(
                rodrigues(forward, up,    d_yaw),
                right, d_pitch,
            ).normalize();

            // Advance + toroidal wrap
            let new_pos = wrap_pos(pos + new_dir * speed, half);

            (new_pos, new_dir, total.min(255) as u8)
        })
        .collect();

    // Write back
    let mut sum_n = 0u64;
    let mut max_n = 0u8;
    for (i, (pos, dir, count)) in results.into_iter().enumerate() {
        world.positions[i]       = pos;
        world.directions[i]      = dir;
        world.neighbor_counts[i] = count;
        sum_n += count as u64;
        if count > max_n { max_n = count; }
    }
    world.mean_neighbors = sum_n as f32 / n as f32;
    world.max_cluster    = max_n;
    world.ticks         += 1;
}

// ─── Rodrigues rotation ────────────────────────────────────────────────────────
/// Rotate vector `v` around unit axis `k` by `angle` radians.
#[inline]
fn rodrigues(v: Vec3, k: Vec3, angle: f32) -> Vec3 {
    let (sin_a, cos_a) = angle.sin_cos();
    v * cos_a + k.cross(v) * sin_a + k * k.dot(v) * (1.0 - cos_a)
}

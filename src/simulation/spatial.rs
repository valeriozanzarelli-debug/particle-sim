//! Spatial hash grid for O(1) average neighbour lookups in a toroidal box.

use bevy::math::Vec3;
use std::collections::HashMap;

pub struct SpatialGrid {
    cells:     HashMap<(i32, i32, i32), Vec<u32>>,
    cell_size: f32,
    box_half:  f32,
    /// Number of cells per axis (for toroidal wrapping).
    n_cells:   i32,
}

impl SpatialGrid {
    /// Build grid from a slice of positions.
    pub fn build(positions: &[Vec3], cell_size: f32, box_half: f32) -> Self {
        let n_cells = ((box_half * 2.0) / cell_size).ceil() as i32;
        let mut cells: HashMap<(i32, i32, i32), Vec<u32>> = HashMap::new();
        for (i, &pos) in positions.iter().enumerate() {
            let key = Self::pos_to_key(pos, cell_size, box_half, n_cells);
            cells.entry(key).or_default().push(i as u32);
        }
        Self { cells, cell_size, box_half, n_cells }
    }

    /// All particle indices in cells within `range_cells` of `pos`.
    /// Uses toroidal wrapping: no particle is "outside" the box.
    pub fn candidates(&self, pos: Vec3, extra_range: i32) -> Vec<u32> {
        let (cx, cy, cz) = Self::pos_to_key(pos, self.cell_size, self.box_half, self.n_cells);
        let nc = self.n_cells;
        let mut out = Vec::new();
        let r = 1 + extra_range;
        for dz in -r..=r {
            for dy in -r..=r {
                for dx in -r..=r {
                    let key = (
                        Self::wrap_key(cx + dx, nc),
                        Self::wrap_key(cy + dy, nc),
                        Self::wrap_key(cz + dz, nc),
                    );
                    if let Some(bucket) = self.cells.get(&key) {
                        out.extend_from_slice(bucket);
                    }
                }
            }
        }
        out
    }

    fn pos_to_key(pos: Vec3, cell_size: f32, box_half: f32, n_cells: i32) -> (i32, i32, i32) {
        let f = |x: f32| {
            let k = ((x + box_half) / cell_size).floor() as i32;
            Self::wrap_key(k, n_cells)
        };
        (f(pos.x), f(pos.y), f(pos.z))
    }

    #[inline]
    fn wrap_key(k: i32, n: i32) -> i32 {
        ((k % n) + n) % n
    }
}

// ─── Toroidal helpers (used by pps.rs too) ────────────────────────────────────

/// Wrap a single coordinate into [−half, +half).
#[inline]
pub fn wrap_coord(x: f32, half: f32) -> f32 {
    let size = half * 2.0;
    let x = x % size;
    if x > half  { x - size }
    else if x < -half { x + size }
    else { x }
}

/// Wrap all three coords.
#[inline]
pub fn wrap_pos(pos: Vec3, half: f32) -> Vec3 {
    Vec3::new(wrap_coord(pos.x, half), wrap_coord(pos.y, half), wrap_coord(pos.z, half))
}

/// Minimum-image displacement (toroidal): shortest vector from a to b.
#[inline]
pub fn min_image(delta: Vec3, half: f32) -> Vec3 {
    Vec3::new(
        min_image_coord(delta.x, half),
        min_image_coord(delta.y, half),
        min_image_coord(delta.z, half),
    )
}

#[inline]
fn min_image_coord(dx: f32, half: f32) -> f32 {
    if dx >  half { dx - 2.0 * half }
    else if dx < -half { dx + 2.0 * half }
    else { dx }
}

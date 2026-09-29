// 3D Primordial Particle System — GPU compute shader (WGSL)
//
// Three entry points, dispatched in sequence each simulation tick:
//   1. clear_grid  — zeroes the spatial-hash cell counters
//   2. build_grid  — inserts each particle into its grid cell (atomic)
//   3. simulate    — reads neighbours, applies PPS rule, writes new state

struct SimParams {
    n:             u32,   // particle count
    n_cells_total: u32,   // total grid cells
    n_cells_axis:  u32,   // cells per axis
    max_per_cell:  u32,   // hard cap on particles stored per cell
    alpha:         f32,   // base turn angle (rad)
    beta:          f32,   // neighbour-influence angle (rad/nbr)
    r_neighbor:    f32,   // detection radius
    speed:         f32,   // forward speed (sim units/tick)
    box_half:      f32,   // half-side of the toroidal cube
    cell_size:     f32,   // spatial grid cell size (= r_neighbor, clamped)
    _pad0:         u32,
    _pad1:         u32,
}

@group(0) @binding(0) var<uniform>             p          : SimParams;
@group(0) @binding(1) var<storage, read>       pos_in     : array<vec4<f32>>;
@group(0) @binding(2) var<storage, read>       dir_in     : array<vec4<f32>>;
@group(0) @binding(3) var<storage, read_write> pos_out    : array<vec4<f32>>;
@group(0) @binding(4) var<storage, read_write> dir_out    : array<vec4<f32>>;
@group(0) @binding(5) var<storage, read_write> nbr_out    : array<u32>;
@group(0) @binding(6) var<storage, read_write> cell_counts: array<atomic<u32>>;
@group(0) @binding(7) var<storage, read_write> cell_data  : array<u32>;

// ── helpers ───────────────────────────────────────────────────────────────────

fn wrap1(x: f32) -> f32 {
    let h = p.box_half;
    let s = h * 2.0;
    var r = x % s;
    if r >  h { r -= s; }
    if r < -h { r += s; }
    return r;
}

fn wrap3(v: vec3<f32>) -> vec3<f32> {
    return vec3(wrap1(v.x), wrap1(v.y), wrap1(v.z));
}

fn min1(dx: f32) -> f32 {
    if dx >  p.box_half { return dx - 2.0 * p.box_half; }
    if dx < -p.box_half { return dx + 2.0 * p.box_half; }
    return dx;
}

fn min_img(d: vec3<f32>) -> vec3<f32> {
    return vec3(min1(d.x), min1(d.y), min1(d.z));
}

fn cell_coord(x: f32) -> i32 {
    let nc = i32(p.n_cells_axis);
    var k = i32(floor((x + p.box_half) / p.cell_size));
    k = ((k % nc) + nc) % nc;
    return k;
}

fn cell_idx(c: vec3<i32>) -> u32 {
    let nc = i32(p.n_cells_axis);
    let cx = ((c.x % nc) + nc) % nc;
    let cy = ((c.y % nc) + nc) % nc;
    let cz = ((c.z % nc) + nc) % nc;
    return u32(cz * nc * nc + cy * nc + cx);
}

fn rodrigues(v: vec3<f32>, k: vec3<f32>, angle: f32) -> vec3<f32> {
    let s = sin(angle);
    let c = cos(angle);
    return v * c + cross(k, v) * s + k * dot(k, v) * (1.0 - c);
}

// ── pass 1: clear grid ────────────────────────────────────────────────────────

@compute @workgroup_size(64)
fn clear_grid(@builtin(global_invocation_id) gid: vec3<u32>) {
    if gid.x < p.n_cells_total {
        atomicStore(&cell_counts[gid.x], 0u);
    }
}

// ── pass 2: build grid ────────────────────────────────────────────────────────

@compute @workgroup_size(64)
fn build_grid(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if i >= p.n { return; }

    let pos  = pos_in[i].xyz;
    let c    = vec3(cell_coord(pos.x), cell_coord(pos.y), cell_coord(pos.z));
    let cidx = cell_idx(c);
    let slot = atomicAdd(&cell_counts[cidx], 1u);
    if slot < p.max_per_cell {
        cell_data[cidx * p.max_per_cell + slot] = i;
    }
}

// ── pass 3: simulate ─────────────────────────────────────────────────────────

@compute @workgroup_size(64)
fn simulate(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if i >= p.n { return; }

    let pos = pos_in[i].xyz;
    let fwd = normalize(dir_in[i].xyz);

    // Build local orthonormal frame
    var right: vec3<f32>;
    if length(cross(fwd, vec3(0.0, 1.0, 0.0))) > 0.001 {
        right = normalize(cross(fwd, vec3(0.0, 1.0, 0.0)));
    } else {
        right = normalize(cross(fwd, vec3(0.0, 0.0, 1.0)));
    }
    let up = normalize(cross(right, fwd));

    let cx  = cell_coord(pos.x);
    let cy  = cell_coord(pos.y);
    let cz  = cell_coord(pos.z);
    let r2  = p.r_neighbor * p.r_neighbor;

    var n_r: i32 = 0;
    var n_l: i32 = 0;
    var n_u: i32 = 0;
    var n_d: i32 = 0;
    var total: u32 = 0u;

    for (var dz = -1; dz <= 1; dz++) {
        for (var dy = -1; dy <= 1; dy++) {
            for (var dx = -1; dx <= 1; dx++) {
                let nc  = vec3(cx + dx, cy + dy, cz + dz);
                let cidx = cell_idx(nc);
                let cnt  = min(atomicLoad(&cell_counts[cidx]), p.max_per_cell);
                for (var s = 0u; s < cnt; s++) {
                    let j = cell_data[cidx * p.max_per_cell + s];
                    if j == i { continue; }
                    let delta = min_img(pos_in[j].xyz - pos);
                    let d2 = dot(delta, delta);
                    if d2 < r2 && d2 > 0.0001 {
                        total += 1u;
                        if dot(delta, right) >= 0.0 { n_r += 1; } else { n_l += 1; }
                        if dot(delta, up)    >= 0.0 { n_u += 1; } else { n_d += 1; }
                    }
                }
            }
        }
    }

    let n_lr    = f32(n_r + n_l);
    let n_ud    = f32(n_u + n_d);
    let sign_lr = f32(sign(f32(n_r - n_l)));
    let sign_ud = f32(sign(f32(n_u - n_d)));

    let d_yaw   = p.alpha + p.beta * n_lr * sign_lr;
    let d_pitch = p.alpha + p.beta * n_ud * sign_ud;

    let new_dir = normalize(rodrigues(rodrigues(fwd, up, d_yaw), right, d_pitch));
    let new_pos = wrap3(pos + new_dir * p.speed);

    pos_out[i] = vec4(new_pos, 0.0);
    dir_out[i] = vec4(new_dir, 0.0);
    nbr_out[i] = min(total, 255u);
}

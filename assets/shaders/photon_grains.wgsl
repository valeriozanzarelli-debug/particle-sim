// Photon grains on the GPU (WGSL) — port of `advance()` in src/photon/physics.rs.
//
// Every grain translates at c along `dir`; its phase is phase0 + k·|pos − origin|.
// Optical elements (per photon, sorted along travel) are met strictly in order:
//   Wall      → Huygens re-emission from a slit (deterministic per-grain hash
//               for the angle, bit-identical to physics.rs `hash3`), else dies
//   Polarizer / BellSplitter → pass (photon-level decision, see `decide`)
//   Screen    → the grain dies and deposits its complex amplitude into a bin.
//
// Three entry points per tick, dispatched in this order in ONE command buffer
// (the host clears `status` first):
//   1. probe   — 1 thread / grain: which elements does the photon's front reach
//                within `dist`, and the Born-rule sums for the front element
//                (workgroup reduction → 1 atomic per workgroup per word)
//   2. decide  — 1 thread / photon slot: Born rule with a pre-drawn uniform,
//                writes photon.op (KILL / SET_JONES) and photon.front
//   3. advance — 1 thread / grain: apply photon.op, move, cross planes, deposit
//                on the screen, count survivors (workgroup reduction)
//
// Nothing is read back per tick. The host reads only `status` (32 B/photon)
// and `screen` (48 B/bin/photon) when it needs them.
//
// Invariant required by the workgroup reductions: each photon's grain range
// starts at a multiple of 64 (WG), so a workgroup never spans two photons.
//
// ── Layouts (std430 / uniform; mirrored with bytemuck in src/photon/gpu.rs) ──
//   Params  96 B   Grain  64 B   Photon  112 B
//   status  8 × u32 per photon slot:
//           [0] reach mask (bit e = element e reached this tick)
//           [1] grains alive after advance
//           [2,3] den  i64 fixed point (lo, hi)  — Σ|J|²        (× prob_scale)
//           [4,5] num  i64 fixed point (lo, hi)  — Σ through/pass (× prob_scale)
//           [6] fate   [7] front (next photon-level element)
//   screen  12 × u32 per (slot, bin):
//           [0,1] Re a_e1  [2,3] Im a_e1  [4,5] Re a_e2  [6,7] Im a_e2
//           (i64 fixed point × amp_scale)
//           [8,9] Σ|a|² self term (i64 × self_scale)   [10] count   [11] pad

const TAU: f32 = 6.283185307179586;
const PI:  f32 = 3.141592653589793;
const WG:  u32 = 64u;

// Grain.state bits
const NEXT_MASK: u32 = 0xFFu;      // bits 0-7  : next element index
const ST_ALIVE:  u32 = 0x100u;     // bit 8
const ST_DIFFR:  u32 = 0x200u;     // bit 9     : re-emitted by a slit
// bits 16-31: photon slot

// Element kinds (Photon.elem_kind low byte; high bits = lab polarizer index)
const K_WALL:      u32 = 1u;
const K_POLARIZER: u32 = 2u;
const K_BELL:      u32 = 3u;
const K_SCREEN:    u32 = 4u;

const OP_NONE:      u32 = 0u;
const OP_KILL:      u32 = 1u;
const OP_SET_JONES: u32 = 2u;

const FATE_IN_FLIGHT: u32 = 0u;
const FATE_WALL:      u32 = 1u;
const FATE_POLARIZER: u32 = 2u;    // | (lab polarizer index << 8)

const STATUS_WORDS: u32 = 8u;
const SCREEN_WORDS: u32 = 12u;
const INVALID: u32 = 0xFFFFFFFFu;

struct Params {
    dist:        f32,       // light path this tick
    k:           f32,       // 2π/λ
    lambda:      f32,
    waist:       f32,       // w0 (Huygens 1/|U_inc| weight)
    theta_max:   f32,       // half-angle of the re-emission fan
    screen_half: f32,
    n_bins:      u32,
    n_grains:    u32,       // grains in the buffer (multiple of 64)
    slit_c0:     f32,
    slit_c1:     f32,
    slit_half:   f32,
    filter_mask: u32,       // bit s = slit s has a polarizer
    filter0:     vec4<f32>, // Jones (re0, im0, re1, im1) of slit 0 polarizer
    filter1:     vec4<f32>,
    amp_scale:   f32,       // fixed point scale of the amplitudes
    self_scale:  f32,       // fixed point scale of |a|²
    prob_scale:  f32,       // fixed point scale of the Born-rule sums
    n_slots:     u32,       // photon slots
}

struct Grain {
    pos:    vec3<f32>,
    amp:    f32,
    dir:    vec3<f32>,
    phase0: f32,
    origin: vec3<f32>,
    state:  u32,
    jones:  vec4<f32>,      // (re e1, im e1, re e2, im e2)
}

struct Photon {
    elem_z:      vec4<f32>, // up to 4 elements, sorted along travel
    elem_kind:   vec4<u32>,
    elem_arg:    vec4<f32>, // polarizer angle
    rand:        vec4<f32>, // pre-drawn uniforms for the Born rule, per element
    op_jones:    vec4<f32>, // written by `decide`
    n_elems:     u32,
    key_lo:      u32,       // seed ^ photon_id·0x9E3779B97F4A7C15 (u64)
    key_hi:      u32,
    first_grain: u32,       // multiple of 64
    n_grains:    u32,
    front:       u32,       // next element for photon-level events
    op:          u32,
    fate:        u32,
}

@group(0) @binding(0) var<uniform>             p       : Params;
@group(0) @binding(1) var<storage, read_write> grains  : array<Grain>;
@group(0) @binding(2) var<storage, read_write> photons : array<Photon>;
@group(0) @binding(3) var<storage, read_write> status  : array<atomic<u32>>;
@group(0) @binding(4) var<storage, read_write> screen  : array<atomic<u32>>;

var<workgroup> sh_f0:   array<f32, 64>;
var<workgroup> sh_f1:   array<f32, 64>;
var<workgroup> sh_u:    array<u32, 64>;
var<workgroup> sh_slot: array<u32, 64>;

// ── complex helpers ───────────────────────────────────────────────────────────

fn cmul(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
    return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
fn conj(a: vec2<f32>) -> vec2<f32> { return vec2(a.x, -a.y); }
/// ⟨a|b⟩
fn inner(a: vec4<f32>, b: vec4<f32>) -> vec2<f32> {
    return cmul(conj(a.xy), b.xy) + cmul(conj(a.zw), b.zw);
}
/// Reduce to [−π, π] before sin/cos.
fn wrap_phase(x: f32) -> f32 { return x - TAU * floor(x / TAU + 0.5); }

/// (sin x, cos x) to ~1e-7 absolute for |x| ≲ 2π. WGSL only guarantees
/// 2^-11 for the built-in sin/cos, and the Intel UHD driver really is that
/// coarse (1.4e-5 measured on the re-emission angle → 1e-3 on positions
/// 70 units later), so the grain physics uses this instead.
fn sincos(x: f32) -> vec2<f32> {
    let q = round(x * 0.6366197723675814);               // x / (π/2)
    let r = (x - q * 1.5707963705062866) + q * 4.371139000186243e-8;
    let r2 = r * r;
    let s = r + r * r2 * (-0.16666666666666666 + r2 * (0.008333333333333333
            + r2 * (-0.0001984126984126984 + r2 * 2.7557319223985893e-6)));
    let c = 1.0 + r2 * (-0.5 + r2 * (0.041666666666666664 + r2 * (-0.001388888888888889
            + r2 * (2.48015873015873e-5 - r2 * 2.755731922398589e-7))));
    switch (i32(q) & 3) {
        case 0: { return vec2(s, c); }
        case 1: { return vec2(c, -s); }
        case 2: { return vec2(-s, -c); }
        default: { return vec2(-c, s); }
    }
}

// ── 64-bit integer emulation (u64 as vec2(lo, hi)) for the grain hash ────────

fn mul_wide(x: u32, y: u32) -> vec2<u32> {
    let xl = x & 0xFFFFu; let xh = x >> 16u;
    let yl = y & 0xFFFFu; let yh = y >> 16u;
    let ll = xl * yl; let lh = xl * yh; let hl = xh * yl; let hh = xh * yh;
    let mid = (ll >> 16u) + (lh & 0xFFFFu) + (hl & 0xFFFFu);
    return vec2((ll & 0xFFFFu) | (mid << 16u), hh + (lh >> 16u) + (hl >> 16u) + (mid >> 16u));
}
fn mul64(a: vec2<u32>, b: vec2<u32>) -> vec2<u32> {
    let ll = mul_wide(a.x, b.x);
    return vec2(ll.x, ll.y + a.x * b.y + a.y * b.x);
}
fn shr64(a: vec2<u32>, n: u32) -> vec2<u32> {   // 0 < n < 32
    return vec2((a.x >> n) | (a.y << (32u - n)), a.y >> n);
}
/// physics.rs `hash3(seed, photon_id, gi)`, with the first two terms folded in `key`.
fn grain_hash(key: vec2<u32>, gi: u32) -> vec2<u32> {
    var x = key ^ mul64(vec2(gi, 0u), vec2(0x27D4EB4Fu, 0xC2B2AE3Du));
    x = x ^ shr64(x, 30u); x = mul64(x, vec2(0x1CE4E5B9u, 0xBF58476Du));
    x = x ^ shr64(x, 27u); x = mul64(x, vec2(0x133111EBu, 0x94D049BBu));
    return x ^ shr64(x, 31u);
}

// ── fixed-point i64 atomics (lo word + carry into hi word) ───────────────────
// Integer adds are associative: the sums do not depend on thread order.

fn to_fixed(x: f32) -> i32 { return i32(round(clamp(x, -2.0e9, 2.0e9))); }

fn screen_add(w: u32, v: i32) {
    if (v == 0) { return; }
    let lo = bitcast<u32>(v);
    let old = atomicAdd(&screen[w], lo);
    var hi = select(0u, 0xFFFFFFFFu, v < 0);
    if (old + lo < old) { hi = hi + 1u; }
    if (hi != 0u) { atomicAdd(&screen[w + 1u], hi); }
}
fn status_add(w: u32, v: i32) {
    if (v == 0) { return; }
    let lo = bitcast<u32>(v);
    let old = atomicAdd(&status[w], lo);
    var hi = select(0u, 0xFFFFFFFFu, v < 0);
    if (old + lo < old) { hi = hi + 1u; }
    if (hi != 0u) { atomicAdd(&status[w + 1u], hi); }
}
fn status_load64(w: u32) -> f32 {   // non-negative sums only
    return f32(atomicLoad(&status[w + 1u])) * 4294967296.0 + f32(atomicLoad(&status[w]));
}

fn grain_index(gid: vec3<u32>, nwg: vec3<u32>) -> u32 {
    return gid.x + gid.y * nwg.x * WG;
}

// ── grain physics (mirrors physics.rs) ───────────────────────────────────────

fn huygens(g_in: Grain, key: vec2<u32>, gi: u32) -> Grain {
    var g = g_in;
    let x = g.pos.x;
    var s = 2u;
    if (abs(x - p.slit_c0) <= p.slit_half) { s = 0u; }
    else if (abs(x - p.slit_c1) <= p.slit_half) { s = 1u; }
    if (s == 2u) { g.state = g.state & ~ST_ALIVE; return g; }

    let phase_here = g.phase0 + p.k * length(g.pos - g.origin);
    if (((p.filter_mask >> s) & 1u) != 0u) {
        let f = select(p.filter0, p.filter1, s == 1u);
        let c = inner(f, g.jones);
        g.jones = vec4(cmul(f.xy, c), cmul(f.zw, c));
    }
    // u = (rnd >> 11) / 2^53, as in physics.rs
    let r = shr64(grain_hash(key, gi), 11u);
    let u = f32(r.y) * 4.76837158203125e-7 + f32(r.x) * 1.1102230246251565e-16;
    let theta = p.theta_max * (2.0 * u - 1.0);
    let sgn = select(-1.0, 1.0, g.dir.z >= 0.0);
    let sc = sincos(theta);
    g.dir = normalize(vec3(sc.x, g.dir.y, sgn * sc.y));
    let inc = exp(-(x * x) / (p.waist * p.waist));
    g.amp = 0.5 * (1.0 + sc.y) / inc;
    g.origin = g.pos;
    g.phase0 = phase_here - 0.25 * PI;
    g.state = g.state | ST_DIFFR;
    return g;
}

fn deposit(g: Grain, slot: u32) {
    let f = (g.pos.x + p.screen_half) / (2.0 * p.screen_half);
    if (!(f >= 0.0 && f < 1.0)) { return; }
    let bin = min(u32(f * f32(p.n_bins)), p.n_bins - 1u);
    let base = (slot * p.n_bins + bin) * SCREEN_WORDS;
    atomicAdd(&screen[base + 10u], 1u);
    if ((g.state & ST_DIFFR) == 0u) { return; }
    // Huygens kernel e^{ikr}/√(λr), Jacobian dθ/dx = cos²θ / L
    let dx = g.pos.x - g.origin.x;
    let dz = g.pos.z - g.origin.z;
    let r = sqrt(dx * dx + dz * dz);
    let cos_t = abs(dz) / r;
    let w = g.amp * abs(dz) / (cos_t * cos_t * sqrt(p.lambda * r));
    let ph = wrap_phase(g.phase0 + p.k * length(g.pos - g.origin));
    let e = sincos(ph).yx * w;
    let a0 = cmul(g.jones.xy, e);
    let a1 = cmul(g.jones.zw, e);
    screen_add(base + 0u, to_fixed(a0.x * p.amp_scale));
    screen_add(base + 2u, to_fixed(a0.y * p.amp_scale));
    screen_add(base + 4u, to_fixed(a1.x * p.amp_scale));
    screen_add(base + 6u, to_fixed(a1.y * p.amp_scale));
    screen_add(base + 8u, to_fixed((dot(a0, a0) + dot(a1, a1)) * p.self_scale));
}

/// Translate by `dist`, meeting the elements strictly in order (physics.rs `advance`).
fn advance_grain(g_in: Grain, ph: Photon, gi: u32, slot: u32) -> Grain {
    var g = g_in;
    var dist = p.dist;
    let key = vec2(ph.key_lo, ph.key_hi);
    loop {
        if ((g.state & ST_ALIVE) == 0u) { break; }
        let next = g.state & NEXT_MASK;
        if (next >= ph.n_elems) { g.pos += g.dir * dist; break; }
        let z = ph.elem_z[next];
        let t = max((z - g.pos.z) / g.dir.z, 0.0);
        if (t > dist) { g.pos += g.dir * dist; break; }
        g.pos += g.dir * t;
        g.pos.z = z;
        dist -= t;
        g.state += 1u;                       // next element
        let kind = ph.elem_kind[next] & 0xFFu;
        if (kind == K_WALL) {
            g = huygens(g, key, gi);
        } else if (kind == K_SCREEN) {
            g.state = g.state & ~ST_ALIVE;
            deposit(g, slot);
            break;
        }
        // K_POLARIZER / K_BELL: photon-level, nothing per grain
    }
    return g;
}

// ── pass 1: probe ────────────────────────────────────────────────────────────

@compute @workgroup_size(64)
fn probe(@builtin(global_invocation_id) gid: vec3<u32>,
         @builtin(local_invocation_index) lid: u32,
         @builtin(num_workgroups) nwg: vec3<u32>) {
    let i = grain_index(gid, nwg);
    var mask = 0u;
    var den = 0.0;
    var num = 0.0;
    var slot = INVALID;
    if (i < p.n_grains) {
        let g = grains[i];
        let s = g.state >> 16u;
        if ((g.state & ST_ALIVE) != 0u && s < p.n_slots) {
            let ph = photons[s];
            if (ph.fate == FATE_IN_FLIGHT && ph.front < ph.n_elems) {
                slot = s;
                for (var e = ph.front; e < ph.n_elems; e++) {
                    let t = (ph.elem_z[e] - g.pos.z) / g.dir.z;
                    if (t >= 0.0 && t <= p.dist) { mask = mask | (1u << e); }
                }
                let kind = ph.elem_kind[ph.front] & 0xFFu;
                let w = dot(g.jones, g.jones);
                if (kind == K_WALL) {
                    let t = (ph.elem_z[ph.front] - g.pos.z) / g.dir.z;
                    let x = g.pos.x + g.dir.x * t;
                    den = w;
                    if (abs(x - p.slit_c0) <= p.slit_half) {
                        num += select(w, dot(inner(p.filter0, g.jones), inner(p.filter0, g.jones)), (p.filter_mask & 1u) != 0u);
                    }
                    if (abs(x - p.slit_c1) <= p.slit_half) {
                        num += select(w, dot(inner(p.filter1, g.jones), inner(p.filter1, g.jones)), (p.filter_mask & 2u) != 0u);
                    }
                } else if (kind == K_POLARIZER) {
                    let a = ph.elem_arg[ph.front];
                    let sa = sincos(a);
                    let c = inner(vec4(sa.y, 0.0, sa.x, 0.0), g.jones);
                    den = w;
                    num = dot(c, c);
                }
            }
        }
    }
    sh_u[lid] = mask; sh_f0[lid] = den; sh_f1[lid] = num; sh_slot[lid] = slot;
    workgroupBarrier();
    for (var st = 32u; st > 0u; st = st >> 1u) {
        if (lid < st) {
            sh_u[lid] = sh_u[lid] | sh_u[lid + st];
            sh_f0[lid] += sh_f0[lid + st];
            sh_f1[lid] += sh_f1[lid + st];
            sh_slot[lid] = min(sh_slot[lid], sh_slot[lid + st]);
        }
        workgroupBarrier();
    }
    if (lid == 0u && sh_slot[0] != INVALID) {
        let base = sh_slot[0] * STATUS_WORDS;
        if (sh_u[0] != 0u) { atomicOr(&status[base], sh_u[0]); }
        status_add(base + 2u, to_fixed(sh_f0[0] * p.prob_scale));
        status_add(base + 4u, to_fixed(sh_f1[0] * p.prob_scale));
    }
}

// ── pass 2: decide (one thread per photon) ───────────────────────────────────

@compute @workgroup_size(64)
fn decide(@builtin(global_invocation_id) gid: vec3<u32>) {
    let s = gid.x;
    if (s >= p.n_slots) { return; }
    var ph = photons[s];
    ph.op = OP_NONE;
    let base = s * STATUS_WORDS;
    if (ph.fate == FATE_IN_FLIGHT) {
        let mask = atomicLoad(&status[base]);
        let den = status_load64(base + 2u);
        let num = status_load64(base + 4u);
        let first = ph.front;
        var known = false;             // jones set by an earlier element this tick
        var kj = vec4(0.0);
        loop {
            if (ph.front >= ph.n_elems) { break; }
            let e = ph.front;
            if (((mask >> e) & 1u) == 0u) { break; }
            let kind = ph.elem_kind[e] & 0xFFu;
            ph.front += 1u;
            if (kind == K_WALL) {
                var pp = 0.0;
                if (den > 0.0 && e == first) { pp = num / den; }
                if (ph.rand[e] >= pp) { ph.fate = FATE_WALL; ph.op = OP_KILL; break; }
            } else if (kind == K_POLARIZER) {
                let sa = sincos(ph.elem_arg[e]);
                let a = vec4(sa.y, 0.0, sa.x, 0.0);
                var pp = 0.0;
                if (known) {
                    let c = inner(a, kj);
                    pp = dot(c, c) / max(dot(kj, kj), 1e-30);
                } else if (den > 0.0 && e == first) {
                    pp = num / den;
                }
                if (ph.rand[e] < pp) {
                    ph.op = OP_SET_JONES; ph.op_jones = a; known = true; kj = a;
                } else {
                    ph.fate = FATE_POLARIZER | (ph.elem_kind[e] & 0xFFFFFF00u);
                    ph.op = OP_KILL;
                    break;
                }
            }
            // K_BELL: the entangled collapse stays on the host; K_SCREEN: nothing
        }
    }
    photons[s] = ph;
    atomicStore(&status[base + 6u], ph.fate);
    atomicStore(&status[base + 7u], ph.front);
}

// ── pass 3: advance ──────────────────────────────────────────────────────────

@compute @workgroup_size(64)
fn advance(@builtin(global_invocation_id) gid: vec3<u32>,
           @builtin(local_invocation_index) lid: u32,
           @builtin(num_workgroups) nwg: vec3<u32>) {
    let i = grain_index(gid, nwg);
    var alive = 0u;
    var slot = INVALID;
    if (i < p.n_grains) {
        var g = grains[i];
        let s = g.state >> 16u;
        if ((g.state & ST_ALIVE) != 0u && s < p.n_slots) {
            slot = s;
            let ph = photons[s];
            if (ph.op == OP_KILL) {
                g.state = g.state & ~ST_ALIVE;
            } else {
                if (ph.op == OP_SET_JONES) { g.jones = ph.op_jones; }
                g = advance_grain(g, ph, i - ph.first_grain, s);
            }
            if ((g.state & ST_ALIVE) != 0u) { alive = 1u; }
            grains[i] = g;
        }
    }
    sh_u[lid] = alive; sh_slot[lid] = slot;
    workgroupBarrier();
    for (var st = 32u; st > 0u; st = st >> 1u) {
        if (lid < st) {
            sh_u[lid] += sh_u[lid + st];
            sh_slot[lid] = min(sh_slot[lid], sh_slot[lid + st]);
        }
        workgroupBarrier();
    }
    if (lid == 0u && sh_slot[0] != INVALID && sh_u[0] != 0u) {
        atomicAdd(&status[sh_slot[0] * STATUS_WORDS + 1u], sh_u[0]);
    }
}

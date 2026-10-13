// Weighted blended order-independent transparency for the photon's energy balls (oit_palline.rs,
// docs/MOTORE_GRAFICO_OIT.md). McGuire & Bavoil 2013: every ball adds c·a·w and a·w to ACCUM and multiplies
// REVEAL by exp(-k a); the composite pass shows the w-weighted mean colour with coverage 1 − REVEAL.
// w = exp(−(z − z_rif)/D): the balls nearest to the eye count most, so the colours do not average to white.
// Same 32-byte instances as istanze.wgsl; the quad corners come from the vertex index (no extra buffer).
// Only balls: no lines, no surfaces.

struct Par {
    clip_from_world: mat4x4<f32>,
    right: vec4<f32>,      // camera right (xyz)
    up: vec4<f32>,         // camera up (xyz)
    eye: vec4<f32>,        // camera position (xyz), w = view alignment |fwd · axis|
    fwd: vec4<f32>,        // camera forward (xyz), w = z_rif
    offset: vec4<f32>,     // layer translation (xyz), w = D (depth of the transparency)
    // size, gain, soft, k (exposure)
    a: vec4<f32>,
    // tail samples ns, visible tail fraction from the side, colour gain of the composite, unused
    b: vec4<f32>,
};

@group(0) @binding(0) var<uniform> par: Par;

struct VIn {
    @builtin(vertex_index) vi: u32,
    @location(0) pos_size: vec4<f32>,
    @location(1) allunga: vec3<f32>,
    @location(2) colore: vec4<f32>,
};

struct VOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) corner: vec2<f32>,
    @location(1) color: vec4<f32>,
    @location(2) w: f32,
};

// Visible weight of a tail ball (j = 0 newest) of `ns`, seen with alignment c = |fwd · axis|.
// From the side only the newest quarter of the tail (as before); along the axis the whole tail (2 λ/c),
// with a flatter taper so the arcs of the swirl read.
fn coda(j: f32, ns: f32, c: f32) -> f32 {
    let q = (j + 1.0) / ns;
    let qmax = mix(par.b.y, 1.0, smoothstep(0.55, 0.95, c));
    let vis = clamp(1.0 + (qmax - q) * ns, 0.0, 1.0);
    let lato = 0.5 * (1.0 - j / ns);
    let dietro = 0.7 * (1.0 - 0.6 * q);
    return vis * mix(lato, dietro, smoothstep(0.55, 0.95, c));
}

@vertex
fn vertex(v: VIn) -> VOut {
    var out: VOut;
    var cs = array<vec2<f32>, 6>(vec2<f32>(-1.0, -1.0), vec2<f32>(1.0, -1.0), vec2<f32>(1.0, 1.0),
                                 vec2<f32>(-1.0, -1.0), vec2<f32>(1.0, 1.0), vec2<f32>(-1.0, 1.0));
    let corner = cs[v.vi % 6u];
    let ws = v.pos_size.w;
    var col = v.colore;
    var size = ws;
    // tails are written with size 0.85 − 0.08 j (fotdisegno::palline and gpu_fotone.wgsl): j from the size,
    // the head alpha from the tail alpha 0.5 (1 − j/ns) a, then the view-dependent tail weight
    let ns = max(par.b.x, 1.0);
    let lungo = smoothstep(0.55, 0.95, par.eye.w);
    var st = v.allunga;
    if (ws < 0.95 && ws > 0.0) {
        let j = clamp(round((0.85 - ws) / 0.08), 0.0, ns - 1.0);
        let a_testa = col.a / max(0.5 * (1.0 - j / ns), 1e-3);
        col.a = clamp(a_testa * coda(j, ns, par.eye.w), 0.0, 1.0);
        // a tail ball reaches half-way to the next newer one (a streak of soft balls) only along the axis
        st = st * lungo;
    } else {
        // along the axis the stretch along E (all directions at once) would hide the swirl: shorter there
        st = st * mix(1.0, 0.35, lungo);
    }
    let world = v.pos_size.xyz + par.offset.xyz;
    let s = par.a.x * size;
    let a = vec2<f32>(dot(st, par.right.xyz), dot(st, par.up.xyz));
    let la = length(a);
    var d = vec2<f32>(1.0, 0.0);
    if (la > 1e-4) { d = a / la; }
    let nrm = vec2<f32>(-d.y, d.x);
    let off = d * (corner.x * (s + la)) + nrm * (corner.y * s);
    out.clip = par.clip_from_world * vec4<f32>(world + par.right.xyz * off.x + par.up.xyz * off.y, 1.0);
    if (col.a < 1e-3) { out.clip = vec4<f32>(0.0, 0.0, -2.0, 1.0); }
    out.corner = corner;
    out.color = col;
    // depth weight (McGuire–Bavoil): nearer balls cover farther ones
    let z = dot(world - par.eye.xyz, par.fwd.xyz);
    out.w = exp(-clamp((z - par.fwd.w) / max(par.offset.w, 1e-3), -3.0, 14.0));
    return out;
}

struct FOut {
    @location(0) accum: vec4<f32>,
    @location(1) reveal: vec4<f32>,
};

@fragment
fn fragment(in: VOut) -> FOut {
    let r2 = dot(in.corner, in.corner);
    if (r2 > 1.0) { discard; }
    // the same soft ball as the alpha path of istanze.wgsl
    let rim = 1.0 - smoothstep(0.55, 1.0, r2);
    let f = mix(rim, exp(-par.a.z * r2), 0.6);
    let a = clamp(in.color.a * par.a.y * f, 0.0, 1.0);
    var o: FOut;
    o.accum = vec4<f32>(in.color.rgb * a * in.w, a * in.w);
    o.reveal = vec4<f32>(exp(-par.a.w * a), 0.0, 0.0, 1.0);
    return o;
}

// ── composite: weighted mean colour × coverage, premultiplied, over the scene ──

@group(0) @binding(1) var t_accum: texture_2d<f32>;
@group(0) @binding(2) var t_reveal: texture_2d<f32>;

struct FsIn {
    @builtin(position) pos: vec4<f32>,
};

@fragment
fn composita(in: FsIn) -> @location(0) vec4<f32> {
    let p = vec2<i32>(in.pos.xy);
    let acc = textureLoad(t_accum, p, 0);
    let rev = textureLoad(t_reveal, p, 0).r;
    let cop = clamp(1.0 - rev, 0.0, 1.0);
    if (cop < 1e-4) { discard; }
    let media = acc.rgb / max(acc.a, 1e-6);
    return vec4<f32>(media * cop * par.b.z, cop);
}

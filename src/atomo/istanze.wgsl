// Instanced soft balls (istanze.rs): one shared quad, one 32-byte instance per ball.
// Same look as balls.wgsl (the mesh path), so the two paths can be compared (exam R2).

#import bevy_pbr::view_transformations::position_world_to_clip
#import bevy_pbr::mesh_view_bindings::view

struct Parametri {
    size: f32,
    gain: f32,
    soft: f32,
    core: f32,
    // xyz: translation of the layer; w: the first `w` instances follow the level of detail
    offset: vec4<f32>,
    // level of detail (grafica.rs): kept fraction f, brightness factor, size factor, smallest radius in pixels
    lod: vec4<f32>,
    // depth cue: camera distance, depth scale, strength, scene units per pixel per unit distance
    prof: vec4<f32>,
    // cut towards the camera: on, front fade width, back fade length (scene units); w: largest radius in pixels (0: none)
    taglio: vec4<f32>,
};

@group(2) @binding(0) var<uniform> params: Parametri;

struct VertexIn {
    @builtin(instance_index) ii: u32,
    @location(0) corner: vec3<f32>,
    @location(3) pos_size: vec4<f32>,
    @location(4) allunga: vec3<f32>,
    @location(5) colore: vec4<f32>,
};

struct VertexOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) corner: vec2<f32>,
    @location(1) color: vec4<f32>,
};

// the stable key of a glow, the same PCG hash as `grafica::chiave` (Rust)
fn chiave(i: u32) -> f32 {
    let s = i * 747796405u + 2891336453u;
    let w = ((s >> ((s >> 28u) + 4u)) ^ s) * 277803737u;
    return f32((w >> 22u) ^ w) / 4294967296.0;
}

@vertex
fn vertex(v: VertexIn) -> VertexOut {
    var out: VertexOut;
    let right = view.world_from_view[0].xyz;
    let up = view.world_from_view[1].xyz;
    var s = params.size * v.pos_size.w;
    var amul = 1.0;
    // level of detail: a stable subset of the first instances, each a little bigger (grafica::lod)
    if (f32(v.ii) < params.offset.w) {
        let f = params.lod.x;
        if (f < 1.0) {
            amul = clamp((f - chiave(v.ii)) / (0.15 * f + 1e-6), 0.0, 1.0);
        }
        amul = amul * params.lod.y;
        s = s * params.lod.z;
    }
    let centro = v.pos_size.xyz + params.offset.xyz;
    let dist = length(centro - view.world_position);
    // never smaller than lod.w pixels
    s = max(s, params.lod.w * params.prof.w * dist);
    // and never bigger than taglio.w pixels: close up a coarse light glow does not swell over the screen; between 1 and 2
    // times that size it fades out, and the finer detail (the sparkles, zzscena::scintille) takes over
    if (params.taglio.w > 0.0) {
        let px = s / max(params.prof.w * dist, 1e-9);
        amul = amul * clamp(2.0 - px / params.taglio.w, 0.0, 1.0);
        s = min(s, params.taglio.w * params.prof.w * dist);
    }
    // depth: nearer glows brighter and a little bigger, farther ones dimmer
    if (params.prof.z > 0.0) {
        let dd = (dist - params.prof.x) / max(params.prof.y, 1e-6);
        amul = amul * clamp(1.0 - params.prof.z * dd, 0.25, 1.6);
        s = s * clamp(1.0 - 0.25 * params.prof.z * dd, 0.7, 1.3);
    }
    // the light cut open towards the camera: in front of the plane through the centre it fades fast, behind it slowly
    if (params.taglio.x > 0.5) {
        let e = normalize(view.world_position - params.offset.xyz);
        let d = dot(v.pos_size.xyz, e);
        if (d > 0.0) {
            amul = amul * exp(-(d * d) / (params.taglio.y * params.taglio.y));
        } else {
            amul = amul * exp(d / params.taglio.z);
        }
    }
    if (amul < 0.003) {
        // nothing to draw: the quad leaves the screen
        out.clip = vec4<f32>(0.0, 0.0, 2.0, 1.0);
        out.corner = vec2<f32>(0.0, 0.0);
        out.color = vec4<f32>(0.0);
        return out;
    }
    let a = vec2<f32>(dot(v.allunga, right), dot(v.allunga, up));
    let la = length(a);
    var d = vec2<f32>(1.0, 0.0);
    if (la > 1e-4) {
        d = a / la;
    }
    let nrm = vec2<f32>(-d.y, d.x);
    let off = d * (v.corner.x * (s + la)) + nrm * (v.corner.y * s);
    out.clip = position_world_to_clip(centro + right * off.x + up * off.y);
    out.corner = v.corner.xy;
    out.color = vec4<f32>(v.colore.rgb, v.colore.a * amul);
    return out;
}

@fragment
fn fragment(in: VertexOut) -> @location(0) vec4<f32> {
    let r2 = dot(in.corner, in.corner);
    if (r2 > 1.0) {
        discard;
    }
#ifdef ALFA
    let rim = 1.0 - smoothstep(0.55, 1.0, r2);
    let f = mix(rim, exp(-params.soft * r2), 0.6);
    return vec4<f32>(in.color.rgb, clamp(in.color.a * params.gain * f, 0.0, 1.0));
#else
    let glow = exp(-params.soft * r2) * (1.0 - r2);
    let core = params.core * exp(-14.0 * r2);
    let a = in.color.a * params.gain * (glow + core);
    return vec4<f32>(in.color.rgb * a, 0.0);
#endif
}

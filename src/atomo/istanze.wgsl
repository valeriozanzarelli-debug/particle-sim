// Instanced soft balls (istanze.rs): one shared quad, one 32-byte instance per ball.
// Same look as balls.wgsl (the mesh path), so the two paths can be compared (exam R2).

#import bevy_pbr::view_transformations::position_world_to_clip
#import bevy_pbr::mesh_view_bindings::view

struct Parametri {
    size: f32,
    gain: f32,
    soft: f32,
    core: f32,
    offset: vec4<f32>,
};

@group(2) @binding(0) var<uniform> params: Parametri;

struct VertexIn {
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

@vertex
fn vertex(v: VertexIn) -> VertexOut {
    var out: VertexOut;
    let right = view.world_from_view[0].xyz;
    let up = view.world_from_view[1].xyz;
    let s = params.size * v.pos_size.w;
    let a = vec2<f32>(dot(v.allunga, right), dot(v.allunga, up));
    let la = length(a);
    var d = vec2<f32>(1.0, 0.0);
    if (la > 1e-4) {
        d = a / la;
    }
    let nrm = vec2<f32>(-d.y, d.x);
    let off = d * (v.corner.x * (s + la)) + nrm * (v.corner.y * s);
    out.clip = position_world_to_clip(v.pos_size.xyz + params.offset.xyz + right * off.x + up * off.y);
    out.corner = v.corner.xy;
    out.color = v.colore;
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

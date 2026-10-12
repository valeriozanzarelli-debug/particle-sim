// Round soft balls of the carbon atom viewer (render.rs): one camera-facing quad per ball, additive glow.
// Vertex data: position, colour (rgb, alpha = brightness), size factor, corner of the quad.

#import bevy_pbr::mesh_functions::get_world_from_local
#import bevy_pbr::view_transformations::position_world_to_clip
#import bevy_pbr::mesh_view_bindings::view

struct BallParams {
    size: f32,
    gain: f32,
    soft: f32,
    core: f32,
};

@group(2) @binding(0) var<uniform> params: BallParams;

struct VertexIn {
    @builtin(instance_index) instance_index: u32,
    @location(0) position: vec3<f32>,
    @location(1) extra: vec4<f32>,
    @location(2) color: vec4<f32>,
    @location(3) corner: vec2<f32>,
};

struct VertexOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) corner: vec2<f32>,
    @location(1) color: vec4<f32>,
};

@vertex
fn vertex(v: VertexIn) -> VertexOut {
    var out: VertexOut;
    let m = get_world_from_local(v.instance_index);
    let world = (m * vec4<f32>(v.position, 1.0)).xyz;
    let right = view.world_from_view[0].xyz;
    let up = view.world_from_view[1].xyz;
    let s = params.size * v.extra.x;
    let q = world + (right * v.corner.x + up * v.corner.y) * s;
    out.clip = position_world_to_clip(q);
    out.corner = v.corner;
    out.color = v.color;
    return out;
}

@fragment
fn fragment(in: VertexOut) -> @location(0) vec4<f32> {
    let r2 = dot(in.corner, in.corner);
    if (r2 > 1.0) {
        discard;
    }
    // soft glow plus a brighter core
    let glow = exp(-params.soft * r2) * (1.0 - r2);
    let core = params.core * exp(-14.0 * r2);
    let a = in.color.a * params.gain * (glow + core);
    return vec4<f32>(in.color.rgb * a, 0.0);
}

// Round soft balls of the Compton viewer (render.rs).
//
// Each ball is a camera-facing quad. A ball flies straight: position = origin + vel · (clock − t0),
// drawn while t0 ≤ clock ≤ t1. A static ball has vel = 0. The glow of a surfer uses its entity's
// transform (use_transform > 0.5) and is drawn as a star.

#import bevy_pbr::mesh_functions::get_world_from_local
#import bevy_pbr::view_transformations::position_world_to_clip
#import bevy_pbr::mesh_view_bindings::view

struct BallParams {
    clock: f32,
    size: f32,
    gain: f32,
    soft: f32,
    star: f32,
    use_transform: f32,
    pad0: f32,
    pad1: f32,
};

@group(2) @binding(0) var<uniform> params: BallParams;

struct VertexIn {
    @builtin(instance_index) instance_index: u32,
    @location(0) position: vec3<f32>,
    // velocity (xyz), unused (w)
    @location(1) vel: vec4<f32>,
    // t0, t1, size factor, unused
    @location(2) life: vec4<f32>,
    @location(3) color: vec4<f32>,
    @location(4) corner: vec2<f32>,
};

struct VertexOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) corner: vec2<f32>,
    @location(1) color: vec4<f32>,
};

@vertex
fn vertex(v: VertexIn) -> VertexOut {
    var out: VertexOut;
    let s = params.clock - v.life.x;
    let hidden = s < 0.0 || params.clock > v.life.y;
    var world = v.position + v.vel.xyz * max(s, 0.0);
    if (params.use_transform > 0.5) {
        world = (get_world_from_local(v.instance_index) * vec4<f32>(world, 1.0)).xyz;
    }
    let right = view.world_from_view[0].xyz;
    let up = view.world_from_view[1].xyz;
    let q = world + (right * v.corner.x + up * v.corner.y) * (params.size * v.life.z);
    out.clip = position_world_to_clip(q);
    if (hidden) {
        out.clip = vec4<f32>(0.0, 0.0, -10.0, 1.0);
    }
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
    let rim = 1.0 - smoothstep(0.55, 1.0, r2);
    var f = mix(rim, exp(-params.soft * r2), 0.6);
    if (params.star > 0.5) {
        let c = abs(in.corner);
        f = exp(-params.soft * r2) + 0.8 * (exp(-260.0 * c.y * c.y) * (1.0 - c.x) + exp(-260.0 * c.x * c.x) * (1.0 - c.y));
    }
    let a = in.color.a * params.gain * f;
    // additive (Bevy's AlphaMode::Add: premultiplied colour, alpha 0)
    return vec4<f32>(in.color.rgb * a, 0.0);
}

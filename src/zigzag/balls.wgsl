// Round soft balls of the zig-zag viewer (render.rs).
//
// Each ball is a camera-facing quad. Colour: either the vertex colour, or (phase.y > 0.5) a hue given by
// the ball's phase minus the global clock phase (the atom's cloud: the phase of the wave turns in time).

#import bevy_pbr::mesh_functions::get_world_from_local
#import bevy_pbr::view_transformations::position_world_to_clip
#import bevy_pbr::mesh_view_bindings::view

struct BallParams {
    clock_phase: f32,
    size: f32,
    gain: f32,
    soft: f32,
};

@group(2) @binding(0) var<uniform> params: BallParams;

struct VertexIn {
    @builtin(instance_index) instance_index: u32,
    @location(0) position: vec3<f32>,
    // phase, use phase colour, size factor, unused
    @location(1) phase: vec4<f32>,
    @location(2) color: vec4<f32>,
    @location(3) corner: vec2<f32>,
};

struct VertexOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) corner: vec2<f32>,
    @location(1) color: vec4<f32>,
};

fn hue(h: f32) -> vec3<f32> {
    let t = 6.2831853 * h;
    return 0.55 + 0.45 * vec3<f32>(cos(t), cos(t - 2.0943951), cos(t + 2.0943951));
}

@vertex
fn vertex(v: VertexIn) -> VertexOut {
    var out: VertexOut;
    let world = (get_world_from_local(v.instance_index) * vec4<f32>(v.position, 1.0)).xyz;
    let right = view.world_from_view[0].xyz;
    let up = view.world_from_view[1].xyz;
    let q = world + (right * v.corner.x + up * v.corner.y) * (params.size * v.phase.z);
    out.clip = position_world_to_clip(q);
    out.corner = v.corner;
    var c = v.color;
    if (v.phase.y > 0.5) {
        c = vec4<f32>(hue((v.phase.x - params.clock_phase) / 6.2831853), v.color.a);
    }
    out.color = c;
    return out;
}

@fragment
fn fragment(in: VertexOut) -> @location(0) vec4<f32> {
    let r2 = dot(in.corner, in.corner);
    if (r2 > 1.0) {
        discard;
    }
    let rim = 1.0 - smoothstep(0.55, 1.0, r2);
    let f = mix(rim, exp(-params.soft * r2), 0.6);
    let a = in.color.a * params.gain * f;
    // additive (Bevy's AlphaMode::Add: premultiplied colour, alpha 0)
    return vec4<f32>(in.color.rgb * a, 0.0);
}

// Round soft balls in 3-D (render3d.rs).
//
// Each ball is a camera-facing quad. A wave ball flies straight at c from its
// origin: position = origin + dir · (clock − t_emit − depth), drawn while it
// has been emitted and has not passed its kill plane. A static dot (a click,
// the photon's glow) has dir = 0.

#import bevy_pbr::mesh_functions::get_world_from_local
#import bevy_pbr::view_transformations::position_world_to_clip
#import bevy_pbr::mesh_view_bindings::view

struct BallParams {
    clock: f32,
    size: f32,
    lim_x: f32,
    lim_y: f32,
    gain: f32,
    soft: f32,
    star: f32,
    additive: f32,
    use_transform: f32,
    pad0: f32,
    pad1: f32,
    pad2: f32,
};

@group(2) @binding(0) var<uniform> params: BallParams;

struct VertexIn {
    @builtin(instance_index) instance_index: u32,
    @location(0) position: vec3<f32>,
    // direction of flight (xyz) and the visual depth behind the front (w)
    @location(1) dir: vec4<f32>,
    // t_emit, z_kill, size factor, unused
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
    let s = params.clock - v.life.x - v.dir.w;
    let p = v.position + v.dir.xyz * max(s, 0.0);
    let hidden = s < 0.0 || p.z > v.life.y || abs(p.x) > params.lim_x || abs(p.y) > params.lim_y;
    var world = p;
    if (params.use_transform > 0.5) {
        world = (get_world_from_local(v.instance_index) * vec4<f32>(p, 1.0)).xyz;
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
    // soft round ball: bright core, smooth rim
    let rim = 1.0 - smoothstep(0.55, 1.0, r2);
    var f = mix(rim, exp(-params.soft * r2), 0.6);
    if (params.star > 0.5) {
        // a star: soft core + four thin spikes
        let c = abs(in.corner);
        f = exp(-params.soft * r2) + 0.8 * (exp(-260.0 * c.y * c.y) * (1.0 - c.x) + exp(-260.0 * c.x * c.x) * (1.0 - c.y));
    }
    let a = in.color.a * params.gain * f;
    if (params.additive > 0.5) {
        // additive (Bevy's AlphaMode::Add: premultiplied colour, alpha 0)
        return vec4<f32>(in.color.rgb * a, 0.0);
    }
    // alpha blending: the ball's own colour, soft edge
    return vec4<f32>(in.color.rgb, clamp(a, 0.0, 1.0));
}

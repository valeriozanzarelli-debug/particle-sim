//! Drawing: EVERYTHING is a ball (no lines, no wireframes). Own small shader `balls.wgsl`: one camera-facing
//! quad per ball, additive glow. The scenes build their ball lists (fotscena.rs, zzscena.rs) and
//! render_zz.rs puts them on screen.
//!
//! Meshes are REPLACED, never modified in place (Bevy 0.15 leaks a modified mesh; see `swap_mesh` in
//! src/photon/render.rs). Physics (x, y, z) → screen (x, z, −y): physics z is up on screen.

use bevy::{
    asset::load_internal_asset,
    pbr::{MaterialPipeline, MaterialPipelineKey},
    prelude::*,
    render::{
        mesh::{Indices, MeshVertexAttribute, MeshVertexBufferLayoutRef, PrimitiveTopology},
        render_asset::RenderAssetUsages,
        render_resource::{AsBindGroup, RenderPipelineDescriptor, ShaderRef, SpecializedMeshPipelineError, VertexFormat},
    },
};
use std::f64::consts::TAU;

pub const BALL_SHADER: Handle<Shader> = Handle::weak_from_u128(0x4c2a_91d3_70be_4f18_a2c6_33e1_0d97_5b21);
pub const ATTR_EXTRA: MeshVertexAttribute = MeshVertexAttribute::new("At_Extra", 0x4c2a_0001, VertexFormat::Float32x4);

#[allow(dead_code)]
mod params {
    use bevy::render::render_resource::ShaderType;
    #[derive(Clone, Copy, Debug, ShaderType)]
    pub struct BallParams {
        pub size: f32,
        pub gain: f32,
        pub soft: f32,
        pub core: f32,
        /// > 0.5: alpha blending (each ball keeps its own colour, as in the photon laboratory), else additive glow
        pub blend: f32,
    }
}
pub use params::BallParams;

#[derive(Asset, TypePath, AsBindGroup, Clone, Debug)]
pub struct BallMaterial {
    #[uniform(0)]
    pub params: BallParams,
    /// alpha blending instead of additive glow
    pub blend: bool,
}

impl Material for BallMaterial {
    fn vertex_shader() -> ShaderRef { BALL_SHADER.into() }
    fn fragment_shader() -> ShaderRef { BALL_SHADER.into() }
    fn alpha_mode(&self) -> AlphaMode { if self.blend { AlphaMode::Blend } else { AlphaMode::Add } }
    fn specialize(
        _pipeline: &MaterialPipeline<Self>,
        descriptor: &mut RenderPipelineDescriptor,
        layout: &MeshVertexBufferLayoutRef,
        _key: MaterialPipelineKey<Self>,
    ) -> Result<(), SpecializedMeshPipelineError> {
        let vl = layout.0.get_layout(&[
            Mesh::ATTRIBUTE_POSITION.at_shader_location(0),
            ATTR_EXTRA.at_shader_location(1),
            Mesh::ATTRIBUTE_COLOR.at_shader_location(2),
            Mesh::ATTRIBUTE_UV_0.at_shader_location(3),
        ])?;
        descriptor.vertex.buffers = vec![vl];
        descriptor.primitive.cull_mode = None;
        Ok(())
    }
}

pub struct BallsPlugin;

impl Plugin for BallsPlugin {
    fn build(&self, app: &mut App) {
        load_internal_asset!(app, BALL_SHADER, "balls.wgsl", Shader::from_wgsl);
        app.add_plugins(MaterialPlugin::<BallMaterial> { prepass_enabled: false, shadows_enabled: false, ..default() });
    }
}

/// One ball: position, size factor, colour (alpha = brightness).
pub type Ball = ([f32; 3], f32, [f32; 4]);

pub fn quads(items: &[Ball]) -> Mesh { quads_s(items, &[]) }

/// Like `quads`, with a STRETCH per ball (screen-space vector, scene units): the ball becomes a soft ellipse
/// elongated along it (half-length added to the radius). Missing entries = round balls.
pub fn quads_s(items: &[Ball], allunga: &[[f32; 3]]) -> Mesh {
    let hidden = [([0.0f32, -1.0e9, 0.0], 0.0f32, [0.0f32; 4])];
    let items = if items.is_empty() { &hidden[..] } else { items };
    let n = items.len();
    let mut pos = Vec::with_capacity(4 * n);
    let mut ex = Vec::with_capacity(4 * n);
    let mut col = Vec::with_capacity(4 * n);
    let mut uv = Vec::with_capacity(4 * n);
    let mut idx: Vec<u32> = Vec::with_capacity(6 * n);
    const CORNERS: [[f32; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    for (q, it) in items.iter().enumerate() {
        let st = allunga.get(q).copied().unwrap_or([0.0; 3]);
        for c in CORNERS {
            pos.push(it.0);
            ex.push([it.1, st[0], st[1], st[2]]);
            col.push(it.2);
            uv.push(c);
        }
        let b = 4 * q as u32;
        idx.extend_from_slice(&[b, b + 1, b + 2, b, b + 2, b + 3]);
    }
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    mesh.insert_attribute(ATTR_EXTRA, ex);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, col);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uv);
    mesh.insert_indices(Indices::U32(idx));
    mesh
}

/// Physics → screen.
#[inline]
pub fn schermo(p: [f64; 3]) -> Vec3 { Vec3::new(p[0] as f32, p[2] as f32, -p[1] as f32) }

/// Cyclic phase palette of the photon laboratory (src/photon/render.rs).
pub fn phase_color(angle: f64, alpha: f32) -> [f32; 4] {
    let t = (angle.rem_euclid(TAU) / TAU) as f32;
    const STOPS: [[f32; 3]; 5] = [[0.20, 0.95, 1.00], [0.25, 0.45, 1.00], [0.70, 0.30, 1.00], [1.00, 0.35, 0.70], [1.00, 0.85, 0.30]];
    let x = t * 5.0;
    let i = (x.floor() as usize) % 5;
    let s = x - x.floor();
    let (a, b) = (STOPS[i], STOPS[(i + 1) % 5]);
    [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s, alpha]
}


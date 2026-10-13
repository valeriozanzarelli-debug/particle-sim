//! GPU INSTANCED BALLS (study MOTORE GRAFICO, docs/MOTORE_GRAFICO.md, exam R2).
//!
//! The old path (render.rs `quads`) rebuilds a mesh of 4 vertices + 6 indices per ball every frame
//! (~208 bytes per ball). Here one ball is ONE INSTANCE of 32 bytes (position + size, stretch, colour as 4 bytes)
//! drawn on a single shared quad: the CPU only fills a flat array, the GPU expands it.
//!
//! Positions are baked in world units by the CPU (translation of the layer added), so the shader needs only
//! the view. Two blend modes: alpha (the energy balls, each keeps its own colour) and additive glow.
//! Only balls: no lines, no surfaces (exam R5).

use bevy::{
    asset::load_internal_asset,
    core_pipeline::core_3d::Transparent3d,
    ecs::{
        query::QueryItem,
        system::{lifetimeless::*, SystemParamItem},
    },
    pbr::{MeshPipeline, MeshPipelineKey, RenderMeshInstances, SetMeshBindGroup, SetMeshViewBindGroup},
    prelude::*,
    render::{
        extract_component::{ExtractComponent, ExtractComponentPlugin},
        mesh::{allocator::MeshAllocator, Indices, MeshVertexBufferLayoutRef, PrimitiveTopology, RenderMesh, RenderMeshBufferInfo},
        render_asset::{RenderAssetUsages, RenderAssets},
        render_phase::{
            AddRenderCommand, DrawFunctions, PhaseItem, PhaseItemExtraIndex, RenderCommand, RenderCommandResult,
            SetItemPipeline, TrackedRenderPass, ViewSortedRenderPhases,
        },
        render_resource::*,
        renderer::RenderDevice,
        sync_world::MainEntity,
        view::ExtractedView,
        Render, RenderApp, RenderSet,
    },
};
use bytemuck::{Pod, Zeroable};
use std::sync::Arc;

use crate::render::Ball;

pub const ISTANZE_SHADER: Handle<Shader> = Handle::weak_from_u128(0x7d13_5a0e_c4b2_4e61_9f08_2b6d_41aa_93c7);

/// One ball on the GPU: 32 bytes.
#[derive(Clone, Copy, Pod, Zeroable, Debug, Default)]
#[repr(C)]
pub struct Istanza {
    /// world position (xyz) and size factor (w)
    pub pos_size: [f32; 4],
    /// stretch vector in world units (the ball becomes an ellipse along its screen projection)
    pub allunga: [f32; 3],
    /// colour rgba, 8 bits each (alpha = brightness)
    pub colore: u32,
}

#[inline]
fn u8c(x: f32) -> u32 { (x.clamp(0.0, 1.0) * 255.0 + 0.5) as u32 }

/// Pack a list of balls (scene units) into instances.
pub fn impacca(items: &[Ball], allunga: &[[f32; 3]]) -> Vec<Istanza> {
    items.iter().enumerate().map(|(i, b)| {
        let c = b.2;
        Istanza {
            pos_size: [b.0[0], b.0[1], b.0[2], b.1],
            allunga: allunga.get(i).copied().unwrap_or([0.0; 3]),
            colore: u8c(c[0]) | (u8c(c[1]) << 8) | (u8c(c[2]) << 16) | (u8c(c[3]) << 24),
        }
    }).collect()
}

/// Instances already on the GPU (written by a compute shader, e.g. the photon's energy balls in gpu_fotone.rs):
/// drawn straight from this buffer, no CPU round trip.
#[derive(Clone)]
pub struct IstanzeGpu { pub buffer: Buffer, pub n: u32 }

/// A layer of instanced balls. `alfa` = alpha blending (else additive glow). The parameters of the old
/// BallMaterial: size, gain, soft, core.
#[derive(Component, Clone)]
pub struct StratoIstanze {
    pub dati: Arc<Vec<Istanza>>,
    pub alfa: bool,
    pub size: f32,
    pub gain: f32,
    pub soft: f32,
    pub core: f32,
    /// translation of the layer (world), applied in the shader so the data need not change when the camera follows
    pub offset: [f32; 3],
    /// when set, the instances come from this GPU buffer instead of `dati`
    pub gpu: Option<IstanzeGpu>,
}

impl Default for StratoIstanze {
    fn default() -> Self { StratoIstanze { dati: Arc::new(vec![]), alfa: false, size: 0.1, gain: 1.0, soft: 3.0, core: 0.3, offset: [0.0; 3], gpu: None } }
}

impl ExtractComponent for StratoIstanze {
    type QueryData = &'static StratoIstanze;
    type QueryFilter = ();
    type Out = Self;
    fn extract_component(item: QueryItem<'_, Self::QueryData>) -> Option<Self> { Some(item.clone()) }
}

/// The shared quad: 4 corners, position only (the corner in x, y).
pub fn quadrato() -> Mesh {
    let mut m = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    m.insert_attribute(Mesh::ATTRIBUTE_POSITION, vec![[-1.0f32, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]]);
    m.insert_indices(Indices::U32(vec![0, 1, 2, 0, 2, 3]));
    m
}

pub struct IstanzePlugin;

impl Plugin for IstanzePlugin {
    fn build(&self, app: &mut App) {
        load_internal_asset!(app, ISTANZE_SHADER, "istanze.wgsl", Shader::from_wgsl);
        app.add_plugins(ExtractComponentPlugin::<StratoIstanze>::default());
        app.sub_app_mut(RenderApp)
            .add_render_command::<Transparent3d, DisegnaIstanze>()
            .init_resource::<SpecializedMeshPipelines<PipelineIstanze>>()
            .add_systems(Render, (
                accoda.in_set(RenderSet::QueueMeshes),
                prepara_buffer.in_set(RenderSet::PrepareResources),
            ));
    }
    fn finish(&self, app: &mut App) {
        app.sub_app_mut(RenderApp).init_resource::<PipelineIstanze>();
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Chiave { mesh: MeshPipelineKey, alfa: bool }

#[allow(clippy::too_many_arguments)]
fn accoda(
    draw_functions: Res<DrawFunctions<Transparent3d>>,
    pipeline: Res<PipelineIstanze>,
    mut pipelines: ResMut<SpecializedMeshPipelines<PipelineIstanze>>,
    cache: Res<PipelineCache>,
    meshes: Res<RenderAssets<RenderMesh>>,
    mesh_instances: Res<RenderMeshInstances>,
    strati: Query<(Entity, &MainEntity, &StratoIstanze)>,
    mut phases: ResMut<ViewSortedRenderPhases<Transparent3d>>,
    views: Query<(Entity, &ExtractedView, &Msaa)>,
) {
    let draw = draw_functions.read().id::<DisegnaIstanze>();
    for (view_entity, view, msaa) in &views {
        let Some(phase) = phases.get_mut(&view_entity) else { continue };
        let view_key = MeshPipelineKey::from_msaa_samples(msaa.samples()) | MeshPipelineKey::from_hdr(view.hdr);
        let rangefinder = view.rangefinder3d();
        for (entity, main_entity, strato) in &strati {
            let Some(mi) = mesh_instances.render_mesh_queue_data(*main_entity) else { continue };
            let Some(mesh) = meshes.get(mi.mesh_asset_id) else { continue };
            let blend = if strato.alfa { MeshPipelineKey::BLEND_ALPHA } else { MeshPipelineKey::BLEND_PREMULTIPLIED_ALPHA };
            let key = Chiave { mesh: view_key | blend | MeshPipelineKey::from_primitive_topology(mesh.primitive_topology()), alfa: strato.alfa };
            let Ok(pid) = pipelines.specialize(&cache, &pipeline, key, &mesh.layout) else { continue };
            // the alpha layer is drawn first (behind the glow), as in the mesh path
            phase.add(Transparent3d {
                entity: (entity, *main_entity),
                pipeline: pid,
                draw_function: draw,
                distance: rangefinder.distance_translation(&mi.translation) + if strato.alfa { 1.0e6 } else { 0.0 },
                batch_range: 0..1,
                extra_index: PhaseItemExtraIndex::NONE,
            });
        }
    }
}

#[derive(Component)]
pub struct BufferIstanze { buffer: Buffer, params: BindGroup, n: u32 }

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Parametri { size: f32, gain: f32, soft: f32, core: f32, offset: [f32; 4] }

fn prepara_buffer(
    mut commands: Commands,
    q: Query<(Entity, &MainEntity, &StratoIstanze)>,
    device: Res<RenderDevice>,
    pipeline: Res<PipelineIstanze>,
    mut cache: Local<std::collections::HashMap<MainEntity, (usize, Buffer, u32)>>,
) {
    for (entity, main, s) in &q {
        // the instance buffer is rebuilt only when the data change (same Arc = same balls)
        let chiave = Arc::as_ptr(&s.dati) as usize;
        let (buffer, n) = if let Some(g) = s.gpu.as_ref().filter(|g| g.n > 0) { (g.buffer.clone(), g.n) } else { match cache.get(main) {
            Some((k, b, n)) if *k == chiave => (b.clone(), *n),
            _ => {
                let vuoto = [Istanza { pos_size: [0.0, -1.0e9, 0.0, 0.0], ..default() }];
                let dati: &[Istanza] = if s.dati.is_empty() { &vuoto } else { &s.dati };
                let b = device.create_buffer_with_data(&BufferInitDescriptor {
                    label: Some("palline a istanze"),
                    contents: bytemuck::cast_slice(dati),
                    usage: BufferUsages::VERTEX,
                });
                cache.insert(*main, (chiave, b.clone(), dati.len() as u32));
                (b, dati.len() as u32)
            }
        } };
        let p = Parametri { size: s.size, gain: s.gain, soft: s.soft, core: s.core, offset: [s.offset[0], s.offset[1], s.offset[2], 0.0] };
        let ub = device.create_buffer_with_data(&BufferInitDescriptor {
            label: Some("parametri palline"),
            contents: bytemuck::bytes_of(&p),
            usage: BufferUsages::UNIFORM,
        });
        let params = device.create_bind_group("parametri palline", &pipeline.layout_params, &BindGroupEntries::single(ub.as_entire_binding()));
        commands.entity(entity).insert(BufferIstanze { buffer, params, n });
    }
}

#[derive(Resource)]
pub struct PipelineIstanze { mesh_pipeline: MeshPipeline, layout_params: BindGroupLayout }

impl FromWorld for PipelineIstanze {
    fn from_world(world: &mut World) -> Self {
        let device = world.resource::<RenderDevice>();
        let layout_params = device.create_bind_group_layout(
            "parametri palline",
            &BindGroupLayoutEntries::single(ShaderStages::VERTEX_FRAGMENT, binding_types::uniform_buffer_sized(false, None)),
        );
        PipelineIstanze { mesh_pipeline: world.resource::<MeshPipeline>().clone(), layout_params }
    }
}

impl SpecializedMeshPipeline for PipelineIstanze {
    type Key = Chiave;
    fn specialize(&self, key: Chiave, layout: &MeshVertexBufferLayoutRef) -> Result<RenderPipelineDescriptor, SpecializedMeshPipelineError> {
        let mut d = self.mesh_pipeline.specialize(key.mesh, layout)?;
        d.label = Some("palline a istanze".into());
        d.vertex.shader = ISTANZE_SHADER;
        d.vertex.buffers.push(VertexBufferLayout {
            array_stride: std::mem::size_of::<Istanza>() as u64,
            step_mode: VertexStepMode::Instance,
            attributes: vec![
                VertexAttribute { format: VertexFormat::Float32x4, offset: 0, shader_location: 3 },
                VertexAttribute { format: VertexFormat::Float32x3, offset: 16, shader_location: 4 },
                VertexAttribute { format: VertexFormat::Unorm8x4, offset: 28, shader_location: 5 },
            ],
        });
        let f = d.fragment.as_mut().unwrap();
        f.shader = ISTANZE_SHADER;
        if key.alfa { f.shader_defs.push("ALFA".into()); }
        d.vertex.shader_defs.retain(|_| true);
        d.layout.push(self.layout_params.clone());
        d.primitive.cull_mode = None;
        if let Some(ds) = d.depth_stencil.as_mut() { ds.depth_write_enabled = false; }
        Ok(d)
    }
}

type DisegnaIstanze = (SetItemPipeline, SetMeshViewBindGroup<0>, SetMeshBindGroup<1>, SetParametri<2>, DisegnaQuadrati);

struct SetParametri<const I: usize>;
impl<P: PhaseItem, const I: usize> RenderCommand<P> for SetParametri<I> {
    type Param = ();
    type ViewQuery = ();
    type ItemQuery = Read<BufferIstanze>;
    fn render<'w>(_: &P, _: (), b: Option<&'w BufferIstanze>, _: SystemParamItem<'w, '_, Self::Param>, pass: &mut TrackedRenderPass<'w>) -> RenderCommandResult {
        let Some(b) = b else { return RenderCommandResult::Skip };
        pass.set_bind_group(I, &b.params, &[]);
        RenderCommandResult::Success
    }
}

struct DisegnaQuadrati;
impl<P: PhaseItem> RenderCommand<P> for DisegnaQuadrati {
    type Param = (SRes<RenderAssets<RenderMesh>>, SRes<RenderMeshInstances>, SRes<MeshAllocator>);
    type ViewQuery = ();
    type ItemQuery = Read<BufferIstanze>;
    fn render<'w>(item: &P, _: (), b: Option<&'w BufferIstanze>, (meshes, mesh_instances, alloc): SystemParamItem<'w, '_, Self::Param>, pass: &mut TrackedRenderPass<'w>) -> RenderCommandResult {
        let alloc = alloc.into_inner();
        let Some(mi) = mesh_instances.render_mesh_queue_data(item.main_entity()) else { return RenderCommandResult::Skip };
        let Some(gpu_mesh) = meshes.into_inner().get(mi.mesh_asset_id) else { return RenderCommandResult::Skip };
        let Some(b) = b else { return RenderCommandResult::Skip };
        let Some(vs) = alloc.mesh_vertex_slice(&mi.mesh_asset_id) else { return RenderCommandResult::Skip };
        pass.set_vertex_buffer(0, vs.buffer.slice(..));
        pass.set_vertex_buffer(1, b.buffer.slice(..));
        match &gpu_mesh.buffer_info {
            RenderMeshBufferInfo::Indexed { index_format, count } => {
                let Some(is) = alloc.mesh_index_slice(&mi.mesh_asset_id) else { return RenderCommandResult::Skip };
                pass.set_index_buffer(is.buffer.slice(..), 0, *index_format);
                pass.draw_indexed(is.range.start..(is.range.start + count), vs.range.start as i32, 0..b.n);
            }
            RenderMeshBufferInfo::NonIndexed => pass.draw(vs.range, 0..b.n),
        }
        RenderCommandResult::Success
    }
}

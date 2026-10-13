//! TRASPARENZA PESATA PER PROFONDITÀ delle palline di energia del fotone (studio mg_oit,
//! docs/MOTORE_GRAFICO_OIT.md, pre-registrazione lab/out/mg_oit_prereg.json).
//!
//! Before: the energy balls were alpha-blended in submission order ("over"). Seen along the flight axis tens of
//! thousands of balls of every hue pile up on each pixel: the result was a white fog that hid the spin swirl.
//! Now: weighted blended order-independent transparency (McGuire & Bavoil 2013), the method of the Python
//! prototype lab/fp_render.py:
//!   pass 0 (offscreen, depth only): the nearest bright ball core on each pixel = the FRONT of the light there;
//!   pass 1 (offscreen, one draw call): every ball adds (c·a·w, a·w) to ACCUM (Rgba16Float) and multiplies
//!          REVEAL (R16Float) by exp(−k a); w = exp(−(z − z_fronte)/D) favours the nearest layer of light;
//!   pass 2 (full screen): colour = Σ c a w / Σ a w, coverage = 1 − exp(−k Σ a), premultiplied over the scene,
//!          before the additive glow layers (surfers, box) of the transparent pass.
//! The exposure k comes from the expected column density of the balls (no read-back, same rule on every preset);
//! the tail of each ball is weighted by the view: from the side only its newest quarter, along the axis all
//! of it (2 λ/c): the arcs of the swirl. Everything drawn is a ball at a computed position: nothing is painted.
//!
//! Why not Bevy 0.15's built-in OIT: it stores per-pixel linked lists in storage buffers (fragment-writable
//! storage, ~8 layers × 8 bytes per pixel) and is wired to the mesh materials; ours are instanced quads, and the
//! weighted blend costs three small render targets and three passes, which the Intel UHD runs at the same cost.

use bevy::{
    asset::load_internal_asset,
    core_pipeline::{
        core_3d::graph::{Core3d, Node3d},
        fullscreen_vertex_shader::fullscreen_shader_vertex_state,
    },
    ecs::query::QueryItem,
    prelude::*,
    render::{
        extract_resource::{ExtractResource, ExtractResourcePlugin},
        render_graph::{NodeRunError, RenderGraphApp, RenderGraphContext, RenderLabel, ViewNode, ViewNodeRunner},
        render_resource::{binding_types::{texture_2d, texture_depth_2d, uniform_buffer_sized}, *},
        renderer::{RenderContext, RenderDevice},
        texture::{CachedTexture, TextureCache},
        view::{ExtractedView, Msaa, ViewTarget},
        Render, RenderApp, RenderSet,
    },
};
use bytemuck::{Pod, Zeroable};
use std::sync::Arc;

use crate::istanze::{Istanza, IstanzeGpu};

pub const OIT_SHADER: Handle<Shader> = Handle::weak_from_u128(0x3b8e_91d2_5c47_4f0a_8e16_d0c2_7a95_e413);

const ACCUM: TextureFormat = TextureFormat::Rgba16Float;
const REVEAL: TextureFormat = TextureFormat::R16Float;
const FRONTE: TextureFormat = TextureFormat::Depth32Float;

/// What the OIT layer draws this frame (main world, extracted every frame). `attivo` false = layer off.
#[derive(Resource, Clone, ExtractResource)]
pub struct OitPalline {
    pub attivo: bool,
    /// the instances on the CPU path (fotdisegno::palline), or …
    pub dati: Arc<Vec<Istanza>>,
    /// … the instances already on the GPU (gpu_fotone.rs)
    pub gpu: Option<IstanzeGpu>,
    pub size: f32,
    pub gain: f32,
    pub soft: f32,
    /// translation of the layer (world)
    pub offset: [f32; 3],
    /// flight axis (world, unit) and centre of the light (world)
    pub asse: [f32; 3],
    pub centro: [f32; 3],
    /// wavelength and waist (world units), pulse length (energy, world units)
    pub lambda: f32,
    pub w0: f32,
    pub sigma: f32,
    /// energy balls (heads) and tail samples per ball
    pub teste: f32,
    pub ns: f32,
    /// knobs (aesthetic, counted in the doc): depth of the transparency (λ), exposure target, visible tail
    /// fraction from the side, colour gain of the composite
    pub profondita: f32,
    pub esposizione: f32,
    pub coda_lato: f32,
    pub guadagno: f32,
}

impl Default for OitPalline {
    fn default() -> Self {
        OitPalline {
            attivo: false, dati: Arc::new(vec![]), gpu: None, size: 0.1, gain: 1.0, soft: 2.0, offset: [0.0; 3], asse: [1.0, 0.0, 0.0],
            centro: [0.0; 3], lambda: 12.0, w0: 12.0, sigma: 8.5, teste: 1.0, ns: 2.0,
            profondita: DEF_PROFONDITA, esposizione: DEF_ESPOSIZIONE, coda_lato: DEF_CODA_LATO, guadagno: DEF_GUADAGNO,
        }
    }
}

pub const DEF_PROFONDITA: f32 = 0.2;
pub const DEF_ESPOSIZIONE: f32 = 1.5;
pub const DEF_CODA_LATO: f32 = 0.25;
pub const DEF_GUADAGNO: f32 = 1.0;
/// from the side (view ⟂ axis) the depth and the exposure are this many times larger: there the colours do not
/// pile up (each line of sight crosses ~one band), and the bands stay as bright as with the old alpha blending
pub const LATO_PROFONDITA: f32 = 2.5;
pub const LATO_ESPOSIZIONE: f32 = 2.0;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Default)]
struct ParGpu {
    clip_from_world: [[f32; 4]; 4],
    right: [f32; 4],
    up: [f32; 4],
    eye: [f32; 4],
    fwd: [f32; 4],
    offset: [f32; 4],
    a: [f32; 4],
    b: [f32; 4],
}

/// Exposure k: the coverage is 1 − exp(−k Σa); k is chosen so that at the expected peak of the column density
/// of the balls Σa·k = `obiettivo`. The light is taken as a Gaussian cloud (σ along the axis, w0/2 across),
/// projected on the screen plane: peak surface density N / (2π √det Σ₂). Same formula on every preset.
pub fn esposizione(o: &OitPalline, right: Vec3, up: Vec3, allineamento: f32) -> f32 {
    let a = Vec3::from(o.asse).normalize_or_zero();
    let (st, sa) = (0.5 * o.w0, o.sigma);
    let cov = |u: Vec3, v: Vec3| st * st * u.dot(v) + (sa * sa - st * st) * u.dot(a) * v.dot(a);
    let det = (cov(right, right) * cov(up, up) - cov(right, up).powi(2)).max(1e-6);
    // instances per head that count: head + tail weights seen with this alignment (as the shader does)
    let ns = o.ns.max(1.0);
    let mut peso = 1.0;
    for j in 0..(ns as usize) {
        peso += peso_coda(j as f32, ns, allineamento, o.coda_lato);
    }
    // a ball: mean alpha ~0.35 × gain, soft profile ~0.5, area π s (s + stretch ~0.15 λ)
    let area = std::f32::consts::PI * o.size * (o.size + 0.15 * o.lambda);
    let alfa = 0.35 * o.gain.min(1.5) * 0.5;
    let colonna = o.teste * peso * alfa * area / (std::f32::consts::TAU * det.sqrt());
    o.esposizione / colonna.max(1e-6)
}

/// The tail weight of the shader (oit_palline.wgsl `coda`), for the exam V2 and the exposure.
pub fn peso_coda(j: f32, ns: f32, c: f32, coda_lato: f32) -> f32 {
    let ss = |e0: f32, e1: f32, x: f32| { let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0); t * t * (3.0 - 2.0 * t) };
    let q = (j + 1.0) / ns;
    let qmax = coda_lato + (1.0 - coda_lato) * ss(0.55, 0.95, c);
    let vis = (1.0 + (qmax - q) * ns).clamp(0.0, 1.0);
    let lato = 0.5 * (1.0 - j / ns);
    let dietro = 0.7 * (1.0 - 0.6 * q);
    vis * (lato + (dietro - lato) * ss(0.55, 0.95, c))
}

/// Exam V2 (docs/MOTORE_GRAFICO_OIT.md): the swirl as DRAWN. The instances exactly as sent to the GPU (CPU path),
/// heads with their tails, seen from the camera: mean signed angle (rad, + = counter-clockwise on the screen)
/// from the oldest visible tail ball to the head, around the projected flight axis, weighted by the visible alpha.
/// Also the sign of the angular momentum about the axis from the same balls in 3-D (Σ (r − c) × (head − tail)).
pub fn vortice_disegnato(dati: &[Istanza], ns: usize, o: &OitPalline, eye: Vec3, fwd: Vec3, right: Vec3, up: Vec3) -> (f64, f64, usize) {
    let off = Vec3::from(o.offset);
    let c3 = Vec3::from(o.centro);
    let asse = Vec3::from(o.asse).normalize_or_zero();
    let proj = |p: Vec3| { let d = p - eye; let z = d.dot(fwd).max(1e-3); (d.dot(right) / z, d.dot(up) / z) };
    let c = fwd.dot(asse).abs();
    // the axis point: the centre of the light, projected
    let (xc, yc) = proj(c3);
    let alfa = |i: &Istanza| ((i.colore >> 24) & 255) as f32 / 255.0;
    let (mut sw, mut swd, mut l3) = (0.0f64, 0.0f64, 0.0f64);
    let mut n = 0usize;
    for g in dati.chunks_exact(ns + 1) {
        let testa = &g[ns];
        let ph = Vec3::new(testa.pos_size[0], testa.pos_size[1], testa.pos_size[2]) + off;
        // the oldest tail ball that the shader still shows (j = ns − 1 … 0, written first)
        let mut scelta = None;
        for (s, t) in g[..ns].iter().enumerate() {
            let j = (ns - 1 - s) as f32;
            if peso_coda(j, ns as f32, c, o.coda_lato) > 0.05 { scelta = Some(t); break; }
        }
        let Some(t) = scelta else { continue };
        let pt = Vec3::new(t.pos_size[0], t.pos_size[1], t.pos_size[2]) + off;
        let (xh, yh) = proj(ph);
        let (xt, yt) = proj(pt);
        let fh = (yh - yc).atan2(xh - xc);
        let ft = (yt - yc).atan2(xt - xc);
        let mut d = (fh - ft) as f64;
        while d > std::f64::consts::PI { d -= std::f64::consts::TAU; }
        while d < -std::f64::consts::PI { d += std::f64::consts::TAU; }
        let w = alfa(testa) as f64;
        sw += w;
        swd += w * d;
        l3 += ((ph - c3).cross(ph - pt)).dot(asse) as f64;
        n += 1;
    }
    (if sw > 0.0 { swd / sw } else { 0.0 }, l3, n)
}

pub struct OitPlugin;

#[derive(Debug, Hash, PartialEq, Eq, Clone, RenderLabel)]
pub struct OitEtichetta;

impl Plugin for OitPlugin {
    fn build(&self, app: &mut App) {
        load_internal_asset!(app, OIT_SHADER, "oit_palline.wgsl", Shader::from_wgsl);
        app.init_resource::<OitPalline>();
        app.add_plugins(ExtractResourcePlugin::<OitPalline>::default());
        let Some(ra) = app.get_sub_app_mut(RenderApp) else { return };
        ra.init_resource::<SpecializedRenderPipelines<PipelineOit>>()
            .add_systems(Render, prepara.in_set(RenderSet::PrepareBindGroups))
            .add_render_graph_node::<ViewNodeRunner<NodoOit>>(Core3d, OitEtichetta)
            // after the opaque scene, before the transparent pass: the glow layers (surfers, box) stay on top
            .add_render_graph_edges(Core3d, (Node3d::MainTransmissivePass, OitEtichetta, Node3d::MainTransparentPass));
    }
    fn finish(&self, app: &mut App) {
        let Some(ra) = app.get_sub_app_mut(RenderApp) else { return };
        ra.init_resource::<PipelineOit>();
    }
}

#[derive(Resource)]
pub struct PipelineOit { layout_fronte: BindGroupLayout, layout_disegno: BindGroupLayout, layout_comp: BindGroupLayout, fronte: CachedRenderPipelineId, disegno: CachedRenderPipelineId }

impl FromWorld for PipelineOit {
    fn from_world(world: &mut World) -> Self {
        let device = world.resource::<RenderDevice>();
        let layout_fronte = device.create_bind_group_layout("oit fronte",
            &BindGroupLayoutEntries::single(ShaderStages::VERTEX_FRAGMENT, uniform_buffer_sized(false, None)));
        let layout_disegno = device.create_bind_group_layout("oit palline",
            &BindGroupLayoutEntries::with_indices(ShaderStages::VERTEX_FRAGMENT, ((0, uniform_buffer_sized(false, None)), (3, texture_depth_2d()))));
        let layout_comp = device.create_bind_group_layout("oit composizione",
            &BindGroupLayoutEntries::sequential(ShaderStages::VERTEX_FRAGMENT, (
                uniform_buffer_sized(false, None),
                texture_2d(TextureSampleType::Float { filterable: false }),
                texture_2d(TextureSampleType::Float { filterable: false }),
            )));
        let somma = BlendComponent { src_factor: BlendFactor::One, dst_factor: BlendFactor::One, operation: BlendOperation::Add };
        let per = BlendComponent { src_factor: BlendFactor::Zero, dst_factor: BlendFactor::Src, operation: BlendOperation::Add };
        let istanze = || vec![VertexBufferLayout {
            array_stride: std::mem::size_of::<Istanza>() as u64,
            step_mode: VertexStepMode::Instance,
            attributes: vec![
                VertexAttribute { format: VertexFormat::Float32x4, offset: 0, shader_location: 0 },
                VertexAttribute { format: VertexFormat::Float32x3, offset: 16, shader_location: 1 },
                VertexAttribute { format: VertexFormat::Unorm8x4, offset: 28, shader_location: 2 },
            ],
        }];
        // pass 0: depth only, the front of the light on each pixel
        let desc_fronte = RenderPipelineDescriptor {
            label: Some("oit fronte".into()),
            layout: vec![layout_fronte.clone()],
            push_constant_ranges: vec![],
            vertex: VertexState { shader: OIT_SHADER, shader_defs: vec![], entry_point: "vertex_fronte".into(), buffers: istanze() },
            primitive: PrimitiveState { cull_mode: None, ..default() },
            depth_stencil: Some(DepthStencilState {
                format: FRONTE, depth_write_enabled: true, depth_compare: CompareFunction::Less,
                stencil: StencilState::default(), bias: DepthBiasState::default(),
            }),
            multisample: MultisampleState::default(),
            fragment: Some(FragmentState { shader: OIT_SHADER, shader_defs: vec![], entry_point: "fragment_fronte".into(), targets: vec![] }),
            zero_initialize_workgroup_memory: false,
        };
        let desc = RenderPipelineDescriptor {
            label: Some("oit palline".into()),
            layout: vec![layout_disegno.clone()],
            push_constant_ranges: vec![],
            vertex: VertexState { shader: OIT_SHADER, shader_defs: vec![], entry_point: "vertex".into(), buffers: istanze() },
            primitive: PrimitiveState { cull_mode: None, ..default() },
            depth_stencil: None,
            multisample: MultisampleState::default(),
            fragment: Some(FragmentState {
                shader: OIT_SHADER, shader_defs: vec![], entry_point: "fragment".into(),
                targets: vec![
                    Some(ColorTargetState { format: ACCUM, blend: Some(BlendState { color: somma, alpha: somma }), write_mask: ColorWrites::ALL }),
                    Some(ColorTargetState { format: REVEAL, blend: Some(BlendState { color: per, alpha: per }), write_mask: ColorWrites::ALL }),
                ],
            }),
            zero_initialize_workgroup_memory: false,
        };
        let cache = world.resource::<PipelineCache>();
        let fronte = cache.queue_render_pipeline(desc_fronte);
        let disegno = cache.queue_render_pipeline(desc);
        PipelineOit { layout_fronte, layout_disegno, layout_comp, fronte, disegno }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct ChiaveComp { formato: TextureFormat, campioni: u32 }

impl SpecializedRenderPipeline for PipelineOit {
    type Key = ChiaveComp;
    fn specialize(&self, key: ChiaveComp) -> RenderPipelineDescriptor {
        // premultiplied "over": the weighted mean colour × coverage covers the scene by the coverage
        let over = BlendComponent { src_factor: BlendFactor::One, dst_factor: BlendFactor::OneMinusSrcAlpha, operation: BlendOperation::Add };
        RenderPipelineDescriptor {
            label: Some("oit composizione".into()),
            layout: vec![self.layout_comp.clone()],
            push_constant_ranges: vec![],
            vertex: fullscreen_shader_vertex_state(),
            primitive: PrimitiveState::default(),
            depth_stencil: None,
            multisample: MultisampleState { count: key.campioni, ..default() },
            fragment: Some(FragmentState {
                shader: OIT_SHADER, shader_defs: vec![], entry_point: "composita".into(),
                targets: vec![Some(ColorTargetState { format: key.formato, blend: Some(BlendState { color: over, alpha: over }), write_mask: ColorWrites::ALL })],
            }),
            zero_initialize_workgroup_memory: false,
        }
    }
}

/// Per view: the two offscreen targets, the bind groups, the composite pipeline.
#[derive(Component)]
pub struct VistaOit {
    fronte: CachedTexture,
    bg_fronte: BindGroup,
    accum: CachedTexture,
    reveal: CachedTexture,
    bg_disegno: BindGroup,
    bg_comp: BindGroup,
    comp: CachedRenderPipelineId,
    buffer: Buffer,
    n: u32,
}

#[allow(clippy::too_many_arguments)]
fn prepara(
    mut commands: Commands,
    o: Option<Res<OitPalline>>,
    device: Res<RenderDevice>,
    pipeline: Res<PipelineOit>,
    mut pipelines: ResMut<SpecializedRenderPipelines<PipelineOit>>,
    cache: Res<PipelineCache>,
    mut texture_cache: ResMut<TextureCache>,
    views: Query<(Entity, &ExtractedView, &ViewTarget, Option<&Msaa>)>,
    mut buf: Local<Option<(usize, Buffer, u32)>>,
) {
    let attivo = o.as_ref().map(|o| o.attivo).unwrap_or(false);
    if !attivo {
        for (e, ..) in &views { commands.entity(e).remove::<VistaOit>(); }
        return;
    }
    let o = o.unwrap();
    // the instance buffer: the GPU one as it is, else rebuilt only when the data change (same Arc = same balls)
    let (buffer, n) = if let Some(g) = o.gpu.as_ref().filter(|g| g.n > 0) { (g.buffer.clone(), g.n) } else {
        let chiave = Arc::as_ptr(&o.dati) as usize;
        match buf.as_ref() {
            Some((k, b, n)) if *k == chiave => (b.clone(), *n),
            _ => {
                let vuoto = [Istanza { pos_size: [0.0, -1.0e9, 0.0, 0.0], ..default() }];
                let dati: &[Istanza] = if o.dati.is_empty() { &vuoto } else { &o.dati };
                let b = device.create_buffer_with_data(&BufferInitDescriptor {
                    label: Some("oit istanze"), contents: bytemuck::cast_slice(dati), usage: BufferUsages::VERTEX,
                });
                *buf = Some((chiave, b.clone(), if o.dati.is_empty() { 0 } else { dati.len() as u32 }));
                (b, if o.dati.is_empty() { 0 } else { dati.len() as u32 })
            }
        }
    };
    for (e, view, target, msaa) in &views {
        let size = target.main_texture().size();
        let tex = |tc: &mut TextureCache, label: &'static str, format: TextureFormat| tc.get(&device, TextureDescriptor {
            label: Some(label), size: Extent3d { width: size.width, height: size.height, depth_or_array_layers: 1 },
            mip_level_count: 1, sample_count: 1, dimension: TextureDimension::D2, format,
            usage: TextureUsages::RENDER_ATTACHMENT | TextureUsages::TEXTURE_BINDING, view_formats: &[],
        });
        let fronte = tex(&mut texture_cache, "oit fronte", FRONTE);
        let accum = tex(&mut texture_cache, "oit accum", ACCUM);
        let reveal = tex(&mut texture_cache, "oit reveal", REVEAL);
        // camera
        let wv = view.world_from_view.compute_matrix();
        let clip_from_world = view.clip_from_world.unwrap_or(view.clip_from_view * wv.inverse());
        let right = wv.x_axis.truncate().normalize_or_zero();
        let up = wv.y_axis.truncate().normalize_or_zero();
        let fwd = -wv.z_axis.truncate().normalize_or_zero();
        let eye = wv.w_axis.truncate();
        let asse = Vec3::from(o.asse).normalize_or_zero();
        let c = fwd.dot(asse).abs();
        // z_rif: the front of the light seen from the eye (its centre minus its half-extent along the view)
        let centro = Vec3::from(o.centro);
        let mezza = ((2.5 * o.sigma * c).powi(2) + (1.5 * o.w0).powi(2) * (1.0 - c * c)).sqrt();
        // depth range of pass 0: everything up to well beyond the light
        let z_max = (centro - eye).length() * 2.0 + 4.0 * mezza + 1.0;
        let ss = |e0: f32, e1: f32, x: f32| { let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0); t * t * (3.0 - 2.0 * t) };
        let lungo = ss(0.55, 0.95, c);
        let k = esposizione(&o, right, up, c) * (LATO_ESPOSIZIONE + (1.0 - LATO_ESPOSIZIONE) * lungo);
        let prof = o.profondita * (LATO_PROFONDITA + (1.0 - LATO_PROFONDITA) * lungo);
        let p = ParGpu {
            clip_from_world: clip_from_world.to_cols_array_2d(),
            right: right.extend(0.0).to_array(),
            up: up.extend(0.0).to_array(),
            eye: [eye.x, eye.y, eye.z, c],
            fwd: [fwd.x, fwd.y, fwd.z, 0.0],
            offset: [o.offset[0], o.offset[1], o.offset[2], prof * o.lambda],
            a: [o.size, o.gain, o.soft, k],
            b: [o.ns, o.coda_lato, o.guadagno, z_max],
        };
        let ub = device.create_buffer_with_data(&BufferInitDescriptor { label: Some("oit parametri"), contents: bytemuck::bytes_of(&p), usage: BufferUsages::UNIFORM });
        let bg_fronte = device.create_bind_group("oit fronte", &pipeline.layout_fronte, &BindGroupEntries::single(ub.as_entire_binding()));
        let bg_disegno = device.create_bind_group("oit palline", &pipeline.layout_disegno,
            &BindGroupEntries::with_indices(((0, ub.as_entire_binding()), (3, &fronte.default_view))));
        let bg_comp = device.create_bind_group("oit composizione", &pipeline.layout_comp,
            &BindGroupEntries::sequential((ub.as_entire_binding(), &accum.default_view, &reveal.default_view)));
        let campioni = msaa.map(|m| m.samples()).unwrap_or(1);
        let formato = if view.hdr { ViewTarget::TEXTURE_FORMAT_HDR } else { target.main_texture_format() };
        let comp = pipelines.specialize(&cache, &pipeline, ChiaveComp { formato, campioni });
        commands.entity(e).insert(VistaOit { fronte, bg_fronte, accum, reveal, bg_disegno, bg_comp, comp, buffer: buffer.clone(), n });
    }
}

#[derive(Default)]
pub struct NodoOit;

impl ViewNode for NodoOit {
    type ViewQuery = (&'static ViewTarget, &'static VistaOit);

    fn run<'w>(&self, _: &mut RenderGraphContext, rc: &mut RenderContext<'w>, (target, v): QueryItem<'w, Self::ViewQuery>, world: &'w World) -> Result<(), NodeRunError> {
        if v.n == 0 { return Ok(()); }
        let cache = world.resource::<PipelineCache>();
        let pip = world.resource::<PipelineOit>();
        let (Some(fronte), Some(disegno), Some(comp)) = (cache.get_render_pipeline(pip.fronte), cache.get_render_pipeline(pip.disegno), cache.get_render_pipeline(v.comp)) else { return Ok(()) };
        {
            let mut pass = rc.begin_tracked_render_pass(RenderPassDescriptor {
                label: Some("oit fronte"),
                color_attachments: &[],
                depth_stencil_attachment: Some(RenderPassDepthStencilAttachment {
                    view: &v.fronte.default_view,
                    depth_ops: Some(Operations { load: LoadOp::Clear(1.0), store: StoreOp::Store }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            pass.set_render_pipeline(fronte);
            pass.set_bind_group(0, &v.bg_fronte, &[]);
            pass.set_vertex_buffer(0, v.buffer.slice(..));
            pass.draw(0..6, 0..v.n);
        }
        {
            let mut pass = rc.begin_tracked_render_pass(RenderPassDescriptor {
                label: Some("oit palline"),
                color_attachments: &[
                    Some(RenderPassColorAttachment { view: &v.accum.default_view, resolve_target: None,
                        ops: Operations { load: LoadOp::Clear(wgpu::Color::TRANSPARENT), store: StoreOp::Store } }),
                    Some(RenderPassColorAttachment { view: &v.reveal.default_view, resolve_target: None,
                        ops: Operations { load: LoadOp::Clear(wgpu::Color::WHITE), store: StoreOp::Store } }),
                ],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            pass.set_render_pipeline(disegno);
            pass.set_bind_group(0, &v.bg_disegno, &[]);
            pass.set_vertex_buffer(0, v.buffer.slice(..));
            pass.draw(0..6, 0..v.n);
        }
        let mut pass = rc.begin_tracked_render_pass(RenderPassDescriptor {
            label: Some("oit composizione"),
            color_attachments: &[Some(target.get_color_attachment())],
            depth_stencil_attachment: None,
            timestamp_writes: None,
            occlusion_query_set: None,
        });
        pass.set_render_pipeline(comp);
        pass.set_bind_group(0, &v.bg_comp, &[]);
        pass.draw(0..3, 0..1);
        Ok(())
    }
}

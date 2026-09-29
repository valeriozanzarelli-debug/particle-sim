//! GPU-accelerated PPS simulation via a dedicated wgpu compute pipeline.
//!
//! Architecture
//! ─────────────────────────────────────────────────────────────────
//! `GpuCompute` is stored as a normal Bevy Resource (Send+Sync are
//! asserted with unsafe because wgpu types require it).  It owns its
//! own `wgpu::Device` + `wgpu::Queue` so it does not compete with
//! Bevy's render device — the physical GPU is shared, only the logical
//! device is separate.
//!
//! Per-frame flow (all N steps batched into ONE command buffer):
//!   • If world was reset / first frame: upload CPU positions/dirs → pos_a/dir_a
//!   • Encode N × (clear_grid → build_grid → simulate) with ping-pong buffers
//!   • Copy final output → staging, submit, poll(Wait), map → CPU
//!   • Write results back into `ParticleWorld`

use bevy::prelude::*;
use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::simulation::world::ParticleWorld;

// ─── public plugin ────────────────────────────────────────────────────────────

pub struct GpuPpsPlugin;

impl Plugin for GpuPpsPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, init_gpu_compute)
            .add_systems(
                Update,
                gpu_step_simulation
                    .run_if(resource_exists::<GpuCompute>)
                    // Must run before the mesh updater so it sees the new positions
                    .before(crate::rendering::point_cloud::update_point_cloud),
            );
    }
}

// ─── GPU parameter uniform ────────────────────────────────────────────────────
// Layout must match the WGSL `SimParams` struct exactly (12 × 4 bytes = 48 B).

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct SimParams {
    n:             u32,
    n_cells_total: u32,
    n_cells_axis:  u32,
    max_per_cell:  u32,
    alpha:         f32,
    beta:          f32,
    r_neighbor:    f32,
    speed:         f32,
    box_half:      f32,
    cell_size:     f32,
    _pad0:         u32,
    _pad1:         u32,
}

// ─── resource ─────────────────────────────────────────────────────────────────

/// Bevy resource wrapping all GPU-compute state.
/// We assert Send+Sync because wgpu types are not Send by default but we
/// guarantee single-threaded access through the ECS system scheduler.
#[derive(Resource)]
pub struct GpuCompute(GpuState);

// SAFETY: we only ever access GpuState from the main Bevy thread (non-parallel systems).
unsafe impl Send for GpuCompute {}
unsafe impl Sync for GpuCompute {}

pub struct GpuState {
    device: wgpu::Device,
    queue:  wgpu::Queue,

    // Particle data — ping-pong pair
    pos_a: wgpu::Buffer,   // vec4<f32> × N
    pos_b: wgpu::Buffer,
    dir_a: wgpu::Buffer,   // vec4<f32> × N
    dir_b: wgpu::Buffer,
    nbr:   wgpu::Buffer,   // u32 × N  (always written to, always the same buffer)

    // Spatial grid
    cell_counts: wgpu::Buffer,   // atomic<u32> × n_cells_total
    cell_data:   wgpu::Buffer,   // u32 × n_cells_total × max_per_cell

    // Uniform buffer
    params_buf: wgpu::Buffer,

    // CPU readback staging
    staging_pos: wgpu::Buffer,
    staging_dir: wgpu::Buffer,
    staging_nbr: wgpu::Buffer,

    // Compute pipelines (reused every frame)
    clear_pipeline:    wgpu::ComputePipeline,
    build_pipeline:    wgpu::ComputePipeline,
    simulate_pipeline: wgpu::ComputePipeline,

    bgl: wgpu::BindGroupLayout,

    // Allocation bookkeeping (to detect when buffers must be recreated)
    alloc_n:             usize,
    alloc_n_cells_total: u32,
    alloc_max_per_cell:  u32,

    // Ping-pong state: false → A is input / B is output this tick
    ping_pong: bool,
    // Last world.gpu_generation we uploaded — if it differs we re-upload
    last_generation: u64,

    pub gpu_name: String,
}

// ─── startup system ───────────────────────────────────────────────────────────

fn init_gpu_compute(mut commands: Commands, world: Res<ParticleWorld>) {
    match GpuState::new(&world) {
        Ok(state) => {
            info!("GPU compute ready — {}", state.gpu_name);
            commands.insert_resource(GpuCompute(state));
        }
        Err(e) => {
            error!("GPU init failed ({e}). Running on CPU.");
        }
    }
}

// ─── per-frame system ─────────────────────────────────────────────────────────

pub fn gpu_step_simulation(
    mut gc:    ResMut<GpuCompute>,
    mut world: ResMut<ParticleWorld>,
) {
    let steps = world.steps_per_frame;
    if steps == 0 { return; }
    gc.0.run_steps(&mut world, steps);
}

// ─── GpuState implementation ─────────────────────────────────────────────────

impl GpuState {
    fn new(world: &ParticleWorld) -> Result<Self, String> {
        // ── adapter + device ──────────────────────────────────────────────
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::VULKAN | wgpu::Backends::DX12,
            ..Default::default()
        });
        let adapter = pollster::block_on(instance.request_adapter(
            &wgpu::RequestAdapterOptions {
                power_preference:       wgpu::PowerPreference::HighPerformance,
                compatible_surface:     None,
                force_fallback_adapter: false,
            },
        ))
        .ok_or("No GPU adapter found")?;

        let gpu_name = adapter.get_info().name.clone();
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label:              Some("pps_compute"),
                required_features:  wgpu::Features::empty(),
                required_limits:    wgpu::Limits::default(),
                memory_hints:       wgpu::MemoryHints::Performance,
            },
            None,
        ))
        .map_err(|e| format!("Device error: {e}"))?;

        // ── shader ────────────────────────────────────────────────────────
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label:  Some("pps_compute"),
            source: wgpu::ShaderSource::Wgsl(
                include_str!("../../assets/shaders/pps_compute.wgsl").into(),
            ),
        });

        // ── bind group layout ─────────────────────────────────────────────
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label:   Some("pps_bgl"),
            entries: &[
                buf_entry(0, wgpu::BufferBindingType::Uniform,                   ),
                buf_entry(1, wgpu::BufferBindingType::Storage { read_only: true  }),
                buf_entry(2, wgpu::BufferBindingType::Storage { read_only: true  }),
                buf_entry(3, wgpu::BufferBindingType::Storage { read_only: false }),
                buf_entry(4, wgpu::BufferBindingType::Storage { read_only: false }),
                buf_entry(5, wgpu::BufferBindingType::Storage { read_only: false }),
                buf_entry(6, wgpu::BufferBindingType::Storage { read_only: false }),
                buf_entry(7, wgpu::BufferBindingType::Storage { read_only: false }),
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(
            &wgpu::PipelineLayoutDescriptor {
                label:                Some("pps_layout"),
                bind_group_layouts:   &[&bgl],
                push_constant_ranges: &[],
            },
        );

        let make = |entry: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label:               Some(entry),
                layout:              Some(&pipeline_layout),
                module:              &shader,
                entry_point:         Some(entry),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
                cache:               None,
            })
        };
        let clear_pipeline    = make("clear_grid");
        let build_pipeline    = make("build_grid");
        let simulate_pipeline = make("simulate");

        // ── initial buffers ───────────────────────────────────────────────
        let (_n_cells_axis, n_cells_total, max_per_cell) =
            grid_params(world.n, world.r_neighbor, world.box_half);

        let bufs = Buffers::alloc(&device, world, n_cells_total, max_per_cell);

        Ok(GpuState {
            device,
            queue,
            pos_a: bufs.pos_a, pos_b: bufs.pos_b,
            dir_a: bufs.dir_a, dir_b: bufs.dir_b,
            nbr:   bufs.nbr,
            cell_counts: bufs.cell_counts,
            cell_data:   bufs.cell_data,
            params_buf:  bufs.params_buf,
            staging_pos: bufs.staging_pos,
            staging_dir: bufs.staging_dir,
            staging_nbr: bufs.staging_nbr,
            clear_pipeline, build_pipeline, simulate_pipeline,
            bgl,
            alloc_n: world.n,
            alloc_n_cells_total: n_cells_total,
            alloc_max_per_cell:  max_per_cell,
            ping_pong: false,
            last_generation: u64::MAX, // force upload on first tick
            gpu_name,
        })
    }

    fn run_steps(&mut self, world: &mut ParticleWorld, steps: u32) {
        let n = world.n;
        let (n_cells_axis, n_cells_total, max_per_cell) =
            grid_params(n, world.r_neighbor, world.box_half);

        // ── reallocate buffers if N or grid size changed ──────────────────
        if n != self.alloc_n
            || n_cells_total != self.alloc_n_cells_total
            || max_per_cell  != self.alloc_max_per_cell
        {
            let bufs = Buffers::alloc(&self.device, world, n_cells_total, max_per_cell);
            self.pos_a = bufs.pos_a; self.pos_b = bufs.pos_b;
            self.dir_a = bufs.dir_a; self.dir_b = bufs.dir_b;
            self.nbr   = bufs.nbr;
            self.cell_counts = bufs.cell_counts;
            self.cell_data   = bufs.cell_data;
            self.params_buf  = bufs.params_buf;
            self.staging_pos = bufs.staging_pos;
            self.staging_dir = bufs.staging_dir;
            self.staging_nbr = bufs.staging_nbr;
            self.alloc_n             = n;
            self.alloc_n_cells_total = n_cells_total;
            self.alloc_max_per_cell  = max_per_cell;
            self.ping_pong       = false;
            self.last_generation = u64::MAX; // force re-upload
        }

        // ── upload CPU state when world was reset ────────────────────────
        if self.last_generation != world.gpu_generation {
            let pos_data: Vec<[f32; 4]> = world.positions.iter()
                .map(|p| [p.x, p.y, p.z, 0.0]).collect();
            let dir_data: Vec<[f32; 4]> = world.directions.iter()
                .map(|d| [d.x, d.y, d.z, 0.0]).collect();
            self.queue.write_buffer(&self.pos_a, 0, bytemuck::cast_slice(&pos_data));
            self.queue.write_buffer(&self.dir_a, 0, bytemuck::cast_slice(&dir_data));
            self.ping_pong       = false;
            self.last_generation = world.gpu_generation;
        }

        // ── write simulation parameters ───────────────────────────────────
        let cell_size = world.r_neighbor.max(2.0_f32);
        let params = SimParams {
            n:             n as u32,
            n_cells_total,
            n_cells_axis,
            max_per_cell,
            alpha:      world.alpha,
            beta:       world.beta,
            r_neighbor: world.r_neighbor,
            speed:      world.speed,
            box_half:   world.box_half,
            cell_size,
            _pad0: 0,
            _pad1: 0,
        };
        self.queue.write_buffer(&self.params_buf, 0, bytemuck::bytes_of(&params));

        // ── pre-build both ping-pong bind groups ──────────────────────────
        // bg_ab: pos_a/dir_a → pos_b/dir_b
        let bg_ab = self.make_bg(&self.pos_a, &self.dir_a, &self.pos_b, &self.dir_b);
        // bg_ba: pos_b/dir_b → pos_a/dir_a
        let bg_ba = self.make_bg(&self.pos_b, &self.dir_b, &self.pos_a, &self.dir_a);

        // ── encode all steps into a single command buffer ─────────────────
        let mut encoder = self.device.create_command_encoder(
            &wgpu::CommandEncoderDescriptor { label: Some("pps") },
        );

        let n_cells_wg = div_ceil(n_cells_total, 64);
        let n_wg       = div_ceil(n as u32, 64);

        for _ in 0..steps {
            let bg = if !self.ping_pong { &bg_ab } else { &bg_ba };

            encode_pass(&mut encoder, &self.clear_pipeline,    bg, n_cells_wg, "clr");
            encode_pass(&mut encoder, &self.build_pipeline,    bg, n_wg,       "bld");
            encode_pass(&mut encoder, &self.simulate_pipeline, bg, n_wg,       "sim");

            self.ping_pong = !self.ping_pong;
        }

        // After `steps` ticks, ping_pong tells us which buffer was last written:
        //   ping_pong is now true  → last write was to B (we started ping_pong=false, wrote B, flipped)
        //   ping_pong is now false → last write was to A
        // More precisely: we wrote to B when ping_pong was false before the tick.
        // After odd number of ticks starting from false: ping_pong = true → last output = B
        // After even number: ping_pong = false → last output = A
        let (last_pos, last_dir) = if self.ping_pong {
            (&self.pos_b, &self.dir_b)
        } else {
            (&self.pos_a, &self.dir_a)
        };

        let byte_v4 = (n * 4 * 4) as u64;
        let byte_u1 = (n * 4)     as u64;
        encoder.copy_buffer_to_buffer(last_pos,    0, &self.staging_pos, 0, byte_v4);
        encoder.copy_buffer_to_buffer(last_dir,    0, &self.staging_dir, 0, byte_v4);
        encoder.copy_buffer_to_buffer(&self.nbr,   0, &self.staging_nbr, 0, byte_u1);

        self.queue.submit(std::iter::once(encoder.finish()));

        // ── synchronous readback ──────────────────────────────────────────
        map_and_wait(
            &self.device,
            &self.staging_pos,
            &self.staging_dir,
            &self.staging_nbr,
        );

        {
            let v = self.staging_pos.slice(..).get_mapped_range();
            let d: &[[f32; 4]] = bytemuck::cast_slice(&v);
            for (i, &[x, y, z, _]) in d.iter().take(n).enumerate() {
                world.positions[i] = Vec3::new(x, y, z);
            }
        }
        self.staging_pos.unmap();

        {
            let v = self.staging_dir.slice(..).get_mapped_range();
            let d: &[[f32; 4]] = bytemuck::cast_slice(&v);
            for (i, &[x, y, z, _]) in d.iter().take(n).enumerate() {
                world.directions[i] = Vec3::new(x, y, z);
            }
        }
        self.staging_dir.unmap();

        {
            let v = self.staging_nbr.slice(..).get_mapped_range();
            let d: &[u32] = bytemuck::cast_slice(&v);
            let mut sum = 0u64;
            let mut max_c = 0u8;
            for (i, &c) in d.iter().take(n).enumerate() {
                let c8 = c.min(255) as u8;
                world.neighbor_counts[i] = c8;
                sum   += c8 as u64;
                if c8 > max_c { max_c = c8; }
            }
            world.mean_neighbors = sum as f32 / n as f32;
            world.max_cluster    = max_c;
        }
        self.staging_nbr.unmap();

        world.ticks += steps as u64;
    }

    fn make_bg(
        &self,
        pos_in:  &wgpu::Buffer,
        dir_in:  &wgpu::Buffer,
        pos_out: &wgpu::Buffer,
        dir_out: &wgpu::Buffer,
    ) -> wgpu::BindGroup {
        self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label:   Some("pps_bg"),
            layout:  &self.bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: self.params_buf .as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: pos_in          .as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: dir_in          .as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: pos_out         .as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: dir_out         .as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: self.nbr        .as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: self.cell_counts.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: self.cell_data  .as_entire_binding() },
            ],
        })
    }
}

// ─── expose dirty flag for world resets ──────────────────────────────────────

impl GpuCompute {
    pub fn gpu_name(&self) -> &str {
        &self.0.gpu_name
    }
}

// ─── helpers ─────────────────────────────────────────────────────────────────

fn grid_params(n: usize, r: f32, box_half: f32) -> (u32, u32, u32) {
    let box_size     = box_half * 2.0;
    let cell_size    = r.max(2.0_f32);
    let n_cells_axis = ((box_size / cell_size).ceil() as u32).clamp(1, 40);
    let n_cells_total = n_cells_axis.pow(3);
    let avg           = (n as f32 / n_cells_total as f32).ceil() as u32;
    let max_per_cell  = (avg * 8).max(64);
    (n_cells_axis, n_cells_total, max_per_cell)
}

fn encode_pass(
    encoder:  &mut wgpu::CommandEncoder,
    pipeline: &wgpu::ComputePipeline,
    bg:       &wgpu::BindGroup,
    wg:       u32,
    label:    &str,
) {
    let mut pass = encoder.begin_compute_pass(
        &wgpu::ComputePassDescriptor { label: Some(label), timestamp_writes: None },
    );
    pass.set_pipeline(pipeline);
    pass.set_bind_group(0, bg, &[]);
    pass.dispatch_workgroups(wg, 1, 1);
}

fn map_and_wait(
    device:      &wgpu::Device,
    staging_pos: &wgpu::Buffer,
    staging_dir: &wgpu::Buffer,
    staging_nbr: &wgpu::Buffer,
) {
    // Use unit-typed channels — we only need the signal that mapping completed.
    let (p_tx, p_rx) = std::sync::mpsc::channel::<()>();
    let (d_tx, d_rx) = std::sync::mpsc::channel::<()>();
    let (n_tx, n_rx) = std::sync::mpsc::channel::<()>();

    staging_pos.slice(..).map_async(wgpu::MapMode::Read, move |_| { let _ = p_tx.send(()); });
    staging_dir.slice(..).map_async(wgpu::MapMode::Read, move |_| { let _ = d_tx.send(()); });
    staging_nbr.slice(..).map_async(wgpu::MapMode::Read, move |_| { let _ = n_tx.send(()); });

    // Wait for the GPU to finish executing all submitted commands,
    // which will also trigger the map callbacks above.
    device.poll(wgpu::Maintain::Wait);

    let _ = p_rx.try_recv();
    let _ = d_rx.try_recv();
    let _ = n_rx.try_recv();
}

fn buf_entry(binding: u32, ty: wgpu::BufferBindingType) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty,
            has_dynamic_offset: false,
            min_binding_size:   None,
        },
        count: None,
    }
}

#[inline]
fn div_ceil(n: u32, d: u32) -> u32 { (n + d - 1) / d }

// ─── buffer allocation ────────────────────────────────────────────────────────

struct Buffers {
    pos_a: wgpu::Buffer, pos_b: wgpu::Buffer,
    dir_a: wgpu::Buffer, dir_b: wgpu::Buffer,
    nbr:   wgpu::Buffer,
    cell_counts: wgpu::Buffer,
    cell_data:   wgpu::Buffer,
    params_buf:  wgpu::Buffer,
    staging_pos: wgpu::Buffer,
    staging_dir: wgpu::Buffer,
    staging_nbr: wgpu::Buffer,
}

impl Buffers {
    fn alloc(
        device:        &wgpu::Device,
        world:         &ParticleWorld,
        n_cells_total: u32,
        max_per_cell:  u32,
    ) -> Self {
        let n = world.n;
        let pos_data: Vec<[f32; 4]> = world.positions.iter()
            .map(|p| [p.x, p.y, p.z, 0.0]).collect();
        let dir_data: Vec<[f32; 4]> = world.directions.iter()
            .map(|d| [d.x, d.y, d.z, 0.0]).collect();

        // All particle buffers need STORAGE + COPY_SRC + COPY_DST:
        //   – STORAGE  : bound as storage buffer in compute shaders
        //   – COPY_SRC : copied to staging buffer for readback
        //   – COPY_DST : written by queue.write_buffer on initial upload
        let rw = wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST;

        let pos_a = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("pos_a"), contents: bytemuck::cast_slice(&pos_data), usage: rw,
        });
        let pos_b = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pos_b"), size: (n * 16) as u64, usage: rw, mapped_at_creation: false,
        });
        let dir_a = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("dir_a"), contents: bytemuck::cast_slice(&dir_data), usage: rw,
        });
        let dir_b = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("dir_b"), size: (n * 16) as u64, usage: rw, mapped_at_creation: false,
        });
        let nbr = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("nbr"), size: (n * 4) as u64, usage: rw, mapped_at_creation: false,
        });

        let grid_rw = wgpu::BufferUsages::STORAGE;
        let cell_counts = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("cell_counts"),
            size:  (n_cells_total * 4) as u64,
            usage: grid_rw,
            mapped_at_creation: false,
        });
        let cell_data_size = ((n_cells_total as u64) * (max_per_cell as u64) * 4).max(4);
        let cell_data = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("cell_data"),
            size:  cell_data_size,
            usage: grid_rw,
            mapped_at_creation: false,
        });

        let params_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("params"),
            size:  std::mem::size_of::<SimParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let staging = wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST;
        let staging_pos = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("stg_pos"), size: (n * 16) as u64, usage: staging, mapped_at_creation: false,
        });
        let staging_dir = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("stg_dir"), size: (n * 16) as u64, usage: staging, mapped_at_creation: false,
        });
        let staging_nbr = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("stg_nbr"), size: (n * 4) as u64,  usage: staging, mapped_at_creation: false,
        });

        Buffers {
            pos_a, pos_b, dir_a, dir_b, nbr,
            cell_counts, cell_data, params_buf,
            staging_pos, staging_dir, staging_nbr,
        }
    }
}

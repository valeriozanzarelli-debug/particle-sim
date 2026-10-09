//! GPU port of the photon grains (wgpu compute) — NOT yet wired into `Sim`.
//!
//! Mirrors `advance()` / `huygens()` / `deposit()` / the photon-level Born
//! decisions of `physics.rs` on a dedicated wgpu device, so a photon can be
//! made of millions of grains. See `docs/GPU_PORT.md` and
//! `assets/shaders/photon_grains.wgsl`.
//!
//! Differences from the old `simulation/gpu_pps.rs`, on purpose:
//!   • no per-frame readback of the grains: they live in `grain_buffer()`
//!     (the renderer is meant to read it directly); per tick nothing is read
//!     back — only `read_status()` (32 B/photon) and `read_screen()`
//!     (48 B/bin/photon) when the host asks;
//!   • no capped per-cell storage: the screen is a histogram of i64
//!     fixed-point atomics, nothing is ever dropped;
//!   • all passes of a tick in one command buffer, no CPU round trip for the
//!     Born-rule decisions (pre-drawn uniforms live in the photon record).
//!
//! Scope: FreeBeam, DoubleSlit, Polarizers are complete. For Bell, the
//! polarizing beam splitter is a pass-through here: the entangled collapse
//! (`entangle::measure`) stays on the host and must be applied with
//! `set_photon_jones` after reading `status` (not yet implemented).

// The module is compiled into `fotone` but not used by the live app yet.
#![allow(dead_code)]

use bytemuck::{Pod, Zeroable};

use crate::physics::{Experiment, Grain, Lab};

/// Threads per workgroup (must match `WG` in the shader). Every photon's
/// grain range starts at a multiple of this.
pub const WG: usize = 64;
pub const MAX_ELEMS: usize = 4;
const STATUS_WORDS: usize = 8;
const SCREEN_WORDS: usize = 12;

const ST_ALIVE: u32 = 0x100;
const ST_DIFFR: u32 = 0x200;

const K_WALL: u32 = 1;
const K_POLARIZER: u32 = 2;
const K_BELL: u32 = 3;
const K_SCREEN: u32 = 4;

pub const FATE_IN_FLIGHT: u32 = 0;
pub const FATE_WALL: u32 = 1;
pub const FATE_POLARIZER: u32 = 2;

// ─── GPU structs (layouts documented in the shader header) ───────────────────

/// WGSL `Params` — uniform, 96 B.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug)]
struct Params {
    dist:        f32,
    k:           f32,
    lambda:      f32,
    waist:       f32,
    theta_max:   f32,
    screen_half: f32,
    n_bins:      u32,
    n_grains:    u32,
    slit_c0:     f32,
    slit_c1:     f32,
    slit_half:   f32,
    filter_mask: u32,
    filter0:     [f32; 4],
    filter1:     [f32; 4],
    amp_scale:   f32,
    self_scale:  f32,
    prob_scale:  f32,
    n_slots:     u32,
}

/// WGSL `Grain` — storage, 64 B, 16-byte aligned. The renderer can bind the
/// grain buffer as an instance vertex buffer: stride 64, `pos` at offset 0
/// (Float32x3), `phase0` at 28, `origin` at 32, `state` at 44, `jones` at 48.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug, PartialEq)]
pub struct GrainGpu {
    pub pos:    [f32; 3],
    pub amp:    f32,
    pub dir:    [f32; 3],
    pub phase0: f32,
    pub origin: [f32; 3],
    /// bits 0-7 next element, 8 alive, 9 diffracted, 16-31 photon slot
    pub state:  u32,
    /// (re e1, im e1, re e2, im e2)
    pub jones:  [f32; 4],
}

/// WGSL `Photon` — storage, 112 B.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug)]
struct PhotonGpu {
    elem_z:      [f32; 4],
    elem_kind:   [u32; 4],
    elem_arg:    [f32; 4],
    rand:        [f32; 4],
    op_jones:    [f32; 4],
    n_elems:     u32,
    key_lo:      u32,
    key_hi:      u32,
    first_grain: u32,
    n_grains:    u32,
    front:       u32,
    op:          u32,
    fate:        u32,
}

const _: () = assert!(std::mem::size_of::<Params>() == 96);
const _: () = assert!(std::mem::size_of::<GrainGpu>() == 64);
const _: () = assert!(std::mem::size_of::<PhotonGpu>() == 112);

impl GrainGpu {
    pub fn from_grain(g: &Grain, slot: u32) -> Self {
        let mut state = (g.next as u32 & 0xFF) | (slot << 16);
        if g.alive { state |= ST_ALIVE; }
        if g.diffracted { state |= ST_DIFFR; }
        Self {
            pos: g.pos.to_array(), amp: g.amp,
            dir: g.dir.to_array(), phase0: g.phase0,
            origin: g.origin.to_array(), state,
            jones: [g.jones[0].re, g.jones[0].im, g.jones[1].re, g.jones[1].im],
        }
    }
    pub fn alive(&self) -> bool { self.state & ST_ALIVE != 0 }
    pub fn diffracted(&self) -> bool { self.state & ST_DIFFR != 0 }
    pub fn next(&self) -> u8 { (self.state & 0xFF) as u8 }
}

// ─── Optical elements of one photon ──────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GpuElement {
    Wall,
    /// `lab_index` into `lab.polarizers`, transmission axis `angle`.
    Polarizer { lab_index: usize, angle: f32 },
    BellSplitter,
    Screen,
}

/// The element list of `Sim::emit` (physics.rs), sorted along travel.
pub fn photon_elements(lab: &Lab, side: usize) -> Vec<(f32, GpuElement)> {
    let sign = if side == 0 { 1.0 } else { -1.0 };
    let mut elems = Vec::new();
    match lab.experiment {
        Experiment::FreeBeam | Experiment::Airy | Experiment::Arago => {}
        Experiment::DoubleSlit => elems.push((lab.wall_z, GpuElement::Wall)),
        Experiment::Polarizers => {
            for (i, p) in lab.polarizers.iter().enumerate() {
                if p.enabled { elems.push((p.z, GpuElement::Polarizer { lab_index: i, angle: p.angle })); }
            }
        }
        Experiment::Bell => {
            let d = if side == 0 { lab.bell.dist_a } else { lab.bell.dist_b };
            elems.push((sign * d, GpuElement::BellSplitter));
        }
    }
    let screen_z = match lab.experiment {
        Experiment::Bell => sign * (if side == 0 { lab.bell.dist_a } else { lab.bell.dist_b } + 15.0),
        _ => lab.screen_z,
    };
    elems.push((screen_z, GpuElement::Screen));
    elems.sort_by(|a, b| (sign * a.0).partial_cmp(&(sign * b.0)).unwrap());
    elems
}

/// physics.rs `hash3(seed, id, gi)` = mix(seed ^ id·C1 ^ gi·C2): the photon
/// key is the part that does not depend on the grain.
fn photon_key(seed: u64, photon_id: u64) -> u64 {
    seed ^ photon_id.wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

// ─── Readbacks ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, Default)]
pub struct PhotonStatus {
    /// Bit e set: the photon's front reached element e during the last tick.
    pub reach_mask: u32,
    /// Grains still alive after the last tick (0 → the photon is finished).
    pub alive: u32,
    /// Born-rule sums of the last probe (front element): Σ|J|², Σ passing.
    pub den: f64,
    pub num: f64,
    /// FATE_* (FATE_POLARIZER | lab_index << 8).
    pub fate: u32,
    /// Next photon-level element.
    pub front: u32,
}

/// One photon's screen, as in physics.rs `ScreenAccum` (without the hits).
#[derive(Clone, Debug)]
pub struct ScreenBins {
    /// Σ a per bin: (re e1, im e1, re e2, im e2).
    pub sum:   Vec<[f64; 4]>,
    /// Σ |a|² per bin.
    pub self_: Vec<f64>,
    pub count: Vec<u32>,
}

impl ScreenBins {
    /// The photon's intensity, exactly as `Sim::finalize` builds it.
    pub fn intensity(&self, diffracted: bool) -> Vec<f64> {
        (0..self.count.len()).map(|b| {
            if diffracted {
                let s = self.sum[b];
                (s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + s[3] * s[3] - self.self_[b]).max(0.0)
            } else {
                self.count[b] as f64
            }
        }).collect()
    }
}

fn i64_at(w: &[u32], i: usize) -> i64 {
    ((w[i + 1] as i64) << 32) | w[i] as i64
}

// ─── The engine ──────────────────────────────────────────────────────────────

pub struct GrainEngine {
    device: wgpu::Device,
    queue:  wgpu::Queue,
    pub adapter_name: String,

    pipe_probe:   wgpu::ComputePipeline,
    pipe_decide:  wgpu::ComputePipeline,
    pipe_advance: wgpu::ComputePipeline,
    bind_group:   wgpu::BindGroup,

    params_buf:  wgpu::Buffer,
    grains_buf:  wgpu::Buffer,
    photons_buf: wgpu::Buffer,
    status_buf:  wgpu::Buffer,
    screen_buf:  wgpu::Buffer,
    staging_status: wgpu::Buffer,
    staging_screen: wgpu::Buffer,

    params:   Params,
    /// Grain capacity (multiple of WG).
    capacity: usize,
    n_slots:  usize,
    bins:     usize,
}

impl GrainEngine {
    /// Allocate room for `capacity` grains shared by `n_slots` photons, each
    /// with a screen of `bins` bins. Uses default limits unless the grain
    /// buffer needs more, and then only what the adapter reports.
    pub fn new(capacity: usize, n_slots: usize, bins: usize) -> Result<Self, String> {
        let capacity = capacity.div_ceil(WG).max(1) * WG;
        let n_slots = n_slots.max(1);
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::VULKAN | wgpu::Backends::DX12 | wgpu::Backends::METAL,
            ..Default::default()
        });
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .ok_or("no GPU adapter found")?;
        let adapter_name = format!("{} ({:?})", adapter.get_info().name, adapter.get_info().backend);

        // Limits: ask for more than the defaults only when needed, and never
        // more than the adapter has.
        let grain_bytes = (capacity * std::mem::size_of::<GrainGpu>()) as u64;
        let screen_bytes = (n_slots * bins * SCREEN_WORDS * 4) as u64;
        let biggest = grain_bytes.max(screen_bytes);
        let have = adapter.limits();
        let mut limits = wgpu::Limits::default();
        if biggest > limits.max_storage_buffer_binding_size as u64 {
            if biggest > have.max_storage_buffer_binding_size as u64 {
                return Err(format!(
                    "{} grains need a {} MiB storage binding; the adapter allows {} MiB",
                    capacity, biggest >> 20, have.max_storage_buffer_binding_size >> 20));
            }
            limits.max_storage_buffer_binding_size = biggest as u32;
        }
        if biggest > limits.max_buffer_size {
            if biggest > have.max_buffer_size {
                return Err(format!("buffer of {} MiB exceeds the adapter's max_buffer_size", biggest >> 20));
            }
            limits.max_buffer_size = biggest;
        }
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("photon_grains"),
            required_features: wgpu::Features::empty(),
            required_limits: limits,
            memory_hints: wgpu::MemoryHints::Performance,
        }, None))
        .map_err(|e| format!("device error: {e}"))?;

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("photon_grains"),
            source: wgpu::ShaderSource::Wgsl(include_str!("../../assets/shaders/photon_grains.wgsl").into()),
        });

        let storage = |binding| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only: false },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("photon_grains_bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                storage(1), storage(2), storage(3), storage(4),
            ],
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("photon_grains_layout"),
            bind_group_layouts: &[&bgl],
            push_constant_ranges: &[],
        });
        let make = |entry: &str| device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(entry),
            layout: Some(&layout),
            module: &shader,
            entry_point: Some(entry),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });
        let (pipe_probe, pipe_decide, pipe_advance) = (make("probe"), make("decide"), make("advance"));

        let buf = |label, size: u64, usage| device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label), size: size.max(16), usage, mapped_at_creation: false,
        });
        use wgpu::BufferUsages as U;
        let params_buf = buf("params", std::mem::size_of::<Params>() as u64, U::UNIFORM | U::COPY_DST);
        // VERTEX: the renderer may bind it as an instance buffer (same device).
        let grains_buf = buf("grains", grain_bytes, U::STORAGE | U::COPY_DST | U::COPY_SRC | U::VERTEX);
        let photons_buf = buf("photons", (n_slots * std::mem::size_of::<PhotonGpu>()) as u64, U::STORAGE | U::COPY_DST);
        let status_bytes = (n_slots * STATUS_WORDS * 4) as u64;
        let status_buf = buf("status", status_bytes, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let screen_buf = buf("screen", screen_bytes, U::STORAGE | U::COPY_DST | U::COPY_SRC);
        let staging_status = buf("staging_status", status_bytes, U::MAP_READ | U::COPY_DST);
        let staging_screen = buf("staging_screen", (bins * SCREEN_WORDS * 4) as u64, U::MAP_READ | U::COPY_DST);

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("photon_grains_bg"),
            layout: &bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: grains_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: photons_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: status_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: screen_buf.as_entire_binding() },
            ],
        });

        // Everything starts dead / empty.
        let mut enc = device.create_command_encoder(&Default::default());
        enc.clear_buffer(&grains_buf, 0, None);
        enc.clear_buffer(&photons_buf, 0, None);
        enc.clear_buffer(&screen_buf, 0, None);
        queue.submit(Some(enc.finish()));

        let params = Params {
            dist: 0.0, k: 1.0, lambda: 1.0, waist: 1.0, theta_max: 0.0, screen_half: 1.0,
            n_bins: bins as u32, n_grains: capacity as u32,
            slit_c0: 0.0, slit_c1: 0.0, slit_half: 0.0, filter_mask: 0,
            filter0: [0.0; 4], filter1: [0.0; 4],
            amp_scale: 1.0, self_scale: 1.0, prob_scale: (1u32 << 24) as f32,
            n_slots: n_slots as u32,
        };
        Ok(Self {
            device, queue, adapter_name,
            pipe_probe, pipe_decide, pipe_advance, bind_group,
            params_buf, grains_buf, photons_buf, status_buf, screen_buf,
            staging_status, staging_screen,
            params, capacity, n_slots, bins,
        })
    }

    pub fn capacity(&self) -> usize { self.capacity }

    /// The grain buffer (`GrainGpu` × capacity) for rendering on this device.
    pub fn grain_buffer(&self) -> &wgpu::Buffer { &self.grains_buf }

    /// Load the laboratory constants (call again when the lab changes).
    pub fn set_lab(&mut self, lab: &Lab) {
        assert_eq!(lab.bins, self.bins, "screen bins fixed at allocation");
        let (c, half) = lab.slits();
        let filters = lab.slit_filters();
        let pack = |f: Option<crate::complex::Jones>| f.map_or([0.0; 4], |j| [j[0].re, j[0].im, j[1].re, j[1].im]);

        // Fixed-point scales: the largest single deposit must fit in 2^30.
        // |a| ≤ w ≤ amp_max · L / (cos²θmax · √(λL)), amp_max = 1/U_inc at the slit edge.
        let edge = c[0].abs().max(c[1].abs()) + half;
        let inc_min = (-(edge * edge) / (lab.waist * lab.waist)).exp().max(1e-30);
        let l = (lab.screen_z - lab.wall_z).abs().max(lab.lambda);
        let cos2 = lab.theta_max.cos().powi(2).max(1e-6);
        let w_max = (l / (cos2 * (lab.lambda * l).sqrt()) / inc_min).max(1e-6);
        let pow2_floor = |x: f32| 2f32.powi(x.log2().floor() as i32);
        let limit = (1u32 << 30) as f32;

        self.params = Params {
            k: lab.k(), lambda: lab.lambda, waist: lab.waist, theta_max: lab.theta_max,
            screen_half: lab.screen_half,
            slit_c0: c[0], slit_c1: c[1], slit_half: half,
            filter_mask: filters[0].is_some() as u32 | (filters[1].is_some() as u32) << 1,
            filter0: pack(filters[0]), filter1: pack(filters[1]),
            amp_scale: pow2_floor(limit / w_max),
            self_scale: pow2_floor(limit / (w_max * w_max)),
            ..self.params
        };
    }

    /// Put one photon into `slot`: its grains go to `first_grain..` (a
    /// multiple of WG), its screen is cleared. `rand[e]` is the uniform the
    /// Born rule of element e will use (draw it from the host RNG).
    pub fn launch(&mut self, slot: usize, first_grain: usize, grains: &[Grain], photon_id: u64,
                  seed: u64, elems: &[(f32, GpuElement)], rand: [f32; MAX_ELEMS]) {
        assert!(slot < self.n_slots);
        assert!(first_grain % WG == 0, "photon ranges must start on a workgroup boundary");
        assert!(elems.len() <= MAX_ELEMS);
        let padded = grains.len().div_ceil(WG) * WG;
        assert!(first_grain + padded <= self.capacity, "grain buffer too small");

        let mut data: Vec<GrainGpu> = grains.iter().map(|g| GrainGpu::from_grain(g, slot as u32)).collect();
        data.resize(padded, GrainGpu { state: (slot as u32) << 16, ..GrainGpu::zeroed() });
        self.queue.write_buffer(&self.grains_buf, (first_grain * 64) as u64, bytemuck::cast_slice(&data));

        let key = photon_key(seed, photon_id);
        let mut ph = PhotonGpu {
            n_elems: elems.len() as u32,
            key_lo: key as u32, key_hi: (key >> 32) as u32,
            first_grain: first_grain as u32, n_grains: grains.len() as u32,
            rand, ..PhotonGpu::zeroed()
        };
        for (e, &(z, el)) in elems.iter().enumerate() {
            ph.elem_z[e] = z;
            (ph.elem_kind[e], ph.elem_arg[e]) = match el {
                GpuElement::Wall => (K_WALL, 0.0),
                GpuElement::Polarizer { lab_index, angle } => (K_POLARIZER | (lab_index as u32) << 8, angle),
                GpuElement::BellSplitter => (K_BELL, 0.0),
                GpuElement::Screen => (K_SCREEN, 0.0),
            };
        }
        self.queue.write_buffer(&self.photons_buf, (slot * std::mem::size_of::<PhotonGpu>()) as u64, bytemuck::bytes_of(&ph));

        let stride = (self.bins * SCREEN_WORDS * 4) as u64;
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&self.screen_buf, slot as u64 * stride, Some(stride));
        self.queue.submit(Some(enc.finish()));
    }

    /// `n` ticks of light path `dist` each, in one command buffer
    /// (clear status → probe → decide → advance, per tick). No readback.
    pub fn tick(&mut self, dist: f32, n: u32) {
        self.params.dist = dist;
        self.queue.write_buffer(&self.params_buf, 0, bytemuck::bytes_of(&self.params));
        let groups = (self.capacity / WG) as u32;
        let (gx, gy) = if groups <= 65_535 { (groups, 1) } else { (65_535, groups.div_ceil(65_535)) };
        let slot_groups = (self.n_slots as u32).div_ceil(WG as u32);

        let mut enc = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("photon_tick") });
        for _ in 0..n {
            enc.clear_buffer(&self.status_buf, 0, None);
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some("photon_tick"), timestamp_writes: None });
            pass.set_bind_group(0, &self.bind_group, &[]);
            pass.set_pipeline(&self.pipe_probe);
            pass.dispatch_workgroups(gx, gy, 1);
            pass.set_pipeline(&self.pipe_decide);
            pass.dispatch_workgroups(slot_groups, 1, 1);
            pass.set_pipeline(&self.pipe_advance);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        self.queue.submit(Some(enc.finish()));
    }

    /// Block until the GPU has finished all submitted work (benchmarks).
    pub fn wait(&self) {
        self.device.poll(wgpu::Maintain::Wait);
    }

    /// Small readback: photon-level state of every slot after the last tick.
    pub fn read_status(&self) -> Vec<PhotonStatus> {
        let w: Vec<u32> = self.read(&self.status_buf, 0, &self.staging_status, (self.n_slots * STATUS_WORDS * 4) as u64);
        w.chunks_exact(STATUS_WORDS).map(|s| PhotonStatus {
            reach_mask: s[0], alive: s[1],
            den: i64_at(s, 2) as f64 / self.params.prob_scale as f64,
            num: i64_at(s, 4) as f64 / self.params.prob_scale as f64,
            fate: s[6], front: s[7],
        }).collect()
    }

    /// Small readback: the screen bins of one photon (never the grains).
    pub fn read_screen(&self, slot: usize) -> ScreenBins {
        let stride = (self.bins * SCREEN_WORDS * 4) as u64;
        let w: Vec<u32> = self.read(&self.screen_buf, slot as u64 * stride, &self.staging_screen, stride);
        let (a, s) = (self.params.amp_scale as f64, self.params.self_scale as f64);
        let mut out = ScreenBins { sum: Vec::new(), self_: Vec::new(), count: Vec::new() };
        for b in w.chunks_exact(SCREEN_WORDS) {
            out.sum.push([i64_at(b, 0) as f64 / a, i64_at(b, 2) as f64 / a, i64_at(b, 4) as f64 / a, i64_at(b, 6) as f64 / a]);
            out.self_.push(i64_at(b, 8) as f64 / s);
            out.count.push(b[10]);
        }
        out
    }

    /// Debug/test only: full copy of a grain range to the CPU.
    #[cfg(test)]
    fn read_grains(&self, first: usize, n: usize) -> Vec<GrainGpu> {
        let bytes = (n * 64) as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("staging_grains"), size: bytes,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        self.read(&self.grains_buf, (first * 64) as u64, &staging, bytes)
    }

    fn read<T: Pod>(&self, src: &wgpu::Buffer, offset: u64, staging: &wgpu::Buffer, bytes: u64) -> Vec<T> {
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(src, offset, staging, 0, bytes);
        self.queue.submit(Some(enc.finish()));
        let slice = staging.slice(..bytes);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx.send(r); });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv().expect("map callback").expect("map failed");
        let out = bytemuck::cast_slice(&slice.get_mapped_range()).to_vec();
        staging.unmap();
        out
    }
}

// ─── CPU vs GPU ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::complex::{inner, jones_norm2};
    use crate::physics::{correlation, Sim};

    /// The CPU's wall decision probability (physics.rs `wall_decision`).
    fn cpu_p_pass(lab: &Lab, grains: &[Grain]) -> f64 {
        let (c, half) = lab.slits();
        let filters = lab.slit_filters();
        let (mut through, mut total) = (0.0f64, 0.0f64);
        for g in grains.iter().filter(|g| g.alive) {
            let w = jones_norm2(g.jones) as f64;
            total += w;
            let t = (lab.wall_z - g.pos.z) / g.dir.z;
            let x = g.pos.x + g.dir.x * t;
            for s in 0..2 {
                if (x - c[s]).abs() <= half {
                    through += match filters[s] { Some(f) => inner(f, g.jones).norm2() as f64, None => w };
                }
            }
        }
        through / total
    }

    fn reaches(grains: &[Grain], z: f32, dist: f32) -> bool {
        grains.iter().any(|g| g.alive && { let t = (z - g.pos.z) / g.dir.z; (0.0..=dist).contains(&t) })
    }

    /// The same photon (same grains, same per-grain hash) on CPU `Sim` and on
    /// the GPU: grain positions mid-flight and the screen intensity.
    #[test]
    fn gpu_matches_cpu_double_slit() {
        let dist = 7.3f32;          // not a divisor of anything: exercises mid-step crossings
        let mid_step = 15;          // clock 109.5: grains between wall and screen
        let mut lab = Lab::preset(Experiment::DoubleSlit);

        // ── CPU: first seed whose photon passes the wall ──
        let mut found = None;
        for seed in 1..=40 {
            lab.seed = seed;
            let mut sim = Sim::new(lab.clone());
            sim.launch();
            let init = sim.photons[0].grains.clone();
            let (mut steps, mut mid, mut wall) = (0, None, None);
            let t0 = std::time::Instant::now();
            while sim.in_flight() && steps < 500 {
                if wall.is_none() && reaches(&sim.photons[0].grains, lab.wall_z, dist) {
                    wall = Some((steps, cpu_p_pass(&lab, &sim.photons[0].grains)));
                }
                sim.step(dist);
                steps += 1;
                if steps == mid_step { mid = sim.photons.first().map(|p| p.grains.clone()); }
            }
            if sim.stats.photons_in_mean == 1 {
                let cpu_ms = t0.elapsed().as_secs_f64() * 1e3;
                found = Some((init, steps, mid.unwrap(), wall.unwrap(), sim.stats.mean_intensity.clone(), cpu_ms));
                break;
            }
        }
        let (init, steps, cpu_mid, (wall_step, p_cpu), cpu_int, cpu_ms) = found.expect("no seed passed the wall");

        let mut eng = match GrainEngine::new(init.len(), 1, lab.bins) {
            Ok(e) => e,
            Err(e) => { println!("no GPU — skipped ({e})"); return; }
        };
        println!("GPU: {}   grains {}   seed {}   steps {} × {dist}", eng.adapter_name, init.len(), lab.seed, steps);
        eng.set_lab(&lab);
        // pid 0 (first photon of the Sim); rand 0 → passes when p > 0, as the CPU did
        eng.launch(0, 0, &init, 0, lab.seed, &photon_elements(&lab, 0), [0.0; MAX_ELEMS]);

        let mut p_gpu = None;
        let mut max_dpos = 0.0f32;
        let t0 = std::time::Instant::now();
        for s in 0..steps {
            eng.tick(dist, 1);
            if s == wall_step {
                let st = eng.read_status()[0];
                assert_eq!(st.reach_mask & 1, 1, "GPU front must reach the wall on the same step");
                p_gpu = Some(st.num / st.den);
            }
            if s + 1 == mid_step {
                let g = eng.read_grains(0, init.len());
                let mut alive_mismatch = 0;
                let mut worst = 0;
                for (i, (c, q)) in cpu_mid.iter().zip(&g).enumerate() {
                    let d = (c.pos - bevy::math::Vec3::from_array(q.pos)).abs().max_element();
                    if d > max_dpos { worst = i; }
                    max_dpos = max_dpos.max(d);
                    if c.alive != q.alive() || c.diffracted != q.diffracted() { alive_mismatch += 1; }
                }
                let (c, q) = (&cpu_mid[worst], &g[worst]);
                println!("worst grain {worst}: CPU pos {:?} dir {:?} origin {:?}\n                 GPU pos {:?} dir {:?} origin {:?}",
                         c.pos, c.dir, c.origin, q.pos, q.dir, q.origin);
                let alive = g.iter().filter(|q| q.alive()).count();
                println!("after {mid_step} steps: alive {alive}, max |Δpos| = {max_dpos:.3e}, alive/diffracted mismatches {alive_mismatch}");
                assert_eq!(alive_mismatch, 0);
            }
        }
        eng.wait();
        let gpu_ms = t0.elapsed().as_secs_f64() * 1e3;
        let st = eng.read_status()[0];
        let screen = eng.read_screen(0);
        let gi = screen.intensity(true);
        let total: f64 = gi.iter().sum();
        let gpu_int: Vec<f32> = gi.iter().map(|v| (v / total) as f32).collect();
        let cpu_int: Vec<f32> = cpu_int.iter().map(|&v| v as f32).collect();
        let corr = correlation(&cpu_int, &gpu_int);
        let max_dint = cpu_int.iter().zip(&gpu_int).map(|(a, b)| (a - b).abs()).fold(0.0f32, f32::max);
        let peak = cpu_int.iter().cloned().fold(0.0f32, f32::max);
        let p_gpu = p_gpu.unwrap();
        println!("wall p_pass: CPU {p_cpu:.6}  GPU {p_gpu:.6}");
        println!("screen: {} grains landed, corr(CPU, GPU) = {corr:.7}, max |ΔI| = {max_dint:.2e} (peak {peak:.2e})",
                 screen.count.iter().sum::<u32>());
        println!("timing (debug build, incl. test readbacks): CPU Sim {cpu_ms:.1} ms, GPU {gpu_ms:.1} ms");

        assert_eq!(st.alive, 0, "all grains must have ended");
        assert!((p_cpu - p_gpu).abs() < 1e-4);
        assert!(max_dpos < 1e-3, "max |Δpos| {max_dpos}");
        assert!(corr > 0.999, "correlation {corr}");
    }

    /// Born rule on the GPU: a photon that draws "absorbed" dies whole at the
    /// wall and leaves nothing on the screen; a polarizer at 90° blocks H light.
    #[test]
    fn gpu_photon_level_decisions() {
        let mut lab = Lab::preset(Experiment::DoubleSlit);
        lab.grains = 5000;
        let mut sim = Sim::new(lab.clone());
        sim.launch();
        let init = sim.photons[0].grains.clone();
        let mut eng = match GrainEngine::new(8192, 2, lab.bins) {
            Ok(e) => e,
            Err(e) => { println!("no GPU — skipped ({e})"); return; }
        };
        eng.set_lab(&lab);
        eng.launch(0, 0, &init, 0, lab.seed, &photon_elements(&lab, 0), [0.999_999; MAX_ELEMS]);

        // Polarizers: H light, one polarizer at 90° → p = 0 → absorbed
        let mut lp = Lab::preset(Experiment::Polarizers);
        lp.grains = 3000;
        lp.polarizers[0].angle = std::f32::consts::FRAC_PI_2;
        let mut sp = Sim::new(lp.clone());
        sp.launch();
        let pinit = sp.photons[0].grains.clone();
        eng.launch(1, 5056, &pinit, 0, lp.seed, &photon_elements(&lp, 0), [0.5; MAX_ELEMS]);

        eng.tick(10.0, 40);
        let st = eng.read_status();
        println!("wall: fate {} alive {}   polarizer: fate {:#x} alive {}", st[0].fate, st[0].alive, st[1].fate, st[1].alive);
        assert_eq!(st[0].fate, FATE_WALL);
        assert_eq!(st[0].alive, 0);
        assert_eq!(eng.read_screen(0).count.iter().sum::<u32>(), 0);
        assert_eq!(st[1].fate, FATE_POLARIZER);
        assert_eq!(st[1].alive, 0);
    }

    /// Throughput with a million grains (prints numbers for docs/GPU_PORT.md).
    #[test]
    fn gpu_million_grains_throughput() {
        let n = 1_000_000;
        let mut lab = Lab::preset(Experiment::DoubleSlit);
        lab.grains = n;
        let mut sim = Sim::new(lab.clone());
        sim.launch();
        let init = sim.photons[0].grains.clone();
        let mut eng = match GrainEngine::new(n, 1, lab.bins) {
            Ok(e) => e,
            Err(e) => { println!("no GPU — skipped ({e})"); return; }
        };
        eng.set_lab(&lab);
        eng.launch(0, 0, &init, 0, lab.seed, &photon_elements(&lab, 0), [0.0; MAX_ELEMS]);
        eng.tick(1.0, 1);
        eng.wait();

        let ticks = 100;
        let t0 = std::time::Instant::now();
        eng.tick(2.0, ticks);
        eng.wait();
        let gpu = t0.elapsed().as_secs_f64() * 1e3 / ticks as f64;

        let t0 = std::time::Instant::now();
        for _ in 0..5 { sim.step(2.0); }
        let cpu = t0.elapsed().as_secs_f64() * 1e3 / 5.0;
        let screen = eng.read_screen(0);
        println!("{} grains, {} MiB of grain buffer: GPU {gpu:.2} ms/tick (probe+decide+advance), CPU Sim {cpu:.1} ms/step; landed so far {}",
                 n, n * 64 >> 20, screen.count.iter().sum::<u32>());
    }
}

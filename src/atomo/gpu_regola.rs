//! LA REGOLA DELLA LUCE SULLA GPU (studio motore_grafico, esame R1).
//!
//! La stessa regola zig-zag di `Fotone3D::tick` (fotone3d.rs), scritta in WGSL
//! (assets/atomo/regola_luce.wgsl): per ogni asse (x,y,z nei tick pari, z,y,x nei dispari) una dispatch che
//! divide F sugli autovettori di S_j e sposta la parte +1 avanti di una cella, la −1 indietro, la 0 ferma.
//! Due terne di buffer di storage a ping-pong (una per componente, vec2<f32> per cella), nessuna copia verso la CPU
//! durante i tick. Collegata alla scena del fotone da gpu_fotone.rs (palline e misure sulla GPU).
//!
//! Funziona in due modi, con lo stesso codice:
//! • dentro Bevy: `GpuRegola::nel_motore(&render_device, &render_queue, n)` (il device di Bevy);
//! • senza finestra (test): `GpuHeadless::nuovo()` crea Instance/Adapter/Device propri.
//! Tutti i metodi prendono `&wgpu::Device` / `&wgpu::Queue`, così la struttura non possiede il device.
//!
//! wgpu: la stessa wgpu 23 di Bevy 0.15 (dipendenza diretta già nel Cargo.toml, una sola versione nel
//! Cargo.lock, quindi i tipi coincidono con quelli di `bevy::render::render_resource`).

use crate::fotone3d::Cf;
use crate::walk::order;

const SORGENTE: &str = include_str!("../../assets/atomo/regola_luce.wgsl");
const WG: [u32; 3] = [8, 8, 4];

/// Device e queue propri, per i test e per l'uso fuori da Bevy.
pub struct GpuHeadless {
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    /// nome e backend dell'adattatore (es. "llvmpipe (LLVM 20.1.2, 256 bits) (Vulkan)")
    pub adattatore: String,
}

impl GpuHeadless {
    /// Cerca un adattatore (anche software, es. llvmpipe). `Err` con il motivo se non c'è.
    pub fn nuovo() -> Result<Self, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .ok_or("nessun adattatore GPU (neanche software) trovato")?;
        let info = adapter.get_info();
        let adattatore = format!("{} ({:?})", info.name, info.backend);
        // chiediamo i limiti dell'adattatore: a n = 192 un buffer pesa 170 MB (> 128 MB di default)
        let limits = adapter.limits();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("regola_luce"),
            required_features: wgpu::Features::empty(),
            required_limits: limits,
            memory_hints: wgpu::MemoryHints::Performance,
        }, None))
        .map_err(|e| format!("errore del device: {e}"))?;
        Ok(GpuHeadless { device, queue, adattatore })
    }
}

/// Il campo F di un fotone sulla GPU e la regola che lo fa avanzare.
///
/// Il campo è diviso in tre buffer, uno per componente (n³ × 8 byte ciascuno), e ce ne sono due terne a
/// ping-pong: a 192³ un buffer pesa 56 MB (dentro i 128 MB minimi di WebGPU per un binding), a 384³ 453 MB.
pub struct GpuRegola {
    pub n: usize,
    /// tick fatti (decide l'ordine degli assi, come `Fotone3D::ticks`)
    pub ticks: u64,
    pipeline: wgpu::ComputePipeline,
    campo: [[wgpu::Buffer; 3]; 2],
    /// bind group [asse][verso: 0 = +1, 1 = −1][sorgente: terna 0 o 1]
    gruppi: Vec<wgpu::BindGroup>,
    /// quale delle due terne contiene il campo attuale
    attuale: usize,
    staging: wgpu::Buffer,
}

impl GpuRegola {
    /// Byte di un buffer di campo (una componente) per lato n.
    pub fn byte_componente(n: usize) -> u64 { (n * n * n * 8) as u64 }
    /// Byte di tutto il campo (tre componenti) per lato n.
    pub fn byte_campo(n: usize) -> u64 { 3 * Self::byte_componente(n) }

    /// Il device regge un buffer di componente a questo n? (limiti di wgpu)
    pub fn entra(device: &wgpu::Device, n: usize) -> Result<(), String> {
        let lim = device.limits();
        let b = Self::byte_componente(n);
        if b > lim.max_storage_buffer_binding_size as u64 || b > lim.max_buffer_size {
            return Err(format!("una componente del campo a {n}³ pesa {} MB, il device ne regge {} MB per binding / {} MB per buffer",
                               b >> 20, lim.max_storage_buffer_binding_size >> 20, lim.max_buffer_size >> 20));
        }
        if lim.max_storage_buffers_per_shader_stage < 8 {
            return Err(format!("servono 8 buffer di storage per shader, il device ne ha {}", lim.max_storage_buffers_per_shader_stage));
        }
        Ok(())
    }

    /// Crea pipeline e buffer (campo a zero) sul device dato.
    pub fn nuovo(device: &wgpu::Device, n: usize) -> Self {
        let modulo = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("regola_luce"),
            source: wgpu::ShaderSource::Wgsl(SORGENTE.into()),
        });
        let voce = |binding, ty| wgpu::BindGroupLayoutEntry {
            binding, visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer { ty, has_dynamic_offset: false, min_binding_size: None },
            count: None,
        };
        let ro = wgpu::BufferBindingType::Storage { read_only: true };
        let rw = wgpu::BufferBindingType::Storage { read_only: false };
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("regola_luce_bgl"),
            entries: &[voce(0, wgpu::BufferBindingType::Uniform), voce(1, ro), voce(2, ro), voce(3, ro), voce(4, rw), voce(5, rw), voce(6, rw)],
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("regola_luce_layout"), bind_group_layouts: &[&bgl], push_constant_ranges: &[],
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("regola_luce"), layout: Some(&layout), module: &modulo, entry_point: Some("main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(), cache: None,
        });
        use wgpu::BufferUsages as U;
        let size = Self::byte_componente(n);
        let mk = |label| device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label), size, usage: U::STORAGE | U::COPY_DST | U::COPY_SRC, mapped_at_creation: false,
        });
        let campo = [[mk("campo_0x"), mk("campo_0y"), mk("campo_0z")], [mk("campo_1x"), mk("campo_1y"), mk("campo_1z")]];
        let staging = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("campo_staging"), size, usage: U::MAP_READ | U::COPY_DST, mapped_at_creation: false,
        });
        let mut gruppi = Vec::with_capacity(12);
        for ax in 0..3u32 {
            for verso in [1i32, -1] {
                let mut dati = [0u8; 16];
                dati[0..4].copy_from_slice(&(n as u32).to_le_bytes());
                dati[4..8].copy_from_slice(&ax.to_le_bytes());
                dati[8..12].copy_from_slice(&verso.to_le_bytes());
                let ub = device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("regola_luce_par"), size: 16, usage: U::UNIFORM, mapped_at_creation: true,
                });
                ub.slice(..).get_mapped_range_mut().copy_from_slice(&dati);
                ub.unmap();
                for s in 0..2 {
                    let (a, b) = (&campo[s], &campo[1 - s]);
                    gruppi.push(device.create_bind_group(&wgpu::BindGroupDescriptor {
                        label: Some("regola_luce_bg"),
                        layout: &bgl,
                        entries: &[
                            wgpu::BindGroupEntry { binding: 0, resource: ub.as_entire_binding() },
                            wgpu::BindGroupEntry { binding: 1, resource: a[0].as_entire_binding() },
                            wgpu::BindGroupEntry { binding: 2, resource: a[1].as_entire_binding() },
                            wgpu::BindGroupEntry { binding: 3, resource: a[2].as_entire_binding() },
                            wgpu::BindGroupEntry { binding: 4, resource: b[0].as_entire_binding() },
                            wgpu::BindGroupEntry { binding: 5, resource: b[1].as_entire_binding() },
                            wgpu::BindGroupEntry { binding: 6, resource: b[2].as_entire_binding() },
                        ],
                    }));
                }
            }
        }
        GpuRegola { n, ticks: 0, pipeline, campo, gruppi, attuale: 0, staging }
    }

    /// Dentro l'app Bevy: usa il device e la queue del render world.
    pub fn nel_motore(device: &bevy::render::renderer::RenderDevice, n: usize) -> Self {
        Self::nuovo(device.wgpu_device(), n)
    }

    /// Carica il campo (es. dopo `Fotone3D::lancia`) e imposta il contatore dei tick.
    pub fn carica(&mut self, queue: &wgpu::Queue, f: &[[Cf; 3]], ticks: u64) {
        assert_eq!(f.len(), self.n * self.n * self.n);
        self.attuale = 0;
        for c in 0..3 {
            let mut b = Vec::with_capacity(f.len() * 8);
            for p in f { b.extend_from_slice(&p[c].re.to_le_bytes()); b.extend_from_slice(&p[c].im.to_le_bytes()); }
            queue.write_buffer(&self.campo[0][c], 0, &b);
        }
        self.ticks = ticks;
    }

    /// Registra `quanti` tick nel command encoder (una compute pass per asse: la barriera tra assi è
    /// garantita dalla fine della pass). Non invia niente: la queue la gestisce chi chiama.
    pub fn registra_tick(&mut self, enc: &mut wgpu::CommandEncoder, quanti: u32) {
        let n = self.n as u32;
        let g = [n.div_ceil(WG[0]), n.div_ceil(WG[1]), n.div_ceil(WG[2])];
        for _ in 0..quanti {
            for ax in order(self.ticks) {
                let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some("regola_luce"), timestamp_writes: None });
                pass.set_pipeline(&self.pipeline);
                pass.set_bind_group(0, &self.gruppi[ax * 4 + self.attuale], &[]);
                pass.dispatch_workgroups(g[0], g[1], g[2]);
                drop(pass);
                self.attuale = 1 - self.attuale;
            }
            self.ticks += 1;
        }
    }

    /// Registra l'inverso esatto di `quanti` tick (come `Fotone3D::indietro`).
    pub fn registra_indietro(&mut self, enc: &mut wgpu::CommandEncoder, quanti: u32) {
        let n = self.n as u32;
        let g = [n.div_ceil(WG[0]), n.div_ceil(WG[1]), n.div_ceil(WG[2])];
        for _ in 0..quanti {
            if self.ticks == 0 { return; }
            self.ticks -= 1;
            for ax in order(self.ticks).into_iter().rev() {
                let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some("regola_luce_inv"), timestamp_writes: None });
                pass.set_pipeline(&self.pipeline);
                pass.set_bind_group(0, &self.gruppi[ax * 4 + 2 + self.attuale], &[]);
                pass.dispatch_workgroups(g[0], g[1], g[2]);
                drop(pass);
                self.attuale = 1 - self.attuale;
            }
        }
    }

    /// Comodità: registra e invia `quanti` tick.
    pub fn tick(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, quanti: u32) {
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("regola_luce") });
        self.registra_tick(&mut enc, quanti);
        queue.submit(Some(enc.finish()));
    }

    /// Quale terna contiene il campo attuale (0 o 1): cambia a ogni asse.
    pub fn attuale(&self) -> usize { self.attuale }

    /// Le due terne di buffer (x, y, z; un vec2<f32> per cella, cella = (i n + j) n + k): da legare ad altri
    /// shader senza passare per la CPU. Quella attuale è `terne()[attuale()]`.
    pub fn terne(&self) -> &[[wgpu::Buffer; 3]; 2] { &self.campo }

    /// La terna che contiene il campo attuale.
    pub fn campo_attuale(&self) -> &[wgpu::Buffer; 3] { &self.campo[self.attuale] }

    /// Rilegge il campo sulla CPU (bloccante: aspetta la GPU). Lento: serve per test, misure e viste di diagnosi.
    pub fn leggi(&self, device: &wgpu::Device, queue: &wgpu::Queue) -> Vec<[Cf; 3]> {
        let size = Self::byte_componente(self.n);
        let nc = self.n * self.n * self.n;
        let mut out = vec![[Cf::ZERO; 3]; nc];
        for c in 0..3 {
            let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("leggi_campo") });
            enc.copy_buffer_to_buffer(&self.campo[self.attuale][c], 0, &self.staging, 0, size);
            queue.submit(Some(enc.finish()));
            let slice = self.staging.slice(..);
            let (tx, rx) = std::sync::mpsc::channel();
            slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx.send(r); });
            device.poll(wgpu::Maintain::Wait);
            rx.recv().expect("map_async").expect("mappatura del buffer");
            {
                let dati = slice.get_mapped_range();
                let fl: &[f32] = bytemuck::cast_slice(&dati);
                for (q, p) in out.iter_mut().enumerate() { p[c] = Cf::new(fl[2 * q], fl[2 * q + 1]); }
            }
            self.staging.unmap();
        }
        out
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::fotone3d::Fotone3D;
    use crate::fotscena::Parametri;
    use std::time::Instant;

    fn diff_rel(a: &[[Cf; 3]], b: &[[Cf; 3]]) -> (f64, f64) {
        let (mut dmax, mut fmax) = (0.0f64, 0.0f64);
        for (p, q) in a.iter().zip(b) {
            for c in 0..3 {
                let (dr, di) = ((p[c].re - q[c].re) as f64, (p[c].im - q[c].im) as f64);
                dmax = dmax.max((dr * dr + di * di).sqrt());
                fmax = fmax.max((q[c].norm2() as f64).sqrt());
            }
        }
        (dmax / fmax.max(1e-300), fmax)
    }

    /// R1 (lab/out/mg_prereg.json): stesso impulso del Laboratorio (Parametri di default, n = 96),
    /// 96 tick su CPU e GPU, max|F_gpu − F_cpu| / max|F_cpu| ≤ 1e-4.
    #[test]
    fn gpu_regola_r1() {
        let gpu = match GpuHeadless::nuovo() {
            Ok(g) => g,
            Err(e) => { eprintln!("R1 SALTATO: {e} (serve Vulkan, anche software llvmpipe)"); return; }
        };
        eprintln!("adattatore: {}", gpu.adattatore);
        let par = Parametri::default();
        let n = par.n();
        let l = par.lancio();
        let mut cpu = Fotone3D::new(n);
        cpu.lancia(&l);
        let iniziale = cpu.f.clone();
        let tick = 96u32;

        // GPU: carica, un giro di riscaldamento (compilazione), poi la corsa misurata
        let mut g = GpuRegola::nuovo(&gpu.device, n);
        g.carica(&gpu.queue, &iniziale, 0);
        g.tick(&gpu.device, &gpu.queue, 1);
        let _ = g.leggi(&gpu.device, &gpu.queue);
        g.carica(&gpu.queue, &iniziale, 0);
        gpu.device.poll(wgpu::Maintain::Wait);
        let t0 = Instant::now();
        g.tick(&gpu.device, &gpu.queue, tick);
        gpu.device.poll(wgpu::Maintain::Wait);
        let ms_gpu = t0.elapsed().as_secs_f64() * 1e3 / tick as f64;
        let f_gpu = g.leggi(&gpu.device, &gpu.queue);

        // CPU
        let t0 = Instant::now();
        for _ in 0..tick { cpu.tick(); }
        let ms_cpu = t0.elapsed().as_secs_f64() * 1e3 / tick as f64;
        let (rel, fmax) = diff_rel(&f_gpu, &cpu.f);
        // controllo che il campo si sia davvero mosso (il confronto non è banale)
        let (mosso, _) = diff_rel(&cpu.f, &iniziale);
        assert!(mosso > 0.5, "il campo non si è mosso: {mosso:e}");

        // controllo in più: l'inverso sulla GPU riporta al lancio
        let mut enc = gpu.device.create_command_encoder(&Default::default());
        g.registra_indietro(&mut enc, tick);
        gpu.queue.submit(Some(enc.finish()));
        let (rel_inv, _) = diff_rel(&g.leggi(&gpu.device, &gpu.queue), &iniziale);

        let passa = rel <= 1e-4;
        eprintln!("R1: n = {n}, {tick} tick, max|F_cpu| = {fmax:.3e}");
        eprintln!("R1: max|F_gpu - F_cpu| / max|F_cpu| = {rel:.3e}  (soglia 1e-4) -> {}", if passa { "PASSA" } else { "FALLISCE" });
        eprintln!("R1: ms/tick CPU (rayon) = {ms_cpu:.2}, GPU ({}) = {ms_gpu:.2}", gpu.adattatore);
        eprintln!("R1: quanto e' cambiato il campo in {tick} tick: max|F_96 - F_0| / max|F| = {mosso:.3}");
        eprintln!("R1: avanti+indietro sulla GPU: max diff / max|F0| = {rel_inv:.3e}");
        let cpu_core = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(0);
        let json = format!(
            "{{\n \"studio\": \"motore_grafico\",\n \"esame\": \"R1_regola_gpu\",\n \"n\": {n},\n \"tick\": {tick},\n \
             \"impulso\": \"Parametri::default() di fotscena (8 celle per lambda, w0 = 1 lambda, circolare +, direzione x)\",\n \
             \"max_F_cpu\": {fmax:e},\n \"diff_rel_max\": {rel:e},\n \"soglia\": 1e-4,\n \"esito\": \"{}\",\n \
             \"ms_tick_cpu\": {ms_cpu:.3},\n \"cpu_core\": {cpu_core},\n \"ms_tick_gpu\": {ms_gpu:.3},\n \"adattatore\": \"{}\",\n \
             \"diff_rel_avanti_indietro_gpu\": {rel_inv:e},\n \"cambio_campo_96_tick\": {mosso:.4},\n \
             \"nota\": \"tempi del cloud: llvmpipe e' Vulkan su CPU, non dice niente di una GPU vera\"\n}}\n",
            if passa { "passa" } else { "fallisce" }, gpu.adattatore);
        let _ = std::fs::write(concat!(env!("CARGO_MANIFEST_DIR"), "/lab/out/mg_r1.json"), json);
        assert!(passa, "R1 fallito: {rel:e} > 1e-4");
    }
}

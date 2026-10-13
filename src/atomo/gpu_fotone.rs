//! IL FOTONE SULLA GPU (studio motore_grafico_gpu, docs/MOTORE_GRAFICO_GPU.md, pre-registrazione
//! lab/out/mg_gpu_prereg.json).
//!
//! Per i preset alta e ultra la scena del fotone tiene TUTTO sulla GPU: il campo (gpu_regola.rs), le palline di
//! energia (stessa regola di `fotdisegno::Traccianti::avanza`), le misure (energia, centro e covarianza, elicità,
//! flusso di Poynting totale, spin letto dalle palline) e le istanze del disegno, che vanno dal buffer della GPU
//! al disegno a istanze (istanze.rs) senza passare dalla CPU. Ogni tick la CPU rilegge solo qualche decina di KB
//! di somme parziali, mai il campo (170 MB a 192³).
//!
//! Regola di Valerio: i preset cambiano solo risoluzione, palline e precisione, mai i fenomeni. Per questo le
//! misure nuove (Poynting, spin dalle palline) sono calcolate anche sul cammino CPU (`poynting_cpu`,
//! `somme_palline_cpu`, `spin_da_somme`), con la stessa formula.
//!
//! Unica differenza dichiarata dalla CPU: una pallina persa (u < 1e-4 del picco) rinasce al posto di una pallina
//! viva scelta a caso (±1 cella), invece che dalla distribuzione |F|² per blocchi (semina) della CPU.

use std::sync::OnceLock;
use std::time::{Duration, Instant};

use bevy::prelude::{Res, Vec3};
use bevy::render::render_resource::Buffer as BBuffer;
use bevy::render::renderer::{RenderDevice, RenderQueue};
use bytemuck::{Pod, Zeroable};
use rayon::prelude::*;

use crate::fotdisegno::{Traccianti, Vista, L_SCIA, SCIA_MAX};
use crate::fotone3d::{base_trasversa, croce, punto, Fotone3D, Lancio, Misura};
use crate::gpu_regola::GpuRegola;
use crate::istanze::{Istanza, IstanzeGpu};
use crate::render::{schermo, Ball};
use crate::rng::Rng;

const SORGENTE: &str = include_str!("../../assets/atomo/gpu_fotone.wgsl");
/// numeri per gruppo di lavoro nelle somme del campo
const NP: usize = 20;
/// gruppi di lavoro per le somme delle palline
const NWG_SPIN: u32 = 64;
/// punti campionati per la CPU (i surfisti)
pub const MAX_SONDE: usize = 8;
const MORTA: u32 = u32::MAX;

/// Tick tra due campioni della scia: λ/3 (stessa formula di fotdisegno::passo_scia: i 4 campioni coprono lo
/// stesso tempo fisico, 4/3 λ/c, a ogni preset).
pub fn passo_scia(lambda_celle: f64) -> u64 { (lambda_celle / 3.0).round().max(1.0) as u64 }

// ───────────────────────── scelta del cammino ─────────────────────────

/// `--regola`: 0 = automatico (GPU per alta e ultra, CPU per bassa e media), 1 = sempre GPU, 2 = sempre CPU.
pub static REGOLA: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);

/// Il preset `qual` (indice di fotscena::QUALITA) usa la GPU?
pub fn vuole_gpu(qual: usize) -> bool {
    match REGOLA.load(std::sync::atomic::Ordering::Relaxed) { 1 => true, 2 => false, _ => qual >= 3 }
}

/// Device e queue del motore (Bevy), presi all'avvio: il thread del fotone li usa per i suoi compute shader.
pub static DISPOSITIVO: OnceLock<(RenderDevice, RenderQueue)> = OnceLock::new();

/// Sistema di avvio (main.rs): copia device e queue del motore per il thread del fotone.
pub fn prendi_dispositivo(d: Res<RenderDevice>, q: Res<RenderQueue>) { let _ = DISPOSITIVO.set((d.clone(), q.clone())); }

/// Aspetta (al massimo `attesa`) che il motore abbia un device.
pub fn dispositivo(attesa: Duration) -> Option<(RenderDevice, RenderQueue)> {
    let t0 = Instant::now();
    loop {
        if let Some(d) = DISPOSITIVO.get() { return Some(d.clone()); }
        if t0.elapsed() > attesa { return None; }
        std::thread::sleep(Duration::from_millis(20));
    }
}

// ───────────────────────── misure comuni ai due cammini ─────────────────────────

/// Flusso di Poynting totale Σ E × B sulle celle (cammino CPU).
pub fn poynting_cpu(f: &Fotone3D) -> [f64; 3] {
    f.f.par_iter().map(|p| {
        let (e, b) = ([p[0].re, p[1].re, p[2].re], [p[0].im, p[1].im, p[2].im]);
        [(e[1] * b[2] - e[2] * b[1]) as f64, (e[2] * b[0] - e[0] * b[2]) as f64, (e[0] * b[1] - e[1] * b[0]) as f64]
    }).reduce(|| [0.0; 3], |a, b| [a[0] + b[0], a[1] + b[1], a[2] + b[2]])
}

/// Somme delle palline (cammino CPU): Σ (p − rif) × v e Σ v.
pub fn somme_palline_cpu(pos: &[[f32; 3]], vel: &[[f32; 3]], rif: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let s = pos.par_iter().zip(vel.par_iter()).map(|(p, v)| {
        let r = [p[0] as f64 - rif[0], p[1] as f64 - rif[1], p[2] as f64 - rif[2]];
        let v = [v[0] as f64, v[1] as f64, v[2] as f64];
        let l = croce(r, v);
        [l[0], l[1], l[2], v[0], v[1], v[2]]
    }).reduce(|| [0.0; 6], |a, b| { let mut o = a; for q in 0..6 { o[q] += b[q]; } o });
    let s = s;
    ([s[0], s[1], s[2]], [s[3], s[4], s[5]])
}

/// Spin letto dalle palline, in ħ per fotone (come FOTONE_PALLINE F2): (2π/λ) ⟨((p − c) × v) · d⟩ con c il centro
/// dell'energia misurato. Dalle somme fatte attorno a `rif`: Σ (p − c) × v = Σ (p − rif) × v − (c − rif) × Σ v.
pub fn spin_da_somme(l_rif: [f64; 3], v_sum: [f64; 3], rif: [f64; 3], centro: [f64; 3], d: [f64; 3], lambda: f64, quante: usize) -> f64 {
    let dc = [centro[0] - rif[0], centro[1] - rif[1], centro[2] - rif[2]];
    let k = croce(dc, v_sum);
    let l = [l_rif[0] - k[0], l_rif[1] - k[1], l_rif[2] - k[2]];
    punto(l, d) / quante.max(1) as f64 * std::f64::consts::TAU / lambda
}

// ───────────────────────── la parte GPU ─────────────────────────

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Default)]
struct Par {
    n: u32, tick: u32, n_palline: u32, ns: u32,
    verso: f32, scia_op: u32, allunga: u32, n_sonde: u32,
    rif: [f32; 4], e1: [f32; 4], e2: [f32; 4], dir: [f32; 4],
}

/// Quello che la CPU rilegge a ogni tick (qualche decina di KB).
#[derive(Clone, Debug)]
pub struct Letture {
    pub misura: Misura,
    /// numeratore e denominatore dell'elicità (Fotone3D::elicita)
    pub elicita_num: f64,
    pub elicita_den: f64,
    pub poynting: [f64; 3],
    pub imax: f32,
    pub emax2: f32,
    /// Σ (p − rif) × v e Σ v sulle palline
    pub l_rif: [f64; 3],
    pub v_sum: [f64; 3],
    /// velocità (xyz) e densità (w) nei punti di sonda
    pub sonde: Vec<[f32; 4]>,
}

/// Campo, palline e misure di un fotone sulla GPU.
pub struct GpuFotone {
    pub regola: GpuRegola,
    pub n: usize,
    pub n_palline: usize,
    pub ns: usize,
    /// ogni quanti tick avanti la storia (scia) delle palline scorre: lo stesso di `Traccianti` (λ/3)
    pub passo_scia: u64,
    passi: u64,
    par: Par,
    par_buf: wgpu::Buffer,
    pos: wgpu::Buffer,
    vel: wgpu::Buffer,
    scia: wgpu::Buffer,
    /// tenuti vivi qui (li usano solo i bind group)
    _eprev: [wgpu::Buffer; 3],
    _stato: wgpu::Buffer,
    eta: wgpu::Buffer,
    sonde_in: wgpu::Buffer,
    staging: wgpu::Buffer,
    parz: wgpu::Buffer,
    parzp: wgpu::Buffer,
    sonde_out: wgpu::Buffer,
    /// le istanze del disegno (STORAGE | VERTEX): istanze.rs le disegna da qui
    pub ist: BBuffer,
    pip: [wgpu::ComputePipeline; 7],
    /// bind group per kernel; quelli che leggono il campo hanno due versioni (terna 0 o 1)
    bg_misura: [wgpu::BindGroup; 2],
    bg_avanza: [wgpu::BindGroup; 2],
    bg_sonde: [wgpu::BindGroup; 2],
    bg_istanze: [wgpu::BindGroup; 2],
    bg_finalizza: wgpu::BindGroup,
    bg_rinasci: wgpu::BindGroup,
    bg_spin: wgpu::BindGroup,
}

const K_MISURA: usize = 0;
const K_FINALIZZA: usize = 1;
const K_AVANZA: usize = 2;
const K_RINASCI: usize = 3;
const K_SONDE: usize = 4;
const K_SPIN: usize = 5;
const K_ISTANZE: usize = 6;

impl GpuFotone {
    /// Il device regge questo preset? (limiti di wgpu: niente sorprese a metà corsa)
    pub fn entra(device: &wgpu::Device, n: usize, n_palline: usize, ns: usize) -> Result<(), String> {
        GpuRegola::entra(device, n)?;
        let lim = device.limits();
        let bi = (n_palline * (ns + 1) * std::mem::size_of::<Istanza>()) as u64;
        if bi > lim.max_storage_buffer_binding_size as u64 || bi > lim.max_buffer_size {
            return Err(format!("le istanze ({n_palline} palline con scia {ns}) pesano {} MB, il device ne regge {} MB", bi >> 20, lim.max_storage_buffer_binding_size >> 20));
        }
        if (n_palline as u64).div_ceil(256) > lim.max_compute_workgroups_per_dimension as u64 { return Err("troppe palline per una dispatch".into()); }
        if ((n * n) as u64).div_ceil(64) > lim.max_compute_workgroups_per_dimension as u64 { return Err("scatola troppo grande per una dispatch".into()); }
        Ok(())
    }

    pub fn nuovo(device: &wgpu::Device, n: usize, n_palline: usize, ns: usize) -> Self {
        assert!(n % 2 == 0 && ns <= SCIA_MAX);
        let regola = GpuRegola::nuovo(device, n);
        // la lunghezza della storia di ogni pallina viene da fotdisegno (stesso passo del buffer caricato)
        let src = SORGENTE.replace("const SCIA_MAX: u32 = 4u;", &format!("const SCIA_MAX: u32 = {}u;", SCIA_MAX));
        let modulo = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("gpu_fotone"), source: wgpu::ShaderSource::Wgsl(src.into()) });
        let pip = ["misura", "finalizza", "avanza", "rinasci", "sonde", "spin", "istanze"].map(|e| device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(e), layout: None, module: &modulo, entry_point: Some(e),
            compilation_options: wgpu::PipelineCompilationOptions::default(), cache: None,
        }));
        use wgpu::BufferUsages as U;
        let mk = |label: &str, size: u64, usage: U| device.create_buffer(&wgpu::BufferDescriptor { label: Some(label), size: size.max(16), usage, mapped_at_creation: false });
        let st = U::STORAGE | U::COPY_DST | U::COPY_SRC;
        let nc = (n * n * n) as u64;
        let ngm = ((n * n) as u64).div_ceil(64);
        let par_buf = mk("fotone_par", std::mem::size_of::<Par>() as u64, U::UNIFORM | U::COPY_DST);
        let eprev = [mk("e_prima_x", nc * 4, st), mk("e_prima_y", nc * 4, st), mk("e_prima_z", nc * 4, st)];
        let parz = mk("parz", ngm * NP as u64 * 4, st);
        let stato = mk("stato", 16, st);
        let np = n_palline as u64;
        let pos = mk("pos", np * 16, st);
        let vel = mk("vel", np * 16, st);
        let scia = mk("scia", np * 16 * SCIA_MAX as u64, st);
        let eta = mk("eta", np * 4, st);
        let ist: BBuffer = mk("istanze_fotone", np * (ns as u64 + 1) * 32, U::STORAGE | U::VERTEX | U::COPY_SRC).into();
        let sonde_in = mk("sonde_in", MAX_SONDE as u64 * 16, U::STORAGE | U::COPY_DST);
        let sonde_out = mk("sonde_out", MAX_SONDE as u64 * 16, st);
        let parzp = mk("parzp", NWG_SPIN as u64 * 8 * 4, st);
        let staging = mk("fotone_staging", ngm * NP as u64 * 4 + NWG_SPIN as u64 * 32 + MAX_SONDE as u64 * 16, U::MAP_READ | U::COPY_DST);

        let bg = |k: usize, voci: &[(u32, &wgpu::Buffer)]| {
            let entries: Vec<wgpu::BindGroupEntry> = voci.iter().map(|(b, buf)| wgpu::BindGroupEntry { binding: *b, resource: buf.as_entire_binding() }).collect();
            device.create_bind_group(&wgpu::BindGroupDescriptor { label: Some("gpu_fotone"), layout: &pip[k].get_bind_group_layout(0), entries: &entries })
        };
        let terne = regola.terne();
        let ist_w: &wgpu::Buffer = &ist;
        let bg_misura = [0, 1].map(|s| bg(K_MISURA, &[(0, &par_buf), (1, &terne[s][0]), (2, &terne[s][1]), (3, &terne[s][2]), (4, &eprev[0]), (5, &eprev[1]), (6, &eprev[2]), (7, &parz)]));
        let bg_avanza = [0, 1].map(|s| bg(K_AVANZA, &[(0, &par_buf), (1, &terne[s][0]), (2, &terne[s][1]), (3, &terne[s][2]), (8, &stato), (9, &pos), (10, &vel), (11, &scia), (12, &eta)]));
        let bg_sonde = [0, 1].map(|s| bg(K_SONDE, &[(0, &par_buf), (1, &terne[s][0]), (2, &terne[s][1]), (3, &terne[s][2]), (14, &sonde_in), (15, &sonde_out)]));
        let bg_istanze = [0, 1].map(|s| bg(K_ISTANZE, &[(0, &par_buf), (1, &terne[s][0]), (2, &terne[s][1]), (3, &terne[s][2]), (8, &stato), (9, &pos), (11, &scia), (12, &eta), (13, ist_w)]));
        let bg_finalizza = bg(K_FINALIZZA, &[(0, &par_buf), (7, &parz), (8, &stato)]);
        let bg_rinasci = bg(K_RINASCI, &[(0, &par_buf), (9, &pos), (10, &vel), (11, &scia), (12, &eta)]);
        let bg_spin = bg(K_SPIN, &[(0, &par_buf), (9, &pos), (10, &vel), (16, &parzp)]);
        GpuFotone {
            regola, n, n_palline, ns, passo_scia: 2, passi: 0, par: Par { n: n as u32, n_palline: n_palline as u32, ns: ns as u32, ..Default::default() }, par_buf,
            pos, vel, scia, _eprev: eprev, _stato: stato, eta, sonde_in, staging, parz, parzp, sonde_out, ist, pip,
            bg_misura, bg_avanza, bg_sonde, bg_istanze, bg_finalizza, bg_rinasci, bg_spin,
        }
    }

    /// Carica campo e palline (dopo il lancio sulla CPU), poi una prima misura (serve il massimo al disegno).
    pub fn carica(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, f: &Fotone3D, tr: &Traccianti, l: &Lancio) {
        assert_eq!(tr.pos.len(), self.n_palline);
        self.regola.carica(queue, &f.f, f.ticks);
        let v4 = |v: &[[f32; 3]]| -> Vec<u8> { let x: Vec<[f32; 4]> = v.iter().map(|p| [p[0], p[1], p[2], 0.0]).collect(); bytemuck::cast_slice(&x).to_vec() };
        queue.write_buffer(&self.pos, 0, &v4(&tr.pos));
        queue.write_buffer(&self.vel, 0, &v4(&tr.vel));
        let eta: Vec<u32> = tr.eta.iter().map(|&e| e as u32).collect();
        queue.write_buffer(&self.eta, 0, bytemuck::cast_slice(&eta));
        let scia: Vec<[f32; 4]> = tr.scia.iter().flat_map(|s| s.iter().map(|p| [p[0], p[1], p[2], 0.0])).collect();
        queue.write_buffer(&self.scia, 0, bytemuck::cast_slice(&scia));
        self.passi = 0;
        let (e1, e2) = base_trasversa(l.dir);
        self.par.e1 = [e1[0] as f32, e1[1] as f32, e1[2] as f32, 0.0];
        self.par.e2 = [e2[0] as f32, e2[1] as f32, e2[2] as f32, 0.0];
        self.par.dir = [l.dir[0] as f32, l.dir[1] as f32, l.dir[2] as f32, 0.0];
        self.par.verso = 0.0;
        self.par.scia_op = 0;
        self.par.tick = f.ticks as u32;
        self.par.n_sonde = 0;
        queue.write_buffer(&self.par_buf, 0, bytemuck::bytes_of(&self.par));
        // una prima misura: scrive e_prima = E del lancio e i massimi (servono al disegno prima del primo tick)
        let mut enc = device.create_command_encoder(&Default::default());
        self.registra_misura(&mut enc);
        queue.submit(Some(enc.finish()));
    }

    fn dispatch(&self, enc: &mut wgpu::CommandEncoder, k: usize, gruppi: u32, terna: usize) {
        let bg = match k {
            K_MISURA => &self.bg_misura[terna],
            K_AVANZA => &self.bg_avanza[terna],
            K_SONDE => &self.bg_sonde[terna],
            K_ISTANZE => &self.bg_istanze[terna],
            K_FINALIZZA => &self.bg_finalizza,
            K_RINASCI => &self.bg_rinasci,
            _ => &self.bg_spin,
        };
        let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some("gpu_fotone"), timestamp_writes: None });
        pass.set_pipeline(&self.pip[k]);
        pass.set_bind_group(0, bg, &[]);
        pass.dispatch_workgroups(gruppi, 1, 1);
    }

    fn registra_misura(&self, enc: &mut wgpu::CommandEncoder) {
        let t = self.regola.attuale();
        self.dispatch(enc, K_MISURA, ((self.n * self.n) as u32).div_ceil(64), t);
        self.dispatch(enc, K_FINALIZZA, 1, t);
    }

    /// Un tick avanti (o indietro) di campo e palline, con le misure. `rif` = centro atteso al tick NUOVO
    /// (le somme dei momenti sono fatte attorno a lui, come `Fotone3D::misura`). `sonde`: punti di cui si vuole
    /// la velocità dell'energia (i surfisti). Bloccante: aspetta la GPU e rilegge solo le somme.
    pub fn passo(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, avanti: bool, rif: [f64; 3], sonde: &[[f32; 3]]) -> Letture {
        self.passi += 1;
        self.par.verso = if avanti { 1.0 } else { -1.0 };
        self.par.scia_op = if avanti && self.passi % self.passo_scia.max(1) == 0 { 1 } else if !avanti { 2 } else { 0 };
        self.par.rif = [rif[0] as f32, rif[1] as f32, rif[2] as f32, 0.0];
        let ns = sonde.len().min(MAX_SONDE);
        self.par.n_sonde = ns as u32;
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("fotone_passo") });
        if avanti { self.regola.registra_tick(&mut enc, 1); } else { self.regola.registra_indietro(&mut enc, 1); }
        self.par.tick = self.regola.ticks as u32;
        queue.write_buffer(&self.par_buf, 0, bytemuck::bytes_of(&self.par));
        if ns > 0 {
            let s: Vec<[f32; 4]> = sonde[..ns].iter().map(|p| [p[0], p[1], p[2], 0.0]).collect();
            queue.write_buffer(&self.sonde_in, 0, bytemuck::cast_slice(&s));
        }
        let t = self.regola.attuale();
        let gp = (self.n_palline as u32).div_ceil(256);
        self.registra_misura(&mut enc);
        self.dispatch(&mut enc, K_AVANZA, gp, t);
        self.dispatch(&mut enc, K_RINASCI, gp, t);
        if ns > 0 { self.dispatch(&mut enc, K_SONDE, 1, t); }
        self.dispatch(&mut enc, K_SPIN, NWG_SPIN, t);
        let ngm = ((self.n * self.n) as u64).div_ceil(64);
        let b0 = ngm * NP as u64 * 4;
        let b1 = NWG_SPIN as u64 * 32;
        enc.copy_buffer_to_buffer(&self.parz, 0, &self.staging, 0, b0);
        enc.copy_buffer_to_buffer(&self.parzp, 0, &self.staging, b0, b1);
        enc.copy_buffer_to_buffer(&self.sonde_out, 0, &self.staging, b0 + b1, MAX_SONDE as u64 * 16);
        queue.submit(Some(enc.finish()));
        let v = leggi_staging(device, &self.staging);
        let f: &[f32] = bytemuck::cast_slice(&v);
        let mut s = [0.0f64; NP];
        let (mut imax, mut emax2) = (0.0f32, 0.0f32);
        for g in 0..ngm as usize {
            let p = &f[g * NP..(g + 1) * NP];
            for q in 0..15 { s[q] += p[q] as f64; }
            imax = imax.max(p[15]);
            emax2 = emax2.max(p[16]);
        }
        let o = ngm as usize * NP;
        let mut sp = [0.0f64; 6];
        for g in 0..NWG_SPIN as usize { for q in 0..6 { sp[q] += f[o + g * 8 + q] as f64; } }
        let o2 = o + NWG_SPIN as usize * 8;
        let sonde_v = (0..ns).map(|i| [f[o2 + 4 * i], f[o2 + 4 * i + 1], f[o2 + 4 * i + 2], f[o2 + 4 * i + 3]]).collect();
        Letture {
            misura: misura_da_somme(&s, rif), elicita_num: s[10], elicita_den: s[11], poynting: [s[12], s[13], s[14]],
            imax, emax2, l_rif: [sp[0], sp[1], sp[2]], v_sum: [sp[3], sp[4], sp[5]], sonde: sonde_v,
        }
    }

    /// Scrive le istanze del disegno dal buffer delle palline (niente CPU). `l_max`: allungamento massimo lungo E
    /// (celle), 0 = palline tonde.
    pub fn disegna(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, l_max: f32) {
        self.par.allunga = (l_max > 0.0) as u32;
        self.par.e1[3] = l_max;
        queue.write_buffer(&self.par_buf, 0, bytemuck::bytes_of(&self.par));
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("fotone_istanze") });
        self.dispatch(&mut enc, K_ISTANZE, (self.n_palline as u32).div_ceil(256), self.regola.attuale());
        queue.submit(Some(enc.finish()));
    }

    /// Le istanze da disegnare (buffer + quante).
    pub fn istanze(&self) -> IstanzeGpu { IstanzeGpu { buffer: self.ist.clone(), n: (self.n_palline * (self.ns + 1)) as u32 } }

    /// Rilegge posizioni, velocità ed età delle palline (test).
    pub fn leggi_palline(&self, device: &wgpu::Device, queue: &wgpu::Queue) -> (Vec<[f32; 3]>, Vec<[f32; 3]>, Vec<u32>) {
        let rd = |b: &wgpu::Buffer| -> Vec<u8> {
            let st = device.create_buffer(&wgpu::BufferDescriptor { label: None, size: b.size(), usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
            let mut enc = device.create_command_encoder(&Default::default());
            enc.copy_buffer_to_buffer(b, 0, &st, 0, b.size());
            queue.submit(Some(enc.finish()));
            leggi_staging(device, &st)
        };
        let p: Vec<[f32; 4]> = bytemuck::cast_slice(&rd(&self.pos)).to_vec();
        let v: Vec<[f32; 4]> = bytemuck::cast_slice(&rd(&self.vel)).to_vec();
        let e: Vec<u32> = bytemuck::cast_slice(&rd(&self.eta)).to_vec();
        (p.iter().map(|q| [q[0], q[1], q[2]]).collect(), v.iter().map(|q| [q[0], q[1], q[2]]).collect(), e)
    }

    /// Rilegge le istanze (test e confronto col disegno della CPU).
    pub fn leggi_istanze(&self, device: &wgpu::Device, queue: &wgpu::Queue) -> Vec<Istanza> {
        let b: &wgpu::Buffer = &self.ist;
        let st = device.create_buffer(&wgpu::BufferDescriptor { label: None, size: b.size(), usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let mut enc = device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(b, 0, &st, 0, b.size());
        queue.submit(Some(enc.finish()));
        bytemuck::cast_slice(&leggi_staging(device, &st)).to_vec()
    }
}

fn leggi_staging(device: &wgpu::Device, b: &wgpu::Buffer) -> Vec<u8> {
    let slice = b.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx.send(r); });
    device.poll(wgpu::Maintain::Wait);
    rx.recv().expect("map_async").expect("mappatura del buffer");
    let v = slice.get_mapped_range().to_vec();
    b.unmap();
    v
}

/// Centro e covarianza dalle somme (stesse formule di `Fotone3D::misura`).
fn misura_da_somme(s: &[f64; NP], rif: [f64; 3]) -> Misura {
    let w = s[0].max(1e-300);
    let m = [s[1] / w, s[2] / w, s[3] / w];
    let c2 = |q: usize, a: usize, b: usize| s[q] / w - m[a] * m[b];
    let cov = [[c2(4, 0, 0), c2(7, 0, 1), c2(8, 0, 2)], [c2(7, 0, 1), c2(5, 1, 1), c2(9, 1, 2)], [c2(8, 0, 2), c2(9, 1, 2), c2(6, 2, 2)]];
    Misura { energia: s[0], centro: [rif[0] + m[0], rif[1] + m[1], rif[2] + m[2]], cov }
}

// ───────────────────────── la scena: il fotone del Laboratorio sulla GPU ─────────────────────────

/// Quello che la scena del fotone (fotscena.rs) tiene quando la regola è sulla GPU.
pub struct FotoneGpu {
    pub g: GpuFotone,
    dev: (RenderDevice, RenderQueue),
    /// il campo in `Fotone3D::f` è quello del tick attuale? (si rilegge solo se una vista di diagnosi lo chiede)
    campo_cpu_ok: bool,
    rng: Rng,
    pub ultime: Option<Letture>,
}

/// Misure di un passo sulla GPU, pronte per la scena.
pub struct PassoGpu {
    pub misura: Misura,
    pub elicita: Option<f64>,
    pub poynting: [f64; 3],
    pub spin: f64,
}

impl FotoneGpu {
    /// Prova a mettere sulla GPU del motore il fotone appena lanciato (campo in `f`, palline in `tr`).
    /// `riuso`: il FotoneGpu del fotone precedente (stessi buffer se il preset è lo stesso).
    pub fn prova(f: &Fotone3D, tr: &Traccianti, l: &Lancio, riuso: Option<FotoneGpu>) -> Result<Self, String> {
        // la stessa scia a ogni preset (regola di Valerio, fotdisegno::palline): SCIA_MAX campioni
        let (n, np, ns) = (f.n, tr.pos.len(), SCIA_MAX);
        let mut fg = match riuso {
            Some(r) if r.g.n == n && r.g.n_palline == np && r.g.ns == ns => r,
            vecchio => {
                drop(vecchio);
                let dev = dispositivo(Duration::from_secs(30)).ok_or("il motore non ha un device GPU")?;
                GpuFotone::entra(dev.0.wgpu_device(), n, np, ns)?;
                let g = GpuFotone::nuovo(dev.0.wgpu_device(), n, np, ns);
                FotoneGpu { g, dev, campo_cpu_ok: true, rng: Rng::new(777), ultime: None }
            }
        };
        let (d, q) = (fg.dev.0.clone(), fg.dev.1.clone());
        fg.g.carica(d.wgpu_device(), &q, f, tr, l);
        fg.g.passo_scia = passo_scia(l.lambda);
        fg.campo_cpu_ok = true;
        fg.ultime = None;
        Ok(fg)
    }

    /// Un tick sulla GPU. Aggiorna `f.ticks` (il campo in `f.f` resta indietro finché non serve) e i surfisti di `tr`.
    pub fn passo(&mut self, f: &mut Fotone3D, tr: &mut Traccianti, avanti: bool, rif: [f64; 3], l: &Lancio) -> PassoGpu {
        let (d, q) = (self.dev.0.clone(), self.dev.1.clone());
        let sonde: Vec<[f32; 3]> = tr.surf.iter().take(MAX_SONDE).map(|s| *s.last().unwrap()).collect();
        let r = self.g.passo(d.wgpu_device(), &q, avanti, rif, &sonde);
        f.ticks = self.g.regola.ticks;
        self.campo_cpu_ok = false;
        // i surfisti: stessa logica di Traccianti::avanza (la velocità viene dalla GPU)
        let verso = if avanti { 1.0f32 } else { -1.0 };
        let imax = r.imax.max(1e-30);
        for (k, s) in r.sonde.iter().enumerate() {
            let p = *tr.surf[k].last().unwrap();
            if s[3] < 0.05 * imax {
                let rr = [self.rng.gauss() as f32 * 2.5, self.rng.gauss() as f32 * 2.5, self.rng.gauss() as f32 * 2.5];
                tr.surf[k] = vec![[rif[0] as f32 + rr[0], rif[1] as f32 + rr[1], rif[2] as f32 + rr[2]]];
                continue;
            }
            let qq = [p[0] + verso * s[0], p[1] + verso * s[1], p[2] + verso * s[2]];
            if avanti { tr.surf[k].push(qq); if tr.surf[k].len() > L_SCIA { tr.surf[k].remove(0); } }
            else if tr.surf[k].len() > 1 { tr.surf[k].pop(); } else { tr.surf[k] = vec![qq]; }
        }
        let elicita = avanti.then(|| r.elicita_num / ((std::f64::consts::TAU / l.lambda).sin() * r.elicita_den.max(1e-300)));
        let spin = spin_da_somme(r.l_rif, r.v_sum, rif, r.misura.centro, l.dir, l.lambda, self.g.n_palline);
        let out = PassoGpu { misura: r.misura, elicita, poynting: r.poynting, spin };
        self.ultime = Some(r);
        out
    }

    /// Rilegge il campo sulla CPU se serve (viste di diagnosi, clic del rivelatore). Lento: 24 n³ byte.
    pub fn campo_cpu(&mut self, f: &mut Fotone3D) {
        if self.campo_cpu_ok { return; }
        let (d, q) = (self.dev.0.clone(), self.dev.1.clone());
        f.f = self.g.regola.leggi(d.wgpu_device(), &q);
        self.campo_cpu_ok = true;
    }

    /// Il disegno: le palline di energia restano sulla GPU (istanze), la CPU fa solo surfisti, bordi della scatola
    /// e, se accese, le viste di diagnosi (che rileggono il campo).
    pub fn palline(&mut self, f: &mut Fotone3D, l: &Lancio, centro: [f64; 3], tr: &Traccianti, vista: Vista, surfisti: bool)
        -> (Vec<Ball>, Vec<Ball>, Vec<[f32; 3]>, Option<IstanzeGpu>) {
        let diagnosi = vista.fronti || vista.e_vettori || vista.b_vettori || vista.poynting || vista.spin || vista.sezione;
        let mut out = if diagnosi {
            self.campo_cpu(f);
            let v = Vista { traccianti: false, ..vista };
            crate::fotdisegno::palline(f, l, centro, tr, None, v, false).0
        } else { bordi(f.n) };
        // no separate surfers with trails (Valerio, 13/10): the energy balls already ride the flow
        if vista.traccianti && surfisti && false { disegna_surfisti(&mut out, tr); }
        if !vista.traccianti { return (out, vec![], vec![], None); }
        let (d, q) = (self.dev.0.clone(), self.dev.1.clone());
        let l_max = if vista.allunga { (0.35 * l.lambda) as f32 } else { 0.0 };
        self.g.disegna(d.wgpu_device(), &q, l_max);
        (out, vec![], vec![], Some(self.g.istanze()))
    }
}

/// I surfisti come in fotdisegno::palline.
fn disegna_surfisti(out: &mut Vec<Ball>, tr: &Traccianti) {
    let sp = |p: &[f32; 3]| schermo([p[0] as f64, p[1] as f64, p[2] as f64]).to_array();
    for s in &tr.surf {
        let k = s.len();
        for (j, p) in s.iter().enumerate() {
            let x = (j + 1) as f32 / k as f32;
            out.push((sp(p), 0.5 + 0.7 * x, [1.0, 0.95, 0.85, 0.12 + 0.6 * x * x]));
        }
        if let Some(p) = s.last() {
            out.push((sp(p), 2.6, [1.0, 0.97, 0.9, 1.0]));
            out.push((sp(p), 6.0, [1.0, 0.9, 0.7, 0.18]));
        }
    }
}

/// La scatola: palline tenui lungo i 12 spigoli (come in fotdisegno::palline).
fn bordi(n: usize) -> Vec<Ball> {
    let h = (n / 2) as f32;
    let (lo, hi) = (-h - 0.5, h - 0.5);
    let steps = 96;
    let mut out = Vec::with_capacity(12 * (steps + 1));
    for e in 0..12 {
        for s in 0..=steps {
            let t = lo + (hi - lo) * s as f32 / steps as f32;
            let (u, v) = ([lo, hi][e & 1], [lo, hi][(e >> 1) & 1]);
            let p = match e / 4 { 0 => [t, u, v], 1 => [u, t, v], _ => [u, v, t] };
            out.push((Vec3::new(p[0], p[2], -p[1]).to_array(), 0.6, [0.45, 0.55, 0.75, 0.09]));
        }
    }
    out
}

/// Una GPU senza finestra per i test (stessa di gpu_regola).
#[cfg(test)]
pub use crate::gpu_regola::GpuHeadless;

//! THE COMPUTED SCENES of the program (docs/SCENE_CALCOLATE.md): states and fields the studies computed, drawn with the
//! core glows (bagliori.rs), the same look as the zig-zag atom. Nothing here adds physics: it reads the exports of
//! assets/scene_calcolate/ and checks that the drawing tells the truth (exams S0-S9 of lab/out/sc_prereg.json).
//!
//!   • helium, H₂ parts A and B: walkers on |Ψ(r₁, r₂)|² inside the program (Metropolis, as lab/h2v_sim.py), each electron a
//!     glow with its short tail. The motion is the sampler, NOT a trajectory (a real ground state has no current).
//!     Part B: a reference electron on a chosen loop; the light of the balls is weighted by g(r₂|r₁) = |Ψ|²/(ρ₁ρ₁).
//!   • H₂ part C: the bond vibrates. Independent photographs of |Ψ_R|², R drawn with the weights of the (v₀+v₁)/√2 packet,
//!     protons at the mean R(t). No tails: the photographs are independent, a tail would be a fake motion.

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use rayon::prelude::*;
use serde_json::Value;

use crate::bagliori::{scia_in, Livello};
use crate::grafica::{self, VistaWorker};
use crate::render::{schermo, Ball};
use crate::rng::Rng;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mostra { Elio, H2 }
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Parte { A, B, C }

impl Mostra {
    pub fn nome(self) -> &'static str { match self { Mostra::Elio => "elio", Mostra::H2 => "h2" } }
    pub fn parti(self) -> &'static [(Parte, &'static str)] {
        match self {
            Mostra::Elio => &[(Parte::A, "la nuvola"), (Parte::B, "buca di Coulomb")],
            Mostra::H2 => &[(Parte::A, "la molecola"), (Parte::B, "dove sta l'altro elettrone"), (Parte::C, "il legame vibra")],
        }
    }
}

/// Screen frames per second the scenes are written for, and frames of one turn of parts B and C (lab/h2v_sim.py, he_v_render.py).
pub const RITMO: f64 = 24.0;
const GIRO_B: [u64; 2] = [240, 288];
const PERIODO_C_S: f64 = 4.0;
/// Pair walkers (two glows each) and the Metropolis of the video: 0.03 a₀ × 3 moves per frame (h2v_sim.py).
pub const N_CAMMINATORI: usize = 2500;
const PASSO: f64 = 0.03;
const MOSSE: usize = 3;
/// Photographs alive in part C and their life (frames).
const N_FOTO_C: usize = 2000;
const VITA_C: u32 = 6;
/// Past positions kept for a tail, and the most tail pieces per frame.
const S_MAX: usize = 16;
const MAX_SCIA: f64 = 60_000.0;
/// How fast a tail's smoothed path follows its glow, per frame.
const LISCIA: f32 = 0.3;
const STOPS: [[f32; 3]; 5] = [[0.20, 0.95, 1.00], [0.25, 0.45, 1.00], [0.70, 0.30, 1.00], [1.00, 0.35, 0.70], [1.00, 0.85, 0.30]];

fn tinta(angle: f32) -> [f32; 3] {
    let t = angle.rem_euclid(std::f32::consts::TAU) / std::f32::consts::TAU;
    let x = t * 5.0;
    let i = (x.floor() as usize) % 5;
    let s = x - x.floor();
    let (p, q) = (STOPS[i], STOPS[(i + 1) % 5]);
    [p[0] + (q[0] - p[0]) * s, p[1] + (q[1] - p[1]) * s, p[2] + (q[2] - p[2]) * s]
}

/// Ψ(r₁, r₂) = Σ w exp(−(G₀ r₁·r₁ + G₁ r₁·r₂ + G₂ r₂·r₂ + G₃ z₁ + G₄ z₂ + G₅)), normalised (h2v_core.WF).
pub struct Stato { pub g: Vec<[f64; 6]>, pub w: Vec<f64> }

impl Stato {
    #[inline]
    pub fn psi(&self, a: &[f64; 3], b: &[f64; 3]) -> f64 {
        let f = [a[0] * a[0] + a[1] * a[1] + a[2] * a[2], a[0] * b[0] + a[1] * b[1] + a[2] * b[2], b[0] * b[0] + b[1] * b[1] + b[2] * b[2], a[2], b[2], 1.0];
        let mut s = 0.0;
        for (g, w) in self.g.iter().zip(&self.w) {
            s += w * (-(f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3] + f[4] * g[4] + f[5] * g[5])).exp();
        }
        s
    }
}

/// One exported state family: the states on the R grid, their photographs, ρ₁ at the R of parts A and B, the packet.
pub struct Dati {
    pub mostra: Mostra,
    pub r: Vec<f64>,
    pub stati: Vec<Stato>,
    foto: Vec<f32>,
    n_foto: usize,
    pub i_rho: usize,
    rho: Vec<f32>,
    rmax: f64,
    nr: usize,
    nz: usize,
    pub rho_max: f64,
    pub z_nuclei: Vec<u32>,
    w0: Vec<f64>,
    w1: Vec<f64>,
    pub r0: f64,
    pub r1: f64,
    pub periodo_fs: f64,
    pub abs_z: Vec<f64>,
    pub rif: Value,
    pub energia: f64,
}

fn leggi_byte(p: &std::path::Path) -> Result<Vec<u8>, String> { std::fs::read(p).map_err(|e| format!("{}: {e}", p.display())) }

impl Dati {
    pub fn carica(mostra: Mostra) -> Result<Dati, String> {
        let d = crate::esportazioni::cartella_assets().join("scene_calcolate").join(mostra.nome());
        let m: Value = serde_json::from_slice(&leggi_byte(&d.join("manifest.json"))?).map_err(|e| e.to_string())?;
        let v = |x: &Value| -> Vec<f64> { x.as_array().map(|a| a.iter().filter_map(|y| y.as_f64()).collect()).unwrap_or_default() };
        let r = v(&m["R"]);
        let n_foto = m["n_foto"].as_u64().ok_or("n_foto mancante")? as usize;
        let termini: Vec<usize> = v(&m["n_termini"]).iter().map(|x| *x as usize).collect();
        let byte = leggi_byte(&d.join("stati.f64"))?;
        if byte.len() != termini.iter().sum::<usize>() * 56 { return Err(format!("{}: stati.f64 di {} byte", d.display(), byte.len())); }
        let num: Vec<f64> = byte.chunks_exact(8).map(|b| f64::from_le_bytes(b.try_into().unwrap())).collect();
        let (mut stati, mut k) = (vec![], 0);
        for t in &termini {
            let mut s = Stato { g: Vec::with_capacity(*t), w: Vec::with_capacity(*t) };
            for q in num[k..k + 7 * t].chunks_exact(7) { s.g.push([q[0], q[1], q[2], q[3], q[4], q[5]]); s.w.push(q[6]); }
            k += 7 * t;
            stati.push(s);
        }
        let f32s = |b: Vec<u8>| -> Vec<f32> { b.chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap())).collect() };
        let foto = f32s(leggi_byte(&d.join("fotografie.f32"))?);
        if foto.len() != r.len() * n_foto * 6 { return Err(format!("{}: fotografie.f32 di {} numeri", d.display(), foto.len())); }
        let rho = f32s(leggi_byte(&d.join("rho1.f32"))?);
        let (nr, nz) = (m["rho"]["nr"].as_u64().unwrap_or(0) as usize, m["rho"]["nz"].as_u64().unwrap_or(0) as usize);
        if rho.len() != nr * nz { return Err(format!("{}: rho1.f32 di {} numeri", d.display(), rho.len())); }
        Ok(Dati {
            mostra, stati, foto, n_foto, rho, nr, nz, rmax: m["rho"]["rmax"].as_f64().unwrap_or(6.0), rho_max: m["rho"]["massimo"].as_f64().unwrap_or(1.0),
            i_rho: m["i_rho"].as_u64().unwrap_or(0) as usize, z_nuclei: v(&m["Z"]).iter().map(|x| *x as u32).collect(),
            w0: v(&m["w0"]), w1: v(&m["w1"]), r0: m["R0"].as_f64().unwrap_or(0.0), r1: m["R1"].as_f64().unwrap_or(0.0),
            periodo_fs: m["periodo_fs"].as_f64().unwrap_or(0.0), abs_z: v(&m["abs_z"]), rif: m["rif"].clone(), energia: m["energia_Eh"].as_f64().unwrap_or(0.0), r,
        })
    }

    /// Photograph i of R index k: the two electrons.
    pub fn fotografia(&self, k: usize, i: usize) -> [[f64; 3]; 2] {
        let o = (k * self.n_foto + i) * 6;
        let f = &self.foto[o..o + 6];
        [[f[0] as f64, f[1] as f64, f[2] as f64], [f[3] as f64, f[4] as f64, f[5] as f64]]
    }

    /// ρ₁ at a point (bilinear in the table, clamped at its edge), R of parts A and B.
    pub fn rho1(&self, p: &[f64; 3]) -> f64 {
        let (dr, dz) = (self.rmax / (self.nr - 1) as f64, 2.0 * self.rmax / (self.nz - 1) as f64);
        let x = ((p[0] * p[0] + p[1] * p[1]).sqrt().min(self.rmax) / dr).min((self.nr - 1) as f64 - 1e-9);
        let y = ((p[2].clamp(-self.rmax, self.rmax) + self.rmax) / dz).min((self.nz - 1) as f64 - 1e-9);
        let (i, j) = (x as usize, y as usize);
        let (s, t) = (x - i as f64, y - j as f64);
        let q = |a: usize, b: usize| self.rho[a * self.nz + b] as f64;
        ((1.0 - s) * ((1.0 - t) * q(i, j) + t * q(i, j + 1)) + s * ((1.0 - t) * q(i + 1, j) + t * q(i + 1, j + 1))).max(1e-300)
    }

    /// g(r₂|r₁) of parts B.
    pub fn g(&self, r1: &[f64; 3], rho_r1: f64, p: &[f64; 3]) -> f64 {
        let s = self.stati[self.i_rho].psi(r1, p);
        s * s / (rho_r1 * self.rho1(p))
    }

    /// Weights of the R grid at phase cos(ωt) of the packet, and the mean R.
    pub fn pesi(&self, cos_wt: f64) -> (Vec<f64>, f64) {
        let mut w: Vec<f64> = self.w0.iter().zip(&self.w1).map(|(a, b)| (a + b * cos_wt).max(0.0)).collect();
        let s: f64 = w.iter().sum();
        w.iter_mut().for_each(|x| *x /= s.max(1e-300));
        (w, self.r0 + self.r1 * cos_wt)
    }
}

/// The reference electron on its chosen loop (declared), u in [0, 1]: the loops of lab/h2v_sim.py and lab/he_v_render.py.
pub fn riferimento(m: Mostra, u: f64) -> [f64; 3] {
    let a = std::f64::consts::TAU * u;
    match m {
        Mostra::H2 => [0.55 * (2.0 * a).sin(), 0.0, -1.25 * a.cos()],
        Mostra::Elio => {
            let v = [a.sin(), 0.35 * (2.0 * a).sin(), a.cos()];
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            [v[0] / n, v[1] / n, v[2] / n]
        }
    }
}

#[derive(Clone, Copy)]
struct Foto { k: usize, i: usize, rot: f64, eta: u32 }

/// The live scene: the walkers of parts A and B, the photographs of part C, the clock.
pub struct Mondo {
    pub dati: Arc<Dati>,
    pub cammina: Vec<[[f64; 3]; 2]>,
    p2: Vec<f64>,
    rngs: Vec<Rng>,
    /// past screen positions of each electron glow, newest at `testa`
    storia: Vec<[[f32; 3]; S_MAX]>,
    testa: usize,
    pieni: usize,
    foto: Vec<Foto>,
    rng: Rng,
    pub parte: Parte,
    /// frames of the present part
    pub f: u64,
    /// mean of g over the balls in the last frame of part B (exam S5)
    pub media_g: f64,
    pub vista: VistaWorker,
    pub ritmo_vero: f64,
}

const BLOCCO: usize = 128;

impl Mondo {
    pub fn new(dati: Arc<Dati>, n: usize) -> Self {
        let k = dati.i_rho;
        let n_foto = dati.n_foto;
        let cammina: Vec<[[f64; 3]; 2]> = (0..n).map(|i| dati.fotografia(k, i % n_foto)).collect();
        let s = &dati.stati[k];
        let p2: Vec<f64> = cammina.par_iter().map(|c| { let v = s.psi(&c[0], &c[1]); v * v }).collect();
        let rngs = (0..n.div_ceil(BLOCCO)).map(|i| Rng::new(1000 + i as u64)).collect();
        let storia = cammina.iter().flat_map(|c| [[schermo(c[0]).to_array(); S_MAX], [schermo(c[1]).to_array(); S_MAX]]).collect();
        let mut rng = Rng::new(77);
        let foto = (0..N_FOTO_C).map(|j| Foto { k, i: (rng.next_u64() % n_foto as u64) as usize, rot: 0.0, eta: (j as u32) % VITA_C }).collect();
        Mondo { dati, cammina, p2, rngs, storia, testa: 0, pieni: 1, foto, rng, parte: Parte::A, f: 0, media_g: 1.0, vista: VistaWorker::default(), ritmo_vero: RITMO }
    }

    pub fn metti_parte(&mut self, p: Parte) { self.parte = p; self.f = 0; self.pieni = 1; }

    /// cos(ωt) of the packet at frame f of part C.
    pub fn cos_wt(&self) -> f64 { (std::f64::consts::TAU * self.f as f64 / (RITMO * PERIODO_C_S)).cos() }

    /// One frame of the sampler. A and B: Metropolis moves on |Ψ|²; C: the photographs age, new ones drawn with the packet's weights.
    pub fn passo(&mut self) {
        self.f += 1;
        if self.parte == Parte::C {
            let (w, _) = self.dati.pesi(self.cos_wt());
            for j in 0..self.foto.len() {
                self.foto[j].eta += 1;
                if self.foto[j].eta >= VITA_C {
                    let mut u = self.rng.uni();
                    let mut k = w.len() - 1;
                    for (q, x) in w.iter().enumerate() { if u < *x { k = q; break; } u -= x; }
                    let i = (self.rng.next_u64() % self.dati.n_foto as u64) as usize;
                    self.foto[j] = Foto { k, i, rot: std::f64::consts::TAU * self.rng.uni(), eta: 0 };
                }
            }
            return;
        }
        let dati = self.dati.clone();
        let s = &dati.stati[dati.i_rho];
        self.cammina.par_chunks_mut(BLOCCO).zip(self.p2.par_chunks_mut(BLOCCO)).zip(self.rngs.par_iter_mut()).for_each(|((cs, ps), rng)| {
            for (c, p) in cs.iter_mut().zip(ps.iter_mut()) {
                for _ in 0..MOSSE {
                    let mut n = *c;
                    for e in n.iter_mut() { for x in e.iter_mut() { *x += PASSO * rng.gauss(); } }
                    let v = s.psi(&n[0], &n[1]);
                    let q = v * v;
                    if rng.uni() * *p < q { *c = n; *p = q; }
                }
            }
        });
        self.testa = (self.testa + 1) % S_MAX;
        self.pieni = (self.pieni + 1).min(S_MAX);
        // the tail follows a smoothed path (drawing choice, declared): the sampler's jagged walk would read as a wire
        let (t, u) = (self.testa, (self.testa + S_MAX - 1) % S_MAX);
        let liscia = |prima: [f32; 3], ora: [f32; 3]| [prima[0] + LISCIA * (ora[0] - prima[0]), prima[1] + LISCIA * (ora[1] - prima[1]), prima[2] + LISCIA * (ora[2] - prima[2])];
        for (i, c) in self.cammina.iter().enumerate() {
            self.storia[2 * i][t] = liscia(self.storia[2 * i][u], schermo(c[0]).to_array());
            self.storia[2 * i + 1][t] = liscia(self.storia[2 * i + 1][u], schermo(c[1]).to_array());
        }
    }

    /// The glows of this frame. Electrons first (their indices are the level-of-detail keys), then the nuclei.
    pub fn disegna(&mut self) -> Livello {
        let mut l = Livello::new("calcolata");
        let d = self.dati.clone();
        // the phase of a nodeless ground state is one, it turns at its energy's frequency: one turn in 8 s here (slowed, declared)
        let ph = std::f32::consts::TAU * (self.f as f32 / RITMO as f32) / 8.0 + 3.6;
        let base = tinta(ph);
        let alfa = |r: f64| (0.03 + 0.30 * (r / d.rho_max).powf(0.35)) as f32;
        let bianco = [1.0, 0.96, 0.9];
        let mut nuclei: Vec<[f64; 3]> = match d.mostra {
            Mostra::Elio => vec![[0.0; 3]],
            Mostra::H2 => vec![[0.0, 0.0, 0.7], [0.0, 0.0, -0.7]],
        };
        match self.parte {
            Parte::A | Parte::B => {
                let b = if self.parte == Parte::B {
                    let u = (self.f % GIRO_B[(d.mostra == Mostra::H2) as usize]) as f64 / GIRO_B[(d.mostra == Mostra::H2) as usize] as f64;
                    Some(riferimento(d.mostra, u))
                } else { None };
                let rho_b = b.map(|r| d.rho1(&r)).unwrap_or(1.0);
                let scala_g = if d.mostra == Mostra::H2 { 0.3 } else { 0.15 };
                let punti: Vec<&[f64; 3]> = self.cammina.iter().flat_map(|c| [&c[0], &c[1]]).collect();
                let luci: Vec<(f64, f64)> = punti.par_iter().map(|p| (d.rho1(p), b.map(|r| d.g(&r, rho_b, p)).unwrap_or(1.0))).collect();
                self.media_g = luci.iter().map(|x| x.1).sum::<f64>() / luci.len().max(1) as f64;
                let mezzo = [0.5 * (STOPS[1][0] + STOPS[4][0]), 0.5 * (STOPS[1][1] + STOPS[4][1]), 0.5 * (STOPS[1][2] + STOPS[4][2])];
                let n = punti.len();
                let f_lod = self.vista.lod_f.clamp(0.0, 1.0);
                let k = ((self.vista.scia_s * self.ritmo_vero).round() as usize).min(S_MAX - 1).min(self.pieni.saturating_sub(1));
                let f_scia = f_lod * (MAX_SCIA / (((f_lod as f64) * n as f64).max(1.0) * k.max(1) as f64)).min(1.0) as f32;
                let mut passati = [[0.0f32; 3]; S_MAX];
                for (i, (p, (r, g))) in punti.iter().zip(&luci).enumerate() {
                    let (col, a) = match b {
                        None => (base, alfa(*r)),
                        Some(_) => {
                            let t = (((*g - 1.0) / scala_g).clamp(-1.0, 1.0)) as f32;
                            let c = if t > 0.0 { STOPS[4] } else { STOPS[1] };
                            let u = t.abs();
                            ([mezzo[0] + (c[0] - mezzo[0]) * u, mezzo[1] + (c[1] - mezzo[1]) * u, mezzo[2] + (c[2] - mezzo[2]) * u], alfa(*r) * g.clamp(0.0, 3.0) as f32)
                        }
                    };
                    l.bagliore(schermo(**p).to_array(), 1.0, col, a);
                    if k > 0 {
                        let w = grafica::peso(i as u32, f_scia);
                        if w > 0.0 {
                            for (j, q) in passati.iter_mut().take(k + 1).enumerate() { *q = self.storia[i][(self.testa + S_MAX - j) % S_MAX]; }
                            passati[0] = schermo(**p).to_array();
                            scia_in(&mut l.scie, &mut l.scie_allunga, &passati[..k + 1], [col[0], col[1], col[2], a], w);
                        }
                    }
                }
                if let Some(r) = b { l.bagliore(schermo(r).to_array(), 7.0, [1.0, 0.2, 0.6], 8.0); }
            }
            Parte::C => {
                let (_, rm) = d.pesi(self.cos_wt());
                nuclei = vec![[0.0, 0.0, 0.5 * rm], [0.0, 0.0, -0.5 * rm]];
                for s in &self.foto {
                    let w = (std::f32::consts::PI * (s.eta as f32 + 0.5) / VITA_C as f32).sin();
                    let c = d.fotografia(s.k, s.i);
                    let (sn, cs) = s.rot.sin_cos();
                    for e in &c {
                        let p = [cs * e[0] - sn * e[1], sn * e[0] + cs * e[1], e[2]];
                        l.bagliore(schermo(p).to_array(), 1.0, base, alfa(d.rho1(&p)) * w);
                    }
                }
            }
        }
        let n_e = l.teste.len();
        for p in &nuclei { l.bagliore(schermo(*p).to_array(), 4.0, bianco, 2.5); }
        // the disc that holds most of the cloud
        l.chiudi(match d.mostra { Mostra::Elio => 1.6, Mostra::H2 => 2.2 });
        if let Some(g) = &mut l.gruppo { g.n = n_e; }
        l
    }

    pub fn righe(&self) -> (String, String) {
        let d = &self.dati;
        let p = d.mostra.parti().iter().find(|x| x.0 == self.parte).map(|x| x.1).unwrap_or("");
        let (nome, en) = match d.mostra { Mostra::Elio => ("ELIO", d.energia), Mostra::H2 => ("H₂", d.energia) };
        let a = format!("{nome} · {p} · E = {en:.6} Eh · stato ECG sulla regola a molti elettroni (risolutore, non la zig-zag che cammina)");
        let b = match self.parte {
            Parte::A => "ogni pallina è una posizione possibile di un elettrone: il moto è il campionatore, non una traiettoria".to_string(),
            Parte::B => format!("elettrone rosa su un giro scelto da noi; luce pesata con g: oro più del solito, blu meno · media di g {:.3}", self.media_g),
            Parte::C => { let (_, r) = d.pesi(self.cos_wt()); format!("(v₀+v₁)/√2: R = {r:.3} a₀, periodo vero {:.3} fs, qui {PERIODO_C_S} s · fotografie indipendenti, pesate col pacchetto (Born–Oppenheimer)", d.periodo_fs) }
        };
        (a, b)
    }
}

// ───────────────────────────── the worker ─────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Comando { Avanti, Pausa }

#[derive(Clone)]
pub struct Uscita {
    pub palline: Arc<Vec<Ball>>,
    pub scie: Arc<Vec<Ball>>,
    pub scie_allunga: Arc<Vec<[f32; 3]>>,
    pub n_testa: usize,
    pub raggio: f64,
    pub scia: f64,
    pub versione: u64,
    pub righe: (String, String),
    pub ms: f64,
    pub errore: Option<String>,
}

pub struct Condiviso {
    pub comando: Comando,
    pub scegli: Option<(Mostra, Parte)>,
    pub vista: VistaWorker,
    pub fine: bool,
    pub out: Uscita,
}

pub fn avvia() -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, scegli: None, vista: VistaWorker::default(), fine: false,
        out: Uscita { palline: Arc::new(vec![]), scie: Arc::new(vec![]), scie_allunga: Arc::new(vec![]), n_testa: 0, raggio: 2.0, scia: 0.0, versione: 0,
                      righe: ("Scene calcolate: scegli cosa vedere".into(), String::new()), ms: 0.0, errore: None },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let mut mondo: Option<Mondo> = None;
        let mut versione = 0u64;
        let mut ms = 0.0;
        let mut ultimo = Instant::now();
        loop {
            let t0 = Instant::now();
            let (cmd, scegli, vista, fine) = { let mut g = cc.lock().unwrap(); (g.comando, g.scegli.take(), g.vista, g.fine) };
            if fine { break; }
            if let Some((m, p)) = scegli {
                if mondo.as_ref().map(|x| x.dati.mostra) != Some(m) {
                    match Dati::carica(m) {
                        Ok(d) => mondo = Some(Mondo::new(Arc::new(d), N_CAMMINATORI)),
                        Err(e) => {
                            mondo = None;
                            let mut g = cc.lock().unwrap();
                            g.out.errore = Some(e.clone());
                            g.out.righe = (format!("{}: dati non trovati, niente disegno finto", m.nome()), e);
                        }
                    }
                }
                if let Some(x) = &mut mondo { x.metti_parte(p); }
            }
            if let Some(x) = &mut mondo {
                let dt = t0.duration_since(ultimo).as_secs_f64();
                ultimo = t0;
                if dt > 0.0 { x.ritmo_vero += 0.1 * (1.0 / dt - x.ritmo_vero); }
                x.vista = vista;
                if cmd == Comando::Avanti { x.passo(); }
                let l = x.disegna();
                let righe = x.righe();
                let gr = l.gruppo.unwrap_or(grafica::Gruppo { n: 1, raggio: 2.0, scia: 0.0 });
                versione += 1;
                ms += 0.15 * (t0.elapsed().as_secs_f64() * 1e3 - ms);
                let mut g = cc.lock().unwrap();
                g.out = Uscita { n_testa: gr.n, raggio: gr.raggio, scia: gr.scia, palline: Arc::new(l.teste), scie: Arc::new(l.scie), scie_allunga: Arc::new(l.scie_allunga), versione, righe, ms, errore: None };
            }
            let want = Duration::from_secs_f64(1.0 / RITMO);
            let el = t0.elapsed();
            if el < want { std::thread::sleep(want - el); }
        }
    });
    c
}

#[cfg(test)]
mod esami {
    //! Exams S0, S3-S6 of lab/out/sc_prereg.json: `cargo test --bin atomo calcolate -- --ignored --nocapture`.
    use super::*;

    fn dati(m: Mostra) -> Arc<Dati> { Arc::new(Dati::carica(m).expect("mancano i dati: python lab/sc_esporta_stati.py")) }

    #[test]
    fn s0_psi() {
        for m in [Mostra::Elio, Mostra::H2] {
            let d = dati(m);
            let dir = crate::esportazioni::cartella_assets().join("scene_calcolate").join(m.nome());
            let j: Value = serde_json::from_slice(&std::fs::read(dir.join("prova_psi.json")).unwrap()).unwrap();
            let x: Vec<f64> = j["X"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
            let psi: Vec<f64> = j["psi"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
            let mut e: f64 = 0.0;
            for (i, v) in psi.iter().enumerate() {
                let q = &x[6 * i..6 * i + 6];
                let r = d.stati[d.i_rho].psi(&[q[0], q[1], q[2]], &[q[3], q[4], q[5]]);
                e = e.max(((r - v) / v.abs().max(1e-300)).abs());
            }
            println!("S0 {}: errore relativo massimo {e:.2e}", m.nome());
            assert!(e <= 1e-10);
        }
    }

    /// S3 and S4: the walkers keep the state (snapshots every 100 frames from 2000 to 6000, 4000 walkers).
    #[test]
    #[ignore]
    fn s3_s4_camminatori() {
        for m in [Mostra::H2, Mostra::Elio] {
            let d = dati(m);
            let mut w = Mondo::new(d.clone(), 4000);
            for _ in 0..2000 { w.passo(); }
            let (mut opp, mut altro, mut cond, mut az, mut r, mut r12, mut vicini, mut rim, mut n) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            let dist = |a: &[f64; 3], b: &[f64; 3]| ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
            let o = [0.0; 3];
            for _ in 0..40 {
                for _ in 0..100 { w.passo(); }
                let c = &w.cammina;
                for (i, q) in c.iter().enumerate() {
                    let (z1, z2) = (q[0][2], q[1][2]);
                    if z1 * z2 < 0.0 { opp += 1.0; }
                    for (a, b) in [(z1, z2), (z2, z1)] { if a < -0.35 { cond += 1.0; if b > 0.0 { altro += 1.0; } } }
                    az += 0.5 * (z1.abs() + z2.abs());
                    r += 0.5 * (dist(&q[0], &o) + dist(&q[1], &o));
                    let d12 = dist(&q[0], &q[1]);
                    r12 += d12;
                    if d12 < 0.5 { vicini += 1.0; }
                    if dist(&q[0], &c[(i + 1) % c.len()][1]) < 0.5 { rim += 1.0; }
                    n += 1.0;
                }
            }
            match m {
                Mostra::H2 => {
                    let (po, pa, zz) = (opp / n, altro / cond, az / n);
                    let (rpo, rpa, rz) = (d.rif["P_opposti"].as_f64().unwrap(), d.rif["P_altro"].as_f64().unwrap(), d.abs_z[d.i_rho]);
                    println!("S3 P_opposti {po:.4} (calcolo {rpo:.4}), P_altro {pa:.4} (calcolo {rpa:.4}), <|z|> {zz:.4} (calcolo {rz:.4})");
                    assert!((po - rpo).abs() <= 0.01 && (pa - rpa).abs() <= 0.02 && (zz - rz).abs() <= 0.01);
                }
                Mostra::Elio => {
                    let (rr, rd, q) = (r / n, r12 / n, vicini / rim);
                    println!("S4 <r> {rr:.5} (0,929472), <r12> {rd:.5} (1,42207), P(r12<0,5) vere/rimescolate {q:.3}");
                    assert!((rr / 0.929472 - 1.0).abs() <= 0.005 && (rd / 1.42207 - 1.0).abs() <= 0.005 && q <= 0.85);
                }
            }
        }
    }

    /// S5: the mean of g over the balls in every frame of a turn of part B.
    #[test]
    #[ignore]
    fn s5_media_g() {
        for (m, giro, soglia) in [(Mostra::Elio, GIRO_B[0], 0.03), (Mostra::H2, GIRO_B[1], 0.05)] {
            let mut w = Mondo::new(dati(m), N_CAMMINATORI);
            w.metti_parte(Parte::B);
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for _ in 0..giro { w.passo(); w.disegna(); lo = lo.min(w.media_g); hi = hi.max(w.media_g); }
            println!("S5 {}: media di g fra {lo:.4} e {hi:.4}", m.nome());
            assert!(lo >= 1.0 - soglia && hi <= 1.0 + soglia);
        }
    }

    /// S6: the period of R(t) on the scene's clock, and the cloud following the protons.
    #[test]
    #[ignore]
    fn s6_vibrazione() {
        let d = dati(Mostra::H2);
        let mut w = Mondo::new(d.clone(), 100);
        w.metti_parte(Parte::C);
        let giro = (RITMO * PERIODO_C_S) as u64;
        // the frames where R is largest and smallest
        let mut rs = vec![];
        for _ in 0..giro { w.passo(); rs.push(d.pesi(w.cos_wt()).1); }
        let imax = rs.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).unwrap().0 as u64 + 1;
        let imin = rs.iter().enumerate().min_by(|a, b| a.1.total_cmp(b.1)).unwrap().0 as u64 + 1;
        // one turn of the scene is one period of cos(ωt), by construction; the real period is read from the levels
        println!("S6 periodo {:.5} fs (h/(E1-E0) = 8,00485), R fra {:.3} e {:.3}", d.periodo_fs, rs[imin as usize - 1], rs[imax as usize - 1]);
        assert!((d.periodo_fs / 8.00485 - 1.0).abs() <= 0.01);
        for (nome, ff) in [("R massimo", imax), ("R minimo", imin)] {
            let (mut s, mut n) = (0.0, 0.0);
            for gg in 0..60u64 {
                w.f = (ff + giro * gg).saturating_sub(VITA_C as u64 + 1);
                for _ in 0..=VITA_C { w.passo(); }
                assert_eq!(w.f, ff + giro * gg);
                for f in &w.foto { let c = d.fotografia(f.k, f.i); s += 0.5 * (c[0][2].abs() + c[1][2].abs()); n += 1.0; }
            }
            let (pw, _) = d.pesi((std::f64::consts::TAU * ff as f64 / giro as f64).cos());
            let calc: f64 = pw.iter().zip(&d.abs_z).map(|(a, b)| a * b).sum();
            println!("S6 {nome}: <|z|> palline {:.4}, calcolo {calc:.4}", s / n);
            assert!((s / n - calc).abs() <= 0.03);
        }
    }
}

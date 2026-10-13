//! The live zig-zag atom on a worker thread (the 3-D walk costs ~30 ms per tick on this PC; the window keeps
//! drawing at full speed while the rule runs at its own pace).
//!
//! What moves (every ball, every tick):
//!   • the cloud: ~20 000 SURFERS, each carried by the current of the walk, v = j/ρ (≤ c): around the spin
//!     axis they circle (the Dirac 1s has a circulating current), after the photon they slosh;
//!     colour = the local phase of the wave, with the rest-mass rotation e^{−imc²t/ħ} removed, so a
//!     stationary state turns its colour at its Bohr frequency while its density stays still;
//!   • six highlighted surfers with short trails of glowing balls;
//!   • the photon (an incoming packet, drawn) and the light the atom sends out: dipole waves whose
//!     brightness is the atom's own oscillating dipole ⟨z⟩(t) measured from the walk, read at the retarded
//!     time t − r/c.

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::rng::Rng;
use crate::misura::{Opzioni, Rivelatore};
use crate::render::{schermo, Ball};
use crate::accoppiato::{self as acc, Esperimento};
use crate::zz::{self, Atomo, Impulso, A};

pub const N_SURF: usize = 20_000;
pub const N_SCIA: usize = 6;
pub const L_SCIA: usize = 40;
pub const N_ONDA_OUT: usize = 7000;
/// ticks of dipole history (sets how far the emitted light is drawn: c·A per tick)
pub const STORIA_DIP: usize = 120;
/// Pulse: duration and field strength (≈ a half transfer 1s → 2p at resonance; calibrated on this lattice)
pub const SIGMA_IMP: f64 = 30.0;
pub const E0_IMP: f64 = 0.010;
/// Pulse of the coupled light: e·E at its peak (weaker: the pulse is a real wave, ~1 λ long)
pub const E0_ACC: f64 = 0.004;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Comando { Avanti, Pausa, Indietro }

#[derive(Clone)]
pub struct Uscita {
    pub palline: Arc<Vec<Ball>>,
    /// balls drawn with alpha blending (the orbital shells)
    pub onde: Arc<Vec<Ball>>,
    pub versione: u64,
    pub ticks: u64,
    pub p1s: f64,
    pub dipolo: f64,
    pub campo: f64,
    pub ms_tick: f64,
    pub omega_res: f64,
    pub legame_1s: f64,
    pub punto: [f64; 3],
    pub pronto: bool,
    /// clicks of the detector: how many, and the last one (λ_C)
    pub n_clic: usize,
    pub ultimo_clic: Option<[f64; 3]>,
    /// coupled light: on, still being prepared, populations (projections on the listened 1s and 2p½ of the 3-D
    /// rule), the 1s → 2p½ difference of the 3-D rule, the frequency of the light sent back (probe), light energy
    pub accoppiato: bool,
    pub preparando: bool,
    pub p_exc: f64,
    pub p2p: f64,
    pub omega21: f64,
    pub omega_emessa: Option<f64>,
    pub e_luce: f64,
}

pub struct Condiviso {
    pub comando: Comando,
    /// ticks per second asked
    pub ritmo: f64,
    /// a photon to send: its frequency (mc²/ħ)
    pub spara: Option<f64>,
    pub ricomincia: bool,
    /// light and atom coupled (the light box of accoppiato.rs) instead of the prescribed pulse
    pub accoppiato: bool,
    pub opz: Opzioni,
    /// a measurement asked (button MISURA) / clear the clicks
    pub misura: bool,
    pub azzera: bool,
    pub fine: bool,
    pub out: Uscita,
}

struct Mondo {
    at: Atomo,
    surf: Vec<[f64; 3]>,
    scie: Vec<Vec<[f64; 3]>>,
    dip: Vec<f64>,
    rng: Rng,
    /// frequency and start tick of the photon packets (for drawing)
    fotoni: Vec<(f64, f64)>,
    dirs: Vec<[f64; 3]>,
    eta: Vec<f64>,
    /// the coupled light + atom (when on, its atom is the one shown)
    esp: Option<Esperimento>,
}

/// The atom shown: the coupled one when the light is coupled.
fn atomo(m: &Mondo) -> &Atomo { match &m.esp { Some(e) => &e.at, None => &m.at } }

fn nuovo_mondo() -> Mondo {
    // one colour of the cube by default (docs/SOTTORETICOLO.md); ATOMO_RETICOLO=pieno for the whole cube
    let at = Atomo::new_in(zz::colore_lab());
    let mut rng = Rng::new(2026);
    let surf = at.nascite(N_SURF, &mut rng);
    let scie = surf.iter().take(N_SCIA).map(|p| vec![*p]).collect();
    // directions of the outgoing light: dipole pattern sin²θ about z; random ages along the retarded time
    let mut dirs = Vec::with_capacity(N_ONDA_OUT);
    while dirs.len() < N_ONDA_OUT {
        let v = [rng.gauss(), rng.gauss(), rng.gauss()];
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let u = [v[0] / n, v[1] / n, v[2] / n];
        if rng.uni() < 1.0 - u[2] * u[2] { dirs.push(u); }
    }
    let eta = (0..N_ONDA_OUT).map(|_| rng.uni()).collect();
    Mondo { at, surf, scie, dip: vec![0.0; STORIA_DIP], rng, fotoni: vec![], dirs, eta, esp: None }
}

/// Hue for a phase (cyclic palette of the photon laboratory).
fn tinta(angle: f32, a: f32) -> [f32; 4] {
    let t = angle.rem_euclid(std::f32::consts::TAU) / std::f32::consts::TAU;
    const STOPS: [[f32; 3]; 5] = [[0.20, 0.95, 1.00], [0.25, 0.45, 1.00], [0.70, 0.30, 1.00], [1.00, 0.35, 0.70], [1.00, 0.85, 0.30]];
    let x = t * 5.0;
    let i = (x.floor() as usize) % 5;
    let s = x - x.floor();
    let (p, q) = (STOPS[i], STOPS[(i + 1) % 5]);
    [p[0] + (q[0] - p[0]) * s, p[1] + (q[1] - p[1]) * s, p[2] + (q[2] - p[2]) * s, a]
}

fn palline(m: &Mondo, surfisti: bool) -> Vec<Ball> {
    use rayon::prelude::*;
    let at = atomo(m);
    let t = at.t();
    let rho_max = at.psi.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).reduce(|| 0.0, f64::max).max(1e-300);
    let mut out: Vec<Ball> = m.surf.par_iter().map(|p| {
        // the rest-mass phase removed: what is left turns at the binding (Bohr) frequency
        let s = at.campiona(*p);
        let rho: f64 = s.iter().map(|c| c.norm2()).sum();
        let mut best = 0;
        for q in 1..4 { if s[q].norm2() > s[best].norm2() { best = q; } }
        let ph = s[best].arg() as f32 + t as f32;
        // brightness from the local density: balls that wandered where the wave is thin are faint
        let a = (0.03 + 0.30 * (rho / rho_max).powf(0.35)) as f32;
        (schermo(*p).to_array(), 1.0, tinta(ph, a))
    }).collect();
    for s in m.scie.iter().filter(|_| surfisti) {
        for (j, p) in s.iter().enumerate() {
            let a = (j + 1) as f32 / s.len() as f32;
            out.push((schermo(*p).to_array(), 0.5 + 0.9 * a, [1.0, 0.97, 0.9, 0.15 + 0.85 * a * a]));
        }
    }
    // incoming photon packets: a drawn packet (phase colours) arriving along +x at c, centred on the atom at
    // the peak of its pulse. Its true wavelength (≈ 55 λ_C) is drawn 6 times shorter.
    for &(w, t0) in m.fotoni.iter().filter(|_| m.esp.is_none()) {
        let x_c = (t - t0) * 1.0; // λ_C, c = 1
        if !(-45.0..30.0).contains(&x_c) { continue; }
        let lam = std::f64::consts::TAU / w / 6.0;
        let kk = std::f64::consts::TAU / lam;
        let mut rng = Rng::new((t0 * 10.0) as u64 + 3);
        for _ in 0..5000 {
            let dz = 1.2 * lam * rng.gauss();
            let (y, z) = (2.5 * lam * rng.gauss(), 2.5 * lam * rng.gauss());
            let env = (-0.5 * (dz / (1.2 * lam)).powi(2)).exp();
            let a = (0.55 * env.sqrt() * (kk * dz).cos().powi(2)) as f32;
            out.push((schermo([x_c + dz, y, z]).to_array(), 1.0, tinta((kk * dz) as f32, a)));
        }
    }
    // outgoing light: balls on shells r = c·age, brightness ∝ |d(t − r/c)| (the dipole the walk had then)
    let dmax = m.dip.iter().fold(0.0f64, |a, b| a.max(b.abs())).max(1e-9);
    let soglia = 0.02;
    if dmax > soglia && m.esp.is_none() {
        let n = m.dip.len();
        for (u, &e) in m.dirs.iter().zip(&m.eta) {
            let age = (e + at.ticks as f64 * 0.003).fract() * (n - 1) as f64;
            let i = n - 1 - age as usize;
            let d = m.dip[i];
            let r = age * A + 1.5;
            let a = (0.22 * (d.abs() / dmax).powf(1.5) * (dmax / 0.6).min(1.0)) as f32;
            if a < 0.01 { continue; }
            let col = if d > 0.0 { [1.0, 0.82, 0.35, a] } else { [0.35, 0.7, 1.0, a] };
            out.push((schermo([u[0] * r, u[1] * r, u[2] * r]).to_array(), 0.8, col));
        }
    }
    out
}

fn passo(m: &mut Mondo, verso: f64, surfisti: bool) {
    match m.esp.as_mut() {
        // coupled: forward only (the light absorbed by the sponge cannot come back)
        Some(e) => { if verso > 0.0 { e.passo(); } else { return; } }
        None => if verso > 0.0 { m.at.avanti(); } else { m.at.indietro(); },
    }
    let at = match &m.esp { Some(e) => &e.at, None => &m.at };
    // with the surfers (Bohm trajectories) the balls ride the current j/ρ; without them the balls are just
    // the wave: they stay where they are and 8 % are drawn again from |ψ|² every tick
    if surfisti { at.muovi_surfisti(&mut m.surf, verso); }
    let k = if surfisti { 40 } else { N_SURF / 12 };
    let nuovi = at.nascite(k, &mut m.rng);
    for (j, p) in nuovi.into_iter().enumerate() {
        let i = N_SCIA + (m.rng.next_u64() as usize + j) % (N_SURF - N_SCIA);
        m.surf[i] = p;
    }
    for (s, p) in m.scie.iter_mut().zip(m.surf.iter()) {
        s.push(*p);
        if s.len() > L_SCIA { s.remove(0); }
    }
    let d = atomo(m).dipolo();
    if verso > 0.0 { m.dip.push(d); m.dip.remove(0); } else { m.dip.insert(0, 0.0); m.dip.pop(); }
}

/// The live light of the coupled box (wave balls, alpha blending): every other cell where E² is high, gold where
/// E_z > 0 and blue where E_z < 0. Brightness relative to the strongest light away from the atom (the near field
/// at the atom would hide the wave that leaves it).
fn palline_luce(e: &Esperimento) -> Vec<Ball> {
    use rayon::prelude::*;
    let f = &e.luce.f;
    let n = f.n;
    let b = acc::B_LUCE;
    let h = (n / 2) as i64;
    let lontano = |i: usize, j: usize, k: usize| { let d = [i as i64 - h, j as i64 - h, k as i64 - h]; d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > 64 };
    let emax = (0..n / 2).into_par_iter().map(|a| {
        let mut m = 0.0f32;
        for bb in 0..n / 2 { for c in 0..n / 2 {
            let (i, j, k) = (2 * a, 2 * bb, 2 * c);
            if !lontano(i, j, k) { continue; }
            let p = &f.f[(i * n + j) * n + k];
            m = m.max(p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re);
        }}
        m
    }).reduce(|| 0.0, f32::max).max(1e-30);
    (0..n / 2).into_par_iter().flat_map_iter(|a| {
        let mut v = Vec::new();
        for bb in 0..n / 2 { for c in 0..n / 2 {
            let (i, j, k) = (2 * a, 2 * bb, 2 * c);
            let p = &f.f[(i * n + j) * n + k];
            let e2 = p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re;
            let x = e2 / emax;
            if x < 0.03 || !lontano(i, j, k) { continue; }
            let q = ((i * n + j) * n + k) as u64 * 3;
            let pos = [0, 1, 2].map(|c2| ([i, j, k][c2] as f64 - h as f64 + 2.0 * (crate::rng::hash01(q + c2 as u64) as f64 - 0.5)) * b);
            let col = if p[2].re > 0.0 { [1.0, 0.78, 0.32] } else { [0.35, 0.72, 1.0] };
            v.push((schermo(pos).to_array(), 4.0, [col[0], col[1], col[2], (0.10 + 0.6 * x.min(1.0).sqrt()).min(0.7)]));
        }}
        v
    }).collect()
}

/// Start the worker. The shared state is read by the window every frame.
pub fn avvia() -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 25.0, spara: None, ricomincia: false, accoppiato: true, opz: Opzioni::default(), misura: false, azzera: false, fine: false,
        out: Uscita { palline: Arc::new(vec![]), onde: Arc::new(vec![]), versione: 0, ticks: 0, p1s: 1.0, dipolo: 0.0, campo: 0.0, ms_tick: 0.0,
                      omega_res: 0.0, legame_1s: 0.0, punto: [0.0; 3], pronto: false, n_clic: 0, ultimo_clic: None,
                      accoppiato: true, preparando: true, p_exc: 0.0, p2p: 0.0, omega21: 0.0, omega_emessa: None, e_luce: 0.0 },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let l1 = zz::livello(-1, 1);
        let l2 = zz::livello(1, 2);
        let omega_res = l2.omega - l1.omega;
        let mut m = nuovo_mondo();
        let mut versione = 0u64;
        let mut ms = 0.0;
        let mut riv = Rivelatore::default();
        // the coupled experiment is prepared in the background (listening to the levels of the 3-D rule, ~30 s)
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || { let _ = tx.send(Esperimento::nuovo(600, zz::colore_lab())); });
        let mut riserva: Option<Esperimento> = None;
        // the clean 1s (an exact level of the 3-D rule, docs/BATTIMENTO_1S.md), taken from the coupled experiment when
        // it is ready: the atom without coupling restarts from it too (the radial start beats with period ~476)
        let mut uno_s: Option<Vec<[crate::cplx::C; 4]>> = None;
        let mut p2p = 0.0;
        let mut opz_prima = Opzioni::default();
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, spara, ric, fine, opz, misura, azzera, accoppiato) = {
                let mut g = cc.lock().unwrap();
                let s = (g.comando, g.ritmo, g.spara.take(), g.ricomincia, g.fine, g.opz, g.misura, g.azzera, g.accoppiato);
                g.ricomincia = false;
                g.misura = false;
                g.azzera = false;
                s
            };
            if fine { break; }
            let mut cambio = false;
            if let Ok(e) = rx.try_recv() {
                uno_s = Some(e.phi1s.clone());
                if m.esp.is_none() { m.at.imposta(e.phi1s.clone()); cambio = true; }
                riserva = Some(e);
            }
            if accoppiato && m.esp.is_none() && riserva.is_some() { m.esp = riserva.take(); cambio = true; }
            if !accoppiato && m.esp.is_some() { riserva = m.esp.take(); cambio = true; }
            if ric {
                match m.esp.as_mut() {
                    Some(e) => e.ricomincia(),
                    None => {
                        let esp = m.esp.take();
                        m = nuovo_mondo();
                        if let Some(v) = &uno_s { m.at.imposta(v.clone()); }
                        m.esp = esp;
                    }
                }
                cambio = true;
            }
            if cambio {
                let nuove = (match &m.esp { Some(e) => &e.at, None => &m.at }).nascite(N_SURF, &mut m.rng);
                m.surf = nuove;
                for (s, q) in m.scie.iter_mut().zip(m.surf.iter()) { *s = vec![*q]; }
            }
            if azzera { riv.azzera(); }
            let mut flash = riv.invecchia();
            if misura {
                // Born rule: one point drawn with |ψ|²
                let p = (match &m.esp { Some(e) => &e.at, None => &m.at }).nascite(1, &mut m.rng)[0];
                riv.registra(p);
                flash = true;
                if opz.collassa {
                    match m.esp.as_mut() { Some(e) => e.at.collassa(p), None => m.at.collassa(p) }
                    m.surf = (match &m.esp { Some(e) => &e.at, None => &m.at }).nascite(N_SURF, &mut m.rng);
                    for (s, q) in m.scie.iter_mut().zip(m.surf.iter()) { *s = vec![*q]; }
                }
            }
            if opz != opz_prima { flash = true; opz_prima = opz; }
            if let Some(w) = spara {
                match m.esp.as_mut() {
                    // a real pulse of light walkers, launched in the light box
                    Some(e) => e.spara(w, E0_ACC),
                    None => {
                        // the pulse peaks when the drawn packet reaches the atom
                        let t_arrivo = m.at.t() + 2.5 * SIGMA_IMP;
                        m.at.impulsi.push(Impulso { e0: E0_IMP, omega: w, t0: t_arrivo, sigma: SIGMA_IMP });
                        m.fotoni.push((w, t_arrivo));
                    }
                }
            }
            match cmd {
                Comando::Avanti => passo(&mut m, 1.0, opz.surfisti),
                Comando::Indietro => if atomo(&m).ticks > 0 { passo(&mut m, -1.0, opz.surfisti) },
                Comando::Pausa => {}
            }
            let moved = cmd != Comando::Pausa || spara.is_some() || ric || flash || misura || azzera || cambio;
            if moved || versione == 0 {
                let mut p = palline(&m, opz.surfisti);
                riv.disegna(&mut p, 3.0, opz.mostra_clic);
                versione += 1;
                let tick_ms = t0.elapsed().as_secs_f64() * 1e3;
                ms += 0.1 * (tick_ms - ms);
                let onde = m.esp.as_ref().map(palline_luce).unwrap_or_default();
                let at = atomo(&m);
                if let Some(e) = &m.esp { if at.ticks % 8 == 0 { p2p = e.p2p(); } }
                let (p_exc, omega21, omega_emessa, e_luce, w_res) = match &m.esp {
                    Some(e) => {
                        let da = e.t_impulso.map(|t| t + 2.0 * acc::SONDA as f64 * acc::B_LUCE + 100.0);
                        (e.serie_p.last().map(|v| v.1).unwrap_or(0.0), e.omega_21(), da.and_then(|d| acc::frequenza(&e.luce.sonda, d)), e.luce.energia(), e.omega_21())
                    }
                    None => (0.0, 0.0, None, 0.0, omega_res),
                };
                let out = Uscita {
                    palline: Arc::new(p), onde: Arc::new(onde), versione, ticks: at.ticks,
                    p1s: if at.ticks % 4 == 0 { at.p_fondamentale() } else { cc.lock().unwrap().out.p1s },
                    dipolo: *m.dip.last().unwrap(), campo: at.campo(at.t()), ms_tick: ms, omega_res: w_res,
                    legame_1s: 1.0 - l1.omega, punto: m.surf[0], pronto: true, n_clic: riv.clic.len(), ultimo_clic: riv.ultimo(),
                    accoppiato: m.esp.is_some(), preparando: accoppiato && m.esp.is_none(), p_exc, p2p, omega21, omega_emessa, e_luce,
                };
                cc.lock().unwrap().out = out;
            }
            // pace: at most `ritmo` ticks per second
            let want = Duration::from_secs_f64(1.0 / ritmo.max(0.5));
            let el = t0.elapsed();
            if el < want { std::thread::sleep(want - el); }
            if cmd == Comando::Pausa { std::thread::sleep(Duration::from_millis(15)); }
        }
    });
    c
}

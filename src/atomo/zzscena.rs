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

use crate::carbonio::Rng;
use crate::render::{schermo, Ball};
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

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Comando { Avanti, Pausa, Indietro }

#[derive(Clone)]
pub struct Uscita {
    pub palline: Arc<Vec<Ball>>,
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
}

pub struct Condiviso {
    pub comando: Comando,
    /// ticks per second asked
    pub ritmo: f64,
    /// a photon to send: its frequency (mc²/ħ)
    pub spara: Option<f64>,
    pub ricomincia: bool,
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
}

fn nuovo_mondo() -> Mondo {
    let at = Atomo::new();
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
    Mondo { at, surf, scie, dip: vec![0.0; STORIA_DIP], rng, fotoni: vec![], dirs, eta }
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

fn palline(m: &Mondo) -> Vec<Ball> {
    use rayon::prelude::*;
    let t = m.at.t();
    let rho_max = m.at.psi.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).reduce(|| 0.0, f64::max).max(1e-300);
    let mut out: Vec<Ball> = m.surf.par_iter().map(|p| {
        // the rest-mass phase removed: what is left turns at the binding (Bohr) frequency
        let s = m.at.campiona(*p);
        let rho: f64 = s.iter().map(|c| c.norm2()).sum();
        let mut best = 0;
        for q in 1..4 { if s[q].norm2() > s[best].norm2() { best = q; } }
        let ph = s[best].arg() as f32 + t as f32;
        // brightness from the local density: balls that wandered where the wave is thin are faint
        let a = (0.03 + 0.30 * (rho / rho_max).powf(0.35)) as f32;
        (schermo(*p).to_array(), 1.0, tinta(ph, a))
    }).collect();
    for s in &m.scie {
        for (j, p) in s.iter().enumerate() {
            let a = (j + 1) as f32 / s.len() as f32;
            out.push((schermo(*p).to_array(), 0.5 + 0.9 * a, [1.0, 0.97, 0.9, 0.15 + 0.85 * a * a]));
        }
    }
    // incoming photon packets: a drawn packet (phase colours) arriving along +x at c, centred on the atom at
    // the peak of its pulse. Its true wavelength (≈ 55 λ_C) is drawn 6 times shorter.
    for &(w, t0) in &m.fotoni {
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
    if dmax > soglia {
        let n = m.dip.len();
        for (u, &e) in m.dirs.iter().zip(&m.eta) {
            let age = (e + m.at.ticks as f64 * 0.003).fract() * (n - 1) as f64;
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

fn passo(m: &mut Mondo, verso: f64) {
    if verso > 0.0 { m.at.avanti(); } else { m.at.indietro(); }
    m.at.muovi_surfisti(&mut m.surf, verso);
    // a few balls are reborn with |ψ|² each tick (keeps the cloud faithful on the coarse lattice)
    let k = 40;
    let nuovi = m.at.nascite(k, &mut m.rng);
    for (j, p) in nuovi.into_iter().enumerate() {
        let i = N_SCIA + (m.rng.next_u64() as usize + j) % (N_SURF - N_SCIA);
        m.surf[i] = p;
    }
    for (s, p) in m.scie.iter_mut().zip(m.surf.iter()) {
        s.push(*p);
        if s.len() > L_SCIA { s.remove(0); }
    }
    let d = m.at.dipolo();
    if verso > 0.0 { m.dip.push(d); m.dip.remove(0); } else { m.dip.insert(0, 0.0); m.dip.pop(); }
}

/// Start the worker. The shared state is read by the window every frame.
pub fn avvia() -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 25.0, spara: None, ricomincia: false, fine: false,
        out: Uscita { palline: Arc::new(vec![]), versione: 0, ticks: 0, p1s: 1.0, dipolo: 0.0, campo: 0.0, ms_tick: 0.0,
                      omega_res: 0.0, legame_1s: 0.0, punto: [0.0; 3], pronto: false },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let l1 = zz::livello(-1, 1);
        let l2 = zz::livello(1, 2);
        let omega_res = l2.omega - l1.omega;
        let mut m = nuovo_mondo();
        let mut versione = 0u64;
        let mut ms = 0.0;
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, spara, ric, fine) = {
                let mut g = cc.lock().unwrap();
                let s = (g.comando, g.ritmo, g.spara.take(), g.ricomincia, g.fine);
                g.ricomincia = false;
                s
            };
            if fine { break; }
            if ric { m = nuovo_mondo(); }
            if let Some(w) = spara {
                // the pulse peaks when the drawn packet reaches the atom
                let t_arrivo = m.at.t() + 2.5 * SIGMA_IMP;
                m.at.impulsi.push(Impulso { e0: E0_IMP, omega: w, t0: t_arrivo, sigma: SIGMA_IMP });
                m.fotoni.push((w, t_arrivo));
            }
            match cmd {
                Comando::Avanti => passo(&mut m, 1.0),
                Comando::Indietro => if m.at.ticks > 0 { passo(&mut m, -1.0) },
                Comando::Pausa => {}
            }
            let moved = cmd != Comando::Pausa || spara.is_some() || ric;
            if moved || versione == 0 {
                let p = palline(&m);
                versione += 1;
                let tick_ms = t0.elapsed().as_secs_f64() * 1e3;
                ms += 0.1 * (tick_ms - ms);
                let out = Uscita {
                    palline: Arc::new(p), versione, ticks: m.at.ticks,
                    p1s: if m.at.ticks % 4 == 0 { m.at.p_fondamentale() } else { cc.lock().unwrap().out.p1s },
                    dipolo: *m.dip.last().unwrap(), campo: m.at.campo(m.at.t()), ms_tick: ms, omega_res,
                    legame_1s: 1.0 - l1.omega, punto: m.surf[0], pronto: true,
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

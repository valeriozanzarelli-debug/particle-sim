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
/// Pulse of the scene (docs/ATOMO.md §3.8, exam C0 of lab/out/diffusa_prereg.json): the π/2 pulse computed from the
/// rule's own pulse area at the atom and its 1s → 2p dipoles, so that the atom really ends half in the 2p and its cloud
/// swings as much as it can. Same duration as E0_ACC (σ = 1 λ: the 6.5 λ box holds no longer pulse).
pub const E0_PI2: f64 = 0.013369;
/// Balls carried by the energy current of the scattered light.
pub const N_FLUSSO: usize = 6000;
/// At most this many light glow balls per light box drawn.
pub const MAX_GRIGLIA: usize = 30_000;

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
    /// scattered view (docs/ATOMO.md §3.8): the declared amplification of the scattered light balls (None: not drawn),
    /// time since the pulse (λ_C/c), energy of the scattered light (same units as e_luce), e·E of the pulse at launch
    pub amp_diffusa: Option<f64>,
    pub t_dal_impulso: Option<f64>,
    pub e_diffusa: f64,
    pub e0_impulso: f64,
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
    /// scattered view: draw the light the atom sends out (F_sc, amplified ×A, declared) and the free pulse dimmed
    pub diffusa: bool,
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
    /// scattered view: balls carried by the energy current of the scattered light (λ_C, box coordinates),
    /// max |E| of the free pulse at launch, smoothed max |E| of the scattered light away from the atom
    flusso: Vec<[f64; 3]>,
    e_ref: f64,
    m_sc: f64,
    /// sampling table of the scattered energy density (cumulative, cells of the 64³ sub-grid), rebuilt every light tick
    cdf: Vec<(f32, u32)>,
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
    Mondo { at, surf, scie, dip: vec![0.0; STORIA_DIP], rng, fotoni: vec![], dirs, eta, esp: None, flusso: vec![], e_ref: 0.0, m_sc: 0.0, cdf: vec![] }
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
    // no separate surfer with a trail of balls (Valerio, 13/10): every ball above is a glow riding the current,
    // and the moving glows together make the shapes
    let _ = surfisti;
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
    if verso > 0.0 && m.esp.is_some() {
        if atomo(m).ticks % acc::SUB as u64 == 0 { aggiorna_diffusa(m); }
        muovi_flusso(m);
    }
}

/// The live light of the coupled box (wave balls, alpha blending): every other cell where E² is high, gold where
/// E_z > 0 and blue where E_z < 0. Brightness relative to the strongest light away from the atom (the near field
/// at the atom would hide the wave that leaves it).
///
/// SCATTERED VIEW (`diffusa`, docs/ATOMO.md §3.8): the free pulse is drawn dimmed (alpha × 0.25); the light the atom
/// sends out (F_sc, its own exact box) is drawn on its own brightness scale, i.e. amplified by A = e_ref / m_sc (declared
/// on screen), rose where E_z > 0 and green where E_z < 0; and the balls of `flusso` (carried by its energy current)
/// glow where its energy density is.
fn palline_luce(m: &Mondo, diffusa: bool) -> Vec<Ball> {
    let e = match &m.esp { Some(e) => e, None => return vec![] };
    match (&e.luce.sc, diffusa) {
        (Some(sc), true) => {
            let mut v = griglia(&e.luce.f.f, None, 0.25, [[1.0, 0.78, 0.32], [0.35, 0.72, 1.0]]);
            if m.e_ref > 0.0 && m.m_sc > 0.0 {
                v.extend(griglia(&sc.f, Some(m.m_sc as f32), 1.0, [[1.0, 0.42, 0.78], [0.35, 1.0, 0.62]]));
                v.extend(palline_flusso(m, sc));
            }
            v
        }
        (Some(sc), false) => {
            // the total light, as before the split
            let tot: Vec<[crate::fotone3d::Cf; 3]> = e.luce.f.f.iter().zip(sc.f.iter()).map(|(a, b)| [0, 1, 2].map(|c| crate::fotone3d::Cf { re: a[c].re + b[c].re, im: a[c].im + b[c].im })).collect();
            griglia(&tot, None, 1.0, [[1.0, 0.78, 0.32], [0.35, 0.72, 1.0]])
        }
        (None, _) => griglia(&e.luce.f.f, None, 1.0, [[1.0, 0.78, 0.32], [0.35, 0.72, 1.0]]),
    }
}

/// Glow balls on every other cell of a light box where E² > 3 % of `scala`² (default: the strongest E away from the
/// atom), colour by the sign of E_z, alpha × `alpha`.
fn griglia(f: &[[crate::fotone3d::Cf; 3]], scala: Option<f32>, alpha: f32, col: [[f32; 3]; 2]) -> Vec<Ball> {
    use rayon::prelude::*;
    let n = acc::N_LUCE;
    let b = acc::B_LUCE;
    let h = (n / 2) as i64;
    let lontano = |i: usize, j: usize, k: usize| { let d = [i as i64 - h, j as i64 - h, k as i64 - h]; d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > 64 };
    let emax = match scala {
        Some(s) => s * s,
        None => (0..n / 2).into_par_iter().map(|a| {
            let mut m = 0.0f32;
            for bb in 0..n / 2 { for c in 0..n / 2 {
                let (i, j, k) = (2 * a, 2 * bb, 2 * c);
                if !lontano(i, j, k) { continue; }
                let p = &f[(i * n + j) * n + k];
                m = m.max(p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re);
            }}
            m
        }).reduce(|| 0.0, f32::max),
    }.max(1e-30);
    let mut out = (0..n / 2).into_par_iter().flat_map_iter(|a| {
        let mut v = Vec::new();
        for bb in 0..n / 2 { for c in 0..n / 2 {
            let (i, j, k) = (2 * a, 2 * bb, 2 * c);
            let p = &f[(i * n + j) * n + k];
            let e2 = p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re;
            let x = e2 / emax;
            if x < 0.03 || !lontano(i, j, k) { continue; }
            let q = ((i * n + j) * n + k) as u64 * 3;
            let pos = [0, 1, 2].map(|c2| ([i, j, k][c2] as f64 - h as f64 + 2.0 * (crate::rng::hash01(q + c2 as u64) as f64 - 0.5)) * b);
            let cc = if p[2].re > 0.0 { col[0] } else { col[1] };
            v.push((schermo(pos).to_array(), 4.0, [cc[0], cc[1], cc[2], alpha * (0.10 + 0.6 * x.min(1.0).sqrt()).min(0.7)]));
        }}
        v
    }).collect::<Vec<Ball>>();
    // at most MAX_GRIGLIA balls (the laptop's Intel): above it, a fixed random subset (hash of the position)
    let tot = out.len();
    if tot > MAX_GRIGLIA {
        let f = MAX_GRIGLIA as f32 / tot as f32;
        out.retain(|b| crate::rng::hash01((b.0[0].to_bits() as u64) ^ ((b.0[1].to_bits() as u64) << 20) ^ ((b.0[2].to_bits() as u64) << 40)) < f);
    }
    out
}

/// Light cell of a point (λ_C, box coordinates centred on the atom), None outside the box.
fn cella(p: [f64; 3]) -> Option<usize> {
    let n = acc::N_LUCE as i64;
    let mut s = 0usize;
    for c in 0..3 {
        let i = (p[c] / acc::B_LUCE).round() as i64 + n / 2;
        if i < 0 || i >= n { return None; }
        s = s * acc::N_LUCE + i as usize;
    }
    Some(s)
}

/// Energy density u = ½(E² + B²) and energy velocity v = S/u = (E × B)/u (|v| ≤ c) of a light cell (F = E + iB).
fn flusso_cella(p: &[crate::fotone3d::Cf; 3]) -> (f64, [f64; 3]) {
    let e = [p[0].re as f64, p[1].re as f64, p[2].re as f64];
    let b = [p[0].im as f64, p[1].im as f64, p[2].im as f64];
    let u = 0.5 * (e.iter().map(|x| x * x).sum::<f64>() + b.iter().map(|x| x * x).sum::<f64>());
    let s = [e[1] * b[2] - e[2] * b[1], e[2] * b[0] - e[0] * b[2], e[0] * b[1] - e[1] * b[0]];
    if u <= 0.0 { return (0.0, [0.0; 3]); }
    (u, [s[0] / u, s[1] / u, s[2] / u])
}

/// The flow balls: one glow ball each, brightness from the scattered energy density where it is.
fn palline_flusso(m: &Mondo, sc: &crate::fotone3d::Fotone3D) -> Vec<Ball> {
    let umax = 0.5 * m.m_sc * m.m_sc;
    m.flusso.iter().filter_map(|p| {
        let s = cella(*p)?;
        let (u, _) = flusso_cella(&sc.f[s]);
        let x = (u / umax.max(1e-300)).min(1.0);
        if x < 0.01 { return None; }
        Some((schermo(*p).to_array(), 2.2, [1.0, 0.93, 0.97, (0.15 + 0.75 * x.sqrt()) as f32]))
    }).collect()
}

/// Every light tick: the scale of the scattered light (A = e_ref / m_sc) and the table to draw new flow balls with
/// probability ∝ u away from the atom (r > 6 light cells; the near field is huge and does not travel).
fn aggiorna_diffusa(m: &mut Mondo) {
    use rayon::prelude::*;
    let e = match &m.esp { Some(e) => e, None => return };
    let sc = match &e.luce.sc { Some(s) => s, None => return };
    if m.e_ref <= 0.0 { return; }
    let n = acc::N_LUCE;
    let h = (n / 2) as i64;
    let interno = (n / 2 - acc::SPUGNA) as i64;
    let celle: Vec<(f32, u32)> = (0..(n / 2).pow(3)).into_par_iter().filter_map(|q| {
        let (a, bb, c) = (q / (n * n / 4), (q / (n / 2)) % (n / 2), q % (n / 2));
        let d = [2 * a as i64 - h, 2 * bb as i64 - h, 2 * c as i64 - h];
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if r2 <= 36 || d.iter().any(|x| x.abs() >= interno) { return None; }
        let s = ((2 * a) * n + 2 * bb) * n + 2 * c;
        let (u, _) = flusso_cella(&sc.f[s]);
        Some((u as f32, s as u32))
    }).collect();
    let emax_far = celle.par_iter().filter(|c| { let s = c.1 as usize; let d = [(s / (n * n)) as i64 - h, ((s / n) % n) as i64 - h, (s % n) as i64 - h]; d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > 64 })
        .map(|c| { let p = &sc.f[c.1 as usize]; ((p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re) as f64).sqrt() }).reduce(|| 0.0, f64::max);
    m.m_sc = emax_far.max(0.97 * m.m_sc);
    let mut acc_u = 0.0f32;
    m.cdf = celle.into_iter().filter(|c| c.0 > 0.0).map(|c| { acc_u += c.0; (acc_u, c.1) }).collect();
}

/// Move the flow balls by one atom tick along v = S/u of the scattered light; balls in empty or absorbing cells, and
/// 0.5 % per tick at random, are drawn again with probability ∝ u.
fn muovi_flusso(m: &mut Mondo) {
    let sc = match m.esp.as_ref().and_then(|e| e.luce.sc.as_ref()) { Some(s) => s, None => return };
    let tot = match m.cdf.last() { Some(c) if c.0 > 0.0 => c.0, _ => return };
    let n = acc::N_LUCE;
    let lim = (n / 2 - acc::SPUGNA) as f64 * acc::B_LUCE;
    let umin = 1e-4 * 0.5 * m.m_sc * m.m_sc;
    if m.flusso.len() < N_FLUSSO { m.flusso.resize(N_FLUSSO, [1e9; 3]); }
    for i in 0..m.flusso.len() {
        let p = m.flusso[i];
        let vivo = match cella(p) {
            Some(s) if p.iter().all(|x| x.abs() < lim) => {
                let (u, v) = flusso_cella(&sc.f[s]);
                if u > umin { m.flusso[i] = [p[0] + v[0] * A, p[1] + v[1] * A, p[2] + v[2] * A]; true } else { false }
            }
            _ => false,
        };
        if !vivo || m.rng.uni() < 0.005 {
            let x = m.rng.uni() as f32 * tot;
            let k = m.cdf.partition_point(|c| c.0 < x).min(m.cdf.len() - 1);
            let s = m.cdf[k].1 as usize;
            let h = (n / 2) as f64;
            let c = [s / (n * n), (s / n) % n, s % n];
            m.flusso[i] = [0, 1, 2].map(|q| (c[q] as f64 - h + 2.0 * (m.rng.uni() - 0.5)) * acc::B_LUCE);
        }
    }
}

/// Start the worker. The shared state is read by the window every frame.
pub fn avvia() -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 25.0, spara: None, ricomincia: false, accoppiato: true, diffusa: true, opz: Opzioni::default(), misura: false, azzera: false, fine: false,
        out: Uscita { palline: Arc::new(vec![]), onde: Arc::new(vec![]), versione: 0, ticks: 0, p1s: 1.0, dipolo: 0.0, campo: 0.0, ms_tick: 0.0,
                      omega_res: 0.0, legame_1s: 0.0, punto: [0.0; 3], pronto: false, n_clic: 0, ultimo_clic: None,
                      accoppiato: true, preparando: true, p_exc: 0.0, p2p: 0.0, omega21: 0.0, omega_emessa: None, e_luce: 0.0,
                      amp_diffusa: None, t_dal_impulso: None, e_diffusa: 0.0, e0_impulso: E0_PI2 },
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
        std::thread::spawn(move || { let _ = tx.send(Esperimento::nuovo(600, zz::colore_lab()).con_diffusa()); });
        let mut riserva: Option<Esperimento> = None;
        // the clean 1s (an exact level of the 3-D rule, docs/BATTIMENTO_1S.md), taken from the coupled experiment when
        // it is ready: the atom without coupling restarts from it too (the radial start beats with period ~476)
        let mut uno_s: Option<Vec<[crate::cplx::C; 4]>> = None;
        let mut p2p = 0.0;
        let mut opz_prima = Opzioni::default();
        let mut in_attesa: Option<f64> = None;
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, spara, ric, fine, opz, misura, azzera, accoppiato, diffusa) = {
                let mut g = cc.lock().unwrap();
                let s = (g.comando, g.ritmo, g.spara.take(), g.ricomincia, g.fine, g.opz, g.misura, g.azzera, g.accoppiato, g.diffusa);
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
                m.e_ref = 0.0;
                m.m_sc = 0.0;
                m.flusso.clear();
                m.cdf.clear();
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
            // a pulse asked while the coupled atom is still being prepared waits for it (it would otherwise go to the
            // uncoupled atom as a prescribed field)
            let spara = match spara { Some(w) if accoppiato && m.esp.is_none() => { in_attesa = Some(w); None } Some(w) => Some(w), None => match m.esp.as_ref() { Some(e) if in_attesa.is_some() => { in_attesa = None; Some(e.omega_21()) } _ => None } };
            if let Some(w) = spara {
                match m.esp.as_mut() {
                    // a real pulse of light walkers, launched in the light box
                    Some(e) => {
                        e.spara(w, E0_PI2);
                        // the reference of the declared amplification: the strongest E of the free pulse at launch
                        m.e_ref = m.e_ref.max(acc::e_max(&e.luce.f.f, false));
                    }
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
                let onde = palline_luce(&m, diffusa);
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
                    amp_diffusa: match &m.esp { Some(e) if diffusa && e.luce.sc.is_some() && m.e_ref > 0.0 && m.m_sc > 0.0 => Some(m.e_ref / m.m_sc), _ => None },
                    t_dal_impulso: m.esp.as_ref().and_then(|e| e.t_impulso.map(|t| e.at.t() - t)),
                    e_diffusa: m.esp.as_ref().map(|e| if at.ticks % 8 == 0 { e.luce.energia_diffusa() } else { cc.lock().unwrap().out.e_diffusa }).unwrap_or(0.0),
                    e0_impulso: E0_PI2,
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

/// Check of the scattered view without a window (test_diffusa::diffusa_vista): run `ticks` atom ticks of `e` (pulse
/// already fired) with the scene's own update of the flow balls; returns (mean radial velocity of the flow balls in
/// units of c, number of flow balls drawn, number of scattered glow balls, amplification A).
#[cfg(test)]
pub fn prova_vista(e: Esperimento, ticks: usize) -> (f64, usize, usize, f64) {
    let mut m = nuovo_mondo();
    m.e_ref = acc::e_max(&e.luce.f.f, false);
    m.esp = Some(e);
    for _ in 0..ticks { passo(&mut m, 1.0, false); }
    let sc = m.esp.as_ref().unwrap().luce.sc.as_ref().unwrap();
    let mut vr = 0.0;
    let mut k = 0;
    for p in &m.flusso {
        if let Some(s) = cella(*p) {
            let (u, v) = flusso_cella(&sc.f[s]);
            if u <= 0.0 { continue; }
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(1e-9);
            vr += (v[0] * p[0] + v[1] * p[1] + v[2] * p[2]) / r;
            k += 1;
        }
    }
    let nf = palline_flusso(&m, sc).len();
    let ng = griglia(&sc.f, Some(m.m_sc as f32), 1.0, [[1.0; 3]; 2]).len();
    (vr / k.max(1) as f64, nf, ng, m.e_ref / m.m_sc.max(1e-300))
}

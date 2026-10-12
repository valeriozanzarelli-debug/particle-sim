//! ORBITALS of the one-electron ion (Nd⁵⁹⁺) from the zig-zag rule: any state n = 1..4, l = 0..3 (s, p, d, f),
//! j = l ± ½, m_j, and superpositions of two of them.
//!
//! Where the numbers come from: the RADIAL form of the same rule (src/zigzag/walk.rs `Radial`: the walk in one
//! partial wave κ, same cell a = 0.5 λ_C, same Coulomb phase). Its eigenstates (ω, g(r), f(r)) are found by
//! "listening" at a frequency (inverse iteration); the n-th state of a partial wave is the (n − l)-th most bound
//! one. The 3-D wave is ψ = (g/r Ω_{κ m}, i f/r Ω_{−κ m}) with the Dirac spin-angle functions Ω built from the
//! spherical harmonics. A superposition evolves with the rule's own frequencies, ψ(t) = Σ c_k ψ_k e^{−i(ω_k − 1)t}
//! (rest-mass phase removed): its density beats at ω₁ − ω₂.
//!
//! Why not on the live 3-D lattice: the live box (48³ cells of 0.5 λ_C, atom scene) is 24 λ_C wide; the 3d
//! of Nd⁵⁹⁺ peaks at ~20 λ_C and the 4f at ~37 λ_C. A live box for n = 4 would need ~250³ cells. So the shapes
//! come from the radial rule (exact for one electron in a central field), the live lattice keeps the 1s.
//!
//! Drawing: soft density SHELLS at three levels of |ψ|² (colour = the phase of the wave, which turns in time and
//! winds around z as e^{imφ}), and SURFERS carried by the current of the wave v = j/ρ (≤ c), a few with trails.

use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::cplx::C;
use crate::render::{phase_color, schermo, Ball};
use crate::rng::{hash01, Rng};
use crate::walk::{surfer_velocity, Radial};
use crate::zz::{dirac_a_chirale, z_alpha, A};

/// One state: principal n, orbital l, 2j and 2m_j.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Stato { pub n: u32, pub l: u32, pub j2: u32, pub m2: i32 }

impl Stato {
    pub fn kappa(&self) -> i32 { if self.j2 == 2 * self.l + 1 { -(self.l as i32 + 1) } else { self.l as i32 } }
    pub fn valido(&self) -> bool {
        self.n >= 1 && self.n <= 4 && self.l < self.n && (self.j2 == 2 * self.l + 1 || (self.l > 0 && self.j2 == 2 * self.l - 1))
            && self.m2.unsigned_abs() <= self.j2 && self.m2.rem_euclid(2) == 1
    }
    /// Make it valid (keep what can be kept).
    pub fn aggiusta(mut self) -> Self {
        self.n = self.n.clamp(1, 4);
        self.l = self.l.min(self.n - 1);
        if !(self.j2 == 2 * self.l + 1 || (self.l > 0 && self.j2 == 2 * self.l - 1)) { self.j2 = 2 * self.l + 1; }
        if self.m2.rem_euclid(2) != 1 { self.m2 = 1; }
        self.m2 = self.m2.clamp(-(self.j2 as i32), self.j2 as i32);
        self
    }
    pub fn nome(&self) -> String {
        let ll = ["s", "p", "d", "f"][self.l as usize];
        let m = if self.m2 > 0 { format!("+{}/2", self.m2) } else { format!("{}/2", self.m2) };
        format!("{}{}{}/2 m={}", self.n, ll, self.j2, m)
    }
    pub fn new(n: u32, l: u32, m2: i32) -> Self { Stato { n, l, j2: 2 * l + 1, m2 }.aggiusta() }
}

/// Exact Dirac energy (units mc²) of level n, κ for Zα (for comparison only).
pub fn dirac_esatto(n: u32, kappa: i32) -> f64 {
    let za = z_alpha();
    let k = kappa.abs() as f64;
    let d = k - (k * k - za * za).sqrt();
    1.0 / (1.0 + (za / (n as f64 - d)).powi(2)).sqrt()
}

/// A radial state of the rule: frequency ω (mc²/ħ), r·g and r·f on sites r_i = (i + ½) a.
#[derive(Clone)]
pub struct Radiale { pub omega: f64, pub gf: Vec<[f64; 2]>, pub a: f64 }

/// The n-th state of partial wave κ of the radial rule: listen at many frequencies, keep the distinct bound
/// states, sorted from the most bound; the n-th of a partial wave with orbital l is number n − l − 1.
pub fn radiale(kappa: i32, n: u32, l: u32) -> Radiale {
    let za = z_alpha();
    let rmax = 30.0 + 9.0 * (n * n) as f64 * (1.0 / za);
    let w = Radial::new(kappa, A, rmax, za);
    let mut trovati: Vec<(f64, Vec<[f64; 2]>)> = Vec::new();
    // bindings from 0.3 down to 0.002 mc², log-spaced
    for k in 0..48 {
        let b = 0.3 * (0.002f64 / 0.3).powf(k as f64 / 47.0);
        let (om, gf) = w.eigen(1.0 - b, 30);
        if om >= 1.0 || om <= 0.0 { continue; }
        if trovati.iter().all(|t| (t.0 - om).abs() > 1e-6) { trovati.push((om, gf)); }
    }
    trovati.sort_by(|a, b| a.0.total_cmp(&b.0));
    let idx = (n - l - 1) as usize;
    let (omega, gf) = trovati.get(idx).cloned().unwrap_or_else(|| trovati.last().cloned().unwrap_or((1.0, vec![[0.0, 0.0]; 2])));
    Radiale { omega, gf, a: A }
}

/// Associated Legendre P_l^m(x) for m ≥ 0 (with the Condon–Shortley phase).
fn legendre(l: i32, m: i32, x: f64) -> f64 {
    let mut pmm = 1.0;
    if m > 0 {
        let s = ((1.0 - x) * (1.0 + x)).max(0.0).sqrt();
        let mut f = 1.0;
        for _ in 0..m { pmm *= -f * s; f += 2.0; }
    }
    if l == m { return pmm; }
    let mut pmmp1 = x * (2 * m + 1) as f64 * pmm;
    if l == m + 1 { return pmmp1; }
    let mut pll = 0.0;
    for ll in (m + 2)..=l {
        pll = (x * (2 * ll - 1) as f64 * pmmp1 - (ll + m - 1) as f64 * pmm) / (ll - m) as f64;
        pmm = pmmp1;
        pmmp1 = pll;
    }
    pll
}

fn fatt(k: i32) -> f64 { (1..=k).map(|v| v as f64).product() }

/// Complex spherical harmonic Y_l^m at the direction (cos θ, φ).
pub fn ylm(l: i32, m: i32, ct: f64, phi: f64) -> C {
    if m.abs() > l || l < 0 { return C::ZERO; }
    let am = m.abs();
    let norm = ((2 * l + 1) as f64 / (4.0 * std::f64::consts::PI) * fatt(l - am) / fatt(l + am)).sqrt();
    let y = C::expi(am as f64 * phi).scale(norm * legendre(l, am, ct));
    if m >= 0 { y } else { let s = if am % 2 == 0 { 1.0 } else { -1.0 }; y.conj().scale(s) }
}

/// Dirac spin-angle function Ω_{κ m} (two components), m = m2/2.
pub fn omega_km(kappa: i32, m2: i32, ct: f64, phi: f64) -> [C; 2] {
    let (l, piu) = if kappa < 0 { (-kappa - 1, true) } else { (kappa, false) };
    let m = m2 as f64 / 2.0;
    let lf = l as f64;
    let (mu, md) = ((m2 - 1) / 2, (m2 + 1) / 2); // m − ½ and m + ½
    let d = 2.0 * lf + 1.0;
    if piu {
        [ylm(l, mu, ct, phi).scale(((lf + m + 0.5) / d).max(0.0).sqrt()), ylm(l, md, ct, phi).scale(((lf - m + 0.5) / d).max(0.0).sqrt())]
    } else {
        [ylm(l, mu, ct, phi).scale(-((lf - m + 0.5) / d).max(0.0).sqrt()), ylm(l, md, ct, phi).scale(((lf + m + 0.5) / d).max(0.0).sqrt())]
    }
}

/// The 3-D Dirac spinor (Dirac representation) of a radial state at p (λ_C).
pub fn spinore(r: &Radiale, kappa: i32, m2: i32, p: [f64; 3]) -> [C; 4] {
    let rr = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(1e-9);
    let x = rr / r.a - 0.5;
    if x >= (r.gf.len() - 1) as f64 { return [C::ZERO; 4]; }
    let i = (x.floor().max(0.0) as usize).min(r.gf.len() - 2);
    let t = (x - i as f64).clamp(0.0, 1.0);
    let g = (r.gf[i][0] * (1.0 - t) + r.gf[i + 1][0] * t) / rr;
    let f = (r.gf[i][1] * (1.0 - t) + r.gf[i + 1][1] * t) / rr;
    let ct = p[2] / rr;
    let phi = p[1].atan2(p[0]);
    let up = omega_km(kappa, m2, ct, phi);
    let lo = omega_km(-kappa, m2, ct, phi);
    [up[0].scale(g), up[1].scale(g), C::I * lo[0].scale(f), C::I * lo[1].scale(f)]
}

/// Radius (λ_C) that holds 99.5 % of the state.
pub fn raggio(r: &Radiale) -> f64 {
    let tot: f64 = r.gf.iter().map(|v| v[0] * v[0] + v[1] * v[1]).sum();
    let mut acc = 0.0;
    for (i, v) in r.gf.iter().enumerate() {
        acc += v[0] * v[0] + v[1] * v[1];
        if acc > 0.995 * tot { return (i as f64 + 0.5) * r.a; }
    }
    r.gf.len() as f64 * r.a
}

pub const G: usize = 64;
pub const N_SURF_ORB: usize = 5000;
pub const N_SCIE: usize = 6;

/// What to show: one or two states with amplitudes (normalised).
#[derive(Clone, Debug, PartialEq)]
pub struct Scelta { pub stati: Vec<(Stato, f64)> }

impl Default for Scelta { fn default() -> Self { Scelta { stati: vec![(Stato::new(2, 1, 1), 0.6f64.sqrt()), (Stato::new(3, 2, 1), 0.4f64.sqrt())] } } }

/// From text like "3d+2p" or "4f" (m = +½, j = l + ½; weights equal).
pub fn scelta_da(s: &str) -> Option<Scelta> {
    let mut stati = vec![];
    for p in s.split('+') {
        let b = p.as_bytes();
        if b.len() < 2 { return None; }
        let n = (b[0] as char).to_digit(10)?;
        let l = "spdf".find(b[1] as char)? as u32;
        stati.push((Stato::new(n, l, 1), 1.0));
    }
    if stati.is_empty() { None } else { Some(Scelta { stati }) }
}

/// The grid of the chosen states: ψ_k on G³ points over [−R, R]³.
pub struct Griglia {
    pub r: f64,
    pub psi: Vec<Vec<[C; 4]>>,
    pub omega: Vec<f64>,
    pub coef: Vec<f64>,
    pub stati: Vec<Stato>,
}

impl Griglia {
    pub fn passo(&self) -> f64 { 2.0 * self.r / G as f64 }
    pub fn x(&self, i: usize) -> f64 { -self.r + (i as f64 + 0.5) * self.passo() }

    pub fn new(s: &Scelta, cache: &mut HashMap<(i32, u32), Radiale>) -> Self {
        let mut rad = Vec::new();
        for (st, _) in &s.stati {
            let k = st.kappa();
            let r = cache.entry((k, st.n)).or_insert_with(|| radiale(k, st.n, st.l)).clone();
            rad.push(r);
        }
        let r = rad.iter().map(raggio).fold(2.0, f64::max) * 1.05;
        let norm: f64 = s.stati.iter().map(|v| v.1 * v.1).sum::<f64>().sqrt().max(1e-12);
        let mut g = Griglia { r, psi: vec![], omega: rad.iter().map(|v| v.omega).collect(), coef: s.stati.iter().map(|v| v.1 / norm).collect(), stati: s.stati.iter().map(|v| v.0).collect() };
        for ((st, _), rr) in s.stati.iter().zip(&rad) {
            let mut v: Vec<[C; 4]> = (0..G * G * G).into_par_iter().map(|q| {
                let p = [g.x(q / (G * G)), g.x((q / G) % G), g.x(q % G)];
                spinore(rr, st.kappa(), st.m2, p)
            }).collect();
            // normalise on the grid
            let nn: f64 = v.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum();
            let k = 1.0 / nn.max(1e-300).sqrt();
            v.par_iter_mut().for_each(|p| for c in p.iter_mut() { *c = c.scale(k); });
            g.psi.push(v);
        }
        g
    }

    /// ψ(t) on the grid (rest-mass phase removed).
    pub fn al_tempo(&self, t: f64) -> Vec<[C; 4]> {
        let fasi: Vec<C> = self.omega.iter().zip(&self.coef).map(|(w, c)| C::expi(-(w - 1.0) * t).scale(*c)).collect();
        (0..G * G * G).into_par_iter().map(|q| {
            let mut o = [C::ZERO; 4];
            for (k, f) in fasi.iter().enumerate() { for c in 0..4 { o[c] += self.psi[k][q][c] * *f; } }
            o
        }).collect()
    }

    /// Slowest period worth seeing (λ_C/c): the beat, or the phase turn of a single state.
    pub fn periodo(&self) -> f64 {
        if self.omega.len() >= 2 && (self.omega[0] - self.omega[1]).abs() > 1e-9 { std::f64::consts::TAU / (self.omega[0] - self.omega[1]).abs() }
        else { std::f64::consts::TAU / (1.0 - self.omega[0]).abs().max(1e-6) }
    }
}

/// ψ at a point from grid values, trilinear (0 outside).
fn campiona(g: &Griglia, v: &[[C; 4]], p: [f64; 3]) -> [C; 4] {
    let h = g.passo();
    let mut i0 = [0usize; 3];
    let mut w = [0.0; 3];
    for c in 0..3 {
        let x = (p[c] + g.r) / h - 0.5;
        if x < 0.0 || x >= (G - 1) as f64 { return [C::ZERO; 4]; }
        let fl = x.floor();
        i0[c] = fl as usize;
        w[c] = x - fl;
    }
    let mut o = [C::ZERO; 4];
    for dx in 0..2 { for dy in 0..2 { for dz in 0..2 {
        let wt = (if dx == 1 { w[0] } else { 1.0 - w[0] }) * (if dy == 1 { w[1] } else { 1.0 - w[1] }) * (if dz == 1 { w[2] } else { 1.0 - w[2] });
        let q = &v[((i0[0] + dx) * G + i0[1] + dy) * G + i0[2] + dz];
        for c in 0..4 { o[c] += q[c] * wt; }
    }}}
    o
}

/// Phase of the dominant upper (large) component: the colour.
fn fase(p: &[C; 4]) -> f64 { if p[0].norm2() >= p[1].norm2() { p[0].arg() } else { p[1].arg() } }

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Comando { Avanti, Pausa, Indietro }

#[derive(Clone)]
pub struct Uscita {
    pub palline: Arc<Vec<Ball>>,
    pub onde: Arc<Vec<Ball>>,
    pub versione: u64,
    pub t: f64,
    pub scelta: Scelta,
    /// (state, ω of the rule, exact Dirac ω)
    pub livelli: Vec<(Stato, f64, f64)>,
    pub r: f64,
    pub periodo: f64,
    pub ms: f64,
    pub pronto: bool,
}

pub struct Condiviso {
    pub comando: Comando,
    pub ritmo: f64,
    pub scegli: Option<Scelta>,
    pub fine: bool,
    pub out: Uscita,
}

struct Mondo {
    g: Griglia,
    scelta: Scelta,
    t: f64,
    dt: f64,
    surf: Vec<[f64; 3]>,
    scie: Vec<Vec<[f64; 3]>>,
    rng: Rng,
}

fn nasci(g: &Griglia, v: &[[C; 4]], k: usize, rng: &mut Rng) -> Vec<[f64; 3]> {
    let w: Vec<f64> = v.iter().map(|p| p.iter().map(|c| c.norm2()).sum()).collect();
    let mut cdf = Vec::with_capacity(w.len());
    let mut s = 0.0;
    for x in &w { s += x; cdf.push(s); }
    let h = g.passo();
    (0..k).map(|_| {
        let u = rng.uni() * s;
        let q = cdf.partition_point(|&c| c < u).min(w.len() - 1);
        [g.x(q / (G * G)) + h * (rng.uni() - 0.5), g.x((q / G) % G) + h * (rng.uni() - 0.5), g.x(q % G) + h * (rng.uni() - 0.5)]
    }).collect()
}

impl Mondo {
    fn new(s: Scelta, cache: &mut HashMap<(i32, u32), Radiale>) -> Self {
        let g = Griglia::new(&s, cache);
        // the slowest motion takes ~200 steps (≈ 7 s at 30 steps per second)
        let dt = g.periodo() / 200.0;
        let v = g.al_tempo(0.0);
        let mut rng = Rng::new(77);
        let surf = nasci(&g, &v, N_SURF_ORB, &mut rng);
        let scie = surf.iter().take(N_SCIE).map(|p| vec![*p]).collect();
        Mondo { g, scelta: s, t: 0.0, dt, surf, scie, rng }
    }

    fn passo(&mut self, verso: f64) -> Vec<[C; 4]> {
        self.t += verso * self.dt;
        let v = self.g.al_tempo(self.t);
        // surfers: v = j/ρ of the wave (≤ c), in 4 midpoint sub-steps
        let g = &self.g;
        let sub = 4;
        let h = verso * self.dt / sub as f64;
        self.surf.par_iter_mut().for_each(|p| {
            for _ in 0..sub {
                let v1 = surfer_velocity(&dirac_a_chirale(campiona(g, &v, *p)));
                let q = [p[0] + 0.5 * h * v1[0], p[1] + 0.5 * h * v1[1], p[2] + 0.5 * h * v1[2]];
                let v2 = surfer_velocity(&dirac_a_chirale(campiona(g, &v, q)));
                for c in 0..3 { p[c] += h * v2[c]; }
            }
        });
        // a few are reborn with |ψ|² every step (they would otherwise drift on the coarse grid)
        let nuovi = nasci(g, &v, 25, &mut self.rng);
        for (j, p) in nuovi.into_iter().enumerate() {
            let i = N_SCIE + (self.rng.next_u64() as usize + j) % (N_SURF_ORB - N_SCIE);
            self.surf[i] = p;
        }
        for (s, p) in self.scie.iter_mut().zip(self.surf.iter()) {
            if verso > 0.0 { s.push(*p); if s.len() > 60 { s.remove(0); } } else if s.len() > 1 { s.pop(); } else { *s = vec![*p]; }
        }
        v
    }
}

/// Shells (glow) and surfers (wave balls) of the wave `v` at time t.
fn disegna(m: &Mondo, v: &[[C; 4]]) -> (Vec<Ball>, Vec<Ball>) {
    let g = &m.g;
    let rho: Vec<f32> = v.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>() as f32).collect();
    let rmax = rho.iter().cloned().fold(0.0f32, f32::max).max(1e-30);
    // three density levels: inner bright shell, middle, faint outer
    const LIV: [(f32, f32, f32); 3] = [(0.50, 0.22, 1.25), (0.18, 0.11, 1.1), (0.05, 0.05, 1.0)];
    let h = g.passo();
    let palline: Vec<Ball> = (0..G * G * G).into_par_iter().flat_map_iter(|q| {
        let mut out = Vec::new();
        let x = rho[q] / rmax;
        for (li, (lv, alfa, sz)) in LIV.iter().enumerate() {
            let d = (x - lv).abs() / lv;
            if d > 0.22 { continue; }
            let w = 1.0 - d / 0.22;
            let col = phase_color(fase(&v[q]), alfa * w);
            for r in 0..3u64 {
                let k = (q as u64 * 9 + li as u64 * 3 + r) * 3;
                let p = [g.x(q / (G * G)) + h * (hash01(k) as f64 - 0.5), g.x((q / G) % G) + h * (hash01(k + 1) as f64 - 0.5), g.x(q % G) + h * (hash01(k + 2) as f64 - 0.5)];
                out.push((schermo(p).to_array(), *sz, col));
            }
        }
        out
    }).collect();
    let mut palline = palline;
    // surfers with trails (bright, glowing)
    for s in &m.scie {
        let k = s.len();
        for (j, p) in s.iter().enumerate() {
            let x = (j + 1) as f32 / k as f32;
            palline.push((schermo(*p).to_array(), 0.5 + 0.8 * x, [1.0, 0.96, 0.88, 0.1 + 0.7 * x * x]));
        }
        if let Some(p) = s.last() { palline.push((schermo(*p).to_array(), 2.2, [1.0, 0.97, 0.9, 1.0])); }
    }
    // the other surfers: wave balls with the phase colour where they are
    let onde: Vec<Ball> = m.surf.par_iter().skip(N_SCIE).map(|p| {
        let s = campiona(g, v, *p);
        let r: f64 = s.iter().map(|c| c.norm2()).sum();
        let a = (0.25 + 0.6 * (r as f32 / rmax).powf(0.35)).min(0.85);
        (schermo(*p).to_array(), 0.7, phase_color(fase(&s), a))
    }).collect();
    (palline, onde)
}

pub fn avvia() -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 30.0, scegli: None, fine: false,
        out: Uscita { palline: Arc::new(vec![]), onde: Arc::new(vec![]), versione: 0, t: 0.0, scelta: Scelta::default(), livelli: vec![],
                      r: 20.0, periodo: 1.0, ms: 0.0, pronto: false },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let mut cache = HashMap::new();
        let mut m = Mondo::new(Scelta::default(), &mut cache);
        let mut v = m.g.al_tempo(0.0);
        let mut versione = 0u64;
        let mut ms = 0.0;
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, scegli, fine) = {
                let mut g = cc.lock().unwrap();
                (g.comando, g.ritmo, g.scegli.take(), g.fine)
            };
            if fine { break; }
            let mut cambiato = false;
            if let Some(s) = scegli { m = Mondo::new(s, &mut cache); v = m.g.al_tempo(0.0); cambiato = true; }
            match cmd {
                Comando::Avanti => { v = m.passo(1.0); cambiato = true; }
                Comando::Indietro => { v = m.passo(-1.0); cambiato = true; }
                Comando::Pausa => {}
            }
            if cambiato || versione == 0 {
                let (p, o) = disegna(&m, &v);
                versione += 1;
                let el = t0.elapsed().as_secs_f64() * 1e3;
                ms += 0.15 * (el.min(500.0) - ms);
                let livelli = m.g.stati.iter().zip(&m.g.omega).map(|(s, w)| (*s, *w, dirac_esatto(s.n, s.kappa()))).collect();
                let out = Uscita { palline: Arc::new(p), onde: Arc::new(o), versione, t: m.t, scelta: m.scelta.clone(), livelli,
                                   r: m.g.r, periodo: m.g.periodo(), ms, pronto: true };
                cc.lock().unwrap().out = out;
            }
            let want = Duration::from_secs_f64(1.0 / ritmo.max(0.5));
            let el = t0.elapsed();
            if el < want { std::thread::sleep(want - el); }
            if cmd == Comando::Pausa { std::thread::sleep(Duration::from_millis(15)); }
        }
    });
    c
}

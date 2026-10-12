//! The FREE PHOTON scene: one photon of zig-zag spin-1 balls flying through an empty 3-D box (fotone3d.rs),
//! on a worker thread (the window keeps drawing while the rule runs at its own pace).
//!
//! SCALE. Without mass the rule has no length of its own: a photon of 1 km and one of 1 pm are the same
//! pattern of balls, only the cell is bigger or smaller. So the box is fixed in WAVELENGTHS (the "quality" =
//! cells per λ) and the wavelength chosen in the panel (radio → gamma) only changes the size of a cell, the
//! scale bar and the clock. The shape is self-similar; what changes physically with λ is only how light
//! couples to matter (not in this scene).
//!
//! Layers, all drawn from the same walker field F = E + icB (one ball per cell, or per block of cells at
//! high quality, at a fixed random point inside it — no lattice pattern on screen):
//!   • E cloud: brightness E² (the crests), colour = direction of E across the flight (a colour wheel):
//!     linear light = two colours alternating crest after crest, circular light = a coloured screw;
//!   • B cloud: brightness c²B², blue;
//!   • E and B vectors: chains of balls along E (wheel colours) or B (blue-white) from the beam axis and 6
//!     lines around it — the E-vector helix of circular light, the swinging comb of linear light;
//!   • Poynting flow S = E × B (energy transport): golden chains along S on a coarse grid;
//!   • spin: chains along the local spin density E(t−1) × E(t) (green = along the flight, magenta = against);
//!   • cross-section: the plane across the flight through the centre of the pulse, drawn denser.
//!
//! Measured live from the field (not drawn): speed of the centre of energy, helicity, beam width and, once
//! the photon has flown, its spreading angle.

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::fotone3d::{base_trasversa, norma, punto, Cf, Fotone3D, Lancio, Misura, Polar};
use crate::render::{phase_color, schermo, Ball};
use crate::rng::hash01;

/// Quality: (cells per wavelength, cells per side of the box). The box is 8–12 λ wide.
pub const QUALITA: [(f64, usize, &str); 4] = [
    (8.0, 96, "veloce: 8 celle per λ, 96³"),
    (12.0, 128, "media: 12 celle per λ, 128³"),
    (16.0, 144, "alta: 16 celle per λ, 144³"),
    (24.0, 192, "massima: 24 celle per λ, 192³ (lenta)"),
];

/// Pulse length σ in wavelengths.
pub const SIGMA_LAMBDA: f64 = 1.0;

pub const DIREZIONI: [(&str, [f64; 3]); 5] = [
    ("x (di lato)", [1.0, 0.0, 0.0]),
    ("y (in profondità)", [0.0, 1.0, 0.0]),
    ("z (in alto)", [0.0, 0.0, 1.0]),
    ("obliqua (1,1,0)", [1.0, 1.0, 0.0]),
    ("diagonale (1,1,1)", [1.0, 1.0, 1.0]),
];

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Parametri {
    /// index into QUALITA
    pub qual: usize,
    /// waist w0 in wavelengths
    pub w0_lambda: f64,
    pub polar: Polar,
    pub dir: usize,
}

impl Default for Parametri {
    fn default() -> Self { Parametri { qual: 0, w0_lambda: 1.0, polar: Polar::CircolarePiu, dir: 0 } }
}

impl Parametri {
    /// Cells per wavelength.
    pub fn cpl(&self) -> f64 { QUALITA[self.qual.min(QUALITA.len() - 1)].0 }
    /// Cells per side.
    pub fn n(&self) -> usize { QUALITA[self.qual.min(QUALITA.len() - 1)].1 }
    pub fn lancio(&self) -> Lancio {
        let n = self.n();
        let lambda = self.cpl();
        let sigma = SIGMA_LAMBDA * lambda;
        let d = norma(DIREZIONI[self.dir % DIREZIONI.len()].1);
        // start far enough from the faces that the tail of the pulse does not wrap around the periodic box
        let t_fuoco = (n as f64 / 2.0 - 2.7 * sigma - 2.0).max(4.0);
        Lancio { lambda, w0: self.w0_lambda * lambda, sigma, dir: d, polar: self.polar, t_fuoco }
    }
}

/// Which layers are drawn.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vista {
    pub e_nuvola: bool,
    pub b_nuvola: bool,
    pub e_vettori: bool,
    pub b_vettori: bool,
    pub poynting: bool,
    pub spin: bool,
    pub sezione: bool,
    /// brightness of the E cloud = E along the launch polarization (as through a polarizer) instead of E²
    pub polarizzatore: bool,
}

impl Default for Vista {
    fn default() -> Self {
        Vista { e_nuvola: true, b_nuvola: false, e_vettori: true, b_vettori: false, poynting: false, spin: false, sezione: false, polarizzatore: false }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Comando { Avanti, Pausa, Indietro }

#[derive(Clone)]
pub struct Uscita {
    pub palline: Arc<Vec<Ball>>,
    pub versione: u64,
    pub ticks: u64,
    pub lanci: u64,
    pub par: Parametri,
    /// centre of energy (cells from the box centre)
    pub centro: [f64; 3],
    /// speed of the centre of energy (c), over the last ticks
    pub v: Option<f64>,
    /// the continuum prediction for this pulse (a focused pulse is slightly slower than c along its axis)
    pub v_continuo: f64,
    pub elicita: Option<f64>,
    /// beam width w = 2√⟨δ⊥²⟩ now (cells), and the spreading angle fitted on the flight so far
    pub w: f64,
    pub theta: Option<f64>,
    pub ms_tick: f64,
    pub pronto: bool,
}

pub struct Condiviso {
    pub comando: Comando,
    pub ritmo: f64,
    /// launch a new photon with these parameters
    pub lancia: Option<Parametri>,
    pub vista: Vista,
    pub fine: bool,
    pub out: Uscita,
}

struct Mondo {
    f: Fotone3D,
    par: Parametri,
    l: Lancio,
    /// measurements per tick of this flight (index = tick)
    storia: Vec<Misura>,
    /// E of the previous tick (for the helicity and the spin layer)
    e_prima: Option<(u64, Vec<[f32; 3]>)>,
    elicita: Option<f64>,
    lanci: u64,
    v_continuo: f64,
}

impl Mondo {
    /// A new flight; the launched field is kept in `cache` (building it costs 0.3–3 s, the next photon is the same).
    fn new(par: Parametri, lanci: u64, cache: &mut Option<(Lancio, Vec<[Cf; 3]>)>) -> Self {
        let n = par.n();
        let mut f = Fotone3D::new(n);
        let l = par.lancio();
        match cache {
            Some((lc, campo)) if *lc == l && campo.len() == f.f.len() => f.f.clone_from(campo),
            _ => { f.lancia(&l); *cache = Some((l, f.f.clone())); }
        }
        let m0 = f.misura(atteso(&l, 0));
        let v_continuo = Fotone3D::velocita_continuo(&l, n);
        Mondo { f, par, l, storia: vec![m0], e_prima: None, elicita: None, lanci, v_continuo }
    }

    /// Advance (or go back) one tick and measure.
    fn passo(&mut self, avanti: bool) {
        let e = if avanti { Some(self.f.campo_e()) } else { None };
        if avanti { self.f.tick(); } else { self.f.indietro(); }
        let t = self.f.ticks as usize;
        let m = self.f.misura(atteso(&self.l, self.f.ticks));
        self.storia.truncate(t);
        self.storia.push(m);
        self.e_prima = e.map(|e| (self.f.ticks, e));
        if let Some((_, e)) = &self.e_prima { self.elicita = Some(self.f.elicita(e, self.l.dir, self.l.lambda)); }
    }

    /// Speed of the centre of energy along the flight: straight-line fit over the last (up to) 16 ticks.
    fn velocita(&self) -> Option<f64> {
        let n = self.storia.len();
        if n < 6 { return None; }
        let pts: Vec<(f64, f64)> = (n.saturating_sub(16)..n).map(|t| (t as f64, punto(self.storia[t].centro, self.l.dir))).collect();
        Some(retta(&pts).1)
    }

    fn larghezza(&self, m: &Misura) -> f64 {
        let (e1, e2) = base_trasversa(self.l.dir);
        2.0 * (0.5 * (m.varianza(e1) + m.varianza(e2))).sqrt()
    }
}

/// Where the centre of the pulse should be at tick t (cells from the box centre).
pub fn atteso(l: &Lancio, t: u64) -> [f64; 3] { let s = t as f64 - l.t_fuoco; [l.dir[0] * s, l.dir[1] * s, l.dir[2] * s] }

/// Least-squares straight line: (intercept, slope).
pub fn retta(p: &[(f64, f64)]) -> (f64, f64) {
    let n = p.len() as f64;
    let (sx, sy) = p.iter().fold((0.0, 0.0), |a, q| (a.0 + q.0, a.1 + q.1));
    let (mx, my) = (sx / n, sy / n);
    let (mut sxx, mut sxy) = (0.0, 0.0);
    for q in p { sxx += (q.0 - mx).powi(2); sxy += (q.0 - mx) * (q.1 - my); }
    let b = sxy / sxx;
    (my - b * mx, b)
}

/// The spreading half-angle from the measured transverse variances of a flight: ⟨δ⊥²⟩(t) = a + b t + c t²
/// (for free light this is exact) → θ = 2√c.
pub fn theta_da(storia: &[Misura], dir: [f64; 3]) -> Option<f64> {
    if storia.len() < 20 { return None; }
    let (e1, e2) = base_trasversa(dir);
    let mut a = [[0.0f64; 3]; 3];
    let mut r = [0.0f64; 3];
    for (t, m) in storia.iter().enumerate() {
        let y = 0.5 * (m.varianza(e1) + m.varianza(e2));
        let t = t as f64;
        let b = [1.0, t, t * t];
        for i in 0..3 { r[i] += b[i] * y; for j in 0..3 { a[i][j] += b[i] * b[j]; } }
    }
    let c = risolvi3(a, r)?;
    if c[2] <= 0.0 { return None; }
    Some(2.0 * c[2].sqrt())
}

fn risolvi3(mut a: [[f64; 3]; 3], mut r: [f64; 3]) -> Option<[f64; 3]> {
    for i in 0..3 {
        let p = (i..3).max_by(|&x, &y| a[x][i].abs().partial_cmp(&a[y][i].abs()).unwrap())?;
        a.swap(i, p);
        r.swap(i, p);
        if a[i][i].abs() < 1e-300 { return None; }
        for j in i + 1..3 {
            let f = a[j][i] / a[i][i];
            for k in i..3 { a[j][k] -= f * a[i][k]; }
            r[j] -= f * r[i];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let mut s = r[i];
        for k in i + 1..3 { s -= a[i][k] * x[k]; }
        x[i] = s / a[i][i];
    }
    Some(x)
}

const BLU: [f32; 3] = [0.40, 0.62, 1.0];

#[inline]
fn d3(a: &[f32; 3], b: &[f32; 3]) -> f32 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
#[inline]
fn x3(a: &[f32; 3], b: &[f32; 3]) -> [f32; 3] { [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]] }

/// A chain of 8 balls from `p` along the unit vector `u`, length `l` (cells).
fn catena(out: &mut Vec<Ball>, p: [f64; 3], u: [f64; 3], l: f64, col: [f32; 3], alfa: f32, size: f32) {
    let k = 7;
    for q in 0..=k {
        let t = l * q as f64 / k as f64;
        let a = alfa * (0.3 + 0.7 * q as f32 / k as f32);
        out.push((schermo([p[0] + u[0] * t, p[1] + u[1] * t, p[2] + u[2] * t]).to_array(), size, [col[0], col[1], col[2], a]));
    }
}

/// The balls of the photon (scene units = cells, box centred on the origin).
pub fn palline(f: &Fotone3D, l: &Lancio, t: u64, e_prima: Option<&[[f32; 3]]>, vista: Vista) -> Vec<Ball> {
    let n = f.n;
    let dir = l.dir;
    let centro = atteso(l, t);
    let (e1d, e2d) = base_trasversa(dir);
    let f32v = |v: [f64; 3]| [v[0] as f32, v[1] as f32, v[2] as f32];
    let nhat = if l.polar == Polar::Verticale { e2d } else { e1d };
    let (e1, e2, nh, df) = (f32v(e1d), f32v(e2d), f32v(nhat), f32v(dir));
    let (imax, emax2, bmax2) = f.f.par_iter()
        .map(|p| (p[0].norm2() + p[1].norm2() + p[2].norm2(), p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re, p[0].im * p[0].im + p[1].im * p[1].im + p[2].im * p[2].im))
        .reduce(|| (0.0, 0.0, 0.0), |a, b| (a.0.max(b.0), a.1.max(b.1), a.2.max(b.2)));
    let imax = imax.max(1e-30);
    let lam = l.lambda;
    // one ball per block of `passo`³ cells: about the same number of balls at every quality
    let passo = ((lam / 8.0).round() as usize).max(1);
    let soglia = 0.006 * imax;
    let h = (n / 2) as f32;
    let nb = n / passo;
    let mut out: Vec<Ball> = Vec::new();
    // ── clouds and cross-section ──
    if vista.e_nuvola || vista.b_nuvola || vista.sezione {
        let mezzo_spessore = 0.5 * passo as f32 + 0.01;
        let cf = f32v(centro);
        let parti: Vec<Vec<Ball>> = (0..nb).into_par_iter().map(|bi| {
            let mut v = Vec::new();
            for bj in 0..nb {
                for bk in 0..nb {
                    let (i, j, k) = (bi * passo, bj * passo, bk * passo);
                    let s = (i * n + j) * n + k;
                    let p = &f.f[s];
                    let i_tot = p[0].norm2() + p[1].norm2() + p[2].norm2();
                    if i_tot < soglia { continue; }
                    let e = [p[0].re, p[1].re, p[2].re];
                    let b = [p[0].im, p[1].im, p[2].im];
                    let q = s as u64 * 3;
                    let jit = passo as f32;
                    let pos = [i as f32 - h + jit * (hash01(q) - 0.5), j as f32 - h + jit * (hash01(q + 1) - 0.5), k as f32 - h + jit * (hash01(q + 2) - 0.5)];
                    let ang = d3(&e, &e2).atan2(d3(&e, &e1)) as f64;
                    let pos64 = [pos[0] as f64, pos[1] as f64, pos[2] as f64];
                    if vista.e_nuvola {
                        let lum = if vista.polarizzatore { d3(&e, &nh).powi(2) } else { d3(&e, &e) };
                        let a = 0.008 + 0.03 * (i_tot / imax).sqrt() + 0.17 * (lum / imax).powf(0.75);
                        v.push((schermo(pos64).to_array(), 1.0, phase_color(ang, a)));
                    }
                    if vista.b_nuvola {
                        let a = 0.006 + 0.15 * (d3(&b, &b) / imax).powf(0.75);
                        // a second ball, a third of a cell away, so E and B can be told apart where they overlap
                        let pb = [pos64[0] + 0.3, pos64[1] - 0.3, pos64[2] + 0.3];
                        v.push((schermo(pb).to_array(), 1.0, [BLU[0], BLU[1], BLU[2], a]));
                    }
                    if vista.sezione {
                        let dp = [pos[0] - cf[0], pos[1] - cf[1], pos[2] - cf[2]];
                        if d3(&dp, &df).abs() < mezzo_spessore {
                            let a = 0.05 + 0.6 * (d3(&e, &e) / imax).powf(0.6);
                            v.push((schermo(pos64).to_array(), 1.6, phase_color(ang, a)));
                        }
                    }
                }
            }
            v
        }).collect();
        for p in parti { out.extend(p); }
    }
    let idx = |v: f64| ((v.round() as i64 + (n / 2) as i64).rem_euclid(n as i64)) as usize;
    let campo = |c: [f64; 3]| f.f[(idx(c[0]) * n + idx(c[1])) * n + idx(c[2])];
    // ── E and B vectors on 7 lines along the flight ──
    if vista.e_vettori || vista.b_vettori {
        let emax = (emax2.max(bmax2) as f64).sqrt().max(1e-30);
        let lung = 1.0 * l.w0;
        let mut linee = vec![[0.0f64; 3]];
        for q in 0..6 {
            let a = std::f64::consts::TAU * q as f64 / 6.0;
            let r = 0.75 * l.w0;
            linee.push([0, 1, 2].map(|c| r * (a.cos() * e1d[c] + a.sin() * e2d[c])));
        }
        for (li, off) in linee.iter().enumerate() {
            let alfa = if li == 0 { 0.9 } else { 0.55 };
            let mut s = -(n as f64) / 2.0;
            while s < n as f64 / 2.0 {
                let c = [0, 1, 2].map(|q| centro[q] + off[q] + dir[q] * s);
                let p = campo(c);
                for (vale, v, colore) in [(vista.e_vettori, [p[0].re, p[1].re, p[2].re], None), (vista.b_vettori, [p[0].im, p[1].im, p[2].im], Some(BLU))] {
                    if !vale { continue; }
                    let v = [v[0] as f64, v[1] as f64, v[2] as f64];
                    let m = punto(v, v).sqrt();
                    if m < 0.08 * emax { continue; }
                    let col = match colore {
                        Some(c) => [0.6 * c[0] + 0.4, 0.6 * c[1] + 0.4, 0.6 * c[2] + 0.4],
                        None => { let pc = phase_color(punto(v, e2d).atan2(punto(v, e1d)), 1.0); [0.7 * pc[0] + 0.3, 0.7 * pc[1] + 0.3, 0.7 * pc[2] + 0.3] }
                    };
                    catena(&mut out, c, [v[0] / m, v[1] / m, v[2] / m], lung * m / emax, col, alfa, 0.8);
                }
                s += (lam / 8.0).max(1.0);
            }
        }
    }
    // ── Poynting flow and spin density on a coarse grid inside the pulse ──
    if vista.poynting || vista.spin {
        let g = ((lam / 2.0).round() as usize).max(3);
        let ng = n / g;
        let sin_w = (std::f64::consts::TAU / lam).sin() as f32;
        let parti: Vec<Vec<Ball>> = (0..ng).into_par_iter().map(|gi| {
            let mut v = Vec::new();
            for gj in 0..ng {
                for gk in 0..ng {
                    let (i, j, k) = (gi * g, gj * g, gk * g);
                    let s = (i * n + j) * n + k;
                    let p = &f.f[s];
                    let i_tot = p[0].norm2() + p[1].norm2() + p[2].norm2();
                    if i_tot < 0.2 * imax { continue; }
                    let pos = [i as f64 - h as f64, j as f64 - h as f64, k as f64 - h as f64];
                    let e = [p[0].re, p[1].re, p[2].re];
                    let b = [p[0].im, p[1].im, p[2].im];
                    if vista.poynting {
                        let sv = x3(&e, &b);
                        let m = d3(&sv, &sv).sqrt();
                        if m > 1e-30 {
                            let u = [(sv[0] / m) as f64, (sv[1] / m) as f64, (sv[2] / m) as f64];
                            catena(&mut v, pos, u, 0.9 * lam * (m / (0.5 * imax)).min(1.5) as f64, [1.0, 0.78, 0.30], 0.5, 0.7);
                        }
                    }
                    if vista.spin {
                        if let Some(ep) = e_prima {
                            let sd = x3(&ep[s], &e);
                            let mm = d3(&sd, &sd).sqrt();
                            let m = mm / sin_w.abs().max(1e-6);
                            if m > 0.05 * imax {
                                let u = [(sd[0] / mm) as f64, (sd[1] / mm) as f64, (sd[2] / mm) as f64];
                                let col = if d3(&sd, &df) >= 0.0 { [0.45, 1.0, 0.55] } else { [1.0, 0.40, 0.85] };
                                catena(&mut v, [pos[0] + 0.5, pos[1] + 0.5, pos[2] + 0.5], u, 0.9 * lam * (m / (0.5 * imax)).min(1.5) as f64, col, 0.5, 0.7);
                            }
                        }
                    }
                }
            }
            v
        }).collect();
        for p in parti { out.extend(p); }
    }
    // ── the box: dim balls along its 12 edges (the 3-D frame while orbiting) ──
    let lo = -h - 0.5;
    let hi = h - 0.5;
    let steps = 64;
    for e in 0..12 {
        for s in 0..=steps {
            let t = lo + (hi - lo) * s as f32 / steps as f32;
            let (u, v) = ([lo, hi][e & 1], [lo, hi][(e >> 1) & 1]);
            let p = match e / 4 { 0 => [t, u, v], 1 => [u, t, v], _ => [u, v, t] };
            out.push((schermo([p[0] as f64, p[1] as f64, p[2] as f64]).to_array(), 0.7, [0.45, 0.55, 0.75, 0.10]));
        }
    }
    out
}

pub fn avvia(par: Parametri) -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 30.0, lancia: None, vista: Vista::default(), fine: false,
        out: Uscita { palline: Arc::new(vec![]), versione: 0, ticks: 0, lanci: 0, par, centro: [0.0; 3], v: None, v_continuo: 1.0,
                      elicita: None, w: 0.0, theta: None, ms_tick: 0.0, pronto: false },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let mut cache = None;
        let mut m = Mondo::new(par, 1, &mut cache);
        let mut versione = 0u64;
        let mut ms = 0.0;
        let mut vista_prima = Vista::default();
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, lancia, vista, fine) = {
                let mut g = cc.lock().unwrap();
                (g.comando, g.ritmo, g.lancia.take(), g.vista, g.fine)
            };
            if fine { break; }
            let mut cambiato = vista != vista_prima;
            vista_prima = vista;
            if let Some(p) = lancia { let l = m.lanci + 1; m = Mondo::new(p, l, &mut cache); cambiato = true; }
            match cmd {
                Comando::Avanti => {
                    // the photon has crossed the box: the next one (same parameters)
                    if m.f.ticks as f64 >= 2.0 * m.l.t_fuoco { let l = m.lanci + 1; m = Mondo::new(m.par, l, &mut cache); }
                    else { m.passo(true); }
                    cambiato = true;
                }
                Comando::Indietro => if m.f.ticks > 0 { m.passo(false); cambiato = true; },
                Comando::Pausa => {}
            }
            if cambiato || versione == 0 {
                let ep = m.e_prima.as_ref().filter(|(t, _)| *t == m.f.ticks).map(|(_, e)| e.as_slice());
                let p = palline(&m.f, &m.l, m.f.ticks, ep, vista);
                versione += 1;
                let tick_ms = t0.elapsed().as_secs_f64() * 1e3;
                ms += 0.15 * (tick_ms.min(500.0) - ms);
                let ult = *m.storia.last().unwrap();
                let out = Uscita {
                    palline: Arc::new(p), versione, ticks: m.f.ticks, lanci: m.lanci, par: m.par, centro: ult.centro,
                    v: m.velocita(), v_continuo: m.v_continuo, elicita: m.elicita, w: m.larghezza(&ult), theta: theta_da(&m.storia, m.l.dir),
                    ms_tick: ms, pronto: true,
                };
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

//! The FREE PHOTON scene: one photon of zig-zag spin-1 balls flying through an empty 3-D box (fotone3d.rs),
//! on a worker thread (the window keeps drawing while the rule runs at its own pace).
//!
//! SCALE. Without mass the rule has no length of its own: a photon of 1 km and one of 1 pm are the same
//! pattern of balls, only the cell is bigger or smaller. So the box is fixed in WAVELENGTHS (the "quality" =
//! cells per λ) and the wavelength chosen in the panel (radio → gamma) only changes the size of a cell, the
//! scale bar and the clock. The shape is self-similar; what changes physically with λ is only how light
//! couples to matter (not in this scene).
//!
//! What is drawn (tracers riding the field, wavefront sheets, the E helix, streaks...) is in fotdisegno.rs.
//!
//! Measured live from the field (not drawn): speed of the centre of energy, helicity, beam width and, once
//! the photon has flown, its spreading angle.

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::fotone3d::{base_trasversa, norma, punto, Cf, Fotone3D, Lancio, Misura, Polar};
use crate::fotdisegno::{self, Traccianti};
use crate::misura::{Opzioni, Rivelatore};
use crate::render::Ball;
use crate::rng::Rng;

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

pub use crate::fotdisegno::Vista;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Comando { Avanti, Pausa, Indietro }

#[derive(Clone)]
pub struct Uscita {
    pub palline: Arc<Vec<Ball>>,
    /// the wave balls (alpha blending, as in the photon laboratory)
    pub onde: Arc<Vec<Ball>>,
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
    pub n_clic: usize,
    pub ultimo_clic: Option<[f64; 3]>,
}

pub struct Condiviso {
    pub comando: Comando,
    pub ritmo: f64,
    /// launch a new photon with these parameters
    pub lancia: Option<Parametri>,
    pub vista: Vista,
    pub opz: Opzioni,
    pub misura: bool,
    pub azzera: bool,
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
    /// the balls that ride the field (drawing only)
    tr: Traccianti,
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
        let tr = Traccianti::new(&f, 1000 + lanci);
        Mondo { f, par, l, storia: vec![m0], e_prima: None, elicita: None, lanci, v_continuo, tr }
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
        self.tr.avanza(&self.f, if avanti { 1.0 } else { -1.0 }, atteso(&self.l, self.f.ticks));
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

pub fn avvia(par: Parametri) -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 30.0, lancia: None, vista: Vista::default(), opz: Opzioni::default(), misura: false, azzera: false, fine: false,
        out: Uscita { palline: Arc::new(vec![]), onde: Arc::new(vec![]), versione: 0, ticks: 0, lanci: 0, par, centro: [0.0; 3], v: None, v_continuo: 1.0,
                      elicita: None, w: 0.0, theta: None, ms_tick: 0.0, pronto: false, n_clic: 0, ultimo_clic: None },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let mut cache = None;
        let mut m = Mondo::new(par, 1, &mut cache);
        let mut versione = 0u64;
        let mut ms = 0.0;
        let mut vista_prima = Vista::default();
        let mut riv = Rivelatore::default();
        let mut opz_prima = Opzioni::default();
        let mut rng = Rng::new(4242);
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, lancia, vista, fine, opz, misura, azzera) = {
                let mut g = cc.lock().unwrap();
                let s = (g.comando, g.ritmo, g.lancia.take(), g.vista, g.fine, g.opz, g.misura, g.azzera);
                g.misura = false;
                g.azzera = false;
                s
            };
            if fine { break; }
            let mut cambiato = vista != vista_prima || opz != opz_prima || riv.invecchia() || azzera;
            vista_prima = vista;
            opz_prima = opz;
            if azzera { riv.azzera(); }
            if misura {
                // a detector click: one point drawn with the energy density |F|² of the photon
                let p = fotdisegno::semina(&m.f, 1, &mut rng)[0];
                riv.registra([p[0] as f64, p[1] as f64, p[2] as f64]);
                cambiato = true;
                // collapse: a photon detected is absorbed; the next one leaves
                if opz.collassa { let l = m.lanci + 1; m = Mondo::new(m.par, l, &mut cache); }
            }
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
                let (mut p, onde) = fotdisegno::palline(&m.f, &m.l, atteso(&m.l, m.f.ticks), &m.tr, ep, vista, opz.surfisti);
                riv.disegna(&mut p, 0.6 * m.l.lambda, opz.mostra_clic);
                versione += 1;
                let tick_ms = t0.elapsed().as_secs_f64() * 1e3;
                ms += 0.15 * (tick_ms.min(500.0) - ms);
                let ult = *m.storia.last().unwrap();
                let out = Uscita {
                    palline: Arc::new(p), onde: Arc::new(onde), versione, ticks: m.f.ticks, lanci: m.lanci, par: m.par, centro: ult.centro,
                    v: m.velocita(), v_continuo: m.v_continuo, elicita: m.elicita, w: m.larghezza(&ult), theta: theta_da(&m.storia, m.l.dir),
                    ms_tick: ms, pronto: true, n_clic: riv.clic.len(), ultimo_clic: riv.ultimo(),
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

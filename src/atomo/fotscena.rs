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
use crate::gpu_fotone::{self, FotoneGpu};
use crate::misura::{Opzioni, Rivelatore};
use crate::render::Ball;
use crate::rng::Rng;

/// Quality: (cells per wavelength, cells per side of the box). The box is 8–12 λ wide (≥ 10,7 λ from 12 cells per λ up to medio-alta).
pub const QUALITA: [(f64, usize, &str); 5] = [
    // floor at 12 cells per λ: below it diffraction changes (ka/tan ka, −21 % at 8), the dynamics would
    // not just blur (docs/MOTORE_OTTIMIZZATO.md, exam EL). The box must stay ≥ ~10 λ too: at 12 cells per λ
    // an 8 λ box spreads the beam −11.5 % against Maxwell, a 10.7 λ box −9.8 % (test fotone_divergenza_scatola).
    (12.0, 128, "bassa: 12 celle per λ, 128³"),
    (16.0, 160, "media: 16 celle per λ, 160³"),
    (20.0, 176, "medio-alta: 20 celle per λ, 176³"),
    (24.0, 192, "alta: 24 celle per λ, 192³"),
    (32.0, 384, "ultra: 32 celle per λ, 384³ (per una GPU da 16 GB; qui lentissima)"),
];

/// Energy balls drawn for each quality (docs/MOTORE_GRAFICO.md; each ball is drawn with a 2-ball tail).
pub const PALLINE_QUAL: [usize; 5] = [20_000, 60_000, 120_000, 250_000, 1_000_000];
/// Tail balls per energy ball (one every 2 ticks) for each quality.
pub const SCIA_QUAL: [usize; 5] = [0, 0, 0, 0, 0];

/// The named presets (`--qualita`): index into QUALITA.
pub fn preset(nome: &str) -> Option<usize> {
    match nome { "bassa" => Some(0), "media" => Some(1), "alta" => Some(3), "ultra" => Some(4), _ => None }
}

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
    /// (cells per λ, cells per side) overriding the preset: tests of the rule pin their own grid, so a change
    /// of the presets does not change what they test
    pub griglia: Option<(f64, usize)>,
}

impl Default for Parametri {
    fn default() -> Self { Parametri { qual: 0, w0_lambda: 1.0, polar: Polar::CircolarePiu, dir: 0, griglia: None } }
}

impl Parametri {
    /// Cells per wavelength.
    pub fn cpl(&self) -> f64 { if let Some((c, _)) = self.griglia { return c; } QUALITA[self.qual.min(QUALITA.len() - 1)].0 }
    /// Cells per side.
    pub fn n(&self) -> usize { if let Some((_, n)) = self.griglia { return n; } QUALITA[self.qual.min(QUALITA.len() - 1)].1 }
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
    /// stretch of each energy ball (along E), parallel to `onde`
    pub allunga: Arc<Vec<[f32; 3]>>,
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
    /// total energy Σ|F|², total Poynting flux Σ E × B and the spin read from the energy balls (ħ per photon):
    /// measured on both paths (CPU and GPU), so every preset shows the same phenomena
    pub energia: f64,
    pub poynting: [f64; 3],
    pub spin_palline: Option<f64>,
    /// where the rule runs ("gpu" or "cpu", with the reason when the GPU was refused)
    pub regola: String,
    /// GPU path: the energy balls as instances already on the GPU (drawn without a CPU round trip)
    pub onde_gpu: Option<crate::istanze::IstanzeGpu>,
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
    /// where the camera is (cells from the box centre): the energy balls are drawn from the farthest to the
    /// nearest, so the alpha blending shows the nearest layer in front (docs/BAGLIORI.md). None = unsorted.
    pub occhio: Option<[f64; 3]>,
    pub out: Uscita,
}

/// Sort the energy balls (and their stretch) from the farthest to the nearest to `occhio`.
pub fn ordina(onde: &mut Vec<Ball>, allunga: &mut Vec<[f32; 3]>, occhio: [f64; 3]) {
    let o = [occhio[0] as f32, occhio[1] as f32, occhio[2] as f32];
    // the balls are already in screen coordinates (render::schermo: (x, z, −y))
    let os = [o[0], o[2], -o[1]];
    let d2 = |b: &Ball| { let p = b.0; (p[0] - os[0]).powi(2) + (p[1] - os[1]).powi(2) + (p[2] - os[2]).powi(2) };
    let mut idx: Vec<(f32, usize)> = onde.iter().enumerate().map(|(i, b)| (-d2(b), i)).collect();
    idx.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let st = allunga.len() == onde.len();
    *onde = idx.iter().map(|&(_, i)| onde[i]).collect();
    if st { *allunga = idx.iter().map(|&(_, i)| allunga[i]).collect(); }
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
    /// alta/ultra: field, balls and measurements on the GPU (gpu_fotone.rs); None = CPU path
    gpu: Option<FotoneGpu>,
    regola: String,
    poynting: [f64; 3],
    spin: Option<f64>,
}

impl Mondo {
    /// A new flight; the launched field is kept in `cache` (building it costs 0.3–3 s, the next photon is the same).
    fn new(par: Parametri, lanci: u64, cache: &mut Option<(Lancio, Vec<[Cf; 3]>)>, riuso: Option<FotoneGpu>) -> Self {
        let n = par.n();
        let mut f = Fotone3D::new(n);
        let l = par.lancio();
        match cache {
            Some((lc, campo)) if *lc == l && campo.len() == f.f.len() => f.f.clone_from(campo),
            _ => { f.lancia(&l); *cache = Some((l, f.f.clone())); }
        }
        let m0 = f.misura(atteso(&l, 0));
        let v_continuo = Fotone3D::velocita_continuo(&l, n);
        let mut tr = Traccianti::new(&f, 1000 + lanci, PALLINE_QUAL[par.qual.min(PALLINE_QUAL.len() - 1)]);
        tr.scia_n = SCIA_QUAL[par.qual.min(SCIA_QUAL.len() - 1)];
        let (gpu, regola) = if gpu_fotone::vuole_gpu(par.qual) {
            match FotoneGpu::prova(&f, &tr, &l, riuso) {
                Ok(g) => (Some(g), "gpu".to_string()),
                Err(e) => { eprintln!("fotone: regola sulla CPU, la GPU non va: {e}"); (None, format!("cpu (gpu: {e})")) }
            }
        } else { (None, "cpu".to_string()) };
        if lanci == 1 { eprintln!("fotone: regola sulla {} ({}³, {} palline)", regola, n, tr.pos.len()); }
        let poynting = gpu_fotone::poynting_cpu(&f);
        Mondo { f, par, l, storia: vec![m0], e_prima: None, elicita: None, lanci, v_continuo, tr, gpu, regola, poynting, spin: None }
    }

    /// Advance (or go back) one tick and measure.
    fn passo(&mut self, avanti: bool) {
        if let Some(g) = &mut self.gpu {
            let t = if avanti { self.f.ticks + 1 } else { self.f.ticks.saturating_sub(1) };
            let r = g.passo(&mut self.f, &mut self.tr, avanti, atteso(&self.l, t), &self.l);
            self.storia.truncate(self.f.ticks as usize);
            self.storia.push(r.misura);
            self.e_prima = None;
            if let Some(h) = r.elicita { self.elicita = Some(h); }
            self.poynting = r.poynting;
            self.spin = Some(r.spin);
            return;
        }
        let e = if avanti { Some(self.f.campo_e()) } else { None };
        if avanti { self.f.tick(); } else { self.f.indietro(); }
        let t = self.f.ticks as usize;
        let m = self.f.misura(atteso(&self.l, self.f.ticks));
        self.storia.truncate(t);
        self.storia.push(m);
        self.e_prima = e.map(|e| (self.f.ticks, e));
        if let Some((_, e)) = &self.e_prima { self.elicita = Some(self.f.elicita(e, self.l.dir, self.l.lambda)); }
        self.tr.avanza(&self.f, if avanti { 1.0 } else { -1.0 }, atteso(&self.l, self.f.ticks));
        // the same measurements as the GPU path (gpu_fotone.rs): Poynting flux and spin from the balls
        self.poynting = gpu_fotone::poynting_cpu(&self.f);
        let rif = atteso(&self.l, self.f.ticks);
        let (lr, vs) = gpu_fotone::somme_palline_cpu(&self.tr.pos, &self.tr.vel, rif);
        self.spin = Some(gpu_fotone::spin_da_somme(lr, vs, rif, m.centro, self.l.dir, self.l.lambda, self.tr.pos.len()));
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

/// For unattended photos: the photon stops at this tick (0 = never), so two runs show the same frame.
pub static TICK_FERMO: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

pub fn avvia(par: Parametri) -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 30.0, lancia: None, vista: Vista::default(), opz: Opzioni::default(), misura: false, azzera: false, fine: false, occhio: None,
        out: Uscita { palline: Arc::new(vec![]), onde: Arc::new(vec![]), allunga: Arc::new(vec![]), versione: 0, ticks: 0, lanci: 0, par, centro: [0.0; 3], v: None, v_continuo: 1.0,
                      elicita: None, w: 0.0, theta: None, ms_tick: 0.0, energia: 0.0, poynting: [0.0; 3], spin_palline: None,
                      regola: String::new(), onde_gpu: None, pronto: false, n_clic: 0, ultimo_clic: None },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let mut cache = None;
        let mut m = Mondo::new(par, 1, &mut cache, None);
        let mut versione = 0u64;
        let mut ms = 0.0;
        let mut vista_prima = Vista::default();
        let mut riv = Rivelatore::default();
        let mut opz_prima = Opzioni::default();
        let mut rng = Rng::new(4242);
        let mut occhio_prima: Option<[f64; 3]> = None;
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, lancia, vista, fine, opz, misura, azzera, occhio) = {
                let mut g = cc.lock().unwrap();
                let s = (g.comando, g.ritmo, g.lancia.take(), g.vista, g.fine, g.opz, g.misura, g.azzera, g.occhio);
                g.misura = false;
                g.azzera = false;
                s
            };
            if fine { break; }
            let fermo = TICK_FERMO.load(std::sync::atomic::Ordering::Relaxed);
            let cmd = if fermo > 0 && m.f.ticks >= fermo && cmd == Comando::Avanti { Comando::Pausa } else { cmd };
            // the camera moved: draw again in the new depth order
            let occhio_mosso = match (occhio, occhio_prima) {
                (Some(a), Some(b)) => { let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
                                        d > 0.01 * (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt() }
                (a, b) => a.is_some() != b.is_some(),
            };
            if occhio_mosso { occhio_prima = occhio; }
            let mut cambiato = vista != vista_prima || opz != opz_prima || riv.invecchia() || azzera || occhio_mosso;
            vista_prima = vista;
            opz_prima = opz;
            if azzera { riv.azzera(); }
            if misura {
                // a detector click: one point drawn with the energy density |F|² of the photon
                if let Some(g) = &mut m.gpu { g.campo_cpu(&mut m.f); }
                let p = fotdisegno::semina(&m.f, 1, &mut rng)[0];
                riv.registra([p[0] as f64, p[1] as f64, p[2] as f64]);
                cambiato = true;
                // collapse: a photon detected is absorbed; the next one leaves
                if opz.collassa { let l = m.lanci + 1; m = Mondo::new(m.par, l, &mut cache, m.gpu.take()); }
            }
            if let Some(p) = lancia { let l = m.lanci + 1; m = Mondo::new(p, l, &mut cache, m.gpu.take()); cambiato = true; }
            match cmd {
                Comando::Avanti => {
                    // the photon has crossed the box: the next one (same parameters)
                    if m.f.ticks as f64 >= 2.0 * m.l.t_fuoco { let l = m.lanci + 1; m = Mondo::new(m.par, l, &mut cache, m.gpu.take()); }
                    else { m.passo(true); }
                    cambiato = true;
                }
                Comando::Indietro => if m.f.ticks > 0 { m.passo(false); cambiato = true; },
                Comando::Pausa => {}
            }
            if cambiato || versione == 0 {
                let ep = m.e_prima.as_ref().filter(|(t, _)| *t == m.f.ticks).map(|(_, e)| e.as_slice());
                let (mut p, onde, allunga, onde_gpu) = match &mut m.gpu {
                    Some(g) => { let c = atteso(&m.l, m.f.ticks); g.palline(&mut m.f, &m.l, c, &m.tr, vista, opz.surfisti, occhio_prima) }
                    None => {
                        let (a, mut b, mut c) = fotdisegno::palline(&m.f, &m.l, atteso(&m.l, m.f.ticks), &m.tr, ep, vista, opz.surfisti);
                        if let Some(o) = occhio_prima { ordina(&mut b, &mut c, o); }
                        (a, b, c, None)
                    }
                };
                riv.disegna(&mut p, 0.6 * m.l.lambda, opz.mostra_clic);
                versione += 1;
                let tick_ms = t0.elapsed().as_secs_f64() * 1e3;
                ms += 0.15 * (tick_ms.min(500.0) - ms);
                let ult = *m.storia.last().unwrap();
                let out = Uscita {
                    palline: Arc::new(p), onde: Arc::new(onde), allunga: Arc::new(allunga), versione, ticks: m.f.ticks, lanci: m.lanci, par: m.par, centro: ult.centro,
                    v: m.velocita(), v_continuo: m.v_continuo, elicita: m.elicita, w: m.larghezza(&ult), theta: theta_da(&m.storia, m.l.dir),
                    ms_tick: ms, energia: ult.energia, poynting: m.poynting, spin_palline: m.spin, regola: m.regola.clone(), onde_gpu,
                    pronto: true, n_clic: riv.clic.len(), ultimo_clic: riv.ultimo(),
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

//! L'ELETTRONE LIBERO A BAGLIORI (studio "bagliori", docs/BAGLIORI.md, pre-registrazione lab/out/gl_prereg.json).
//!
//! Un pacchetto gaussiano dell'elettrone (o del muone) fermo nel vuoto, con la stessa regola zig-zag di tutto il
//! resto (src/zigzag/walk.rs, `Dirac3D`: corri a c, poi ruota di θ = m c a / ħ). Niente nucleo, niente luce.
//!
//! Cosa si vede: ~40 000 BAGLIORI, ognuno un surfista che nasce con |ψ|² e poi cavalca la corrente della regola,
//! v = j/ρ (≤ c), letta dallo spinore nel punto. Colore = fase locale dell'onda (tolta la rotazione di riposo
//! e^{−imc²t/ħ}); luminosità dalla densità locale. Nessuna scia, nessuna pallina speciale, nessuna forma disegnata:
//! le forme sono quelle che fanno i bagliori muovendosi insieme.
//!
//! Il pacchetto di partenza è quello dell'esame F3 di docs/ZIGZAG.md (σ = 3 λ_C, θ = 0,25, spinore
//! (L, R) = ((0, i), (1, 0))/√2): mescola energie positive e negative, quindi TREMA (Zitterbewegung): il suo centro
//! gira su un cerchio di raggio ~ħ/2mc a frequenza ~2mc²/ħ. Lì era misurato da ⟨r⟩ dell'onda; qui si rimisura dai
//! bagliori disegnati (esami E1–E3).
//!
//! Il muone: stessa scena con la massa 206,768 volte più grande (la massa è un INGRESSO, non esce dalla regola).
//! Senza nucleo la regola non ha una lunghezza sua: il muone libero è lo stesso disegno con la cella 207 volte più
//! piccola e il tick 207 volte più corto. Cambiano solo il righello e l'orologio (come il fotone e la sua λ).

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::cplx::C;
use crate::render::{phase_color, schermo, Ball};
use crate::rng::Rng;
use crate::walk::{order, surfer_velocity, Dirac3D, HALF};

/// Lato della scatola (celle), flip θ = m c a / ħ (= passo del reticolo in λ_C), larghezza del pacchetto (λ_C).
pub const N_EL: usize = 96;
pub const THETA: f64 = 0.25;
pub const SIGMA: f64 = 3.0;
/// Bagliori disegnati.
pub const N_BAGLIORI: usize = 40_000;
/// Dopo questo tick il pacchetto si è allargato troppo per la scatola periodica: si ricomincia.
pub const TICK_MAX: u64 = 200;
/// Massa del muone in masse dell'elettrone (CODATA 2018): un ingresso.
pub const MASSA_MUONE: f64 = 206.768_283;

/// Per le foto: la scena si ferma a questo tick (0 = mai), come fotscena::TICK_FERMO.
pub static TICK_FERMO: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Particella { Elettrone, Muone }

impl Particella {
    pub fn nome(self) -> &'static str { match self { Particella::Elettrone => "elettrone", Particella::Muone => "muone" } }
    /// m / mₑ
    pub fn massa(self) -> f64 { match self { Particella::Elettrone => 1.0, Particella::Muone => MASSA_MUONE } }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Comando { Avanti, Pausa }

#[derive(Clone)]
pub struct Uscita {
    pub palline: Arc<Vec<Ball>>,
    pub versione: u64,
    pub ticks: u64,
    /// tempo in ħ/mc² della particella mostrata
    pub t: f64,
    /// centro dei bagliori e centro dell'onda ⟨r⟩ (in λ_C della particella)
    pub centro_palline: [f64; 3],
    pub centro_onda: [f64; 3],
    pub particella: Particella,
    pub ms_tick: f64,
    pub pronto: bool,
}

pub struct Condiviso {
    pub comando: Comando,
    pub ritmo: f64,
    pub ricomincia: bool,
    /// cambia particella (ricomincia)
    pub particella: Option<Particella>,
    pub fine: bool,
    pub out: Uscita,
}

/// Lo spinore iniziale dell'esame F3 (chirale: L₀, L₁, R₀, R₁).
pub fn spinore_f3() -> [C; 4] {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    [C::ZERO, C::new(0.0, s), C::real(s), C::ZERO]
}

/// Il pacchetto: gaussiana exp(−r²/4σ²) per lo spinore, normalizzato (come zz_core.gaussian3).
pub fn pacchetto(n: usize, theta: f64, sigma: f64, spin: [C; 4]) -> Dirac3D {
    let mut d = Dirac3D::new(n, theta, theta);
    let xs: Vec<f64> = (0..n).map(|i| d.x(i)).collect();
    d.psi.par_iter_mut().enumerate().for_each(|(s, p)| {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        let r2 = xs[i] * xs[i] + xs[j] * xs[j] + xs[k] * xs[k];
        let g = (-r2 / (4.0 * sigma * sigma)).exp();
        for q in 0..4 { p[q] = spin[q].scale(g); }
    });
    let nn = d.norm2().sqrt();
    d.psi.par_iter_mut().for_each(|p| for c in p.iter_mut() { *c = c.scale(1.0 / nn); });
    d
}

/// `k` punti nati con |ψ|² (una cella scelta con il suo peso, poi un punto a caso dentro la cella).
pub fn nascite(d: &Dirac3D, k: usize, rng: &mut Rng) -> Vec<[f64; 3]> {
    let n = d.n;
    let mut cdf = Vec::with_capacity(d.psi.len());
    let mut s = 0.0;
    for p in &d.psi { s += p.iter().map(|c| c.norm2()).sum::<f64>(); cdf.push(s); }
    (0..k).map(|_| {
        let u = rng.uni() * s;
        let i = cdf.partition_point(|&c| c < u).min(cdf.len() - 1);
        let c = [i / (n * n), (i / n) % n, i % n];
        [0, 1, 2].map(|q| d.x(c[q]) + d.a * (rng.uni() - 0.5))
    }).collect()
}

/// ⟨r⟩ dell'onda (λ_C).
pub fn centro_onda(d: &Dirac3D) -> [f64; 3] {
    let n = d.n;
    let xs: Vec<f64> = (0..n).map(|i| d.x(i)).collect();
    let (m, w) = d.psi.par_iter().enumerate().map(|(s, p)| {
        let w: f64 = p.iter().map(|c| c.norm2()).sum();
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        ([w * xs[i], w * xs[j], w * xs[k]], w)
    }).reduce(|| ([0.0; 3], 0.0), |a, b| ([a.0[0] + b.0[0], a.0[1] + b.0[1], a.0[2] + b.0[2]], a.1 + b.1));
    [m[0] / w, m[1] / w, m[2] / w]
}

/// Centro dei bagliori (λ_C).
pub fn centro_palline(p: &[[f64; 3]]) -> [f64; 3] {
    let s = p.par_iter().map(|q| *q).reduce(|| [0.0; 3], |a, b| [a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
    let k = p.len().max(1) as f64;
    [s[0] / k, s[1] / k, s[2] / k]
}

/// Un tick della regola e dei bagliori. Ogni bagliore cavalca v = j/ρ con la regola del trapezio nel tempo:
/// v₁ letto al tick t nel punto p, v₂ letto al tick t+1 in p + a v₁, passo a (v₁ + v₂)/2 (secondo ordine; il
/// tremito ha un periodo di ~12 tick, un passo di Eulero lo sfaserebbe).
pub fn passo(d: &mut Dirac3D, surf: &mut [[f64; 3]]) {
    let a = d.a;
    let v1: Vec<[f64; 3]> = surf.par_iter().map(|p| surfer_velocity(&d.sample(*p))).collect();
    d.tick();
    let d = &*d;
    surf.par_iter_mut().zip(v1.par_iter()).for_each(|(p, v1)| {
        let q = [p[0] + a * v1[0], p[1] + a * v1[1], p[2] + a * v1[2]];
        let v2 = surfer_velocity(&d.sample(q));
        for c in 0..3 { p[c] += 0.5 * a * (v1[c] + v2[c]); }
    });
}

/// LA CORRENTE ESATTA DELLA REGOLA (docs/SURFISTI_ESATTI.md): la probabilità che passa da ogni link (s, s + e_j)
/// durante il prossimo tick, la parte che va avanti da s meno quella che torna indietro da s + e_j, sotto-passo per
/// sotto-passo nell'ordine della regola. Lettore soltanto: rifà le proiezioni di `Dirac3D::shift` (walk.rs, privata)
/// con le costanti pubbliche HALF e order su una copia di ψ; l'evoluzione resta `Dirac3D::tick`. Debito dichiarato
/// (docs/BAGLIORI.md): con docs/patch/walk_tick_flusso.patch la regola darebbe i flussi da sola, senza la copia.
pub fn flussi(d: &Dirac3D) -> [Vec<f64>; 3] {
    let n = d.n;
    let nb = |s: usize, j: usize, dd: isize| -> usize {
        let (i, k, l) = (s / (n * n), (s / n) % n, s % n);
        let w = |x: usize| ((x as isize + dd).rem_euclid(n as isize)) as usize;
        match j { 0 => (w(i) * n + k) * n + l, 1 => (i * n + w(k)) * n + l, _ => (i * n + k) * n + w(l) }
    };
    let mut psi = d.psi.clone();
    let mut tmp = vec![[C::ZERO; 4]; psi.len()];
    let mut f = [vec![0.0; psi.len()], vec![0.0; psi.len()], vec![0.0; psi.len()]];
    for j in order(d.ticks) {
        let [u, dn] = HALF[j];
        tmp.par_iter_mut().zip(psi.par_iter()).for_each(|(t, p)| {
            t[0] = u[0].conj() * p[0] + u[1].conj() * p[1];
            t[1] = dn[0].conj() * p[0] + dn[1].conj() * p[1];
            t[2] = u[0].conj() * p[2] + u[1].conj() * p[3];
            t[3] = dn[0].conj() * p[2] + dn[1].conj() * p[3];
        });
        let t = &tmp;
        f[j].par_iter_mut().enumerate().for_each(|(s, v)| {
            let q = &t[nb(s, j, 1)];
            *v = t[s][1].norm2() + t[s][2].norm2() - q[0].norm2() - q[3].norm2();
        });
        psi.par_iter_mut().enumerate().for_each(|(s, p)| {
            let lp = t[nb(s, j, 1)][0];
            let lm = t[nb(s, j, -1)][1];
            let rp = t[nb(s, j, -1)][2];
            let rm = t[nb(s, j, 1)][3];
            p[0] = u[0] * lp + dn[0] * lm;
            p[1] = u[1] * lp + dn[1] * lm;
            p[2] = u[0] * rp + dn[0] * rm;
            p[3] = u[1] * rp + dn[1] * rm;
        });
    }
    f
}

/// Un tick con i bagliori portati dalla CORRENTE ESATTA (`flussi`): dentro la cella del sito s la velocità lungo j
/// va linearmente dal flusso della faccia di sinistra a quello della faccia di destra, diviso la probabilità della
/// cella (come gli elementi di Raviart–Thomas); spostamento = v·a, al massimo una cella per asse. La media dei
/// bagliori segue così ⟨r⟩ della regola (Σ flussi · a = d⟨r⟩ esatto), non la corrente ψ†αψ del continuo.
pub fn passo_flusso(d: &mut Dirac3D, surf: &mut [[f64; 3]]) {
    let n = d.n;
    let a = d.a;
    let f = flussi(d);
    let rho: Vec<f64> = d.psi.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum()).collect();
    let x0 = d.x(0);
    surf.par_iter_mut().for_each(|p| {
        let mut ii = [0usize; 3];
        let mut fr = [0f64; 3];
        for c in 0..3 {
            let g = (p[c] - x0) / a;
            let r = g.round();
            ii[c] = (r as i64).rem_euclid(n as i64) as usize;
            fr[c] = g - r + 0.5; // 0 = faccia di sinistra, 1 = di destra
        }
        let s = (ii[0] * n + ii[1]) * n + ii[2];
        let r = rho[s].max(1e-300);
        for c in 0..3 {
            let mut im = ii;
            im[c] = (ii[c] + n - 1) % n;
            let sm = (im[0] * n + im[1]) * n + im[2];
            let v = ((1.0 - fr[c]) * f[c][sm] + fr[c] * f[c][s]) / r;
            p[c] += a * v.clamp(-1.0, 1.0);
        }
    });
    d.tick();
}

/// Il lettore della corrente usato dalla scena: true = corrente esatta della regola (`passo_flusso`, emendamento
/// dopo l'esame E3, docs/BAGLIORI.md), false = spinore interpolato (`passo`, quello pre-registrato).
pub const LETTORE_ESATTO: bool = true;

/// Il disegno: un bagliore tondo per surfista, colore = fase (tolta la rotazione di riposo), luce dalla densità.
pub fn palline(d: &Dirac3D, surf: &[[f64; 3]]) -> Vec<Ball> {
    let t = d.ticks as f64 * d.a;
    let rho_max = d.psi.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).reduce(|| 0.0, f64::max).max(1e-300);
    surf.par_iter().map(|p| {
        let s = d.sample(*p);
        let rho: f64 = s.iter().map(|c| c.norm2()).sum();
        let mut best = 0;
        for q in 1..4 { if s[q].norm2() > s[best].norm2() { best = q; } }
        let ph = s[best].arg() + t;
        let a = (0.03 + 0.30 * (rho / rho_max).powf(0.35)) as f32;
        (schermo(*p).to_array(), 1.0, phase_color(ph, a))
    }).collect()
}

struct Mondo { d: Dirac3D, surf: Vec<[f64; 3]>, part: Particella }

fn nuovo(part: Particella, rng: &mut Rng) -> Mondo {
    let d = pacchetto(N_EL, THETA, SIGMA, spinore_f3());
    let surf = nascite(&d, N_BAGLIORI, rng);
    Mondo { d, surf, part }
}

pub fn avvia(part: Particella) -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 12.0, ricomincia: false, particella: None, fine: false,
        out: Uscita { palline: Arc::new(vec![]), versione: 0, ticks: 0, t: 0.0, centro_palline: [0.0; 3], centro_onda: [0.0; 3],
                      particella: part, ms_tick: 0.0, pronto: false },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let mut rng = Rng::new(2027);
        let mut m = nuovo(part, &mut rng);
        let mut versione = 0u64;
        let mut ms = 0.0;
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, ric, nuova, fine) = {
                let mut g = cc.lock().unwrap();
                let s = (g.comando, g.ritmo, g.ricomincia, g.particella.take(), g.fine);
                g.ricomincia = false;
                s
            };
            if fine { break; }
            let mut cambio = false;
            if let Some(p) = nuova { m = nuovo(p, &mut rng); cambio = true; }
            if ric || m.d.ticks >= TICK_MAX { let p = m.part; m = nuovo(p, &mut rng); cambio = true; }
            let fermo = TICK_FERMO.load(std::sync::atomic::Ordering::Relaxed);
            let cmd = if fermo > 0 && m.d.ticks >= fermo { Comando::Pausa } else { cmd };
            if cmd == Comando::Avanti {
                if LETTORE_ESATTO { passo_flusso(&mut m.d, &mut m.surf); } else { passo(&mut m.d, &mut m.surf); }
                cambio = true;
            }
            if cambio || versione == 0 {
                let p = palline(&m.d, &m.surf);
                versione += 1;
                ms += 0.15 * (t0.elapsed().as_secs_f64() * 1e3 - ms);
                let out = Uscita {
                    palline: Arc::new(p), versione, ticks: m.d.ticks, t: m.d.ticks as f64 * m.d.a,
                    centro_palline: centro_palline(&m.surf), centro_onda: centro_onda(&m.d), particella: m.part, ms_tick: ms, pronto: true,
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

/// Fit di x(t) = c₀ + c₁ t + A sin(ω t + φ) (come zbw_fit di lab/zz_free.py): per ogni ω i quattro coefficienti
/// lineari ai minimi quadrati, ω cercato su una griglia fine e poi raffinato. Restituisce (A, ω).
pub fn fit_tremito(t: &[f64], x: &[f64], w_lo: f64, w_hi: f64) -> (f64, f64) {
    let resid = |w: f64| -> (f64, f64) {
        // base: 1, t, sin ωt, cos ωt
        let mut ata = [[0.0f64; 4]; 4];
        let mut atb = [0.0f64; 4];
        for (ti, xi) in t.iter().zip(x) {
            let b = [1.0, *ti, (w * ti).sin(), (w * ti).cos()];
            for i in 0..4 { atb[i] += b[i] * xi; for j in 0..4 { ata[i][j] += b[i] * b[j]; } }
        }
        let c = risolvi4(ata, atb);
        let mut r = 0.0;
        for (ti, xi) in t.iter().zip(x) {
            let f = c[0] + c[1] * ti + c[2] * (w * ti).sin() + c[3] * (w * ti).cos();
            r += (f - xi) * (f - xi);
        }
        (r, (c[2] * c[2] + c[3] * c[3]).sqrt())
    };
    let k = 2000;
    let mut best = (f64::MAX, w_lo);
    for i in 0..=k {
        let w = w_lo + (w_hi - w_lo) * i as f64 / k as f64;
        let r = resid(w).0;
        if r < best.0 { best = (r, w); }
    }
    // sezione aurea attorno al minimo
    let h = (w_hi - w_lo) / k as f64;
    let (mut a, mut b) = (best.1 - h, best.1 + h);
    let g = 0.618_033_988_75;
    for _ in 0..80 {
        let c1 = b - g * (b - a);
        let c2 = a + g * (b - a);
        if resid(c1).0 < resid(c2).0 { b = c2; } else { a = c1; }
    }
    let w = 0.5 * (a + b);
    (resid(w).1, w)
}

fn risolvi4(mut a: [[f64; 4]; 4], mut r: [f64; 4]) -> [f64; 4] {
    for i in 0..4 {
        let p = (i..4).max_by(|&x, &y| a[x][i].abs().partial_cmp(&a[y][i].abs()).unwrap()).unwrap();
        a.swap(i, p);
        r.swap(i, p);
        if a[i][i].abs() < 1e-300 { return [0.0; 4]; }
        for j in i + 1..4 {
            let f = a[j][i] / a[i][i];
            for k in i..4 { a[j][k] -= f * a[i][k]; }
            r[j] -= f * r[i];
        }
    }
    let mut x = [0.0; 4];
    for i in (0..4).rev() {
        let mut s = r[i];
        for k in i + 1..4 { s -= a[i][k] * x[k]; }
        x[i] = s / a[i][i];
    }
    x
}

//! AN ATOM DRAWN FROM ITS COMPUTED ORBITALS (lab/out/orbatomo_prereg.json): one group of glows for each spin-orbital of the
//! determinant in `assets/esportazioni/<nome>/orbitali.json`, drawn from |R Y_lm|² and carried by the Pauli current
//! v_φ = m/ρ − s ½ ∂ln|ψ|²/∂ρ, the small-velocity limit of the Dirac current the hydrogen surfers ride. The flow is all
//! azimuthal, so each glow turns on its own circle at a constant rate and the density stays exact: no sampler, no drift.

use std::path::Path;

use serde_json::Value;

use crate::bagliori::Storie;
use crate::render::{schermo, Ball};
use crate::rng::Rng;

/// Screen seconds per atomic unit of time (slowed, declared on screen). A glow never goes faster than c.
pub const S_PER_AU: f64 = 4.0;
const C_AU: f64 = 137.036;
const STOPS: [[f32; 3]; 5] = [[0.20, 0.95, 1.00], [0.25, 0.45, 1.00], [0.70, 0.30, 1.00], [1.00, 0.35, 0.70], [1.00, 0.85, 0.30]];

#[derive(Clone, Debug)]
pub struct Orb { pub nome: String, pub l: u32, pub p: Vec<f64> }
#[derive(Clone, Copy, Debug)]
pub struct SpinOrb { pub orb: usize, pub m: i32, pub s: f64 }
#[derive(Clone, Debug)]
pub struct Dati { pub r: Vec<f64>, pub orb: Vec<Orb>, pub so: Vec<SpinOrb> }

impl Dati {
    pub fn leggi(cartella: &Path) -> Option<Dati> {
        let v: Value = serde_json::from_slice(&std::fs::read(cartella.join("orbitali.json")).ok()?).ok()?;
        let nums = |x: &Value| -> Vec<f64> { x.as_array().map(|a| a.iter().filter_map(|y| y.as_f64()).collect()).unwrap_or_default() };
        let orb = v["orbitali"].as_array()?.iter().map(|o| Orb { nome: o["nome"].as_str().unwrap_or("?").into(), l: o["l"].as_u64().unwrap_or(0) as u32, p: nums(&o["P"]) }).collect();
        let so = v["spin_orbitali"].as_array()?.iter().map(|o| SpinOrb { orb: o["orb"].as_u64().unwrap_or(0) as usize, m: o["m"].as_i64().unwrap_or(0) as i32, s: o["s"].as_f64().unwrap_or(1.0) }).collect();
        Some(Dati { r: nums(&v["r"]), orb, so })
    }

    /// P(r) and dP/dr on the table (linear between its points; 0 outside).
    fn p_dp(&self, o: usize, r: f64) -> (f64, f64) {
        let (x, p) = (&self.r, &self.orb[o].p);
        if r <= x[0] || r >= x[x.len() - 1] { return (0.0, 0.0); }
        let i = x.partition_point(|a| *a < r).clamp(1, x.len() - 1);
        let d = (p[i] - p[i - 1]) / (x[i] - x[i - 1]);
        (p[i - 1] + d * (r - x[i - 1]), d)
    }

    /// v_φ of spin-orbital `q` at (ρ, z), limited to c.
    pub fn v_phi(&self, q: SpinOrb, rho: f64, z: f64) -> f64 {
        let r = (rho * rho + z * z).sqrt().max(1e-12);
        let (p, dp) = self.p_dp(q.orb, r);
        // ln|ψ|² = 2 ln|P/r| + ln|Θ|²: Θ² ∝ 1 (l = 0), z²/r² (m = 0), ρ²/r² (|m| = 1)
        let d_rad = if p.abs() > 1e-300 { 2.0 * (dp / p - 1.0 / r) * rho / r } else { 0.0 };
        let d_ang = match (self.orb[q.orb].l, q.m.abs()) { (0, _) => 0.0, (_, 0) => -2.0 * rho / (r * r), _ => 2.0 / rho.max(1e-12) - 2.0 * rho / (r * r) };
        (q.m as f64 / rho.max(1e-12) - q.s * 0.5 * (d_rad + d_ang)).clamp(-C_AU, C_AU)
    }

    /// A point (ρ, z) of spin-orbital q drawn from |R Y_lm|².
    fn pesca(&self, q: SpinOrb, cdf: &[f64], rng: &mut Rng) -> (f64, f64) {
        let u = rng.uni() * cdf[cdf.len() - 1];
        let i = cdf.partition_point(|c| *c < u).clamp(1, cdf.len() - 1);
        let t = (u - cdf[i - 1]) / (cdf[i] - cdf[i - 1]).max(1e-300);
        let r = self.r[i - 1] + t * (self.r[i] - self.r[i - 1]);
        let ct = loop {
            let c = 2.0 * rng.uni() - 1.0;
            let w = match (self.orb[q.orb].l, q.m.abs()) { (0, _) => 1.0, (_, 0) => c * c, _ => 1.0 - c * c };
            if rng.uni() < w { break c; }
        };
        (r * (1.0 - ct * ct).max(0.0).sqrt(), r * ct)
    }
}

/// One glow: where it is on its circle, how fast it turns, what it looks like.
#[derive(Clone, Copy, Debug)]
pub struct Bagl { pub rho: f64, pub z: f64, pub phi0: f64, pub om: f64, pub m: i32, pub cl: u32, pub alfa: f32 }

pub struct Atomo { pub dati: Dati, pub g: Vec<Bagl>, pub t: f64, pub raggio: f64, storie: Storie, z_nucleo: u32 }

impl Atomo {
    pub fn new(dati: Dati, z_nucleo: u32, n_tot: usize) -> Self {
        let mut rng = Rng::new(808);
        let k = n_tot / dati.so.len().max(1);
        let mut g = Vec::with_capacity(n_tot);
        for (j, q) in dati.so.iter().enumerate() {
            let p = &dati.orb[q.orb].p;
            let mut cdf = vec![0.0; dati.r.len()];
            for i in 1..cdf.len() { cdf[i] = cdf[i - 1] + 0.5 * (p[i] * p[i] + p[i - 1] * p[i - 1]) * (dati.r[i] - dati.r[i - 1]); }
            for _ in 0..k {
                let (rho, z) = dati.pesca(*q, &cdf, &mut rng);
                let om = dati.v_phi(*q, rho, z) / rho.max(1e-6);
                g.push(Bagl { rho, z, phi0: std::f64::consts::TAU * rng.uni(), om, m: q.m, cl: q.orb as u32 + j as u32 * 0, alfa: 0.25 });
            }
        }
        // the order of the glows is random, so the level of detail thins all the orbitals alike
        for i in (1..g.len()).rev() { let j = (rng.next_u64() % (i as u64 + 1)) as usize; g.swap(i, j); }
        let mut rr: Vec<f64> = g.iter().map(|b| (b.rho * b.rho + b.z * b.z).sqrt()).collect();
        rr.sort_by(|a, b| a.total_cmp(b));
        let raggio = rr.get((rr.len() as f64 * 0.95) as usize).copied().unwrap_or(3.0).max(0.5);
        let storie = Storie::new(g.iter().map(|b| schermo([b.rho * b.phi0.cos(), b.rho * b.phi0.sin(), b.z]).to_array()), 1.0);
        Atomo { dati, g, t: 0.0, raggio, storie, z_nucleo }
    }

    fn pos(&self, b: &Bagl) -> [f64; 3] { let f = b.phi0 + b.om * self.t; [b.rho * f.cos(), b.rho * f.sin(), b.z] }

    /// `dt` screen seconds.
    pub fn passo(&mut self, dt: f64) {
        self.t += dt / S_PER_AU;
        let posizioni: Vec<[f32; 3]> = self.g.iter().map(|b| schermo(self.pos(b)).to_array()).collect();
        self.storie.registra(posizioni.into_iter());
    }

    pub fn disegna(&self, v: &crate::grafica::VistaWorker, ritmo: f64) -> (Vec<Ball>, Vec<Ball>, Vec<[f32; 3]>) {
        let mut p: Vec<Ball> = self.g.iter().map(|b| {
            let f = b.phi0 + b.om * self.t;
            // colour = the orbital's phase m φ, one starting hue per shell
            let t = ((b.m as f64 * f) as f32 + 1.3 * b.cl as f32).rem_euclid(std::f32::consts::TAU) / std::f32::consts::TAU * 5.0;
            let i = (t.floor() as usize) % 5;
            let s = t - t.floor();
            let (a, c) = (STOPS[i], STOPS[(i + 1) % 5]);
            (schermo([b.rho * f.cos(), b.rho * f.sin(), b.z]).to_array(), 1.0, [a[0] + (c[0] - a[0]) * s, a[1] + (c[1] - a[1]) * s, a[2] + (c[2] - a[2]) * s, b.alfa])
        }).collect();
        let (sc, sa) = self.storie.scie(&p, v.lod_f, v.scia_s, ritmo);
        p.push((schermo([0.0; 3]).to_array(), 1.5 + 0.3 * (self.z_nucleo as f32).sqrt(), [1.0, 0.6, 0.25, 0.9]));
        (p, sc, sa)
    }

    pub fn righe(&self, nome: &str, simbolo: &str, energia: Option<f64>, livello: &str, np: usize, nf: usize) -> (String, String) {
        let conf: Vec<String> = self.dati.so.iter().map(|q| format!("{}{}{:+}", self.dati.orb[q.orb].nome, if q.s > 0.0 { "↑" } else { "↓" }, q.m)).collect();
        (format!("{nome} ({simbolo}) · livello {} · orbitali calcolati: {} · esami {np} PASS / {nf} FAIL", livello.to_uppercase(), conf.join(" ")),
         format!("{}palline su |ψ|² di ogni orbitale, portate dalla corrente di Pauli (orbitale + spin) · elettroni indipendenti: la correlazione è solo nelle energie · 1 u.a. di tempo = {S_PER_AU} s", energia.map(|e| format!("E = {e:.6} Eh · ")).unwrap_or_default()))
    }
}

#[cfg(test)]
mod esami {
    use super::*;

    fn dati(n: &str) -> Dati { Dati::leggi(&crate::esportazioni::cartella_assets().join("esportazioni").join(n)).expect("python lab/tavg_f5_orbitali.py") }

    /// O1, O2, O4.
    #[test]
    fn orb_palline() {
        let d = dati("O");
        let integra = |f: &dyn Fn(usize) -> f64| -> f64 { (1..d.r.len()).map(|i| 0.5 * (f(i) + f(i - 1)) * (d.r[i] - d.r[i - 1])).sum() };
        for o in &d.orb {
            let n = integra(&|i| o.p[i] * o.p[i]);
            println!("O1 {} norma {n:.6}", o.nome);
            assert!((n - 1.0).abs() <= 1e-3);
        }
        let a = Atomo::new(d.clone(), 8, 200_000);
        let n = a.g.len() as f64;
        let rm: f64 = a.g.iter().map(|b| (b.rho * b.rho + b.z * b.z).sqrt()).sum::<f64>() / n;
        let atteso: f64 = d.so.iter().map(|q| { let o = &d.orb[q.orb]; integra(&|i| o.p[i] * o.p[i] * d.r[i]) }).sum::<f64>() / d.so.len() as f64;
        println!("O2 <r> palline {rm:.4}, orbitali {atteso:.4}");
        assert!((rm / atteso - 1.0).abs() <= 0.02);
        // O4: ⟨z²⟩ / ⟨x²⟩ from the determinant: ⟨cos²θ⟩ = 1/3 (l = 0), 3/5 (m = 0), 1/5 (|m| = 1)
        let (mut z2, mut r2) = (0.0, 0.0);
        for q in &d.so {
            let o = &d.orb[q.orb];
            let rr = integra(&|i| o.p[i] * o.p[i] * d.r[i] * d.r[i]);
            let c = match (o.l, q.m.abs()) { (0, _) => 1.0 / 3.0, (_, 0) => 0.6, _ => 0.2 };
            z2 += rr * c;
            r2 += rr;
        }
        let atteso = z2 / (0.5 * (r2 - z2));
        let (bz, bx): (f64, f64) = a.g.iter().fold((0.0, 0.0), |s, b| (s.0 + b.z * b.z, s.1 + 0.5 * b.rho * b.rho));
        println!("O4 <z^2>/<x^2> palline {:.4}, determinante {atteso:.4}", bz / bx);
        assert!((bz / bx / atteso - 1.0).abs() <= 0.05);
    }

    /// O3: hydrogen-like 1s: v_φ = Z sinθ.
    #[test]
    fn orb_corrente_1s() {
        let z = 8.0f64;
        let r: Vec<f64> = (1..20000).map(|i| i as f64 * 2e-4).collect();
        let p: Vec<f64> = r.iter().map(|x| 2.0 * z.powf(1.5) * x * (-z * x).exp()).collect();
        let d = Dati { r, orb: vec![Orb { nome: "1s".into(), l: 0, p }], so: vec![] };
        let v = d.v_phi(SpinOrb { orb: 0, m: 0, s: 1.0 }, 0.1, 0.1);
        let atteso = z * std::f64::consts::FRAC_1_SQRT_2;
        println!("O3 v_phi {v:.4} contro {atteso:.4}");
        assert!((v / atteso - 1.0).abs() <= 0.01);
    }
}

//! The numbers the viewer uses, all produced by our own runs (`atomo --misura`, lab/atomo_misura.sh) and
//! stored as `chiave = valore` lines in assets/atomo/risultati.txt. A copy is compiled into the program,
//! so the executable works even without the assets folder.

use crate::carbonio::Param;

pub const FILE: &str = "assets/atomo/risultati.txt";
const DENTRO: &str = include_str!("../../assets/atomo/risultati.txt");

#[derive(Clone, Debug)]
pub struct Risultati {
    pub righe: Vec<(String, f64)>,
}

impl Risultati {
    pub fn carica() -> Self {
        let txt = std::fs::read_to_string(FILE).unwrap_or_else(|_| DENTRO.to_string());
        Self::da_testo(&txt)
    }
    pub fn da_testo(txt: &str) -> Self {
        let righe = txt.lines().filter_map(|l| {
            let l = l.split('#').next()?.trim();
            let (k, v) = l.split_once('=')?;
            Some((k.trim().to_string(), v.trim().parse().ok()?))
        }).collect();
        Risultati { righe }
    }
    pub fn get(&self, k: &str) -> Option<f64> { self.righe.iter().rev().find(|(a, _)| a == k).map(|(_, v)| *v) }
    pub fn o(&self, k: &str, d: f64) -> f64 { self.get(k).unwrap_or(d) }
    pub fn metti(&mut self, k: &str, v: f64) {
        self.righe.retain(|(a, _)| a != k);
        self.righe.push((k.to_string(), v));
    }
    pub fn salva(&self, intestazione: &str) {
        let mut s = String::new();
        for l in intestazione.lines() { s.push_str(&format!("# {}\n", l)); }
        for (k, v) in &self.righe { s.push_str(&format!("{} = {:.10}\n", k, v)); }
        let _ = std::fs::create_dir_all("assets/atomo");
        std::fs::write(FILE, s).expect("cannot write assets/atomo/risultati.txt");
    }

    /// Optimised guide parameters of a state ("c", "ione", "ecc").
    pub fn param(&self, s: &str) -> Param {
        let d = Param { z1: 5.58, zs: 2.02, zp: 1.78, z3: 0.30, d3: 0.53, b: 1.0 };
        Param {
            z1: self.o(&format!("{s}_z1"), d.z1),
            zs: self.o(&format!("{s}_zs"), d.zs),
            zp: self.o(&format!("{s}_zp"), d.zp),
            z3: self.o(&format!("{s}_z3"), d.z3),
            d3: self.o(&format!("{s}_d3"), d.d3),
            b: 1.0,
        }
    }
    pub fn metti_param(&mut self, s: &str, p: &Param) {
        for (k, v) in [("z1", p.z1), ("zs", p.zs), ("zp", p.zp), ("z3", p.z3), ("d3", p.d3)] {
            self.metti(&format!("{s}_{k}"), v);
        }
    }
    /// DMC energies (hartree) of the three states.
    pub fn e_c(&self) -> f64 { self.o("dmc_c", -37.825) }
    pub fn e_ione(&self) -> f64 { self.o("dmc_ione", -37.412) }
    pub fn e_ecc(&self) -> f64 { self.o("dmc_ecc", -37.542) }
    /// Errors of ΔE (hartree): sum in quadrature.
    pub fn err_ecc(&self) -> f64 { (self.o("dmc_c_err", 0.004).powi(2) + self.o("dmc_ecc_err", 0.004).powi(2)).sqrt() }
    pub fn err_ione(&self) -> f64 { (self.o("dmc_c_err", 0.004).powi(2) + self.o("dmc_ione_err", 0.004).powi(2)).sqrt() }
}

//! The world of the viewer: one carbon atom (its balls), the nucleus moving with the temperature,
//! photons and the photo-electron. Deterministic: the state after a frame depends only on the state before
//! and the frame's `Input` (all random draws come from the seeded generators inside), which is what makes
//! reversible time possible (tempo.rs).
//!
//! Units: positions in bohr, physical time in seconds.

use crate::carbonio::{Dmc, Guida, Rng, Stato};
use crate::costanti as k;
use crate::risultati::Risultati;

/// Carbon: nuclear charge (input).
pub const Z: f64 = 6.0;
/// Mass of the ¹²C nucleus (kg): 12 u (definition of u) minus 6 electron masses (the electrons' binding,
/// ~1 keV, changes it by 10⁻⁷: neglected).
pub fn massa_nucleo_kg() -> f64 { 12.0 * 1.660_539_066_60e-27 - 6.0 * k::M_E }
pub const K_B: f64 = 1.380_649e-23;
/// Gas pressure around the atom for the collision time (1 atm).
pub const PRESSIONE_PA: f64 = 101_325.0;
/// Speed of light in bohr / s.
pub fn c_bohr_s() -> f64 { k::C_L / k::bohr_m() }
/// Photon energy (hartree) at wavelength λ (nm, vacuum): E = hc/λ.
pub fn hartree_da_nm(l: f64) -> f64 { 2.0 * std::f64::consts::PI * k::HBAR * k::C_L / (l * 1e-9) / k::hartree_j() }

/// Live balls and tick (the same tick as the validation runs).
pub const N_BALLS: usize = 1500;
pub const TAU: f64 = 0.005;
/// Photons start this far from the nucleus (bohr) and are forgotten beyond `R_FINE`.
pub const R_PARTENZA: f64 = 40.0;
pub const R_FINE: f64 = 70.0;
/// Waits that in reality last nanoseconds (spontaneous emission, recombination) are shown compressed:
/// mean screen seconds. Their true durations are not computed (debts, docs/ATOMO.md).
pub const ATTESA_ECC_S: f64 = 3.0;
pub const ATTESA_IONE_S: f64 = 5.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Esito { InVolo, Assorbito, Ionizza, Passa }

#[derive(Clone, Debug)]
pub struct Fotone {
    /// centre of the packet (bohr) and direction
    pub pos: [f64; 3],
    pub dir: [f64; 3],
    pub lambda_nm: f64,
    /// outgoing photon: a spherical dipole wave around `pos`, radius `raggio`, axis `dir`
    pub uscente: bool,
    pub raggio: f64,
    pub esito: Esito,
    pub id: u64,
}

#[derive(Clone, Debug)]
pub struct Fotoelettrone { pub pos: [f64; 3], pub vel: [f64; 3], pub ev: f64 }

/// Everything except the balls (cheap to copy every frame for the history).
#[derive(Clone, Debug)]
pub struct Piccolo {
    pub t_s: f64,
    pub t_schermo: f64,
    pub stato: Stato,
    pub nucleo: [f64; 3],
    pub v_nucleo: [f64; 3],
    pub fotoni: Vec<Fotone>,
    pub elettroni: Vec<Fotoelettrone>,
    /// screen time at which the excited atom / the ion returns to the ground state
    pub attesa: Option<f64>,
    pub rng: Rng,
    pub tick: u64,
    pub prossimo_id: u64,
    pub registro: Vec<String>,
    pub assorbiti: u32,
    pub ionizzati: u32,
    pub passati: u32,
    pub emessi: u32,
}

#[derive(Clone)]
pub struct Mondo {
    pub p: Piccolo,
    pub dmc: Dmc,
}

/// What the user (and the clock) did in one frame.
#[derive(Clone, Copy, Debug, Default)]
pub struct Input {
    pub dt_s: f64,
    pub dt_schermo: f64,
    pub ticks: u32,
    pub t_kelvin: f64,
    /// shoot a photon: wavelength (nm) and direction (unit)
    pub spara: Option<(f64, [f64; 3])>,
}

/// The constants of the three states, from our own DMC runs.
#[derive(Clone, Debug)]
pub struct Livelli {
    pub e: [f64; 3],
    pub de_ecc: f64,
    pub de_ecc_err: f64,
    pub ip: f64,
    pub ip_err: f64,
}

impl Livelli {
    pub fn da(r: &Risultati) -> Self {
        Livelli { e: [r.e_c(), r.e_ecc(), r.e_ione()], de_ecc: r.e_ecc() - r.e_c(), de_ecc_err: r.err_ecc(),
                  ip: r.e_ione() - r.e_c(), ip_err: r.err_ione() }
    }
    pub fn energia(&self, s: Stato) -> f64 {
        match s { Stato::Fondamentale => self.e[0], Stato::Eccitato => self.e[1], Stato::Ione => self.e[2] }
    }
    pub fn lambda_ecc_nm(&self) -> f64 { k::nm_from_hartree(self.de_ecc) }
    pub fn lambda_ip_nm(&self) -> f64 { k::nm_from_hartree(self.ip) }
}

pub fn guida(r: &Risultati, s: Stato) -> Guida {
    let p = match s { Stato::Fondamentale => r.param("c"), Stato::Eccitato => r.param("ecc"), Stato::Ione => r.param("ione") };
    Guida { z: Z, stato: s, p }
}

fn norm(v: [f64; 3]) -> f64 { (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() }

/// Direction drawn with the dipole pattern sin²θ about `asse` (rejection sampling).
fn dir_dipolo(rng: &mut Rng, asse: [f64; 3]) -> [f64; 3] {
    loop {
        let v = [rng.gauss(), rng.gauss(), rng.gauss()];
        let n = norm(v);
        let u = [v[0] / n, v[1] / n, v[2] / n];
        let c = u[0] * asse[0] + u[1] * asse[1] + u[2] * asse[2];
        if rng.uni() < 1.0 - c * c { return u; }
    }
}

impl Mondo {
    pub fn new(r: &Risultati, seed: u64) -> Self {
        let mut dmc = Dmc::new(guida(r, Stato::Fondamentale), TAU, N_BALLS, seed);
        // a short start of the rule, so the first picture is already the cloud
        for _ in 0..60 { dmc.passo(); }
        Mondo {
            p: Piccolo { t_s: 0.0, t_schermo: 0.0, stato: Stato::Fondamentale, nucleo: [0.0; 3], v_nucleo: [0.0; 3],
                         fotoni: vec![], elettroni: vec![], attesa: None, rng: Rng::new(seed ^ 0x7EA), tick: 0,
                         prossimo_id: 1, registro: vec![], assorbiti: 0, ionizzati: 0, passati: 0, emessi: 0 },
            dmc,
        }
    }

    fn log(&mut self, s: String) {
        self.p.registro.push(s);
        if self.p.registro.len() > 6 { self.p.registro.remove(0); }
    }

    /// Change the state keeping the balls (the cloud relaxes into the new shape under the rule).
    fn passa_a(&mut self, s: Stato, r: &Risultati) {
        let old = self.p.stato;
        let (ne0, nu0) = (old.ne(), old.up().len());
        let n = self.dmc.n();
        let mut x = Vec::with_capacity(n * 3 * s.ne());
        for i in 0..n {
            let w = &self.dmc.x[i * 3 * ne0..(i + 1) * 3 * ne0];
            match (old.ne(), s.ne()) {
                (a, b) if a == b => x.extend_from_slice(w),
                (6, 5) => {
                    // the 4th up electron leaves (it is the one the photon took away)
                    x.extend_from_slice(&w[..3 * (nu0 - 1)]);
                    x.extend_from_slice(&w[3 * nu0..]);
                }
                _ => {
                    // an electron arrives: placed at random around the ion (≈ 2 bohr), then the rule acts
                    x.extend_from_slice(&w[..3 * nu0]);
                    for _ in 0..3 { x.push(2.0 * self.p.rng.gauss()); }
                    x.extend_from_slice(&w[3 * nu0..]);
                }
            }
        }
        self.dmc.da_configurazioni(guida(r, s), x);
        self.p.stato = s;
    }

    pub fn passo(&mut self, inp: &Input, r: &Risultati, lv: &Livelli) {
        let dt = inp.dt_s;
        self.p.t_s += dt;
        self.p.t_schermo += inp.dt_schermo;
        // ── the balls ──
        for _ in 0..inp.ticks {
            self.dmc.passo();
            self.p.tick += 1;
        }
        // ── the nucleus in a gas at temperature T (Maxwell–Boltzmann; collisions at the mean free time) ──
        let t = inp.t_kelvin.max(0.0);
        if t > 0.0 && dt > 0.0 {
            let m = massa_nucleo_kg();
            let sig_v = (K_B * t / m).sqrt() / k::bohr_m(); // bohr/s per component
            let v_mean = (8.0 * K_B * t / (std::f64::consts::PI * m)).sqrt();
            let dens = PRESSIONE_PA / (K_B * t);
            // collision diameter: twice the radius of OUR cloud (⟨r⟩ of the valence balls ≈ 1.7 bohr)
            let d = 2.0 * 1.7 * k::bohr_m();
            let tau_c = 1.0 / (dens * std::f64::consts::PI * d * d * v_mean * 2f64.sqrt());
            let l_box = (K_B * t / PRESSIONE_PA).powf(1.0 / 3.0) / k::bohr_m();
            // move in sub-steps of at most a tenth of the collision time
            let nsub = ((dt / tau_c * 10.0).ceil() as usize).clamp(1, 200);
            let h = dt / nsub as f64;
            for _ in 0..nsub {
                if self.p.v_nucleo == [0.0; 3] || self.p.rng.uni() < 1.0 - (-h / tau_c).exp() {
                    for c in 0..3 { self.p.v_nucleo[c] = sig_v * self.p.rng.gauss(); }
                }
                for c in 0..3 {
                    self.p.nucleo[c] += self.p.v_nucleo[c] * h;
                    // the neighbours (mean spacing at 1 atm) as walls
                    if self.p.nucleo[c].abs() > 0.5 * l_box {
                        self.p.nucleo[c] = self.p.nucleo[c].signum() * 0.5 * l_box;
                        self.p.v_nucleo[c] = -self.p.v_nucleo[c];
                    }
                }
            }
        } else if t == 0.0 {
            self.p.v_nucleo = [0.0; 3];
        }
        // ── a new photon ──
        if let Some((lam, dir)) = inp.spara {
            let n = self.p.nucleo;
            let id = self.p.prossimo_id;
            self.p.prossimo_id += 1;
            self.p.fotoni.push(Fotone { pos: [n[0] - R_PARTENZA * dir[0], n[1] - R_PARTENZA * dir[1], n[2] - R_PARTENZA * dir[2]],
                                        dir, lambda_nm: lam, uscente: false, raggio: 0.0, esito: Esito::InVolo, id });
        }
        // ── photons fly at c; the incoming one meets the atom when it passes the nucleus ──
        let c = c_bohr_s();
        let mut decisioni = vec![];
        for (i, f) in self.p.fotoni.iter_mut().enumerate() {
            if f.uscente { f.raggio += c * dt; continue; }
            let rel0: f64 = (0..3).map(|q| (f.pos[q] - self.p.nucleo[q]) * f.dir[q]).sum();
            if f.esito == Esito::Assorbito || f.esito == Esito::Ionizza { continue; }
            for q in 0..3 { f.pos[q] += f.dir[q] * c * dt; }
            let rel1: f64 = (0..3).map(|q| (f.pos[q] - self.p.nucleo[q]) * f.dir[q]).sum();
            if f.esito == Esito::InVolo && rel0 < 0.0 && rel1 >= 0.0 { decisioni.push(i); }
        }
        for i in decisioni { self.incontro(i, r, lv); }
        // ── photo-electrons fly away ──
        for e in self.p.elettroni.iter_mut() { for q in 0..3 { e.pos[q] += e.vel[q] * dt; } }
        let nn = self.p.nucleo;
        self.p.elettroni.retain(|e| norm([e.pos[0] - nn[0], e.pos[1] - nn[1], e.pos[2] - nn[2]]) < 5.0 * R_FINE);
        self.p.fotoni.retain(|f| {
            if f.uscente { f.raggio < R_FINE } else {
                f.esito == Esito::Assorbito || f.esito == Esito::Ionizza
                    || norm([f.pos[0] - nn[0], f.pos[1] - nn[1], f.pos[2] - nn[2]]) < R_FINE
            }
        });
        // absorbed photons vanish once the atom has taken them
        self.p.fotoni.retain(|f| f.uscente || f.esito == Esito::InVolo || f.esito == Esito::Passa);
        // ── the compressed wait ends: back to the ground state, one photon out ──
        if let Some(ts) = self.p.attesa {
            if self.p.t_schermo >= ts {
                self.p.attesa = None;
                let de = match self.p.stato {
                    Stato::Eccitato => lv.de_ecc,
                    Stato::Ione => {
                        // recombination with an electron at rest: the photon carries the binding energy
                        self.p.elettroni.clear();
                        lv.ip
                    }
                    Stato::Fondamentale => 0.0,
                };
                if de > 0.0 {
                    let was = self.p.stato;
                    self.passa_a(Stato::Fondamentale, r);
                    let asse = [0.0, 0.0, 1.0];
                    let dir = dir_dipolo(&mut self.p.rng, asse);
                    let lam = k::nm_from_hartree(de);
                    let id = self.p.prossimo_id;
                    self.p.prossimo_id += 1;
                    self.p.fotoni.push(Fotone { pos: self.p.nucleo, dir, lambda_nm: lam, uscente: true, raggio: 0.0, esito: Esito::Passa, id });
                    self.p.emessi += 1;
                    self.log(if was == Stato::Ione {
                        format!("ricombinazione: esce un fotone di {:.1} nm", lam)
                    } else {
                        format!("diseccitazione 3s → 2p: esce un fotone di {:.1} nm", lam)
                    });
                }
            }
        }
    }

    /// The photon reaches the atom: what happens is decided by OUR computed levels only.
    fn incontro(&mut self, i: usize, r: &Risultati, lv: &Livelli) {
        let (lam, dir) = (self.p.fotoni[i].lambda_nm, self.p.fotoni[i].dir);
        // photon energy in the atom's frame (first-order Doppler from the thermal velocity)
        let vpar = (0..3).map(|q| self.p.v_nucleo[q] * dir[q]).sum::<f64>() / c_bohr_s();
        let e = hartree_da_nm(lam) * (1.0 - vpar);
        let stato = self.p.stato;
        let ip_qui = lv.e[2] - lv.energia(stato);
        if stato != Stato::Ione && e >= ip_qui {
            self.p.fotoni[i].esito = Esito::Ionizza;
            let ke = e - ip_qui;
            // speed of the photo-electron from its kinetic energy (non-relativistic, a.u.: v = √(2E))
            let v_au = (2.0 * ke).sqrt();
            let v = v_au / k::au_time_s(); // bohr / s
            let u = dir_dipolo(&mut self.p.rng, [0.0, 0.0, 1.0]);
            let n = self.p.nucleo;
            self.p.elettroni.push(Fotoelettrone { pos: n, vel: [u[0] * v, u[1] * v, u[2] * v], ev: ke * k::hartree_ev() });
            self.passa_a(Stato::Ione, r);
            self.p.ionizzati += 1;
            let w = self.p.t_schermo - ATTESA_IONE_S * (1.0 - self.p.rng.uni()).ln();
            self.p.attesa = Some(w);
            self.log(format!("ionizzazione: l'elettrone esce con {:.2} eV", ke * k::hartree_ev()));
        } else if stato == Stato::Fondamentale && (e - lv.de_ecc).abs() <= lv.de_ecc_err.max(1e-6) {
            self.p.fotoni[i].esito = Esito::Assorbito;
            self.passa_a(Stato::Eccitato, r);
            self.p.assorbiti += 1;
            let w = self.p.t_schermo - ATTESA_ECC_S * (1.0 - self.p.rng.uni()).ln();
            self.p.attesa = Some(w);
            self.log(format!("assorbito ({:.1} nm): un elettrone 2p sale sul 3s", lam));
        } else {
            self.p.fotoni[i].esito = Esito::Passa;
            self.p.passati += 1;
            self.log(format!("{:.1} nm: nessun livello calcolato a questa energia, passa", lam));
        }
    }
}

//! The event: one photon meets one helium atom.
//!
//! Every number comes from the balls (assets/molecole/risultati.txt, `molecole --misura elio1p|eliocampo`):
//! the energies of 1¹S and 2¹P, the line λ, the dipole matrix element → f, A, lifetime, the polarizability.
//! Standard quantum mechanics decides what happens:
//!   • one photon (a Gaussian beam focused on the atom, waist w0 = λ) is taken with probability
//!     P = σ(λ)/A_eff, σ = (8π/3)k⁴|α(ω)|² (= 3λ²/2π at the line centre, Rayleigh far from it);
//!   • a taken resonant photon is ABSORBED: given that outcome, the atom's (conditional) wave function is
//!     a ψ_g + b e^{−iΔE t} ψ_e, with |b|² rising while the packet passes — its density oscillates at the
//!     transition frequency (dipole "sloshing"); then the atom is in 2¹P;
//!   • it re-emits after an exponentially distributed time with the computed lifetime, the new photon
//!     going out as a dipole wave (sin²θ about the polarization), its photon ball along a direction drawn
//!     from sin²θ; the cloud goes back to 1¹S through the same kind of superposition;
//!   • a taken photon far from resonance is Rayleigh scattered (no real excitation).

use crate::config::*;
use crate::dmc::{Dmc, Rng, Sistema};
use crate::elio::{self, LuceElio};
use crate::misura::Risultati;
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Fase { Attesa, Volo, Eccitato, Emissione }

impl Fase {
    pub fn nome(self) -> &'static str {
        match self {
            Fase::Attesa => "attesa del prossimo fotone",
            Fase::Volo => "fotone in volo",
            Fase::Eccitato => "atomo eccitato 1s2p ¹P",
            Fase::Emissione => "emissione spontanea",
        }
    }
}

/// The incoming photon.
#[derive(Clone, Debug)]
pub struct FotoneIn {
    pub id: u64,
    /// time (s) when it was sent
    pub t0: f64,
    pub lambda_nm: f64,
    pub lam_draw: f64,
    /// taken by the atom (probability P), and whether that is a real absorption (resonance)
    pub preso: bool,
    pub risonante: bool,
    /// start of its photon ball (screen x, y at z = Z_START)
    pub surf0: [f64; 2],
    /// the atom has already reacted to it
    pub fatto: bool,
}

/// The outgoing photon (emitted or scattered): a dipole wave leaving the atom.
#[derive(Clone, Debug)]
pub struct FotoneOut {
    pub id: u64,
    /// time (s) when the centre of its envelope leaves the atom
    pub t0: f64,
    pub lam_draw: f64,
    /// direction of its photon ball (screen coordinates, unit)
    pub dir: [f64; 3],
}

#[derive(Default, Clone, Copy, Debug)]
pub struct Contatori { pub inviati: u64, pub assorbiti: u64, pub passati: u64, pub emessi: u64, pub diffusi: u64 }

pub struct Scontro {
    pub luce: LuceElio,
    pub lambda_nm: f64,
    pub t: f64,
    pub fase: Fase,
    pub t_fase: f64,
    pub fin: Option<FotoneIn>,
    pub fout: Option<FotoneOut>,
    /// |b|² and the relative phase φ = ΔE t of the superposition (on screen: slowed)
    pub b2: f64,
    pub phi: f64,
    /// drawn waiting time of the excited atom (screen s) and the physical one (s)
    pub attesa: f64,
    pub attesa_fisica: f64,
    /// random numbers of the EVENT (photon taken, waiting time, directions) — independent of the frame rate
    pub rng: Rng,
    /// random numbers of the DRAWING (which balls are drawn, the surfer's start)
    pub rng_vis: Rng,
    pub cont: Contatori,
    pub g: Dmc,
    pub e: Dmc,
    pub sg: Sistema,
    pub se: Sistema,
    /// √(∫ψ_g²/∫ψ_e²): turns the guide ratio ψ_e/ψ_g into the ratio of the normalised functions
    pub norma: f64,
    /// the surfer: one configuration point of the two electrons (bohr)
    pub surfer: [f64; 6],
    pub surfer_in_moto: bool,
    pub pausa: bool,
    /// dipole ⟨z1 + z2⟩ (bohr) of the cloud drawn in the last frame
    pub dipolo: f64,
    next_id: u64,
}

/// Cumulative normal distribution (Abramowitz–Stegun 7.1.26 for erf, |error| < 1.5e-7).
pub fn phi_norm(x: f64) -> f64 {
    let z = x.abs() / 2f64.sqrt();
    let t = 1.0 / (1.0 + 0.327_591_1 * z);
    let y = 1.0 - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t + 0.254_829_592) * t * (-z * z).exp();
    if x >= 0.0 { 0.5 * (1.0 + y) } else { 0.5 * (1.0 - y) }
}

/// Drawn wavelength (screen units) of a photon of λ nm — not to scale, see `config.rs`.
pub fn lambda_draw(lambda_nm: f64) -> f64 { LAMBDA_DRAW_RES * (lambda_nm / 58.4334).sqrt() }

/// A unit vector drawn from the dipole pattern sin²θ about the screen Y axis (the polarization).
pub fn dir_dipolo(rng: &mut Rng) -> [f64; 3] {
    loop {
        let ct = 2.0 * rng.uni() - 1.0;
        let ph = TAU * rng.uni();
        let st = (1.0 - ct * ct).sqrt();
        // sin²θ about Y: θ measured from Y, cos θ = ct
        if rng.uni() < st * st {
            return [st * ph.cos(), ct, st * ph.sin()];
        }
    }
}

impl Scontro {
    pub fn new(res: &Risultati, seed: u64) -> Self {
        let luce = LuceElio::da(res).expect("helium light not computed: run `molecole --misura elio1p` and `--misura eliocampo`");
        let sg = elio::sys_fondamentale(0.0);
        let se = elio::sys_1p(luce.kb);
        let mut g = Dmc::new(sg.clone(), TAU_LIVE, N_POOL, seed ^ 11);
        let mut e = Dmc::new(se.clone(), TAU_LIVE, N_POOL, seed ^ 12);
        g.metropolis = true;
        e.metropolis = true;
        // equilibrate the two clouds (the excited one needs its 2p electron to spread out)
        for _ in 0..(5.0 / TAU_LIVE) as usize { g.passo(); }
        for _ in 0..(25.0 / TAU_LIVE) as usize { e.passo(); }
        let mut s = Scontro {
            lambda_nm: luce.lambda_nm, luce, t: 0.0, fase: Fase::Attesa, t_fase: PAUSA - 0.5, fin: None, fout: None,
            b2: 0.0, phi: 0.0, attesa: 0.0, attesa_fisica: 0.0, rng: Rng::new(seed), rng_vis: Rng::new(seed ^ 0x5EED), cont: Contatori::default(),
            g, e, sg, se, norma: (1.0 / luce.norma_e_g).sqrt(), surfer: [0.0; 6], surfer_in_moto: false, pausa: false,
            dipolo: 0.0, next_id: 1,
        };
        s.nuovo_surfer_fondamentale();
        s
    }

    /// Envelope of the packets: σ_z = 1 drawn wavelength.
    pub fn sigma_z(lam_draw: f64) -> f64 { lam_draw }

    /// Angular frequency (rad per screen second) at which the drawn superposition oscillates: the same as
    /// the drawn wave of the photon passing the atom (c on screen / drawn λ).
    pub fn omega_scr(&self) -> f64 { TAU * C_SCR / lambda_draw(self.lambda_nm) }

    /// Slow-motion factor of the oscillation: real transition frequency / on-screen frequency.
    pub fn rallentamento_osc(&self) -> f64 { self.luce.omega0_s() / self.omega_scr() }

    pub fn prob(&self) -> f64 { self.luce.prob(self.lambda_nm) }

    fn nuovo_surfer_fondamentale(&mut self) {
        let n = self.g.n();
        let i = (self.rng_vis.next_u64() % n as u64) as usize;
        self.surfer.copy_from_slice(&self.g.x[6 * i..6 * i + 6]);
        self.surfer_in_moto = false;
    }

    fn invia(&mut self) {
        let p = self.prob();
        let lam_draw = lambda_draw(self.lambda_nm);
        // its photon ball starts where the beam's intensity puts it (Gaussian beam, waist w0 = λ at the atom)
        let w0 = lam_draw;
        let zr = PI * w0 * w0 / lam_draw;
        let w = w0 * (1.0 + (Z_START / zr).powi(2)).sqrt();
        let (gx, gy) = (self.rng.gauss(), self.rng.gauss());
        let f = FotoneIn {
            id: self.next_id, t0: self.t, lambda_nm: self.lambda_nm, lam_draw,
            preso: self.rng.uni() < p, risonante: self.luce.risonante(self.lambda_nm),
            surf0: [0.5 * w * gx, 0.5 * w * gy], fatto: false,
        };
        self.next_id += 1;
        self.cont.inviati += 1;
        eprintln!("FASE volo t={:.2} lambda={:.4} P={:.3e} preso={}", self.t, self.lambda_nm, p, f.preso);
        self.fin = Some(f);
        self.fase = Fase::Volo;
        self.t_fase = 0.0;
    }

    fn emetti(&mut self, t0: f64, lam_draw: f64) {
        let dir = dir_dipolo(&mut self.rng);
        self.fout = Some(FotoneOut { id: self.next_id, t0, lam_draw, dir });
        self.next_id += 1;
    }

    /// Time (s after sending) at which the centre of the incoming packet is at the atom.
    pub fn t_arrivo() -> f64 { -Z_START / C_SCR }

    /// One frame of the event (dt in screen seconds).
    pub fn avanza(&mut self, dt: f64) {
        for _ in 0..TICKS_PER_FRAME {
            self.g.passo();
            self.e.passo();
        }
        if self.pausa { return; }
        self.t += dt;
        self.t_fase += dt;
        let w = self.omega_scr();
        match self.fase {
            Fase::Attesa => {
                self.b2 = 0.0;
                if self.t_fase > PAUSA { self.invia(); }
            }
            Fase::Volo => {
                let f = self.fin.clone().unwrap();
                let st = Self::sigma_z(f.lam_draw) / C_SCR;
                let s = (self.t - f.t0 - Self::t_arrivo()) / st;
                // marks for an external screen capture (stderr)
                let zc = Z_START + C_SCR * (self.t - f.t0);
                if zc >= -38.0 && zc - C_SCR * dt < -38.0 { eprintln!("MOMENTO arrivo t={:.2}", self.t); }
                if self.b2 < 0.4 && phi_norm(s) >= 0.4 && f.preso && f.risonante { eprintln!("MOMENTO oscillazione t={:.2}", self.t); }
                if f.preso && f.risonante {
                    if self.b2 == 0.0 && s > -3.0 {
                        self.nuovo_surfer_fondamentale();
                        self.surfer_in_moto = true;
                        self.phi = 0.0;
                        eprintln!("FASE assorbimento t={:.2}", self.t);
                    }
                    if s > -3.0 {
                        self.b2 = phi_norm(s).max(1e-9);
                        self.phi += w * dt;
                        self.muovi_surfer(dt);
                    }
                    if s > 3.5 {
                        self.b2 = 1.0;
                        self.surfer_in_moto = false;
                        self.cont.assorbiti += 1;
                        // spontaneous emission after an exponential time with the computed lifetime
                        self.attesa_fisica = -(1.0 - self.rng.uni()).ln() * self.luce.vita_s();
                        self.attesa = self.attesa_fisica * S_DECAY;
                        eprintln!("FASE eccitato t={:.2} attesa={:.3} ns", self.t, self.attesa_fisica * 1e9);
                        self.fin = None;
                        self.fase = Fase::Eccitato;
                        self.t_fase = 0.0;
                    }
                } else {
                    if f.preso && !f.fatto && s > 0.0 {
                        // Rayleigh scattering: the photon leaves as a dipole wave, the atom stays in 1¹S
                        self.cont.diffusi += 1;
                        self.emetti(self.t, f.lam_draw);
                        if let Some(ff) = self.fin.as_mut() { ff.fatto = true; }
                        eprintln!("FASE diffusione t={:.2}", self.t);
                    }
                    let zc = Z_START + C_SCR * (self.t - f.t0);
                    let out_gone = self.fout.as_ref().map_or(true, |o| C_SCR * (self.t - o.t0) > Z_END + 10.0);
                    if zc > Z_END + 3.0 * Self::sigma_z(f.lam_draw) && out_gone {
                        if !f.preso { self.cont.passati += 1; }
                        self.fin = None;
                        self.fout = None;
                        self.fase = Fase::Attesa;
                        self.t_fase = 0.0;
                        eprintln!("FASE attesa t={:.2}", self.t);
                    }
                }
            }
            Fase::Eccitato => {
                self.b2 = 1.0;
                if self.t_fase > self.attesa {
                    let lam = lambda_draw(self.lambda_nm);
                    let st = Self::sigma_z(lam) / C_SCR;
                    self.emetti(self.t + 3.0 * st, lam);
                    self.fase = Fase::Emissione;
                    self.t_fase = 0.0;
                    self.surfer_in_moto = true;
                    self.phi = 0.0;
                    eprintln!("FASE emissione t={:.2}", self.t);
                }
            }
            Fase::Emissione => {
                let o = self.fout.clone().unwrap();
                let st = Self::sigma_z(o.lam_draw) / C_SCR;
                let s = (self.t - o.t0) / st;
                let r = C_SCR * (self.t - o.t0);
                if r >= 28.0 && r - C_SCR * dt < 28.0 { eprintln!("MOMENTO emissione t={:.2}", self.t); }
                if s < 3.5 {
                    self.b2 = (1.0 - phi_norm(s)).max(1e-9);
                    self.phi += w * dt;
                    self.muovi_surfer(dt);
                } else {
                    if self.surfer_in_moto { self.cont.emessi += 1; }
                    self.b2 = 0.0;
                    self.surfer_in_moto = false;
                }
                if C_SCR * (self.t - o.t0) > Z_END + 10.0 {
                    self.fout = None;
                    self.fase = Fase::Attesa;
                    self.t_fase = 0.0;
                    eprintln!("FASE attesa t={:.2}", self.t);
                }
            }
        }
    }

    /// ψ_e/ψ_g of the NORMALISED guides at configuration x (sign included).
    pub fn rapporto(&self, x: &[f64]) -> f64 {
        self.se.segno(x) * (self.se.ln_guida(x) - self.sg.ln_guida(x)).exp() * self.norma
    }

    /// Bohm velocity (a.u.) of the two electrons for Ψ = a ψ_g + b e^{−iφ} ψ_e (a, b real ≥ 0):
    /// v = Im(∇Ψ/Ψ) = a b r sinφ (∇ln ψ_g − ∇ln ψ_e) / (a² + b² r² + 2 a b r cosφ),  r = ψ_e/ψ_g.
    pub fn velocita_bohm(&self, x: &[f64], a: f64, b: f64, phi: f64, v: &mut [f64; 6]) {
        let r = self.rapporto(x);
        let (mut dg, mut de) = ([0.0; 6], [0.0; 6]);
        self.sg.valuta(x, &mut dg);
        self.se.valuta(x, &mut de);
        let den = (a * a + b * b * r * r + 2.0 * a * b * r * phi.cos()).max(1e-12);
        let c = a * b * r * phi.sin() / den;
        for i in 0..6 { v[i] = c * (dg[i] - de[i]); }
    }

    /// Move the surfer with the Bohm velocity of the drawn superposition for one frame.
    /// One screen second is ω_scr/ΔE atomic units of time (the oscillation is slowed by `rallentamento_osc`).
    fn muovi_surfer(&mut self, dt: f64) {
        let t_au = dt * self.omega_scr() / self.luce.de;
        let n = (t_au / 0.01).ceil().max(1.0) as usize;
        let h = t_au / n as f64;
        let (a, b) = ((1.0 - self.b2).max(0.0).sqrt(), self.b2.sqrt());
        let mut v = [0.0; 6];
        for k in 0..n {
            // φ advances with the physical time inside the frame
            let ph = self.phi - (n - k) as f64 * h * self.luce.de;
            let x = self.surfer;
            self.velocita_bohm(&x, a, b, ph, &mut v);
            // cap the step near the nodes of Ψ (there |v| → ∞): at most 0.05 bohr per electron per sub-step
            for e in 0..2 {
                let vv = (v[3 * e] * v[3 * e] + v[3 * e + 1] * v[3 * e + 1] + v[3 * e + 2] * v[3 * e + 2]).sqrt() * h;
                let f = if vv > 0.05 { 0.05 / vv } else { 1.0 };
                for c in 0..3 { self.surfer[3 * e + c] += f * h * v[3 * e + c]; }
            }
        }
    }

    /// The balls to draw: every electron of every drawn ball (physics coordinates, bohr) and its spin index.
    /// Stationary phases: the ground or the excited cloud as they are. Superposition: candidates from the
    /// ground cloud (probability a²) or the excited one (b²), kept with probability w/2,
    /// w = |aψ_g + b e^{−iφ}ψ_e|² / (a²ψ_g² + b²ψ_e²) = 1 + 2ab r cosφ / (a² + b²r²) ∈ [0, 2].
    pub fn punti(&mut self) -> Vec<([f64; 3], usize)> {
        let mut out = Vec::with_capacity(2 * N_POOL);
        let push = |out: &mut Vec<([f64; 3], usize)>, x: &[f64]| {
            out.push(([x[0], x[1], x[2]], 0));
            out.push(([x[3], x[4], x[5]], 1));
        };
        if self.b2 <= 1e-6 {
            for x in self.g.x.chunks(6) { push(&mut out, x); }
        } else if self.b2 >= 1.0 - 1e-6 {
            for x in self.e.x.chunks(6) { push(&mut out, x); }
        } else {
            let (a, b) = ((1.0 - self.b2).sqrt(), self.b2.sqrt());
            let cphi = self.phi.cos();
            for _ in 0..N_POOL {
                let from_g = self.rng_vis.uni() < a * a;
                let pool = if from_g { &self.g } else { &self.e };
                let i = (self.rng_vis.next_u64() % pool.n() as u64) as usize;
                let x = &pool.x[6 * i..6 * i + 6];
                let r = self.rapporto(x);
                let w = 1.0 + 2.0 * a * b * r * cphi / (a * a + b * b * r * r).max(1e-300);
                if self.rng_vis.uni() * 2.0 < w {
                    let xx: [f64; 6] = x.try_into().unwrap();
                    push(&mut out, &xx);
                }
            }
        }
        self.dipolo = 2.0 * out.iter().map(|(p, _)| p[2]).sum::<f64>() / out.len().max(1) as f64;
        out
    }
}

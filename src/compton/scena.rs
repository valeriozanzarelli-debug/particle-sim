//! One shot of the viewer (scene units: 1 = λ of the incoming photon; screen seconds), and the
//! statistics of many shots (the same surfer rule, run without drawing).

use bevy::prelude::*;

use crate::config::*;
use crate::fisica::*;

/// physics (x = polarization, y, z = beam) → screen (x → up, y → −X, z → Z): right-handed both.
pub fn v3(p: [f64; 3]) -> Vec3 { Vec3::new(-p[1] as f32, p[0] as f32, p[2] as f32) }

/// One shot: where the two surfers are born and what happens when they meet.
#[derive(Clone, Debug)]
pub struct Colpo {
    pub id: u64,
    /// electron surfer at birth (drawn units: inside the drawn electron packet)
    pub e0: Vec3,
    /// photon surfer at birth
    pub g0: Vec3,
    /// relative transverse position b at the meeting (m, real scale)
    pub b_m: [f64; 2],
    /// screen time when the photon surfer reaches the plane of the electron surfer
    pub t_c: f32,
    /// screen time when the centre of the photon packet reaches the centre of the electron packet
    pub t_w: f32,
    pub t_end: f32,
    pub esito: Esito,
    /// azimuth of the event: the drawing is turned by −rot about the beam, so the plane of the
    /// event (beam, photon, electron) is the plane facing the camera (a pure change of viewpoint)
    pub rot: f64,
}

impl Colpo {
    pub fn nuovo(id: u64, rng: &mut Rng, kn: &Kn, pk: &Pacchetti, pol: Pol, solo_diffusi: bool) -> Self {
        let e0 = Vec3::new(rng.gauss() as f32, rng.gauss() as f32, rng.gauss() as f32) * S_E_DRAW;
        let p = pk.p_scatter(kn);
        let b_m = nascita_b(rng, pk, p, solo_diffusi);
        let esito = regola(kn, pk, b_m, pol);
        let rot = match esito { Esito::Diffuso(e) => e.phi, Esito::Passa => 0.0 };
        let turn = |p: [f64; 3]| { let (s, c) = (-rot).sin_cos(); v3([c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]]) };
        // drawn offset of the photon surfer: b in units of the relative width, times the drawn width
        let s = pk.s_rel();
        let bd = turn([b_m[0] / s, b_m[1] / s, 0.0]) * S_GAMMA_DRAW;
        let z0 = Z_START + N_CYCLES * (rng.uni() as f32 - 0.5);
        let g0 = Vec3::new(e0.x + bd.x, e0.y + bd.y, z0);
        let t_c = (e0.z - z0) / C_SCENE;
        let t_w = -Z_START / C_SCENE;
        Colpo { id, e0, g0, b_m, t_c, t_w, t_end: t_w + R_END / C_SCENE, esito, rot }
    }

    /// A physics vector, turned by −rot about the beam, on screen.
    pub fn ruota(&self, p: [f64; 3]) -> Vec3 {
        let (s, c) = (-self.rot).sin_cos();
        v3([c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]])
    }

    /// The direction (θ, φ) on screen.
    pub fn ruota_dir(&self, theta: f64, phi: f64) -> Vec3 { v3(dir(theta, phi - self.rot)) }

    pub fn evento(&self) -> Option<Evento> {
        match self.esito { Esito::Diffuso(e) => Some(e), Esito::Passa => None }
    }

    /// Photon surfer at screen time t.
    pub fn fotone(&self, t: f32) -> Vec3 {
        let straight = self.g0 + Vec3::Z * C_SCENE * t;
        match self.esito {
            Esito::Diffuso(e) if t > self.t_c => {
                let pc = self.g0 + Vec3::Z * C_SCENE * self.t_c;
                pc + self.ruota_dir(e.theta, e.phi) * C_SCENE * (t - self.t_c)
            }
            _ => straight,
        }
    }

    /// Electron surfer at screen time t (at rest until hit: the surfer of a packet at rest).
    pub fn elettrone(&self, t: f32) -> Vec3 {
        match self.esito {
            Esito::Diffuso(e) if t > self.t_c => {
                let n = self.ruota(e.pe).normalize_or_zero();
                self.e0 + n * e.beta_e as f32 * C_SCENE * (t - self.t_c)
            }
            _ => self.e0,
        }
    }
}

/// The statistics of many shots, all scattering (born inside the scattering disk).
#[derive(Clone, Debug)]
pub struct Statistica {
    pub e_kev: f64,
    pub pol: Pol,
    pub eventi: Vec<Evento>,
    /// histogram of cos θ (N_BIN_MU bins on [−1, 1]) and the Klein–Nishina expectation
    pub h_mu: Vec<f64>,
    pub kn_mu: Vec<f64>,
    pub chi2_mu: f64,
    /// histogram of the electron kinetic energy (N_BIN_T bins on [0, 1.05 edge])
    pub h_t: Vec<f64>,
    pub kn_t: Vec<f64>,
    pub chi2_t: f64,
    pub t_max_seen: f64,
    /// largest |p_γ' + p_e − p_in| / p_in over all events
    pub max_resid: f64,
}

pub const N_BIN_MU: usize = 40;
pub const N_BIN_T: usize = 50;

/// Pearson χ² of counts vs expectations (bins with expectation ≥ 5 only); returns (χ², dof).
pub fn chi2(obs: &[f64], exp: &[f64]) -> (f64, usize) {
    let mut c = 0.0;
    let mut n: usize = 0;
    for (o, e) in obs.iter().zip(exp) {
        if *e >= 5.0 { c += (o - e).powi(2) / e; n += 1; }
    }
    (c, n.saturating_sub(1))
}

/// Run `n` shots born inside the scattering disk through the surfer rule.
pub fn statistica(e_kev: f64, pol: Pol, n: usize, seed: u64) -> Statistica {
    let kn = Kn::new(e_kev);
    let pk = crate::pacchetti(&kn);
    let p = pk.p_scatter(&kn);
    let mut rng = Rng::new(seed);
    let mut eventi = Vec::with_capacity(n);
    let mut max_resid: f64 = 0.0;
    while eventi.len() < n {
        let b = nascita_b(&mut rng, &pk, p, true);
        if let Esito::Diffuso(e) = regola(&kn, &pk, b, pol) {
            let r = [e.pg[0] + e.pe[0], e.pg[1] + e.pe[1], e.pg[2] + e.pe[2] - e_kev];
            max_resid = max_resid.max((r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt() / e_kev);
            eventi.push(e);
        }
    }
    let nf = n as f64;
    let mut h_mu = vec![0.0; N_BIN_MU];
    for e in &eventi {
        let i = (((e.theta.cos() + 1.0) * 0.5 * N_BIN_MU as f64) as usize).min(N_BIN_MU - 1);
        h_mu[i] += 1.0;
    }
    let kn_mu: Vec<f64> = (0..N_BIN_MU).map(|i| {
        let a = -1.0 + 2.0 * i as f64 / N_BIN_MU as f64;
        let b = -1.0 + 2.0 * (i + 1) as f64 / N_BIN_MU as f64;
        nf * (kn.cdf_mu(b) - kn.cdf_mu(a))
    }).collect();
    // electron energy: T(μ) is monotone (T grows as μ falls), so the expected counts of a T bin are
    // the KN probability of the matching μ interval
    let edge = compton_edge(e_kev);
    let tmax = 1.05 * edge;
    let mut h_t = vec![0.0; N_BIN_T];
    let mut t_max_seen: f64 = 0.0;
    for e in &eventi {
        t_max_seen = t_max_seen.max(e.t_e);
        let i = ((e.t_e / tmax * N_BIN_T as f64) as usize).min(N_BIN_T - 1);
        h_t[i] += 1.0;
    }
    let mu_of_t = |t: f64| -> f64 {
        if t >= edge { return -1.0; }
        // T = E − E'  →  E' = E − T  →  1 − μ = (E/E' − 1)/k
        let e2 = e_kev - t;
        (1.0 - (e_kev / e2 - 1.0) / kn.k).clamp(-1.0, 1.0)
    };
    let kn_t: Vec<f64> = (0..N_BIN_T).map(|i| {
        let t0 = tmax * i as f64 / N_BIN_T as f64;
        let t1 = tmax * (i + 1) as f64 / N_BIN_T as f64;
        nf * (kn.cdf_mu(mu_of_t(t0)) - kn.cdf_mu(mu_of_t(t1)))
    }).collect();
    let chi2_mu = chi2(&h_mu, &kn_mu).0;
    let chi2_t = chi2(&h_t, &kn_t).0;
    Statistica { e_kev, pol, eventi, h_mu, kn_mu, chi2_mu, h_t, kn_t, chi2_t, t_max_seen, max_resid }
}

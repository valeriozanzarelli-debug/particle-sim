//! What the viewer shows: one system at a time, its balls running LIVE with the rule of dmc.rs.
//! The numbers in the panel come from `molecole --misura` (assets/molecole/risultati.txt).

use crate::dmc::*;
use crate::misura::{a_cuspide, biatomica, Risultati};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Scelta { H1s, H2p, H2s, H3d, H3p, H2piu, H2, He }

impl Scelta {
    pub const TUTTE: [Scelta; 8] = [Scelta::H1s, Scelta::H2p, Scelta::H2s, Scelta::H3d, Scelta::H3p, Scelta::H2piu, Scelta::H2, Scelta::He];
    pub fn nome(self) -> &'static str {
        match self {
            Scelta::H1s => "H 1s", Scelta::H2p => "H 2p", Scelta::H2s => "H 2s", Scelta::H3d => "H 3d",
            Scelta::H3p => "H 3p", Scelta::H2piu => "H2+", Scelta::H2 => "H2", Scelta::He => "He",
        }
    }
    /// Electrons carried by every ball.
    pub fn elettroni(self) -> usize { match self { Scelta::H2 | Scelta::He => 2, _ => 1 } }
    pub fn biatomica(self) -> bool { matches!(self, Scelta::H2piu | Scelta::H2) }
    /// Equilibrium distance of the nuclei (bohr): the computed one if available.
    pub fn r_eq(self, res: &Risultati) -> f64 {
        let ang = crate::costanti::bohr_ang();
        match self {
            Scelta::H2piu => res.f("h2p_mol_re").map_or(2.0, |r| r / ang),
            Scelta::H2 => res.f("h2_mol_re").map_or(1.4, |r| r / ang),
            _ => 0.0,
        }
    }
    /// Harmonic frequency (cm⁻¹) computed for the molecule (for the vibration animation).
    pub fn omega(self, res: &Risultati) -> Option<f64> {
        match self {
            Scelta::H2piu => res.f("h2p_mol_we"),
            Scelta::H2 => res.f("h2_mol_we"),
            _ => None,
        }
    }
    /// The system with its node rules and guide (guide parameters from the measurements, else exact defaults).
    pub fn sistema(self, res: &Risultati, r: f64) -> Sistema {
        let h = |vincoli: Vec<Vincolo>, guida: Guida| Sistema {
            nuclei: vec![Nucleo { pos: [0.0; 3], z: 1.0 }], spin: vec![1], guida, vincoli, campo_z: 0.0,
        };
        let a = |k: &str, d: f64| res.f(k).unwrap_or(d);
        match self {
            Scelta::H1s => h(vec![], Guida::Orbitali { a: a("h_a", 1.0), jastrow: false }),
            Scelta::H2p => h(vec![Vincolo::PianoZ], Guida::Eccitato { ang: Angolare::Z, nodi: [0.0; 2], nn: 0, a: a("h2p_a", 0.5) }),
            Scelta::H3d => h(vec![Vincolo::PianoX, Vincolo::PianoY],
                             Guida::Eccitato { ang: Angolare::XY, nodi: [0.0; 2], nn: 0, a: a("h3d_a", 1.0 / 3.0) }),
            // radial-node states: the whole state (both pockets), node where the balls found it
            Scelta::H2s => {
                let rn = a("h2s_node", 2.0);
                h(vec![], Guida::Eccitato { ang: Angolare::S, nodi: [rn, 0.0], nn: 1, a: 0.5 })
            }
            Scelta::H3p => {
                let rn = a("h3p_node", 6.0);
                h(vec![], Guida::Eccitato { ang: Angolare::Z, nodi: [rn, 0.0], nn: 1, a: 1.0 / 3.0 })
            }
            Scelta::H2piu => biatomica(r, 1.0, 1),
            Scelta::H2 => biatomica(r, 1.0, 2),
            Scelta::He => Sistema {
                nuclei: vec![Nucleo { pos: [0.0; 3], z: 2.0 }], spin: vec![1, -1],
                guida: Guida::Orbitali { a: a_cuspide(2.0, 1e9), jastrow: true }, vincoli: vec![], campo_z: 0.0,
            },
        }
    }
    /// Tick used live (a.u. of time).
    pub fn tau(self) -> f64 {
        match self { Scelta::H3d | Scelta::H3p => 0.04, Scelta::H2p | Scelta::H2s => 0.02, _ => 0.01 }
    }
    /// States whose guide has a radial node: the balls must not cross it (fixed-node rule via the guide sign).
    pub fn nodo_radiale(self) -> Option<&'static str> {
        match self { Scelta::H2s => Some("h2s_node"), Scelta::H3p => Some("h3p_node"), _ => None }
    }
}

/// A live population of balls for one system.
pub fn nuova(s: Scelta, res: &Risultati, n: usize, seed: u64) -> Dmc {
    let r = s.r_eq(res);
    let mut sys = s.sistema(res, r);
    // the full 2s / 3p state: two pockets, each kept by the node of the guide; the viewer starts the balls
    // on both sides by using the node only through the guide (balls do not cross: drift → ∞ at the node,
    // the accept/reject rule refuses any step that changes the sign of the guide)
    if s.nodo_radiale().is_some() { sys.vincoli.clear(); }
    let mut d = Dmc::new(sys, s.tau(), n, seed);
    d.metropolis = true;
    d
}

/// Every electron of every ball, as points to draw (ne points per ball).
pub fn punti(d: &Dmc) -> Vec<([f32; 3], usize)> {
    let ne = d.sys.ne();
    let mut v = Vec::with_capacity(d.n() * ne);
    for e in 0..ne {
        for p in d.elettrone(e) { v.push(([p[0] as f32, p[1] as f32, p[2] as f32], e)); }
    }
    v
}

//! The measurements (no window): `atomo --misura veloce|tutto`.
//! 1. our variational exponents (VMC) for C, C⁺ and C* (2p3s);
//! 2. fixed-node DMC energies of the three states → total energy, ionization energy, 2p→3s excitation;
//! 3. the shell boundary from the radial distribution of the balls;
//! 4. the proton knot.
//! Numbers → assets/atomo/risultati.txt (read by the viewer), summary on stdout.
//! Pre-registered tolerances (docs/ATOMO.md §3): E(C) within 1 % of −37.8450; IP within 3 % of 11.2603 eV;
//! 2p→3s within 3 % of 7.482 eV.

use std::time::Instant;

use crate::carbonio::{confini, istogramma_radiale, misura_dmc, ottimizza, Guida, Param, Stato};
use crate::costanti as k;
use crate::risultati::Risultati;
use crate::skyrme::Profilo;
use crate::ui::{E_ESATTA_C, ECC_NIST_EV, IP_NIST_EV};

pub struct Corsa { pub n_conf: usize, pub giri: usize, pub n0: usize, pub eq: usize, pub k: usize, pub tau: f64 }

pub fn corsa(modo: &str) -> Corsa {
    match modo {
        "veloce" => Corsa { n_conf: 12_000, giri: 2, n0: 800, eq: 1500, k: 3000, tau: 0.005 },
        _ => Corsa { n_conf: 30_000, giri: 3, n0: 2000, eq: 4000, k: 40_000, tau: 0.005 },
    }
}

pub fn main(modo: &str) {
    let t0 = Instant::now();
    let c = corsa(modo);
    let z = crate::mondo::Z;
    let mut r = Risultati::carica();
    eprintln!("== atomo --misura {modo}: n0 {} eq {} k {} tau {}", c.n0, c.eq, c.k, c.tau);
    // 1. variational exponents (start: the bare hydrogen-like guess screened by hand-free coordinate descent)
    let mut p0 = Param::idrogenoidi(z);
    p0.z1 = 5.7;
    p0.zs = 1.6;
    p0.zp = 1.6;
    p0.z3 = 0.35;
    p0.d3 = 0.8;
    let (pc, vc, vce) = ottimizza(z, Stato::Fondamentale, p0, c.n_conf, c.giri, 11, true);
    eprintln!("VMC C   {:.4} ± {:.4}  {:?}  [{:.0} s]", vc, vce, pc, t0.elapsed().as_secs_f64());
    let (pi, vi, vie) = ottimizza(z, Stato::Ione, pc, c.n_conf, c.giri, 12, true);
    eprintln!("VMC C+  {:.4} ± {:.4}  {:?}  [{:.0} s]", vi, vie, pi, t0.elapsed().as_secs_f64());
    let mut pe0 = pc;
    pe0.z3 = 0.35;
    pe0.d3 = 0.8;
    let (pe, ve, vee) = ottimizza(z, Stato::Eccitato, pe0, c.n_conf, c.giri, 13, true);
    eprintln!("VMC C*  {:.4} ± {:.4}  {:?}  [{:.0} s]", ve, vee, pe, t0.elapsed().as_secs_f64());
    r.metti_param("c", &pc);
    r.metti_param("ione", &pi);
    r.metti_param("ecc", &pe);
    r.metti("vmc_c", vc);
    r.metti("vmc_ione", vi);
    r.metti("vmc_ecc", ve);
    // 2. fixed-node DMC
    let mut e = [0.0; 3];
    let mut err = [0.0; 3];
    for (j, (s, p, nome)) in [(Stato::Fondamentale, pc, "c"), (Stato::Ione, pi, "ione"), (Stato::Eccitato, pe, "ecc")].into_iter().enumerate() {
        let (en, er, d) = misura_dmc(Guida { z, stato: s, p }, c.tau, c.n0, c.eq, c.k, 100 + j as u64);
        e[j] = en;
        err[j] = er;
        r.metti(&format!("dmc_{nome}"), en);
        r.metti(&format!("dmc_{nome}_err"), er);
        let h = istogramma_radiale(&d, crate::DR_ISTO, crate::R_ISTO);
        let cf = confini(&h, crate::DR_ISTO);
        for (i, b) in cf.iter().enumerate() { r.metti(&format!("guscio_{nome}_{i}"), *b); }
        eprintln!("DMC {:?}  {:.5} ± {:.5}  gusci {:?} bohr  [{:.0} s]", s, en, er, cf, t0.elapsed().as_secs_f64());
    }
    // 4. the proton
    let pr = Profilo::calcola();
    r.metti("skyrme_a", pr.a);
    r.metti("skyrme_b", pr.numero_barionico());
    r.metti("skyrme_e12pi2", pr.energia_12pi2());
    r.metti("skyrme_rms_fm", pr.raggio_rms_fm());
    r.metti("skyrme_massa_mev", pr.massa_mev());
    let ev = k::hartree_ev();
    let ip = (e[1] - e[0]) * ev;
    let ip_e = (err[0].powi(2) + err[1].powi(2)).sqrt() * ev;
    let ex = (e[2] - e[0]) * ev;
    let ex_e = (err[0].powi(2) + err[2].powi(2)).sqrt() * ev;
    let dev_e = (e[0] - E_ESATTA_C) / E_ESATTA_C.abs();
    let pf = |b: bool| if b { "PASS" } else { "FAIL" };
    let righe = [
        format!("V1 energia C   {:.5} ± {:.5} Eh  vs esatta {:.4}: {:+.3} %  -> {}", e[0], err[0], E_ESATTA_C, 100.0 * dev_e, pf(dev_e.abs() < 0.01)),
        format!("V2 ionizzazione {:.3} ± {:.3} eV vs {:.4}: {:+.2} %  -> {}", ip, ip_e, IP_NIST_EV, 100.0 * (ip - IP_NIST_EV) / IP_NIST_EV, pf(((ip - IP_NIST_EV) / IP_NIST_EV).abs() < 0.03)),
        format!("V3 2p->3s       {:.3} ± {:.3} eV vs {:.3}: {:+.2} %  -> {}  (lambda {:.2} nm)", ex, ex_e, ECC_NIST_EV, 100.0 * (ex - ECC_NIST_EV) / ECC_NIST_EV,
                pf(((ex - ECC_NIST_EV) / ECC_NIST_EV).abs() < 0.03), k::nm_from_hartree(e[2] - e[0])),
        format!("V4 protone B={:.5} E/12pi2={:.5} rms={:.3} fm M={:.1} MeV", pr.numero_barionico(), pr.energia_12pi2(), pr.raggio_rms_fm(), pr.massa_mev()),
        format!("corsa '{}' in {:.0} s", modo, t0.elapsed().as_secs_f64()),
    ];
    for l in &righe { println!("{l}"); }
    if modo != "veloce" {
        r.salva(&format!("atomo --misura {modo} (src/atomo/misura.rs)\n{}", righe.join("\n")));
        eprintln!("scritto {}", crate::risultati::FILE);
    }
}

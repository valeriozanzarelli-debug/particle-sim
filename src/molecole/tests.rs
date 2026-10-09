//! Tests of the molecule level: the rule itself (short runs), the guide algebra, the helpers,
//! and the stored measurements against experiment.

use crate::dmc::*;
use crate::misura::{a_cuspide, biatomica, corri, livelli, Risultati, FILE};
use crate::scena::{nuova, punti, Scelta};

fn h(guida: Guida, vincoli: Vec<Vincolo>) -> Sistema {
    Sistema { nuclei: vec![Nucleo { pos: [0.0; 3], z: 1.0 }], spin: vec![1], guida, vincoli, campo_z: 0.0 }
}

/// The plain rule (no guide) finds the hydrogen ground state: E = −½ hartree within 1 %.
#[test]
fn rule_without_guide_finds_hydrogen() {
    let mut d = Dmc::new(h(Guida::Nessuna, vec![]), 0.01, 3000, 1);
    let r = corri(&mut d, 5.0, 25.0, |_| {});
    assert!((r.e.0 + 0.5).abs() < 0.005, "E = {}", r.e.0);
    assert!((d.n() as f64 - 3000.0).abs() < 600.0, "population not stable: {}", d.n());
}

/// With the exact orbital as guide the local energy is −½ everywhere (zero variance): checks the algebra.
#[test]
fn exact_guide_has_constant_local_energy() {
    let s = h(Guida::Orbitali { a: 1.0, jastrow: false }, vec![]);
    let mut d = vec![0.0; 3];
    for p in [[0.3, -0.2, 0.9], [2.0, 1.0, -3.0], [0.01, 0.0, 0.02]] {
        assert!((s.valuta(&p, &mut d) + 0.5).abs() < 1e-10);
    }
    // 2p_z with a = 1/2 and 3d_xy with a = 1/3 are exact too
    let s = h(Guida::Eccitato { ang: Angolare::Z, nodi: [0.0; 2], nn: 0, a: 0.5 }, vec![]);
    assert!((s.valuta(&[0.4, 0.7, 1.3], &mut d) + 0.125).abs() < 1e-10);
    let s = h(Guida::Eccitato { ang: Angolare::XY, nodi: [0.0; 2], nn: 0, a: 1.0 / 3.0 }, vec![]);
    assert!((s.valuta(&[1.4, 0.7, -2.3], &mut d) + 1.0 / 18.0).abs() < 1e-10);
    // 2s: (r − 2) e^{−r/2}
    let s = h(Guida::Eccitato { ang: Angolare::S, nodi: [2.0, 0.0], nn: 1, a: 0.5 }, vec![]);
    assert!((s.valuta(&[1.4, 0.7, -2.3], &mut d) + 0.125).abs() < 1e-10);
}

/// Drift and local energy of the H2 guide agree with finite differences of ln ψ_g.
#[test]
fn h2_guide_derivatives() {
    let s = biatomica(1.4, 1.0, 2);
    let x = [0.3, -0.1, 0.5, -0.4, 0.2, -0.6];
    let mut g = vec![0.0; 6];
    let el = s.valuta(&x, &mut g);
    let hh = 1e-4;
    let mut lap = 0.0;
    for i in 0..6 {
        let (mut xp, mut xm) = (x, x);
        xp[i] += hh;
        xm[i] -= hh;
        let (lp, l0, lm) = (s.ln_guida(&xp), s.ln_guida(&x), s.ln_guida(&xm));
        let gi = (lp - lm) / (2.0 * hh);
        assert!((gi - g[i]).abs() < 1e-6, "drift {i}: {gi} vs {}", g[i]);
        lap += (lp - 2.0 * l0 + lm) / (hh * hh);
    }
    let g2: f64 = g.iter().map(|v| v * v).sum();
    let el_fd = -0.5 * (lap + g2) + s.potenziale(&x);
    assert!((el - el_fd).abs() < 1e-4, "E_loc {el} vs {el_fd}");
}

/// The node rule: balls in the half space z > 0 find the 2p level (−1/8 hartree) within 1 %.
#[test]
fn node_rule_gives_2p() {
    let s = h(Guida::Eccitato { ang: Angolare::Z, nodi: [0.0; 2], nn: 0, a: 0.49 }, vec![Vincolo::PianoZ]);
    let mut d = Dmc::new(s, 0.02, 2000, 3);
    let r = corri(&mut d, 30.0, 30.0, |_| {});
    assert!((r.e.0 + 0.125).abs() < 0.00125, "E = {}", r.e.0);
}

/// The cusp rule gives a = 1 for one nucleus and a > 1 for two close nuclei.
#[test]
fn cusp_exponent() {
    assert!((a_cuspide(1.0, 1e9) - 1.0).abs() < 1e-12);
    let a = a_cuspide(1.0, 2.0);
    assert!((a - (1.0 + (-2.0 * a).exp())).abs() < 1e-12 && a > 1.1);
}

/// Every ball carries the right number of electrons: the cloud of H2 and He integrates to 2 electrons
/// per ball, H and H2+ to 1 (counted on the points the viewer draws).
#[test]
fn electrons_per_ball() {
    let res = Risultati::default();
    for (s, ne) in [(Scelta::H1s, 1), (Scelta::H2piu, 1), (Scelta::H2, 2), (Scelta::He, 2)] {
        let d = nuova(s, &res, 500, 7);
        let p = punti(&d);
        let per_ball = p.len() as f64 / d.n() as f64;
        assert!((per_ball - ne as f64).abs() / (ne as f64) < 0.005, "{s:?}: {per_ball}");
        assert_eq!(d.sys.ne(), ne);
        let nuc: f64 = d.sys.nuclei.iter().map(|n| n.z).sum();
        assert_eq!(nuc as usize, if s == Scelta::H2piu { 2 } else { ne }, "neutral systems except H2+");
    }
}

/// Levels of −½μ⁻¹ d²/dR² + ½kR²: (n + ½) ω.
#[test]
fn nuclear_levels_of_harmonic_well() {
    let (mu, kk) = (918.0, 0.37);
    let w = (kk / mu as f64).sqrt();
    let lv = livelli(|r| 0.5 * kk * (r - 3.0) * (r - 3.0), mu, 2.0, 4.0, 4000, 3);
    for (n, e) in lv.iter().enumerate() {
        assert!((e / w - (n as f64 + 0.5)).abs() < 1e-3, "level {n}: {}", e / w);
    }
}

/// The stored measurements (assets/molecole/risultati.txt, produced by `--misura`) against experiment.
#[test]
fn stored_measurements_match_experiment() {
    let path = format!("{}/{}", env!("CARGO_MANIFEST_DIR"), FILE);
    let r = Risultati::load(&path);
    let chk = |k: &str, exp: f64, tol_rel: f64| {
        let v = r.f(k).unwrap_or_else(|| panic!("{k} missing in {path}"));
        assert!(((v - exp) / exp).abs() < tol_rel, "{k} = {v}, experiment {exp}");
    };
    chk("h_E_guided", -0.5, 1e-4);
    chk("line_lyman_alpha_nm", 121.567, 1e-4);
    chk("line_balmer_alpha_nm", 656.461, 5e-4);
    chk("h_alpha", 4.5, 0.06);
    chk("h2_mol_re", 0.74144, 0.02);
    chk("h2_mol_d0", 4.47813, 0.02);
    chk("h2_mol_we", 4401.21, 0.03);
    chk("he_E_ev_mu", -79.0052, 0.002);
    // helium and light (`--misura elio1p`, `--misura eliocampo`)
    chk("he1p_E", -2.123843, 1e-3);           // fixed-node 1s2p ¹P vs exact nonrelativistic (error stated: +5.5e-5 Eh)
    chk("he1p_lambda_nm", 58.4334, 0.01);     // resonance line 1¹S → 2¹P, NIST, within 1 %
    chk("he1p_de_ev", 21.2180, 0.01);
    chk("he1p_f", 0.2762, 0.15);              // oscillator strength, within 15 %
    chk("he1p_A", 1.799e9, 0.15);             // Einstein A, NIST, within 15 %
    chk("he_alpha", 1.3832, 0.10);            // static polarizability, exact 1.3832 a0³
}

/// Drift and local energy of the helium 1s2p ¹P guide agree with finite differences of ln ψ_g; the guide is
/// symmetric under exchange of the electrons (singlet) and odd under z → −z of both (¹P, M = 0).
#[test]
fn he_1s2p_guide_derivatives_and_symmetry() {
    let s = crate::elio::sys_1p(0.5);
    let x = [0.3, -0.1, 0.5, -0.9, 1.2, -1.6];
    let mut g = vec![0.0; 6];
    let el = s.valuta(&x, &mut g);
    let hh = 1e-4;
    let mut lap = 0.0;
    for i in 0..6 {
        let (mut xp, mut xm) = (x, x);
        xp[i] += hh;
        xm[i] -= hh;
        let (lp, l0, lm) = (s.ln_guida(&xp), s.ln_guida(&x), s.ln_guida(&xm));
        assert!(((lp - lm) / (2.0 * hh) - g[i]).abs() < 1e-6, "drift {i}");
        lap += (lp - 2.0 * l0 + lm) / (hh * hh);
    }
    let g2: f64 = g.iter().map(|v| v * v).sum();
    let el_fd = -0.5 * (lap + g2) + s.potenziale(&x);
    assert!((el - el_fd).abs() < 1e-4, "E_loc {el} vs {el_fd}");
    let p = crate::elio::psi(&s, &x);
    let swap = [x[3], x[4], x[5], x[0], x[1], x[2]];
    let mirror = [x[0], x[1], -x[2], x[3], x[4], -x[5]];
    assert!((crate::elio::psi(&s, &swap) - p).abs() < 1e-12 * p.abs().max(1e-300));
    assert!((crate::elio::psi(&s, &mirror) + p).abs() < 1e-12 * p.abs().max(1e-300));
}

/// The screened helium guide: drift and local energy agree with finite differences; the cusp is exact
/// (no 1/r divergence of the local energy at the nucleus).
#[test]
fn he_screened_guide_derivatives() {
    let s = crate::elio::sys_schermato(1.7, 0.0);
    let x = [0.3, -0.1, 0.5, -0.9, 1.2, -1.6];
    let mut g = vec![0.0; 6];
    let el = s.valuta(&x, &mut g);
    let hh = 1e-4;
    let mut lap = 0.0;
    for i in 0..6 {
        let (mut xp, mut xm) = (x, x);
        xp[i] += hh;
        xm[i] -= hh;
        let (lp, l0, lm) = (s.ln_guida(&xp), s.ln_guida(&x), s.ln_guida(&xm));
        assert!(((lp - lm) / (2.0 * hh) - g[i]).abs() < 1e-6, "drift {i}");
        lap += (lp - 2.0 * l0 + lm) / (hh * hh);
    }
    let g2: f64 = g.iter().map(|v| v * v).sum();
    assert!((el - (-0.5 * (lap + g2) + s.potenziale(&x))).abs() < 1e-4);
    let (e1, e2) = (s.valuta(&[1e-6, 0.0, 0.0, 0.0, 1.0, 0.0], &mut g), s.valuta(&[1e-4, 0.0, 0.0, 0.0, 1.0, 0.0], &mut g));
    assert!((e1 - e2).abs() < 1e-2, "cusp: {e1} vs {e2}");
}

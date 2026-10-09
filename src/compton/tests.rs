//! Tests of the Compton binary: kinematics to 1e-9, Klein–Nishina statistics of the surfer rule,
//! the Compton edge, conservation per event, polarization, the unconditioned scattered fraction.

use std::f64::consts::PI;

use crate::fisica::*;
use crate::misura::{errori_cinematica, test_frazione, test_polarizzazione, ENERGIE};
use crate::scena::{chi2, statistica, Colpo};

/// χ² threshold: dof + 4√(2 dof) (≈ p < 10⁻⁴ for a correct sampler; fixed before running).
fn soglia(dof: usize) -> f64 { dof as f64 + 4.0 * (2.0 * dof as f64).sqrt() }

#[test]
fn cinematica_esatta() {
    let m = errori_cinematica();
    assert!(m[0] < 1e-9, "Compton formula vs momenta: {:e}", m[0]);
    assert!(m[1] < 1e-9, "energy conservation: {:e}", m[1]);
    assert!(m[2] < 1e-9, "momentum conservation: {:e}", m[2]);
    assert!(m[3] < 1e-9, "electron angle: {:e}", m[3]);
}

#[test]
fn costanti() {
    // h/mₑc from hc and mₑc²; σ_T; Mo Kα wavelength
    assert!((HC_KEV_PM / ME_C2_KEV - LAMBDA_C_PM).abs() < 1e-9);
    assert!((sigma_thomson() - 6.652_458_7321e-29).abs() < 1e-37);
    assert!((lambda_pm(17.48) - 70.93).abs() < 0.005);
}

#[test]
fn sezione_totale() {
    for &e in &[2.48e-3, 0.1, 1.0, 17.48, 100.0, 662.0, 1332.5, 1.0e4] {
        let kn = Kn::new(e);
        let f = kn_total_formula(kn.k);
        assert!((kn.sigma_re2 - f).abs() / f < 1e-9, "E={e}: {} vs {}", kn.sigma_re2, f);
    }
    // Thomson limit
    let kn = Kn::new(2.48e-3);
    assert!((kn.sigma_m2() / sigma_thomson() - 1.0).abs() < 2e-5);
}

#[test]
fn inversione_tabella() {
    let kn = Kn::new(662.0);
    for i in 0..=200 {
        let v = i as f64 / 200.0;
        let mu = kn.mu_of_rank(v);
        assert!((kn.cdf_mu(mu) - v).abs() < 1e-12, "v={v} mu={mu} cdf={}", kn.cdf_mu(mu));
    }
}

#[test]
fn bordo_compton() {
    assert!((compton_edge(662.0) - 477.65).abs() < 0.005, "{}", compton_edge(662.0));
    let s = statistica(662.0, Pol::NonPolarizzata, 20_000, 3);
    let edge = compton_edge(662.0);
    assert!(s.t_max_seen <= edge * (1.0 + 1e-12));
    assert!(s.t_max_seen > 0.99 * edge, "highest T {} keV", s.t_max_seen);
    let (c, d) = chi2(&s.h_t, &s.kn_t);
    assert!(c < soglia(d), "χ²(T) = {c} / {d}");
}

#[test]
fn distribuzione_angolare_klein_nishina() {
    for (i, (nome, e)) in ENERGIE.iter().enumerate() {
        let s = statistica(*e, Pol::NonPolarizzata, 20_000, 100 + i as u64);
        let (c, d) = chi2(&s.h_mu, &s.kn_mu);
        println!("{nome} {e} keV: χ² = {c:.1} / {d}");
        assert!(c < soglia(d), "{nome}: χ² = {c} / {d}");
        assert!(s.max_resid < 1e-9, "{nome}: momentum residual {}", s.max_resid);
    }
}

#[test]
fn polarizzazione() {
    let (c, d, r) = test_polarizzazione(17.48, 40_000, 7);
    assert!(c < soglia(d), "χ²(φ) = {c} / {d}");
    // near θ = 90° almost nothing goes along the electric field
    assert!(r < 0.25, "ratio φ≈0 / φ≈90° = {r}");
}

#[test]
fn frazione_senza_condizionamento() {
    let (f, p) = test_frazione(17.48, 2.0e-14, 200_000, 11);
    let sd = (p * (1.0 - p) / 200_000.0).sqrt();
    assert!((f - p).abs() < 4.0 * sd, "fraction {f} vs p {p} (sd {sd})");
}

#[test]
fn colpo_coerente() {
    // the surfers of a scattering shot: photon on its line at c, electron at β c along p_e
    let kn = Kn::new(662.0);
    let pk = crate::pacchetti(&kn);
    let mut rng = Rng::new(5);
    for id in 0..200 {
        let c = Colpo::nuovo(id, &mut rng, &kn, &pk, Pol::NonPolarizzata, true);
        let e = c.evento().expect("conditioned shots scatter");
        let t1 = c.t_c + 2.0;
        let dg = (c.fotone(t1) - c.fotone(c.t_c)).length();
        let de = (c.elettrone(t1) - c.elettrone(c.t_c)).length();
        assert!((dg - 2.0 * crate::config::C_SCENE).abs() < 1e-3);
        assert!((de - 2.0 * crate::config::C_SCENE * e.beta_e as f32).abs() < 1e-3);
        assert!(e.theta >= 0.0 && e.theta <= PI);
    }
}


//! THE LIVE TESTBED (lab/out/banco_prereg.json): the rule's own wave packet against the live zig-zag atom, true charge,
//! amplitude from N photons of energy. Nothing is drawn or piloted: it reports what happens.
//!   cargo test --bin atomo banco_corsa -- --ignored --nocapture

use crate::accoppiato::{self as acc, Esperimento};

fn scrivi(righe: &[serde_json::Value], e_vero: f64, n_pi2: f64) {
    let v = serde_json::json!({ "e_vero": e_vero, "N_pi2": n_pi2, "corse": righe });
    std::fs::write("lab/out/banco_risultati.json", serde_json::to_string_pretty(&v).unwrap()).unwrap();
}

#[test]
#[ignore]
fn banco_corsa() {
    let t0 = std::time::Instant::now();
    let mut x = Esperimento::nuovo(600, None).con_diffusa();
    x.k_acc = 1.0;
    let e_vero = (4.0 * std::f64::consts::PI * acc::alpha()).sqrt();
    assert!((acc::e_eff() * x.k_acc / acc::K_ACC / e_vero - 1.0).abs() < 1e-12, "B2");
    let w = x.omega_21();
    println!("preparato in {:.0} s, omega_21 {w:.5}", t0.elapsed().as_secs_f64());
    // N_pi2: the N whose peak e·E is the π/2 pulse's; e·E ∝ √N
    let mut prova = Esperimento::nuovo(30, None);
    prova.k_acc = 1.0;
    let ee1 = prova.spara_fotoni(w, 1.0);
    drop(prova);
    let n_pi2 = (crate::zzscena::E0_PI2 / ee1).powi(2);
    println!("e·E di picco con 1 fotone {ee1:.4e}; N_pi2 = {n_pi2:.0}");
    let mut righe = vec![];
    for n in [0.0, 1.0, 1e2, 1e4, 1e5, n_pi2] {
        let t1 = std::time::Instant::now();
        x.ricomincia();
        for _ in 0..200 { x.passo(); }
        let p0 = x.serie_p.last().map(|v| v.1).unwrap_or(0.0);
        let (e_rapp, ee) = if n > 0.0 {
            let ee = x.spara_fotoni(w, n);
            (x.luce.energia() / (n * w), ee)
        } else { (0.0, 0.0) };
        let t_spara = x.at.t();
        for _ in 0..1600 { x.passo(); }
        let t1f = x.at.t();
        let coda: Vec<f64> = x.serie_p.iter().filter(|v| v.0 > t1f - 100.0).map(|v| v.1).collect();
        let dp = coda.iter().sum::<f64>() / coda.len().max(1) as f64 - p0;
        let d = acc::ampiezza_a(&x.serie_d, w, t1f - 400.0, t1f);
        let r = serde_json::json!({ "N": n, "p0": p0, "dp": dp, "D": d, "e_E_picco": ee, "energia_su_N_omega": e_rapp, "e_diffusa": x.luce.energia_diffusa(), "t_spara": t_spara, "sec": t1.elapsed().as_secs_f64() });
        println!("{r}");
        righe.push(r);
        scrivi(&righe, e_vero, n_pi2);
    }
    println!("fatto in {:.0} s", t0.elapsed().as_secs_f64());
}

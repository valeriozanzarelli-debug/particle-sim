//! THE LIVE TESTBED (lab/out/banco_prereg.json): the rule's own wave packet against the live zig-zag atom, true charge,
//! amplitude from N photons of energy. Nothing is drawn or piloted: it reports what happens.
//!   cargo test --bin atomo banco_corsa -- --ignored --nocapture

use crate::accoppiato::{self as acc, Esperimento};
use crate::cplx::C;

fn scrivi(righe: &[serde_json::Value], e_vero: f64, n_pi2: f64) {
    let v = serde_json::json!({ "e_vero": e_vero, "N_pi2": n_pi2, "corse": righe });
    std::fs::write("lab/out/banco_risultati.json", serde_json::to_string_pretty(&v).unwrap()).unwrap();
}

fn sopra(a: &[[C; 4]], b: &[[C; 4]]) -> C {
    use rayon::prelude::*;
    a.par_iter().zip(b.par_iter()).map(|(x, y)| { let mut s = C::ZERO; for q in 0..4 { s += x[q].conj() * y[q]; } s }).reduce(|| C::ZERO, |u, v| u + v)
}

/// ∫ E(t) e^{iωt} dt of a series sampled every light tick: its modulus.
fn trasformata(serie: &[(f64, f64)], w: f64) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    for (t, e) in serie { re += e * (w * t).cos(); im += e * (w * t).sin(); }
    (re * re + im * im).sqrt() * acc::B_LUCE
}

/// Phase 2 (lab/out/banco_livelli_prereg.json): where the absorbed light goes, level by level.
#[test]
#[ignore]
fn banco_livelli() {
    use crate::orbitali::{dirac_esatto, radiale};
    let t0 = std::time::Instant::now();
    let mut x = Esperimento::nuovo(600, None).con_diffusa();
    x.k_acc = 1.0;
    let e = (4.0 * std::f64::consts::PI * acc::alpha()).sqrt();
    let w21 = x.omega_21();
    let w1s = x.omega1s;
    let lv = x.at.livelli_di_dipolo(&x.phi1s, 2, x.omega1s + 0.05, x.omega1s + 0.095, 0.0015, 6000, 1e-3);
    println!("omega1s {:.6}, omega_21 {w21:.6}, {} livelli in {:.0} s", x.omega1s, lv.len(), t0.elapsed().as_secs_f64());
    // ⟨z²⟩ − ⟨z⟩² of the rule's 1s
    let (mut z1, mut z2) = (0.0, 0.0);
    for (s, p) in x.phi1s.iter().enumerate() {
        let c = crate::zz::sito(s, None);
        let z = crate::zz::x_of(c[2]);
        let r: f64 = p.iter().map(|q| q.norm2()).sum();
        z1 += r * z;
        z2 += r * z * z;
    }
    let var_z = z2 - z1 * z1;
    // Dirac from the radial rule on the zoom mesh
    let (r1, r32, r12) = (radiale(-1, 1, 0), radiale(-2, 2, 1), radiale(1, 2, 1));
    let raggio = |a: &crate::orbitali::Radiale, b: &crate::orbitali::Radiale| -> f64 {
        a.gf.iter().zip(&b.gf).zip(a.r.iter().zip(&a.h)).map(|((p, q), (r, h))| (p[0] * q[0] + p[1] * q[1]) * r * h).sum()
    };
    let (d32, d12) = (2f64.sqrt() / 3.0 * raggio(&r1, &r32), raggio(&r1, &r12) / 3.0);
    let (wd32, wd12) = (dirac_esatto(2, -2) - dirac_esatto(1, -1), dirac_esatto(2, 1) - dirac_esatto(1, -1));
    println!("Dirac: d3/2 {d32:.4}, d1/2 {d12:.4}, somma d^2 {:.4}, omega {wd32:.5} {wd12:.5}; <z^2>-<z>^2 regola {var_z:.4}, Dirac {:.4}", d32 * d32 + d12 * d12, raggio_z2(&r1));
    // the runs
    let h = acc::N_LUCE / 2;
    let centro = ((h * acc::N_LUCE + h) * acc::N_LUCE) + h;
    let mut corsa = |n: f64| {
        x.ricomincia();
        for _ in 0..200 { x.passo(); }
        let p0 = x.serie_p.last().map(|v| v.1).unwrap_or(0.0);
        if n > 0.0 { x.spara_fotoni(w21, n); }
        let (mut eg, mut ep) = (vec![], vec![]);
        for k in 0..1600 {
            x.passo();
            if k % acc::SUB == 0 {
                let t = x.at.t();
                eg.push((t, x.luce.campo_libero_atomo()[2]));
                ep.push((t, x.luce.f.f[centro][2].re as f64));
            }
        }
        let tf = x.at.t();
        let coda: Vec<f64> = x.serie_p.iter().filter(|v| v.0 > tf - 100.0).map(|v| v.1).collect();
        let dp = coda.iter().sum::<f64>() / coda.len().max(1) as f64 - p0;
        (x.at.psi.clone(), dp, eg, ep)
    };
    let (psi0, dp0, _, _) = corsa(0.0);
    let mut righe = vec![];
    let mut somme = vec![];
    for n in [1.0, 100.0] {
        let (psi, dp, eg, ep) = corsa(n);
        let delta: Vec<[C; 4]> = psi.iter().zip(&psi0).map(|(a, b)| [0, 1, 2, 3].map(|q| a[q] - b[q])).collect();
        let p_tot: f64 = delta.iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum();
        let (mut somma, mut prev) = (0.0, 0.0);
        let mut livelli = vec![];
        for l in &lv {
            let wj = l.lv.omega - w1s;
            let pj = sopra(&l.lv.v, &delta).norm2();
            let dd = l.d[2].norm2();
            let pr = (e * trasformata(&eg, wj)).powi(2) * dd;
            somma += pj;
            prev += pr;
            livelli.push(serde_json::json!({ "omega": l.lv.omega, "resid": l.lv.resid, "d2": dd, "P": pj, "previsione": pr }));
        }
        let dirac = (e * trasformata(&eg, wd32)).powi(2) * d32 * d32 + (e * trasformata(&eg, wd12)).powi(2) * d12 * d12;
        let chiusura = (e * trasformata(&eg, w21)).powi(2) * var_z;
        let r = serde_json::json!({ "N": n, "dp_meno_controllo": dp - dp0, "P_tot": p_tot, "somma_P": somma, "previsione_regola": prev, "previsione_Dirac": dirac, "chiusura": chiusura,
            "L2": somma / prev, "L3": somma / dirac, "L4": p_tot / chiusura, "L5": (dp - dp0) / p_tot, "L6_resto": 1.0 - somma / p_tot,
            "E_gauss_su_punto": trasformata(&eg, w21) / trasformata(&ep, w21), "livelli": livelli });
        println!("{r}");
        somme.push(somma);
        righe.push(r);
    }
    println!("L1 = {:.3}; fatto in {:.0} s", somme[1] / somme[0], t0.elapsed().as_secs_f64());
    let v = serde_json::json!({ "omega_1s": x.omega1s, "omega_21": w21, "var_z_regola": var_z, "var_z_Dirac": raggio_z2(&r1), "d32_Dirac": d32, "d12_Dirac": d12, "L1": somme[1] / somme[0], "corse": righe });
    std::fs::write("lab/out/banco_livelli_risultati.json", serde_json::to_string_pretty(&v).unwrap()).unwrap();
}

/// ⟨z²⟩ = ⟨r²⟩/3 of a radial s state (r·g, r·f on the mesh).
fn raggio_z2(a: &crate::orbitali::Radiale) -> f64 {
    a.gf.iter().zip(a.r.iter().zip(&a.h)).map(|(p, (r, h))| (p[0] * p[0] + p[1] * p[1]) * r * r * h).sum::<f64>() / 3.0
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

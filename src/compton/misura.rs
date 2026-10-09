//! Measurements without the window: the validation table (`--misura`) and the event files for the
//! independent Python check (`--eventi N --out file.csv`, read by lab/compton_check.py).

use std::f64::consts::PI;
use std::io::Write;

use crate::fisica::*;
use crate::scena::{chi2, statistica};

/// Largest relative error of the kinematics over a grid of energies and angles:
/// (λ' from the momenta vs Compton's formula, energy conservation, momentum conservation,
///  electron angle vs tan φ = cot(θ/2)/(1+k)).
pub fn errori_cinematica() -> [f64; 4] {
    let mut m = [0.0f64; 4];
    for &e in &[2.5e-3, 1.0, 17.48, 100.0, 662.0, 1173.2, 1332.5, 1.0e4] {
        for i in 0..=720 {
            let th = PI * i as f64 / 720.0;
            for &ph in &[0.0, 0.7, 2.0, 4.5] {
                let ev = cinematica(e, th, ph);
                let lam = lambda_pm(e);
                // λ' from the photon momentum |p'| = h/λ'
                let p2 = (ev.pg[0].powi(2) + ev.pg[1].powi(2) + ev.pg[2].powi(2)).sqrt();
                let lam2 = HC_KEV_PM / p2;
                let shift_formula = LAMBDA_C_PM * (1.0 - th.cos());
                m[0] = m[0].max(((lam2 - lam) - shift_formula).abs() / lam2);
                // energy: E_in + mc² = E' + sqrt(p_e² + m²)
                let pe2 = ev.pe[0].powi(2) + ev.pe[1].powi(2) + ev.pe[2].powi(2);
                let ee = (pe2 + ME_C2_KEV * ME_C2_KEV).sqrt();
                m[1] = m[1].max(((e + ME_C2_KEV) - (p2 + ee)).abs() / (e + ME_C2_KEV));
                // momentum per event
                let r = [ev.pg[0] + ev.pe[0], ev.pg[1] + ev.pe[1], ev.pg[2] + ev.pe[2] - e];
                m[2] = m[2].max((r[0].powi(2) + r[1].powi(2) + r[2].powi(2)).sqrt() / e);
                if th > 1e-6 && th < PI - 1e-6 {
                    m[3] = m[3].max((ev.theta_e - electron_angle_formula(e, th)).abs());
                }
            }
        }
    }
    m
}

/// Azimuth test with linear polarization: events with 80° < θ < 100°, φ in 12 bins vs Klein–Nishina.
/// Returns (χ², dof, counts near φ = 0 / counts near φ = 90°).
pub fn test_polarizzazione(e_kev: f64, n: usize, seed: u64) -> (f64, usize, f64) {
    let kn = Kn::new(e_kev);
    let pk = crate::pacchetti(&kn);
    let p = pk.p_scatter(&kn);
    let mut rng = Rng::new(seed);
    let nb = 12;
    let mut h = vec![0.0; nb];
    let (a, b) = (80f64.to_radians(), 100f64.to_radians());
    let mut tot = 0usize;
    let mut done = 0usize;
    while done < n {
        let bb = nascita_b(&mut rng, &pk, p, true);
        if let Esito::Diffuso(e) = regola(&kn, &pk, bb, Pol::LineareX) {
            done += 1;
            if e.theta > a && e.theta < b {
                let i = ((e.phi.rem_euclid(2.0 * PI)) / (2.0 * PI) * nb as f64) as usize;
                h[i.min(nb - 1)] += 1.0;
                tot += 1;
            }
        }
    }
    // expected: ∫∫ dσ/dΩ sinθ dθ dφ over the band and bin (2-D Gauss–Legendre via fine midpoint grid)
    let mut w = vec![0.0; nb];
    let (nt, np) = (200, 120);
    for i in 0..nb {
        let f0 = 2.0 * PI * i as f64 / nb as f64;
        let f1 = 2.0 * PI * (i + 1) as f64 / nb as f64;
        let mut s = 0.0;
        for it in 0..nt {
            let th = a + (b - a) * (it as f64 + 0.5) / nt as f64;
            for ip in 0..np {
                let ph = f0 + (f1 - f0) * (ip as f64 + 0.5) / np as f64;
                s += kn_domega(kn.k, th, ph, Pol::LineareX) * th.sin();
            }
        }
        w[i] = s;
    }
    let ws: f64 = w.iter().sum();
    let exp: Vec<f64> = w.iter().map(|x| x / ws * tot as f64).collect();
    let (c, dof) = chi2(&h, &exp);
    if std::env::var("DEBUG_PHI").is_ok() { println!("h={h:?}
exp={exp:?}"); }
    let near0 = h[0] + h[nb - 1] + h[nb / 2 - 1] + h[nb / 2];
    let near90 = h[nb / 4 - 1] + h[nb / 4] + h[3 * nb / 4 - 1] + h[3 * nb / 4];
    (c, dof, near0 / near90)
}

/// Fraction of shots that scatter when the surfers are born with the true |ψ|² (no conditioning),
/// with an artificially narrow relative width so that p is large enough to count.
pub fn test_frazione(e_kev: f64, s_rel_m: f64, n: usize, seed: u64) -> (f64, f64) {
    let kn = Kn::new(e_kev);
    let pk = Pacchetti { s_gamma_m: s_rel_m, s_e_m: 0.0 };
    let p = pk.p_scatter(&kn);
    let mut rng = Rng::new(seed);
    let mut k = 0usize;
    for _ in 0..n {
        let b = nascita_b(&mut rng, &pk, p, false);
        if matches!(regola(&kn, &pk, b, Pol::NonPolarizzata), Esito::Diffuso(_)) { k += 1; }
    }
    (k as f64 / n as f64, p)
}

pub const ENERGIE: [(&str, f64); 3] = [("Mo Kα", 17.48), ("Cs-137", 662.0), ("Co-60", 1332.5)];

pub fn tabella() {
    println!("Compton — validazione (seed fissi)\n");
    let m = errori_cinematica();
    println!("1. cinematica (8 energie × 721 angoli × 4 azimut):");
    println!("   |Δλ − λc(1−cosθ)|/λ'   max = {:.2e}", m[0]);
    println!("   energia conservata      max = {:.2e}", m[1]);
    println!("   impulso conservato      max = {:.2e}", m[2]);
    println!("   angolo elettrone vs tanφ = cot(θ/2)/(1+k)   max = {:.2e} rad", m[3]);
    println!("\n   sezione d'urto totale (r_e²): tabella numerica vs formula chiusa di Klein–Nishina");
    for &e in &[2.48e-3, 1.0, 17.48, 662.0, 1173.2, 1332.5] {
        let kn = Kn::new(e);
        let f = kn_total_formula(kn.k);
        println!("   E = {:>9.4} keV  σ = {:.10} (tab)  {:.10} (formula)  σ/σ_T = {:.6}  rel.diff {:.1e}",
                 e, kn.sigma_re2, f, kn.sigma_re2 / (8.0 * PI / 3.0), (kn.sigma_re2 - f).abs() / f);
    }
    println!("\n2.-3. distribuzione angolare e spettro dell'elettrone (20 000 eventi, regola dei surfisti):");
    for (i, (nome, e)) in ENERGIE.iter().enumerate() {
        let s = statistica(*e, Pol::NonPolarizzata, 20_000, 100 + i as u64);
        let (c1, d1) = chi2(&s.h_mu, &s.kn_mu);
        let (c2, d2) = chi2(&s.h_t, &s.kn_t);
        println!("   {nome:>7} {e:>8.2} keV: χ²(cosθ) = {c1:.1}/{d1}   χ²(T_e) = {c2:.1}/{d2}   T_e max visto {:.3} keV, \
                  bordo {:.3} keV   max|Δp|/p = {:.1e}", s.t_max_seen, compton_edge(*e), s.max_resid);
    }
    println!("   bordo Compton a 662 keV: {:.3} keV (atteso 477.65)", compton_edge(662.0));
    let (c, d, r) = test_polarizzazione(17.48, 40_000, 7);
    println!("\n   polarizzazione lineare, 80°<θ<100°: χ²(φ) = {c:.1}/{d}, conteggi vicino a φ=0 / φ=90° = {r:.3}");
    let (f, p) = test_frazione(17.48, 2.0e-14, 200_000, 11);
    println!("   frazione diffusa senza condizionamento (s = 20 fm, solo test): {f:.5} vs p = {p:.5}");
    println!("\n4. Compton 1923 (Mo Kα su grafite):");
    for &(ang, lam, dl) in COMPTON_1923.iter() {
        let ours = lam_out(lam, ang.to_radians()) - lam;
        println!("   θ = {ang:>5.1}°: Δλ misurato {dl:.2} pm   nostro {ours:.4} pm   diff {:+.1} %", 100.0 * (dl - ours) / ours);
    }
    for &ang in &[45.0f64, 90.0, 135.0] {
        println!("   θ = {ang:>5.1}°: λ' = {:.3} pm (λ = 70.93 pm)", lam_out(70.93, ang.to_radians()));
    }
    println!("\n   probabilità di diffusione per fotone (pacchetto largo 10 λ, elettrone 1 nm):");
    for &e in &[2.48e-3, 17.48, 662.0, 1332.5] {
        let kn = Kn::new(e);
        let pk = crate::pacchetti(&kn);
        println!("   E = {e:>9.4} keV: σ = {:.4e} m², p = {:.3e}, raggio del disco √(σ/π) = {:.3e} m",
                 kn.sigma_m2(), pk.p_scatter(&kn), (kn.sigma_m2() / PI).sqrt());
    }
}

/// Write `n` scattering events (surfer rule, born inside the scattering disk) to a CSV file.
pub fn csv(e_kev: f64, pol: Pol, n: usize, seed: u64, out: &str) {
    let kn = Kn::new(e_kev);
    let pk = crate::pacchetti(&kn);
    let p = pk.p_scatter(&kn);
    let mut rng = Rng::new(seed);
    if let Some(dir) = std::path::Path::new(out).parent() { let _ = std::fs::create_dir_all(dir); }
    let f = std::fs::File::create(out).expect("cannot create the output file");
    let mut w = std::io::BufWriter::new(f);
    writeln!(w, "# e_kev={e_kev:.17e} lam_pm={:.17e} seed={seed} pol={:?} p_scatter={p:.6e} sigma_m2={:.10e}",
             kn.lam_pm, pol, kn.sigma_m2()).unwrap();
    writeln!(w, "theta,phi,e_out,lam_out,pgx,pgy,pgz,pex,pey,pez,t_e,theta_e,beta_e,bx,by").unwrap();
    let mut k = 0;
    while k < n {
        let b = nascita_b(&mut rng, &pk, p, true);
        if let Esito::Diffuso(e) = regola(&kn, &pk, b, pol) {
            writeln!(w, "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.6e},{:.6e}",
                     e.theta, e.phi, e.e_out, e.lam_out, e.pg[0], e.pg[1], e.pg[2], e.pe[0], e.pe[1], e.pe[2],
                     e.t_e, e.theta_e, e.beta_e, b[0], b[1]).unwrap();
            k += 1;
        }
    }
    println!("{n} eventi -> {out}");
}

//! Tests of the collision: the cross-section per photon, the dipole pattern, the drawn superposition,
//! and the default event (one absorption and one emission within 15 s of opening).

use crate::elio::LuceElio;
use crate::fisica::{dir_dipolo, Fase, Scontro};
use crate::misura::{Risultati, FILE};
use std::f64::consts::PI;

fn risultati() -> Risultati { Risultati::load(&format!("{}/{}", env!("CARGO_MANIFEST_DIR"), FILE)) }

/// At the line centre the cross-section is the resonant one, 3λ²/2π, and P per photon is 3/π².
#[test]
fn resonant_cross_section_at_line_centre() {
    let l = LuceElio::da(&risultati()).expect("helium light in the results file");
    let lam = l.lambda_nm;
    let s0 = 3.0 * (lam * 1e-9).powi(2) / (2.0 * PI);
    assert!((l.sigma(lam) / s0 - 1.0).abs() < 1e-6, "σ = {} vs 3λ²/2π = {s0}", l.sigma(lam));
    assert!((l.prob(lam) - 3.0 / (PI * PI)).abs() < 1e-6);
    // half width: at Δω = Γ/2 the Lorentzian is half
    let w = LuceElio::omega_s(lam) + 0.5 * l.a_s;
    let lam_half = 2.0 * PI * crate::costanti::C_L / w * 1e9;
    assert!((l.sigma(lam_half) / s0 - 0.5).abs() < 1e-3, "{}", l.sigma(lam_half) / s0);
    // far off resonance: Rayleigh, tiny
    let p500 = l.prob(500.0);
    assert!(p500 > 0.0 && p500 < 1e-15, "P(500 nm) = {p500}");
    assert!(!l.risonante(500.0) && l.risonante(lam));
}

/// Directions of the emitted photon: ⟨cos²θ⟩ = 1/5 for the dipole pattern sin²θ.
#[test]
fn dipole_pattern() {
    let mut rng = crate::dmc::Rng::new(5);
    let n = 200_000;
    let c2: f64 = (0..n).map(|_| dir_dipolo(&mut rng)[1].powi(2)).sum::<f64>() / n as f64;
    assert!((c2 - 0.2).abs() < 0.003, "<cos²θ> = {c2}");
}

/// The drawn superposition (a = b = 1/√2) carries a dipole that follows cos φ: the cloud sloshes.
#[test]
fn superposition_sloshes() {
    let mut s = Scontro::new(&risultati(), 9);
    s.b2 = 0.5;
    let zm = |s: &mut Scontro, phi: f64| {
        s.phi = phi;
        let (mut z, mut n) = (0.0, 0.0);
        for _ in 0..20 {
            for (p, _) in s.punti() { z += p[2]; n += 1.0; }
        }
        2.0 * z / n // per ball: z1 + z2
    };
    let (zp, zm_) = (zm(&mut s, 0.0), zm(&mut s, PI));
    assert!(zp > 0.15 && zm_ < -0.15, "<z1+z2> at φ = 0: {zp}, at φ = π: {zm_}");
}

/// The default view: within 15 s of opening the first photon is absorbed and re-emitted.
#[test]
fn default_event_within_15_s() {
    let mut s = Scontro::new(&risultati(), crate::SEME);
    let mut assorbito = None;
    let mut emesso = None;
    let mut t = 0.0;
    while t < 15.0 {
        s.avanza(1.0 / 60.0);
        t += 1.0 / 60.0;
        if assorbito.is_none() && s.fase == Fase::Eccitato { assorbito = Some(t); }
        if emesso.is_none() && s.cont.emessi > 0 { emesso = Some(t); }
    }
    assert!(assorbito.is_some() && emesso.is_some(), "absorbed at {assorbito:?}, emitted at {emesso:?}");
    assert_eq!(s.cont.inviati, s.cont.assorbiti + s.cont.passati + s.cont.diffusi + (s.fin.is_some() as u64));
}


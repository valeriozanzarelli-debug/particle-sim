//! Tests of the laboratory. `cargo test --bin atomo` (fast ones); `cargo test --bin atomo -- --ignored`
//! for the long ones. The photon tests F1–F4 are the pre-registered checks of docs/ATOMO.md §6.

use crate::fotone3d::{punto, Fotone3D, Lancio, Misura, Polar};
use crate::fotscena::{retta, theta_da, Parametri};
use crate::skyrme::Profilo;

#[test]
fn nodo_b1() {
    let p = Profilo::calcola();
    assert!((p.numero_barionico() - 1.0).abs() < 1e-3);
    assert!((p.energia_12pi2() - 1.2315).abs() < 2e-3, "{}", p.energia_12pi2());
    assert!((p.raggio_rms_fm() - 0.59).abs() < 0.02);
}

/// The zig-zag atom runs backwards exactly: 60 ticks forward with a photon pulse, 60 back → the start.
#[test]
fn zigzag_indietro_esatto() {
    let mut at = crate::zz::Atomo::new();
    at.impulsi.push(crate::zz::Impulso { e0: 0.01, omega: 0.11, t0: 15.0, sigma: 10.0 });
    for _ in 0..60 { at.avanti(); }
    for _ in 0..60 { at.indietro(); }
    let err: f64 = at.psi.iter().zip(at.psi0.iter()).map(|(a, b)| (0..4).map(|q| (a[q] - b[q]).norm2()).sum::<f64>()).sum();
    eprintln!("|Δψ|² after forward+back = {err:e}");
    assert!(err < 1e-9);
}

/// The resonance comes out of the walk itself: the same pulse moves the dipole much more at the 1s→2p
/// frequency of the rule than 30 % off it (several minutes in a dev build: ignored by default).
#[test]
#[ignore]
fn zigzag_risonanza() {
    let w = crate::zz::livello(1, 2).omega - crate::zz::livello(-1, 1).omega;
    let amp = |omega: f64| {
        let mut at = crate::zz::Atomo::new();
        at.impulsi.push(crate::zz::Impulso { e0: 0.008, omega, t0: 125.0, sigma: 50.0 });
        let mut m: f64 = 0.0;
        for k in 0..600 { at.avanti(); if k > 450 { m = m.max(at.dipolo().abs()); } }
        m
    };
    let (a1, a2) = (amp(w), amp(1.3 * w));
    eprintln!("dipole after the pulse: on resonance {a1:.3}, 1.3× off {a2:.3}");
    assert!(a1 > 2.0 * a2);
}

// ─────────────────────────────── the free photon (pre-registered F1–F4) ───────────────────────────────

/// One flight of the photon: measurements at every tick, and the helicity at the focus.
struct Volo { n: usize, l: Lancio, storia: Vec<Misura>, elicita: f64, energia0: f64 }

fn vola(par: Parametri) -> Volo {
    let mut f = Fotone3D::new(par.n());
    let l = par.lancio();
    f.lancia(&l);
    let atteso = |t: u64| { let s = t as f64 - l.t_fuoco; [l.dir[0] * s, l.dir[1] * s, l.dir[2] * s] };
    let mut storia = vec![f.misura(atteso(0))];
    let t_f = l.t_fuoco.round() as u64;
    let mut elicita = f64::NAN;
    let fine = (2.0 * l.t_fuoco) as u64;
    for t in 1..=fine {
        let prima = if t == t_f { Some(f.campo_e()) } else { None };
        f.tick();
        if let Some(e) = prima { elicita = f.elicita(&e, l.dir, l.lambda); }
        storia.push(f.misura(atteso(t)));
    }
    let energia0 = storia[0].energia;
    Volo { n: par.n(), l, storia, elicita, energia0 }
}

impl Volo {
    /// Speed of the centre of energy along the flight: straight line over the whole flight.
    fn velocita(&self) -> f64 {
        let p: Vec<(f64, f64)> = self.storia.iter().enumerate().map(|(t, m)| (t as f64, punto(m.centro, self.l.dir))).collect();
        retta(&p).1
    }
}

/// F1a: a wide beam (w0 = 3λ) flies along a lattice axis at c within 0.5 %.
#[test]
fn fotone_f1_velocita_c() {
    let v = vola(Parametri { qual: 0, w0_lambda: 3.0, polar: Polar::Orizzontale, dir: 0 });
    let vel = v.velocita();
    let vc = Fotone3D::velocita_continuo(&v.l, v.n);
    eprintln!("F1a: v = {vel:.5} c (continuum, same pulse: {vc:.5} c); energy kept to {:e}", v.storia.last().unwrap().energia / v.energia0 - 1.0);
    assert!((vel - 1.0).abs() < 0.005, "v = {vel}");
}

/// F1b: the default beam (w0 = λ) against the continuum value for the same pulse, within 0.5 %
/// (quality 8 and 16 cells per λ).
#[test]
fn fotone_f1_velocita_continuo() {
    for qual in [0, 2] {
        let v = vola(Parametri { qual, ..Parametri::default() });
        let vel = v.velocita();
        let vc = Fotone3D::velocita_continuo(&v.l, v.n);
        eprintln!("F1b q{qual}: v = {vel:.5} c, continuum {vc:.5} c ({:+.3} %)", 100.0 * (vel / vc - 1.0));
        assert!((vel / vc - 1.0).abs() < 0.005);
    }
}

/// F2: helicity ±1 for circular light, 0 for linear, within 0.02 (8 and 16 cells per λ).
#[test]
fn fotone_f2_elicita() {
    for qual in [0, 2] {
        for (p, atteso) in [(Polar::CircolarePiu, 1.0), (Polar::CircolareMeno, -1.0), (Polar::Orizzontale, 0.0), (Polar::Verticale, 0.0)] {
            let v = vola(Parametri { qual, w0_lambda: 1.0, polar: p, dir: 0 });
            eprintln!("F2 q{qual}: {:?} helicity {:+.4} (expected {atteso:+})", p, v.elicita);
            assert!((v.elicita - atteso).abs() < 0.02, "{p:?}: {}", v.elicita);
        }
    }
}

/// F3: the beam spreads with half-angle λ/(π w0) within 5 % (w0 = 1.5 λ; 8 and 16 cells per λ).
#[test]
fn fotone_f3_divergenza() {
    for qual in [0, 2] {
        let par = Parametri { qual, w0_lambda: 1.5, polar: Polar::Orizzontale, dir: 0 };
        let v = vola(par);
        let th = theta_da(&v.storia, v.l.dir).unwrap();
        let att = 1.0 / (std::f64::consts::PI * par.w0_lambda);
        eprintln!("F3 q{qual}: θ = {th:.4} rad, λ/(π w0) = {att:.4} ({:+.2} %)", 100.0 * (th / att - 1.0));
        assert!((th / att - 1.0).abs() < 0.05);
    }
}

/// F4: the photon flies back exactly (every step of the rule is unitary).
#[test]
fn fotone_f4_indietro_esatto() {
    let mut f = Fotone3D::new(48);
    let l = Lancio { lambda: 8.0, w0: 8.0, sigma: 8.0, dir: [1.0, 0.0, 0.0], polar: Polar::CircolarePiu, t_fuoco: 10.0 };
    f.lancia(&l);
    let f0 = f.f.clone();
    for _ in 0..31 { f.tick(); }
    for _ in 0..31 { f.indietro(); }
    let err: f64 = f.f.iter().zip(&f0).map(|(a, b)| (0..3).map(|c| { let (x, y) = (a[c].re - b[c].re, a[c].im - b[c].im); (x * x + y * y) as f64 }).sum::<f64>()).sum();
    eprintln!("F4: Σ|ΔF|² after 31 ticks forward and back = {err:e} (Σ|F|² = 1)");
    assert!(err < 1e-8);
}

/// Informative (not pre-registered as pass/fail): the lattice in oblique directions.
#[test]
#[ignore]
fn fotone_direzioni_oblique() {
    for dir in 0..5 {
        let v = vola(Parametri { qual: 0, w0_lambda: 1.5, polar: Polar::CircolarePiu, dir });
        let vc = Fotone3D::velocita_continuo(&v.l, v.n);
        let th = theta_da(&v.storia, v.l.dir).unwrap_or(f64::NAN);
        eprintln!("dir {dir}: v = {:.4} c (continuum {vc:.4}), helicity {:+.3}, θ = {th:.4}", v.velocita(), v.elicita);
    }
}

#[test]
#[ignore]
fn bench_fotone() {
    use std::time::Instant;
    for q in 0..4 {
        let par = Parametri { qual: q, ..Parametri::default() };
        let mut f = Fotone3D::new(par.n());
        let l = par.lancio();
        let t = Instant::now(); f.lancia(&l); let tl = t.elapsed();
        let t = Instant::now(); for _ in 0..5 { f.tick(); } let tt = t.elapsed() / 5;
        let t = Instant::now(); for _ in 0..5 { f.misura([0.0; 3]); } let tm = t.elapsed() / 5;
        let e = f.campo_e();
        let t = Instant::now(); let mut nb = 0; for _ in 0..5 { nb = crate::fotscena::palline(&f, &l, 5, Some(&e), Default::default()).len(); } let tp = t.elapsed() / 5;
        eprintln!("qual {q} (n {}): lancia {tl:?}, tick {tt:?}, misura {tm:?}, palline {tp:?} ({nb} palline)", par.n());
    }
}

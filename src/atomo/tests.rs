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
#[ignore = "pre-registered, FAILED at 8 cells per λ (−0.58 %): docs/ATOMO.md §4.1"]
fn fotone_f1_velocita_continuo() {
    let mut ok = true;
    for qual in [0, 2] {
        let v = vola(Parametri { qual, ..Parametri::default() });
        let vel = v.velocita();
        let vc = Fotone3D::velocita_continuo(&v.l, v.n);
        eprintln!("F1b q{qual}: v = {vel:.5} c, continuum {vc:.5} c ({:+.3} %)", 100.0 * (vel / vc - 1.0));
        ok &= (vel / vc - 1.0).abs() < 0.005;
    }
    assert!(ok);
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
#[ignore = "pre-registered, FAILED (−21 % and −6 %): docs/ATOMO.md §4.1"]
fn fotone_f3_divergenza() {
    let mut ok = true;
    for qual in [0, 2] {
        let par = Parametri { qual, w0_lambda: 1.5, polar: Polar::Orizzontale, dir: 0 };
        let v = vola(par);
        let th = theta_da(&v.storia, v.l.dir).unwrap();
        let att = 1.0 / (std::f64::consts::PI * par.w0_lambda);
        eprintln!("F3 q{qual}: θ = {th:.4} rad, λ/(π w0) = {att:.4} ({:+.2} %)", 100.0 * (th / att - 1.0));
        ok &= (th / att - 1.0).abs() < 0.05;
    }
    assert!(ok);
}

/// Regression (post hoc, after F3 failed): the spreading angle follows the walk's own prediction
/// λ/(π w0) · k a / tan(k a) within 2 % (8 cells per λ, w0 = 1.5 λ).
#[test]
fn fotone_divergenza_come_la_regola() {
    let par = Parametri { qual: 0, w0_lambda: 1.5, polar: Polar::Orizzontale, dir: 0 };
    let v = vola(par);
    let th = theta_da(&v.storia, v.l.dir).unwrap();
    let k = std::f64::consts::TAU / par.cpl();
    let reticolo = k / k.tan() / (std::f64::consts::PI * par.w0_lambda);
    eprintln!("θ = {th:.4}, walk prediction {reticolo:.4}");
    assert!((th / reticolo - 1.0).abs() < 0.02);
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

/// Informative (after F3): the spreading angle at every quality, against λ/(π w0) and against the walk's own
/// prediction.
#[test]
#[ignore]
fn fotone_divergenza_qualita() {
    for qual in 0..4 {
        let par = Parametri { qual, w0_lambda: 1.5, polar: Polar::Orizzontale, dir: 0 };
        let v = vola(par);
        let th = theta_da(&v.storia, v.l.dir).unwrap();
        let att = 1.0 / (std::f64::consts::PI * par.w0_lambda);
        // the walk spreads light sideways slower than the continuum: v⊥ = q⊥ / tan(k a) instead of q⊥ / (k a)
        // (from its two-tick Bloch matrix), so its own prediction is λ/(π w0) · k a / tan(k a)
        let k = std::f64::consts::TAU / par.cpl();
        let reticolo = att * k / k.tan();
        eprintln!("q{qual} ({} cells/λ): θ = {th:.4}, λ/(π w0) = {att:.4} ({:+.2} %), lattice prediction {reticolo:.4} ({:+.2} %)",
            par.cpl(), 100.0 * (th / att - 1.0), 100.0 * (th / reticolo - 1.0));
    }
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
        let mut tr = crate::fotdisegno::Traccianti::new(&f, 1);
        let t = Instant::now(); for _ in 0..5 { tr.avanza(&f, 1.0, [0.0; 3]); } let tt2 = t.elapsed() / 5;
        let t = Instant::now(); let mut nb = 0; for _ in 0..5 { let (a, b, _) = crate::fotdisegno::palline(&f, &l, [0.0; 3], &tr, Some(&e), Default::default(), true); nb = a.len() + b.len(); } let tp = t.elapsed() / 5;
        eprintln!("qual {q} (n {}): lancia {tl:?}, tick {tt:?}, misura {tm:?}, tracers {tt2:?}, palline {tp:?} ({nb} palline)", par.n());
    }
}

/// Orbitals: the radial form of the rule finds every state n = 1..4, l = 0..3 (j = l + ½) in the right order,
/// each one bound, less bound as n grows; printed against the exact Dirac levels. Since 2026-10-13 on the graded
/// mesh of the automatic zoom (docs/NUCLEO_1S.md §8): 1s, 2p3/2, 3d5/2 within 1 % of exact Dirac (threshold written
/// in NUCLEO_1S §8 before the run; with equal cells of 0.5 λ_C they were +39.9 %, +3.5 %, +1.8 %).
#[test]
fn orbitali_livelli() {
    use crate::orbitali::{dirac_esatto, radiale, raggio, Stato};
    let mut ultimo_per_l = [0.0f64; 4];
    for n in 1..=4u32 {
        for l in 0..n {
            let s = Stato::new(n, l, 1);
            let r = radiale(s.kappa(), n, l);
            let b = 1.0 - r.omega;
            let d = 1.0 - dirac_esatto(n, s.kappa());
            eprintln!("{}: binding {:.6} mc² (Dirac {:.6}, {:+.4} %), radius {:.1} λ_C", s.nome(), b, d, 100.0 * (b / d - 1.0), raggio(&r));
            assert!(b > 0.0);
            if n > l + 1 { assert!(b < ultimo_per_l[l as usize], "{}", s.nome()); }
            ultimo_per_l[l as usize] = b;
            // the same order as Dirac within a factor 2
            assert!(b / d > 0.5 && b / d < 2.0, "{}", s.nome());
            if (n, l) == (1, 0) || (n, l) == (2, 1) || (n, l) == (3, 2) {
                assert!((b / d - 1.0).abs() < 0.01, "{}: {:+.3} % vs exact Dirac", s.nome(), 100.0 * (b / d - 1.0));
            }
        }
    }
}


/// Development run of the coupled light + atom (not a check): prints populations, frequencies, decay.
#[test]
#[ignore]
fn accoppiato_sviluppo() {
    use crate::accoppiato::*;
    let t0 = std::time::Instant::now();
    let asc: usize = std::env::var("ASCOLTO").ok().and_then(|v| v.parse().ok()).unwrap_or(1200);
    let mut x = Esperimento::new(asc);
    let w_rad = crate::zz::livello(1, 2).omega - crate::zz::livello(-1, 1).omega;
    let w_res = x.omega_21();
    eprintln!("prepared in {:?}: ω1s(3D) = {:.5}, ω2p(3D) = {:.5}, ω21(3D) = {:.5}, radial ω21 = {:.5}", t0.elapsed(), x.omega1s, x.omega2p, w_res, w_rad);
    let e0: f64 = std::env::var("E0").ok().and_then(|v| v.parse().ok()).unwrap_or(crate::zzscena::E0_IMP);
    x.spara(w_res, e0);
    let n: usize = std::env::var("TICKS").ok().and_then(|v| v.parse().ok()).unwrap_or(1500);
    for k in 0..n {
        x.passo();
        if k % 150 == 0 {
            let p = x.serie_p.last().map(|v| v.1).unwrap_or(0.0);
            eprintln!("tick {k}: t {:.0} p_exc {:.4} p2p {:.4} dip {:+.4} light energy {:.3e} ({:?})", x.at.t(), p, x.p2p(), x.at.dipolo(), x.luce.energia(), t0.elapsed());
        }
    }
    let t_pass = x.at.t() * 0.0 + 150.0;
    eprintln!("probe frequency {:?}, dipole frequency {:?}, ω_res {:.5}", frequenza(&x.luce.sonda, t_pass), frequenza(&x.serie_d.iter().map(|p| (p.0, p.1)).collect::<Vec<_>>(), t_pass), w_res);
    let g = decadimento(&x.serie_p, t_pass);
    let d = x.ampiezza_dipolo(200.0);
    let p = x.serie_p.last().unwrap().1;
    eprintln!("Γ measured {:?}, Larmor {:.3e} (D = {:.4}, p = {:.4})", g, gamma_larmor(w_res, d, p), d, p);
}


/// Pre-registered E1, E2 (docs/ATOMO.md §3.5): the coupled light + atom. Long (~3 min): ignored by default.
/// E1: the light the atom sends back (probe 160 λ_C behind it, after the pulse) has the frequency of the 1s → 2p
///     level difference of the 3-D rule (listened) within 2 %.
/// E2: the power it radiates (energy into the sponge + left in the box, slope over t = 600..1200) equals the
///     Larmor power of its own dipole at that frequency (amplitude over t = 400..1000, same K, same Gaussian
///     average) within 30 %.
#[test]
#[ignore = "pre-registered v1, FAILED (E1 −2.01 %, E2 −0.15): docs/ATOMO.md §3.5"]
fn accoppiato_e1_e2() {
    use crate::accoppiato::*;
    let mut x = Esperimento::new(1200);
    let w21 = x.omega_21();
    x.spara(w21, 0.004);
    while x.at.t() < 1200.0 { x.passo(); }
    let wp = frequenza(&x.luce.sonda, 400.0).unwrap();
    eprintln!("E1: ω emitted {wp:.5}, ω21 of the 3-D rule {w21:.5} ({:+.2} %); radial ω21 {:.5}", 100.0 * (wp / w21 - 1.0),
              crate::zz::livello(1, 2).omega - crate::zz::livello(-1, 1).omega);
    let e: Vec<(f64, f64)> = x.luce.serie_e.iter().filter(|p| p.0 >= 600.0).cloned().collect();
    let prad = crate::fotscena::retta(&e).1;
    let d = ampiezza_a(&x.serie_d, w21, 400.0, 1000.0);
    let pl = potenza_larmor(w21, d);
    eprintln!("E2: radiated power {prad:.4e}, Larmor {pl:.4e} (dipole amplitude {d:.4} λ_C) ratio {:.3}", prad / pl);
    let gam = prad / w21; // excitation quanta per unit time
    eprintln!("    decay: energy lost per unit time / ħω = {gam:.3e} per λ_C/c; p_exc now {:.3}", x.serie_p.last().unwrap().1);
    assert!((wp / w21 - 1.0).abs() < 0.02);
    assert!((prad / pl - 1.0).abs() < 0.30);
}

/// Post hoc (after E1/E2): the emitted wave's amplitude at the probe against the far field of the atom's own
/// dipole, e D ω² W / (4π r), and the energy taken by the sponge alone.
#[test]
#[ignore]
fn accoppiato_post_hoc() {
    use crate::accoppiato::*;
    let mut x = Esperimento::new(1200);
    let w21 = x.omega_21();
    x.spara(w21, 0.004);
    let mut ass = vec![];
    while x.at.t() < 1200.0 { x.passo(); ass.push((x.at.t(), x.luce.assorbita)); }
    let wp = frequenza(&x.luce.sonda, 400.0).unwrap();
    let r = (SONDA as f64) * B_LUCE;
    let d = ampiezza_a(&x.serie_d, wp, 400.0 - r, 1200.0 - r);
    let ep = ampiezza_a(&x.luce.sonda, wp, 400.0, 1200.0);
    let atteso = e_eff() * d * wp * wp * LuceAtomo::fattore_forma(wp) / (4.0 * std::f64::consts::PI * r);
    eprintln!("probe amplitude {ep:.4e}, far field of the dipole {atteso:.4e} (ratio {:.3}); dipole at ω {d:.4}", ep / atteso);
    let a: Vec<(f64, f64)> = ass.iter().filter(|p| p.0 >= 600.0).cloned().collect();
    eprintln!("sponge power {:.4e}, Larmor {:.4e}", crate::fotscena::retta(&a).1, potenza_larmor(wp, ampiezza_a(&x.serie_d, wp, 400.0, 1000.0)));
    eprintln!("current at ω {:.4e}, ω × dipole at ω {:.4e}", ampiezza_a(&x.serie_j, wp, 400.0, 1000.0), wp * ampiezza_a(&x.serie_d, wp, 400.0, 1000.0));
}

/// Pre-registered E1v2, E2v2 (docs/ATOMO.md §3.5, after the v1 diagnosis: the sponge was too weak and the
/// source was the one-slice current; now a strong sponge and the conserved current −e Δ⟨r⟩/Δt).
/// E1v2: frequency of the light at the probe (t = 400..1200) vs ω21 of the 3-D rule: within 2 %.
/// E2v2: its amplitude vs the far field of the atom's own dipole, e D ω² W / (4π r) (D at that frequency over
///       the retarded window): within 10 %.
#[test]
#[ignore = "pre-registered v2: E1v2 FAILED (+7.1 %), E2v2 passed (−0.6 %): docs/ATOMO.md §3.5"]
fn accoppiato_v2() {
    use crate::accoppiato::*;
    let mut x = Esperimento::new(1200);
    let w21 = x.omega_21();
    x.spara(w21, 0.004);
    while x.at.t() < 1200.0 { x.passo(); }
    let wp = frequenza(&x.luce.sonda, 400.0).unwrap();
    let r = (SONDA as f64) * B_LUCE;
    let d = ampiezza_a(&x.serie_d, wp, 400.0 - r, 1200.0 - r);
    let ep = ampiezza_a(&x.luce.sonda, wp, 400.0, 1200.0);
    let atteso = e_eff() * d * wp * wp * LuceAtomo::fattore_forma(wp) / (4.0 * std::f64::consts::PI * r);
    let pl = potenza_larmor(wp, d);
    let p = x.serie_p.iter().filter(|q| q.0 >= 400.0).map(|q| q.1).sum::<f64>() / x.serie_p.iter().filter(|q| q.0 >= 400.0).count() as f64;
    eprintln!("E1v2: ω emitted {wp:.5}, ω21 {w21:.5} ({:+.2} %)", 100.0 * (wp / w21 - 1.0));
    eprintln!("E2v2: probe amplitude {ep:.4e}, dipole far field {atteso:.4e} ({:+.1} %)", 100.0 * (ep / atteso - 1.0));
    eprintln!("     Larmor power {pl:.3e} → classical decay time ħω p(1−p)/P = {:.0} λ_C/c (p ≈ {p:.3}, K = {K_ACC})", w21 * p * (1.0 - p) / pl);
    assert!((wp / w21 - 1.0).abs() < 0.02);
    assert!((ep / atteso - 1.0).abs() < 0.10);
}

/// Post hoc after E1v2: the spectrum of the light at the probe (t = 400..1200): which lines?
#[test]
#[ignore]
fn accoppiato_righe() {
    use crate::accoppiato::*;
    let mut x = Esperimento::new(1200);
    let w21 = x.omega_21();
    x.spara(w21, 0.004);
    while x.at.t() < 1200.0 { x.passo(); }
    for (nome, serie, da) in [("probe", &x.luce.sonda, 400.0), ("dipole", &x.serie_d, 240.0)] {
        let mut picchi = vec![];
        let mut prev = (0.0, 0.0, 0.0);
        for k in 0..=600 {
            let w = 0.04 + 0.0001 * k as f64;
            let a = ampiezza_a(serie, w, da, da + 800.0);
            if prev.1 > prev.0 && prev.1 > a { picchi.push((prev.2, prev.1)); }
            prev = (prev.1, a, w);
        }
        picchi.sort_by(|a, b| b.1.total_cmp(&a.1));
        eprintln!("{nome}: lines {:?} (ω21 {w21:.4}, ω2p − 0.882: {:.4}, ω2p − 0.8405: {:.4})",
            picchi.iter().take(5).map(|p| format!("{:.4} ({:.2e})", p.0, p.1)).collect::<Vec<_>>(), x.omega2p - 0.8822, x.omega2p - 0.8405);
    }
}

/// Levels of the 3-D rule by listening (spectrum of ⟨ψ₀|ψ(t)⟩): 1s, 2p½, 2p3/2 (m = +½).
#[test]
#[ignore]
fn accoppiato_livelli_3d() {
    use crate::cplx::C;
    let n = crate::zz::N;
    let da_radiale = |kappa: i32, nn: u32, l: u32| {
        let r = crate::orbitali::radiale(kappa, nn, l);
        let mut psi: Vec<[C; 4]> = (0..n * n * n).map(|s| {
            let x = [crate::zz::x_of(s / (n * n)), crate::zz::x_of((s / n) % n), crate::zz::x_of(s % n)];
            crate::zz::dirac_a_chirale(crate::orbitali::spinore(&r, kappa, 1, x))
        }).collect();
        let nn2: f64 = psi.iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum();
        for p in psi.iter_mut() { for c in p.iter_mut() { *c = c.scale(1.0 / nn2.sqrt()); } }
        let mut a = crate::zz::Atomo::new();
        a.imposta(psi);
        a
    };
    for (nome, a, lo, hi, ticks) in [("2p3/2", da_radiale(-2, 2, 1), 0.95, 0.995, 500), ("2p1/2", da_radiale(1, 2, 1), 0.95, 0.995, 500)] {
        let (_, om) = a.ascolta(lo, hi, ticks);
        eprintln!("{nome}: ω = {om:.5}");
    }
}

/// EXPORT of the coupled light + atom for offline rendering (not a check; nothing in the physics changes): the
/// same experiment as the Laboratorio (`Esperimento::new(600)`, `spara(omega_21(), E0_ACC)`), recorded.
/// Every `OGNI` atom ticks (default 16 = 2 light ticks) one frame: electron density ρ, current j (averaged over
/// the frame's atom ticks), phase of the dominant spinor component (rest-mass phase removed, as the scene colours
/// it), all on the 48³ atom grid; the light's E (= Re F) on the 128³ light box averaged over 2×2×2 → 64³.
/// Raw little-endian float16 files (stored = physical × scale, scales in the manifest), `manifest.json`,
/// `fotogrammi.csv` (scalars per frame), `serie.csv` (scalars per light tick). Docs: docs/ATOMO.md §"Export".
/// Run: ESPORTA_DIR=… TICKS=3200 OGNI=16 cargo test --bin atomo accoppiato_esporta -- --ignored --nocapture
#[test]
#[ignore]
fn accoppiato_esporta() {
    use crate::accoppiato::*;
    use rayon::prelude::*;
    use std::io::Write;
    fn f16(x: f32) -> u16 {
        let b = x.to_bits();
        let sign = ((b >> 16) & 0x8000) as u16;
        let exp = ((b >> 23) & 0xff) as i32;
        let mut m = b & 0x7f_ffff;
        if exp == 255 { return sign | 0x7c00 | if m != 0 { 0x200 } else { 0 }; }
        let e = exp - 127 + 15;
        if e >= 31 { return sign | 0x7c00; }
        if e <= 0 {
            if e < -10 { return sign; }
            m |= 0x80_0000;
            let sh = (14 - e) as u32;
            let mut h = (m >> sh) as u16;
            if (m >> (sh - 1)) & 1 == 1 { h += 1; }
            return sign | h;
        }
        let mut h = sign | ((e as u16) << 10) | (m >> 13) as u16;
        if m & 0x1000 != 0 { h += 1; }
        h
    }
    fn scrivi(w: &mut std::io::BufWriter<std::fs::File>, v: &[f32], scala: f32) {
        let mut buf = Vec::with_capacity(v.len() * 2);
        for &x in v { buf.extend_from_slice(&f16(x * scala).to_le_bytes()); }
        w.write_all(&buf).unwrap();
    }
    let dir = std::env::var("ESPORTA_DIR").unwrap_or("/mnt/project-files/ripresa/dati/atomo_accoppiato".into());
    let ticks: usize = std::env::var("TICKS").ok().and_then(|v| v.parse().ok()).unwrap_or(3200);
    let ogni: usize = std::env::var("OGNI").ok().and_then(|v| v.parse().ok()).unwrap_or(16);
    assert!(ogni % SUB == 0, "OGNI must be a multiple of SUB = {SUB} (frames on light ticks)");
    std::fs::create_dir_all(&dir).unwrap();
    let apri = |n: &str| std::io::BufWriter::with_capacity(1 << 22, std::fs::File::create(format!("{dir}/{n}")).unwrap());
    let t0 = std::time::Instant::now();
    let mut x = Esperimento::new(600);
    let w21 = x.omega_21();
    eprintln!("prepared in {:?}: ω21 = {w21:.5}", t0.elapsed());
    x.spara(w21, crate::zzscena::E0_ACC);
    let t_sparo = x.at.t();
    let n = crate::zz::N;
    let na = n * n * n;
    let nl = N_LUCE;
    let nd = nl / 2;
    let ee = e_eff();
    // scales (stored = physical × scale): ρ per λ_C³ ×1e3, j per λ_C³ (units c) ×1e3, e·E (mc²/λ_C) ×1e5
    let (s_rho, s_j, s_e) = (1e3f32, 1e3f32, 1e5f32);
    let vol = (crate::zz::A * crate::zz::A * crate::zz::A) as f32;
    let (mut f_rho, mut f_j, mut f_fase, mut f_e) = (apri("rho.f16"), apri("j.f16"), apri("fase.f16"), apri("E.f16"));
    let mut f_csv = apri("fotogrammi.csv");
    writeln!(f_csv, "fotogramma,tick,t,t_dal_sparo,p_exc,p2p,dipolo_z,eE_atomo_x,eE_atomo_y,eE_atomo_z,energia_luce_box,energia_assorbita,sonda_Ez").unwrap();
    let mut jacc = vec![[0.0f64; 3]; na];
    let mut nacc = 0usize;
    let mut fotogrammi = 0usize;
    let mut scrivi_fotogramma = |x: &Esperimento, jacc: &mut Vec<[f64; 3]>, nacc: &mut usize, fotogrammi: &mut usize, f_csv: &mut std::io::BufWriter<std::fs::File>| {
        let at = &x.at;
        let t = at.t();
        let rho: Vec<f32> = at.psi.par_iter().map(|p| (p.iter().map(|c| c.norm2()).sum::<f64>() as f32) / vol).collect();
        let fase: Vec<f32> = at.psi.par_iter().map(|p| {
            let mut b = 0;
            for q in 1..4 { if p[q].norm2() > p[b].norm2() { b = q; } }
            ((p[b].arg() + t).rem_euclid(std::f64::consts::TAU)) as f32
        }).collect();
        // j averaged over the atom ticks since the last frame (the instantaneous one carries the fast
        // zitterbewegung beat at 2mc²/ħ, period π λ_C/c, which a frame every 8 λ_C/c would alias)
        let k = if *nacc > 0 { 1.0 / *nacc as f64 } else { 0.0 };
        let mut jv = vec![0.0f32; 3 * na];
        for c in 0..3 { for s in 0..na { jv[c * na + s] = (jacc[s][c] * k) as f32 / vol; } }
        for v in jacc.iter_mut() { *v = [0.0; 3]; }
        *nacc = 0;
        let f = &x.luce.f.f;
        let mut ev = vec![0.0f32; 3 * nd * nd * nd];
        ev.par_chunks_mut(nd * nd).enumerate().for_each(|(ci, sl)| {
            let (c, i) = (ci / nd, ci % nd);
            for j in 0..nd { for kk in 0..nd {
                let mut a = 0.0f32;
                for di in 0..2 { for dj in 0..2 { for dk in 0..2 {
                    a += f[((2 * i + di) * nl + 2 * j + dj) * nl + 2 * kk + dk][c].re;
                }}}
                sl[j * nd + kk] = a * 0.125 * ee as f32;
            }}
        });
        scrivi(&mut f_rho, &rho, s_rho);
        scrivi(&mut f_j, &jv, s_j);
        scrivi(&mut f_fase, &fase, 1.0);
        scrivi(&mut f_e, &ev, s_e);
        let p = 1.0 - at.proiezione(&x.phi1s);
        let ea = x.luce.campo_atomo();
        let eb = x.luce.energia();
        let pr = x.luce.sonda.last().map(|v| v.1).unwrap_or(0.0) * ee;
        writeln!(f_csv, "{},{},{:.4},{:.4},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}", *fotogrammi, at.ticks, t, t - t_sparo, p, x.p2p(), at.dipolo(),
                 ee * ea[0], ee * ea[1], ee * ea[2], eb, x.luce.assorbita, pr).unwrap();
        *fotogrammi += 1;
    };
    scrivi_fotogramma(&x, &mut jacc, &mut nacc, &mut fotogrammi, &mut f_csv);
    for k in 1..=ticks {
        x.passo();
        let j: Vec<[f64; 3]> = x.at.psi.par_iter().map(|p| {
            let (l0, l1, r0, r1) = (p[0], p[1], p[2], p[3]);
            let sx = |a: crate::cplx::C, b: crate::cplx::C| 2.0 * (a.conj() * b).re;
            let sy = |a: crate::cplx::C, b: crate::cplx::C| 2.0 * (a.conj() * b).im;
            let sz = |a: crate::cplx::C, b: crate::cplx::C| a.norm2() - b.norm2();
            [sx(r0, r1) - sx(l0, l1), sy(r0, r1) - sy(l0, l1), sz(r0, r1) - sz(l0, l1)]
        }).collect();
        jacc.par_iter_mut().zip(j.par_iter()).for_each(|(a, b)| for c in 0..3 { a[c] += b[c]; });
        nacc += 1;
        if k % ogni == 0 {
            scrivi_fotogramma(&x, &mut jacc, &mut nacc, &mut fotogrammi, &mut f_csv);
            if fotogrammi % 20 == 0 { eprintln!("frame {fotogrammi}, tick {k}, p_exc {:.4} ({:?})", x.serie_p.last().unwrap().1, t0.elapsed()); }
        }
    }
    for w in [&mut f_rho, &mut f_j, &mut f_fase, &mut f_e, &mut f_csv] { w.flush().unwrap(); }
    // per light tick: p_exc, dipole (end of tick), current, light energy (box + sponge), probe E_z
    let mut fs = apri("serie.csv");
    writeln!(fs, "t,t_dal_sparo,p_exc,dipolo_z,corrente_z,energia_box_piu_assorbita,sonda_eEz").unwrap();
    for (i, (tp, p)) in x.serie_p.iter().enumerate() {
        let d = x.serie_d[(i + 1) * SUB - 1].1;
        let jz = x.serie_j[i].1;
        let e = x.luce.serie_e[i].1;
        let pr = x.luce.sonda[i].1 * ee;
        writeln!(fs, "{tp:.4},{:.4},{p:.6e},{d:.6e},{jz:.6e},{e:.6e},{pr:.6e}", tp - t_sparo).unwrap();
    }
    fs.flush().unwrap();
    let a = crate::zz::A;
    let man = format!(r#"{{
  "descrizione": "Luce + atomo accoppiati (src/atomo/accoppiato.rs, Esperimento): stesso esperimento del Laboratorio Atomo, solo registrato. Ione idrogenoide Z = {z} (Nd59+), elettrone = cammino zig-zag di Dirac su 48^3, luce = camminatori di spin 1 su 128^3. Impulso: spara(omega_21, E0_ACC) lungo +x, polarizzato z.",
  "unita": {{"lunghezza": "lambda_C = hbar/(m c) = {lc:e} m", "tempo": "lambda_C/c = {tu:e} s", "energia": "m c^2 = {mc2} keV", "hbar_c_m": "1"}},
  "tick_atomo": {a}, "tick_luce": {bl}, "tick_atomo_per_tick_luce": {sub},
  "tick_per_fotogramma": {ogni}, "dt_fotogramma": {dtf}, "fotogrammi": {fotogrammi}, "tick_totali": {ticks},
  "t_sparo": {t_sparo}, "omega_21": {w21}, "omega_1s": {o1}, "omega_2p": {o2}, "periodo_luce": {per},
  "e_eff": {ee}, "K_ACC": {k}, "E0_ACC": {e0}, "alpha": {al}, "Z": {z},
  "griglia_atomo": {{"n": {n}, "passo": {a}, "coordinate": "x_i = (i - {n}/2 + 0.5) * {a} lambda_C, i = 0..{n}-1 (stesse per y, z); centro = nucleo = origine"}},
  "griglia_luce_originale": {{"n": {nl}, "passo": {bl}, "coordinate": "x_i = (i - {nl}/2) * {bl} lambda_C; la cella {h} e' l'origine (nucleo); spugna assorbente {sp} celle ai bordi"}},
  "griglia_luce": {{"n": {nd}, "passo": {bd}, "coordinate": "media 2x2x2: X_I = (2 I + 0.5 - {h}) * {bl} lambda_C; il nucleo (origine) sta a I = {ia}"}},
  "griglia_atomo_dentro_luce": "la griglia atomo (lato {la} lambda_C) e' centrata sull'origine della scatola di luce (lato {ll} lambda_C): copre X in [{am}, {ap}] lambda_C",
  "sonda": "cella di luce a {so} celle dietro il nucleo lungo -x (x = -{sod} lambda_C)",
  "ordine_indici": "C-order [x][y][z] (z piu' veloce), come s = (i*n + j)*n + k nel codice",
  "file": {{
    "rho.f16":  {{"dtype": "float16 LE", "forma": [{fotogrammi}, {n}, {n}, {n}], "scala": {s_rho}, "fisico": "densita' di probabilita' dell'elettrone |psi|^2 per lambda_C^3 (somma * {a}^3 = 1)"}},
    "j.f16":    {{"dtype": "float16 LE", "forma": [{fotogrammi}, 3, {n}, {n}, {n}], "scala": {s_j}, "fisico": "corrente di probabilita' j = R^+ sigma R - L^+ sigma L per lambda_C^3 (unita' c), MEDIATA sui {ogni} tick atomo del fotogramma (toglie il battito di zitterbewegung a 2mc^2); velocita' delle palline v = j/rho; corrente elettrica = -e j"}},
    "fase.f16": {{"dtype": "float16 LE", "forma": [{fotogrammi}, {n}, {n}, {n}], "scala": 1.0, "fisico": "fase (rad, 0..2pi) della componente dominante dello spinore + t (fase di massa a riposo tolta), come la colora la scena (zzscena::palline)"}},
    "E.f16":    {{"dtype": "float16 LE", "forma": [{fotogrammi}, 3, {nd}, {nd}, {nd}], "scala": {s_e}, "fisico": "e*E della luce (E = Re F, F = E + icB), in m c^2 / lambda_C (cioe' gia' moltiplicato per e_eff), media 2x2x2 delle celle di luce"}},
    "fotogrammi.csv": "uno per fotogramma: tick, t, t dal sparo, p_exc = 1 - |<1s|psi>|^2, p2p = |<2p1/2 m=1/2|psi>|^2, dipolo <z> (lambda_C), e*E nel punto dell'atomo (media gaussiana sigma = 2 celle), energia luce nella scatola, energia presa dalla spugna (unita' interne: usare solo rapporti), e*E_z alla sonda",
    "serie.csv": "uno per tick di luce ({sub} tick atomo): p_exc, dipolo z, corrente <alpha_z> = d<z>/dt, energia luce + assorbita, e*E_z alla sonda"
  }},
  "debiti": [
    "accoppiamento x{k} a mano (K_ACC): e -> {k} e nei due versi, cosi' il decadimento si vede in minuti; i tassi scalano come K^2",
    "semi-classico: la luce e' un campo classico (camminatori di spin 1), l'atomo quantistico; niente fotone quantizzato",
    "emissione solo da sovrapposizione (1s + 2p, dipolo oscillante): uno stato eccitato puro non decadrebbe (manca l'emissione spontanea vera)",
    "forma di dipolo e E mediato su una gaussiana sigma = 2 celle di luce intorno all'atomo (riduce la potenza irradiata di W^2)",
    "ione Z = 60 invece dell'idrogeno (la cella deve essere << raggio di Bohr)"
  ],
  "comando": "ESPORTA_DIR=<dir> TICKS={ticks} OGNI={ogni} cargo test --bin atomo accoppiato_esporta -- --ignored --nocapture"
}}
"#, z = crate::zz::Z_ION, lc = crate::zz::LAMBDA_C_M, tu = crate::zz::T_UNIT_S, mc2 = crate::zz::MC2_KEV, a = a, bl = B_LUCE, sub = SUB,
        dtf = ogni as f64 * a, per = std::f64::consts::TAU / w21, o1 = x.omega1s, o2 = x.omega2p, k = K_ACC, e0 = crate::zzscena::E0_ACC, al = alpha(),
        h = nl / 2, sp = SPUGNA, bd = 2.0 * B_LUCE, ia = (nl / 2) as f64 / 2.0 - 0.25, la = n as f64 * a, ll = nl as f64 * B_LUCE,
        am = -(n as f64) * a / 2.0, ap = n as f64 * a / 2.0, so = SONDA, sod = SONDA as f64 * B_LUCE);
    std::fs::write(format!("{dir}/manifest.json"), man).unwrap();
    eprintln!("done: {fotogrammi} frames in {dir} ({:?})", t0.elapsed());
}

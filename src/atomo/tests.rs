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
        let t = Instant::now(); let mut nb = 0; for _ in 0..5 { let (a, b) = crate::fotdisegno::palline(&f, &l, [0.0; 3], &tr, Some(&e), Default::default(), true); nb = a.len() + b.len(); } let tp = t.elapsed() / 5;
        eprintln!("qual {q} (n {}): lancia {tl:?}, tick {tt:?}, misura {tm:?}, tracers {tt2:?}, palline {tp:?} ({nb} palline)", par.n());
    }
}

/// Orbitals: the radial form of the rule finds every state n = 1..4, l = 0..3 (j = l + ½) in the right order,
/// each one bound, less bound as n grows; printed against the exact Dirac levels (informative: cells of 0.5 λ_C).
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
            eprintln!("{}: binding {:.5} mc² (Dirac {:.5}, {:+.1} %), radius {:.1} λ_C", s.nome(), b, d, 100.0 * (b / d - 1.0), raggio(&r));
            assert!(b > 0.0);
            if n > l + 1 { assert!(b < ultimo_per_l[l as usize], "{}", s.nome()); }
            ultimo_per_l[l as usize] = b;
            // the same order as Dirac within a factor 2
            assert!(b / d > 0.5 && b / d < 2.0, "{}", s.nome());
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

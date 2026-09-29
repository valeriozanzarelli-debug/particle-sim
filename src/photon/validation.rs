//! Validation of the photon against its real properties.
//!
//!   cargo test --release --bin fotone -- --nocapture --test-threads=1

use bevy::math::Vec3;

use crate::entangle::{correlation_local, correlation_quantum};
use crate::physics::*;

const DEG: f32 = std::f32::consts::PI / 180.0;

fn sim(e: Experiment, f: impl FnOnce(&mut Lab)) -> Sim {
    let mut lab = Lab::preset(e);
    f(&mut lab);
    Sim::new(lab)
}

/// Every grain moves at exactly c; its phase turns 2π per wavelength —
/// for any wavelength (no dispersion in vacuum).
#[test]
fn speed_of_light_and_phase_rotation() {
    for lambda in [0.5f32, 1.0, 2.0] {
        let mut s = sim(Experiment::FreeBeam, |l| { l.lambda = lambda; l.grains = 2000; });
        s.launch();
        let k = s.lab.k();
        let before: Vec<(Vec3, f32)> = s.photons[0].grains.iter().map(|g| (g.pos, g.phase(k))).collect();
        s.step(10.0);
        let mut worst_speed = 0.0f32;
        let mut worst_turns = 0.0f32;
        for (g, (p0, ph0)) in s.photons[0].grains.iter().zip(&before) {
            worst_speed = worst_speed.max(((g.pos - *p0).length() - 10.0).abs());
            let turns = (g.phase(k) - ph0) / TAU;
            worst_turns = worst_turns.max((turns - 10.0 / lambda).abs());
        }
        println!("λ={lambda}: |Δd − c·t| ≤ {worst_speed:.1e}, |turns − d/λ| ≤ {worst_turns:.1e}");
        assert!(worst_speed < 1e-3);
        assert!(worst_turns < 1e-3);
    }
}

/// Circular light: E turns about the direction of flight, one full turn per
/// period; the sense of rotation is the helicity (+1 or −1).
#[test]
fn circular_polarization_helicity() {
    for (pol, sign) in [(Polarization::CircPlus, 1.0f32), (Polarization::CircMinus, -1.0)] {
        let mut s = sim(Experiment::FreeBeam, |l| { l.polarization = pol; l.grains = 100; });
        s.launch();
        let k = s.lab.k();
        let lam = s.lab.lambda;
        for g in &s.photons[0].grains {
            // At a fixed point in space, time runs: clock advances, phase stays
            let e0 = g.e_field(k, 0.0);
            let e8 = g.e_field(k, lam / 8.0);
            let e4 = g.e_field(k, lam / 4.0);
            let e1 = g.e_field(k, lam);
            assert!(sign * e0.cross(e8).dot(g.dir) > 0.0, "wrong sense of rotation");
            assert!(e0.dot(e4).abs() < 1e-3, "quarter period must be 90°");
            assert!((e0 - e1).length() < 1e-3, "one period must be a full turn");
            assert!((e0.length() - e4.length()).abs() < 1e-3, "circular: constant |E|");
        }
        println!("{pol:?}: E turns with helicity {sign:+}, 360° per period, |E| constant");
    }
}

/// E is always perpendicular to the direction of flight — also after the
/// slits bend the grains.
#[test]
fn field_is_transverse() {
    let mut s = sim(Experiment::DoubleSlit, |l| { l.grains = 20_000; l.seed = 3; });
    let mut checked = 0;
    for _ in 0..50 {
        s.launch();
        s.step(60.0); // past the wall
        if let Some(ph) = s.photons.first() {
            for g in ph.grains.iter().filter(|g| g.alive && g.diffracted) {
                let e = g.e_field(s.lab.k(), s.clock);
                assert!(e.dot(g.dir).abs() < 1e-4, "longitudinal E");
                checked += 1;
            }
        }
        while s.in_flight() { s.step(40.0); }
        if checked > 1000 { break; }
    }
    println!("E·u = 0 on {checked} diffracted grains");
    assert!(checked > 1000);
}

/// A Gaussian beam spreads exactly as w(z) = w0 √(1 + (z/z_R)²), and at the
/// waist Δx·Δθ = λ/4π: the minimum allowed by Heisenberg (Δx·Δp = ħ/2).
#[test]
fn gaussian_beam_diffraction_and_uncertainty() {
    let mut s = sim(Experiment::FreeBeam, |l| { l.grains = 60_000; l.packet_len = 0.5; });
    s.launch();
    let zr = s.lab.rayleigh();

    let sigma = |v: &mut dyn Iterator<Item = f32>| {
        let xs: Vec<f32> = v.collect();
        let m = xs.iter().sum::<f32>() / xs.len() as f32;
        (xs.iter().map(|x| (x - m).powi(2)).sum::<f32>() / xs.len() as f32).sqrt()
    };
    let g0 = &s.photons[0].grains;
    let sx = sigma(&mut g0.iter().map(|g| g.pos.x));
    let st = sigma(&mut g0.iter().map(|g| g.dir.x));
    let product = sx * st;
    let bound = s.lab.lambda / (4.0 * std::f32::consts::PI);
    println!("waist: Δx·Δθ = {product:.5}  vs  λ/4π = {bound:.5}");
    assert!((product / bound - 1.0).abs() < 0.03);

    let mut z_done = 0.0;
    for z in [zr, 2.0 * zr, 5.0 * zr] {
        s.step(z - z_done);
        z_done = z;
        let g = &s.photons[0].grains;
        let zm = g.iter().map(|g| g.pos.z).sum::<f32>() / g.len() as f32;
        let w = 2.0 * sigma(&mut g.iter().map(|g| g.pos.x));
        let theory = s.lab.beam_radius(zm);
        println!("z = {zm:6.1}: w = {w:6.2}   theory {theory:6.2}");
        assert!((w / theory - 1.0).abs() < 0.03);
    }
}

fn run_until_detected(s: &mut Sim, n: u64) {
    let mut guard = 0;
    while s.stats.detected < n && guard < 200_000 {
        s.run_one();
        guard += 1;
    }
}

fn mean_profile(s: &Sim) -> Vec<f32> {
    let n = s.stats.photons_in_mean.max(1) as f64;
    s.stats.mean_intensity.iter().map(|v| (v / n) as f32).collect()
}

/// Young: fringes at spacing λL/d emerge from the grains' phase comparison,
/// and match the exact Huygens–Fresnel integral.
#[test]
fn double_slit_interference() {
    let mut s = sim(Experiment::DoubleSlit, |l| { l.grains = 40_000; l.seed = 7; });
    run_until_detected(&mut s, 40);
    let lab = &s.lab;
    let sim_p = mean_profile(&s);
    let ref_p = double_slit_reference(lab);
    let expected = lab.lambda * (lab.screen_z - lab.wall_z) / lab.slit_sep;
    let sp_sim = fringe_spacing(lab, &sim_p).unwrap();
    let sp_ref = fringe_spacing(lab, &ref_p).unwrap();
    let corr = correlation(&sim_p, &ref_p);
    let vis = visibility(lab, &sim_p);
    let p_pass = s.stats.detected as f32 / s.stats.emitted as f32;
    println!("fringe spacing: grains {sp_sim:.2}, exact integral {sp_ref:.2}, λL/d {expected:.2}");
    println!("correlation with exact pattern {corr:.4}, visibility {vis:.3}, photons through the slits {:.1}%", 100.0 * p_pass);
    assert!((sp_sim / expected - 1.0).abs() < 0.03);
    assert!(corr > 0.97);
    assert!(vis > 0.9);
}

/// Marking which slit the photon used (orthogonal polarizers on the slits)
/// erases the fringes; parallel polarizers keep them.
#[test]
fn which_path_information_erases_fringes() {
    let mut v = [0.0f32; 2];
    for (i, wp) in [WhichPath::Parallel, WhichPath::Orthogonal].into_iter().enumerate() {
        let mut s = sim(Experiment::DoubleSlit, |l| { l.grains = 40_000; l.which_path = wp; l.seed = 11; });
        run_until_detected(&mut s, 30);
        v[i] = visibility(&s.lab, &mean_profile(&s));
        println!("{wp:?}: visibility {:.3}", v[i]);
    }
    assert!(v[0] > 0.9);
    assert!(v[1] < 0.15);
}

/// One photon at a time: each lands at a single point, yet the dots
/// build up the interference pattern.
#[test]
fn single_photons_build_the_pattern() {
    let mut s = sim(Experiment::DoubleSlit, |l| { l.grains = 8_000; l.seed = 5; });
    run_until_detected(&mut s, 1500);
    let group = 8;
    let coarse = |v: Vec<f32>| v.chunks(group).map(|c| c.iter().sum::<f32>()).collect::<Vec<f32>>();
    let hist = coarse(s.stats.hist.iter().map(|&h| h as f32).collect());
    let reference = coarse(double_slit_reference(&s.lab));
    let corr = correlation(&hist, &reference);
    println!("{} single clicks vs exact pattern: correlation {corr:.3}", s.stats.detected);
    assert!(corr > 0.9);
}

/// Malus: a single photon crosses a polarizer at angle θ with probability
/// cos²θ. Three polarizers 0°/45°/90° let 1/4 through (with none in the
/// middle: zero).
#[test]
fn malus_law_and_three_polarizers() {
    let n = 6000;
    let run = |f: &dyn Fn(&mut Lab)| {
        let mut s = sim(Experiment::Polarizers, |l| { l.grains = 50; l.seed = 13; f(l); });
        for _ in 0..n { s.run_one(); }
        (s.stats.detected as f32 / n as f32, malus_theory(&s.lab))
    };
    let tol = 4.0 * (0.25f32 / n as f32).sqrt();

    for deg in [0.0f32, 30.0, 60.0, 90.0] {
        let (m, t) = run(&|l: &mut Lab| { l.polarizers[0].angle = deg * DEG; });
        println!("θ = {deg:4}°: passed {m:.3}   cos²θ = {t:.3}");
        assert!((m - t).abs() < tol.max(1e-6));
    }
    let three = |mid: bool| run(&move |l: &mut Lab| {
        l.polarizers[0].angle = 0.0;
        l.polarizers[1] = PolarizerCfg { z: 60.0, angle: 45.0 * DEG, enabled: mid };
        l.polarizers[2] = PolarizerCfg { z: 100.0, angle: 90.0 * DEG, enabled: true };
    });
    let (m0, _) = three(false);
    let (m1, t1) = three(true);
    println!("0°/90°: {m0:.3} (theory 0)   0°/45°/90°: {m1:.3} (theory {t1:.3})");
    assert_eq!(m0, 0.0);
    assert!((m1 - t1).abs() < tol);
}

fn bell(hidden: bool, pairs: usize) -> Sim {
    let mut s = sim(Experiment::Bell, |l| { l.grains = 20; l.bell.hidden_variables = hidden; l.seed = 17; });
    for _ in 0..pairs { s.run_one(); }
    s
}

/// Entanglement, the real test: CHSH. Quantum mechanics predicts S = 2√2;
/// any local model is bounded by 2. Each side alone stays 50/50 whatever the
/// other side measures: no faster-than-light signalling.
#[test]
fn bell_chsh_violation_and_no_signalling() {
    let q = bell(false, 40_000);
    let c = bell(true, 40_000);
    let sq = q.stats.bell.chsh_s().unwrap();
    let sc = c.stats.bell.chsh_s().unwrap();
    let (m0, m1) = (q.stats.bell.bob_plus_given_a(0).unwrap(), q.stats.bell.bob_plus_given_a(1).unwrap());
    println!("CHSH  quantum S = {sq:.3}  (theory 2√2 = 2.828)   local model S = {sc:.3}  (bound 2)");
    println!("no signalling: P(Bob +) = {m0:.3} if Alice uses 0°, {m1:.3} if 45°");
    for i in 0..2 { for j in 0..2 {
        let e = BellStats::correlation(&q.stats.bell.chsh[i][j]).unwrap();
        let t = correlation_quantum(CHSH_A[i], CHSH_B[j]);
        let l = correlation_local(CHSH_A[i], CHSH_B[j]);
        println!("  E({:4.1}°, {:4.1}°) = {e:+.3}   quantum {t:+.3}   local {l:+.3}", CHSH_A[i] / DEG, CHSH_B[j] / DEG);
    }}
    assert!((sq - 2.0 * 2f32.sqrt()).abs() < 0.06);
    assert!(sc <= 2.06);
    assert!((m0 - m1).abs() < 0.02);
    assert!((m0 - 0.5).abs() < 0.02);
}

/// Measuring Alice fixes Bob's polarization at that instant, while Bob is
/// still in flight and far away.
#[test]
fn collapse_is_instantaneous() {
    let mut s = sim(Experiment::Bell, |l| {
        l.grains = 200; l.bell.random_chsh = false; l.bell.angle_a = 30.0 * DEG; l.seed = 19;
    });
    s.launch();
    assert!(s.photons.iter().all(|p| !p.pol_defined));
    // Alice's splitter is at 40, Bob's at 70: step until Alice is measured
    while s.photons.iter().any(|p| !p.pol_defined) { s.step(1.0); }
    let a = s.photons.iter().find(|p| p.side == 0).unwrap();
    let b = s.photons.iter().find(|p| p.side == 1).unwrap();
    let bz = b.grains.iter().map(|g| g.pos.z).sum::<f32>() / b.grains.len() as f32;
    let (ja, jb) = (a.grains[0].jones, b.grains[0].jones);
    let overlap = (ja[0].conj() * jb[0] + ja[1].conj() * jb[1]).norm2();
    println!("Alice measured at 30°; Bob, {:.0} units away and not yet measured, now has |⟨A|B⟩|² = {overlap:.3}", 40.0 - bz);
    assert!(bz > -70.0, "Bob must not have reached his splitter yet");
    assert!((overlap - 1.0).abs() < 1e-4);
}


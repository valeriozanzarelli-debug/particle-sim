//! Validation of the base rule v2 in 3-D (`base3d.rs`) against the Python
//! reference (lab/base_exp3d.py, numbers in docs/REGOLA_BASE.md §4, rule v2).
//! Same geometries and tolerances as the Python scripts (written there before
//! the runs); fewer photons, so the statistical tolerance max(3σ, …) is wider.
//! Seeds are fixed; the random streams are not numpy's, so the numbers agree
//! with the Python ones within statistics, not bit for bit.
//!
//!   cargo test --bin fotone base3d -- --nocapture --test-threads=1

use std::sync::Arc;
use std::time::Instant;

use rayon::prelude::*;

use crate::base3d::*;
use crate::pilot::{hash2, uniform, Z};

const LAM: f64 = 1.0;

fn src_obj(name: &str, z: f64, mask: Mask, w0: f64, zw: f64) -> Obj {
    let mut o = Obj::new(name, Kind::Source, z, mask);
    o.src = Some(SrcBeam { w0, w0y: None, zw, jones: [Z::new(1.0, 0.0), Z::ZERO] });
    o
}

fn run_photons(b: &Arc<Bench3>, seed: u64, n: u64) -> Vec<Ball3> {
    (0..n).into_par_iter().map(|id| {
        let mut p = Ball3::new(b.clone(), seed, id);
        p.run();
        p.trail.clear();
        p
    }).collect()
}

/// Standard deviation of an estimator over bootstrap resamples.
fn bootstrap(x: &[f64], nb: usize, est: impl Fn(&[f64]) -> Option<f64> + Sync) -> f64 {
    let v: Vec<f64> = (0..nb).into_par_iter().filter_map(|k| {
        let key = hash2(0xB007, k as u64);
        let s: Vec<f64> = (0..x.len()).map(|i| x[((uniform(key, i as u64) * x.len() as f64) as usize).min(x.len() - 1)]).collect();
        est(&s)
    }).collect();
    let m = v.iter().sum::<f64>() / v.len() as f64;
    (v.iter().map(|q| (q - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
}

/// The boundary values of a far plane read by FFT (lattice convolution) are
/// the direct sum of the wavelets (Python: `Scene3D.field` on the points).
#[test]
fn base3d_lattice_fft_equals_direct_sum() {
    let mut objs = vec![src_obj("source", -20.0, Mask::Full, 3.0, 20.0)];
    let mut pol = Obj::new("polarizer", Kind::Polarizer(0), 20.0, Mask::Full);
    let a = [Z::new(0.6f64.cos(), 0.0), Z::new(0.6f64.sin(), 0.0)];
    pol.t = Some(projector(a));
    objs.push(pol);
    objs.push(Obj::new("screen", Kind::Screen, 60.0, Mask::Full));
    let t = Instant::now();
    let b = Bench3::build("t".into(), LAM, objs, [-150.0, 150.0, -150.0, 150.0, -21.0, 61.0]);
    println!("bench built in {:.2} s: {}", t.elapsed().as_secs_f32(),
        b.ball_counts().iter().map(|c| format!("{} {}×{}", c.obj, c.n0, c.ndir)).collect::<Vec<_>>().join(", "));
    let (x0, y0, step, nx, ny) = (-7.3, -6.1, 0.0625, 61, 47);
    let fft = b.lattice_tangential(0, 20.0, x0, y0, step, nx, ny);
    let mut worst = 0.0f64;
    let mut scale = 0.0f64;
    for q in (0..nx * ny).step_by(37) {
        let (e, _) = b.field([x0 + (q % nx) as f64 * step, y0 + (q / nx) as f64 * step, 20.0], 0);
        for c in 0..2 {
            worst = worst.max((fft[q][c] - e[c]).abs());
            scale = scale.max(e[c].abs());
        }
    }
    println!("FFT lattice vs direct sum: max |ΔE| = {worst:.2e} (max |E| = {scale:.3})");
    assert!(worst < 1e-9 * scale.max(1e-300));
}

/// Free Gaussian beam in 3-D (Python `exp_free`): the photon balls' width at
/// zR, 2zR, 5zR follows w(z), and the same photons expand as w(z)/w(1λ).
#[test]
fn base3d_free_beam_width() {
    let w0 = 3.0;
    let zr = std::f64::consts::PI * w0 * w0;
    let zs = [zr, 2.0 * zr, 5.0 * zr];
    let mut objs = vec![src_obj("source", 0.0, Mask::Full, w0, 0.0)];
    for z in zs { objs.push(Obj::new("plane", Kind::Plane, z, Mask::Full)); }
    objs.push(Obj::new("plane 1.0", Kind::Plane, 1.0, Mask::Full));
    objs.push(Obj::new("screen", Kind::Screen, 5.0 * zr + 1.0, Mask::Full));
    let b = Arc::new(Bench3::build("free".into(), LAM, objs, [-150.0, 150.0, -150.0, 150.0, -1.0, 5.0 * zr + 2.0]));
    let c = &b.ball_counts()[0];
    println!("free 3-D: {} origins × {} directions = {:.4e} balls (level 0), {:.2e} all levels [Python: 9 743 750 · 2.1e11]",
        c.n0, c.ndir, c.balls0, c.balls_all);
    let n = 2500u64;
    let t = Instant::now();
    let ph = run_photons(&b, 12, n);
    let secs = t.elapsed().as_secs_f32();
    let clicks = ph.iter().filter(|p| p.end == Some(End::Click)).count();
    println!("{n} photons, {clicks} clicks [{secs:.1} s]");
    let width = |k: usize| -> (f64, usize) {
        let pts: Vec<[f64; 3]> = ph.iter().flat_map(|p| p.rec.iter().filter(|r| r.0 == k).map(|r| r.1)).collect();
        let m = pts.len() as f64;
        let var = |a: usize| {
            let mu = pts.iter().map(|p| p[a]).sum::<f64>() / m;
            pts.iter().map(|p| (p[a] - mu).powi(2)).sum::<f64>() / m
        };
        (2.0 * (0.5 * (var(0) + var(1))).sqrt(), pts.len())
    };
    // Python v2 (7000 photons, 6982 clicks): widths and same-photon expansions
    let py_w = [4.3102, 6.8208, 15.5621];
    let py_r = [1.4154, 2.2398, 5.1103];
    let py_n = 6982.0f64;
    let (w1, _) = width(4);
    let mut ok = true;
    for (i, &z) in zs.iter().enumerate() {
        let (wm, cnt) = width(1 + i);
        let wth = w0 * (1.0 + (z / zr).powi(2)).sqrt();
        let sig = 1.0 / (4.0 * cnt as f64).sqrt();
        let err = wm / wth - 1.0;
        let pass = err.abs() <= (3.0 * sig).max(0.02);
        // agreement with the Python number within the two statistical errors
        let sp = 1.0 / (4.0 * py_n).sqrt();
        let agree = (wm / py_w[i] - 1.0).abs() <= 3.0 * (sig * sig + sp * sp).sqrt();
        let rth = (1.0 + (z / zr).powi(2)).sqrt() / (1.0 + (1.0 / zr).powi(2)).sqrt();
        let rpass = (wm / w1 / rth - 1.0).abs() <= 0.02;
        println!("  z = {:.0} zR: width {wm:.4} theory {wth:.4} err {:+.2}% (tol max(3σ={:.3}, 0.02)) {} | Python {:.4}: {} | \
                  expansion w/w(1λ) {:.4} theory {rth:.4} Python {:.4}: {}",
            z / zr, 100.0 * err, 3.0 * sig, if pass { "PASS" } else { "FAIL" }, py_w[i], if agree { "agrees" } else { "DIFFERS" },
            wm / w1, py_r[i], if rpass { "PASS" } else { "FAIL" });
        ok &= pass && agree && rpass;
    }
    assert!(ok);
    assert!(secs < 90.0, "too slow: {secs:.0} s");
}

/// Airy disk (Python `exp_airy`): D = 6 λ, screen at 200 λ, plane wave; the
/// first zero of the click density at 41.52 (L·tan(asin(1.2197 λ/D)));
/// Python v2: 41.63 ± 0.24.
#[test]
fn base3d_airy_first_zero() {
    let (d, la) = (6.0, 200.0);
    let objs = vec![
        src_obj("aperture", 0.0, Mask::CircIn(0.0, 0.0, d / 2.0), f64::INFINITY, 0.0),
        Obj::new("screen", Kind::Screen, la, Mask::Full),
    ];
    let b = Arc::new(Bench3::build("airy".into(), LAM, objs, [-400.0, 400.0, -400.0, 400.0, -1.0, la + 1.0]));
    let c = &b.ball_counts()[0];
    println!("airy: {} origins × {} directions = {:.4e} balls (level 0) [Python 219 252 · 4.8e9]", c.n0, c.ndir, c.balls0);
    let n: u64 = std::env::var("BASE3D_N").ok().and_then(|v| v.parse().ok()).unwrap_or(5000);
    let t = Instant::now();
    let ph = run_photons(&b, 61, n);
    let secs = t.elapsed().as_secs_f32();
    let r: Vec<f64> = ph.iter().filter(|p| p.end == Some(End::Click)).map(|p| p.x[0].hypot(p.x[1])).collect();
    let lost = ph.iter().filter(|p| p.end == Some(End::SourceCoin)).count();
    let rf = airy_first_zero(LAM, d, la);
    let est = |x: &[f64]| min_position(x, rf, 0.3 * rf, |q| q);
    let rm = est(&r).expect("a minimum");
    let sr = bootstrap(&r, 200, est);
    let tol = (3.0 * sr).max(0.02 * rf);
    let inside = r.iter().filter(|&&q| q < 2.2 * rf).count() as f64;
    let frac = r.iter().filter(|&&q| q < rf).count() as f64 / inside;
    println!("{n} photons, {} clicks, {lost} removed by the source coin [{secs:.1} s]", r.len());
    println!("first zero {rm:.2} ± {sr:.2}  theory {rf:.2} (exact RS 41.52)  tol {tol:.2}  err {:+.2}%  | Python 41.63 ± 0.24",
        100.0 * (rm / rf - 1.0));
    println!("fraction inside the first zero, of the clicks within 2.2 r1: {frac:.3}");
    assert!((rm - rf).abs() <= tol);
    assert!((rm - 41.63).abs() <= 3.0 * (sr * sr + 0.24 * 0.24).sqrt(), "differs from the Python number");
    assert!(secs < 90.0, "too slow: {secs:.0} s");
}

/// The fast sin/cos of the wavelet sums agrees with the library's.
#[test]
fn base3d_fast_sincos() {
    let mut worst = 0.0f64;
    for i in 0..2_000_000 {
        let x = -50.0 + 1e-4 * i as f64 * 37.0;        // −50 … 7350: k·R in every experiment
        let (s, c) = sincos(x);
        let (s0, c0) = x.sin_cos();
        worst = worst.max((s - s0).abs()).max((c - c0).abs());
    }
    println!("max |fast − std| = {worst:.2e}");
    assert!(worst < 1e-14);
}

/// Cost of the guidance sums (diagnostic, prints only).
#[test]
fn base3d_cost_of_a_field_sum() {
    use crate::physics::Experiment as E;
    for e in [E::FreeBeam, E::DoubleSlit, E::Polarizers, E::Airy, E::Arago] {
        let t = Instant::now();
        let b = Bench3::from_lab(&crate::physics::Lab::default_for(e));
        let c = b.ball_counts();
        println!("{e:?}: built in {:.2} s; {}; balls per photon {:.3e} (level 0), {:.2e} (all levels)", t.elapsed().as_secs_f32(),
            c.iter().map(|q| format!("{} {}x{}", q.obj, q.n0, q.ndir)).collect::<Vec<_>>().join(", "),
            c.iter().map(|q| q.balls0).sum::<f64>(), c.iter().map(|q| q.balls_all).sum::<f64>());
    }
    let lab = crate::physics::Lab::default_for(crate::physics::Experiment::DoubleSlit);
    let t = Instant::now();
    let b = Bench3::from_lab(&lab);
    println!("double slit bench built in {:.2} s: {}", t.elapsed().as_secs_f32(),
        b.ball_counts().iter().map(|c| format!("{} {}×{}", c.obj, c.n0, c.ndir)).collect::<Vec<_>>().join(", "));
    let n0 = b.fans[0].l0.len();
    let t = Instant::now();
    let mut acc = 0.0;
    let reps = 200;
    for i in 0..reps {
        let (j, _) = b.current([0.1 * i as f64 % 3.0, 0.5, -20.0], Some(0));
        acc += j[2];
    }
    let dt = t.elapsed().as_secs_f64() / reps as f64;
    let (e1, d1) = b.field([1.3, -0.7, -17.0], 0);
    let (e2, d2) = b.field_l0_scalar([1.3, -0.7, -17.0], 0);
    let scale = e2.iter().chain(d2.iter()).map(|v| v.abs()).fold(0.0, f64::max);
    let diff = e1.iter().zip(&e2).chain(d1.iter().zip(&d2)).map(|(a, b)| (*a - *b).abs()).fold(0.0, f64::max);
    println!("4-lane far-field sum vs scalar add3: max diff {diff:.2e} (scale {scale:.3})");
    let t = Instant::now();
    for i in 0..reps { let _ = std::hint::black_box(b.field([0.01 * i as f64, 0.5, -20.0], 0)); }
    let t1 = t.elapsed().as_secs_f64() / reps as f64;
    let t = Instant::now();
    for i in 0..reps { let _ = std::hint::black_box(b.field_l0_scalar([0.01 * i as f64, 0.5, -20.0], 0)); }
    let t2 = t.elapsed().as_secs_f64() / reps as f64;
    println!("per origin: field {:.1} ns, scalar add3 {:.1} ns", 1e9 * t1 / n0 as f64, 1e9 * t2 / n0 as f64);
    assert!(diff < 1e-12 * scale);
    println!("source fan: {n0} origins, {:.3} ms per field evaluation (all threads), {:.1} ns per origin [{acc:.3e}]", 1e3 * dt, 1e9 * dt / n0 as f64);
    let t = Instant::now();
    let mut p = Ball3::new(Arc::new(b), 1, 3);
    p.run();
    println!("one photon: {:?} in {:.0} ms, {} steps", p.end, 1e3 * t.elapsed().as_secs_f64(), p.steps);
}

/// The program's own polarizer laboratory (P1 at 30°, light polarized at 0°):
/// the local coin lets cos²30° = 0.75 of the photons through, whole.
#[test]
fn base3d_malus_in_the_program_lab() {
    let lab = crate::physics::Lab::default_for(crate::physics::Experiment::Polarizers);
    let b = Arc::new(Bench3::from_lab(&lab));
    let n = 800u64;
    let t = Instant::now();
    let ph = run_photons(&b, 5, n);
    let passed = ph.iter().filter(|p| p.end == Some(End::Click)).count() as f64;
    let coin = ph.iter().filter(|p| matches!(p.end, Some(End::Coin(i)) if b.objs[i].kind == Kind::Polarizer(0))).count() as f64;
    let arrived = passed + coin;
    let p = passed / arrived;
    let th = (30f64.to_radians()).cos().powi(2);
    let sig = (th * (1.0 - th) / arrived).sqrt();
    println!("P(pass) = {p:.4}  cos²30° = {th:.4}  (σ {sig:.4}, {arrived} photons at the polarizer) [{:.1} s]", t.elapsed().as_secs_f32());
    assert!((p - th).abs() <= (3.0 * sig).max(0.01));
}

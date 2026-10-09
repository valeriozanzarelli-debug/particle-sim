//! Validation of the pilot rule (`pilot.rs`) — the Rust port of
//! `lab/photon_rule.py` / `lab/validate.py`. Reference numbers of the Python
//! lab: docs/REGOLA_FOTONE.md §5–§6. All runs are seeded.
//!
//!   cargo test --bin fotone pilot -- --nocapture --test-threads=1

use std::sync::Arc;
use std::time::Instant;

use rayon::prelude::*;

use crate::physics::*;
use crate::pilot::*;

/// scipy.special.hankel1(0|1, u), generated once:
/// (u, Re H0, Im H0, Re H1, Im H1)
const HANKEL_SCIPY: [(f64, f64, f64, f64, f64); 32] = [
    (0.01, 0.9999750001562496, -3.0054556370836463, 0.004999937500264314, -63.67859628206066),
    (0.05, 0.9993750976494685, -1.97931100081721, 0.02499218831376048, -12.789855171174974),
    (0.1, 0.9975015620660401, -1.5342386513503667, 0.04993752603624239, -6.4589510947020266),
    (0.25, 0.9844359292958527, -0.9315730249300587, 0.12402597732272712, -2.7041052293152825),
    (0.5, 0.9384698072408126, -0.44451873350670656, 0.24226845767487393, -1.4714723926702433),
    (0.75, 0.8642422751666486, -0.1371727693857724, 0.3492436021748622, -1.0375945507692854),
    (1.0, 0.7651976865579664, 0.088256964215677, 0.44005058574493355, -0.7812128213002889),
    (1.5, 0.5118276717359181, 0.38244892379775886, 0.5579365079100995, -0.4123086269739113),
    (2.0, 0.22389077914123567, 0.5103756726497453, 0.5767248077568734, -0.1070324315409375),
    (2.5, -0.048383776468197956, 0.49807035961523194, 0.49709410246427405, 0.1459181379667858),
    (2.9, -0.22431154579196805, 0.40791176923625017, 0.375427481813096, 0.2959400546076747),
    (2.999, -0.25971270943712194, 0.3771745500619621, 0.33943194156540585, 0.3244056105372726),
    (3.0, -0.2600519549019336, 0.3768500100127904, 0.3390589585259364, 0.32467442479180014),
    (3.3, -0.3442962603988846, 0.26909199505453385, 0.2206634529852411, 0.3878529310237099),
    (4.0, -0.3971498098638474, -0.016940739325064957, -0.06604332802354913, 0.3979257105571),
    (5.0, -0.1775967713143384, -0.3085176252490338, -0.3275791375914653, 0.14786314339122694),
    (6.0, 0.150645257250997, -0.2881946839815792, -0.27668385812756563, -0.17501034430039833),
    (7.5, 0.26633965788037844, 0.11731328614820861, 0.1352484275797055, -0.25912851048611635),
    (9.0, -0.09033361118287618, 0.24993669828502468, 0.2453117865733253, 0.10431457519671594),
    (10.5, -0.23664819446234722, -0.06753037249787641, -0.07885001422733152, 0.23370422835726865),
    (11.5, -0.06765394811166525, -0.22523211169118795, -0.22837862066532355, 0.057942547143000823),
    (11.999, 0.04746583057345669, -0.22529430160059621, -0.223513306194832, -0.05687870171368411),
    (12.0, 0.04768931079683358, -0.22523731263436142, -0.2234471044906276, -0.057099218260896575),
    (12.5, 0.14688405470042115, -0.1712143068446693, -0.16548380461475973, -0.15383825653750122),
    (14.0, 0.17107347611045867, 0.12719256858218367, 0.13337515469879324, -0.16664484185617232),
    (17.0, -0.16985425215118363, -0.09263719844232368, -0.09766849275778064, 0.16720503607723375),
    (20.0, 0.16702466434058325, 0.06264059680938383, 0.06683312417585006, -0.16551161436252138),
    (30.0, -0.08636798358104025, -0.11729573168666407, -0.11875106261662301, 0.08442557066174726),
    (50.0, 0.05581232766925183, -0.0980649954700771, -0.09751182812517517, -0.05679566856201478),
    (100.0, 0.019985850304223125, -0.07724431336508318, -0.07714535201411218, -0.0203723120027598),
    (250.0, -0.02605337342520425, -0.04321684544036628, -0.04326903841033076, 0.025966992185484594),
    (1000.0, 0.024786686152420176, 0.004715917977622814, 0.004728311907089527, -0.02478433129235178),
];

/// H0 and H1 (asymptotic for u ≥ 12, Abramowitz–Stegun below) against scipy.
#[test]
fn pilot_hankel_matches_scipy() {
    let (mut worst0, mut worst1, mut at0, mut at1) = (0.0f64, 0.0f64, 0.0, 0.0);
    for &(u, r0, i0, r1, i1) in &HANKEL_SCIPY {
        let (h0, h1) = hankel01(u);
        let e0 = (h0 - Z::new(r0, i0)).abs() / Z::new(r0, i0).abs();
        let e1 = (h1 - Z::new(r1, i1)).abs() / Z::new(r1, i1).abs();
        if e0 > worst0 { worst0 = e0; at0 = u; }
        if e1 > worst1 { worst1 = e1; at1 = u; }
    }
    println!("max relative error vs scipy: H0 {worst0:.2e} (u = {at0}), H1 {worst1:.2e} (u = {at1})");
    assert!(worst0 < 1e-6 && worst1 < 1e-6);
}

fn lab(e: Experiment, f: impl FnOnce(&mut Lab)) -> Lab {
    let mut l = Lab::default_for(e);
    l.rule = Rule::Pilot;      // the program opens with the base 3-D rule; these tests are the 2-D pilot's
    f(&mut l);
    l
}

/// Photons `first..first+n` of `seed`, each run alone (headless), in parallel.
fn run_many(b: &Arc<Bench>, seed: u64, first: u64, n: u64) -> Vec<Outcome> {
    (first..first + n).into_par_iter().map(|id| run_headless(b, seed, id)).collect()
}

/// 1 — an isolated photon ball stays on a straight ray.
///  (a) it hears nobody: it keeps its own direction;
///  (b) its wave is ONE origin (N = 1): it moves along the wavelet's ray.
#[test]
fn pilot_isolated_ball_is_straight() {
    let l = lab(Experiment::DoubleSlit, |_| {});
    let b = Arc::new(Bench::new(&l));

    let mut p = PilotPhoton::new(b.clone(), 1, 0);
    (p.started, p.fan, p.x, p.z, p.u) = (true, 1, 4.3, 0.1, 0.07);
    let deaf = |_: usize, _: usize, _: usize, _: f64| false;
    p.advance(120.0, Some(&deaf), false, None);
    let dev_a = (p.x - (4.3 + 0.07 * (p.z - 0.1))).abs();

    let mut one = (*b).clone();
    one.fans[1].xs = vec![4.0];
    one.fans[1].amp = vec![[Z::new(1.0, 0.0), Z::ZERO]];
    let one = Arc::new(one);
    let mut q = PilotPhoton::new(one, 1, 0);
    (q.started, q.fan, q.x, q.z, q.u) = (true, 1, 4.5, 1.0, 0.0);
    q.advance(120.0, None, false, None);
    let dev_b = (q.x - (4.0 + 0.5 * q.z)).abs();

    println!("no neighbours: |x − straight ray| = {dev_a:.1e} λ after Δz = {:.1}", p.z - 0.1);
    println!("one origin:    |x − wavelet ray| = {dev_b:.1e} λ after Δz = {:.1}", q.z - 1.0);
    assert!(p.outcome.is_none() && q.outcome.is_none());
    assert!(dev_a < 1e-6 && dev_b < 1e-6);
}

/// 2 — free Gaussian beam: the photon balls' landing width follows w(z).
#[test]
fn pilot_free_beam_width() {
    let l = lab(Experiment::FreeBeam, |_| {});
    let b = Arc::new(Bench::new(&l));
    let zr = l.rayleigh() as f64;
    let zs = [0.5, zr, 2.0 * zr, 5.0 * zr];      // 0.5 = just behind the source (δ)
    let n = 20_000u64;
    let t = Instant::now();
    let runs: Vec<Option<[f64; 4]>> = (0..n).into_par_iter().map(|id| {
        let mut p = PilotPhoton::new(b.clone(), 11, id);
        let mut out = [0.0; 4];
        for (o, &z) in out.iter_mut().zip(&zs) {
            p.advance(f64::INFINITY, None, false, Some(z));
            if p.outcome.is_some() { return None; }
            *o = p.x;
        }
        Some(out)
    }).collect();
    let ok: Vec<[f64; 4]> = runs.iter().flatten().copied().collect();
    println!("free beam w0 = {}, z_R = {zr:.2}: {n} photons, lost at the source {} [{:.1} s]",
        l.waist, n as usize - ok.len(), t.elapsed().as_secs_f32());
    for (i, &z) in zs.iter().enumerate() {
        let xs: Vec<f64> = ok.iter().map(|o| o[i]).collect();
        let m = xs.iter().sum::<f64>() / xs.len() as f64;
        let w = 2.0 * (xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / xs.len() as f64).sqrt();
        let th = l.waist as f64 * (1.0 + (z / zr).powi(2)).sqrt();
        let err = w / th - 1.0;
        println!("  z = {z:6.1}: width {w:6.3}  theory {th:6.3}  error {:+.2}%  (stat. ±{:.1}%)",
            100.0 * err, 100.0 / (2.0 * xs.len() as f64).sqrt());
        assert!(err.abs() < 0.03);
    }
}

// ─── Double-slit scoring (same as lab/optics.py) ──────────────────────────────

const BINS: usize = 240;
const HALF: f64 = 60.0;
const COARSE: usize = 8;

fn xs_bins() -> Vec<f64> { (0..BINS).map(|i| -HALF + (i as f64 + 0.5) * 2.0 * HALF / BINS as f64).collect() }

/// Exact 2-D RS flux j_z on the screen behind the slits, incident Gaussian
/// exp(−x²/w0²) on the wall (as `optics.screen_reference("flux")`).
fn ds_reference(l: &Lab, jones: [J2; 2]) -> Vec<f64> {
    let k = std::f64::consts::TAU / l.lambda as f64;
    let (c, h) = l.slits();
    let (w0, lz) = (l.waist as f64, (l.screen_z - l.wall_z) as f64);
    let m = 400;
    let ds = 2.0 * h as f64 / m as f64;
    let mut src: Vec<(f64, J2)> = Vec::new();
    for s in 0..2 {
        for i in 0..m {
            let x = c[s] as f64 - h as f64 + (i as f64 + 0.5) * ds;
            let w = (-(x * x) / (w0 * w0)).exp() * ds;
            src.push((x, [jones[s][0].s(w), jones[s][1].s(w)]));
        }
    }
    let sub = 8;
    let fine: Vec<f64> = (0..BINS * sub).into_par_iter().map(|i| {
        let x = -HALF + (i as f64 + 0.5) * 2.0 * HALF / (BINS * sub) as f64;
        let (mut u, mut uz) = ([Z::ZERO; 2], [Z::ZERO; 2]);
        for (s, a) in &src {
            let (g, _, gz) = kernel(k, x - s, lz);
            for p in 0..2 { u[p] += g * a[p]; uz[p] += gz * a[p]; }
        }
        (0..2).map(|p| (u[p].conj() * uz[p]).im).sum::<f64>()
    }).collect();
    let p: Vec<f64> = fine.chunks(sub).map(|c| c.iter().sum()).collect();
    let s: f64 = p.iter().sum();
    p.iter().map(|v| v / s).collect()
}

fn hist(clicks: &[f64]) -> Vec<f64> {
    let mut h = vec![0.0; BINS];
    for &x in clicks {
        let f = (x + HALF) / (2.0 * HALF);
        if (0.0..1.0).contains(&f) { h[(f * BINS as f64) as usize] += 1.0; }
    }
    h
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut ab, mut aa, mut bb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) { ab += (x - ma) * (y - mb); aa += (x - ma).powi(2); bb += (y - mb).powi(2); }
    ab / (aa * bb).sqrt()
}

fn coarse_corr(clicks: &[f64], r: &[f64]) -> f64 {
    let c = |v: &[f64]| v.chunks(COARSE).map(|c| c.iter().sum()).collect::<Vec<f64>>();
    pearson(&c(&hist(clicks)), &c(r))
}

/// Relative error of the click fraction in each bright fringe (|x| < 45).
fn fringe_errors(clicks: &[f64], r: &[f64], period: f64) -> (Vec<f64>, f64, f64) {
    let xs = xs_bins();
    let fine: Vec<f64> = (0..BINS).map(|i| {
        let g = |j: isize| if j < 0 || j >= BINS as isize { 0.0 } else { r[j as usize] };
        (g(i as isize - 1) + g(i as isize) + g(i as isize + 1)) / 3.0
    }).collect();
    let mx = fine.iter().cloned().fold(0.0, f64::max);
    let n = clicks.iter().filter(|x| x.abs() < HALF).count().max(1) as f64;
    let mut errs = Vec::new();
    for i in 1..BINS - 1 {
        if fine[i] >= fine[i - 1] && fine[i] > fine[i + 1] && xs[i].abs() < 45.0 && fine[i] > 0.05 * mx {
            let (a, b) = (xs[i] - period / 4.0, xs[i] + period / 4.0);
            let pth: f64 = (0..BINS).filter(|&j| xs[j] >= a && xs[j] < b).map(|j| r[j]).sum();
            let pm = clicks.iter().filter(|&&x| x >= a && x < b).count() as f64 / n;
            errs.push(pm / pth - 1.0);
        }
    }
    let rms = (errs.iter().map(|e| e * e).sum::<f64>() / errs.len() as f64).sqrt();
    let max = errs.iter().map(|e| e.abs()).fold(0.0, f64::max);
    (errs, rms, max)
}

/// (I0 − Id)/(I0 + Id): counts within ±2.5 of the centre vs the first dark fringes.
fn visibility_counts(clicks: &[f64], period: f64) -> f64 {
    let cnt = |c: f64| clicks.iter().filter(|&&x| (x - c).abs() < 2.5).count() as f64;
    let (i0, id) = (cnt(0.0), 0.5 * (cnt(-period / 2.0) + cnt(period / 2.0)));
    (i0 - id) / (i0 + id).max(1e-30)
}

fn visibility_prob(r: &[f64], period: f64) -> f64 {
    let xs = xs_bins();
    let cnt = |c: f64| (0..BINS).filter(|&i| (xs[i] - c).abs() < 2.5).map(|i| r[i]).sum::<f64>();
    let (i0, id) = (cnt(0.0), 0.5 * (cnt(-period / 2.0) + cnt(period / 2.0)));
    (i0 - id) / (i0 + id)
}

struct SlitRun { landings: Vec<f64>, clicks: usize, emitted: usize, wall: usize, source: usize, secs: f32 }

/// Emit photons (headless, in parallel) until `want` clicks on the screen.
fn run_slits(b: &Arc<Bench>, seed: u64, want: usize) -> SlitRun {
    let t = Instant::now();
    let mut r = SlitRun { landings: Vec::new(), clicks: 0, emitted: 0, wall: 0, source: 0, secs: 0.0 };
    while r.clicks < want {
        let batch = 2048;
        for o in run_many(b, seed, r.emitted as u64, batch) {
            match o {
                Outcome::Click(x) => { r.landings.push(x); r.clicks += 1; }
                Outcome::Missed(x) => r.landings.push(x),
                Outcome::Wall => r.wall += 1,
                Outcome::Source => r.source += 1,
                Outcome::Polarizer(_) => unreachable!(),
            }
        }
        r.emitted += batch as usize;
    }
    r.secs = t.elapsed().as_secs_f32();
    r
}

fn plain() -> [J2; 2] {
    let j = Polarization::D.jones();
    let z = [Z::from_c(j[0]), Z::from_c(j[1])];
    [z, z]
}

/// 3 — double slit: clicks = landing points of the photon balls.
/// Python (docs §5, variant I, 5858 clicks): corr 0.9949, fringe rms 3.9 %,
/// max 4.7 %, absorbed 78.9 %.
#[test]
fn pilot_double_slit() {
    let l = lab(Experiment::DoubleSlit, |_| {});
    let b = Arc::new(Bench::new(&l));
    let reference = ds_reference(&l, plain());
    let period = (l.lambda * (l.screen_z - l.wall_z) / l.slit_sep) as f64;
    let r = run_slits(&b, 7, 3000);
    let corr = coarse_corr(&r.landings, &reference);
    let (errs, rms, max) = fringe_errors(&r.landings, &reference, period);
    let absorbed = r.wall as f64 / r.emitted as f64;
    let own: Vec<f64> = b.screen_flux.iter().map(|&v| v as f64).collect();
    println!("double slit: {} photons, {} clicks [{:.1} s]", r.emitted, r.clicks, r.secs);
    println!("  corr (8-bin) vs exact RS flux {corr:.4}  (vs the program's own wave {:.4}; wave vs exact {:.5})",
        coarse_corr(&r.landings, &own), pearson(&own, &reference));
    println!("  fringe error rms {:.1}%  max {:.1}%  per fringe {:?}", 100.0 * rms, 100.0 * max,
        errs.iter().map(|e| (100.0 * e).round() as i32).collect::<Vec<_>>());
    println!("  absorbed at the wall {:.1}% (Python 78.9%), lost at the source {} ({:.2}%)",
        100.0 * absorbed, r.source, 100.0 * r.source as f64 / r.emitted as f64);
    println!("  visibility of the clicks {:.3}, of the exact flux {:.3}",
        visibility_counts(&r.landings, period), visibility_prob(&reference, period));
    assert!(corr >= 0.98);
    assert!(rms <= 0.08);
    assert!((absorbed - 0.789).abs() <= 0.02);
}

/// 4 — which-path polarizers on the slits (light at 45°): parallel keep the
/// fringes, orthogonal erase them. Python: V = 0.925 / −0.023.
#[test]
fn pilot_which_path() {
    let mut v = [0.0; 2];
    for (n, wp) in [WhichPath::Parallel, WhichPath::Orthogonal].into_iter().enumerate() {
        let l = lab(Experiment::DoubleSlit, |l| l.which_path = wp);
        let b = Arc::new(Bench::new(&l));
        let f = l.slit_filters();
        let din = plain()[0];
        let jones: [J2; 2] = std::array::from_fn(|s| {
            let fj = f[s].unwrap();
            let fz = [Z::from_c(fj[0]), Z::from_c(fj[1])];
            let c = fz[0].conj() * din[0] + fz[1].conj() * din[1];
            [fz[0] * c, fz[1] * c]
        });
        let reference = ds_reference(&l, jones);
        let period = (l.lambda * (l.screen_z - l.wall_z) / l.slit_sep) as f64;
        let r = run_slits(&b, 21 + n as u64, 1500);
        v[n] = visibility_counts(&r.landings, period);
        println!("which-path {wp:?}: {} photons, {} clicks, absorbed {:.1}% [{:.1} s]: visibility {:.3} (exact flux {:.3}), corr {:.4}",
            r.emitted, r.clicks, 100.0 * r.wall as f64 / r.emitted as f64, r.secs, v[n],
            visibility_prob(&reference, period), coarse_corr(&r.landings, &reference));
    }
    assert!(v[0] > 0.85);
    assert!(v[1] < 0.15);
}

/// 5 — Malus with an absorbing polarizer (local flux coin of the photon ball).
/// Python M1: 0.9962 / 0.4885 / 0.0000 at 0 / 45 / 90°.
#[test]
fn pilot_malus() {
    let n = 2000u64;
    for deg in [0.0f32, 45.0, 90.0] {
        let l = lab(Experiment::Polarizers, |l| {
            l.polarizers[0].angle = deg.to_radians();
            l.polarizers[0].enabled = true;
            l.polarizers[1].enabled = false;
            l.polarizers[2].enabled = false;
        });
        let b = Arc::new(Bench::new(&l));
        let t = Instant::now();
        let out = run_many(&b, 31, 0, n);
        let passed = out.iter().filter(|o| matches!(o, Outcome::Click(_) | Outcome::Missed(_))).count() as f64;
        let absorbed = out.iter().filter(|o| matches!(o, Outcome::Polarizer(0))).count() as f64;
        let lost = out.iter().filter(|o| matches!(o, Outcome::Source)).count();
        let arrived = passed + absorbed;
        let p = passed / arrived;
        let th = (deg.to_radians() as f64).cos().powi(2);
        let sigma = ((th * (1.0 - th)).max(1.0 / arrived) / arrived).sqrt();
        println!("θ = {deg:4.1}°: passed {p:.4}  cos²θ {th:.4}  (σ {sigma:.4}, {arrived} photons at the polarizer, {lost} lost at the source) [{:.1} s]",
            t.elapsed().as_secs_f32());
        assert!((p - th).abs() < (3.0 * sigma).max(0.01));
    }
}

/// The live run (wave balls stored, the photon ball checks which rays are
/// open) and the headless run (no balls) give the same photon, bit for bit.
#[test]
fn pilot_live_equals_headless() {
    let l = lab(Experiment::DoubleSlit, |l| l.grains = 6000);
    let run = |headless: bool| {
        let mut s = Sim::new(l.clone());
        s.headless = headless;
        let mut trail = Vec::new();
        while s.stats.emitted < 20 || s.in_flight() {
            if !s.in_flight() { s.launch(); }
            s.step(1.7);
            if let Some(pp) = s.photons.first().and_then(|p| p.pilot.as_ref()) { trail.push((pp.x, pp.z)); }
        }
        (trail, s.stats.clicks.clone(), s.stats.absorbed_wall)
    };
    let (a, ca, wa) = run(false);
    let (b, cb, wb) = run(true);
    println!("{} steps compared; clicks {:?}; absorbed {wa}", a.len(), ca.iter().map(|c| c.x).collect::<Vec<_>>());
    assert_eq!(a, b);
    assert_eq!(ca, cb);
    assert_eq!(wa, wb);
}

/// Reversible time with the pilot rule: rewinding and replaying gives the
/// same wave balls and the same photon ball, exactly.
#[test]
fn pilot_rewind_and_replay_are_exact() {
    use crate::timeline::{apply, History, Op};
    let l = lab(Experiment::DoubleSlit, |l| l.grains = 3000);
    let mut sim = Sim::new(l);
    let mut h = History::new(&sim);
    let snap = |s: &Sim| -> (Vec<[f32; 3]>, Vec<(u64, u64)>) {
        (s.photons.iter().flat_map(|p| p.grains.iter().map(|g| g.pos.to_array())).collect(),
         s.photons.iter().filter_map(|p| p.pilot.as_ref()).map(|pp| (pp.x.to_bits(), pp.z.to_bits())).collect())
    };
    let mut trail = Vec::new();
    for i in 0..300 {
        let op = if !sim.in_flight() { Op::Launch } else { Op::Step(0.7 + (i % 5) as f32 * 0.3) };
        apply(&mut sim, op);
        h.record(op, &sim);
        trail.push((h.cursor(), snap(&sim), sim.stats.emitted));
    }
    let mut shown = sim.clone();
    let base = h.base();
    for (frame, s, emitted) in trail.iter().rev().filter(|t| t.0 >= base) {
        h.move_to(&mut shown, *frame);
        assert!(&snap(&shown) == s, "frame {frame}");
        assert_eq!(shown.stats.emitted, *emitted);
    }
    for (frame, s, _) in trail.iter().filter(|t| t.0 >= base) {
        h.move_to(&mut shown, *frame);
        assert!(&snap(&shown) == s, "frame {frame}");
    }
    println!("{} frames replayed exactly ({} photons emitted)", trail.len() as u64 - base, sim.stats.emitted);
}

/// Cost of the animation: one 62 800-ball double-slit photon at speed ×1
/// (60 frames per second of light-clock). Simulation only (no drawing).
#[test]
fn pilot_frame_cost() {
    let l = lab(Experiment::DoubleSlit, |_| {});
    let mut s = Sim::new(l);
    let d = crate::timeline::base_speed(&s.lab) / 60.0;
    s.bench();
    // time a photon that crosses the slits and clicks (most die at the wall)
    let (mut n, mut worst, mut total, mut max_balls) = (0, 0.0f32, 0.0f32, 0);
    while s.stats.detected == 0 {
        s.launch();
        (n, worst, total, max_balls) = (0, 0.0, 0.0, 0);
        while s.in_flight() {
            let t1 = Instant::now();
            s.step(d);
            let dt = t1.elapsed().as_secs_f32();
            worst = worst.max(dt);
            total += dt;
            n += 1;
            max_balls = max_balls.max(s.photons.first().map_or(0, |p| p.grains.len()));
        }
    }
    println!("photon #{} clicked at x = {:.2}", s.stats.emitted, s.stats.clicks[0].x);
    let mean = total / n as f32;
    let b = s.bench();
    println!("fans: {}", b.fans.iter().map(|f| format!("{}×{}", f.xs.len(), f.nu)).collect::<Vec<_>>().join(" → "));
    println!("{n} frames of {d:.2} λ: mean {:.2} ms, worst {:.2} ms per step ({max_balls} wave balls at most)",
        1e3 * mean, 1e3 * worst);
    assert!(mean < 0.05, "a step must fit a 20 fps frame");
}

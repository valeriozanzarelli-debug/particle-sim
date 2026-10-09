//! Tests of the zig-zag binary: unitarity, the lattice dispersion relation, speed ≤ c, the Bloch matrix
//! against the walk itself, hydrogen from the radial walk vs the Python lab and the Dirac formula.

use crate::atom::{self, Hydrogen, ALPHA, MC2_EV};
use crate::cplx::C;
use crate::walk::*;
use crate::Rng;

fn random_state(n: usize, rng: &mut Rng) -> Vec<[C; 4]> {
    (0..n * n * n).map(|_| [0; 4].map(|_| C::new(rng.gauss(), rng.gauss()))).collect()
}

#[test]
fn unitarity_dirac_with_potential() {
    let mut rng = Rng::new(1);
    let n = 12;
    let mut d = Dirac3D::new(n, 0.3, 0.3);
    d.psi = random_state(n, &mut rng);
    let v: Vec<f64> = (0..n * n * n).map(|_| rng.uni() * 5.0).collect();
    d.set_potential(&v);
    let n0 = d.norm2();
    for _ in 0..50 { d.tick(); }
    let e = (d.norm2() - n0).abs() / n0;
    assert!(e < 1e-12, "norm drift {e:e}");
}

#[test]
fn unitarity_light() {
    let mut rng = Rng::new(2);
    let n = 12;
    let mut l = Light3D::new(n);
    for p in l.f.iter_mut() { *p = [0; 3].map(|_| C::new(rng.gauss(), rng.gauss())); }
    let n0 = l.norm2();
    for _ in 0..50 { l.tick(); }
    let e = (l.norm2() - n0).abs() / n0;
    assert!(e < 1e-12, "norm drift {e:e}");
}

#[test]
fn unitarity_radial() {
    let w = Radial::new(-1, 0.5, 300.0, ALPHA);
    let mut rng = Rng::new(3);
    let x: Vec<C> = (0..2 * w.n).map(|_| C::new(rng.gauss(), rng.gauss())).collect();
    let mut y = vec![C::ZERO; 2 * w.n];
    w.apply(&x, &mut y);
    let (a, b): (f64, f64) = (x.iter().map(|c| c.norm2()).sum(), y.iter().map(|c| c.norm2()).sum());
    assert!((a - b).abs() / a < 1e-12);
}

/// The Bloch matrix is what the walk does to a plane wave (two ticks, any spinor).
#[test]
fn bloch_matrix_is_the_walk() {
    let n = 8;
    let theta = 0.37;
    let kk = [1usize, 2, 3].map(|m| 2.0 * std::f64::consts::PI * m as f64 / n as f64);
    let spin = [C::new(0.3, -0.2), C::new(-0.5, 0.1), C::new(0.2, 0.7), C::new(0.4, 0.0)];
    let mut d = Dirac3D::new(n, 1.0, theta);
    for s in 0..n * n * n {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        let ph = C::expi(kk[0] * i as f64 + kk[1] * j as f64 + kk[2] * k as f64);
        d.psi[s] = spin.map(|c| c * ph);
    }
    d.tick();
    d.tick();
    let w = bloch_dirac(kk, theta);
    for q in 0..4 {
        let pred = (0..4).fold(C::ZERO, |acc, r| acc + w[q][r] * spin[r]);
        assert!((d.psi[0][q] - pred).abs() < 1e-12, "component {q}: {:?} vs {:?}", d.psi[0][q], pred);
    }
}

/// Lattice dispersion: along an axis cos ω = cos θ cos(k a) exactly; rest energy = m c² exactly;
/// E → √(p² + m²) with an error ∝ θ² (second-order rule).
#[test]
fn dispersion_relation() {
    for theta in [0.4, 0.2, 0.1] {
        for p in [0.0, 0.3, 1.0, 2.0] {
            let e = walk_energy(p, theta);
            let exact_lattice = (theta.cos() * (p * theta).cos()).acos() / theta;
            assert!((e - exact_lattice).abs() < 1e-9, "θ={theta} p={p}: {e} vs {exact_lattice}");
        }
        assert!((walk_energy(0.0, theta) - 1.0).abs() < 1e-12);
    }
    let err = |th: f64| (walk_energy(1.0, th) / 2f64.sqrt() - 1.0).abs();
    let r = err(0.2) / err(0.1);
    assert!((3.0..5.0).contains(&r), "error ratio {r}");
    assert!(err(0.1) < 2e-3);
}

/// The surfer never exceeds c, for any spinor, and the walk moves every eigen-part by exactly one cell.
#[test]
fn speed_le_c() {
    let mut rng = Rng::new(4);
    for _ in 0..100_000 {
        let s = [0; 4].map(|_| C::new(rng.gauss(), rng.gauss()));
        let v = surfer_velocity(&s);
        assert!((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() <= 1.0 + 1e-12);
    }
    // a right-handed ball with spin along +x moves exactly one cell along +x per tick (θ = 0): speed c.
    // Bloch view: a plane wave along x with that spinor gets the phase e^{−ika} per tick, exactly.
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let w = bloch_dirac([0.4, 0.0, 0.0], 0.0);
    // eigenvalue of the R spin-+x state = e^{-2i k a} (two ticks, one cell per tick)
    let v = [C::ZERO, C::ZERO, C::real(h), C::real(h)];
    let wv: Vec<C> = (0..4).map(|q| (0..4).fold(C::ZERO, |a, r| a + w[q][r] * v[r])).collect();
    let lam = C::expi(-0.8);
    for q in 0..4 { assert!((wv[q] - v[q] * lam).abs() < 1e-12); }
}

/// Hydrogen 1s on the θ = 0.2 lattice equals the Python lab value (lab/out/zz_hydrogen.log:
/// 13.7826491 eV) and the extrapolated fine structure equals Dirac within 1e-3.
#[test]
fn hydrogen_radial_walk() {
    let lv = atom::level(1, -1, 0.2);
    let e = lv.binding * MC2_EV;
    assert!((e / 13.782_649_1 - 1.0).abs() < 1e-7, "1s at θ=0.2: {e}");
    assert_eq!(atom::nodes(&lv.gf), 0);
    let y = |n, k| {
        let mut out = [0.0; 3];
        for (i, &th) in atom::THETAS.iter().enumerate() { out[i] = atom::level(n, k, th).binding * MC2_EV; }
        atom::extrapolate(out)
    };
    let fs = y(2, 1) - y(2, -2);
    let fs_d = (atom::dirac_formula(2, -2, ALPHA) - atom::dirac_formula(2, 1, ALPHA)) * MC2_EV;
    assert!((fs / fs_d - 1.0).abs() < 1e-3, "fine structure {fs} vs {fs_d}");
}

#[test]
fn spinor_harmonics_normalized() {
    // ∫ |Ω_κm|² dΩ = 1 for the states drawn
    for (kappa, tm) in [(-1, 1), (1, 1), (-2, 1), (2, 1), (-2, 3)] {
        let (nt, np) = (200, 200);
        let mut s = 0.0;
        for i in 0..nt {
            let th = std::f64::consts::PI * (i as f64 + 0.5) / nt as f64;
            for j in 0..np {
                let ph = 2.0 * std::f64::consts::PI * j as f64 / np as f64;
                let o = atom::omega(kappa, tm, th, ph);
                s += (o[0].norm2() + o[1].norm2()) * th.sin();
            }
        }
        s *= std::f64::consts::PI / nt as f64 * 2.0 * std::f64::consts::PI / np as f64;
        assert!((s - 1.0).abs() < 1e-3, "κ={kappa}: {s}");
    }
    let _ = Hydrogen::mu();
}


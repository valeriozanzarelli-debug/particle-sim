//! Headless runs of the Bohr-scale hydrogen (src/atomo/idrogeno.rs): exams, development runs, exports.
//! Usage: RAYON_NUM_THREADS=2 cargo run --bin idrogeno -- <command> [args]
#[path = "../zigzag/cplx.rs"]
#[allow(dead_code)]
mod cplx;
mod idrogeno;

use cplx::C;
use idrogeno::*;
use std::time::Instant;
use rayon::prelude::*;

const ASSET_BASE: &str = "assets/atomo/agglomerato_zz.json";
/// the learned rule: IDRO_ASSET overrides (e.g. the parity/cubic-symmetric assets/atomo/agglomerato_zz_simm_v1.json)
fn asset() -> String { std::env::var("IDRO_ASSET").unwrap_or_else(|_| ASSET_BASE.to_string()) }

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("help");
    match cmd {
        "fft" => test_fft(),
        "debug" => dev_debug(),
        "massa" => { for th in [0.26765f64, 0.13382, 0.0334560545] { for kp in [0.02f64, 0.05, 0.1] {
            let w = |kc: f64| { let m = bloch_due_tick([kc, 0.0, 0.0], th); let tr = m[0][0] + m[1][1] + m[2][2] + m[3][3]; (tr.re / 4.0).acos() };
            let kc = kp * th; let e = (w(kc) - w(0.0)) / (2.0 * th);
            println!("theta {th}: k_phys {kp}: (E(k)-E(0)) / (k^2/2) = {:.6}  (relativistic sqrt(1+k^2)-1 over k^2/2: {:.6}; theta*cot theta {:.6})",
                     e / (kp * kp / 2.0), ((1.0 + kp * kp).sqrt() - 1.0) / (kp * kp / 2.0), th / th.tan()); } } }
        "testfine" => dev_test_fine(),
        "attraversa" => { let a: Vec<f64> = args[2..].iter().map(|s| s.parse().unwrap()).collect();
                          dev_attraversa(a[0] as usize, a[1] as usize, a[2], a[3], a[4], a[5] as usize); }
        "blocco" => dev_blocco(args[2].parse().unwrap(), args[3].parse().unwrap(), args[4].parse().unwrap()),
        "ibrido" => dev_ibrido(args[2].parse().unwrap(), args[3].parse().unwrap(), args[4].parse().unwrap(), args.get(5).map(|s| s.as_str()).unwrap_or("fine")),
        "livello" => {
            // livello <B_fine> <N> <state> : eigen-solve on one learned level
            let b: usize = args[2].parse().unwrap();
            let n: usize = args[3].parse().unwrap();
            let kind = args[4].as_str();
            let off: f64 = args.get(5).map(|s| s.parse().unwrap()).unwrap_or(0.0);
            dev_livello(b, n, kind, off);
        }
        "esame_ibrido" => esame_ibrido(args[2].parse().unwrap(), args[3].parse().unwrap(), args[4].parse().unwrap(), args[5].parse().unwrap(), args.get(6).cloned()),
        "esame_risonanza" => esame_risonanza(args[2].parse().unwrap(), args[3].parse().unwrap(), &args[4], args[5].parse().unwrap(),
                                             args[6].parse().unwrap(), args[7].parse().unwrap(), args.get(8).cloned()),
        "dispersione" => {
            // learned coarse G(k) against the fine walk (b ticks per coarse step, derived angle) and against relativity
            let reg = carica_regola(&asset(), "chain_A_3d", 256);
            let n = 64; let lv = Livello::con_nucleo(&reg, n, 0.0, 0.0);
            for b in [32usize, 16] {
                let a = lv.h / b as f64; let ang = angolo_fine(&lv, a);
                let w = |kc: [f64; 3]| { let m = bloch_due_tick([kc[0] / b as f64, kc[1] / b as f64, kc[2] / b as f64], ang);
                    let tr = m[0][0] + m[1][1] + m[2][2] + m[3][3]; (tr.re / 4.0).acos() };
                for (mm, dir) in [(2usize, 0usize), (4, 0), (8, 0), (16, 0), (4, 3), (8, 3)] {
                    let k = 2.0 * std::f64::consts::PI * mm as f64 / n as f64;
                    let (kv, s) = if dir == 0 { ([k, 0.0, 0.0], lv.idx(mm, 0, 0)) } else { ([k, k, k], lv.idx(mm, mm, mm)) };
                    let k2 = kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2];
                    let gw = b as f64 / 2.0 * (w(kv) - w([0.0; 3]));
                    let kp2 = k2 / (lv.h * lv.h);
                    let rel = lv.tau * ((1.0 + kp2).sqrt() - 1.0);
                    println!("b={b} k_cell={:.3} {}: G_rule/rel {:.6}  walk/rel {:.6}  rule/walk-1 {:+.3e}", k2.sqrt(), if dir == 0 { "asse" } else { "diag" },
                             lv.gk[s] / rel, gw / rel, lv.gk[s] / gw - 1.0);
                }
            }
        }
        "parita" => { let (b, n): (usize, usize) = (args[2].parse().unwrap(), args[3].parse().unwrap());
                      let reg = carica_regola(&asset(), "chain_A_3d", b); let lv = Livello::con_nucleo(&reg, n, 1.0, 0.0);
                      for st in ["1s", "2pz", "2px", "2py", "2s"] { let v = carica_stato(&format!("{}/{st}_B{b}_N{n}.c32", args[4]));
                          println!("RISULTATO parita {st}: 1 - |P| = {:.3e}", 1.0 - parita(&lv, &v).abs()); } }
        "spettro" => esame_spettro(args[2].parse().unwrap(), args[3].parse().unwrap(), args[4].parse().unwrap(), args.get(5).cloned()),
        _ => eprintln!("commands: fft | livello B N 1s|2s|2p"),
    }
}

fn test_fft() {
    for &n in &[8usize, 12, 16, 24, 48, 64] {
        let f = Fft::new(n);
        let mut a: Vec<C> = (0..n * n * n).map(|s| C::new(((s * 7919) % 101) as f64 / 101.0 - 0.5, ((s * 104729) % 37) as f64 / 37.0 - 0.5)).collect();
        let orig = a.clone();
        f.tre(&mut a, false);
        // naive check of a few coefficients
        let mut err: f64 = 0.0;
        for &(kx, ky, kz) in &[(1usize, 2usize, 3usize), (0, 0, 0), (n - 1, 1, 0), (n / 2, n / 3, 1)] {
            let mut s = C::ZERO;
            for x in 0..n { for y in 0..n { for z in 0..n {
                let ang = -2.0 * std::f64::consts::PI * ((kx * x + ky * y + kz * z) % n) as f64 / n as f64;
                s += orig[(x * n + y) * n + z] * C::expi(ang);
            }}}
            err = err.max((s - a[(kx * n + ky) * n + kz]).abs());
        }
        f.tre(&mut a, true);
        let back = a.iter().zip(orig.iter()).map(|(p, q)| (*p - *q).abs()).fold(0.0, f64::max);
        println!("n = {n}: forward vs naive {err:.2e}, round trip {back:.2e}");
    }
    for &n in &[128usize, 192, 256] {
        let f = Fft::new(n);
        let mut a = vec![C::new(1.0, 0.5); n * n * n];
        let t = Instant::now();
        f.tre(&mut a, false);
        f.tre(&mut a, true);
        println!("n = {n}: forward+inverse {:.3} s", t.elapsed().as_secs_f64());
    }
}

fn dev_livello(b: usize, n: usize, kind: &str, off: f64) {
    let t = Instant::now();
    let reg = carica_regola(&asset(), "chain_A_3d", b);
    let lv = Livello::con_nucleo(&reg, n, 1.0, off);
    eprintln!("level B = {b}: h = {:.4} (a0/{:.2}), phi0 = {:.6}, box ±{:.2} a0, max|Im g| = {:.2e}, spin = {:.1e} ({:.1} s)",
              lv.h, a0() / lv.h, lv.phi0, n as f64 / 2.0 * lv.h / a0(), lv.im_g, lv.spin, t.elapsed().as_secs_f64());
    let (e1, e2s, e2p) = riferimenti_ev();
    let s1 = lobpcg(&lv, lv.guess("1s"), &[], 0.5 * ALPHA * ALPHA, 1e-9, 200, true);
    println!("B={b} N={n} off={off} 1s: E_rule = {:.6} eV (gen {:.6}) vs Dirac×mu {:.6}: {:+.4} % | resid_rule {:.1e} it {} ({:.0} s)",
             s1.e_regola * ev(), s1.e_gen * ev(), e1, 100.0 * (s1.e_regola * ev() / e1 - 1.0), s1.resid_regola, s1.iter, t.elapsed().as_secs_f64());
    if kind == "1s" { return; }
    let s2p = lobpcg(&lv, lv.guess("2p"), &[&s1.v], ALPHA * ALPHA / 8.0, 1e-9, 400, true);
    println!("B={b} N={n} off={off} 2p: E_rule = {:.6} eV vs Dirac×mu centroid {:.6}: {:+.4} % | resid_rule {:.1e} it {} ({:.0} s)",
             s2p.e_regola * ev(), e2p, 100.0 * (s2p.e_regola * ev() / e2p - 1.0), s2p.resid_regola, s2p.iter, t.elapsed().as_secs_f64());
    if kind == "2p" { return; }
    let px = permuta(&s2p.v, n, 2);
    let py = permuta(&s2p.v, n, 1);
    let s2 = lobpcg(&lv, lv.guess("2s"), &[&s1.v, &s2p.v, &px, &py], ALPHA * ALPHA / 8.0, 1e-9, 400, true);
    println!("B={b} N={n} off={off} 2s: E_rule = {:.6} eV vs Dirac×mu {:.6}: {:+.4} % | resid_rule {:.1e} it {} ({:.0} s)",
             s2.e_regola * ev(), e2s, 100.0 * (s2.e_regola * ev() / e2s - 1.0), s2.resid_regola, s2.iter, t.elapsed().as_secs_f64());
}

/// the state with axes permuted: ax = 2 swaps x↔z, ax = 1 swaps y↔z (2p_z → 2p_x, 2p_y)
pub fn permuta(v: &[C], n: usize, ax: usize) -> Vec<C> {
    let mut o = vec![C::ZERO; v.len()];
    for i in 0..n { for j in 0..n { for k in 0..n {
        let src = if ax == 2 { (k * n + j) * n + i } else { (i * n + k) * n + j };
        o[(i * n + j) * n + k] = v[src];
    }}}
    o
}

#[allow(dead_code)]
pub fn dev_debug() {
    let reg = carica_regola(&asset(), "chain_A_3d", 256);
    let n = 64;
    let lv = Livello::con_nucleo(&reg, n, 1.0, 0.5);
    for m in [(1usize, 0usize, 0usize), (0, 1, 0), (2, 0, 0), (1, 1, 1)] {
        let s = (m.0 * n + m.1) * n + m.2;
        let k = 2.0 * std::f64::consts::PI / n as f64;
        let k2 = k * k * ((m.0 * m.0 + m.1 * m.1 + m.2 * m.2) as f64);
        println!("G at m={:?}: {:.6e}  k^2/(2h) = {:.6e}", m, lv.gk[s], k2 / (2.0 * lv.h));
    }
    let tp = 2.0 * std::f64::consts::PI;
    let mis = (-lv.tau - lv.phi0 + lv.gk[0]).rem_euclid(tp);
    println!("G(0) = {:.6e}, tau = {:.10}, phi0 = {:.10}, theta_f = {:.10}, b*theta_f = {:.10}; fine-vs-coarse phase per step at k=0: {:.6e} (wrapped to (-pi,pi]: {:.6e})",
             lv.gk[0], lv.tau, lv.phi0, reg.theta_fine, reg.theta_fine * reg.b_fine as f64, mis, if mis > std::f64::consts::PI { mis - tp } else { mis });
    let c = lv.c0;
    for d in 0..4 { println!("V[c0-1+{d}] = {:.6e}  -alpha/|x| = {:.6e}", lv.v[lv.idx(c + d, c + d, c + d)], -ALPHA / (3f64.sqrt() * lv.x(c + d)).abs()); }
}

/// checks of the fine block: G_j, R∘P = I, walk vs Bloch matrix on a plane wave, reconstruction accuracy
pub fn dev_test_fine() {
    let reg = carica_regola(&asset(), "chain_A_3d", 256);
    let lv = Livello::new(&reg, 16, 0.0);
    let th = lv.h / 32.0;
    let g = g_positiva(th);
    for j in 0..3 { println!("G_{j} = [[{:.5},{:.5}],[{:.5},{:.5}]] (re/im of 00: {:.5}/{:.5})", g[j][0][0].abs(), g[j][0][1].abs(), g[j][1][0].abs(), g[j][1][1].abs(), g[j][0][0].re, g[j][0][0].im); }
    // free fine block, prolong a smooth coarse field and restrict back
    let mut f = Fine::new(&lv, 1, 1, 32, 0.0);
    let n = lv.n;
    let coarse: Vec<C> = (0..n * n * n).map(|s| {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        let (x, y, z) = (lv.x(i), lv.x(j), lv.x(k));
        C::expi(0.02 * x - 0.01 * z) * (-(x * x + y * y + z * z) / (2.0 * 30.0f64.powi(2))).exp()
    }).collect();
    f.prolunga(&coarse, true);
    let back = f.restringi(false);
    for (s, v) in &back { println!("R(P(c)) at {s}: {:.12e} vs {:.12e}", v.abs(), coarse[*s].abs()); }
    // walk on a plane wave inside the patch vs the Bloch matrix
    let p = f.p;
    let k = [0.3, -0.2, 0.1];
    let s0 = [C::new(0.3, 0.1), C::new(-0.2, 0.5), C::new(0.7, 0.0), C::new(0.1, -0.4)];
    for c in 0..4 { for i in 0..p { for j in 0..p { for l in 0..p {
        f.psi[c][(i * p + j) * p + l] = s0[c] * C::expi(k[0] * i as f64 + k[1] * j as f64 + k[2] * l as f64);
    }}}}
    f.tick = 0;
    f.tick(1);
    f.tick(2);
    // the block walks with the derived angle and the uniform rest phase e^{-i(a_f - angle)} per tick
    let w = bloch_due_tick(k, f.ang);
    let (i, j, l) = (p / 2, p / 2 + 3, p / 2 - 5);
    let s = (i * p + j) * p + l;
    let ph = C::expi(k[0] * i as f64 + k[1] * j as f64 + k[2] * l as f64) * C::expi(-2.0 * (f.th - f.ang));
    println!("block: a_f = {:.6}, derived angle {:.6}", f.th, f.ang);
    let mut err: f64 = 0.0;
    for a in 0..4 { let mut v = C::ZERO; for bq in 0..4 { v += w[a][bq] * s0[bq]; } err = err.max((v * ph - f.psi[a][s]).abs()); }
    println!("walk vs Bloch (2 ticks): {err:.2e}");
    // positive band check: a positive-energy plane wave from the relation, its energy
    let t = relazione_positiva([0.01, 0.0, 0.0], th);
    println!("T(k=0.01 x) = {:?}", t);
}

/// phase-slope energy from overlaps C(n) (n = 0..), fit on n >= from: E = −(slope − φ₀)/τ
pub fn energia_da_fase(cs: &[C], from: usize, phi0: f64, tau: f64) -> f64 {
    let mut ph = vec![0.0; cs.len()];
    for n in 1..cs.len() {
        let mut d = (cs[n] * cs[n - 1].conj()).arg() - phi0;
        while d > std::f64::consts::PI { d -= 2.0 * std::f64::consts::PI; }
        while d < -std::f64::consts::PI { d += 2.0 * std::f64::consts::PI; }
        ph[n] = ph[n - 1] + d;
    }
    let xs: Vec<f64> = (from..cs.len()).map(|n| n as f64).collect();
    let ys: Vec<f64> = (from..cs.len()).map(|n| ph[n]).collect();
    let m = xs.len() as f64;
    let (sx, sy) = (xs.iter().sum::<f64>() / m, ys.iter().sum::<f64>() / m);
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..xs.len() { num += (xs[i] - sx) * (ys[i] - sy); den += (xs[i] - sx).powi(2); }
    -(num / den) / tau
}

pub fn dev_ibrido(n: usize, ni: usize, passi: usize, calcio: &str) {
    let t0 = Instant::now();
    let reg = carica_regola(&asset(), "chain_A_3d", 256);
    let off = if ni % 2 == 1 || ni == 0 { 0.0 } else { 0.5 };
    let lv = Livello::con_nucleo(&reg, n, 1.0, off);
    let s1 = lobpcg(&lv, lv.guess("1s"), &[], 0.5 * ALPHA * ALPHA, 1e-10, 200, false);
    let (e1, _, _) = riferimenti_ev();
    println!("coarse 1s (N={n}, off={off}): {:.6} eV ({:+.4} %) [{:.0} s]", s1.e_regola * ev(), 100.0 * (s1.e_regola * ev() / e1 - 1.0), t0.elapsed().as_secs_f64());
    let v0 = s1.v.clone();
    let lv2 = Livello::con_nucleo(&reg, n, 1.0, off);
    let fine = Fine::new(&lv2, ni.max(1), 1, 32, 1.0);
    println!("fine block: patch {}^3, a_f = {:.5}, angle {:.6} (kappa_rule {:.6}), G_x = {:?}", fine.p, fine.th, fine.ang, kappa_regola(&lv), fine.gj[0][0][1]);
    let cal = if calcio == "grosso" { Calcio::Grosso } else { Calcio::DalFine };
    let mut at = Atomo::new(lv2, v0.clone(), if ni > 0 { Some(fine) } else { None }, cal);
    at.covarianza = std::env::var("IDRO_COV").is_ok();
    let mut cs = vec![dot(&v0, &at.psi)];
    let mut cc = vec![C::ONE];
    let mut co = v0.clone();
    for k in 1..=passi {
        let t = Instant::now();
        at.passo();
        let dt = t.elapsed().as_secs_f64();
        lv.passo(&mut co);
        cs.push(dot(&v0, &at.psi));
        cc.push(dot(&v0, &co));
        if k % 4 == 0 || k < 4 {
            let (ni_, bt, gs) = match at.fine.as_ref() { Some(f) => { let (a, b) = f.interno_norma_beta(); (a, b, f.giu_su) } None => (0.0, 0.0, 0.0) };
            println!("step {k}: |C| {:.9} norm {:.9} arg(C_h/C_c) {:+.4e} interior p {:.4e} beta {:.4e} down/up {:.1e} ({:.2} s/step)",
                     cs[k].abs(), norm2(&at.psi), (cs[k] * cc[k].conj()).arg(), ni_, bt, gs, dt);
        }
    }
    for from in [4usize, 8, 16] {
        if from + 4 > passi { continue; }
        let eh = energia_da_fase(&cs, from, at.lv.phi0, at.lv.tau);
        let ec = energia_da_fase(&cc, from, lv.phi0, lv.tau);
        println!("from step {from}: hybrid {:.6} eV ({:+.4} %), coarse {:.6} eV ({:+.4} %), hybrid/coarse - 1 = {:+.4e}",
                 eh * ev(), 100.0 * (eh * ev() / e1 - 1.0), ec * ev(), 100.0 * (ec * ev() / e1 - 1.0), eh / ec - 1.0);
    }
    // near-nucleus shape: hybrid and coarse values vs the exact 1s block means (continuum Schrödinger, same norm far out)
    let c0 = at.lv.c0;
    let exact = |d: [i64; 3]| -> f64 {
        let ns = 16;
        let mut s = 0.0;
        for a in 0..ns { for b in 0..ns { for c in 0..ns {
            let x = [(d[0] as f64 + off + (a as f64 + 0.5) / ns as f64 - 0.5) * at.lv.h, (d[1] as f64 + off + (b as f64 + 0.5) / ns as f64 - 0.5) * at.lv.h,
                     (d[2] as f64 + off + (c as f64 + 0.5) / ns as f64 - 0.5) * at.lv.h];
            s += (-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt() / a0()).exp();
        }}}
        s / (ns * ns * ns) as f64
    };
    let ph = C::expi(-cs[passi].arg());
    let refd = [4i64, 0, 0];
    let id = |d: [i64; 3]| at.lv.idx((c0 as i64 + d[0]) as usize, (c0 as i64 + d[1]) as usize, (c0 as i64 + d[2]) as usize);
    let (eh0, ec0, ex0) = ((at.psi[id(refd)] * ph).re, v0[id(refd)].re, exact(refd));
    for d in [[0i64, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1], [2, 0, 0], [3, 0, 0]] {
        let e = exact(d) / ex0;
        println!("cell {:?}: exact {:.6}  coarse/exact {:.6}  hybrid/exact {:.6}", d, e, v0[id(d)].re / ec0 / e, (at.psi[id(d)] * ph).re / eh0 / e);
    }
    println!("total {:.0} s", t0.elapsed().as_secs_f64());
}

/// PASS-THROUGH (no nucleus, z = 0): a packet crosses the fine block; hybrid vs the coarse rule alone.
/// returns (fidelity |<c|h>|^2, reflected fraction (k_x < 0 part of h − c), norm change)
pub fn dev_attraversa(n: usize, ni: usize, kcell: f64, sig: f64, x0: f64, passi: usize) -> (f64, f64, f64) {
    let t0 = Instant::now();
    let reg = carica_regola(&asset(), "chain_A_3d", 256);
    let off = if ni % 2 == 1 { 0.0 } else { 0.5 };
    let lv = Livello::con_nucleo(&reg, n, 0.0, off);
    let mut c: Vec<C> = (0..n * n * n).map(|s| {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        let (x, y, z) = (lv.x(i) / lv.h - x0, lv.x(j) / lv.h, lv.x(k) / lv.h);
        C::expi(kcell * x) * (-(x * x + y * y + z * z) / (4.0 * sig * sig)).exp()
    }).collect();
    let k = 1.0 / norm2(&c).sqrt();
    scal(&mut c, C::real(k));
    let lv2 = Livello::con_nucleo(&reg, n, 0.0, off);
    let fine = Fine::new(&lv2, ni, 1, 32, 0.0);
    let mut at = Atomo::new(lv2, c.clone(), Some(fine), Calcio::DalFine);
    for s in 1..=passi {
        at.passo();
        lv.passo(&mut c);
        if s % 50 == 0 || passi <= 32 { let o = dot(&c, &at.psi); eprintln!("  step {s}: norm h {:.10}, |<c|h>|^2 {:.10}, arg<c|h> {:+.3e} ({:.0} s)", norm2(&at.psi), o.norm2(), o.arg(), t0.elapsed().as_secs_f64()); }
    }
    let fid = dot(&c, &at.psi).norm2();
    let mut d = at.psi.clone();
    axpy(&mut d, C::real(-1.0), &c);
    lv.fft.tre(&mut d, false);
    let mut back = 0.0;
    let mut tot = 0.0;
    for i in 0..n { for j in 0..n { for l in 0..n {
        let v = d[(i * n + j) * n + l].norm2();
        tot += v;
        if i > n / 2 { back += v; }
    }}}
    let nn = (n * n * n) as f64;
    let refl = back / nn;
    let dn = norm2(&at.psi) - 1.0;
    println!("attraversa N={n} ni={ni} k={kcell} sigma={sig}: fidelity {fid:.8}, reflected {refl:.3e}, |h-c|^2 {:.3e}, norm change {dn:+.3e} ({:.0} s)",
             tot / nn, t0.elapsed().as_secs_f64());
    (fid, refl, dn)
}

/// time blocking: m steps at once, e^{−iV mτ/2} e^{−i m G} e^{−iV mτ/2}; energy of the 1s (and 2p) read from it
pub fn dev_blocco(b: usize, n: usize, off: f64) {
    let reg = carica_regola(&asset(), "chain_A_3d", b);
    let lv = Livello::con_nucleo(&reg, n, 1.0, off);
    let s1 = lobpcg(&lv, lv.guess("1s"), &[], 0.5 * ALPHA * ALPHA, 1e-10, 200, false);
    let s2 = lobpcg(&lv, lv.guess("2p"), &[&s1.v], ALPHA * ALPHA / 8.0, 1e-10, 400, false);
    for m in [1usize, 2, 4, 8, 16, 32, 64, 128] {
        let mm = lv.molt_n(m as f64);
        let k: Vec<C> = lv.v.iter().map(|v| C::expi(-0.5 * v * lv.tau * m as f64)).collect();
        let mut es = vec![];
        for v in [&s1.v, &s2.v] {
            let mut f = v.clone();
            for (a, kk) in f.iter_mut().zip(k.iter()) { *a = *a * *kk; }
            lv.fft.tre(&mut f, false);
            for (a, q) in f.iter_mut().zip(mm.iter()) { *a = *a * *q; }
            lv.fft.tre(&mut f, true);
            for (a, kk) in f.iter_mut().zip(k.iter()) { *a = *a * *kk; }
            let lam = dot(v, &f) * C::expi(-(m as f64) * lv.phi0);
            es.push(-lam.arg() / (m as f64 * lv.tau));
        }
        println!("blocco B={b} N={n} m={m} (dt = {:.4} a.u.): 1s {:.7} eV, 2p {:.7} eV, rel. to m=1: {:+.3e} {:+.3e}",
                 m as f64 * lv.tau * ALPHA * ALPHA, es[0] * ev(), es[1] * ev(), es[0] / s1.e_regola - 1.0, es[1] / s2.e_regola - 1.0);
    }
}

/// spettro B N off [dir]: 1s, 2p_z, 2p_x, 2p_y, 2s of one learned level (each a separate LOBPCG with deflation of the
/// lower states already found), with parity and the 2p triplet splitting; optionally saves the states (f32) in dir
fn esame_spettro(b: usize, n: usize, off: f64, dir: Option<String>) {
    let t = Instant::now();
    let reg = carica_regola(&asset(), "chain_A_3d", b);
    let lv = Livello::con_nucleo(&reg, n, 1.0, off);
    println!("level B = {b}: h = {:.4} (a0/{:.2}), tau = {:.4}, phi0 = {:.6}, box ±{:.2} a0, N = {n}, off = {off}, max|Im g| = {:.2e}",
             lv.h, a0() / lv.h, lv.tau, lv.phi0, n as f64 / 2.0 * lv.h / a0(), lv.im_g);
    let (e1, e2s, e2p) = riferimenti_ev();
    let riga = |nome: &str, s: &Autostato, rif: f64, lv: &Livello| {
        println!("{nome}: E_rule = {:.6} eV (gen {:.6}) vs Dirac×mu {:.6}: {:+.4} % | parity {:+.6} | resid_rule {:.1e} resid {:.1e} it {} ({:.0} s)",
                 s.e_regola * ev(), s.e_gen * ev(), rif, 100.0 * (s.e_regola * ev() / rif - 1.0), parita(lv, &s.v),
                 s.resid_regola, s.resid, s.iter, t.elapsed().as_secs_f64());
    };
    let tol = 3e-9;
    let s1 = lobpcg(&lv, lv.guess("1s"), &[], 0.5 * ALPHA * ALPHA, tol, 300, true);
    riga("1s ", &s1, e1, &lv);
    let pz = lobpcg(&lv, lv.guess("2p"), &[&s1.v], ALPHA * ALPHA / 8.0, tol, 300, true);
    riga("2pz", &pz, e2p, &lv);
    let px = lobpcg(&lv, permuta(&pz.v, n, 2), &[&s1.v, &pz.v], ALPHA * ALPHA / 8.0, tol, 300, true);
    riga("2px", &px, e2p, &lv);
    let py = lobpcg(&lv, permuta(&pz.v, n, 1), &[&s1.v, &pz.v, &px.v], ALPHA * ALPHA / 8.0, tol, 300, true);
    riga("2py", &py, e2p, &lv);
    let s2 = lobpcg(&lv, lv.guess("2s"), &[&s1.v, &pz.v, &px.v, &py.v], ALPHA * ALPHA / 8.0, tol, 300, true);
    riga("2s ", &s2, e2s, &lv);
    let ep = [pz.e_regola, px.e_regola, py.e_regola];
    let em = (ep[0] + ep[1] + ep[2]) / 3.0;
    let spl = ep.iter().map(|e| (e - em).abs()).fold(0.0, f64::max);
    println!("2p mean {:.6} eV ({:+.4} %), triplet spread {:.2e} relative; 2s-2p {:+.6} eV; 1s->2p {:.6} eV (Dirac×mu {:.6}: {:+.4} %)",
             em * ev(), 100.0 * (em * ev() / e2p - 1.0), spl / em.abs(), (s2.e_regola - em) * ev(),
             (em - s1.e_regola) * ev(), e2p - e1, 100.0 * ((em - s1.e_regola) * ev() / (e2p - e1) - 1.0));
    if let Some(d) = dir {
        std::fs::create_dir_all(&d).unwrap();
        for (nome, s) in [("1s", &s1), ("2pz", &pz), ("2px", &px), ("2py", &py), ("2s", &s2)] {
            salva_stato(&format!("{d}/{nome}_B{b}_N{n}.c32"), &s.v);
        }
    }
    println!("total {:.0} s", t.elapsed().as_secs_f64());
}

fn scrivi(path: &str, data: &[f64], scala: f64, append: bool) {
    use std::io::Write;
    let f = std::fs::OpenOptions::new().create(true).append(append).write(true).truncate(!append).open(path).unwrap();
    let mut w = std::io::BufWriter::new(f);
    scrivi_f16(&mut w, data, scala);
    w.flush().unwrap();
}
/// central cube of side m of a coarse field (nucleus cell c0 at the centre)
fn ritaglio(lv: &Livello, f: &[f64], m: usize) -> Vec<f64> {
    let s0 = lv.c0 - m / 2;
    let mut o = Vec::with_capacity(m * m * m);
    for i in 0..m { for j in 0..m { for k in 0..m { o.push(f[lv.idx(s0 + i, s0 + j, s0 + k)]); } } }
    o
}
fn dens(f: &[C]) -> Vec<f64> { f.par_iter().map(|c| c.norm2()).collect() }
fn h_media(lv: &Livello, f: &[C]) -> f64 { let mut o = vec![C::ZERO; f.len()]; lv.ham(f, &mut o); dot(f, &o).re / norm2(f) }

/// E1/E2(1s)/E3/E5/E6/E8 — the hybrid atom: coarse a0/16 learned rule + fine zig-zag block on the nucleus.
/// esame_ibrido B N ni passi [export dir]   (B = 256: h = a0/16; B = 512: h = a0/8)
fn esame_ibrido(bl: usize, n: usize, ni: usize, passi: usize, esporta: Option<String>) {
    let t0 = Instant::now();
    let reg = carica_regola(&asset(), "chain_A_3d", bl);
    let off = if ni % 2 == 1 || ni == 0 { 0.0 } else { 0.5 };
    let lv = Livello::con_nucleo(&reg, n, 1.0, off);
    let s1 = lobpcg(&lv, lv.guess("1s"), &[], 0.5 * ALPHA * ALPHA, 3e-9, 300, false);
    let (e1, _, _) = riferimenti_ev();
    println!("coarse 1s (B={bl} N={n} box ±{:.1} a0 off={off}): E_rule {:.6} eV ({:+.4} %) resid_rule {:.1e} [{:.0} s]", n as f64 / 2.0 * lv.h / a0(),
             s1.e_regola * ev(), 100.0 * (s1.e_regola * ev() / e1 - 1.0), s1.resid_regola, t0.elapsed().as_secs_f64());
    let v0 = s1.v.clone();
    let lv2 = Livello::con_nucleo(&reg, n, 1.0, off);
    let bf: usize = std::env::var("IDRO_B").ok().and_then(|v| v.parse().ok()).unwrap_or(32);
    let fine = if ni > 0 { Some(Fine::new(&lv2, ni, 1, bf, 1.0)) } else { None };
    if let Some(f) = fine.as_ref() {
        println!("fine block: ni = {ni}, patch {}^3, a_f = {:.5} (a0/{:.0}), walk angle {:.6} (kappa_rule {:.6})", f.p, f.th, a0() / f.th, f.ang, kappa_regola(&lv));
    }
    let mut at = Atomo::new(lv2, v0.clone(), fine, Calcio::DalFine);
    at.covarianza = std::env::var("IDRO_COV").is_ok();
    let q0 = norm2(&at.psi);
    let h0 = h_media(&lv, &at.psi);
    let mut cs = vec![dot(&v0, &at.psi)];
    let mut cc = vec![C::ONE];
    let mut co = v0.clone();
    let meta = passi / 2;
    let (mut rho_meta, mut fine_meta, mut q_meta, mut h_meta) = (vec![], vec![], 0.0, 0.0);
    let interno = |f: &Fine| -> Vec<f64> {
        let (lo, hi) = (f.nh * f.b, (f.nh + f.ni) * f.b);
        let mut o = vec![];
        for i in lo..hi { for j in lo..hi { for k in lo..hi { let q = f.idx(i, j, k); o.push((0..4).map(|c| f.psi[c][q].norm2()).sum::<f64>()); } } }
        o
    };
    let mut tempi = vec![];
    let vk = lv.velocita_k();
    let m_g = 64usize.min(n);
    let a0c = (a0() / lv.h).powi(3);
    for k in 1..=passi {
        let t = Instant::now();
        let ultimo = k == passi && esporta.is_some() && at.fine.is_some();
        if ultimo {
            let d = esporta.clone().unwrap();
            let (pf, ff) = (format!("{d}/rho_fine.f16"), format!("{d}/chir_fine.f16"));
            let _ = std::fs::remove_file(&pf); let _ = std::fs::remove_file(&ff);
            let mut oss = |f: &Fine, _t: usize| {
                let (lo, hi) = (f.nh * f.b, (f.nh + f.ni) * f.b);
                let a3 = (a0() / f.th).powi(3);
                let mut r = vec![]; let mut c = vec![];
                for i in lo..hi { for j in lo..hi { for kk in lo..hi { let q = f.idx(i, j, kk);
                    let (l, rr) = (f.psi[0][q].norm2() + f.psi[1][q].norm2(), f.psi[2][q].norm2() + f.psi[3][q].norm2());
                    r.push((l + rr) * a3); c.push((rr - l) * a3); } } }
                scrivi(&pf, &r, 1.0, true); scrivi(&ff, &c, 1.0, true);
            };
            at.passo_osserva(&mut oss);
        } else { at.passo(); }
        tempi.push(t.elapsed().as_secs_f64());
        lv.passo(&mut co);
        cs.push(dot(&v0, &at.psi));
        cc.push(dot(&v0, &co));
        if k == meta {
            rho_meta = dens(&at.psi); q_meta = norm2(&at.psi); h_meta = h_media(&lv, &at.psi);
            if let Some(f) = at.fine.as_ref() { fine_meta = interno(f); }
        }
        if k % 8 == 0 { eprintln!("  step {k}/{passi}: norm {:.10} ({:.2} s/step)", norm2(&at.psi), tempi[k - 1]); }
    }
    let q1 = norm2(&at.psi);
    let h1 = h_media(&lv, &at.psi);
    let from = 8.min(passi / 2);
    if passi >= 48 {
        for f0 in [8usize, passi / 3, passi / 2, 2 * passi / 3] {
            let e = energia_da_fase(&cs, f0, at.lv.phi0, at.lv.tau);
            println!("RISULTATO finestra [{f0},{passi}]: E_hybrid {:.6} eV ({:+.4} %)", e * ev(), 100.0 * (e * ev() / e1 - 1.0));
        }
    }
    let eh = energia_da_fase(&cs, from, at.lv.phi0, at.lv.tau);
    let ec = energia_da_fase(&cc, from, lv.phi0, lv.tau);
    let rho1 = dens(&at.psi);
    let rmax = rho_meta.iter().cloned().fold(0.0, f64::max);
    let staz_g = rho1.iter().zip(rho_meta.iter()).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max) / rmax;
    let staz_f = match at.fine.as_ref() { Some(f) => { let r = interno(f); let m = fine_meta.iter().cloned().fold(0.0, f64::max);
        r.iter().zip(fine_meta.iter()).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max) / m } None => 0.0 };
    let tm = tempi.iter().skip(1).sum::<f64>() / (tempi.len() - 1).max(1) as f64;
    let ticks = at.fine.as_ref().map(|f| f.b).unwrap_or(1);
    println!("RISULTATO ibrido N={n} ni={ni} passi={passi}: E_hybrid {:.6} eV ({:+.4} %), E_coarse_only {:.6} eV ({:+.4} %), hybrid/coarse-1 {:+.4e}",
             eh * ev(), 100.0 * (eh * ev() / e1 - 1.0), ec * ev(), 100.0 * (ec * ev() / e1 - 1.0), eh / ec - 1.0);
    println!("RISULTATO conservazione: charge drift per step {:+.3e} (whole run), {:+.3e} (2nd half); <H_coarse> drift {:+.3e} (whole), {:+.3e} (2nd half)",
             (q1 - q0) / passi as f64, (q1 - q_meta) / (passi - meta) as f64, (h1 - h0) / h0.abs(), (h1 - h_meta) / h0.abs());
    println!("reflux: {} (last factor-1 {:+.3e}), asset {}", at.reflux, at.ultimo_reflux, asset());
    println!("RISULTATO stazionarieta: coarse density max|drho|/max rho {:.3e}, fine block {:.3e} (step {meta} -> {passi})", staz_g, staz_f);
    println!("RISULTATO costo: {:.3} s per coarse step ({} fine ticks inside), 2 threads; coarse-only step {:.3} s",
             tm, ticks, { let t = Instant::now(); lv.passo(&mut co); t.elapsed().as_secs_f64() });
    if let Some(d) = esporta.as_ref() {
        let r = ritaglio(&lv, &rho1, m_g).iter().map(|x| x * a0c).collect::<Vec<_>>();
        scrivi(&format!("{d}/rho_grosso.f16"), &r, 1.0, false);
        let j = lv.corrente(&at.psi, &vk);
        let mut jj = vec![];
        for a in 0..3 { jj.extend(ritaglio(&lv, &j[a], m_g).iter().map(|x| x * a0c / ALPHA)); }
        scrivi(&format!("{d}/j_grosso.f16"), &jj, 1.0, false);
        println!("RISULTATO corrente spuria (1s vero: 0): max|j| = {:.3e} alpha c/a0^3 (max rho {:.4} a0^-3)",
                 jj.iter().fold(0.0f64, |m, x| m.max(x.abs())), r.iter().cloned().fold(0.0, f64::max));
        println!("export: coarse crop {m_g}^3 (h = a0/{:.0}), fine interior {}^3 x {} ticks (a_f = a0/{:.0})", a0() / lv.h,
                 at.fine.as_ref().map(|f| f.ni * f.b).unwrap_or(0), ticks, at.fine.as_ref().map(|f| a0() / f.th).unwrap_or(0.0));
    }
    println!("total {:.0} s", t0.elapsed().as_secs_f64());
}

/// E4 — the 1s→2p resonance from the time evolution of (1s + 2p_z)/√2 with the learned rule (level B, box N):
/// `passi` true steps of the rule (m = 1): ω from the phase slopes of the projections, dipole ⟨z⟩(t).
/// Export (far view): `fot` frames over one Lyman-α period built from the two rule eigenstates with their rule
/// energies (the "postcard" evolution — exact for the rule up to the eigen residuals; checked against the true steps).
fn esame_risonanza(b: usize, n: usize, dir: &str, passi: usize, _m: usize, fot: usize, esporta: Option<String>) {
    let t0 = Instant::now();
    let reg = carica_regola(&asset(), "chain_A_3d", b);
    let lv = Livello::con_nucleo(&reg, n, 1.0, 0.0);
    let s1 = carica_stato(&format!("{dir}/1s_B{b}_N{n}.c32"));
    let p2 = carica_stato(&format!("{dir}/2pz_B{b}_N{n}.c32"));
    let (e1r, r1) = lv.energia_regola(&s1);
    let (e2r, r2) = lv.energia_regola(&p2);
    let (e1, _, e2p) = riferimenti_ev();
    let dref = e2p - e1;
    println!("level B={b} N={n} (h = a0/{:.0}, box ±{:.1} a0): E1s {:.6} eV (resid {:.1e}), E2p {:.6} eV (resid {:.1e}); eigen difference {:.6} eV vs Dirac×mu {:.6} ({:+.4} %); <1s|2p> = {:.1e}",
             a0() / lv.h, n as f64 / 2.0 * lv.h / a0(), e1r * ev(), r1, e2r * ev(), r2, (e2r - e1r) * ev(), dref, 100.0 * ((e2r - e1r) * ev() / dref - 1.0), dot(&s1, &p2).abs());
    let mut psi = s1.clone();
    axpy(&mut psi, C::ONE, &p2);
    let k = 1.0 / norm2(&psi).sqrt();
    scal(&mut psi, C::real(k));
    let zc: Vec<f64> = (0..n * n * n).map(|s| lv.x(s % n) / a0()).collect();
    let dip = |f: &[C]| -> f64 { f.par_iter().zip(zc.par_iter()).map(|(c, z)| c.norm2() * z).sum::<f64>() / norm2(f) };
    // the postcard: ψ(t) = c1 e^{−iE1 t} 1s + c2 e^{−iE2 t} 2p (rest phase removed, as in the energies)
    let (c10, c20) = (dot(&s1, &psi), dot(&p2, &psi));
    let cartolina = |t: f64| -> Vec<C> {
        let (a, bq) = (c10 * C::expi(-e1r * t), c20 * C::expi(-e2r * t));
        s1.par_iter().zip(p2.par_iter()).map(|(x, y)| a * *x + bq * *y).collect()
    };
    // for the far view: the three 2p states found by the solver are mixed (the triplet is degenerate to 3e-5), so the
    // p state along z is rebuilt as Σ_i d_i* p_i / |d| with d_i = ⟨1s|z|p_i⟩, each p_i keeping its own rule energy
    let px = carica_stato(&format!("{dir}/2px_B{b}_N{n}.c32"));
    let py = carica_stato(&format!("{dir}/2py_B{b}_N{n}.c32"));
    let ps = [&p2, &px, &py];
    let es: Vec<f64> = ps.iter().map(|q| lv.energia_regola(q).0).collect();
    let dz: Vec<C> = ps.iter().map(|q| s1.par_iter().zip(q.par_iter()).zip(zc.par_iter()).map(|((a, c), z)| a.conj() * *c * *z).reduce(|| C::ZERO, |x, y| x + y)).collect();
    let dn = dz.iter().map(|d| d.norm2()).sum::<f64>().sqrt();
    println!("2p triplet: E = {:.6} {:.6} {:.6} eV, <1s|z|p_i> = {:.4} {:.4} {:.4} a0 -> |d| = {:.5} a0 (exact 128√2/243 = 0.74494)",
             es[0] * ev(), es[1] * ev(), es[2] * ev(), dz[0].abs(), dz[1].abs(), dz[2].abs(), dn);
    let cz: Vec<C> = dz.iter().map(|d| d.conj().scale(1.0 / dn)).collect();
    let r2 = 0.5f64.sqrt();
    let cartolina_z = |t: f64| -> Vec<C> {
        let a = C::expi(-e1r * t).scale(r2);
        let w: Vec<C> = (0..3).map(|i| cz[i] * C::expi(-es[i] * t).scale(r2)).collect();
        (0..s1.len()).into_par_iter().map(|q| a * s1[q] + w[0] * p2[q] + w[1] * px[q] + w[2] * py[q]).collect()
    };
    if passi > 0 {
    let mut f = psi.clone();
    let (mut c1, mut c2) = (vec![dot(&s1, &f)], vec![dot(&p2, &f)]);
    let mut ds = vec![dip(&f)];
    for k in 1..=passi {
        lv.passo(&mut f);
        c1.push(dot(&s1, &f)); c2.push(dot(&p2, &f));
        if k % 10 == 0 { ds.push(dip(&f)); }
        if k % 25 == 0 { eprintln!("  step {k}/{passi} ({:.0} s)", t0.elapsed().as_secs_f64()); }
    }
    let ea = energia_da_fase(&c1, 0, lv.phi0, lv.tau);
    let eb = energia_da_fase(&c2, 0, lv.phi0, lv.tau);
    let tfin = passi as f64 * lv.tau;
    // compare the true-step state with the postcard at the same time (rest phase e^{i φ0 n} removed)
    let mut g = cartolina(tfin);
    scal(&mut g, C::expi(lv.phi0 * passi as f64));
    let fid = dot(&g, &f).norm2() / (norm2(&g) * norm2(&f));
    let dpred = dip(&g);
    println!("RISULTATO E4 (m=1, {passi} steps = {:.4} a.u.): E1s(t) {:.6} eV, E2p(t) {:.6} eV, omega {:.6} eV vs Dirac×mu {:.6}: {:+.4} %; vs eigen difference {:+.2e}; |c1|^2+|c2|^2 = {:.9}; dipole {:.6} -> {:.6} a0 (postcard {:.6}); fidelity true vs postcard {:.10}",
             tfin * ALPHA * ALPHA, ea * ev(), eb * ev(), (eb - ea) * ev(), dref, 100.0 * ((eb - ea) * ev() / dref - 1.0),
             (eb - ea) / (e2r - e1r) - 1.0, c1[passi].norm2() + c2[passi].norm2(), ds[0], dip(&f), dpred, fid);
    let lam_nm = 1239.84198 / ((e2r - e1r) * ev());
    println!("RISULTATO E4: Lyman-alpha of the rule {:.4} nm (Dirac×mu centroid {:.4} nm); period {:.4} a.u. = {:.4e} s = {:.0} steps",
             lam_nm, 1239.84198 / dref, 2.0 * std::f64::consts::PI / (e2r - e1r) * ALPHA * ALPHA, 2.0 * std::f64::consts::PI / (e2r - e1r) * ALPHA * ALPHA * 2.4188843265857e-17 * (1.0 / mu()),
             2.0 * std::f64::consts::PI / (e2r - e1r) / lv.tau);
    }
    if let Some(d) = esporta.as_ref() {
        let vk = lv.velocita_k();
        let mg = 96usize.min(n);
        let a0c = (a0() / lv.h).powi(3);
        let per = 2.0 * std::f64::consts::PI / ((es[0] + es[1] + es[2]) / 3.0 - e1r);
        let mut dtxt = String::new();
        for q in 0..fot {
            let t = per * q as f64 / fot as f64;
            let g = cartolina_z(t);
            let r: Vec<f64> = ritaglio(&lv, &dens(&g), mg).iter().map(|x| x * a0c).collect();
            scrivi(&format!("{d}/rho_lontano.f16"), &r, 1.0, q > 0);
            let j = lv.corrente(&g, &vk);
            let mut jj = vec![];
            for a in 0..3 { jj.extend(ritaglio(&lv, &j[a], mg).iter().map(|x| x * a0c / ALPHA)); }
            scrivi(&format!("{d}/j_lontano.f16"), &jj, 1.0, q > 0);
            dtxt += &format!("{:.6} {:.8}\n", t * ALPHA * ALPHA, dip(&g));
        }
        std::fs::write(format!("{d}/dipolo_lontano.txt"), dtxt).unwrap();
        println!("export: {fot} far frames over one period, crop {mg}^3 at h = a0/{:.0}", a0() / lv.h);
    }
    println!("total {:.0} s", t0.elapsed().as_secs_f64());
}

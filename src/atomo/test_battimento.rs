//! BATTIMENTO DELL'1s (docs/BATTIMENTO_1S.md): why, in the coupled run, p_exc = 1 − |⟨1s|ψ⟩|² keeps swinging
//! between ~0.2 and ~0.88 after the pulse has left. These are diagnostics that dump time series (CSV) for the
//! offline analysis in lab/battimento_analisi.py; nothing in the physics changes. All ignored (long).
//! Run: BATT_DIR=<dir> cargo test --bin atomo battimento -- --ignored --nocapture --test-threads 1

use std::io::Write;

use crate::cplx::C;
use crate::zz::Atomo;

fn dir() -> String { std::env::var("BATT_DIR").unwrap_or("/tmp/battimento".into()) }
fn env_usize(k: &str, d: usize) -> usize { std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d) }

pub fn sovrapp(a: &[[C; 4]], b: &[[C; 4]]) -> C {
    use rayon::prelude::*;
    a.par_iter().zip(b.par_iter()).map(|(x, y)| { let mut s = C::ZERO; for q in 0..4 { s += x[q].conj() * y[q]; } s })
        .reduce(|| C::ZERO, |u, v| u + v)
}

fn csv(nome: &str, testa: &str, righe: &[Vec<f64>]) {
    std::fs::create_dir_all(dir()).unwrap();
    let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{}/{nome}", dir())).unwrap());
    writeln!(f, "{testa}").unwrap();
    for r in righe { writeln!(f, "{}", r.iter().map(|v| format!("{v:.10e}")).collect::<Vec<_>>().join(",")).unwrap(); }
}

fn salva_psi(nome: &str, psi: &[[C; 4]]) {
    std::fs::create_dir_all(dir()).unwrap();
    let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{}/{nome}", dir())).unwrap());
    for p in psi { for c in p { f.write_all(&c.re.to_le_bytes()).unwrap(); f.write_all(&c.im.to_le_bytes()).unwrap(); } }
}

/// Free evolution of `start` (no light) for `ticks` ticks: per tick t, ⟨refs_k|ψ(t)⟩ (complex) and ⟨z⟩.
fn libera(start: &Atomo, refs: &[&[[C; 4]]], ticks: usize) -> Vec<Vec<f64>> {
    let mut a = start.clone();
    a.impulsi.clear();
    a.campo_est = [0.0; 3];
    let mut out = vec![];
    for _ in 0..ticks {
        let mut r = vec![a.t()];
        for f in refs { let c = sovrapp(f, &a.psi); r.push(c.re); r.push(c.im); }
        r.push(a.dipolo());
        out.push(r);
        a.avanti();
    }
    out
}


/// The Laboratorio's 1s as it was before the cure (zz::Atomo::ascolta(0.885, 0.915, asc), every tick): the
/// diagnosis is about this one.
fn vecchio_1s(x: &mut crate::accoppiato::Esperimento, asc: usize) {
    let (phi, om) = Atomo::new().ascolta(0.885, 0.915, asc);
    x.phi1s = phi.clone();
    x.omega1s = om;
    x.at.imposta(phi);
}

/// Phase 1: the Laboratorio's preparation (Esperimento::new(ASCOLTO), default 600) and the free evolution of
/// (a) the radial 1s put on the lattice (the start state of the listening) and (b) the listened 1s, no light.
#[test]
#[ignore]
fn battimento_libero() {
    let t0 = std::time::Instant::now();
    let asc = env_usize("ASCOLTO", 600);
    let ticks = env_usize("TICKS", 6000);
    let mut x = crate::accoppiato::Esperimento::new(asc);
    vecchio_1s(&mut x, asc);
    eprintln!("ascolto {asc}: ω1s {:.5}, ω2p½ {:.5}, ω21 {:.5} ({:?})", x.omega1s, x.omega2p, x.omega_21(), t0.elapsed());
    let rad = Atomo::new();
    let r = libera(&rad, &[&rad.psi, &x.phi1s], ticks);
    csv(&format!("libero_radiale.csv"), "t,c00re,c00im,c1sre,c1sim,z", &r);
    let mut a = Atomo::new();
    a.imposta(x.phi1s.clone());
    let r = libera(&a, &[&x.phi1s], ticks);
    csv(&format!("libero_1s_a{asc}.csv"), "t,cre,cim,z", &r);
    let pmin = r.iter().map(|v| v[1] * v[1] + v[2] * v[2]).fold(1.0, f64::min);
    eprintln!("free listened 1s: max p_exc {:.4} over {} ticks ({:?})", 1.0 - pmin, ticks, t0.elapsed());
    salva_psi(&format!("phi1s_a{asc}.bin"), &x.phi1s);
}

/// Phase 2: the coupled run exactly as in the export (Esperimento::new(ASCOLTO), spara(ω21, E0_ACC)): per atom
/// tick ⟨φ1s|ψ⟩, ⟨z⟩, field at the atom, light energy. At t = T_DOPO (default 500, pulse gone) the state is saved
/// and then evolved WITHOUT light (atom alone) for TICKS_DOPO ticks, recording ⟨ψ_dopo|ψ(t)⟩ and ⟨φ1s|ψ(t)⟩:
/// the spectrum of these says which levels of the 3-D rule the post-pulse state contains.
#[test]
#[ignore]
fn battimento_accoppiato() {
    use crate::accoppiato::*;
    let t0 = std::time::Instant::now();
    let asc = env_usize("ASCOLTO", 600);
    let ticks = env_usize("TICKS", 3200);
    let t_dopo = env_usize("T_DOPO", 500) as f64;
    let ticks_dopo = env_usize("TICKS_DOPO", 8000);
    let mut x = Esperimento::new(asc);
    vecchio_1s(&mut x, asc);
    let w21 = x.omega_21();
    let e00 = x.luce.energia();
    x.spara(w21, crate::zzscena::E0_ACC);
    let e0 = x.luce.energia();
    eprintln!("ω1s {:.5} ω2p {:.5} ω21 {w21:.5}; light energy {e0:.3e} (before {e00:.1e}) ({:?})", x.omega1s, x.omega2p, t0.elapsed());
    let mut rec = vec![];
    let mut dopo: Option<Atomo> = None;
    for k in 0..ticks {
        let c = sovrapp(&x.phi1s, &x.at.psi);
        let ce = x.at.campo_est;
        let el = if k % SUB == 0 { x.luce.energia() / e0 } else { f64::NAN };
        rec.push(vec![x.at.t(), c.re, c.im, x.at.dipolo(), ce[0], ce[1], ce[2], el]);
        if dopo.is_none() && x.at.t() >= t_dopo { dopo = Some(x.at.clone()); }
        x.passo();
        if k % 400 == 0 { eprintln!("tick {k} t {:.0} p_exc {:.4} E_luce/E0 {:.2e} ({:?})", x.at.t(), 1.0 - c.norm2(), el, t0.elapsed()); }
    }
    csv(&format!("accoppiato_a{asc}.csv"), "t,cre,cim,z,eEx,eEy,eEz,eluce", &rec);
    let d = dopo.unwrap();
    salva_psi(&format!("psi_dopo_a{asc}.bin"), &d.psi);
    let pd = d.psi.clone();
    let r = libera(&d, &[&pd, &x.phi1s], ticks_dopo);
    csv(&format!("dopo_libero_a{asc}.csv"), "t,cdre,cdim,c1sre,c1sim,z", &r);
    eprintln!("done ({:?})", t0.elapsed());
}

pub use crate::zz::Livello3D;
fn norma(v: &[[C; 4]]) -> f64 { sovrapp(v, v).re.sqrt() }
/// Filter diagonalisation of the 3-D rule: now `Atomo::livelli_esatti` in zz.rs (moved there after the exam).
pub fn ritz(start: &Atomo, lo: f64, hi: f64, dw: f64, ticks: usize) -> Vec<Livello3D> { start.livelli_esatti(lo, hi, dw, ticks) }

/// Sublattice of a site: the walk moves every coordinate by ±1 per tick, so (i − j) mod 2 and (j − k) mod 2 never
/// change: the cube splits into 4 sublattices that the rule never connects (the 3-D "doubling").
pub fn sottoreticolo(s: usize) -> usize {
    let n = crate::zz::N;
    let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
    ((i + j) % 2) * 2 + (j + k) % 2
}

/// Overlap of v with w restricted to each sublattice.
fn per_sottoreticolo(w: &[[C; 4]], v: &[[C; 4]]) -> [C; 4] {
    let mut o = [C::ZERO; 4];
    for (s, (a, b)) in w.iter().zip(v.iter()).enumerate() {
        let mut acc = C::ZERO;
        for q in 0..4 { acc += a[q].conj() * b[q]; }
        o[sottoreticolo(s)] += acc;
    }
    o
}

/// ⟨a| z |b⟩
fn dip(a: &[[C; 4]], b: &[[C; 4]]) -> C {
    let n = crate::zz::N;
    let mut s = C::ZERO;
    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        let z = crate::zz::x_of(i % n);
        for q in 0..4 { s += x[q].conj() * y[q] * z; }
    }
    s
}

/// Phase 3: the levels of the 3-D rule near the 1s by filter diagonalisation, their make-up (weight in the
/// radial start, in the listened 1s of the Laboratorio, per sublattice, z-coupling to the 1s), and the free
/// evolution of the clean 1s.
#[test]
#[ignore]
fn battimento_ritz() {
    let t0 = std::time::Instant::now();
    let ticks = env_usize("RITZ_TICKS", 1200);
    let rad = Atomo::new();
    let lv = ritz(&rad, 0.82, 0.94, 0.005, ticks);
    eprintln!("ritz on the radial start, {ticks} ticks: {} vectors ({:?})", lv.len(), t0.elapsed());
    let mut x = crate::accoppiato::Esperimento::new(600);
    vecchio_1s(&mut x, 600);
    let s1 = lv.iter().filter(|l| l.resid < 0.05).min_by(|a, b| (a.omega - 0.9).abs().total_cmp(&(b.omega - 0.9).abs())).unwrap();
    let mut righe = vec![];
    for l in lv.iter().filter(|l| l.resid < 0.2) {
        let w0 = sovrapp(&l.v, &rad.psi).norm2();
        let wl = sovrapp(&l.v, &x.phi1s).norm2();
        let sub = per_sottoreticolo(&rad.psi, &l.v);
        let zz = dip(&s1.v, &l.v).norm2().sqrt();
        eprintln!("  ω {:.5} resid {:.1e}  |⟨rad|v⟩|² {:.4}  |⟨φ1s_lab|v⟩|² {:.4}  |⟨1s|z|v⟩| {:.3}  sub ⟨rad|v⟩ {:?}",
            l.omega, l.resid, w0, wl, zz, sub.iter().map(|c| format!("{:+.3}{:+.3}i", c.re, c.im)).collect::<Vec<_>>());
        righe.push(vec![l.omega, l.resid, w0, wl, zz]);
    }
    csv("ritz_livelli.csv", "omega,resid,w_rad,w_lab,z_1s", &righe);
    // free evolution of the clean 1s
    let mut a = Atomo::new();
    a.imposta(s1.v.clone());
    let r = libera(&a, &[&s1.v], 3000);
    let pmax = r.iter().map(|v| 1.0 - (v[1] * v[1] + v[2] * v[2])).fold(0.0, f64::max);
    eprintln!("clean 1s ω {:.5} resid {:.2e}: free max p_exc over 3000 ticks {:.2e} ({:?})", s1.omega, s1.resid, pmax, t0.elapsed());
    csv("libero_1s_ritz.csv", "t,cre,cim,z", &r);
    salva_psi("phi1s_ritz.bin", &s1.v);
}

fn carica_psi(nome: &str) -> Vec<[C; 4]> {
    let b = std::fs::read(format!("{}/{nome}", dir())).unwrap();
    let f: Vec<f64> = b.chunks(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect();
    f.chunks(8).map(|c| [C::new(c[0], c[1]), C::new(c[2], c[3]), C::new(c[4], c[5]), C::new(c[6], c[7])]).collect()
}

/// Shape of a wave: ⟨r⟩, weight inside r < 1 λ_C, "roughness" Σ|v(s) − v(s + x̂)|² / Σ|v|² (≈ 0 smooth, 4 staggered).
fn forma(v: &[[C; 4]]) -> (f64, f64, f64) {
    let n = crate::zz::N;
    let (mut r1, mut w1, mut rug, mut tot) = (0.0, 0.0, 0.0, 0.0);
    for (s, p) in v.iter().enumerate() {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        let x = [crate::zz::x_of(i), crate::zz::x_of(j), crate::zz::x_of(k)];
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        let w: f64 = p.iter().map(|c| c.norm2()).sum();
        tot += w; r1 += w * r; if r < 1.0 { w1 += w; }
        if i + 2 < n {
            // same sublattice neighbour two cells along x
            let q = &v[s + 2 * n * n];
            rug += (0..4).map(|c| (p[c] - q[c]).norm2()).sum::<f64>();
        }
    }
    (r1 / tot, w1 / tot, rug / tot)
}

/// Phase 3b: what the levels near the 1s are (shape), and how the Laboratorio's listened 1s and the post-pulse
/// state decompose on the exact levels of the rule.
#[test]
#[ignore]
fn battimento_carattere() {
    let t0 = std::time::Instant::now();
    let rad = Atomo::new();
    let mut x = crate::accoppiato::Esperimento::new(600);
    vecchio_1s(&mut x, 600);
    let mut lab = Atomo::new();
    lab.imposta(x.phi1s.clone());
    let dopo = carica_psi("psi_dopo_a600.bin");
    let mut ad = Atomo::new();
    ad.imposta(dopo.clone());
    for (nome, st) in [("radiale", &rad), ("phi1s_lab", &lab), ("psi_dopo", &ad)] {
        let lv = ritz(st, 0.80, 1.00, 0.005, 1600);
        let mut cap = 0.0;
        eprintln!("== start {nome} ({:?})", t0.elapsed());
        for l in lv.iter().filter(|l| l.resid < 1e-3) {
            let w = sovrapp(&l.v, &st.psi).norm2();
            cap += w;
            if w < 1e-3 { continue; }
            let (r, w1, rug) = forma(&l.v);
            eprintln!("  ω {:.5} resid {:.1e} weight {:.4}  ⟨r⟩ {:.3} w(r<1) {:.3} rough {:.3}", l.omega, l.resid, w, r, w1, rug);
        }
        eprintln!("  total weight on exact levels in [0.80, 1.00]: {cap:.4}");
    }
    let (r, w1, rug) = forma(&rad.psi);
    eprintln!("radial start: ⟨r⟩ {r:.3} w(r<1) {w1:.3} rough {rug:.3}");
}

#[test]
#[ignore]
fn battimento_ruvido() {
    let mut a = Atomo::new();
    for t in 0..8 { let f = forma(&a.psi); eprintln!("tick {t}: ⟨r⟩ {:.3} rough {:.3}", f.0, f.2); a.avanti(); }
    let mut ph = a.psi.clone();
    // per-component roughness
    for q in 0..4 {
        let n = crate::zz::N;
        let (mut r, mut t) = (0.0, 0.0);
        for s in 0..ph.len() { let i = s / (n * n); t += ph[s][q].norm2(); if i + 2 < n { r += (ph[s][q] - ph[s + 2 * n * n][q]).norm2(); } }
        eprintln!("component {q}: weight {t:.3} rough {:.3}", r / t);
    }
    ph.clear();
}

#[test]
#[ignore]
fn battimento_linea() {
    let rad = Atomo::new();
    let lv = ritz(&rad, 0.82, 0.94, 0.005, 1200);
    let n = crate::zz::N;
    let mut sel: Vec<&Livello3D> = lv.iter().filter(|l| l.resid < 1e-6 && sovrapp(&l.v, &rad.psi).norm2() > 0.01).collect();
    sel.sort_by(|a, b| a.omega.total_cmp(&b.omega));
    // projection of the start on the 1s pair (the two levels near 0.8995)
    let mut pair = vec![[C::ZERO; 4]; rad.psi.len()];
    for l in sel.iter().filter(|l| (l.omega - 0.8995).abs() < 0.001) {
        let c = sovrapp(&l.v, &rad.psi);
        for (x, y) in pair.iter_mut().zip(l.v.iter()) { for q in 0..4 { x[q] += c * y[q]; } }
    }
    let np = norma(&pair);
    for x in pair.iter_mut() { for c in x.iter_mut() { *c = c.scale(1.0 / np); } }
    let fp = forma(&pair);
    eprintln!("1s pair projection: ⟨r⟩ {:.3} rough {:.3}", fp.0, fp.2);
    for l in &sel {
        let c = sovrapp(&l.v, &rad.psi);
        let ph = c.scale(1.0 / c.norm2().sqrt()).conj();
        let line: Vec<String> = (16..32).map(|i| { let s = (i * n + 24) * n + 24; let a = l.v[s][0] * ph; let r: f64 = l.v[s].iter().map(|c| c.norm2()).sum();
            format!("{:+.0}°/{:.2e}", a.arg().to_degrees(), r) }).collect();
        eprintln!("ω {:.5}: {}", l.omega, line.join(" "));
    }
    let line: Vec<String> = (16..32).map(|i| { let s = (i * n + 24) * n + 24; let r: f64 = rad.psi[s].iter().map(|c| c.norm2()).sum(); format!("{:+.0}°/{:.2e}", rad.psi[s][0].arg().to_degrees(), r) }).collect();
    eprintln!("radial : {}", line.join(" "));
    // free evolution of the pair projection (even ticks)
    let mut a = Atomo::new();
    a.imposta(pair.clone());
    let r = libera(&a, &[&pair], 3200);
    let pmax = r.iter().step_by(2).map(|v| 1.0 - (v[1] * v[1] + v[2] * v[2])).fold(0.0, f64::max);
    eprintln!("pair projection: free max p_exc (even ticks, 3200 ticks) {pmax:.3e}");
}

/// Weight of a wave on each parity class (i mod 2, j mod 2, k mod 2) and on each residue i mod 4 (x only).
fn classi(v: &[[C; 4]]) -> ([f64; 8], [f64; 4]) {
    let n = crate::zz::N;
    let (mut c8, mut c4) = ([0.0; 8], [0.0; 4]);
    for (s, p) in v.iter().enumerate() {
        let (i, j, k) = (s / (n * n), (s / n) % n, s % n);
        let w: f64 = p.iter().map(|c| c.norm2()).sum();
        c8[(i % 2) * 4 + (j % 2) * 2 + k % 2] += w;
        c4[i % 4] += w;
    }
    (c8, c4)
}

#[test]
#[ignore]
fn battimento_classi() {
    let rad = Atomo::new();
    let lv = ritz(&rad, 0.82, 0.94, 0.005, 1200);
    let mut sel: Vec<&Livello3D> = lv.iter().filter(|l| l.resid < 1e-6 && sovrapp(&l.v, &rad.psi).norm2() > 0.01).collect();
    sel.sort_by(|a, b| a.omega.total_cmp(&b.omega));
    for l in &sel {
        let (c8, c4) = classi(&l.v);
        // spinor make-up: weight on each component (L0, L1, R0, R1)
        let mut sp = [0.0; 4];
        for p in &l.v { for q in 0..4 { sp[q] += p[q].norm2(); } }
        eprintln!("ω {:.5}: parity classes {:?}  i mod 4 {:?}  spinor {:?}", l.omega,
            c8.iter().map(|x| format!("{x:.3}")).collect::<Vec<_>>(), c4.iter().map(|x| format!("{x:.3}")).collect::<Vec<_>>(),
            sp.iter().map(|x| format!("{x:.3}")).collect::<Vec<_>>());
    }
}

#[test]
#[ignore]
fn battimento_salva_livelli() {
    let rad = Atomo::new();
    let lv = ritz(&rad, 0.82, 0.94, 0.005, 1200);
    let mut sel: Vec<&Livello3D> = lv.iter().filter(|l| l.resid < 1e-6 && sovrapp(&l.v, &rad.psi).norm2() > 0.01).collect();
    sel.sort_by(|a, b| a.omega.total_cmp(&b.omega));
    for l in &sel { salva_psi(&format!("liv_{:.5}.bin", l.omega), &l.v); }
    salva_psi("radiale.bin", &rad.psi);
}

/// The pre-registered choice (lab/out/battimento_prereg.json): among the exact levels near the 1s (residual
/// < 1e-6, weight > 1 % in the radial start), the one with the least weight inside r < 1 λ_C (least touched
/// by the first cells around the nucleus). Returns (chosen, all of them sorted by ω).
pub fn scegli_1s(rad: &Atomo) -> (Livello3D, Vec<Livello3D>) { rad.uno_s_esatto(1200) }

/// PRE-REGISTERED EXAM B1–B4 (lab/out/battimento_prereg.json, committed before this was run).
#[test]
#[ignore]
fn battimento_esame() {
    use crate::accoppiato::*;
    let t0 = std::time::Instant::now();
    let rad = Atomo::new();
    let (s1, altri) = scegli_1s(&rad);
    eprintln!("chosen 1s ω {:.5} (resid {:.1e}, w(r<1) {:.3}); others {:?} ({:?})", s1.omega, s1.resid, forma(&s1.v).1,
        altri.iter().map(|l| format!("{:.5}", l.omega)).collect::<Vec<_>>(), t0.elapsed());
    // B1
    let mut a = Atomo::new();
    a.imposta(s1.v.clone());
    let r = libera(&a, &[&s1.v], 3200);
    let b1 = r.iter().step_by(crate::accoppiato::SUB).map(|v| 1.0 - (v[1] * v[1] + v[2] * v[2])).fold(0.0, f64::max);
    eprintln!("B1: free max p_exc {b1:.2e} (threshold 1e-3) -> {}", if b1 < 1e-3 { "PASS" } else { "FAIL" });
    // the coupled run with the new 1s
    let mut x = Esperimento::new(600);
    x.phi1s = s1.v.clone();
    x.omega1s = s1.omega;
    x.at.imposta(s1.v.clone());
    let w21 = x.omega_21();
    x.spara(w21, crate::zzscena::E0_ACC);
    eprintln!("ω1s {:.5} ω2p½ {:.5} ω21 {w21:.5} ({:?})", x.omega1s, x.omega2p, t0.elapsed());
    for _ in 0..3200 { x.passo(); }
    let fin: Vec<f64> = x.serie_p.iter().filter(|p| p.0 >= 450.0 && p.0 <= 1600.0).map(|p| p.1).collect();
    let (mx, mn) = (fin.iter().cloned().fold(0.0, f64::max), fin.iter().cloned().fold(1.0, f64::min));
    let media = fin.iter().sum::<f64>() / fin.len() as f64;
    eprintln!("B2: p_exc in [450,1600] max {mx:.4} min {mn:.4} swing {:.4} (threshold 0.05) -> {}", mx - mn, if mx - mn < 0.05 { "PASS" } else { "FAIL" });
    let b3: f64 = altri.iter().map(|l| sovrapp(&l.v, &x.at.psi).norm2()).sum();
    eprintln!("B3: weight on the other multiplet levels at t = {:.0}: {b3:.2e} (threshold 0.01) -> {}", x.at.t(), if b3 < 0.01 { "PASS" } else { "FAIL" });
    eprintln!("B4: mean p_exc {media:.4} (threshold > 0.01) -> {}", if media > 0.01 { "PASS" } else { "FAIL" });
    let fd = frequenza(&x.serie_d.iter().filter(|p| p.0 >= 450.0).cloned().collect::<Vec<_>>(), 450.0);
    eprintln!("info: dipole frequency after the pulse {:?}; 2p3/2 (0.97880, listened) − ω1s = {:.5}; ω1s vs Dirac 0.89907: {:+.3} %",
        fd, 0.97880 - s1.omega, 100.0 * (s1.omega / 0.89907 - 1.0));
    let rows: Vec<Vec<f64>> = x.serie_p.iter().map(|p| vec![p.0, p.1]).collect();
    csv("esame_p.csv", "t,p_exc", &rows);
    let rows: Vec<Vec<f64>> = x.serie_d.iter().map(|p| vec![p.0, p.1]).collect();
    csv("esame_d.csv", "t,z", &rows);
    eprintln!("done ({:?})", t0.elapsed());
}

/// B1-type check for the scene WITHOUT coupling (zzscena): the atom `Atomo::new()` restarted from the clean 1s of
/// `Esperimento::new(600)` (as the scene does when the experiment is ready), no light, 3200 ticks: max p_exc on the
/// ticks the scene shows it (every 4) < 1e-3 (threshold of B1, lab/out/battimento_prereg.json).
#[test]
#[ignore]
fn battimento_scena_libera() {
    let x = crate::accoppiato::Esperimento::new(600);
    let mut at = Atomo::new();
    at.imposta(x.phi1s.clone());
    let mut pmax: f64 = 0.0;
    for _ in 0..3200 {
        if at.ticks % 4 == 0 { pmax = pmax.max(1.0 - at.p_fondamentale()); }
        at.avanti();
    }
    eprintln!("scene without coupling, clean 1s ω {:.5}: max p_exc {pmax:.2e} over 3200 ticks", x.omega1s);
    assert!(pmax < 1e-3);
}

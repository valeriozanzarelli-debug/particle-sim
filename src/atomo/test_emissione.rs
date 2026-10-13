//! EMISSIONE NEL LABORATORIO (docs/ATOMO.md §3.7, from docs/EMISSIONE.md §3.2 and §7): the timing of the light ↔
//! atom coupling of `accoppiato.rs`. The original coupling reads the light field before the light tick and deposits
//! after it, not transverse: a fake friction (×~2.2 in the Python reference). The patch: centred field + transverse
//! deposit (`Accoppiamento::CENTRATO`). Pre-registration: lab/out/emissione_lab_prereg.json.
//! Fast tests: `cargo test --bin atomo emissione`. Exam: EMLAB_DIR=<dir> cargo test --bin atomo emissione_esame
//! -- --ignored --nocapture --test-threads 1

use crate::accoppiato::*;
use crate::cplx::C;

fn env_usize(k: &str, d: usize) -> usize { std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d) }

/// The light rule's two-tick rotation: its axis is the longitudinal mode (→ k̂ for small k), P_⊥ kills it.
#[test]
fn emissione_proiettore_trasverso() {
    for k in [[0.01, 0.02, -0.015], [0.3, -0.1, 0.2], [1.0, 0.5, -2.0], [0.0, 0.0, 0.2]] {
        let r = rotazione_luce(k);
        let p = proiettore_trasverso(k);
        // P is a projector of rank 2: P² = P, trace 2
        let tr = p[0][0] + p[1][1] + p[2][2];
        assert!((tr - 2.0).abs() < 1e-12);
        for i in 0..3 { for l in 0..3 {
            let pp: f64 = (0..3).map(|m| p[i][m] * p[m][l]).sum();
            assert!((pp - p[i][l]).abs() < 1e-12);
        }}
        // the axis (1 − P) is left fixed by R
        let n2: Vec<[f64; 3]> = (0..3).map(|c| [0, 1, 2].map(|i| if i == c { 1.0 } else { 0.0 } - p[i][c])).collect();
        for v in n2 { let rv = [0, 1, 2].map(|i| (0..3).map(|m| r[i][m] * v[m]).sum::<f64>()); for i in 0..3 { assert!((rv[i] - v[i]).abs() < 1e-9); } }
        let kn = (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]).sqrt();
        if kn < 0.05 {
            for i in 0..3 { for l in 0..3 {
                let cont = if i == l { 1.0 } else { 0.0 } - k[i] * k[l] / (kn * kn);
                assert!((p[i][l] - cont).abs() < 0.02, "small k: P_perp → 1 − k̂k̂ᵀ");
            }}
        }
    }
}

fn dir() -> String { std::env::var("EMLAB_DIR").unwrap_or("/tmp/emissione_lab".into()) }
fn sovrapp(a: &[[C; 4]], b: &[[C; 4]]) -> C {
    use rayon::prelude::*;
    a.par_iter().zip(b.par_iter()).map(|(x, y)| { let mut s = C::ZERO; for q in 0..4 { s += x[q].conj() * y[q]; } s }).reduce(|| C::ZERO, |u, v| u + v)
}

/// The exam's states: the exact 1s of `Esperimento::new(600)`, and the exact levels the 1s reaches with light along z
/// (`Atomo::livelli_di_dipolo` from z·1s, [0.966, 0.986], dw 0.0015, 6000 ticks): 2p3/2* = the strongest |d_z|² in
/// [0.977, 0.981], 2p½* = the strongest in [0.9745, 0.9775]. Cached in EMLAB_DIR/stati.bin (not in the repo).
pub struct Stati { pub s1: Vec<[C; 4]>, pub w1: f64, pub p32: Vec<[C; 4]>, pub w32: f64, pub d32: [C; 3], pub p12: Vec<[C; 4]>, pub w12: f64, pub d12: [C; 3], pub res32: f64, pub res12: f64 }

fn scrivi_v(f: &mut Vec<u8>, v: &[[C; 4]]) { for p in v { for c in p { f.extend_from_slice(&c.re.to_le_bytes()); f.extend_from_slice(&c.im.to_le_bytes()); } } }
fn leggi_v(b: &[u8], off: &mut usize, n: usize) -> Vec<[C; 4]> {
    let mut f = || { let x = f64::from_le_bytes(b[*off..*off + 8].try_into().unwrap()); *off += 8; x };
    (0..n).map(|_| { let mut p = [C::ZERO; 4]; for c in p.iter_mut() { let re = f(); let im = f(); *c = C::new(re, im); } p }).collect()
}

fn nucleo_env() -> crate::zz::Nucleo { if std::env::var("NUCLEO").as_deref() == Ok("imparato") { crate::zz::Nucleo::Imparato } else { crate::zz::Nucleo::Nudo } }

pub fn stati() -> Stati {
    let nu = nucleo_env();
    let n = crate::zz::n_siti(None);
    let file = if nu == crate::zz::Nucleo::Nudo { format!("{}/stati.bin", dir()) } else { format!("{}/stati_imparato.bin", dir()) };
    if let Ok(b) = std::fs::read(&file) {
        let mut off = 0;
        let mut f = |o: &mut usize| { let x = f64::from_le_bytes(b[*o..*o + 8].try_into().unwrap()); *o += 8; x };
        let (w1, w32, w12, res32, res12) = (f(&mut off), f(&mut off), f(&mut off), f(&mut off), f(&mut off));
        let s1 = leggi_v(&b, &mut off, n);
        let p32 = leggi_v(&b, &mut off, n);
        let p12 = leggi_v(&b, &mut off, n);
        let d32 = crate::zz::dipolo_tra(&s1, &p32, None);
        let d12 = crate::zz::dipolo_tra(&s1, &p12, None);
        return Stati { s1, w1, p32, w32, d32, p12, w12, d12, res32, res12 };
    }
    let x = Esperimento::nuovo_nucleo(600, None, nu);
    let lv = x.at.livelli_di_dipolo(&x.phi1s, 2, 0.966, 0.986, 0.0015, 6000, 1e-3);
    for l in &lv { eprintln!("  level ω {:.6} resid {:.1e} |d_z|² {:.4}", l.lv.omega, l.lv.resid, l.d[2].norm2()); }
    let p32 = lv.iter().find(|l| l.lv.omega > 0.977 && l.lv.omega < 0.981).expect("no 2p3/2 level");
    let p12 = lv.iter().find(|l| l.lv.omega > 0.9745 && l.lv.omega < 0.9775).expect("no 2p1/2 level");
    let mut b = vec![];
    for v in [x.omega1s, p32.lv.omega, p12.lv.omega, p32.lv.resid, p12.lv.resid] { b.extend_from_slice(&v.to_le_bytes()); }
    scrivi_v(&mut b, &x.phi1s);
    scrivi_v(&mut b, &p32.lv.v);
    scrivi_v(&mut b, &p12.lv.v);
    std::fs::create_dir_all(dir()).unwrap();
    std::fs::write(&file, b).unwrap();
    Stati { s1: x.phi1s.clone(), w1: x.omega1s, p32: p32.lv.v.clone(), w32: p32.lv.omega, d32: p32.d, p12: p12.lv.v.clone(), w12: p12.lv.omega, d12: p12.d, res32: p32.lv.resid, res12: p12.lv.resid }
}

/// The coupled experiment on the exam's states (1s exact, 2p½* as the Laboratorio's "2p": the pulse frequency).
fn esperimento(st: &Stati, acc: Accoppiamento) -> Esperimento {
    Esperimento::da_stati(st.s1.clone(), st.w1, st.p12.clone(), st.w12, nucleo_env()).con_accoppiamento(acc)
}

/// (1s + 2p3/2*)/√2.
fn sovrapposizione(st: &Stati) -> Vec<[C; 4]> {
    let k = std::f64::consts::FRAC_1_SQRT_2;
    st.s1.iter().zip(st.p32.iter()).map(|(a, b)| [0, 1, 2, 3].map(|q| (a[q] + b[q]).scale(k))).collect()
}

/// Linear fit (slope) of y(t).
fn pendenza(v: &[(f64, f64)]) -> f64 { crate::fotscena::retta(v).1 }

/// Damping in windows: (t1, t2, A measured = −d ln(p/(1−p))/dt, A Larmor of the measured dipole at ω, ratio,
/// light energy gained, atom energy lost ω Δp, their ratio).
pub fn smorzamento(x: &Esperimento, w: f64, k: f64, finestre: &[(f64, f64)]) -> Vec<[f64; 8]> {
    finestre.iter().map(|&(t1, t2)| {
        let lp: Vec<(f64, f64)> = x.serie_p.iter().filter(|q| q.0 >= t1 && q.0 <= t2).map(|q| (q.0, (q.1 / (1.0 - q.1)).ln())).collect();
        let am = -pendenza(&lp);
        let ps: Vec<f64> = x.serie_p.iter().filter(|q| q.0 >= t1 && q.0 <= t2).map(|q| q.1).collect();
        let pm = ps.iter().sum::<f64>() / ps.len() as f64;
        let d = ampiezza_a(&x.serie_d, w, t1, t2);
        let ap = gamma_larmor(w, d, pm) * (k / K_ACC).powi(2);
        let at = |t: f64| x.serie_p.iter().min_by(|a, b| (a.0 - t).abs().total_cmp(&(b.0 - t).abs())).unwrap().1;
        let el = |t: f64| x.luce.serie_e.iter().min_by(|a, b| (a.0 - t).abs().total_cmp(&(b.0 - t).abs())).unwrap().1;
        let (de_l, de_a) = (el(t2) - el(t1), w * (at(t1) - at(t2)));
        [t1, t2, am, ap, am / ap, de_l, de_a, de_l / de_a]
    }).collect()
}

fn json_righe(r: &[[f64; 8]]) -> String {
    let v: Vec<String> = r.iter().map(|q| format!("{{\"t1\": {}, \"t2\": {}, \"A_misurato\": {:.6e}, \"A_larmor\": {:.6e}, \"rapporto\": {:.5}, \"dE_luce\": {:.6e}, \"dE_atomo\": {:.6e}, \"luce_su_atomo\": {:.5}}}", q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7])).collect();
    format!("[{}]", v.join(", "))
}

fn salva(nome: &str, testo: &str) {
    std::fs::create_dir_all(dir()).unwrap();
    let nome = &if nucleo_env() == crate::zz::Nucleo::Imparato { format!("imparato_{nome}") } else { nome.to_string() };
    std::fs::write(format!("{}/{nome}", dir()), testo).unwrap();
    eprintln!("{nome}: {testo}");
}

fn serie_csv(nome: &str, x: &Esperimento) {
    use std::io::Write;
    let nome = &if nucleo_env() == crate::zz::Nucleo::Imparato { format!("imparato_{nome}") } else { nome.to_string() };
    let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{}/{nome}", dir())).unwrap());
    writeln!(f, "t,p_exc,z_fine_tick,energia_luce_piu_assorbita,sonda_Ez").unwrap();
    for (i, (t, p)) in x.serie_p.iter().enumerate() {
        writeln!(f, "{t:.3},{p:.9e},{:.9e},{:.9e},{:.9e}", x.serie_d[(i + 1) * SUB - 1].1, x.luce.serie_e[i].1, x.luce.sonda[i].1).unwrap();
    }
}

/// Exam D (pre-registered, lab/out/emissione_lab_prereg.json): damping of (1s + 2p3/2*)/√2 with K = 10, T = 2000,
/// windows [200, 800], [800, 1400], [1400, 2000], for the coupling ACC (ritardato | centrato | centrato+trasverso).
#[test]
#[ignore]
fn emissione_esame_d() {
    let st = stati();
    let w = st.w32 - st.w1;
    let acc = match std::env::var("ACC").as_deref() {
        Ok("ritardato") => Accoppiamento::RITARDATO,
        Ok("centrato") => Accoppiamento { centrato: true, trasverso: false },
        _ => Accoppiamento::CENTRATO,
    };
    let tmax: f64 = std::env::var("TMAX").ok().and_then(|v| v.parse().ok()).unwrap_or(2000.0);
    let t0 = std::time::Instant::now();
    let mut x = esperimento(&st, acc);
    x.parti_da(sovrapposizione(&st));
    while x.at.t() < tmax - 1e-9 { x.passo(); }
    let fin: Vec<(f64, f64)> = if tmax >= 2000.0 { vec![(200.0, 800.0), (800.0, 1400.0), (1400.0, 2000.0)] } else { vec![(100.0, tmax)] };
    let r = smorzamento(&x, w, K_ACC, &fin);
    let m = r.iter().map(|q| q[4]).sum::<f64>() / r.len() as f64;
    let tot = smorzamento(&x, w, K_ACC, &[(fin[0].0, fin[fin.len() - 1].1)]);
    serie_csv(&format!("d_{}.csv", acc.nome()), &x);
    salva(&format!("d_{}.json", acc.nome()), &format!("{{\"accoppiamento\": \"{}\", \"omega\": {w:.6}, \"d32\": [{:.5}, {:.5}, {:.5}], \"W2\": {:.5}, \"finestre\": {}, \"rapporto_medio\": {m:.5}, \"intero\": {}, \"p_fine\": {:.6}, \"secondi\": {:.0}}}",
        acc.nome(), st.d32[0].abs(), st.d32[1].abs(), st.d32[2].abs(), LuceAtomo::fattore_forma(w).powi(2), json_righe(&r), json_righe(&tot), x.serie_p.last().unwrap().1, t0.elapsed().as_secs_f64()));
}

/// Exam F (still 1s): `Esperimento::new(600)` with coupling ACC, no pulse, 3200 atom ticks: max p_exc.
#[test]
#[ignore]
fn emissione_esame_fermo() {
    for acc in [Accoppiamento::RITARDATO, Accoppiamento::CENTRATO] {
        let st = stati();
        let mut x = esperimento(&st, acc);
        x.parti_da(st.s1.clone());
        let mut pmax: f64 = 0.0;
        for _ in 0..3200 { x.passo(); if let Some(p) = x.serie_p.last() { pmax = pmax.max(p.1); } }
        salva(&format!("fermo_{}.json", acc.nome()), &format!("{{\"accoppiamento\": \"{}\", \"p_exc_max\": {pmax:.3e}, \"energia_luce\": {:.3e}}}", acc.nome(), x.luce.serie_e.last().unwrap().1));
    }
}

/// Exam L (line): the Laboratorio experiment (1s exact, spara(omega_21(), E0_ACC)) with coupling ACC, listened at the
/// probe for t in [400, TL] (TL = 8000): strongest and second line in [0.070, 0.090] (Hann, step 1e-5, parabolic peak)
/// against ω(2p3/2*) − ω1s and ω(2p½*) − ω1s.
#[test]
#[ignore]
fn emissione_esame_riga() {
    let st = stati();
    let acc = if std::env::var("ACC").as_deref() == Ok("ritardato") { Accoppiamento::RITARDATO } else { Accoppiamento::CENTRATO };
    let tl: f64 = std::env::var("TL").ok().and_then(|v| v.parse().ok()).unwrap_or(8000.0);
    let t0 = std::time::Instant::now();
    let mut x = esperimento(&st, acc);
    x.parti_da(st.s1.clone());
    let w21 = x.omega_21();
    x.spara(w21, crate::zzscena::E0_ACC);
    let mut pmax: f64 = 0.0;
    while x.at.t() < tl - 1e-9 { x.passo(); if let Some(p) = x.serie_p.last() { pmax = pmax.max(p.1); } }
    let spet: Vec<(f64, f64)> = (0..=2000).map(|k| { let om = 0.070 + 1e-5 * k as f64; (om, ampiezza_a(&x.luce.sonda, om, 400.0, tl)) }).collect();
    let mut picchi: Vec<(f64, f64)> = vec![];
    for i in 1..spet.len() - 1 {
        if spet[i].1 > spet[i - 1].1 && spet[i].1 >= spet[i + 1].1 {
            let (y0, y1, y2) = (spet[i - 1].1, spet[i].1, spet[i + 1].1);
            let den = y0 - 2.0 * y1 + y2;
            let om = spet[i].0 + if den.abs() > 0.0 { 1e-5 * 0.5 * (y0 - y2) / den } else { 0.0 };
            picchi.push((om, y1));
        }
    }
    picchi.sort_by(|a, b| b.1.total_cmp(&a.1));
    let primo = picchi[0];
    let secondo = picchi.iter().find(|p| (p.0 / primo.0 - 1.0).abs() > 0.01).cloned().unwrap_or((0.0, 0.0));
    let (a32, a12) = (st.w32 - st.w1, st.w12 - st.w1);
    serie_csv(&format!("riga_{}.csv", acc.nome()), &x);
    let tp = x.serie_p.iter().max_by(|a, b| a.1.total_cmp(&b.1)).unwrap().0;
    let pf = x.serie_p.last().unwrap().1;
    salva(&format!("riga_{}.json", acc.nome()), &format!("{{\"accoppiamento\": \"{}\", \"w21_lab\": {w21:.6}, \"atteso_32\": {a32:.6}, \"atteso_12\": {a12:.6}, \"primo\": [{:.6}, {:.4e}], \"secondo\": [{:.6}, {:.4e}], \"scarto_primo\": {:.5}, \"scarto_secondo\": {:.5}, \"rapporto_ampiezze\": {:.4}, \"picchi\": {:?}, \"p_exc_max\": {pmax:.5}, \"t_picco\": {tp:.1}, \"p_exc_fine\": {pf:.5}, \"secondi\": {:.0}}}",
        acc.nome(), primo.0, primo.1, secondo.0, secondo.1, primo.0 / a32 - 1.0, secondo.0 / a12 - 1.0, primo.1 / secondo.1.max(1e-300),
        picchi.iter().take(6).map(|p| (format!("{:.5}", p.0), format!("{:.2e}", p.1))).collect::<Vec<_>>(), t0.elapsed().as_secs_f64()));
}

/// Exam K1 (true charge): (1s + 2p3/2*)/√2 with K = 1, centred+transverse coupling, T = TK (4000); and the same state
/// without light (free walk) to measure the drift of p that the light has nothing to do with.
#[test]
#[ignore]
fn emissione_esame_k1() {
    let st = stati();
    let w = st.w32 - st.w1;
    let tk: f64 = std::env::var("TK").ok().and_then(|v| v.parse().ok()).unwrap_or(4000.0);
    let t0 = std::time::Instant::now();
    // free drift
    let mut a = crate::zz::Atomo::new();
    a.con_nucleo(nucleo_env());
    a.imposta(sovrapposizione(&st));
    let mut lib = vec![];
    while a.t() < tk - 1e-9 { a.avanti(); if a.ticks % SUB as u64 == 0 { let p = 1.0 - sovrapp(&st.s1, &a.psi).norm2(); lib.push((a.t(), (p / (1.0 - p)).ln())); } }
    let deriva = -pendenza(&lib.iter().filter(|q| q.0 >= 200.0).cloned().collect::<Vec<_>>());
    let mut x = esperimento(&st, Accoppiamento::CENTRATO);
    x.k_acc = 1.0;
    x.parti_da(sovrapposizione(&st));
    while x.at.t() < tk - 1e-9 { x.passo(); }
    let r = smorzamento(&x, w, 1.0, &[(200.0, tk)]);
    serie_csv("k1.csv", &x);
    salva("k1.json", &format!("{{\"omega\": {w:.6}, \"intero\": {}, \"deriva_libera\": {deriva:.4e}, \"deriva_su_atteso\": {:.4}, \"p_fine\": {:.7}, \"secondi\": {:.0}}}",
        json_righe(&r), deriva / r[0][3], x.serie_p.last().unwrap().1, t0.elapsed().as_secs_f64()));
}

/// Period-4 checkerboard of the density (the 4 valleys): Σ |ρ̂(k)|² over the wave vectors with every component in
/// {0, ±12, 24} (period 4 cells on 48) except k = 0, over |ρ̂(0)|².
pub fn scacchiera(v: &[[C; 4]]) -> f64 {
    use rayon::prelude::*;
    let n = crate::zz::N;
    let rho: Vec<f64> = v.iter().map(|p| p.iter().map(|c| c.norm2()).sum()).collect();
    let ks = [0usize, 12, 24, 36];
    let tw = std::f64::consts::TAU / n as f64;
    let hat = |k: [usize; 3]| -> C {
        rho.par_iter().enumerate().map(|(s, r)| { let c = crate::zz::sito(s, None); C::expi(-tw * (k[0] * c[0] + k[1] * c[1] + k[2] * c[2]) as f64).scale(*r) })
            .reduce(|| C::ZERO, |a, b| a + b)
    };
    let mut s = 0.0;
    for &a in &ks { for &b in &ks { for &c in &ks {
        if a + b + c == 0 { continue; }
        s += hat([a, b, c]).norm2();
    }}}
    s / hat([0, 0, 0]).norm2()
}

/// Development (declared, amendment 1): the learned nucleus in the Laboratorio's lattice: exact 1s levels, which one
/// uno_s_esatto takes, checkerboard numbers.
#[test]
#[ignore]
fn emissione_sviluppo_nucleo() {
    use crate::zz::{Atomo, Nucleo};
    for nu in [Nucleo::Nudo, Nucleo::Imparato] {
        let t0 = std::time::Instant::now();
        let mut a = Atomo::new();
        a.con_nucleo(nu);
        eprintln!("{}: radial 1s checkerboard {:.3e}", nu.nome(), scacchiera(&a.psi));
        let lv = a.livelli_esatti(0.82, 0.94, 0.005, env_usize("TICKS1S", 1200));
        for l in lv.iter().take(8) { eprintln!("  ω {:.6} resid {:.1e} weight {:.4} inner {:.4} checker {:.3e}", l.omega, l.resid, sovrapp(&l.v, &a.psi).norm2(), crate::zz::peso_interno(&l.v, 1.0, None), scacchiera(&l.v)); }
        eprintln!("  ({:?})", t0.elapsed());
    }
}

/// Exam NC (amendment 1): the period-4 checkerboard of the 1s used, bare (exact level) vs learned (smooth group), and of
/// the radial start state (the floor of the measure).
#[test]
#[ignore]
fn emissione_esame_scacchiera() {
    std::env::set_var("NUCLEO", "nudo");
    let b = stati();
    std::env::set_var("NUCLEO", "imparato");
    let i = stati();
    let r = crate::zz::Atomo::new();
    let (cb, ci, cr) = (scacchiera(&b.s1), scacchiera(&i.s1), scacchiera(&r.psi));
    let (pb, pi) = (scacchiera(&b.p32), scacchiera(&i.p32));
    salva("scacchiera.json", &format!("{{\"C4_1s_nudo\": {cb:.4e}, \"C4_1s_imparato\": {ci:.4e}, \"C4_radiale\": {cr:.4e}, \"rapporto\": {:.4e}, \"C4_2p32_nudo\": {pb:.4e}, \"C4_2p32_imparato\": {pi:.4e}, \"w1s_imparato\": {:.6}, \"w32\": {:.6}, \"w12\": {:.6}, \"d32\": {:.4}, \"d12\": {:.4}}}", ci / cb, i.w1, i.w32, i.w12, i.d32[2].norm2(), i.d12[2].norm2()));
}

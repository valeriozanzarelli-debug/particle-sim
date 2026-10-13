//! ONE SUBLATTICE ("colore") of the 3-D zig-zag atom (docs/SOTTORETICOLO.md). The rule moves every amplitude one
//! cell along x, y and z each tick, so (i+j) mod 2 and (j+k) mod 2 are conserved: the cube is 4 lattices that never
//! exchange amplitude. `Atomo` with `colore = Some(c)` stores and evolves only one of them (N³/4 sites).
//! Fast test: `cargo test --bin atomo sottoreticolo`. Exam (pre-registered, lab/out/sottoreticolo_prereg.json):
//! SOTTO_DIR=<dir> cargo test --bin atomo sottoreticolo_esame -- --ignored --nocapture --test-threads 1

use std::time::Instant;

use crate::cplx::C;
use crate::zz::{self, Atomo, Impulso, N};

fn sovrapp(a: &[[C; 4]], b: &[[C; 4]]) -> C {
    a.iter().zip(b.iter()).map(|(x, y)| { let mut s = C::ZERO; for q in 0..4 { s += x[q].conj() * y[q]; } s })
        .fold(C::ZERO, |u, v| u + v)
}
fn norma2(a: &[[C; 4]]) -> f64 { a.iter().map(|p| p.iter().map(|c| c.norm2()).sum::<f64>()).sum() }
fn normalizza(mut v: Vec<[C; 4]>) -> Vec<[C; 4]> {
    let k = 1.0 / norma2(&v).sqrt();
    for p in v.iter_mut() { for c in p.iter_mut() { *c = c.scale(k); } }
    v
}
pub const COLORI: [[u8; 2]; 4] = [[0, 0], [0, 1], [1, 0], [1, 1]];

/// Weight of a whole-cube wave on each colour.
pub fn pesi_colori(v: &[[C; 4]]) -> [f64; 4] {
    let mut w = [0.0; 4];
    for (s, p) in v.iter().enumerate() {
        let c = zz::colore_di(zz::sito(s, None));
        w[(2 * c[0] + c[1]) as usize] += p.iter().map(|c| c.norm2()).sum::<f64>();
    }
    w
}

/// Largest |difference|² between a one-colour wave (expanded) and a whole-cube wave, and the whole-cube weight off
/// that colour.
fn scarto(compatto: &[[C; 4]], pieno: &[[C; 4]], col: [u8; 2]) -> (f64, f64) {
    let e = zz::espandi_colore(compatto, col);
    let d: f64 = e.iter().zip(pieno.iter()).map(|(a, b)| (0..4).map(|q| (a[q] - b[q]).norm2()).sum::<f64>()).sum();
    let w = pesi_colori(pieno);
    let fuori: f64 = w.iter().sum::<f64>() - w[(2 * col[0] + col[1]) as usize];
    (d, fuori)
}

/// The one-colour rule IS the whole-cube rule restricted to a colour: the same start (radial 1s on colour c, zero
/// elsewhere) evolved 24 ticks with a pulse along z and a transverse live field, both ways; the waves agree to
/// rounding, nothing leaks to the other colours, and the observables agree. All 4 colours.
#[test]
fn sottoreticolo_equivalente() {
    let rad = Atomo::new();
    let w = pesi_colori(&rad.psi);
    eprintln!("radial 1s on the 4 colours: {w:?}");
    for col in COLORI {
        let phi = normalizza(zz::estrai_colore(&rad.psi, col));
        let mut a = Atomo::new_in(Some(col));
        a.imposta(phi.clone());
        let mut b = Atomo::new();
        b.imposta(zz::espandi_colore(&phi, col));
        for at in [&mut a, &mut b] {
            at.impulsi.push(Impulso { e0: 0.02, omega: 0.07, t0: 6.0, sigma: 4.0 });
            at.campo_est = [0.003, -0.002, 0.001];
        }
        for _ in 0..24 { a.avanti(); b.avanti(); }
        let (d, fuori) = scarto(&a.psi, &b.psi, col);
        let (pa, pb) = (a.posizione(), b.posizione());
        let dp = (0..3).map(|k| (pa[k] - pb[k]).abs()).fold(0.0, f64::max);
        let dd = (a.dipolo() - b.dipolo()).abs();
        let dg = (a.p_fondamentale() - b.p_fondamentale()).abs();
        let dn = (norma2(&a.psi) - 1.0).abs();
        // the Born sampling and the trilinear sample see the same wave
        let q = [0.3, -0.7, 1.1];
        let (sa, sb) = (a.campiona(q), b.campiona(q));
        eprintln!("colour {col:?}: |Δψ|² {d:.2e}, off-colour weight {fuori:.2e}, Δ⟨r⟩ {dp:.2e}, Δ⟨z⟩ {dd:.2e}, Δp1s {dg:.2e}, |norm−1| {dn:.2e}, ρ sample {:.3e} / {:.3e}",
                  sa.iter().map(|c| c.norm2()).sum::<f64>(), sb.iter().map(|c| c.norm2()).sum::<f64>());
        assert!(d < 1e-24 && fuori == 0.0 && dp < 1e-13 && dd < 1e-13 && dg < 1e-13 && dn < 1e-12);
        for _ in 0..24 { a.indietro(); b.indietro(); }
        let back: f64 = a.psi.iter().zip(phi.iter()).map(|(x, y)| (0..4).map(|q| (x[q] - y[q]).norm2()).sum::<f64>()).sum();
        assert!(back < 1e-24, "one colour backwards: {back:e}");
    }
    // the whole cube is the sum of its 4 colours, evolved separately
    let mut full = Atomo::new();
    full.impulsi.push(Impulso { e0: 0.02, omega: 0.07, t0: 6.0, sigma: 4.0 });
    let mut parti: Vec<Atomo> = COLORI.iter().map(|&c| {
        let mut a = Atomo::new_in(Some(c));
        a.imposta(zz::estrai_colore(&rad.psi, c)); // not renormalised: the colour's share
        a.impulsi = full.impulsi.clone();
        a
    }).collect();
    for _ in 0..12 { full.avanti(); for p in parti.iter_mut() { p.avanti(); } }
    let mut somma = vec![[C::ZERO; 4]; N * N * N];
    for (p, &c) in parti.iter().zip(COLORI.iter()) {
        for (x, y) in somma.iter_mut().zip(zz::espandi_colore(&p.psi, c)) { for q in 0..4 { x[q] += y[q]; } }
    }
    let d: f64 = somma.iter().zip(full.psi.iter()).map(|(a, b)| (0..4).map(|q| (a[q] - b[q]).norm2()).sum::<f64>()).sum();
    eprintln!("whole cube − Σ colours after 12 ticks: |Δψ|² = {d:.2e}");
    assert!(d < 1e-24);
}

/// Born sampling on one colour: a site stands for its Wigner–Seitz cell (4 cells of the cube: |δ|∞ ≤ a, |δ|₁ ≤ 3a/2,
/// ⟨δ²⟩ = 19/32 a²). The balls of the radial 1s on one colour give the wave's own ⟨r²⟩ (on the sites) plus that
/// smoothing, within 5 standard errors (200 000 balls); the whole cube the same with its cell (⟨δ²⟩ = a²/4).
#[test]
fn sottoreticolo_palline() {
    let rad = Atomo::new();
    let mut one = Atomo::new_in(Some(zz::COLORE_LAB));
    one.imposta(normalizza(zz::estrai_colore(&rad.psi, zz::COLORE_LAB)));
    let siti = |a: &Atomo| a.psi.iter().enumerate().map(|(s, p)| {
        let q = zz::sito(s, a.colore);
        let r2: f64 = q.iter().map(|&i| zz::x_of(i).powi(2)).sum();
        r2 * p.iter().map(|c| c.norm2()).sum::<f64>()
    }).sum::<f64>();
    let palle = |a: &Atomo, seme: u64| {
        let p = a.nascite(200_000, &mut crate::rng::Rng::new(seme));
        let n = p.len() as f64;
        let r2: Vec<f64> = p.iter().map(|x| x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).collect();
        let m = r2.iter().sum::<f64>() / n;
        (m, (r2.iter().map(|x| (x - m).powi(2)).sum::<f64>() / n).sqrt() / n.sqrt())
    };
    let a2 = zz::A * zz::A;
    for (at, k, seme) in [(&rad, 0.25 * a2, 7), (&one, 19.0 / 32.0 * a2, 8)] {
        let (m, sd) = palle(at, seme);
        let z = (m - k - siti(at)) / sd;
        eprintln!("colour {:?}: ⟨r²⟩ sites {:.4}, balls {:.4} − smoothing {:.4}: {:+.2} σ", at.colore, siti(at), m, k, z);
        assert!(z.abs() < 5.0);
    }
}

/// Development (declared, before the pre-registration): which colours the Laboratorio's states populate, whether the
/// 4 colours behave alike, and the cost of one tick.
#[test]
#[ignore]
fn sottoreticolo_sviluppo() {
    let rad = Atomo::new();
    eprintln!("radial 1s (Atomo::new) per colour: {:?}", pesi_colori(&rad.psi));
    let r2p = Atomo::con_stato(1, 2);
    eprintln!("radial 2p½ per colour: {:?}", pesi_colori(&r2p.psi));
    // timing, interleaved
    let mut full = Atomo::new();
    let mut one = Atomo::new_in(Some(zz::COLORE_LAB));
    let (mut tf, mut to) = (0.0, 0.0);
    for _ in 0..5 {
        let t0 = Instant::now(); for _ in 0..10 { full.avanti(); } tf += t0.elapsed().as_secs_f64();
        let t0 = Instant::now(); for _ in 0..10 { one.avanti(); } to += t0.elapsed().as_secs_f64();
    }
    eprintln!("avanti: whole cube {:.2} ms/tick, one colour {:.2} ms/tick, ×{:.2}", tf / 50.0 * 1e3, to / 50.0 * 1e3, tf / to);
    // the clean 1s of the Laboratorio (whole cube) and its colours
    let t0 = Instant::now();
    let (s1, altri) = rad.uno_s_esatto(1200);
    eprintln!("uno_s_esatto whole cube: ω {:.8} resid {:.1e} ({:?}); per colour {:?}", s1.omega, s1.resid, t0.elapsed(), pesi_colori(&s1.v));
    for l in &altri { eprintln!("  other level ω {:.8} resid {:.1e} colours {:?}", l.omega, l.resid, pesi_colori(&l.v)); }
    let t0 = Instant::now();
    let one = Atomo::new_in(Some(zz::COLORE_LAB));
    let (u1, altri1) = one.uno_s_esatto(1200);
    eprintln!("uno_s_esatto one colour: ω {:.8} resid {:.1e} ({:?})", u1.omega, u1.resid, t0.elapsed());
    for l in &altri1 { eprintln!("  other level ω {:.8} resid {:.1e}", l.omega, l.resid); }
    // each colour of the whole-cube 1s, renormalised: is it the same level? same response to a pulse?
    let w21 = zz::livello(1, 2).omega - zz::livello(-1, 1).omega;
    let mut serie: Vec<Vec<(f64, f64)>> = vec![];
    for col in COLORI {
        let v = zz::estrai_colore(&s1.v, col);
        if norma2(&v) < 1e-6 { eprintln!("colour {col:?}: empty"); serie.push(vec![]); continue; }
        let phi = normalizza(v);
        let mut a = Atomo::new_in(Some(col));
        a.imposta(phi.clone());
        let mut b = a.clone();
        b.avanti(); b.avanti();
        let lam = sovrapp(&phi, &b.psi);
        a.impulsi.push(Impulso { e0: 0.010, omega: w21, t0: 75.0, sigma: 30.0 });
        let mut s = vec![];
        for _ in 0..800 { a.avanti(); s.push((1.0 - a.p_fondamentale(), a.dipolo())); }
        eprintln!("colour {col:?}: ω {:.10} |λ| {:.12}, p_exc(800) {:.6e}, max|⟨z⟩| {:.5e}", -lam.arg(), lam.norm2().sqrt(), s[799].0,
                  s.iter().map(|x| x.1.abs()).fold(0.0, f64::max));
        serie.push(s);
    }
    for i in 0..4 { for j in i + 1..4 {
        if serie[i].is_empty() || serie[j].is_empty() { continue; }
        let dp = serie[i].iter().zip(&serie[j]).map(|(a, b)| (a.0 - b.0).abs()).fold(0.0, f64::max);
        let dd = serie[i].iter().zip(&serie[j]).map(|(a, b)| (a.1 - b.1).abs()).fold(0.0, f64::max);
        eprintln!("colours {i} vs {j}: max Δp_exc {dp:.2e}, max Δ⟨z⟩ {dd:.2e}");
    }}
}

/// Development (declared): cost of one light tick of the coupled box (128³, once every 8 atom ticks).
#[test]
#[ignore]
fn sottoreticolo_sviluppo_luce() {
    let mut l = crate::accoppiato::LuceAtomo::new();
    l.passo([0.0; 3]);
    let t0 = Instant::now();
    for _ in 0..10 { l.passo([0.0; 3]); }
    eprintln!("light tick {:.1} ms (per atom tick {:.2} ms)", t0.elapsed().as_secs_f64() / 10.0 * 1e3, t0.elapsed().as_secs_f64() / 80.0 * 1e3);
    let full = Atomo::new();
    let one = Atomo::new_in(Some(zz::COLORE_LAB));
    let t0 = Instant::now(); for _ in 0..20 { let _ = full.posizione(); let _ = full.proiezione(&full.psi0); } let tf = t0.elapsed().as_secs_f64() / 20.0;
    let t0 = Instant::now(); for _ in 0..20 { let _ = one.posizione(); let _ = one.proiezione(&one.psi0); } let to = t0.elapsed().as_secs_f64() / 20.0;
    eprintln!("posizione+proiezione: {:.2} / {:.2} ms", tf * 1e3, to * 1e3);
}

fn omega_u2(a: &Atomo, phi: &[[C; 4]]) -> f64 {
    let mut b = a.clone();
    b.impulsi.clear();
    b.campo_est = [0.0; 3];
    b.imposta(phi.to_vec());
    b.avanti();
    b.avanti();
    -sovrapp(phi, &b.psi).arg()
}

/// The series of one coupled run: per tick ⟨z⟩; per light tick p_exc, p2p, ⟨x⟩, ⟨y⟩, probe, light energy.
#[derive(Default)]
struct Serie { d: Vec<f64>, p: Vec<f64>, p2: Vec<f64>, x: Vec<f64>, y: Vec<f64>, sonda: Vec<f64>, en: Vec<f64>, norma: f64, secondi: f64 }

fn passo_misurato(x: &mut crate::accoppiato::Esperimento, s: &mut Serie) {
    let t0 = Instant::now();
    x.passo();
    s.secondi += t0.elapsed().as_secs_f64();
    s.d.push(x.serie_d.last().unwrap().1);
    if x.at.ticks % crate::accoppiato::SUB as u64 == 0 {
        s.p.push(x.serie_p.last().unwrap().1);
        s.p2.push(x.p2p());
        let r = x.at.posizione();
        s.x.push(r[0]);
        s.y.push(r[1]);
        s.sonda.push(x.luce.sonda.last().unwrap().1);
        s.en.push(x.luce.serie_e.last().unwrap().1);
    }
}

fn rel(a: &[f64], b: &[f64]) -> f64 {
    let m = b.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max) / m
}
fn abs_max(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max) }

fn num(v: f64) -> String { if v.is_finite() { format!("{v:e}") } else { "null".into() } }
fn arr(v: &[f64]) -> String { format!("[{}]", v.iter().map(|x| num(*x)).collect::<Vec<_>>().join(", ")) }
fn obj(kv: &[(&str, String)]) -> String { format!("{{{}}}", kv.iter().map(|(k, v)| format!("\"{k}\": {v}")).collect::<Vec<_>>().join(", ")) }
fn maxabs(v: &[f64]) -> f64 { v.iter().fold(0.0f64, |a, b| a.max(b.abs())) }

/// THE EXAM (lab/out/sottoreticolo_prereg.json). Writes sottoreticolo_results.json in SOTTO_DIR (default lab/out).
#[test]
#[ignore]
fn sottoreticolo_esame() {
    use crate::accoppiato::Esperimento;
    let dir = std::env::var("SOTTO_DIR").unwrap_or("lab/out".into());
    let ticks: usize = std::env::var("SOTTO_TICKS").ok().and_then(|v| v.parse().ok()).unwrap_or(3200);
    let col = zz::COLORE_LAB;
    let e0 = crate::zzscena::E0_ACC;
    let mut out: Vec<(&str, String)> = vec![];
    let mut esiti: Vec<(&str, bool)> = vec![];
    // preparation (S3)
    let t0 = Instant::now();
    let mut xf = Esperimento::nuovo(600, None);
    let t_prep_f = t0.elapsed().as_secs_f64();
    let t0 = Instant::now();
    let mut xo = Esperimento::nuovo(600, Some(col));
    let t_prep_o = t0.elapsed().as_secs_f64();
    eprintln!("prepared: whole {t_prep_f:.1} s, one colour {t_prep_o:.1} s; ω1s {:.10} / {:.10}, ω21 {:.8} / {:.8}", xf.omega1s, xo.omega1s, xf.omega_21(), xo.omega_21());
    out.push(("S3_prepara_s", obj(&[("intero", num(t_prep_f)), ("colore", num(t_prep_o)), ("rapporto", num(t_prep_f / t_prep_o))])));
    out.push(("pesi_colori", obj(&[("ordine", "\"[00, 01, 10, 11] = [(i+j)%2, (j+k)%2]\"".into()),
        ("radiale_1s", arr(&pesi_colori(&Atomo::new().psi))), ("radiale_2p", arr(&pesi_colori(&Atomo::con_stato(1, 2).psi))),
        ("uno_s_lab_intero", arr(&pesi_colori(&xf.phi1s))), ("due_p_lab_intero", arr(&pesi_colori(&xf.phi2p)))])));
    esiti.push(("S3", t_prep_f / t_prep_o >= 2.5));

    // D1, D2: density
    let momenti = |v: &[[C; 4]], c: zz::Colore| {
        let (mut m1, mut m2, mut w) = (0.0, 0.0, 0.0);
        for (s, p) in v.iter().enumerate() {
            let q = zz::sito(s, c);
            let x = [zz::x_of(q[0]), zz::x_of(q[1]), zz::x_of(q[2])];
            let r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
            let rho: f64 = p.iter().map(|c| c.norm2()).sum();
            m1 += rho * r2.sqrt(); m2 += rho * r2; w += rho;
        }
        (m1 / w, m2 / w)
    };
    let (sf, so) = (momenti(&xf.phi1s, None), momenti(&xo.phi1s, Some(col)));
    let d1 = ((so.0 / sf.0 - 1.0).abs()).max((so.1 / sf.1 - 1.0).abs());
    let palle = |a: &Atomo, seme: u64| {
        let p = a.nascite(400_000, &mut crate::rng::Rng::new(seme));
        let n = p.len() as f64;
        let r2: Vec<f64> = p.iter().map(|x| x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).collect();
        let m2 = r2.iter().sum::<f64>() / n;
        let m1 = r2.iter().map(|x| x.sqrt()).sum::<f64>() / n;
        let sd = (r2.iter().map(|x| (x - m2).powi(2)).sum::<f64>() / n).sqrt() / n.sqrt();
        (m1, m2, sd)
    };
    let a2 = zz::A * zz::A;
    let (bf, bo) = (palle(&xf.at, 11), palle(&xo.at, 12));
    let zf = (bf.1 - 0.25 * a2 - sf.1) / bf.2;
    let zo = (bo.1 - 19.0 / 32.0 * a2 - so.1) / bo.2;
    eprintln!("D1 sites ⟨r⟩ {:.5}/{:.5} ⟨r²⟩ {:.5}/{:.5} (max rel {d1:.2e}); D2 balls ⟨r²⟩−kernel z-scores {zf:.2}/{zo:.2}, ⟨r⟩ {:.4}/{:.4}", sf.0, so.0, sf.1, so.1, bf.0, bo.0);
    out.push(("densita", obj(&[("siti_r_intero", num(sf.0)), ("siti_r_colore", num(so.0)), ("siti_r2_intero", num(sf.1)), ("siti_r2_colore", num(so.1)), ("D1_rel", num(d1)),
        ("palline_r_intero", num(bf.0)), ("palline_r_colore", num(bo.0)), ("palline_r2_intero", num(bf.1)), ("palline_r2_colore", num(bo.1)),
        ("kernel_d2_intero", num(0.25 * a2)), ("kernel_d2_colore", num(19.0 / 32.0 * a2)), ("z_intero", num(zf)), ("z_colore", num(zo))])));
    esiti.push(("D1", d1 <= 5e-3));
    esiti.push(("D2", zf.abs() <= 5.0 && zo.abs() <= 5.0 && (bo.0 / bf.0 - 1.0).abs() <= 0.01));

    // S1: the rule alone
    let (mut af, mut ao) = (xf.at.clone(), xo.at.clone());
    let (mut tf, mut to) = (vec![], vec![]);
    for _ in 0..5 {
        let t = Instant::now(); for _ in 0..20 { af.avanti(); } tf.push(t.elapsed().as_secs_f64() / 20.0 * 1e3);
        let t = Instant::now(); for _ in 0..20 { ao.avanti(); } to.push(t.elapsed().as_secs_f64() / 20.0 * 1e3);
    }
    let med = |v: &[f64]| { let mut v = v.to_vec(); v.sort_by(|a, b| a.total_cmp(b)); v[v.len() / 2] };
    let (mf, mo) = (med(&tf), med(&to));
    eprintln!("S1 avanti: {mf:.2} / {mo:.2} ms → ×{:.2}", mf / mo);
    out.push(("S1_avanti", obj(&[("ms_intero", num(mf)), ("ms_colore", num(mo)), ("rapporto", num(mf / mo)), ("blocchi_intero_ms", arr(&tf)), ("blocchi_colore_ms", arr(&to))])));
    esiti.push(("S1", mf / mo >= 3.0));
    drop(af); drop(ao);

    // M1
    let byte = |n: usize| (n * (3 * 64 + 16 + 8 + 2 * 64)) as f64 / 1e6;
    out.push(("M1_memoria", obj(&[("MB_intero", num(byte(zz::n_siti(None)))), ("MB_colore", num(byte(zz::n_siti(Some(col))))),
        ("rapporto", num(zz::n_siti(None) as f64 / zz::n_siti(Some(col)) as f64)),
        ("filtro_uno_s_MB_intero", num(25.0 * 64.0 * zz::n_siti(None) as f64 / 1e6)), ("filtro_uno_s_MB_colore", num(25.0 * 64.0 * zz::n_siti(Some(col)) as f64 / 1e6)),
        ("luce_MB", num((crate::accoppiato::N_LUCE.pow(3) * 3 * 8) as f64 / 1e6))])));
    esiti.push(("M1", zz::n_siti(None) == 4 * zz::n_siti(Some(col))));

    // E2: the Laboratorio as it was vs one colour
    let w21f = xf.omega_21();
    let w21o = xo.omega_21();
    xf.ricomincia(); xo.ricomincia();
    xf.spara(w21f, e0); xo.spara(w21o, e0);
    let (mut sf2, mut so2) = (Serie::default(), Serie::default());
    for k in 0..ticks {
        passo_misurato(&mut xf, &mut sf2);
        passo_misurato(&mut xo, &mut so2);
        if k % 800 == 799 { eprintln!("E2 tick {}: p_exc {:.5} / {:.5}", k + 1, sf2.p.last().unwrap(), so2.p.last().unwrap()); }
    }
    let e2 = [("d_omega1s", (xf.omega1s - xo.omega1s).abs()), ("d_omega21_rel", (w21f - w21o).abs() / w21f),
        ("omega1s_intero", xf.omega1s), ("omega1s_colore", xo.omega1s), ("omega21_intero", w21f), ("omega21_colore", w21o),
        ("dipolo_rel", rel(&so2.d, &sf2.d)), ("sonda_rel", rel(&so2.sonda, &sf2.sonda)), ("energia_luce_rel", rel(&so2.en, &sf2.en)),
        ("p_exc_abs", abs_max(&so2.p, &sf2.p)), ("p_exc_max_intero", maxabs(&sf2.p)), ("p_exc_max_colore", maxabs(&so2.p)),
        ("p2p_abs", abs_max(&so2.p2, &sf2.p2)), ("p_exc_fine_intero", *sf2.p.last().unwrap()), ("p_exc_fine_colore", *so2.p.last().unwrap()),
        ("dipolo_max", maxabs(&sf2.d)), ("sonda_max", maxabs(&sf2.sonda))];
    eprintln!("E2: {e2:?}");
    let g = |v: &[(&str, f64)], k: &str| v.iter().find(|x| x.0 == k).unwrap().1;
    let e2ok = g(&e2, "d_omega1s") <= 1e-8 && g(&e2, "d_omega21_rel") <= 1e-5 && g(&e2, "dipolo_rel") <= 1e-6 && g(&e2, "sonda_rel") <= 1e-6
        && g(&e2, "energia_luce_rel") <= 1e-6 && g(&e2, "p_exc_abs") <= 5e-3;
    out.push(("E2", obj(&e2.iter().map(|(k, v)| (*k, num(*v))).collect::<Vec<_>>())));
    out.push(("serie_E2_per_tick_luce", obj(&[("p_exc_intero", arr(&sf2.p)), ("p_exc_colore", arr(&so2.p)), ("sonda_intero", arr(&sf2.sonda)), ("sonda_colore", arr(&so2.sonda))])));
    esiti.push(("E2", e2ok));

    // E1: the same one-colour start on both storages
    let phi = xo.phi1s.clone();
    xf.phi1s = zz::espandi_colore(&phi, col);
    xf.phi2p = zz::espandi_colore(&xo.phi2p, col);
    xf.omega1s = xo.omega1s;
    xf.omega2p = xo.omega2p;
    xf.ricomincia(); xo.ricomincia();
    let w21 = xo.omega_21();
    xf.spara(w21, e0); xo.spara(w21, e0);
    let (mut rf, mut ro) = (Serie::default(), Serie::default());
    for k in 0..ticks {
        passo_misurato(&mut xf, &mut rf);
        passo_misurato(&mut xo, &mut ro);
        if k % 800 == 799 { eprintln!("E1 tick {}: p_exc {:.6e} / {:.6e}", k + 1, rf.p.last().unwrap(), ro.p.last().unwrap()); }
    }
    rf.norma = norma2(&xf.at.psi);
    ro.norma = norma2(&xo.at.psi);
    let w = pesi_colori(&xf.at.psi);
    let fuori = w.iter().sum::<f64>() - w[(2 * col[0] + col[1]) as usize];
    let (of, oo) = (omega_u2(&xf.at, &xf.phi1s), omega_u2(&xo.at, &xo.phi1s));
    let e1 = [("p_exc_rel", rel(&ro.p, &rf.p)), ("p2p_rel", rel(&ro.p2, &rf.p2)), ("dipolo_rel", rel(&ro.d, &rf.d)), ("x_rel", rel(&ro.x, &rf.x)), ("y_rel", rel(&ro.y, &rf.y)),
        ("sonda_rel", rel(&ro.sonda, &rf.sonda)), ("energia_luce_rel", rel(&ro.en, &rf.en)), ("norma_intero", rf.norma), ("norma_colore", ro.norma), ("d_norma", (rf.norma - ro.norma).abs()),
        ("omega1s_intero", of), ("omega1s_colore", oo), ("d_omega1s", (of - oo).abs()), ("peso_fuori_colore", fuori),
        ("p_exc_max", maxabs(&rf.p)), ("dipolo_max", maxabs(&rf.d)), ("sonda_max", maxabs(&rf.sonda))];
    eprintln!("E1: {e1:?}");
    let e1ok = ["p_exc_rel", "p2p_rel", "dipolo_rel", "x_rel", "y_rel", "sonda_rel", "energia_luce_rel", "d_norma", "d_omega1s"].iter().all(|k| g(&e1, k) <= 1e-12) && fuori <= 1e-30;
    out.push(("E1", obj(&e1.iter().map(|(k, v)| (*k, num(*v))).collect::<Vec<_>>())));
    esiti.push(("E1", e1ok));
    // S2: the coupled step, timed over the E1 runs (and E2)
    let s2 = rf.secondi / ro.secondi;
    let ms = |x: f64| num(x / ticks as f64 * 1e3);
    out.push(("S2_passo_accoppiato", obj(&[("ms_intero", ms(rf.secondi)), ("ms_colore", ms(ro.secondi)), ("rapporto", num(s2)), ("E2_ms_intero", ms(sf2.secondi)), ("E2_ms_colore", ms(so2.secondi))])));
    esiti.push(("S2", s2 >= 1.5));
    eprintln!("S2 passo: {} / {} ms → ×{s2:.2}", ms(rf.secondi), ms(ro.secondi));
    out.push(("ticks", ticks.to_string()));
    out.push(("rayon_threads", rayon::current_num_threads().to_string()));
    let es = obj(&esiti.iter().map(|(k, v)| (*k, if *v { "\"PASS\"".to_string() } else { "\"FAIL\"".to_string() })).collect::<Vec<_>>());
    eprintln!("ESITI: {es}");
    out.insert(0, ("esiti", es));
    std::fs::write(format!("{dir}/sottoreticolo_results.json"), obj(&out)).unwrap();
}

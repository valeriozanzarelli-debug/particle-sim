//! EXACT SCATTERED LIGHT (docs/ATOMO.md §3.8, pre-registration lab/out/diffusa_prereg.json).
//! The light rule, the sponge and the deposit are linear: evolving F_free (the pulse, no source) and F_sc (only the
//! atom's current) and letting the atom feel F_free + F_sc must give the single-box run (X1, X2). Then F_sc is the light
//! the atom scatters and emits, with no subtraction noise. C0/C1: a π/2 pulse computed from the rule's own dipoles.

use crate::accoppiato::{self as acc, Esperimento};
use rayon::prelude::*;

fn norme(a: &acc::LuceAtomo, b: &acc::LuceAtomo) -> (f64, f64, f64) {
    // ||F_A − (F_free + F_sc)||², ||F_A||², ||F_sc||²
    let sc = b.sc.as_ref().unwrap();
    (0..a.f.f.len()).into_par_iter().map(|s| {
        let (p, q, r) = (a.f.f[s], b.f.f[s], sc.f[s]);
        let mut d = 0.0;
        let mut na = 0.0;
        let mut ns = 0.0;
        for c in 0..3 {
            let (dr, di) = (p[c].re as f64 - q[c].re as f64 - r[c].re as f64, p[c].im as f64 - q[c].im as f64 - r[c].im as f64);
            d += dr * dr + di * di;
            na += p[c].norm2() as f64;
            ns += r[c].norm2() as f64;
        }
        (d, na, ns)
    }).reduce(|| (0.0, 0.0, 0.0), |x, y| (x.0 + y.0, x.1 + y.1, x.2 + y.2))
}

/// max |E| (real part) of a box, all cells or only farther than 8 light cells from the atom.
pub fn emax(f: &[[crate::fotone3d::Cf; 3]], lontano: bool) -> f64 {
    let n = acc::N_LUCE;
    let h = (n / 2) as i64;
    f.par_iter().enumerate().map(|(s, p)| {
        if lontano {
            let d = [(s / (n * n)) as i64 - h, ((s / n) % n) as i64 - h, (s % n) as i64 - h];
            if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= 64 { return 0.0; }
        }
        ((p[0].re * p[0].re + p[1].re * p[1].re + p[2].re * p[2].re) as f64).sqrt()
    }).reduce(|| 0.0, f64::max)
}

const W32: f64 = 0.079365;
const W12: f64 = 0.076614;
const D32: f64 = 0.97868;
const D12: f64 = 0.6029;

struct Corsa { pmax: f64, s32: f64, s12: f64, a_amp: f64, dip: f64, e_libera_max: f64, e_sc_max: f64 }

/// Run a split experiment for `ticks` atom ticks; optional single-box twin `a` in lockstep (exam X1/X2).
fn corri(b: &mut Esperimento, mut a: Option<&mut Esperimento>, ticks: usize, rep: &mut serde_like::Mappa) -> Corsa {
    let e = acc::e_eff();
    let (mut s32, mut s12) = ((0.0, 0.0), (0.0, 0.0));
    let mut e_libera_max: f64 = emax(&b.luce.f.f, false);
    let mut e_sc_max: f64 = 0.0;
    let (mut x1, mut x1_sc, mut x2): (f64, f64, f64) = (0.0, 0.0, 0.0);
    let (mut ms_a, mut ms_b) = (0.0, 0.0);
    let mut nl = 0;
    for _ in 0..ticks {
        let t0 = std::time::Instant::now();
        b.passo();
        ms_b += t0.elapsed().as_secs_f64() * 1e3;
        if let Some(a) = a.as_deref_mut() {
            let t0 = std::time::Instant::now();
            a.passo();
            ms_a += t0.elapsed().as_secs_f64() * 1e3;
        }
        let t = b.at.t();
        let ez = e * b.luce.campo_libero_atomo()[2] * crate::zz::A;
        s32.0 += ez * (W32 * t).cos(); s32.1 += ez * (W32 * t).sin();
        s12.0 += ez * (W12 * t).cos(); s12.1 += ez * (W12 * t).sin();
        if b.serie_p.len() > nl {
            nl = b.serie_p.len();
            if let Some(a) = a.as_deref() {
                let (d, na, ns) = norme(&a.luce, &b.luce);
                x1 = x1.max((d / na.max(1e-300)).sqrt());
                if ns > 0.0 { x1_sc = x1_sc.max((d / ns).sqrt()); }
                x2 = x2.max((a.serie_p.last().unwrap().1 - b.serie_p.last().unwrap().1).abs());
            }
            if nl % 4 == 0 {
                e_libera_max = e_libera_max.max(emax(&b.luce.f.f, false));
                e_sc_max = e_sc_max.max(emax(&b.luce.sc.as_ref().unwrap().f, true));
            }
            if nl % 25 == 0 { eprintln!("  t {:.0}  p_exc {:.4}  x1 {:.2e}  x2 {:.2e}", t, b.serie_p.last().unwrap().1, x1, x2); }
        }
    }
    if a.is_some() {
        rep.num("X1_max_rel_L2", x1);
        rep.num("info_err_su_diffusa", x1_sc);
        rep.num("X2_max_dp_exc", x2);
        rep.num("X3_ms_tick_A", ms_a / ticks as f64);
        rep.num("X3_ms_tick_B", ms_b / ticks as f64);
        rep.num("X3_rapporto", ms_b / ms_a.max(1e-9));
    }
    let pmax = b.serie_p.iter().map(|p| p.1).fold(0.0, f64::max);
    Corsa { pmax, s32: (s32.0 * s32.0 + s32.1 * s32.1).sqrt(), s12: (s12.0 * s12.0 + s12.1 * s12.1).sqrt(),
            a_amp: e_libera_max / e_sc_max.max(1e-300), dip: b.ampiezza_dipolo(200.0), e_libera_max, e_sc_max }
}

/// A tiny JSON writer (no serde in this crate's tests).
mod serde_like {
    pub struct Mappa(pub Vec<(String, String)>);
    impl Mappa {
        pub fn num(&mut self, k: &str, v: f64) { self.0.push((k.into(), format!("{v:.6e}"))); }
        pub fn txt(&mut self, k: &str, v: &str) { self.0.push((k.into(), format!("{v:?}"))); }
        pub fn json(&self) -> String {
            let r: Vec<String> = self.0.iter().map(|(k, v)| format!("  {k:?}: {v}")).collect();
            format!("{{\n{}\n}}\n", r.join(",\n"))
        }
    }
}

#[test]
#[ignore]
fn diffusa_esame() {
    let t0 = std::time::Instant::now();
    let x = Esperimento::nuovo(600, None);
    eprintln!("preparato in {:.0} s: omega_21 {:.5}", t0.elapsed().as_secs_f64(), x.omega_21());
    let nuovo = || Esperimento::da_stati(x.phi1s.clone(), x.omega1s, x.phi2p.clone(), x.omega2p, crate::zz::NUCLEO);
    let mut rep = serde_like::Mappa(vec![]);
    rep.num("omega_21", x.omega_21());
    let ticks = 1200;
    // X1/X2/X3 and C0
    let mut a = nuovo();
    let mut b = nuovo().con_diffusa();
    a.spara(x.omega_21(), crate::zzscena::E0_ACC);
    b.spara(x.omega_21(), crate::zzscena::E0_ACC);
    let r = corri(&mut b, Some(&mut a), ticks, &mut rep);
    drop(a);
    let meta = (D32 * D32 * r.s32 * r.s32 + D12 * D12 * r.s12 * r.s12).sqrt();
    let e0_pi2 = crate::zzscena::E0_ACC * (std::f64::consts::FRAC_PI_4) / meta;
    rep.num("C0_S32", r.s32);
    rep.num("C0_S12", r.s12);
    rep.num("C0_theta_mezzi_e0_0.004", meta);
    rep.num("C0_e0_pi2", e0_pi2);
    rep.num("info_p_previsto_e0_0.004", meta.sin().powi(2));
    rep.num("info_p_max_misurato_B", r.pmax);
    rep.num("info_A_B", r.a_amp);
    rep.num("info_Elibera_max_B", r.e_libera_max);
    rep.num("info_Esc_max_lontano_B", r.e_sc_max);
    rep.num("info_dipolo_B", r.dip);
    eprintln!("B: p_max {:.4}  previsto {:.4}  e0_pi2 {:.5}  A {:.1}", r.pmax, meta.sin().powi(2), e0_pi2, r.a_amp);
    let x1 = rep.0.iter().find(|p| p.0 == "X1_max_rel_L2").unwrap().1.parse::<f64>().unwrap();
    let x2 = rep.0.iter().find(|p| p.0 == "X2_max_dp_exc").unwrap().1.parse::<f64>().unwrap();
    rep.txt("X1", if x1 <= 1e-6 { "PASS" } else { "FAIL" });
    rep.txt("X2", if x2 <= 1e-6 { "PASS" } else { "FAIL" });
    drop(b);
    // C1: the π/2 pulse
    let mut c = nuovo().con_diffusa();
    c.spara(x.omega_21(), e0_pi2);
    let rc = corri(&mut c, None, ticks, &mut rep);
    rep.num("C1_p_max", rc.pmax);
    rep.txt("C1", if (0.35..=0.65).contains(&rc.pmax) { "PASS" } else { "FAIL" });
    rep.num("info_A_C", rc.a_amp);
    rep.num("info_dipolo_C", rc.dip);
    rep.num("info_Esc_max_lontano_C", rc.e_sc_max);
    rep.num("secondi", t0.elapsed().as_secs_f64());
    eprintln!("C: p_max {:.4}  A {:.1}  dipolo {:.3}", rc.pmax, rc.a_amp, rc.dip);
    let out = std::env::var("DIFFUSA_OUT").unwrap_or("lab/out/diffusa_results.json".into());
    std::fs::write(&out, rep.json()).unwrap();
    eprintln!("{}", rep.json());
}

/// Cheap check (runs by default): with no atom current, splitting does not change the free light, and a source
/// deposited in the split box equals the same source in the single box (to f32 rounding).
#[test]
fn diffusa_linearita_luce() {
    let mut a = acc::LuceAtomo::new();
    let mut b = acc::LuceAtomo::new();
    b.con_diffusa();
    a.impulso(0.075, 0.004);
    b.impulso(0.075, 0.004);
    for k in 0..6 {
        let j = [0.0, 0.0, 1e-3 * (k as f64 * 0.3).sin()];
        a.passo(j);
        b.passo(j);
    }
    let (d, na, _) = norme(&a, &b);
    assert!((d / na).sqrt() < 1e-6, "split light differs: {:.2e}", (d / na).sqrt());
    assert!((a.energia() / b.energia() - 1.0).abs() < 1e-5);
}

/// The scattered view's data (no window): after the π/2 pulse the flow balls ride the energy current of the scattered
/// light outwards (mean radial velocity > 0.5 c), and glow balls exist. Reported, not a pre-registered exam.
#[test]
#[ignore]
fn diffusa_vista() {
    let mut x = Esperimento::nuovo(600, None).con_diffusa();
    x.spara(x.omega_21(), crate::zzscena::E0_PI2);
    let (vr, nf, ng, a) = crate::zzscena::prova_vista(x, 800);
    eprintln!("vista diffusa: v_r medio {:.3} c, palline di flusso {}, palline di luce diffusa {}, A = {:.1}", vr, nf, ng, a);
    assert!(vr > 0.5 && nf > 100 && ng > 100);
}

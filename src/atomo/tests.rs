//! Tests of the carbon viewer. `cargo test --bin atomo` (fast ones);
//! `cargo test --bin atomo -- --ignored` for the long energy validation (several minutes).

use crate::carbonio::{confini, istogramma_radiale, misura_dmc, Guida, Stato};
use crate::mondo::{Input, Livelli, Mondo, Z};
use crate::risultati::Risultati;
use crate::skyrme::Profilo;
use crate::tempo::Storia;
use crate::ui::E_ESATTA_C;

/// Pre-registered: the fixed-node DMC energy of carbon within 1 % of the exact non-relativistic −37.8450 Eₕ
/// (short run: a few seconds).
#[test]
fn energia_carbonio_entro_1_percento() {
    let r = Risultati::carica();
    let g = Guida { z: Z, stato: Stato::Fondamentale, p: r.param("c") };
    let (e, err, d) = misura_dmc(g, 0.005, 600, 1200, 2000, 5);
    let dev = (e - E_ESATTA_C) / E_ESATTA_C.abs();
    eprintln!("E(C) = {e:.4} ± {err:.4}  ({:+.3} %)", 100.0 * dev);
    assert!(dev.abs() < 0.01, "E = {e}");
    // the balls find the K/L shell boundary by themselves, between 0.3 and 1.0 bohr
    let h = istogramma_radiale(&d, crate::DR_ISTO, crate::R_ISTO);
    let c = confini(&h, crate::DR_ISTO);
    eprintln!("shell boundaries {c:?}");
    assert!(!c.is_empty() && c[0] > 0.3 && c[0] < 1.0, "{c:?}");
}

/// The long version (as in lab/atomo_misura.sh): 2000 balls, 40 000 ticks, total energy within 1 %.
#[test]
#[ignore]
fn energia_carbonio_lunga() {
    let r = Risultati::carica();
    let g = Guida { z: Z, stato: Stato::Fondamentale, p: r.param("c") };
    let (e, err, _) = misura_dmc(g, 0.005, 2000, 4000, 40_000, 100);
    eprintln!("E(C) = {e:.5} ± {err:.5}");
    assert!(((e - E_ESATTA_C) / E_ESATTA_C).abs() < 0.01);
}

#[test]
fn nodo_b1() {
    let p = Profilo::calcola();
    assert!((p.numero_barionico() - 1.0).abs() < 1e-3);
    assert!((p.energia_12pi2() - 1.2315).abs() < 2e-3, "{}", p.energia_12pi2());
    assert!((p.raggio_rms_fm() - 0.59).abs() < 0.02);
}

/// Going back and playing forward again rebuilds exactly the same balls (reversible time).
#[test]
fn tempo_reversibile_esatto() {
    let r = Risultati::carica();
    let lv = Livelli::da(&r);
    let mut m = Mondo::new(&r, 3);
    let mut s = Storia::new(&m);
    let mut inputs = vec![];
    for f in 0..30 {
        let inp = Input { dt_s: 1e-18, dt_schermo: 0.016, ticks: 1, t_kelvin: 300.0,
                          spara: if f == 3 { Some((lv.lambda_ecc_nm(), [1.0, 0.0, 0.0])) } else { None } };
        m.passo(&inp, &r, &lv);
        s.registra(&m, inp);
        inputs.push(inp);
    }
    let finale = m.dmc.x.clone();
    let nucleo = m.p.nucleo;
    let mut m2 = s.riprendi_da(14, &r, &lv);
    for inp in &inputs[14..] { m2.passo(inp, &r, &lv); }
    assert_eq!(m2.p.nucleo, nucleo);
    assert_eq!(m2.dmc.x, finale);
}

/// The photon rule: at the computed 2p→3s energy it is absorbed; in the visible it passes; in the far UV
/// it ionizes.
#[test]
fn fotoni() {
    let r = Risultati::carica();
    let lv = Livelli::da(&r);
    for (lam, atteso) in [(lv.lambda_ecc_nm(), Stato::Eccitato), (500.0, Stato::Fondamentale), (60.0, Stato::Ione)] {
        let mut m = Mondo::new(&r, 4);
        m.passo(&Input { dt_s: 0.0, dt_schermo: 0.0, ticks: 0, t_kelvin: 0.0, spara: Some((lam, [1.0, 0.0, 0.0])) }, &r, &lv);
        // 45 bohr at c
        let dt = 45.0 / crate::mondo::c_bohr_s();
        m.passo(&Input { dt_s: dt, dt_schermo: 0.0, ticks: 1, t_kelvin: 0.0, spara: None }, &r, &lv);
        assert_eq!(m.p.stato, atteso, "lambda {lam}");
        assert_eq!(m.dmc.g.ne(), atteso.ne());
    }
}

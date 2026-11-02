//! Esami dello studio "bagliori" (docs/BAGLIORI.md, pre-registrazione lab/out/gl_prereg.json).
//!
//! S1–S3: il vortice al bordo del fascio, misurato dalle palline di energia DISEGNATE (le stesse posizioni che
//! vanno sullo schermo, fotdisegno::Traccianti), per luce circolare +, circolare − e lineare.
//! E1–E3: il tremito (Zitterbewegung) dell'elettrone libero, misurato dai BAGLIORI disegnati (elscena.rs).
//!
//! cargo test --release --bin atomo bagliori -- --ignored --nocapture   (scrive lab/out/gl_*.json)

use crate::elscena::{self, centro_onda, centro_palline, fit_tremito, nascite, pacchetto, spinore_f3};
use crate::fotdisegno::Traccianti;
use crate::fotone3d::{base_trasversa, Fotone3D, Polar};
use crate::fotscena::{atteso, Parametri, PALLINE_QUAL};
use crate::rng::Rng;

const OUT: &str = "lab/out";

/// Velocità angolare media (rad/tick) delle palline disegnate nella corona 0,7–1,3 w0 attorno all'asse del fascio,
/// tra il tick t_fuoco − 4 e t_fuoco + 4 (attorno al fuoco), dallo spostamento vero di ogni pallina in un tick.
/// Restituisce (Ω, palline usate, Ω atteso = 2σ c/(k w0²)).
fn vortice(polar: Polar) -> (f64, usize, f64) {
    let par = Parametri { qual: 0, polar, ..Default::default() };
    let l = par.lancio();
    let mut f = Fotone3D::new(par.n());
    f.lancia(&l);
    let mut tr = Traccianti::new(&f, 1000 + 1, PALLINE_QUAL[0]);
    let (e1, e2) = base_trasversa(l.dir);
    let t_f = l.t_fuoco.round() as u64;
    let (mut somma, mut quante) = (0.0f64, 0usize);
    for t in 1..=t_f + 4 {
        let prima = tr.pos.clone();
        let eta_prima = tr.eta.clone();
        f.tick();
        tr.avanza(&f, 1.0, atteso(&l, f.ticks));
        if t + 4 < t_f { continue; }
        let c = f.misura(atteso(&l, f.ticks)).centro;
        for i in 0..tr.pos.len() {
            if tr.eta[i] <= eta_prima[i] { continue; } // rinata in questo tick
            let (p, q) = (prima[i], tr.pos[i]);
            let r = |x: [f32; 3]| { let d = [x[0] as f64 - c[0], x[1] as f64 - c[1], x[2] as f64 - c[2]];
                [d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2], d[0] * e2[0] + d[1] * e2[1] + d[2] * e2[2]] };
            let (a, b) = (r(p), r(q));
            let ra = (a[0] * a[0] + a[1] * a[1]).sqrt();
            if ra < 0.7 * l.w0 || ra > 1.3 * l.w0 { continue; }
            let mut dphi = b[1].atan2(b[0]) - a[1].atan2(a[0]);
            if dphi > std::f64::consts::PI { dphi -= std::f64::consts::TAU; }
            if dphi < -std::f64::consts::PI { dphi += std::f64::consts::TAU; }
            somma += dphi;
            quante += 1;
        }
    }
    let k = std::f64::consts::TAU / l.lambda;
    let segno = match polar { Polar::CircolarePiu => 1.0, Polar::CircolareMeno => -1.0, _ => 0.0 };
    (somma / quante.max(1) as f64, quante, segno * 2.0 / (k * l.w0 * l.w0))
}

#[test]
#[ignore]
fn bagliori_vortice_fotone() {
    let mut righe = vec![];
    let mut om = vec![];
    for (p, nome) in [(Polar::CircolarePiu, "circ+"), (Polar::CircolareMeno, "circ-"), (Polar::Orizzontale, "lineare")] {
        let (w, n, th) = vortice(p);
        eprintln!("{nome}: Ω bordo = {w:+.5} rad/tick da {n} spostamenti; atteso (fascio gaussiano, spin) {th:+.5}");
        righe.push(format!("  \"{nome}\": {{\"omega_rad_tick\": {w:e}, \"spostamenti\": {n}, \"omega_atteso\": {th:e}}}"));
        om.push((w, th));
    }
    let s1 = om[0].0 > 0.0 && om[1].0 < 0.0 && om[0].0 * om[1].0 < 0.0;
    let s2 = om[2].0.abs() < 0.15 * 0.5 * (om[0].0.abs() + om[1].0.abs());
    let s3 = [0, 1].iter().all(|&i| { let r = om[i].0 / om[i].1; (0.6..=1.4).contains(&r) });
    eprintln!("S1 segni opposti: {s1}; S2 lineare ~0: {s2}; S3 Ω/atteso in [0,6; 1,4]: {s3} ({:.3}, {:.3})", om[0].0 / om[0].1, om[1].0 / om[1].1);
    let json = format!("{{\n \"studio\": \"bagliori\",\n \"esami\": \"S1 S2 S3\",\n{},\n \"S1\": {s1},\n \"S2\": {s2},\n \"S3\": {s3},\n \"rapporti_S3\": [{}, {}]\n}}\n",
                       righe.join(",\n"), om[0].0 / om[0].1, om[1].0 / om[1].1);
    std::fs::write(format!("{OUT}/gl_vortice.json"), json).unwrap();
}

/// Il tremito dai bagliori: la stessa corsa di F3 (zz_free.py: 96³, θ = 0,25, σ = 3, t ≤ 4π) con 40 000 bagliori.
fn tremito(nome: &str, lettore: u8) -> String {
    let mut d = pacchetto(elscena::N_EL, elscena::THETA, elscena::SIGMA, spinore_f3());
    let mut rng = Rng::new(2027);
    let mut s = nascite(&d, elscena::N_BAGLIORI, &mut rng);
    let nt = (4.0 * std::f64::consts::PI / elscena::THETA).round() as usize;
    let (mut ts, mut bx, mut bz, mut ox, mut oz) = (vec![], vec![], vec![], vec![], vec![]);
    let mut dmax = 0.0f64;
    let mut ymax = 0.0f64;
    for t in 0..=nt {
        let (b, o) = (centro_palline(&s), centro_onda(&d));
        ts.push(t as f64 * d.a);
        bx.push(b[0]); bz.push(b[2]); ox.push(o[0]); oz.push(o[2]);
        let dd = ((b[0] - o[0]).powi(2) + (b[1] - o[1]).powi(2) + (b[2] - o[2]).powi(2)).sqrt();
        dmax = dmax.max(dd);
        ymax = ymax.max(b[1].abs());
        if t < nt {
            match lettore {
                0 => elscena::passo(&mut d, &mut s),
                1 => elscena::passo_flusso(&mut d, &mut s),
                _ => {
                    // diagnosi: Eulero in avanti con lo spinore interpolato al tick t
                    let v: Vec<[f64; 3]> = s.iter().map(|p| crate::walk::surfer_velocity(&d.sample(*p))).collect();
                    for (p, v) in s.iter_mut().zip(&v) { for c in 0..3 { p[c] += d.a * v[c]; } }
                    d.tick();
                }
            }
        }
    }
    let (abx, wbx) = fit_tremito(&ts, &bx, 1.5, 2.6);
    let (abz, wbz) = fit_tremito(&ts, &bz, 1.5, 2.6);
    let (aox, wox) = fit_tremito(&ts, &ox, 1.5, 2.6);
    let (aoz, woz) = fit_tremito(&ts, &oz, 1.5, 2.6);
    eprintln!("{nome}: bagliori x A {abx:.4} ω {wbx:.4}, z A {abz:.4} ω {wbz:.4}; onda x A {aox:.4} ω {wox:.4}, z A {aoz:.4} ω {woz:.4}; |bagliori − onda| max {dmax:.4}, |y| max {ymax:.4}");
    // riferimenti di docs/ZIGZAG.md F3 (regola, Python): raggio 0,433/0,438, ω 2,071
    let e1 = [wbx, wbz].iter().all(|w| (w / 2.071 - 1.0).abs() <= 0.03);
    let e2 = (abx / 0.433 - 1.0).abs() <= 0.10 && (abz / 0.438 - 1.0).abs() <= 0.10;
    let e3 = dmax <= 0.06;
    eprintln!("E1 frequenza {e1}; E2 raggio {e2}; E3 bagliori = onda {e3}");
    let serie: Vec<String> = (0..ts.len()).map(|i| format!("[{:.4}, {:.5}, {:.5}, {:.5}, {:.5}]", ts[i], bx[i], bz[i], ox[i], oz[i])).collect();
    format!("{{\n \"studio\": \"bagliori\",\n \"esami\": \"E1 E2 E3\",\n \"bagliori\": {{\"x\": {{\"A\": {abx}, \"omega\": {wbx}}}, \"z\": {{\"A\": {abz}, \"omega\": {wbz}}}}},\n \
        \"onda\": {{\"x\": {{\"A\": {aox}, \"omega\": {wox}}}, \"z\": {{\"A\": {aoz}, \"omega\": {woz}}}}},\n \"rif_ZIGZAG_F3\": {{\"A_x\": 0.433, \"A_z\": 0.438, \"omega\": 2.071}},\n \
        \"scarto_max_bagliori_onda\": {dmax},\n \"y_max\": {ymax},\n \"E1\": {e1},\n \"E2\": {e2},\n \"E3\": {e3},\n \"serie_t_bx_bz_ox_oz\": [{}]\n}}\n", serie.join(", "))
}

#[test]
#[ignore]
fn bagliori_tremito_elettrone() {
    // pre-registrato: spinore interpolato, trapezio nel tempo
    let j = tremito("elettrone (spinore interpolato)", 0);
    std::fs::write(format!("{OUT}/gl_tremito.json"), j).unwrap();
}

/// Dopo il FAIL di E3 (emendamento 1, docs/BAGLIORI.md): la stessa corsa con i bagliori sulla corrente esatta della
/// regola (link per link), e con Eulero sullo spinore interpolato (solo diagnosi).
#[test]
#[ignore]
fn bagliori_tremito_elettrone_flusso() {
    let j = tremito("elettrone (corrente esatta)", 1);
    std::fs::write(format!("{OUT}/gl_tremito_flusso.json"), j).unwrap();
    let j = tremito("elettrone (Eulero, spinore interpolato)", 2);
    std::fs::write(format!("{OUT}/gl_tremito_eulero.json"), j).unwrap();
}

/// I flussi letti (elscena::flussi) muovono esattamente ⟨r⟩ della regola: Σ_link F a = ⟨r⟩(t+1) − ⟨r⟩(t).
#[test]
fn bagliori_flussi_esatti() {
    let mut d = pacchetto(64, 0.25, 1.0, spinore_f3());
    for _ in 0..3 { d.tick(); }
    for _ in 0..4 {
        let f = elscena::flussi(&d);
        let c0 = centro_onda(&d);
        d.tick();
        let c1 = centro_onda(&d);
        for j in 0..3 {
            let s: f64 = f[j].iter().sum::<f64>() * d.a;
            assert!((s - (c1[j] - c0[j])).abs() < 1e-9, "asse {j}: flussi {s} contro d<r> {}", c1[j] - c0[j]);
        }
    }
}

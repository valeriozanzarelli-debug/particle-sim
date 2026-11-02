//! Esami G1–G5 dello studio motore_grafico_gpu (lab/out/mg_gpu_prereg.json + mg_gpu_amend.json).
//! Girano con `cargo test --release --bin atomo gpu_g -- --nocapture --test-threads 1` (serve Vulkan, anche
//! llvmpipe). Scrivono i risultati in lab/out/mg_gpu_*.json; il giudizio di G5 lo fa lab/mg_gpu_g5.py.

use std::time::Instant;

use crate::fotdisegno::Traccianti;
use crate::fotone3d::{Cf, Fotone3D, Lancio, Polar};
use crate::fotscena::{atteso, theta_da, Parametri, PALLINE_QUAL};
use crate::gpu_fotone::{poynting_cpu, somme_palline_cpu, spin_da_somme, GpuFotone};
use crate::gpu_regola::GpuHeadless;

const OUT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/lab/out");

fn gpu() -> Option<GpuHeadless> {
    match GpuHeadless::nuovo() {
        Ok(g) => { eprintln!("adattatore: {}", g.adattatore); Some(g) }
        Err(e) => { eprintln!("SALTATO: {e} (serve Vulkan, anche software llvmpipe)"); None }
    }
}

/// Una riga di misure (stesse sui due cammini).
#[derive(Clone, Copy, Debug, Default)]
struct Riga { energia: f64, centro: [f64; 3], poynting: [f64; 3], elicita: f64, spin: f64 }

fn json_riga(r: &Riga) -> String {
    format!("{{\"energia\": {:e}, \"centro\": [{:.6}, {:.6}, {:.6}], \"poynting\": [{:e}, {:e}, {:e}], \"elicita\": {:.6}, \"spin\": {:.6}}}",
            r.energia, r.centro[0], r.centro[1], r.centro[2], r.poynting[0], r.poynting[1], r.poynting[2], r.elicita, r.spin)
}

/// Un tick del cammino CPU con le sue misure (come fotscena::Mondo::passo).
fn passo_cpu(f: &mut Fotone3D, tr: &mut Traccianti, l: &Lancio) -> (Riga, crate::fotone3d::Misura) {
    let e = f.campo_e();
    f.tick();
    let rif = atteso(l, f.ticks);
    let m = f.misura(rif);
    let h = f.elicita(&e, l.dir, l.lambda);
    tr.avanza(f, 1.0, rif);
    let p = poynting_cpu(f);
    let (lr, vs) = somme_palline_cpu(&tr.pos, &tr.vel, rif);
    let s = spin_da_somme(lr, vs, rif, m.centro, l.dir, l.lambda, tr.pos.len());
    (Riga { energia: m.energia, centro: m.centro, poynting: p, elicita: h, spin: s }, m)
}

/// Un tick del cammino GPU con le sue misure (come gpu_fotone::FotoneGpu::passo).
fn passo_gpu(g: &mut GpuFotone, d: &GpuHeadless, l: &Lancio) -> (Riga, crate::fotone3d::Misura) {
    let rif = atteso(l, g.regola.ticks + 1);
    let r = g.passo(&d.device, &d.queue, true, rif, &[]);
    let h = r.elicita_num / ((std::f64::consts::TAU / l.lambda).sin() * r.elicita_den.max(1e-300));
    let s = spin_da_somme(r.l_rif, r.v_sum, rif, r.misura.centro, l.dir, l.lambda, g.n_palline);
    (Riga { energia: r.misura.energia, centro: r.misura.centro, poynting: r.poynting, elicita: h, spin: s }, r.misura)
}

fn nuovo_fotone(par: &Parametri) -> (Fotone3D, Lancio, Traccianti) {
    let n = par.n();
    let l = par.lancio();
    let mut f = Fotone3D::new(n);
    f.lancia(&l);
    let mut tr = Traccianti::new(&f, 1001, PALLINE_QUAL[par.qual]);
    tr.scia_n = crate::fotdisegno::SCIA_MAX;
    (f, l, tr)
}

fn norma3(v: [f64; 3]) -> f64 { (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() }

/// G1 (campo dopo 96 tick), G2 (palline dopo 48 tick), G3 (misure tick per tick) al preset alta.
#[test]
fn gpu_g123_alta() {
    let Some(d) = gpu() else { return };
    let par = Parametri { qual: 3, ..Default::default() };
    let (mut f, l, mut tr) = nuovo_fotone(&par);
    let n = f.n;
    let mut g = GpuFotone::nuovo(&d.device, n, tr.pos.len(), tr.scia_n);
    g.carica(&d.device, &d.queue, &f, &tr, &l);
    // riscaldamento della GPU (compilazione) su una copia, poi ricarica
    let _ = g.passo(&d.device, &d.queue, true, atteso(&l, 1), &[]);
    g.carica(&d.device, &d.queue, &f, &tr, &l);
    let f0 = f.f.clone();
    let (mut ms_cpu, mut ms_gpu) = (0.0, 0.0);
    let (mut d_e, mut d_p, mut d_s_48, mut d_s_max, mut d_c, mut d_h) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut righe = vec![];
    let (mut g2_media, mut g2_max, mut g2_conf, mut rinate_cpu, mut rinate_gpu) = (0.0, 0.0f64, 0usize, 0usize, 0usize);
    for t in 1..=96u64 {
        let t0 = Instant::now();
        let (rc, _) = passo_cpu(&mut f, &mut tr, &l);
        ms_cpu += t0.elapsed().as_secs_f64() * 1e3;
        let t0 = Instant::now();
        let (rg, _) = passo_gpu(&mut g, &d, &l);
        ms_gpu += t0.elapsed().as_secs_f64() * 1e3;
        if t <= 48 {
            d_e = d_e.max((rg.energia / rc.energia - 1.0).abs());
            let dp = [rg.poynting[0] - rc.poynting[0], rg.poynting[1] - rc.poynting[1], rg.poynting[2] - rc.poynting[2]];
            d_p = d_p.max(norma3(dp) / norma3(rc.poynting));
            d_c = d_c.max(norma3([rg.centro[0] - rc.centro[0], rg.centro[1] - rc.centro[1], rg.centro[2] - rc.centro[2]]));
            d_h = d_h.max((rg.elicita - rc.elicita).abs());
            let ds = ((rg.spin - rc.spin) / rc.spin).abs();
            if t >= 24 { d_s_max = d_s_max.max(ds); }
            if t == 48 { d_s_48 = ds; }
            righe.push(format!("  {{\"tick\": {t}, \"cpu\": {}, \"gpu\": {}}}", json_riga(&rc), json_riga(&rg)));
        }
        if t == 48 {
            let (pg, _, eg) = g.leggi_palline(&d.device, &d.queue);
            let mut somma = 0.0;
            for i in 0..pg.len() {
                let rc_ = tr.eta[i] < 30_000;
                let rg_ = eg[i] < 30_000;
                if rc_ { rinate_cpu += 1; }
                if rg_ { rinate_gpu += 1; }
                if rc_ || rg_ { continue; }
                let p = tr.pos[i];
                let q = pg[i];
                let dd = (((p[0] - q[0]) as f64).powi(2) + ((p[1] - q[1]) as f64).powi(2) + ((p[2] - q[2]) as f64).powi(2)).sqrt();
                somma += dd;
                g2_max = g2_max.max(dd);
                g2_conf += 1;
            }
            g2_media = somma / g2_conf.max(1) as f64;
            eprintln!("G2: {g2_conf} palline confrontate, scarto medio {:.3e} celle = {:.3e} λ, massimo {:.3e} celle; rinate CPU {rinate_cpu}, GPU {rinate_gpu}",
                      g2_media, g2_media / l.lambda, g2_max);
        }
        if t % 12 == 0 { eprintln!("tick {t}: ms CPU {:.0}, GPU {:.0} (medie)", ms_cpu / t as f64, ms_gpu / t as f64); }
    }
    // G1: il campo dopo 96 tick
    let fg = g.regola.leggi(&d.device, &d.queue);
    let rel = diff_rel(&fg, &f.f);
    let mosso = diff_rel(&f.f, &f0);
    let np = tr.pos.len();
    let g1 = rel <= 1e-4;
    let g2_info = (rinate_cpu.max(rinate_gpu) as f64) <= 0.01 * np as f64;
    let g2 = g2_media / l.lambda <= 0.01;
    let g3 = d_e <= 1e-4 && d_p <= 1e-4 && d_s_48 <= 5e-3 && d_s_max <= 5e-3;
    eprintln!("G1: max|F_gpu - F_cpu|/max|F_cpu| = {rel:.3e} (campo mosso {mosso:.3}) -> {}", if g1 { "PASSA" } else { "FALLISCE" });
    eprintln!("G2: scarto medio {:.3e} λ -> {}{}", g2_media / l.lambda, if g2 { "PASSA" } else { "FALLISCE" }, if g2_info { "" } else { " (non informativo: troppe rinate)" });
    eprintln!("G3: energia {d_e:.2e}, Poynting {d_p:.2e}, spin t48 {d_s_48:.2e}, spin max 24..48 {d_s_max:.2e}; centro {d_c:.2e} celle, elicità {d_h:.2e} -> {}",
              if g3 { "PASSA" } else { "FALLISCE" });
    let json = format!(
        "{{\n \"studio\": \"motore_grafico_gpu\",\n \"esami\": \"G1 G2 G3\",\n \"adattatore\": \"{}\",\n \"preset\": \"alta\",\n \"n\": {n},\n \"celle_per_lambda\": {},\n \"palline\": {np},\n \"scia\": {},\n \
         \"G1\": {{\"tick\": 96, \"diff_rel_max\": {rel:e}, \"soglia\": 1e-4, \"campo_mosso\": {mosso:.4}, \"esito\": \"{}\"}},\n \
         \"G2\": {{\"tick\": 48, \"palline_confrontate\": {g2_conf}, \"scarto_medio_celle\": {g2_media:e}, \"scarto_medio_lambda\": {:e}, \"scarto_max_celle\": {g2_max:e}, \"rinate_cpu\": {rinate_cpu}, \"rinate_gpu\": {rinate_gpu}, \"informativo\": {g2_info}, \"soglia_lambda\": 0.01, \"esito\": \"{}\"}},\n \
         \"G3\": {{\"tick\": \"1..48\", \"energia_rel_max\": {d_e:e}, \"poynting_rel_max\": {d_p:e}, \"spin_rel_t48\": {d_s_48:e}, \"spin_rel_max_24_48\": {d_s_max:e}, \"centro_max_celle\": {d_c:e}, \"elicita_max\": {d_h:e}, \"soglie\": \"energia 1e-4, Poynting 1e-4, spin 5e-3\", \"esito\": \"{}\"}},\n \
         \"tempi_ms_per_tick\": {{\"cpu_rayon\": {:.1}, \"gpu\": {:.1}, \"cpu_core\": {}, \"nota\": \"llvmpipe e' Vulkan su CPU: non dice niente di una GPU vera\"}},\n \"righe\": [\n{}\n ]\n}}\n",
        d.adattatore, par.cpl(), tr.scia_n, if g1 { "passa" } else { "fallisce" }, g2_media / l.lambda, if g2 { "passa" } else { "fallisce" },
        if g3 { "passa" } else { "fallisce" }, ms_cpu / 96.0, ms_gpu / 96.0, std::thread::available_parallelism().map(|x| x.get()).unwrap_or(0), righe.join(",\n"));
    std::fs::write(format!("{OUT}/mg_gpu_g123.json"), json).unwrap();
    assert!(g1 && g2 && g3, "G1 {g1} G2 {g2} G3 {g3}");
}

fn diff_rel(a: &[[Cf; 3]], b: &[[Cf; 3]]) -> f64 {
    let (mut dmax, mut fmax) = (0.0f64, 0.0f64);
    for (p, q) in a.iter().zip(b) {
        for c in 0..3 {
            let (dr, di) = ((p[c].re - q[c].re) as f64, (p[c].im - q[c].im) as f64);
            dmax = dmax.max((dr * dr + di * di).sqrt());
            fmax = fmax.max((q[c].norm2() as f64).sqrt());
        }
    }
    dmax / fmax.max(1e-300)
}

/// G5: la stessa scena a bassa (cammino CPU) e ad alta (cammino GPU), 96 tick; scrive le serie, il giudizio
/// è in lab/mg_gpu_g5.py.
fn g5(qual: usize, polar: Polar, nome: &str) {
    let par = Parametri { qual, polar, ..Default::default() };
    let (mut f, l, mut tr) = nuovo_fotone(&par);
    let mut storia = vec![f.misura(atteso(&l, 0))];
    let mut righe = vec![];
    let cammino;
    let t0 = Instant::now();
    if qual >= 3 {
        let Some(d) = gpu() else { return };
        cammino = "gpu";
        let mut g = GpuFotone::nuovo(&d.device, f.n, tr.pos.len(), tr.scia_n);
        g.carica(&d.device, &d.queue, &f, &tr, &l);
        for _ in 1..=96 { let (r, m) = passo_gpu(&mut g, &d, &l); storia.push(m); righe.push(r); }
    } else {
        cammino = "cpu";
        for _ in 1..=96 { let (r, m) = passo_cpu(&mut f, &mut tr, &l); storia.push(m); righe.push(r); }
    }
    let ms = t0.elapsed().as_secs_f64() * 1e3 / 96.0;
    let theta = theta_da(&storia, l.dir);
    let th_att = l.lambda / (std::f64::consts::PI * l.w0);
    eprintln!("G5 {nome}: cammino {cammino}, theta {theta:?} (Maxwell {th_att:.4}), spin al tick 96 {:.4}, {ms:.0} ms/tick", righe[95].spin);
    let r: Vec<String> = righe.iter().enumerate().map(|(t, r)| format!("  {{\"tick\": {}, \"m\": {}}}", t + 1, json_riga(r))).collect();
    let json = format!(
        "{{\n \"studio\": \"motore_grafico_gpu\",\n \"esame\": \"G5\",\n \"preset_qual\": {qual},\n \"celle_per_lambda\": {},\n \"n\": {},\n \"palline\": {},\n \"polar\": \"{nome}\",\n \"cammino\": \"{cammino}\",\n \
         \"dir\": [{}, {}, {}],\n \"lambda_celle\": {},\n \"w0_celle\": {},\n \"theta\": {},\n \"theta_maxwell\": {th_att},\n \"ms_tick\": {ms:.1},\n \"righe\": [\n{}\n ]\n}}\n",
        par.cpl(), f.n, tr.pos.len(), l.dir[0], l.dir[1], l.dir[2], l.lambda, l.w0, theta.map(|t| t.to_string()).unwrap_or("null".into()), r.join(",\n"));
    std::fs::write(format!("{OUT}/mg_gpu_g5_{}_{nome}.json", if qual >= 3 { "alta" } else { "bassa" }), json).unwrap();
}

#[test] fn gpu_g5_bassa_circp() { g5(0, Polar::CircolarePiu, "circp"); }
#[test] fn gpu_g5_bassa_circm() { g5(0, Polar::CircolareMeno, "circm"); }
#[test] fn gpu_g5_bassa_lin() { g5(0, Polar::Orizzontale, "lin"); }
#[test] fn gpu_g5_alta_circp() { g5(3, Polar::CircolarePiu, "circp"); }
#[test] fn gpu_g5_alta_circm() { g5(3, Polar::CircolareMeno, "circm"); }
#[test] fn gpu_g5_alta_lin() { g5(3, Polar::Orizzontale, "lin"); }

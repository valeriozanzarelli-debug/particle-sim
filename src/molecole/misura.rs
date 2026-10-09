//! The measurements: every number of docs/MOLECOLE.md comes from here, with the ball rule of dmc.rs.
//!
//! `molecole --misura <bersaglio>` with bersaglio = idrogeno | eccitati | luce | h2piu | h2 | elio | tutto.
//! Each target writes its numbers into assets/molecole/risultati.txt (read by the viewer and the tests).
//! Seeds are fixed: rerunning gives the same numbers.

use crate::costanti as k;
use crate::dmc::*;
use std::collections::BTreeMap;
use std::time::Instant;

pub const FILE: &str = "assets/molecole/risultati.txt";

// ───────────────────────── results file ─────────────────────────
#[derive(Default, Clone, Debug)]
pub struct Risultati(pub BTreeMap<String, String>);

impl Risultati {
    pub fn parse(text: &str) -> Self {
        let mut m = BTreeMap::new();
        for l in text.lines() {
            let l = l.split('#').next().unwrap_or("").trim();
            if let Some((a, b)) = l.split_once('=') {
                m.insert(a.trim().to_string(), b.trim().to_string());
            }
        }
        Risultati(m)
    }
    pub fn load(path: &str) -> Self {
        std::fs::read_to_string(path).map(|t| Self::parse(&t)).unwrap_or_default()
    }
    pub fn save(&self, path: &str) {
        let mut s = String::from("# Numbers computed by `molecole --misura` (ball rule, dmc.rs). Do not edit by hand.\n");
        for (a, b) in &self.0 { s += &format!("{a} = {b}\n"); }
        if let Some(dir) = std::path::Path::new(path).parent() { let _ = std::fs::create_dir_all(dir); }
        std::fs::write(path, s).expect("cannot write results");
    }
    pub fn f(&self, key: &str) -> Option<f64> { self.0.get(key).and_then(|v| v.parse().ok()) }
    pub fn put(&mut self, key: &str, v: f64) { self.0.insert(key.to_string(), format!("{v:.9e}")); }
    pub fn put_s(&mut self, key: &str, v: &str) { self.0.insert(key.to_string(), v.to_string()); }
}

// ───────────────────────── running the rule ─────────────────────────
pub struct Esito {
    /// mixed estimator (⟨V⟩ without guide, ⟨E_loc⟩ with guide), error
    pub e: (f64, f64),
    /// growth estimator, error
    pub g: (f64, f64),
}

/// Equilibrate for `t_eq`, then measure for `t_run` (atomic units of time); `obs` sees the balls every 10 ticks.
pub fn corri(d: &mut Dmc, t_eq: f64, t_run: f64, mut obs: impl FnMut(&Dmc)) -> Esito {
    let neq = (t_eq / d.tau).round() as usize;
    let nrun = (t_run / d.tau).round() as usize;
    for _ in 0..neq { d.passo(); }
    let (mut em, mut eg) = (Vec::with_capacity(nrun), Vec::with_capacity(nrun));
    for i in 0..nrun {
        let (a, b) = d.passo();
        em.push(a);
        eg.push(b);
        if i % 10 == 0 { obs(d); }
    }
    Esito { e: media_errore(&em), g: media_errore(&eg) }
}

/// Weighted straight-line fit E(τ) = E0 + b τ; returns (E0, error of E0).
pub fn estrapola(taus: &[f64], e: &[(f64, f64)]) -> (f64, f64) {
    let (mut s, mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (t, (y, dy)) in taus.iter().zip(e) {
        let w = 1.0 / (dy * dy).max(1e-30);
        s += w; sx += w * t; sy += w * y; sxx += w * t * t; sxy += w * t * y;
    }
    let det = s * sxx - sx * sx;
    let a = (sxx * sy - sx * sxy) / det;
    (a, (sxx / det).sqrt())
}

/// Radial histogram of electron 0 (relative to the origin), bin width `dr`, up to `rmax`.
pub struct Istogramma { pub dr: f64, pub h: Vec<f64>, pub n: f64 }

impl Istogramma {
    pub fn new(dr: f64, rmax: f64) -> Self { Istogramma { dr, h: vec![0.0; (rmax / dr) as usize], n: 0.0 } }
    pub fn add(&mut self, d: &Dmc) {
        for p in d.elettrone(0) {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let i = (r / self.dr) as usize;
            if i < self.h.len() { self.h[i] += 1.0; }
            self.n += 1.0;
        }
    }
    /// Ball density per unit volume (normalised to the balls counted), at bin centres.
    pub fn densita(&self) -> Vec<(f64, f64)> {
        self.h.iter().enumerate().map(|(i, c)| {
            let r0 = i as f64 * self.dr;
            let vol = 4.0 / 3.0 * std::f64::consts::PI * ((r0 + self.dr).powi(3) - r0.powi(3));
            (r0 + 0.5 * self.dr, c / (self.n * vol))
        }).collect()
    }
    /// Fit ln(ρ / f(r)) = c − a r on [r0, r1] (weights = counts); `f` removes the known polynomial part.
    pub fn decadimento(&self, r0: f64, r1: f64, f: impl Fn(f64) -> f64) -> f64 {
        let (mut s, mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for ((r, rho), c) in self.densita().into_iter().zip(&self.h) {
            if r < r0 || r > r1 || *c < 20.0 { continue; }
            let y = (rho / f(r)).ln();
            let w = *c;
            s += w; sx += w * r; sy += w * y; sxx += w * r * r; sxy += w * r * y;
        }
        -(s * sxy - sx * sy) / (s * sxx - sx * sx)
    }
}

fn idrogeno_sys(guida: Guida, vincoli: Vec<Vincolo>, campo: f64) -> Sistema {
    Sistema { nuclei: vec![Nucleo { pos: [0.0; 3], z: 1.0 }], spin: vec![1], guida, vincoli, campo_z: campo }
}

fn ev(e: f64) -> f64 { e * k::hartree_ev() }

// ───────────────────────── 1. hydrogen atom ─────────────────────────
pub fn idrogeno(res: &mut Risultati) {
    println!("\n=== 1. Hydrogen atom: one electron, one proton (rule without guide) ===");
    let taus = [0.02, 0.01, 0.005];
    let mut es = vec![];
    let mut hist = Istogramma::new(0.05, 12.0);
    for (i, &tau) in taus.iter().enumerate() {
        let t = Instant::now();
        let mut d = Dmc::new(idrogeno_sys(Guida::Nessuna, vec![], 0.0), tau, 20000, 11 + i as u64);
        let last = i == taus.len() - 1;
        let r = corri(&mut d, 20.0, 150.0, |d| if last { hist.add(d) });
        println!("  tau {tau:<6} E = {:.6} ± {:.6} hartree (growth {:.6})   [{:.0} s]", r.e.0, r.e.1, r.g.0, t.elapsed().as_secs_f64());
        res.put(&format!("h_plain_E_tau{tau}"), r.e.0);
        es.push(r.e);
    }
    let (e0, de0) = estrapola(&taus, &es);
    let a = hist.decadimento(0.3, 5.0, |r| r.powi(0));
    println!("  tau -> 0: E = {e0:.6} ± {de0:.6} hartree = {:.4} ± {:.4} eV   (exact -13.6057 eV)", ev(e0), ev(de0));
    println!("  ball density ∝ exp(-a r): a = {a:.4} /a0   (the cloud of unguided balls is ψ itself: exact ψ ∝ e^(-r/a0))");
    res.put("h_E_plain", e0);
    res.put("h_E_plain_err", de0);
    res.put("h_a", a);

    println!("  second generation: balls guided by OUR measured orbital exp(-a r)");
    let gt = [0.02, 0.01, 0.005];
    let mut ge = vec![];
    let mut hist2 = Istogramma::new(0.05, 12.0);
    let mut d = Dmc::new(idrogeno_sys(Guida::Orbitali { a, jastrow: false }, vec![], 0.0), 0.01, 20000, 21);
    for (i, &tau) in gt.iter().enumerate() {
        d = Dmc::new(idrogeno_sys(Guida::Orbitali { a, jastrow: false }, vec![], 0.0), tau, 20000, 21 + i as u64);
        let last = i == gt.len() - 1;
        let r = corri(&mut d, 10.0, 60.0, |d| if last { hist2.add(d) });
        println!("  guided tau {tau:<6} E = {:.7} ± {:.7}", r.e.0, r.e.1);
        ge.push(r.e);
    }
    let r = Esito { e: estrapola(&gt, &ge), g: (0.0, 0.0) };
    let a2 = hist2.decadimento(0.3, 5.0, |r| r.powi(0));
    println!("  guided, tau->0: E = {:.7} ± {:.7} hartree = {:.5} eV; ball density ∝ exp(-{a2:.4} r) (|ψ|² exact: e^(-2r/a0))",
             r.e.0, r.e.1, ev(r.e.0));
    // radial profile versus the exact |ψ|² (normalised): largest deviation for r < 4 a0
    let mut dev: f64 = 0.0;
    for ((rr, rho), c) in hist2.densita().into_iter().zip(&hist2.h) {
        if rr > 0.2 && rr < 4.0 && *c > 100.0 {
            let exact = (-2.0 * rr).exp() / std::f64::consts::PI;
            dev = dev.max((rho - exact).abs() / exact);
        }
    }
    println!("  profile check: max |ρ_balls − e^(-2r)/π| / exact = {:.2} % for 0.2 < r < 4 a0 (closer to the proton the
   finite tick distorts the cloud: first bin -8 % at tau 0.005, shrinking with tau)", 100.0 * dev);
    res.put("h_E_guided", r.e.0);
    res.put("h_E_guided_err", r.e.1);
    res.put("h_rho_decay", a2);
    res.put("h_profile_dev", dev);
    let p: Vec<[f64; 3]> = d.elettrone(0).take(1).collect();
    println!("  surfer (real stationary ψ: Bohm velocity ∇S = 0, it stays where it was born): at {:?} a0", p[0]);
    res.put_s("h_surfer", &format!("{:.4} {:.4} {:.4}", p[0][0], p[0][1], p[0][2]));
}

// ───────────────────────── 2. excited levels and emitted light ─────────────────────────
struct Stato { nome: &'static str, n: u32, ang: Angolare, piani: Vec<Vincolo>, l: i32 }

fn stato_guidato(st: &Stato, a: f64, nodo: Option<f64>, taus: &[f64], t_run: f64, seed: u64) -> (f64, f64) {
    let mut es = vec![];
    for (i, &tau) in taus.iter().enumerate() {
        let (nodi, nn) = match nodo { Some(r) => ([r, 0.0], 1), None => ([0.0; 2], 0) };
        let mut v = st.piani.clone();
        if let Some(r) = nodo { v.push(Vincolo::SferaFuori(r)); }
        let sys = idrogeno_sys(Guida::Eccitato { ang: st.ang, nodi, nn, a }, v, 0.0);
        let mut d = Dmc::new(sys, tau, 10000, seed + i as u64);
        let r = corri(&mut d, 40.0 * (st.n * st.n) as f64, t_run, |_| {});
        println!("     {} guided tau {tau:<5} E = {:.7} ± {:.7}", st.nome, r.e.0, r.e.1);
        es.push(r.e);
    }
    estrapola(taus, &es)
}

/// Unguided balls in one pocket; returns (growth energy, error, fitted decay of the cloud).
fn tasca(st: &Stato, extra: Option<Vincolo>, tau: f64, t_run: f64, seed: u64, nodo: Option<f64>) -> (f64, f64, f64) {
    let mut v = st.piani.clone();
    if let Some(e) = extra { v.push(e); }
    let mut d = Dmc::new(idrogeno_sys(Guida::Nessuna, v, 0.0), tau, 10000, seed);
    let mut h = Istogramma::new(0.1, 60.0);
    let r = corri(&mut d, 40.0 * (st.n * st.n) as f64, t_run, |d| h.add(d));
    let l = st.l;
    let n = st.n as f64;
    let a = h.decadimento(n * 1.0, n * 5.0, |r| r.powi(l) * nodo.map_or(1.0, |rn| (r - rn).abs()));
    (r.g.0, r.g.1, a)
}

pub fn eccitati(res: &mut Risultati, quale: &str) {
    println!("\n=== 2. Excited levels of hydrogen ({}) ===",
             if quale == "simm" { "symmetry nodes: 2p, 3d" } else { "radial node, found by the balls" });
    let st2p = Stato { nome: "2p", n: 2, ang: Angolare::Z, piani: vec![Vincolo::PianoZ], l: 1 };
    let st3d = Stato { nome: "3d", n: 3, ang: Angolare::XY, piani: vec![Vincolo::PianoX, Vincolo::PianoY], l: 2 };
    let st2s = Stato { nome: "2s", n: 2, ang: Angolare::S, piani: vec![], l: 0 };
    let st3p = Stato { nome: "3p", n: 3, ang: Angolare::Z, piani: vec![Vincolo::PianoZ], l: 1 };
    // symmetry nodes: 2p (plane z = 0), 3d_xy (planes x = 0 and y = 0)
    let simm: &[(&Stato, u64)] = if quale != "simm" { &[] } else { &[(&st2p, 100u64), (&st3d, 200)] };
    for &(st, seed) in simm {
        let t = Instant::now();
        let (eg, deg, a) = tasca(st, None, 0.02, 300.0, seed, None);
        println!("  {}: unguided balls  E = {eg:.6} ± {deg:.6}  cloud decay a = {a:.4}  [{:.0} s]", st.nome, t.elapsed().as_secs_f64());
        let (eq, deq) = stato_guidato(st, a, None, &[0.04, 0.02, 0.01], 200.0, seed + 10);
        println!("  {}: guided, tau->0  E = {eq:.7} ± {deq:.7} hartree   (exact -1/{}) [{:.0} s]",
                 st.nome, st.n * st.n * 2, t.elapsed().as_secs_f64());
        res.put(&format!("h{}_E_plain", st.nome), eg);
        res.put(&format!("h{}_E_plain_err", st.nome), deg);
        res.put(&format!("h{}_a", st.nome), a);
        res.put(&format!("h{}_E", st.nome), eq);
        res.put(&format!("h{}_E_err", st.nome), deq);
    }
    // radial nodes: the node sits where the two pockets (inside / outside) have the SAME energy
    let rad: &[(&Stato, [f64; 3], u64)] = match quale {
        "2s" => &[(&st2s, [1.6, 2.0, 2.4], 300u64)],
        "3p" => &[(&st3p, [5.4, 6.0, 6.6], 400)],
        _ => &[],
    };
    for &(st, cands, seed) in rad {
        let t = Instant::now();
        // two rounds: coarse candidates, then ±5 % around the first crossing; parabola through 3 points
        let mut cs = cands.to_vec();
        let (mut rn, mut ecross, mut drn, mut so) = (0.0, 0.0, 0.0, 0.0);
        for round in 0..2u64 {
            let mut pts = vec![];
            let mut errs = vec![];
            let tr = if round == 0 { 0.5 } else { 1.0 };
            for (i, &r) in cs.iter().enumerate() {
                let sd = seed + 10 * round + 2 * i as u64;
                let (ei, dei, _) = tasca(st, Some(Vincolo::SferaDentro(r)), 0.0025 * (st.n * st.n) as f64, 150.0 * tr, sd, None);
                let (eo, deo, _) = tasca(st, Some(Vincolo::SferaFuori(r)), 0.02, 300.0 * tr, sd + 1, None);
                println!("  {}: node sphere r = {r:.3}: inside E = {ei:.5} ± {dei:.5}, outside E = {eo:.5} ± {deo:.5}", st.nome);
                pts.push((r, ei, eo));
                errs.push((dei, deo));
            }
            let par = |f: &dyn Fn(&(f64, f64, f64)) -> f64, x: f64| {
                // Lagrange parabola through the three points
                let mut y = 0.0;
                for i in 0..3 {
                    let mut l = 1.0;
                    for j in 0..3 { if i != j { l *= (x - pts[j].0) / (pts[i].0 - pts[j].0); } }
                    y += l * f(&pts[i]);
                }
                y
            };
            // bisection for E_in(r) = E_out(r) on the candidate interval
            let (mut lo, mut hi) = (cs[0], cs[2]);
            let g = |x: f64| par(&|p| p.1, x) - par(&|p| p.2, x);
            for _ in 0..60 { let m = 0.5 * (lo + hi); if g(lo) * g(m) <= 0.0 { hi = m } else { lo = m } }
            rn = 0.5 * (lo + hi);
            ecross = par(&|p| p.2, rn);
            // error of the node position: errors of the two pocket energies / difference of their slopes
            let h = 1e-4;
            let si = (par(&|p| p.1, rn + h) - par(&|p| p.1, rn - h)) / (2.0 * h);
            so = (par(&|p| p.2, rn + h) - par(&|p| p.2, rn - h)) / (2.0 * h);
            let e2 = errs.iter().map(|e| e.0 * e.0 + e.1 * e.1).sum::<f64>() / 3.0;
            drn = e2.sqrt() / (si - so).abs();
            cs = vec![rn * 0.95, rn, rn * 1.05];
        }
        // second generation at the crossing node: guided, outer pocket
        let (_, _, a) = tasca(st, Some(Vincolo::SferaFuori(rn)), 0.02, 200.0, seed + 50, Some(rn));
        let (eq, deq) = stato_guidato(st, a, Some(rn), &[0.04, 0.02], 150.0, seed + 60);
        let de_node = so.abs() * drn;
        let derr = (deq * deq + de_node * de_node).sqrt();
        println!("  {}: node found at r = {rn:.4} ± {drn:.4} a0 (the two pockets agree), E(cross) = {ecross:.6}", st.nome);
        println!("  {}: guided E = {eq:.7} ± {derr:.7} hartree (node uncertainty {de_node:.1e} included)  (exact -1/{}) [{:.0} s]",
                 st.nome, 2 * st.n * st.n, t.elapsed().as_secs_f64());
        res.put(&format!("h{}_node", st.nome), rn);
        res.put(&format!("h{}_node_err", st.nome), drn);
        res.put(&format!("h{}_E", st.nome), eq);
        res.put(&format!("h{}_E_err", st.nome), derr);
    }
}

/// Photons emitted between the computed levels (whatever levels are available in the results).
pub fn righe(res: &mut Risultati) {
    let mut e = BTreeMap::new();
    for st in ["1s", "2p", "3d", "2s", "3p"] {
        let key = if st == "1s" { "h_E_guided".to_string() } else { format!("h{st}_E") };
        if let (Some(v), Some(dv)) = (res.f(&key), res.f(&format!("{key}_err"))) { e.insert(st, (v, dv)); }
    }
    // emitted light: ΔE = h c / λ, with the finite proton mass (two-body problem: energies × μ/m_e)
    let mu = k::proton_me() / (k::proton_me() + 1.0);
    let lines = [("Lyman-alpha", "2p", "1s", 121.567), ("Lyman-beta", "3p", "1s", 102.572), ("Balmer-alpha", "3d", "2p", 656.461)];
    println!("  emitted photons (vacuum wavelengths; NIST ASD):");
    for (name, up, lo, nist) in lines {
        let (Some(&(eu, deu)), Some(&(el, del))) = (e.get(up), e.get(lo)) else { continue };
        let de = (eu - el) * mu;
        let dde = (deu * deu + del * del).sqrt() * mu;
        let lam = k::nm_from_hartree(de);
        let dlam = lam * dde / de;
        println!("    {name:13} {up}->{lo}: λ = {lam:.4} ± {dlam:.4} nm   NIST {nist} nm   error {:+.3} nm", lam - nist);
        let key = name.to_lowercase().replace('-', "_");
        res.put(&format!("line_{key}_nm"), lam);
        res.put(&format!("line_{key}_err"), dlam);
        res.put(&format!("line_{key}_nist"), nist);
    }
    for (a, b) in [("2s", "2p"), ("3p", "3d")] {
        if let (Some(x), Some(y)) = (e.get(a), e.get(b)) {
            println!("  degeneracy (same n, same energy): E({a}) - E({b}) = {:.7} ± {:.7} hartree", x.0 - y.0, (x.1 * x.1 + y.1 * y.1).sqrt());
            res.put(&format!("deg_{a}_{b}"), x.0 - y.0);
        }
    }
}

// ───────────────────────── 2b. hydrogen and light: polarizability, line strength ─────────────────────────
pub fn luce(res: &mut Risultati, parte: &str) {
    println!("\n=== 2b. Hydrogen and light: {} ===",
             if parte == "campo" { "polarizability" } else { "Lyman-alpha strength, cross-sections" });
    let a = res.f("h_a").expect("run `--misura idrogeno` first");
    if parte == "campo" { polarizzabilita(res, a); } else { riga(res, a); }
}

fn polarizzabilita(res: &mut Risultati, a: f64) {
    // (a) energy shift in a weak uniform field: E(F) = E(0) − ½ α F² − γ F⁴/24
    let fs = [0.0, 0.01, 0.02, 0.03, 0.04];
    let mut pts = vec![];
    for (i, &f) in fs.iter().enumerate() {
        let t = Instant::now();
        let mut d = Dmc::new(idrogeno_sys(Guida::Orbitali { a, jastrow: false }, vec![], f), 0.01, 20000, 500 + i as u64);
        let mut zsum = (0.0, 0.0);
        let r = corri(&mut d, 20.0, 150.0, |d| { for p in d.elettrone(0) { zsum.0 += p[2]; zsum.1 += 1.0; } });
        let zm = zsum.0 / zsum.1;
        println!("  F = {f:.2} a.u.: E = {:.8} ± {:.8}   <z>_balls = {zm:+.5}  [{:.0} s]", r.e.0, r.e.1, t.elapsed().as_secs_f64());
        pts.push((f, r.e.0, r.e.1, zm));
    }
    // least squares ΔE = −½α F² + c F⁴ (E(0) as a free parameter too)
    let (mut ata, mut atb) = ([[0.0f64; 3]; 3], [0.0f64; 3]);
    for &(f, e, de, _) in &pts {
        let row = [1.0, f * f, f.powi(4)];
        let w = 1.0 / (de * de).max(1e-20);
        for i in 0..3 { atb[i] += w * row[i] * e; for j in 0..3 { ata[i][j] += w * row[i] * row[j]; } }
    }
    let inv = inv3(ata);
    let c1: f64 = (0..3).map(|j| inv[1][j] * atb[j]).sum();
    let alpha = -2.0 * c1;
    let dalpha = 2.0 * inv[1][1].sqrt();
    // induced dipole of the cloud (mixed estimate, extrapolated: 2·mixed − guide, guide is symmetric)
    let alpha_dip: Vec<f64> = pts.iter().skip(1).map(|p| -2.0 * p.3 / p.0).collect();
    println!("  α (energy, F→0) = {alpha:.3} ± {dalpha:.3} a0³   (exact 4.5)");
    println!("  α from the shift of the cloud, 2<z>/F at F = 0.01..0.04: {:?}", alpha_dip.iter().map(|x| (x * 1000.0).round() / 1000.0).collect::<Vec<_>>());
    res.put("h_alpha", alpha);
    res.put("h_alpha_err", dalpha);
    res.put("h_alpha_dip", alpha_dip[0]);
}

fn riga(res: &mut Risultati, a: f64) {
    let alpha = res.f("h_alpha").expect("run `--misura campo` first");

    // (b) dipole matrix element <1s|z|2p0> from the two clouds of guided balls (ball density = ψ_g·ψ)
    let a2 = res.f("h2p_a").expect("run `--misura eccitati` first");
    // ψ = (ball density) / guide holds exactly only for tau -> 0: two pairs of ticks, straight-line extrapolation
    let mut zt = vec![];
    for (ip, &(t1, t2)) in [(0.01, 0.02), (0.005, 0.01)].iter().enumerate() {
        let (z, dz) = elemento_z12(a, a2, t1, t2, 600 + 10 * ip as u64);
        println!("  <1s|z|2p0> with ticks tau(1s) = {t1}, tau(2p) = {t2}: {z:.5} ± {dz:.5} a0");
        zt.push((t2, (z, dz)));
    }
    let (zm, dz) = estrapola(&[zt[0].0, zt[1].0], &[zt[0].1, zt[1].1]);
    let e1 = res.f("h_E_guided").unwrap();
    let e2 = res.f("h2p_E").unwrap();
    let mu = k::proton_me() / (k::proton_me() + 1.0);
    let de = (e2 - e1) * mu;
    // f(1s→2p) = (2/3) ΔE Σ_m |<1s|r|2p_m>|² = 2 ΔE |<1s|z|2p0>|² (a.u., z in reduced-mass bohr)
    let z_mu = zm / mu;
    let fosc = 2.0 * de * z_mu * z_mu;
    // A(2p→1s) = 4 ω³ |d|² / (3 c³), c = 1/α_fs in a.u.
    let c_au = 1.0 / k::alpha_fs();
    let a_au = 4.0 * de.powi(3) * z_mu * z_mu / (3.0 * c_au.powi(3));
    let a_s = a_au / k::au_time_s();
    let da_s = a_s * 2.0 * dz / zm;
    println!("  <1s|z|2p0> from the clouds = {zm:.4} ± {dz:.4} a0   (exact 0.7449)");
    println!("  oscillator strength f(1s→2p) = {fosc:.4} ± {:.4}  (exact 0.4162)", fosc * 2.0 * dz / zm);
    println!("  Einstein A(2p→1s) = {a_s:.4e} ± {da_s:.1e} 1/s  (exact 6.2649e8)   lifetime {:.4} ns (exact 1.596)", 1e9 / a_s);
    res.put("h_z12", zm);
    res.put("h_z12_err", dz);
    res.put("h_f12", fosc);
    res.put("h_A21", a_s);
    res.put("h_A21_err", da_s);
    // (c) the "rules for light" of the hydrogen atom
    let a0 = k::bohr_m();
    let lam = 500e-9;
    let kk = 2.0 * std::f64::consts::PI / lam;
    let ray = |al: f64| 8.0 * std::f64::consts::PI / 3.0 * kk.powi(4) * (al * a0.powi(3)).powi(2);
    let lam_ly = res.f("line_lyman_alpha_nm").unwrap_or(121.567) * 1e-9;
    let sres = 3.0 * lam_ly * lam_ly / (2.0 * std::f64::consts::PI);
    println!("  Rayleigh σ(500 nm, static α) = {:.3e} m²  (with exact α = 4.5: {:.3e} m²)", ray(alpha), ray(4.5));
    println!("  resonant σ at Lyman-alpha = 3λ²/2π = {sres:.3e} m²   (exact λ: {:.3e})", 3.0 * 121.567e-9f64.powi(2) / (2.0 * std::f64::consts::PI));
    println!("  max resonant scattering rate Γ/2 = {:.3e} photons/s  (exact 3.132e8)", a_s / 2.0);
    res.put("h_rayleigh500", ray(alpha));
    res.put("h_sigma_res", sres);
    res.put("h_rate_max", a_s / 2.0);
}

/// <1s|z|2p0> from two clouds of guided balls, 4 time blocks.
/// The ball density is ψ_g·ψ, and ψ/ψ_g depends on r only (the guide carries the exact angular factor,
/// fixed by symmetry), so each cloud measures the radial correction R(r) = ψ/ψ_g:
///   R(r) ∝ (balls in the shell) / ∫_shell ψ_g² dV,  ψ = ψ_g R.
/// Shells with fewer than 30 balls keep the last measured R (their weight in the integrals is < 1e-4).
fn elemento_z12(a: f64, a2: f64, t1: f64, t2: f64, seed: u64) -> (f64, f64) {
    use std::f64::consts::PI;
    let nblk = 4;
    let mut zs = vec![];
    let g1 = Guida::Orbitali { a, jastrow: false };
    let g2 = Guida::Eccitato { ang: Angolare::Z, nodi: [0.0; 2], nn: 0, a: a2 };
    let mut d1 = Dmc::new(idrogeno_sys(g1, vec![], 0.0), t1, 20000, seed);
    let mut d2 = Dmc::new(idrogeno_sys(g2, vec![Vincolo::PianoZ], 0.0), t2, 20000, seed + 1);
    corri(&mut d1, 10.0, 0.0, |_| {});
    corri(&mut d2, 60.0, 0.0, |_| {});
    let (dr, rmax) = (0.05, 40.0);
    // ∫ over a shell of f(r) dr by 20-point midpoint
    let shell = |r0: f64, f: &dyn Fn(f64) -> f64| (0..20).map(|i| f(r0 + (i as f64 + 0.5) * dr / 20.0)).sum::<f64>() * dr / 20.0;
    let radial = |h: &Istogramma, w: &dyn Fn(f64) -> f64| -> Vec<f64> {
        let mut last = 0.0;
        h.h.iter().enumerate().map(|(i, &c)| {
            if c >= 30.0 { last = c / h.n / shell(i as f64 * dr, w); }
            last
        }).collect()
    };
    for _ in 0..nblk {
        let mut h1 = Istogramma::new(dr, rmax);
        let mut h2 = Istogramma::new(dr, rmax);
        corri(&mut d1, 0.0, 40.0, |d| h1.add(d));
        corri(&mut d2, 0.0, 80.0, |d| h2.add(d));
        // shell integrals of ψ_g²: 1s over the sphere, 2p (r cosθ e^{-a2 r}) over the upper half space
        let w1 = |r: f64| 4.0 * PI * r * r * (-2.0 * a * r).exp();
        let w2 = |r: f64| 2.0 * PI / 3.0 * r.powi(4) * (-2.0 * a2 * r).exp();
        let r1 = radial(&h1, &w1);
        let r2 = radial(&h2, &w2);
        // norms over all space and the matrix element (angular integrals done exactly: ∫cos²θ dΩ = 4π/3)
        let (mut n1, mut n2, mut zz) = (0.0, 0.0, 0.0);
        for i in 0..r1.len() {
            let r0 = i as f64 * dr;
            n1 += r1[i] * r1[i] * shell(r0, &|r| 4.0 * PI * r * r * (-2.0 * a * r).exp());
            n2 += r2[i] * r2[i] * shell(r0, &|r| 4.0 * PI / 3.0 * r.powi(4) * (-2.0 * a2 * r).exp());
            zz += r1[i] * r2[i] * shell(r0, &|r| 4.0 * PI / 3.0 * r.powi(4) * (-(a + a2) * r).exp());
        }
        zs.push(zz / (n1 * n2).sqrt());
    }
    let zm = zs.iter().sum::<f64>() / nblk as f64;
    let dz = (zs.iter().map(|z| (z - zm).powi(2)).sum::<f64>() / (nblk * (nblk - 1)) as f64).sqrt();
    (zm, dz)
}

fn inv3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (a, b) = ((j + 1) % 3, (j + 2) % 3);
            let (c, d) = ((i + 1) % 3, (i + 2) % 3);
            r[i][j] = (m[a][c] * m[b][d] - m[a][d] * m[b][c]) / det;
        }
    }
    r
}

// ───────────────────────── 3. molecules: H2+ and H2 ─────────────────────────
/// Decay of our 1s orbital on each of two equal nuclei fixed by the electron–nucleus cusp (Kato):
/// at a nucleus ∂lnφ/∂r = −Z  ⇒  a = 1 + exp(−a Z R) for φ = e^{−aZ r_A} + e^{−aZ r_B}.
/// (For one nucleus this gives a = 1, which the unguided hydrogen cloud confirms: a = 0.9987.)
pub fn a_cuspide(z: f64, r: f64) -> f64 {
    let mut a = 1.0;
    for _ in 0..200 { a = 1.0 + (-a * z * r).exp(); }
    a
}

pub fn biatomica(nuclei_r: f64, z: f64, n_el: usize) -> Sistema {
    let a = a_cuspide(z, nuclei_r);
    Sistema {
        nuclei: vec![Nucleo { pos: [0.0, 0.0, -0.5 * nuclei_r], z }, Nucleo { pos: [0.0, 0.0, 0.5 * nuclei_r], z }],
        spin: if n_el == 1 { vec![1] } else { vec![1, -1] },
        guida: Guida::Orbitali { a, jastrow: n_el > 1 },
        vincoli: vec![],
        campo_z: 0.0,
    }
}

/// Fit V(R) = E∞ + Σ_{k=1..5} c_k e^{−kβR} (E∞ fixed: the separated atoms, computed by us), β chosen by χ².
pub struct Curva { pub e_inf: f64, pub beta: f64, pub c: Vec<f64> }

impl Curva {
    pub fn v(&self, r: f64) -> f64 { self.e_inf + self.c.iter().enumerate().map(|(i, c)| c * (-(i as f64 + 1.0) * self.beta * r).exp()).sum::<f64>() }
    pub fn fit(pts: &[(f64, f64, f64)], e_inf: f64) -> (Curva, f64) {
        let mut best: Option<(Curva, f64)> = None;
        for ib in 0..30 {
            let beta = 0.2 + 0.04 * ib as f64;
            let nb = 5;
            let mut ata = vec![vec![0.0; nb]; nb];
            let mut atb = vec![0.0; nb];
            for &(r, e, de) in pts {
                let row: Vec<f64> = (0..nb).map(|i| (-(i as f64 + 1.0) * beta * r).exp()).collect();
                let w = 1.0 / (de * de);
                for i in 0..nb { atb[i] += w * row[i] * (e - e_inf); for j in 0..nb { ata[i][j] += w * row[i] * row[j]; } }
            }
            let c = solve(ata, atb);
            let cv = Curva { e_inf, beta, c };
            let chi: f64 = pts.iter().map(|&(r, e, de)| ((cv.v(r) - e) / de).powi(2)).sum();
            if best.as_ref().is_none_or(|b| chi < b.1) { best = Some((cv, chi)); }
        }
        best.unwrap()
    }
    /// Minimum (Re, V(Re)) and curvature k = V''(Re).
    pub fn minimo(&self) -> (f64, f64, f64) {
        let mut r = 0.5;
        let mut best = (r, self.v(r));
        while r < 8.0 { let v = self.v(r); if v < best.1 { best = (r, v); } r += 1e-4; }
        let h = 1e-3;
        let re = best.0;
        let kk = (self.v(re + h) - 2.0 * self.v(re) + self.v(re - h)) / (h * h);
        (re, best.1, kk)
    }
}

fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for i in 0..n {
        let p = (i..n).max_by(|&x, &y| a[x][i].abs().total_cmp(&a[y][i].abs())).unwrap();
        a.swap(i, p);
        b.swap(i, p);
        for r in i + 1..n {
            let f = a[r][i] / a[i][i];
            for c in i..n { a[r][c] -= f * a[i][c]; }
            b[r] -= f * b[i];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() { x[i] = (b[i] - (i + 1..n).map(|c| a[i][c] * x[c]).sum::<f64>()) / a[i][i]; }
    x
}

/// Lowest `nlev` levels of −1/(2μ) d²/dR² + V(R) (J = 0) by finite differences + Sturm bisection.
pub fn livelli(v: impl Fn(f64) -> f64, mu: f64, r0: f64, r1: f64, npt: usize, nlev: usize) -> Vec<f64> {
    let h = (r1 - r0) / (npt + 1) as f64;
    let off = -0.5 / (mu * h * h);
    let d: Vec<f64> = (1..=npt).map(|i| v(r0 + i as f64 * h) + 1.0 / (mu * h * h)).collect();
    let count = |x: f64| { // eigenvalues below x
        let mut c = 0;
        let mut q = d[0] - x;
        if q < 0.0 { c += 1; }
        for &di in d.iter().skip(1) {
            q = di - x - off * off / if q.abs() < 1e-300 { 1e-300 } else { q };
            if q < 0.0 { c += 1; }
        }
        c
    };
    let lo0 = d.iter().cloned().fold(f64::INFINITY, f64::min) - 2.0 * off.abs();
    let hi0 = d.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 2.0 * off.abs();
    (0..nlev).map(|kk| {
        let (mut lo, mut hi) = (lo0, hi0);
        for _ in 0..200 { let m = 0.5 * (lo + hi); if count(m) > kk { hi = m } else { lo = m } }
        0.5 * (lo + hi)
    }).collect()
}

/// Properties of a diatomic curve: (Re Å, De eV, ωe cm⁻¹, ν01 cm⁻¹, ZPE cm⁻¹, D0 eV).
pub fn proprieta(c: &Curva, mu: f64) -> [f64; 6] {
    let (re, vmin, kk) = c.minimo();
    let we = (kk / mu).sqrt() * k::hartree_cm1();
    let lv = livelli(|r| c.v(r), mu, 0.4, 14.0, 6000, 2);
    [re * k::bohr_ang(), (c.e_inf - vmin) * k::hartree_ev(), we, (lv[1] - lv[0]) * k::hartree_cm1(),
     (lv[0] - vmin) * k::hartree_cm1(), (c.e_inf - lv[0]) * k::hartree_ev()]
}

pub fn molecola(res: &mut Risultati, n_el: usize) {
    let e_h = res.f("h_E_guided").expect("run `--misura idrogeno` first");
    let (name, key, grid, e_inf, mu, exp) = if n_el == 1 {
        ("H2+", "h2p_mol", vec![1.2, 1.5, 1.7, 1.85, 2.0, 2.15, 2.3, 2.6, 3.0, 3.5, 4.0, 5.0, 6.5, 8.0],
         e_h, { let mh = k::proton_me() + 1.0; mh * k::proton_me() / (mh + k::proton_me()) },
         [1.052, 2.6507, 2321.7, 2191.13])
    } else {
        ("H2", "h2_mol", vec![0.9, 1.1, 1.2, 1.25, 1.3, 1.33, 1.4, 1.47, 1.5, 1.55, 1.6, 1.7, 1.9, 2.2, 2.6, 3.2, 4.0, 5.5],
         2.0 * e_h, (k::proton_me() + 1.0) / 2.0, [0.74144, 4.47813, 4401.21, 4161.17])
    };
    println!("
=== 3. {name}: two protons, {n_el} electron(s); guide = our 1s orbital on each proton{} ===",
             if n_el > 1 { " × e–e cusp" } else { "" });
    println!("  (with the per-ball accept/reject of the move; points cached in the results file: rerun to continue)");
    let t0 = Instant::now();
    // tick 0.02 a.u.; checked at half the tick at the minimum
    let tau = 0.02;
    let t_run = if n_el == 1 { 300.0 } else { 200.0 };
    let mut pts = vec![];
    // near the minimum the curvature needs less noise: there each R gets extra independent runs
    // (different seeds), combined by weighted mean
    let extra = |r: f64| if n_el == 2 && (1.2..=1.6).contains(&r) { 3 } else { 1 };
    for (i, &r) in grid.iter().enumerate() {
        let mut runs = vec![];
        for rep in 0..extra(r) {
            let sfx = if rep == 0 { String::new() } else { format!("_run{rep}") };
            let (ke, kd) = (format!("{key}_E_R{r}{sfx}"), format!("{key}_Eerr_R{r}{sfx}"));
            if let (Some(e), Some(de)) = (res.f(&ke), res.f(&kd)) { runs.push((e, de)); continue; }
            if t0.elapsed().as_secs_f64() > 300.0 { println!("  budget reached: rerun to continue"); return; }
            let t = Instant::now();
            let mut d = Dmc::new(biatomica(r, 1.0, n_el), tau, 10000, 1000 * n_el as u64 + i as u64 + 100 * rep as u64);
            d.metropolis = true;
            let e = corri(&mut d, 10.0, t_run, |_| {});
            println!("  R = {r:4.2} bohr{sfx}: E = {:.6} ± {:.6} hartree  [{:.0} s]", e.e.0, e.e.1, t.elapsed().as_secs_f64());
            res.put(&ke, e.e.0);
            res.put(&kd, e.e.1);
            res.save(FILE);
            runs.push(e.e);
        }
        let w: f64 = runs.iter().map(|r| 1.0 / (r.1 * r.1)).sum();
        let e = runs.iter().map(|r| r.0 / (r.1 * r.1)).sum::<f64>() / w;
        pts.push((r, e, 1.0 / w.sqrt()));
    }
    // time-step check at the minimum: half the tick
    let rmin = if n_el == 1 { 2.0 } else { 1.4 };
    let kt = format!("{key}_E_halftau_R{rmin}");
    if res.f(&kt).is_none() {
        if t0.elapsed().as_secs_f64() > 200.0 { println!("  budget reached: rerun to continue"); return; }
        let mut d = Dmc::new(biatomica(rmin, 1.0, n_el), 0.5 * tau, 10000, 5000 + 10 * n_el as u64);
        d.metropolis = true;
        let e = corri(&mut d, 10.0, t_run, |_| {});
        res.put(&kt, e.e.0);
        res.put(&format!("{key}_Eerr_halftau_R{rmin}"), e.e.1);
        res.save(FILE);
    }
    let e01 = pts.iter().find(|p| p.0 == rmin).unwrap();
    println!("  time-step check at R = {rmin}: tau {tau} {:.6} ± {:.6}, tau {} {:.6} ± {:.6}",
             e01.1, e01.2, 0.5 * tau, res.f(&kt).unwrap(), res.f(&format!("{key}_Eerr_halftau_R{rmin}")).unwrap());
    res.put(&format!("{key}_tau"), tau);
    let (c, chi) = Curva::fit(&pts, e_inf);
    let p = proprieta(&c, mu);
    // bootstrap: refit with Gaussian noise of the size of each error bar
    let mut rng = Rng::new(77);
    let mut boot = vec![];
    // the point errors are scaled by √(χ²/point) when the fit says they are too small
    let scale = (chi / (pts.len() as f64 - 6.0)).sqrt().max(1.0);
    println!("  error bars of the points scaled by {scale:.2} in the re-fits");
    for _ in 0..200 {
        let q: Vec<(f64, f64, f64)> = pts.iter().map(|&(r, e, de)| (r, e + scale * de * rng.gauss(), de)).collect();
        boot.push(proprieta(&Curva::fit(&q, e_inf).0, mu));
    }
    // robust spread: half the 16–84 % interval (a few re-fits pick a bad β and land far away)
    let sd: Vec<f64> = (0..6).map(|j| {
        let mut v: Vec<f64> = boot.iter().map(|b| b[j]).filter(|x| x.is_finite()).collect();
        v.sort_by(|a, b| a.total_cmp(b));
        let q = |f: f64| v[((v.len() - 1) as f64 * f).round() as usize];
        0.5 * (q(0.84) - q(0.16))
    }).collect();
    println!("  fit: χ²/point = {:.2}, β = {:.2}", chi / pts.len() as f64, c.beta);
    // cross-check of the curvature: local quartic through the points within 0.35 bohr of the minimum
    let re0 = c.minimo().0;
    let loc: Vec<&(f64, f64, f64)> = pts.iter().filter(|p| (p.0 - re0).abs() <= 0.35).collect();
    if loc.len() >= 6 {
        let nb = 5;
        let mut ata = vec![vec![0.0; nb]; nb];
        let mut atb = vec![0.0; nb];
        for p in &loc {
            let row: Vec<f64> = (0..nb).map(|i| (p.0 - re0).powi(i as i32)).collect();
            let w = 1.0 / (p.2 * p.2);
            for i in 0..nb { atb[i] += w * row[i] * p.1; for j in 0..nb { ata[i][j] += w * row[i] * row[j]; } }
        }
        let cc = solve(ata, atb);
        let w_loc = (2.0 * cc[2] / mu).sqrt() * k::hartree_cm1();
        println!("  cross-check: local quartic on {} points near the minimum -> ωe = {w_loc:.1} cm⁻¹", loc.len());
        res.put(&format!("{key}_we_local"), w_loc);
    }
    let names = ["Re (Å)", "De (eV)", "ωe (cm⁻¹)", "ν(0→1) (cm⁻¹)", "ZPE (cm⁻¹)", "D0 (eV)"];
    let expv = [Some(exp[0]), None, Some(exp[2]), Some(exp[3]), None, Some(exp[1])];
    for j in 0..6 {
        let e = expv[j].map_or("".into(), |x| format!("   experiment {x}   error {:+.4}", p[j] - x));
        println!("  {:14} = {:.5} ± {:.5}{e}", names[j], p[j], sd[j]);
    }
    let ks = ["re", "de", "we", "nu01", "zpe", "d0"];
    for j in 0..6 {
        res.put(&format!("{key}_{}", ks[j]), p[j]);
        res.put(&format!("{key}_{}_err", ks[j]), sd[j]);
    }
    res.put(&format!("{key}_line_um"), 1e4 / p[3]);
    println!("  its vibrational photon: λ = 1/ν = {:.4} µm (infrared)", 1e4 / p[3]);
}

// ───────────────────────── 4. helium ─────────────────────────
pub fn elio(res: &mut Risultati) {
    println!("\n=== 4. Helium atom: 2 electrons (opposite spins: no sign problem, exact) ===");
    // guide: our 1s for Z = 2 (decay fixed by the cusp: a = 1 → e^{-2r}) for each electron × e–e cusp
    let taus = [0.02, 0.01, 0.005];
    let mut es = vec![];
    for (i, &tau) in taus.iter().enumerate() {
        let t = Instant::now();
        let sys = Sistema { nuclei: vec![Nucleo { pos: [0.0; 3], z: 2.0 }], spin: vec![1, -1],
                            guida: Guida::Orbitali { a: a_cuspide(2.0, 1e9), jastrow: true }, vincoli: vec![], campo_z: 0.0 };
        let mut d = Dmc::new(sys, tau, 10000, 900 + i as u64);
        d.metropolis = true;
        let e = corri(&mut d, 5.0, 150.0, |_| {});
        println!("  tau {tau:<6} E = {:.6} ± {:.6} hartree  [{:.0} s]", e.e.0, e.e.1, t.elapsed().as_secs_f64());
        es.push(e.e);
    }
    let (e0, de0) = estrapola(&taus, &es);
    let mu = k::alpha_me() / (k::alpha_me() + 1.0);
    println!("  tau -> 0: E = {e0:.6} ± {de0:.6} hartree = {:.4} ± {:.4} eV (fixed nucleus)", ev(e0), ev(de0));
    println!("  with the finite nuclear mass (μ/m_e, mass polarization neglected): {:.4} eV   experiment -79.0052 eV (NIST: 24.5874 + 54.4178)",
             ev(e0) * mu);
    res.put("he_E", e0);
    res.put("he_E_err", de0);
    res.put("he_E_ev_mu", ev(e0) * mu);
}

pub fn main(which: &str) {
    let t = Instant::now();
    let mut res = Risultati::load(FILE);
    println!("constants -> atomic units: a0 = {:.10} Å, E_h = {:.8} eV, m_p/m_e = {:.5}",
             k::bohr_ang(), k::hartree_ev(), k::proton_me());
    let all = which == "tutto";
    if all || which == "idrogeno" { idrogeno(&mut res); res.save(FILE); }
    if all || which == "eccitati" { eccitati(&mut res, "simm"); righe(&mut res); res.save(FILE); }
    if all || which == "nodo2s" { eccitati(&mut res, "2s"); righe(&mut res); res.save(FILE); }
    if all || which == "nodo3p" { eccitati(&mut res, "3p"); righe(&mut res); res.save(FILE); }
    if all || which == "campo" { luce(&mut res, "campo"); res.save(FILE); }
    if all || which == "riga" { luce(&mut res, "riga"); res.save(FILE); }
    if all || which == "h2piu" { molecola(&mut res, 1); res.save(FILE); }
    if all || which == "h2" { molecola(&mut res, 2); res.save(FILE); }
    if all || which == "elio" { elio(&mut res); res.save(FILE); }
    println!("\n[{which}: {:.0} s] -> {FILE}", t.elapsed().as_secs_f64());
}

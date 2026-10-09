//! HELIUM AND LIGHT: the first excited singlet 1s2p ¹P, the resonance line 1s² ¹S → 1s2p ¹P,
//! its strength (dipole matrix element from the two ball clouds), the static polarizability,
//! and the photon–atom cross-section used by the collision viewer (`scontro`).
//!
//! Inputs: the same as everywhere (ħ, e, mₑ, ε0, c, Z = 2, mass of the α particle). Orbitals: ONLY our
//! hydrogen-like ones — 1s = e^{−k_a r}, 2p_z = z e^{−k_b r}. k_a = Z = 2 is fixed by the electron–nucleus
//! cusp (Coulomb law, as for every guide of this level); k_b is chosen by minimizing OUR variational energy
//! of the 1s2p guide (`vmc_kb`). No external orbital, no Hartree–Fock.
//!
//! The excited state is a fixed-node calculation: the balls live in the pockets of the guide
//! ψ_g = φ_a(r1)φ_b(r2) + φ_b(r1)φ_a(r2) (× e–e cusp) and never cross its node (accept/reject rule
//! refuses such steps). The node is NOT exact: symmetry only says ψ(z1,z2) = −ψ(−z1,−z2), which fixes
//! ψ = 0 on the 4-D set z1 = z2 = 0, not on a 5-D surface. The node of our orbitals is an approximation,
//! and its error is checked against the exact (variational, Drake) energy of 2¹P in the doc.

use crate::costanti as k;
use crate::dmc::*;
use crate::misura::{corri, estrapola, Risultati};
use rayon::prelude::*;
use std::f64::consts::PI;
use std::time::Instant;

/// Exponent of our 1s for Z = 2, fixed by the electron–nucleus cusp: ∂lnφ/∂r = −Z.
pub const KA: f64 = 2.0;

/// Helium ground state 1s² ¹S: our 1s for Z = 2 on both electrons × e–e cusp (the guide of `misura::elio`).
pub fn sys_fondamentale(campo: f64) -> Sistema {
    Sistema { nuclei: vec![Nucleo { pos: [0.0; 3], z: 2.0 }], spin: vec![1, -1],
              guida: Guida::Orbitali { a: crate::misura::a_cuspide(2.0, 1e9), jastrow: true }, vincoli: vec![], campo_z: campo }
}

/// Helium ground state with the screened guide (used for the field runs: much less noise).
pub fn sys_schermato(zeta: f64, campo: f64) -> Sistema {
    Sistema { nuclei: vec![Nucleo { pos: [0.0; 3], z: 2.0 }], spin: vec![1, -1],
              guida: Guida::Schermato { zeta }, vincoli: vec![], campo_z: campo }
}

/// Variational energy ⟨ψ|H|ψ⟩/⟨ψ|ψ⟩ and the variance of E_loc of any two-electron guide (same sampling).
pub fn vmc(s: &Sistema, n: usize, seed: u64) -> (f64, f64) {
    let nchunk = 64;
    let per = n / nchunk;
    let parts: Vec<[f64; 3]> = (0..nchunk).into_par_iter().map(|c| {
        let mut rng = Rng::new(seed ^ (c as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let (mut x, mut d, mut a) = ([0.0; 6], [0.0; 6], [0.0; 3]);
        for _ in 0..per {
            let p = campiona(&mut rng, &mut x);
            let w = (2.0 * s.ln_guida(&x)).exp() / p;
            let e = s.valuta(&x, &mut d);
            a[0] += w;
            a[1] += w * e;
            a[2] += w * e * e;
        }
        a
    }).collect();
    let t = |i: usize| parts.iter().map(|a| a[i]).sum::<f64>();
    let e = t(1) / t(0);
    (e, t(2) / t(0) - e * e)
}

/// ζ of the screened guide: the minimum of our variational energy (scan + parabola).
pub fn vmc_zeta() -> f64 {
    let zs = [1.55, 1.65, 1.75, 1.85, 1.95];
    let es: Vec<f64> = zs.iter().map(|&z| vmc(&sys_schermato(z, 0.0), 2_000_000, 91).0).collect();
    let (mut ata, mut atb) = ([[0.0f64; 3]; 3], [0.0f64; 3]);
    for (&x, &y) in zs.iter().zip(&es) {
        let row = [1.0, x, x * x];
        for i in 0..3 { atb[i] += row[i] * y; for j in 0..3 { ata[i][j] += row[i] * row[j]; } }
    }
    let c = solve3(ata, atb);
    if c[2] > 0.0 { (-c[1] / (2.0 * c[2])).clamp(zs[0], zs[4]) } else { 1.75 }
}

/// Helium 1s2p ¹P (M = 0) with the guide built from our 1s (k_a = 2) and 2p_z (k_b).
pub fn sys_1p(kb: f64) -> Sistema {
    Sistema { nuclei: vec![Nucleo { pos: [0.0; 3], z: 2.0 }], spin: vec![1, -1],
              guida: Guida::S1P { ka: KA, kb }, vincoli: vec![], campo_z: 0.0 }
}

/// Signed guide value (unnormalised) ψ_g(x) = sign · e^{ln|ψ_g|}.
pub fn psi(s: &Sistema, x: &[f64]) -> f64 { s.segno(x) * s.ln_guida(x).exp() }

// ───────────────────────── integrals of the guides (plain Monte Carlo, importance sampled) ─────────────────────────
/// Integrals of the two guides over the 6-D space, by Monte Carlo with a known normalised density p
/// (each electron: ½ of e^{−4r} + ½ of e^{−0.8r}, uniform direction — broader than both guides, so every
/// estimator has finite variance).
#[derive(Clone, Copy, Debug, Default)]
pub struct Integrali {
    /// ∫ψ_g² (ground guide), ∫ψ_e² (1s2p guide)
    pub ig: f64,
    pub ie: f64,
    /// ⟨ψ_g|z1+z2|ψ_e⟩/√(ig·ie): the dipole matrix element of the two guides (variational estimate)
    pub dvar: f64,
    /// variational energy of the 1s2p guide ⟨ψ_e|H|ψ_e⟩/⟨ψ_e|ψ_e⟩ and of the ground guide
    pub evar_e: f64,
    pub evar_g: f64,
    /// statistical errors of dvar and evar_e
    pub dvar_err: f64,
    pub evar_e_err: f64,
}

const BETAS: [f64; 2] = [4.0, 0.8];

fn campiona(rng: &mut Rng, out: &mut [f64]) -> f64 {
    let mut p = 1.0;
    for e in 0..2 {
        let b = if rng.uni() < 0.5 { BETAS[0] } else { BETAS[1] };
        let r = -((1.0 - rng.uni()) * (1.0 - rng.uni()) * (1.0 - rng.uni())).ln() / b;
        let ct = 2.0 * rng.uni() - 1.0;
        let ph = 2.0 * PI * rng.uni();
        let st = (1.0 - ct * ct).sqrt();
        out[3 * e] = r * st * ph.cos();
        out[3 * e + 1] = r * st * ph.sin();
        out[3 * e + 2] = r * ct;
        let pe: f64 = BETAS.iter().map(|b| 0.5 * b.powi(3) / (8.0 * PI) * (-b * r).exp()).sum();
        p *= pe;
    }
    p
}

/// The integrals with `n` samples (deterministic for a given seed).
pub fn integrali(kb: f64, n: usize, seed: u64) -> Integrali {
    let sg = sys_fondamentale(0.0);
    let se = sys_1p(kb);
    let nchunk = 64;
    let per = n / nchunk;
    // per chunk: Σ wg, Σ we, Σ cross, Σ we·E_e, Σ wg·E_g
    let parts: Vec<[f64; 5]> = (0..nchunk).into_par_iter().map(|c| {
        let mut rng = Rng::new(seed ^ (c as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let mut x = [0.0; 6];
        let mut d = [0.0; 6];
        let mut s = [0.0; 5];
        for _ in 0..per {
            let p = campiona(&mut rng, &mut x);
            let (g, e) = (psi(&sg, &x), psi(&se, &x));
            s[0] += g * g / p;
            s[1] += e * e / p;
            s[2] += g * (x[2] + x[5]) * e / p;
            if e != 0.0 { s[3] += e * e / p * se.valuta(&x, &mut d); }
            s[4] += g * g / p * sg.valuta(&x, &mut d);
        }
        s
    }).collect();
    let tot = |i: usize| parts.iter().map(|s| s[i]).sum::<f64>() / (per * nchunk) as f64;
    let (ig, ie) = (tot(0), tot(1));
    // errors from the spread of the chunk estimates
    let spread = |f: &dyn Fn(&[f64; 5]) -> f64| {
        let v: Vec<f64> = parts.iter().map(f).collect();
        let m = v.iter().sum::<f64>() / v.len() as f64;
        (v.iter().map(|y| (y - m).powi(2)).sum::<f64>() / (v.len() * (v.len() - 1)) as f64).sqrt()
    };
    Integrali {
        ig, ie,
        dvar: tot(2) / (ig * ie).sqrt(),
        evar_e: tot(3) / ie,
        evar_g: tot(4) / ig,
        dvar_err: spread(&|s| s[2] / (s[0] * s[1]).sqrt()),
        evar_e_err: spread(&|s| s[3] / s[1]),
    }
}

/// k_b of our 2p_z: the minimum of the variational energy of the 1s2p guide (parabola through a scan).
pub fn vmc_kb() -> (f64, Vec<(f64, f64, f64)>) {
    let kbs = [0.40, 0.45, 0.50, 0.55, 0.60];
    let scan: Vec<(f64, f64, f64)> = kbs.iter().map(|&kb| {
        let i = integrali(kb, 2_000_000, 77);
        (kb, i.evar_e, i.evar_e_err)
    }).collect();
    // least-squares parabola through the scan
    let (mut ata, mut atb) = ([[0.0f64; 3]; 3], [0.0f64; 3]);
    for &(x, y, _) in &scan {
        let row = [1.0, x, x * x];
        for i in 0..3 { atb[i] += row[i] * y; for j in 0..3 { ata[i][j] += row[i] * row[j]; } }
    }
    let c = solve3(ata, atb);
    let kb = if c[2] > 0.0 { (-c[1] / (2.0 * c[2])).clamp(kbs[0], kbs[kbs.len() - 1]) } else { 0.5 };
    (kb, scan)
}

fn solve3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> [f64; 3] {
    for i in 0..3 {
        let p = (i..3).max_by(|&x, &y| a[x][i].abs().partial_cmp(&a[y][i].abs()).unwrap()).unwrap();
        a.swap(i, p);
        b.swap(i, p);
        for j in i + 1..3 {
            let f = a[j][i] / a[i][i];
            for kk in i..3 { a[j][kk] -= f * a[i][kk]; }
            b[j] -= f * b[i];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        x[i] = (b[i] - (i + 1..3).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}

// ───────────────────────── the light of helium: cross-section per photon ─────────────────────────
/// Everything the collision needs, computed by the balls (read from assets/molecole/risultati.txt).
/// (Used by the `scontro` binary and the tests; the `molecole` viewer only shows the stored numbers.)
#[allow(dead_code)]
#[derive(Clone, Copy, Debug)]
pub struct LuceElio {
    /// transition energy (hartree, reduced mass) and its wavelength (nm, vacuum)
    pub de: f64,
    pub lambda_nm: f64,
    /// oscillator strength, Einstein A (1/s), static polarizability (a0³)
    pub f: f64,
    pub a_s: f64,
    pub alpha0: f64,
    /// energies (hartree, fixed nucleus) of 1¹S and 2¹P, k_b of our 2p, norm ratio ∫ψ_e²/∫ψ_g² of the guides
    pub e_g: f64,
    pub e_e: f64,
    pub kb: f64,
    pub norma_e_g: f64,
}

#[allow(dead_code)]
impl LuceElio {
    pub fn da(res: &Risultati) -> Option<Self> {
        Some(LuceElio {
            de: res.f("he1p_de")?, lambda_nm: res.f("he1p_lambda_nm")?, f: res.f("he1p_f")?, a_s: res.f("he1p_A")?,
            alpha0: res.f("he_alpha")?, e_g: res.f("he1p_Eg_comb")?, e_e: res.f("he1p_E")?, kb: res.f("he1p_kb")?,
            norma_e_g: res.f("he1p_norma_e_g")?,
        })
    }
    /// Angular frequency (rad/s) of a photon of wavelength λ (nm).
    pub fn omega_s(lambda_nm: f64) -> f64 { 2.0 * PI * k::C_L / (lambda_nm * 1e-9) }
    /// Line centre (rad/s).
    pub fn omega0_s(&self) -> f64 { Self::omega_s(self.lambda_nm) }
    /// Complex polarizability (a0³) at angular frequency ω (rad/s): the computed line as a damped oscillator
    /// with strength f, frequency ΔE and width Γ = A (all from the balls), plus the rest of the static
    /// polarizability (α0 − f/ΔE², i.e. all other transitions) taken as frequency independent.
    pub fn alpha(&self, w_s: f64) -> (f64, f64) {
        let w = w_s * k::au_time_s(); // a.u.
        let g = self.a_s * k::au_time_s();
        let w0 = self.de;
        let rest = self.alpha0 - self.f / (w0 * w0);
        // f / (w0² − w² − i w Γ)
        let (re, im) = (w0 * w0 - w * w, -w * g);
        let den = re * re + im * im;
        (rest + self.f * re / den, -self.f * im / den)
    }
    /// Cross-section (m²) for scattering a photon of wavelength λ: σ = (8π/3) k⁴ |α(ω)|².
    /// At the line centre this is exactly the resonant cross-section 3λ²/2π (because Γ = A = 2ω²f/3c³ comes
    /// from the same f and ΔE); far from it, it is Rayleigh scattering with the (dynamic) polarizability.
    pub fn sigma(&self, lambda_nm: f64) -> f64 {
        let a0 = k::bohr_m();
        let kk = 2.0 * PI / (lambda_nm * 1e-9);
        let (re, im) = self.alpha(Self::omega_s(lambda_nm));
        8.0 * PI / 3.0 * kk.powi(4) * (re * re + im * im) * a0.powi(6)
    }
    /// Probability that ONE photon, a Gaussian beam focused on the atom with waist w0 = λ (a tight but
    /// paraxial focus), is taken by the atom: σ / A_eff with A_eff = π w0²/2 (the peak fluence of one
    /// photon is 1/A_eff). At resonance: (3λ²/2π)/(πλ²/2) = 3/π² = 0.304.
    pub fn prob(&self, lambda_nm: f64) -> f64 {
        let w0 = lambda_nm * 1e-9;
        (self.sigma(lambda_nm) / (PI * w0 * w0 / 2.0)).min(1.0)
    }
    /// Is the photon within the resonance (|ω − ω0| < 50 Γ)? Then a taken photon is a real absorption
    /// (the atom goes to 1s2p and later re-emits); otherwise it is Rayleigh scattering (no real excitation).
    pub fn risonante(&self, lambda_nm: f64) -> bool {
        (Self::omega_s(lambda_nm) - self.omega0_s()).abs() < 50.0 * self.a_s
    }
    pub fn vita_s(&self) -> f64 { 1.0 / self.a_s }
}

// ───────────────────────── the measurements ─────────────────────────
fn ev(e: f64) -> f64 { e * k::hartree_ev() }

/// Taus of the excited-state and matrix-element runs.
const TAUS: [f64; 3] = [0.02, 0.01, 0.005];
const N_BALLS: usize = 10000;
/// independent ground-state runs per tick
const NREP_G: usize = 3;

/// One DMC run of the ground or the excited state, collecting the mixed estimate of the dipole matrix
/// element. On the ground balls (density ψ_g·Φ_g):  M_g = ⟨Z ψ_e/ψ_g⟩ · √(I_g/I_e) = ⟨Φ_g|Z|ψ_e⟩/⟨Φ_g|ψ_g⟩;
/// on the excited balls: M_e = ⟨Z ψ_g/ψ_e⟩ · √(I_e/I_g). Z = z1 + z2. Returns (E, dE, M, dM).
fn corsa(eccitato: bool, kb: f64, tau: f64, ints: &Integrali, seed: u64) -> (f64, f64, f64, f64) {
    let sg = sys_fondamentale(0.0);
    let se = sys_1p(kb);
    let sys = if eccitato { se.clone() } else { sg.clone() };
    let mut d = Dmc::new(sys, tau, N_BALLS, seed);
    d.metropolis = true;
    let norm = if eccitato { (ints.ie / ints.ig).sqrt() } else { (ints.ig / ints.ie).sqrt() };
    let mut ms = vec![];
    let t_eq = if eccitato { 25.0 } else { 5.0 };
    let r = corri(&mut d, t_eq, 120.0, |d| {
        let m = 6;
        let (mut s, mut n) = (0.0, 0.0);
        for w in d.x.chunks(m) {
            let z = w[2] + w[5];
            let (lg, le) = (sg.ln_guida(w), se.ln_guida(w));
            let ratio = if eccitato { (lg - le).exp() * se.segno(w) } else { (le - lg).exp() * se.segno(w) };
            s += z * ratio;
            n += 1.0;
        }
        ms.push(s / n * norm);
    });
    let (m, dm) = media_errore(&ms);
    (r.e.0, r.e.1, m, dm)
}

/// `--misura elio1p`: k_b by VMC, the 1s2p ¹P level (fixed node), the dipole matrix element, f, A, lifetime.
/// Resumable: every (state, τ) run is stored as soon as it ends and skipped next time (same k_b).
pub fn elio1p(res: &mut Risultati) {
    println!("\n=== 5. Helium 1s2p 1P: fixed node from our orbitals, resonance line, its strength ===");
    let t0 = Instant::now();
    let kb = match res.f("he1p_kb") {
        Some(kb) => { println!("  k_b of our 2p (stored): {kb:.4}"); kb }
        None => {
            let (kb, scan) = vmc_kb();
            for (x, e, de) in &scan { println!("  VMC of the 1s2p guide: k_b = {x:.2}  E = {e:.5} ± {de:.5} hartree"); }
            println!("  minimum of our variational energy: k_b = {kb:.4} (effective charge of the 2p electron 2·k_b = {:.3})", 2.0 * kb);
            res.put("he1p_kb", kb);
            for (x, e, _) in &scan { res.put(&format!("he1p_vmc_E_kb{x:.2}"), *e); }
            res.save(crate::misura::FILE);
            kb
        }
    };
    let ints = integrali(kb, 8_000_000, 4242);
    println!("  guides: <E>_var(1s2p) = {:.5} ± {:.5}, <E>_var(1s2) = {:.5}; d_var = <g|z1+z2|e> = {:.5} ± {:.5} a0  [{:.0} s]",
             ints.evar_e, ints.evar_e_err, ints.evar_g, ints.dvar, ints.dvar_err, t0.elapsed().as_secs_f64());
    res.put("he1p_Evar", ints.evar_e);
    res.put("he_Evar", ints.evar_g);
    res.put("he1p_dvar", ints.dvar);
    res.put("he1p_norma_e_g", ints.ie / ints.ig);
    // the ground-state mixed estimate M_g is the noisiest number (ψ_e/ψ_g grows at large r): 3 independent
    // runs per tick for the ground state, 1 for the excited state
    for rep in 0..NREP_G {
    for (i, &tau) in TAUS.iter().enumerate() {
        for ecc in [false, true] {
            if ecc && rep > 0 { continue; }
            let tag = if ecc { "he1p".to_string() } else if rep == 0 { "heg".to_string() } else { format!("heg_r{rep}") };
            let key = format!("{tag}_E_tau{tau}");
            if res.f(&key).is_some() { continue; }
            if t0.elapsed().as_secs_f64() > 420.0 {
                println!("  (time budget reached: run `--misura elio1p` again to continue)");
                res.save(crate::misura::FILE);
                return;
            }
            let t = Instant::now();
            let (e, de, m, dm) = corsa(ecc, kb, tau, &ints, 1300 + 10 * i as u64 + ecc as u64 + 100 * rep as u64);
            println!("  {} tau {tau:<6} E = {e:.6} ± {de:.6}   M_{} = {m:.5} ± {dm:.5}  [{:.0} s]",
                     if ecc { "1s2p" } else { "1s2 " }, if ecc { "e" } else { "g" }, t.elapsed().as_secs_f64());
            res.put(&key, e);
            res.put(&format!("{tag}_Eerr_tau{tau}"), de);
            res.put(&format!("{tag}_M_tau{tau}"), m);
            res.put(&format!("{tag}_Merr_tau{tau}"), dm);
            res.save(crate::misura::FILE);
        }
    }
    }
    riassunto_1p(res, &ints);
}

fn riassunto_1p(res: &mut Risultati, ints: &Integrali) {
    // weighted mean over the independent runs of one quantity (value key q, error key qerr)
    let g = |tag: &str, q: &str, tau: f64| -> f64 {
        let (val, err) = if q.ends_with("err") { (q.trim_end_matches("err"), q) } else { (q, "") };
        let qe = if err.is_empty() { format!("{q}err") } else { err.to_string() };
        let tags: Vec<String> = if tag == "heg" { (0..NREP_G).map(|r| if r == 0 { "heg".into() } else { format!("heg_r{r}") }).collect() }
                                else { vec![tag.to_string()] };
        let (mut sw, mut swx) = (0.0, 0.0);
        for t in &tags {
            if let (Some(v), Some(e)) = (res.f(&format!("{t}_{val}_tau{tau}")), res.f(&format!("{t}_{qe}_tau{tau}"))) {
                sw += 1.0 / (e * e);
                swx += v / (e * e);
            }
        }
        if err.is_empty() { swx / sw } else { (1.0 / sw).sqrt() }
    };
    let es: Vec<(f64, f64)> = TAUS.iter().map(|&t| (g("he1p", "E", t), g("he1p", "Eerr", t))).collect();
    let (e1p, de1p) = estrapola(&TAUS, &es);
    let eg_here: Vec<(f64, f64)> = TAUS.iter().map(|&t| (g("heg", "E", t), g("heg", "Eerr", t))).collect();
    let (eg2, deg2) = estrapola(&TAUS, &eg_here);
    // ground state: combine with the earlier independent `--misura elio` runs (same rule, same guide)
    let (eg, deg) = match (res.f("he_E"), res.f("he_E_err")) {
        (Some(a), Some(da)) => {
            let (wa, wb) = (1.0 / (da * da), 1.0 / (deg2 * deg2));
            ((a * wa + eg2 * wb) / (wa + wb), (1.0 / (wa + wb)).sqrt())
        }
        _ => (eg2, deg2),
    };
    println!("  1s2p 1P, tau -> 0: E = {e1p:.6} ± {de1p:.6} hartree  (exact nonrelativistic -2.123843; fixed-node error {:+.6})",
             e1p + 2.123_843_086);
    println!("  1s2 1S, these runs tau -> 0: {eg2:.6} ± {deg2:.6}; combined with `--misura elio`: {eg:.6} ± {deg:.6} (exact -2.903724)");
    let mu = k::alpha_me() / (k::alpha_me() + 1.0);
    let de = (e1p - eg) * mu;
    let dde = (de1p * de1p + deg * deg).sqrt() * mu;
    let lam = k::nm_from_hartree(de);
    let dlam = lam * dde / de;
    println!("  resonance line: ΔE = {:.4} ± {:.4} eV, λ = {lam:.3} ± {dlam:.3} nm   (NIST 21.2180 eV, 58.4334 nm)",
             ev(de), ev(dde));
    // matrix element: extrapolated estimator M_g + M_e − d_var at each tau, then tau -> 0 (straight line)
    let ds: Vec<(f64, f64)> = TAUS.iter().map(|&t| {
        let v = g("heg", "M", t) + g("he1p", "M", t) - ints.dvar;
        let e = (g("heg", "Merr", t).powi(2) + g("he1p", "Merr", t).powi(2) + ints.dvar_err.powi(2)).sqrt();
        println!("    tau {t:<6} M_g = {:.5}  M_e = {:.5}  d_var = {:.5}  ->  d = M_g + M_e - d_var = {v:.5} ± {e:.5}",
                 g("heg", "M", t), g("he1p", "M", t), ints.dvar);
        (v, e)
    }).collect();
    let (dz, ddz) = estrapola(&TAUS, &ds);
    let dmean = ds.iter().map(|x| x.0).sum::<f64>() / ds.len() as f64;
    let ddz = ddz.max((ds[2].0 - dz).abs());
    let d_mu = dz / mu;
    let fosc = 2.0 * de * d_mu * d_mu;
    let c_au = 1.0 / k::alpha_fs();
    let a_au = 4.0 * de.powi(3) * d_mu * d_mu / (3.0 * c_au.powi(3));
    let a_s = a_au / k::au_time_s();
    let rel = 2.0 * ddz / dz;
    println!("  <1S|z1+z2|1P0> = {dz:.4} ± {ddz:.4} a0 (tau -> 0; mean over taus {dmean:.4}; guides alone {:.4})", ints.dvar);
    println!("  f = {fosc:.4} ± {:.4} (exp./theory 0.2762), A = {a_s:.4e} ± {:.1e} 1/s (NIST 1.799e9), lifetime {:.4} ns (0.555)",
             fosc * rel, a_s * rel, 1e9 / a_s);
    res.put("he1p_E", e1p);
    res.put("he1p_E_err", de1p);
    res.put("he1p_Eg_comb", eg);
    res.put("he1p_Eg_comb_err", deg);
    res.put("he1p_de", de);
    res.put("he1p_de_err", dde);
    res.put("he1p_de_ev", ev(de));
    res.put("he1p_lambda_nm", lam);
    res.put("he1p_lambda_err", dlam);
    res.put("he1p_d", dz);
    res.put("he1p_d_err", ddz);
    res.put("he1p_f", fosc);
    res.put("he1p_f_err", fosc * rel);
    res.put("he1p_A", a_s);
    res.put("he1p_A_err", a_s * rel);
    res.put("he1p_tau_ns", 1e9 / a_s);
}

/// `--misura eliocampo`: static polarizability of helium, E(F) = E(0) − ½αF² − γF⁴/24 (fit, F → 0).
/// Resumable per field value.
pub fn eliocampo(res: &mut Risultati) {
    println!("\n=== 6. Helium in a uniform electric field: static polarizability ===");
    let t0 = Instant::now();
    let fs = [0.0, 0.05, 0.08, 0.11, 0.14];
    let tau = 0.01;
    // the guide: our 1s with an effective charge ζ (cusp kept exact), ζ from the minimum of our VMC energy
    let zeta = match res.f("he_zeta") {
        Some(z) => z,
        None => {
            let z = vmc_zeta();
            let (e1, v1) = vmc(&sys_fondamentale(0.0), 4_000_000, 5);
            let (e2, v2) = vmc(&sys_schermato(z, 0.0), 4_000_000, 5);
            println!("  screened guide: ζ = {z:.4}; VMC E = {e2:.5} (variance of E_loc {v2:.4}) vs the unscreened guide {e1:.5} ({v1:.4})");
            res.put("he_zeta", z);
            res.put("he_zeta_Evar", e2);
            res.put("he_zeta_var", v2);
            res.put("he_Evar_var", v1);
            z
        }
    };
    for (i, &f) in fs.iter().enumerate() {
        let key = format!("he_campof_E_F{f}");
        if res.f(&key).is_some() { continue; }
        if t0.elapsed().as_secs_f64() > 450.0 {
            println!("  (time budget reached: run `--misura eliocampo` again to continue)");
            return;
        }
        let t = Instant::now();
        let mut d = Dmc::new(sys_schermato(zeta, f), tau, 20000, 1500 + i as u64);
        d.metropolis = true;
        let mut zs = vec![];
        let r = corri(&mut d, 5.0, 80.0, |d| {
            let (s, n) = d.x.chunks(6).fold((0.0, 0.0), |a, w| (a.0 + w[2] + w[5], a.1 + 1.0));
            zs.push(s / n);
        });
        let (zm, dzm) = media_errore(&zs);
        println!("  F = {f:.2} a.u.: E = {:.6} ± {:.6}   <z1+z2>_balls = {zm:+.5} ± {dzm:.5}  [{:.0} s]",
                 r.e.0, r.e.1, t.elapsed().as_secs_f64());
        res.put(&format!("he_campof_Z_F{f}"), zm);
        res.put(&key, r.e.0);
        res.put(&format!("he_campof_Eerr_F{f}"), r.e.1);
        res.save(crate::misura::FILE);
    }
    let pts: Vec<(f64, f64, f64)> = fs.iter().map(|&f| (f, res.f(&format!("he_campof_E_F{f}")).unwrap(),
                                                         res.f(&format!("he_campof_Eerr_F{f}")).unwrap())).collect();
    let (mut ata, mut atb) = ([[0.0f64; 3]; 3], [0.0f64; 3]);
    for &(f, e, de) in &pts {
        let row = [1.0, f * f, f.powi(4)];
        let w = 1.0 / (de * de).max(1e-20);
        for i in 0..3 { atb[i] += w * row[i] * e; for j in 0..3 { ata[i][j] += w * row[i] * row[j]; } }
    }
    let c = solve3(ata, atb);
    // error of α: from the inverse of the normal matrix (column 1)
    let mut e1 = [0.0; 3];
    e1[1] = 1.0;
    let col = solve3(ata, e1);
    let alpha = -2.0 * c[1];
    let dalpha = 2.0 * col[1].sqrt();
    let chi2: f64 = pts.iter().map(|&(f, e, de)| ((c[0] + c[1] * f * f + c[2] * f.powi(4) - e) / de).powi(2)).sum();
    println!("  α = {alpha:.3} ± {dalpha:.3} a0³   (exact 1.3832; χ² = {chi2:.1} for 2 degrees of freedom); γ/24 term {:.1}",
             c[2]);
    // check: the induced dipole of the cloud, extrapolated estimate 2·<Z>_mixed − <Z>_guide (the guide is
    // symmetric: <Z>_guide = 0)  ->  α = −2<Z>_mixed / F
    let adip: Vec<f64> = fs.iter().skip(1).map(|&f| -2.0 * res.f(&format!("he_campof_Z_F{f}")).unwrap_or(0.0) / f).collect();
    println!("  check, α from the shift of the cloud (−2<z1+z2>/F) at F = {:?}: {:?}", &fs[1..],
             adip.iter().map(|x| (x * 1000.0).round() / 1000.0).collect::<Vec<_>>());
    res.put("he_alpha_dip", adip[0]);
    res.put("he_alpha", alpha);
    res.put("he_alpha_err", dalpha);
}

//! Hydrogen from the zig-zag rule (radial form, see lab/zz_hydrogen.py and docs/ZIGZAG.md).
//!
//! For every Dirac partial wave κ the electron walker (flip θ = a/λ_C) runs in the Coulomb phase of a
//! point proton; the stationary states are the frequencies at which the walk rings, found by "infinite
//! listening" (inverse iteration on the one-tick walk). Three lattices θ = 0.2, 0.1, 0.05 are extrapolated
//! to a → 0 (the rule is second order in a). The 3-D cloud of the viewer is built from the θ = 0.05
//! radial functions (g, f) and the spinor spherical harmonics.

use crate::cplx::C;
use crate::walk::Radial;

pub const ALPHA: f64 = 7.297_352_5693e-3;
pub const MC2_EV: f64 = 510_998.950_00;
pub const H_EVS: f64 = 4.135_667_696e-15;
pub const HC_EVNM: f64 = 1239.841_984_33;
pub const MP_ME: f64 = 1836.152_673_43;
pub const LAMBDA_C_M: f64 = 3.861_592_6796e-13;
pub const T_UNIT_S: f64 = 1.288_088_668_19e-21;
pub const C_MS: f64 = 299_792_458.0;
pub const THETAS: [f64; 3] = [0.2, 0.1, 0.05];

/// The levels computed by the viewer (name, n, κ).
pub const LEVELS: [(&str, u32, i32); 9] = [
    ("1s1/2", 1, -1), ("2s1/2", 2, -1), ("2p1/2", 2, 1), ("2p3/2", 2, -2), ("3s1/2", 3, -1),
    ("3p1/2", 3, 1), ("3p3/2", 3, -2), ("3d3/2", 3, 2), ("3d5/2", 3, -3),
];

/// Sommerfeld–Dirac energy (units m c², rest energy included). CHECK ONLY, never an input.
pub fn dirac_formula(n: u32, kappa: i32, za: f64) -> f64 {
    let k = kappa.unsigned_abs() as f64;
    let d = k - (k * k - za * za).sqrt();
    1.0 / (1.0 + (za / (n as f64 - d)).powi(2)).sqrt()
}

pub fn l_of(kappa: i32) -> u32 { if kappa > 0 { kappa as u32 } else { (-kappa - 1) as u32 } }

/// One level on one lattice: binding energy (m c²) and radial functions (g, f) on r_i = (i+½)a.
pub struct Level {
    pub binding: f64,
    pub a: f64,
    pub gf: Vec<[f64; 2]>,
}

/// Number of sign changes of g (ignoring the tiny tail).
pub fn nodes(gf: &[[f64; 2]]) -> usize {
    let gmax = gf.iter().map(|v| v[0].abs()).fold(0.0, f64::max);
    let mut last = 0.0f64;
    let mut n = 0;
    for v in gf {
        if v[0].abs() < 1e-6 * gmax { continue; }
        if last != 0.0 && v[0].signum() != last { n += 1; }
        last = v[0].signum();
    }
    n
}

/// Listen to the walk of partial wave κ for level n on lattice θ. The guess only chooses which ring to
/// listen to (the Bohr scale α²/2n² m c², from the inputs); the energy is the walk's own frequency.
pub fn level(n: u32, kappa: i32, theta: f64) -> Level {
    let rmax = (15.0 + 8.0 * (n * n) as f64) / ALPHA;
    let w = Radial::new(kappa, theta, rmax, ALPHA);
    let guess = 1.0 - 0.5 * ALPHA * ALPHA / (n * n) as f64 * 1.003;
    let (om, gf) = w.eigen(guess, 6);
    Level { binding: 1.0 - om, a: theta, gf }
}

/// Fit y = y0 + c1 θ² + c2 θ⁴ through three lattices → y0.
pub fn extrapolate(y: [f64; 3]) -> f64 {
    let t: Vec<f64> = THETAS.iter().map(|t| t * t).collect();
    // Lagrange interpolation in s = θ² evaluated at s = 0
    let mut y0 = 0.0;
    for i in 0..3 {
        let mut l = 1.0;
        for j in 0..3 { if j != i { l *= (0.0 - t[j]) / (t[i] - t[j]); } }
        y0 += y[i] * l;
    }
    y0
}

/// Everything the panel shows about hydrogen.
pub struct Hydrogen {
    /// binding energies (eV) per level per lattice, and extrapolated
    pub e_lat: Vec<[f64; 3]>,
    pub e_ext: Vec<f64>,
    pub e_dirac: Vec<f64>,
    /// radial functions at θ = 0.05 for 1s, 2s, 2p½, 2p3/2
    pub states: Vec<Level>,
    pub nodes_ok: bool,
    pub seconds: f64,
}

impl Hydrogen {
    pub fn compute() -> Self {
        let t0 = std::time::Instant::now();
        let mut e_lat = vec![];
        let mut e_ext = vec![];
        let mut e_dirac = vec![];
        let mut states = vec![];
        let mut nodes_ok = true;
        for (k, &(_, n, kappa)) in LEVELS.iter().enumerate() {
            let mut y = [0.0; 3];
            for (i, &th) in THETAS.iter().enumerate() {
                let lv = level(n, kappa, th);
                y[i] = lv.binding * MC2_EV;
                if nodes(&lv.gf) as u32 != n - l_of(kappa) - 1 { nodes_ok = false; }
                if i == 2 && k < 4 { states.push(lv); }
            }
            e_lat.push(y);
            e_ext.push(extrapolate(y));
            e_dirac.push((1.0 - dirac_formula(n, kappa, ALPHA)) * MC2_EV);
        }
        Self { e_lat, e_ext, e_dirac, states, nodes_ok, seconds: t0.elapsed().as_secs_f64() }
    }

    pub fn mu() -> f64 { 1.0 / (1.0 + 1.0 / MP_ME) }

    pub fn idx(name: &str) -> usize { LEVELS.iter().position(|l| l.0 == name).unwrap() }

    /// 2p3/2 − 2p1/2 (Hz), reduced mass.
    pub fn fine_structure_hz(&self) -> f64 {
        (self.e_ext[Self::idx("2p1/2")] - self.e_ext[Self::idx("2p3/2")]) * Self::mu() / H_EVS
    }

    /// Lyman-α centroid (2p₃/₂ : 2p₁/₂ = 2 : 1), nm.
    pub fn lyman_alpha_nm(&self) -> f64 {
        let e1 = self.e_ext[0];
        let e = (2.0 * (e1 - self.e_ext[Self::idx("2p3/2")]) + (e1 - self.e_ext[Self::idx("2p1/2")])) / 3.0 * Self::mu();
        HC_EVNM / e
    }

    /// Balmer-α centroid, nm: E1 components weighted by statistical population × Einstein A of the nl
    /// multiplet (hydrogen values 2.245e7, 6.31e6, 6.47e7 s⁻¹ — the lab computes them from the walk) ×
    /// relative line strength inside the multiplet.
    pub fn balmer_alpha_nm(&self) -> f64 {
        let e = |a: &str, b: &str| self.e_ext[Self::idx(b)] - self.e_ext[Self::idx(a)];
        let comps = [
            ("3p3/2", "2s1/2", 6.0 * 2.245e7 * 2.0 / 3.0), ("3p1/2", "2s1/2", 6.0 * 2.245e7 / 3.0),
            ("3s1/2", "2p3/2", 2.0 * 6.31e6 * 2.0 / 3.0), ("3s1/2", "2p1/2", 2.0 * 6.31e6 / 3.0),
            ("3d5/2", "2p3/2", 10.0 * 6.47e7 * 9.0 / 15.0), ("3d3/2", "2p3/2", 10.0 * 6.47e7 / 15.0),
            ("3d3/2", "2p1/2", 10.0 * 6.47e7 * 5.0 / 15.0),
        ];
        let ws: f64 = comps.iter().map(|c| c.2).sum();
        let em: f64 = comps.iter().map(|c| c.2 * e(c.0, c.1)).sum::<f64>() / ws;
        HC_EVNM / (em * Self::mu())
    }

    /// density-weighted |v|/(αc) of the 1s surfer at the equator (θ = 0.05 lattice).
    pub fn surfer_1s(&self) -> f64 {
        let s = &self.states[0];
        let (mut num, mut den) = (0.0, 0.0);
        for (i, v) in s.gf.iter().enumerate() {
            let r = (i as f64 + 0.5) * s.a;
            if r < 0.01 / ALPHA || r > 8.0 / ALPHA { continue; }
            let rho = v[0] * v[0] + v[1] * v[1];
            num += 2.0 * (v[0] * v[1]).abs();
            den += rho;
        }
        num / den / ALPHA
    }
}

// ─────────────── the 3-D state from (g, f): ψ = (g/r Ω_κm, i f/r Ω_−κm) (Dirac representation) ───────────────

/// Spherical harmonic Y_lm(θ, φ) for l ≤ 2.
pub fn ylm(l: i32, m: i32, th: f64, ph: f64) -> C {
    use std::f64::consts::PI;
    if m.abs() > l { return C::ZERO; }
    let (st, ct) = th.sin_cos();
    let e = C::expi(m as f64 * ph);
    let v = match (l, m.abs()) {
        (0, 0) => 0.5 * (1.0 / PI).sqrt(),
        (1, 0) => 0.5 * (3.0 / PI).sqrt() * ct,
        (1, 1) => -0.5 * (1.5 / PI).sqrt() * st,
        (2, 0) => 0.25 * (5.0 / PI).sqrt() * (3.0 * ct * ct - 1.0),
        (2, 1) => -0.5 * (7.5 / PI).sqrt() * st * ct,
        (2, 2) => 0.25 * (7.5 / PI).sqrt() * st * st,
        _ => 0.0,
    };
    // Condon–Shortley: Y_l,−m = (−1)^m conj(Y_lm)
    if m < 0 {
        let s = if m.abs() % 2 == 1 { -1.0 } else { 1.0 };
        (e * v).scale(s) * C::ONE
    } else {
        e * v
    }
}

/// Spinor spherical harmonic Ω_κ,m (m = m_j, as 2·m_j to stay integer).
pub fn omega(kappa: i32, two_m: i32, th: f64, ph: f64) -> [C; 2] {
    let l = l_of(kappa) as i32;
    let lf = l as f64;
    let m = two_m as f64 / 2.0;
    let (mlo, mhi) = ((two_m - 1) / 2, (two_m + 1) / 2);
    if kappa < 0 {
        // j = l + ½
        let a = ((lf + m + 0.5) / (2.0 * lf + 1.0)).max(0.0).sqrt();
        let b = ((lf - m + 0.5) / (2.0 * lf + 1.0)).max(0.0).sqrt();
        [ylm(l, mlo, th, ph) * a, ylm(l, mhi, th, ph) * b]
    } else {
        // j = l − ½
        let a = ((lf - m + 0.5) / (2.0 * lf + 1.0)).max(0.0).sqrt();
        let b = ((lf + m + 0.5) / (2.0 * lf + 1.0)).max(0.0).sqrt();
        [ylm(l, mlo, th, ph) * (-a), ylm(l, mhi, th, ph) * b]
    }
}

/// Dirac spinor (Dirac representation) at the point p (units λ_C) for the radial state `lv` of κ, m_j.
pub fn spinor(lv: &Level, kappa: i32, two_m: i32, p: [f64; 3]) -> [C; 4] {
    let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(1e-12);
    let th = (p[2] / r).clamp(-1.0, 1.0).acos();
    let ph = p[1].atan2(p[0]);
    let x = r / lv.a - 0.5;
    let i = (x.floor().max(0.0) as usize).min(lv.gf.len() - 2);
    let t = (x - i as f64).clamp(0.0, 1.0);
    let g = lv.gf[i][0] * (1.0 - t) + lv.gf[i + 1][0] * t;
    let f = lv.gf[i][1] * (1.0 - t) + lv.gf[i + 1][1] * t;
    let up = omega(kappa, two_m, th, ph);
    let lo = omega(-kappa, two_m, th, ph);
    [up[0] * (g / r), up[1] * (g / r), C::I * lo[0] * (f / r), C::I * lo[1] * (f / r)]
}

/// Surfer velocity in the Dirac representation: v = 2 Re(φ† σ χ)/(|φ|² + |χ|²).
pub fn velocity_dirac(s: &[C; 4]) -> [f64; 3] {
    let rho: f64 = s.iter().map(|c| c.norm2()).sum();
    if rho < 1e-300 { return [0.0; 3]; }
    let (p0, p1, c0, c1) = (s[0], s[1], s[2], s[3]);
    // φ†σχ components
    let sx = p0.conj() * c1 + p1.conj() * c0;
    let sy = (p0.conj() * c1).scale(1.0) * C::new(0.0, -1.0) + p1.conj() * c0 * C::I;
    let sz = p0.conj() * c0 - p1.conj() * c1;
    [2.0 * sx.re / rho, 2.0 * sy.re / rho, 2.0 * sz.re / rho]
}

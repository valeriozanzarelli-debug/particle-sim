//! The proton as a KNOT of the pion field: the B = 1 hedgehog Skyrmion (docs/NODO_3D.md), computed
//! here from its ordinary differential equation (no table, no file).
//!
//! Skyrme units (energy F_π/4e, length 2/(e F_π)); static energy of the hedgehog U = exp(i f(r) r̂·τ):
//!   E = 4π ∫ [ r² f′² + 2 sin²f (1 + f′²) + sin⁴f / r² ] dr ,
//! whose minimum obeys
//!   (r² + 2 sin²f) f″ + 2 r f′ + sin 2f (f′² − 1 − sin²f / r²) = 0 ,   f(0) = π, f(∞) = 0
//! (massless pions). We find the slope a = −f′(0) by shooting (bisection). Baryon density:
//!   b(r) = −sin²f f′ / (2π² r²)   (∫ b d³x = 1: one baryon, the winding number).
//! To turn the pure numbers into fm and MeV the model needs 2 measured inputs (Adkins–Nappi–Witten 1983:
//! F_π = 129 MeV and e = 5.45, fitted to the nucleon and Δ masses). That is a DEBT (see docs/ATOMO.md).

/// ANW 1983 fitted constants (MeV, dimensionless): the two inputs of the Skyrme model.
pub const F_PI_MEV: f64 = 129.0;
pub const E_SKYRME: f64 = 5.45;
/// ħc in MeV·fm (from ħ, c, e of costanti.rs: 197.3269804).
pub const HBARC_MEV_FM: f64 = 197.326_980_4;

/// Length unit 2/(e F_π) in fm.
pub fn unita_fm() -> f64 { 2.0 * HBARC_MEV_FM / (E_SKYRME * F_PI_MEV) }
/// Energy unit F_π/(4e) in MeV.
pub fn unita_mev() -> f64 { F_PI_MEV / (4.0 * E_SKYRME) }

pub struct Profilo {
    /// grid step (Skyrme units) and f on the grid r_i = i·h
    pub h: f64,
    pub f: Vec<f64>,
    pub df: Vec<f64>,
    /// −f′(0)
    pub a: f64,
}

fn rhs(r: f64, f: f64, g: f64) -> f64 {
    let s = f.sin();
    let s2f = (2.0 * f).sin();
    (s2f * (1.0 - g * g + s * s / (r * r)) - 2.0 * r * g) / (r * r + 2.0 * s * s)
}

/// Integrate from r ≈ 0 with slope −a. Returns (+1 overshoot: f < 0; −1 undershoot: f turns up; 0) and the
/// path up to the decision point.
fn spara(a: f64, h: f64, rmax: f64) -> (i32, Vec<f64>, Vec<f64>) {
    let r0 = 1e-6;
    let (mut f, mut g) = (std::f64::consts::PI - a * r0, -a);
    let mut r = r0;
    let mut fs = vec![std::f64::consts::PI];
    let mut gs = vec![-a];
    let n = (rmax / h) as usize;
    for i in 1..=n {
        let target = i as f64 * h;
        let dt = target - r;
        // RK4 on (f, g)
        let k1f = g;
        let k1g = rhs(r, f, g);
        let k2f = g + 0.5 * dt * k1g;
        let k2g = rhs(r + 0.5 * dt, f + 0.5 * dt * k1f, g + 0.5 * dt * k1g);
        let k3f = g + 0.5 * dt * k2g;
        let k3g = rhs(r + 0.5 * dt, f + 0.5 * dt * k2f, g + 0.5 * dt * k2g);
        let k4f = g + dt * k3g;
        let k4g = rhs(r + dt, f + dt * k3f, g + dt * k3g);
        f += dt / 6.0 * (k1f + 2.0 * k2f + 2.0 * k3f + k4f);
        g += dt / 6.0 * (k1g + 2.0 * k2g + 2.0 * k3g + k4g);
        r = target;
        fs.push(f);
        gs.push(g);
        if f < 0.0 { return (1, fs, gs); }
        if g > 0.0 { return (-1, fs, gs); }
    }
    (0, fs, gs)
}

impl Profilo {
    /// Solve the hedgehog profile on a grid of step h up to rmax (Skyrme units).
    pub fn calcola() -> Self {
        let h = 0.002;
        let rmax = 12.0;
        let (mut lo, mut hi) = (0.5, 5.0);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            let (s, _, _) = spara(mid, h, rmax);
            if s > 0 { hi = mid; } else { lo = mid; }
        }
        let a = 0.5 * (lo + hi);
        let (_, mut f, mut df) = spara(a, h, rmax);
        // the shot is trusted while f is clearly positive and falling; beyond, the exact massless tail
        // f = C / r² continued from the last trusted point
        let mut cut = f.len() - 1;
        for i in 1..f.len() {
            if f[i] < 0.02 || df[i] > -1e-6 { cut = i; break; }
        }
        f.truncate(cut + 1);
        df.truncate(cut + 1);
        let rc = cut as f64 * h;
        let c = f[cut] * rc * rc;
        let n = (rmax / h) as usize;
        for i in cut + 1..=n {
            let r = i as f64 * h;
            f.push(c / (r * r));
            df.push(-2.0 * c / (r * r * r));
        }
        Profilo { h, f, df, a }
    }

    pub fn r(&self, i: usize) -> f64 { i as f64 * self.h }

    /// Baryon density b(r) (per Skyrme-unit volume).
    pub fn barioni(&self, i: usize) -> f64 {
        if i == 0 {
            // limit r → 0: sin f ≈ a r, f′ = −a → b = a³ / 2π²
            return self.a.powi(3) / (2.0 * std::f64::consts::PI.powi(2));
        }
        let r = self.r(i);
        let s = self.f[i].sin();
        -s * s * self.df[i] / (2.0 * std::f64::consts::PI.powi(2) * r * r)
    }

    /// ∫ b d³x (should be 1).
    pub fn numero_barionico(&self) -> f64 {
        let mut s = 0.0;
        for i in 1..self.f.len() { s += 4.0 * std::f64::consts::PI * self.r(i).powi(2) * self.barioni(i) * self.h; }
        s
    }

    /// Static energy in units of 12π² (Faddeev–Bogomolny bound = 1; known value 1.2315).
    pub fn energia_12pi2(&self) -> f64 {
        let mut s = 0.0;
        for i in 1..self.f.len() {
            let r = self.r(i);
            let (sf, g) = (self.f[i].sin(), self.df[i]);
            s += (r * r * g * g + 2.0 * sf * sf * (1.0 + g * g) + sf.powi(4) / (r * r)) * self.h;
        }
        4.0 * std::f64::consts::PI * s / (12.0 * std::f64::consts::PI.powi(2))
    }

    /// Root-mean-square radius of the baryon density, fm.
    pub fn raggio_rms_fm(&self) -> f64 {
        let mut s = 0.0;
        for i in 1..self.f.len() { s += 4.0 * std::f64::consts::PI * self.r(i).powi(4) * self.barioni(i) * self.h; }
        s.sqrt() * unita_fm()
    }

    /// Mass of the classical knot in MeV (E × F_π/4e).
    pub fn massa_mev(&self) -> f64 { self.energia_12pi2() * 12.0 * std::f64::consts::PI.powi(2) * unita_mev() }


    /// Isospin moment of inertia Λ (Skyrme units, ħ_sk = 2e² in these units), T = ½ Λ ω²:
    /// Λ = (16π/3) ∫ r² sin²f (1 + f′² + sin²f / r²) dr, plus the massless tail f = C/r² beyond the grid.
    pub fn inerzia(&self) -> f64 {
        let mut s = 0.0;
        for i in 1..self.f.len() {
            let r = self.r(i);
            let (sf, g) = (self.f[i].sin(), self.df[i]);
            s += r * r * sf * sf * (1.0 + g * g + sf * sf / (r * r)) * self.h;
        }
        let (rl, c) = self.coda();
        s += c * c / rl;
        16.0 * std::f64::consts::PI / 3.0 * s
    }
    /// End of the grid and the tail constant C (f = C/r²).
    pub fn coda(&self) -> (f64, f64) {
        let i = self.f.len() - 1;
        let r = self.r(i);
        (r, self.f[i] * r * r)
    }
    /// ħ in Skyrme units: ħc / (energy unit × length unit) = 2 e².
    pub fn hbar_sk() -> f64 { 2.0 * E_SKYRME * E_SKYRME }
    /// M_Δ − M_N = 3ħ²/(2Λ) in MeV (ANW fixed e, F_π with it: a self-check of Λ, should be ≈ 293 MeV).
    pub fn delta_n_mev(&self) -> f64 { 3.0 * Self::hbar_sk().powi(2) / (2.0 * self.inerzia()) * unita_mev() }
    /// Angular velocity of the spinning nucleon (J = ½): ω = ħ √(J(J+1)) / Λ, in 1/s.
    pub fn omega_rot_s(&self) -> f64 {
        let w_sk = Self::hbar_sk() * (0.75f64).sqrt() / self.inerzia(); // per Skyrme time unit (length unit / c)
        w_sk / (unita_fm() * 1e-15 / 299_792_458.0)
    }
    /// Breathing (scaling) mode from Derrick's scaling: V(λ) = λE₂ + E₄/λ, kinetic ½ I λ̇² with
    /// I = 2·4π ∫ r⁴ f′² (1 + 2 sin²f / r²) dr. ω² = V″(1)/I = 2E₄/I. Returns (ω in 1/s, ħω in MeV).
    pub fn respiro(&self) -> (f64, f64) {
        let (mut e4, mut ii) = (0.0, 0.0);
        for i in 1..self.f.len() {
            let r = self.r(i);
            let (sf, g) = (self.f[i].sin(), self.df[i]);
            e4 += (2.0 * sf * sf * g * g + sf.powi(4) / (r * r)) * self.h;
            ii += r.powi(4) * g * g * (1.0 + 2.0 * sf * sf / (r * r)) * self.h;
        }
        let (rl, c) = self.coda();
        ii += 4.0 * c * c / rl;
        let e4 = 4.0 * std::f64::consts::PI * e4;
        let ii = 2.0 * 4.0 * std::f64::consts::PI * ii;
        let w_sk = (2.0 * e4 / ii).sqrt();
        let t_unit = unita_fm() * 1e-15 / 299_792_458.0;
        (w_sk / t_unit, w_sk * Self::hbar_sk() * unita_mev())
    }

    /// Balls drawn with the baryon density: (position in fm, relative density 0..1).
    pub fn palline(&self, n: usize, seed: u64) -> Vec<([f64; 3], f32)> {
        let mut rng = crate::carbonio::Rng::new(seed);
        let w: Vec<f64> = (0..self.f.len()).map(|i| self.r(i).powi(2) * self.barioni(i)).collect();
        let mut cdf = Vec::with_capacity(w.len());
        let mut s = 0.0;
        for x in &w { s += x; cdf.push(s); }
        let b0 = self.barioni(0);
        (0..n).map(|_| {
            let u = rng.uni() * s;
            let i = cdf.partition_point(|&c| c < u).min(w.len() - 1);
            let r = (i as f64 + rng.uni() - 0.5).max(0.0) * self.h;
            let ct = 2.0 * rng.uni() - 1.0;
            let ph = std::f64::consts::TAU * rng.uni();
            let st = (1.0 - ct * ct).sqrt();
            let rf = r * unita_fm();
            ([rf * st * ph.cos(), rf * st * ph.sin(), rf * ct], (self.barioni(i) / b0) as f32)
        }).collect()
    }
}

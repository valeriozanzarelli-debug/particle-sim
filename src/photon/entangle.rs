//! Polarization entanglement of two photons.
//!
//! Quantum (the real thing): the pair has ONE joint polarization state
//!   |Φ+⟩ = (|H⟩_A|H⟩_B + |V⟩_A|V⟩_B) / √2
//! and neither photon has a polarization of its own. Measuring A at angle a
//! gives + with the Born probability ‖(⟨a| ⊗ 1)|ψ⟩‖² = 1/2, and at that
//! instant B's state becomes the corresponding conditional state — wherever
//! B is. Result: E(a,b) = cos 2(a − b), and the CHSH value reaches 2√2,
//! while each side alone still sees 50/50 whatever the other does
//! (no signalling: the collapse cannot carry a message).
//!
//! Local hidden variables (for comparison): the pair leaves with a shared
//! polarization angle λ, and each side answers from λ and its own setting
//! only. No local model can exceed |S| = 2 (Bell's theorem).

use crate::complex::{linear, Jones, C};

const HALF_PI: f32 = std::f32::consts::FRAC_PI_2;

pub struct Pair {
    /// Joint state, index 2·iA + iB over the basis {H, V} of each beam.
    psi: [C; 4],
    /// Local hidden-variable model: the shared angle λ.
    hidden: Option<f32>,
    /// Per side: (CHSH setting index, angle, outcome +).
    results: [Option<(Option<usize>, f32, bool)>; 2],
}

impl Pair {
    pub fn phi_plus() -> Self {
        let s = std::f32::consts::FRAC_1_SQRT_2;
        Self { psi: [C::real(s), C::ZERO, C::ZERO, C::real(s)], hidden: None, results: [None; 2] }
    }

    pub fn hidden(lambda: f32) -> Self {
        Self { psi: [C::ZERO; 4], hidden: Some(lambda), results: [None; 2] }
    }

    /// True when each photon carries a definite polarization of its own.
    pub fn is_local(&self) -> bool { self.hidden.is_some() }

    /// What each photon carries before any measurement.
    pub fn local_jones(&self) -> Jones {
        match self.hidden {
            Some(l) => linear(l),
            None => [C::ZERO; 2],   // undefined: only the pair has a state
        }
    }

    pub fn record(&mut self, side: usize, idx: Option<usize>, angle: f32, plus: bool) {
        self.results[side] = Some((idx, angle, plus));
    }

    /// Both sides measured: (idx_a, a, +a, idx_b, b, +b).
    pub fn completed(&self) -> Option<(Option<usize>, f32, bool, Option<usize>, f32, bool)> {
        let (ia, a, pa) = self.results[0]?;
        let (ib, b, pb) = self.results[1]?;
        Some((ia, a, pa, ib, b, pb))
    }
}

/// Measure photon `side` (0 = A, 1 = B) with a polarizing beam splitter at
/// `angle`, using the uniform random number `r`. Returns the outcome and,
/// if the other photon is still unmeasured, its collapsed polarization.
pub fn measure(pair: &mut Pair, side: usize, angle: f32, r: f32) -> (bool, Option<Jones>) {
    let other_unmeasured = pair.results[1 - side].is_none();

    if let Some(lambda) = pair.hidden {
        // Deterministic local answer: + if λ is within 45° of the axis
        let plus = (angle - lambda).cos().powi(2) >= 0.5;
        return (plus, None);
    }

    if !other_unmeasured {
        // The pair already collapsed: this photon is in a definite state,
        // stored in `psi` as a product state by the first measurement.
        let j = factor(&pair.psi, side);
        let a = linear(angle);
        let p = (a[0].conj() * j[0] + a[1].conj() * j[1]).norm2();
        return (r < p, None);
    }

    // Conditional states of the other photon for the two outcomes
    let plus_v  = project(&pair.psi, side, angle);
    let minus_v = project(&pair.psi, side, angle + HALF_PI);
    let p_plus = plus_v[0].norm2() + plus_v[1].norm2();
    let plus = r < p_plus;
    let v = if plus { plus_v } else { minus_v };
    let n = (v[0].norm2() + v[1].norm2()).sqrt();
    let other = [v[0].scale(1.0 / n), v[1].scale(1.0 / n)];

    // Collapse to the product state |mine⟩ ⊗ |other⟩
    let mine = if plus { linear(angle) } else { linear(angle + HALF_PI) };
    let (ja, jb) = if side == 0 { (mine, other) } else { (other, mine) };
    for ia in 0..2 {
        for ib in 0..2 {
            pair.psi[2 * ia + ib] = ja[ia] * jb[ib];
        }
    }
    (plus, Some(other))
}

/// (⟨angle| on `side`) |ψ⟩ → unnormalized state of the other photon.
fn project(psi: &[C; 4], side: usize, angle: f32) -> Jones {
    let (c, s) = (angle.cos(), angle.sin());
    if side == 0 {
        [psi[0] * c + psi[2] * s, psi[1] * c + psi[3] * s]
    } else {
        [psi[0] * c + psi[1] * s, psi[2] * c + psi[3] * s]
    }
}

/// Factor of a product state belonging to `side`.
fn factor(psi: &[C; 4], side: usize) -> Jones {
    // For a product state ja ⊗ jb, pick the largest component to divide out.
    let (i, _) = psi.iter().enumerate()
        .max_by(|a, b| a.1.norm2().partial_cmp(&b.1.norm2()).unwrap()).unwrap();
    let (ia, ib) = (i / 2, i % 2);
    let v = if side == 0 { [psi[ib], psi[2 + ib]] } else { [psi[2 * ia], psi[2 * ia + 1]] };
    let n = (v[0].norm2() + v[1].norm2()).sqrt();
    [v[0].scale(1.0 / n), v[1].scale(1.0 / n)]
}

/// Quantum prediction for |Φ+⟩.
pub fn correlation_quantum(a: f32, b: f32) -> f32 { (2.0 * (a - b)).cos() }

/// Best local model (shared λ, sign rule): triangle wave in a − b.
pub fn correlation_local(a: f32, b: f32) -> f32 {
    let d = (a - b).rem_euclid(std::f32::consts::PI);
    let d = d.min(std::f32::consts::PI - d);
    1.0 - 4.0 * d / std::f32::consts::PI
}

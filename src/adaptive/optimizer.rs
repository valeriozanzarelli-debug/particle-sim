/// Self-improving adaptive optimizer.
///
/// Every OPTIMIZER_INTERVAL steps the system measures a "stability score"
/// and nudges the physical parameters toward configurations that produce
/// more stable (lower-energy, more coherent) structures.
///
/// This is analogous to gradient descent on the stability landscape:
/// parameters that improve stability are reinforced; parameters that reduce
/// it are corrected.

use bevy::prelude::*;
use std::collections::VecDeque;
use crate::config::{
    OPTIMIZER_INTERVAL, OPTIMIZER_LR, OPTIMIZER_MAX_DELTA, STABILITY_WINDOW,
    STRONG_A, STRONG_B, STRONG_R0, STRONG_R1, NUCLEAR_D, NUCLEAR_R_NUC,
    CHAIN_ATTENUATION, H_EFF,
};
use crate::physics::components::{Mass, ParticleKind, Velocity};
use crate::physics::integrator::SimClock;

// ─── Adaptive parameter resource ──────────────────────────────────────────────

#[derive(Resource, Clone)]
pub struct AdaptiveParams {
    // Strong force
    pub strong_a: f32,
    pub strong_b: f32,
    pub strong_r0: f32,
    pub strong_r1: f32,
    // Nuclear residual
    pub nuclear_d: f32,
    pub nuclear_r_nuc: f32,
    // EM chain
    pub chain_attenuation: f32,
    // Orbital barrier
    pub h_eff: f32,
    // Previous best score (for gradient direction).
    pub last_score: f32,
    // Which parameter to perturb next (round-robin).
    pub perturb_index: usize,
    // Last perturbation direction (+1 or -1) for each parameter.
    pub perturb_signs: [f32; 8],
}

impl Default for AdaptiveParams {
    fn default() -> Self {
        Self {
            strong_a: STRONG_A,
            strong_b: STRONG_B,
            strong_r0: STRONG_R0,
            strong_r1: STRONG_R1,
            nuclear_d: NUCLEAR_D,
            nuclear_r_nuc: NUCLEAR_R_NUC,
            chain_attenuation: CHAIN_ATTENUATION,
            h_eff: H_EFF,
            last_score: 0.0,
            perturb_index: 0,
            perturb_signs: [1.0; 8],
        }
    }
}

impl AdaptiveParams {
    pub fn as_slice_mut(&mut self) -> [&mut f32; 8] {
        [
            &mut self.strong_a,
            &mut self.strong_b,
            &mut self.strong_r0,
            &mut self.strong_r1,
            &mut self.nuclear_d,
            &mut self.nuclear_r_nuc,
            &mut self.chain_attenuation,
            &mut self.h_eff,
        ]
    }

    /// Names of each parameter for HUD display.
    pub fn param_names() -> [&'static str; 8] {
        [
            "strong_A", "strong_B", "r0_strong", "r1_strong",
            "nuclear_D", "r_nuc", "chain_att", "h_eff",
        ]
    }
}

// ─── Score history resource ───────────────────────────────────────────────────

#[derive(Resource, Default)]
pub struct StabilityHistory {
    pub scores: VecDeque<f32>,
    pub current_score: f32,
    /// Whether the system has reached a stable hydrogen atom.
    pub hydrogen_stable: bool,
    /// Count of consecutive stable steps.
    pub stable_streak: u32,
}

impl StabilityHistory {
    pub fn push(&mut self, score: f32) {
        if self.scores.len() >= STABILITY_WINDOW {
            self.scores.pop_front();
        }
        self.scores.push_back(score);
        self.current_score = score;
    }

    pub fn mean(&self) -> f32 {
        if self.scores.is_empty() {
            return 0.0;
        }
        self.scores.iter().sum::<f32>() / self.scores.len() as f32
    }

    pub fn trend(&self) -> f32 {
        if self.scores.len() < 2 {
            return 0.0;
        }
        let n = self.scores.len();
        let recent: f32 = self.scores.iter().rev().take(n / 2).sum::<f32>()
            / (n / 2) as f32;
        let older: f32 = self.scores.iter().take(n / 2).sum::<f32>() / (n / 2) as f32;
        recent - older
    }
}

// ─── Stability scoring ────────────────────────────────────────────────────────

/// Computes a scalar stability score from the current particle state.
/// Higher = more stable.
///
/// Components:
///   1. Low kinetic energy (particles not flying around) → +score
///   2. Low velocity variance (coherent motion) → +score
///   3. Proton formation (3 quarks close together) → +bonus
///   4. Electron in Bohr-radius range → +bonus
pub fn compute_stability_score(
    particles: &Query<(&Transform, &Velocity, &Mass, &ParticleKind)>,
) -> f32 {
    use crate::config::{BOHR_RADIUS_FM, FM_TO_WORLD};

    let mut total_ke: f32 = 0.0;
    let mut vel_sum = Vec3::ZERO;
    let mut vel_sq_sum: f32 = 0.0;
    let mut count = 0u32;

    let mut quark_positions: Vec<Vec3> = Vec::new();
    let mut electron_pos: Option<Vec3> = None;
    let mut proton_pos: Option<Vec3> = None;

    for (t, vel, mass, kind) in particles.iter() {
        if *kind == ParticleKind::Mediator {
            continue;
        }
        let ke = 0.5 * mass.0 * vel.0.length_squared();
        total_ke += ke;
        vel_sum += vel.0;
        vel_sq_sum += vel.0.length_squared();
        count += 1;

        match kind {
            ParticleKind::Quark => quark_positions.push(t.translation),
            ParticleKind::Electron => electron_pos = Some(t.translation),
            ParticleKind::Proton => proton_pos = Some(t.translation),
            _ => {}
        }
    }

    if count == 0 {
        return 0.0;
    }

    // Low KE → high score component.
    let ke_score = 1.0 / (1.0 + total_ke * 0.001);

    // Low velocity variance → coherence score.
    let mean_v_sq = (vel_sum / count as f32).length_squared();
    let variance = (vel_sq_sum / count as f32 - mean_v_sq).max(0.0);
    let coherence_score = 1.0 / (1.0 + variance * 0.01);

    // Quark confinement bonus: are the quarks close together?
    let quark_bonus = if quark_positions.len() >= 3 {
        let centroid = quark_positions.iter().copied().fold(Vec3::ZERO, |a, b| a + b)
            / quark_positions.len() as f32;
        let spread: f32 = quark_positions
            .iter()
            .map(|p| (*p - centroid).length())
            .sum::<f32>()
            / quark_positions.len() as f32;
        1.0 / (1.0 + spread * 10.0)
    } else {
        0.0
    };

    // Electron orbital bonus: is the electron near the Bohr radius?
    let nucleus = proton_pos
        .or_else(|| {
            if quark_positions.len() >= 3 {
                Some(
                    quark_positions.iter().copied().fold(Vec3::ZERO, |a, b| a + b)
                        / quark_positions.len() as f32,
                )
            } else {
                None
            }
        });

    let electron_bonus = if let (Some(e_pos), Some(n_pos)) = (electron_pos, nucleus) {
        let r_world = (e_pos - n_pos).length();
        let r_fm = r_world / FM_TO_WORLD;
        let bohr = BOHR_RADIUS_FM;
        // Gaussian centred on Bohr radius.
        let diff = (r_fm - bohr) / (bohr * 0.2);
        (-diff * diff * 0.5).exp()
    } else {
        0.0
    };

    ke_score * 0.3 + coherence_score * 0.3 + quark_bonus * 0.2 + electron_bonus * 0.2
}

// ─── Optimiser system ─────────────────────────────────────────────────────────

pub fn run_optimizer(
    clock: Res<SimClock>,
    particles: Query<(&Transform, &Velocity, &Mass, &ParticleKind)>,
    mut params: ResMut<AdaptiveParams>,
    mut history: ResMut<StabilityHistory>,
) {
    if clock.step_count % OPTIMIZER_INTERVAL != 0 {
        return;
    }

    let score = compute_stability_score(&particles);
    history.push(score);

    // Check hydrogen stability: sustained high score.
    if score > 0.75 {
        history.stable_streak += 1;
        if history.stable_streak > 50 {
            history.hydrogen_stable = true;
        }
    } else {
        history.stable_streak = 0;
    }

    // Gradient direction: did the score improve since last optimisation?
    let improved = score > params.last_score;
    params.last_score = score;

    // Round-robin: adjust one parameter at a time.
    let idx = params.perturb_index % 8;
    let sign = if improved {
        params.perturb_signs[idx]
    } else {
        -params.perturb_signs[idx]
    };

    let lr = OPTIMIZER_LR;
    let max_delta = OPTIMIZER_MAX_DELTA;

    {
        let slices = params.as_slice_mut();
        let val = &mut *slices[idx];
        let delta = (*val * lr * sign).clamp(-*val * max_delta, *val * max_delta);
        *val += delta;
        // Keep parameters in physically sensible ranges.
        if *val < 0.001 {
            *val = 0.001;
        }
    }

    params.perturb_signs[idx] = sign;
    params.perturb_index += 1;
}

//! The only numbers that enter: physical constants (CODATA 2018), nuclear charges and nuclear masses.
//! Atomic units (bohr, hartree) are DERIVED here from ħ, e, m_e, ε0, c — nothing else is tabulated.

pub const HBAR: f64 = 1.054_571_817e-34; // J s
pub const E_CH: f64 = 1.602_176_634e-19; // C
pub const M_E: f64 = 9.109_383_7015e-31; // kg
pub const EPS0: f64 = 8.854_187_8128e-12; // F/m
pub const C_L: f64 = 299_792_458.0; // m/s
/// Nuclear masses (the nuclei only), kg: proton and alpha particle.
pub const M_PROTON: f64 = 1.672_621_923_69e-27;
pub const M_ALPHA: f64 = 6.644_657_3357e-27;

use std::f64::consts::PI;

/// Bohr radius a0 = 4π ε0 ħ² / (m_e e²), metres.
pub fn bohr_m() -> f64 { 4.0 * PI * EPS0 * HBAR * HBAR / (M_E * E_CH * E_CH) }
/// Hartree = ħ² / (m_e a0²), joules.
pub fn hartree_j() -> f64 { HBAR * HBAR / (M_E * bohr_m().powi(2)) }
pub fn hartree_ev() -> f64 { hartree_j() / E_CH }
pub fn bohr_ang() -> f64 { bohr_m() * 1e10 }
/// hartree -> cm⁻¹ (E / (h c))
pub fn hartree_cm1() -> f64 { hartree_j() / (2.0 * PI * HBAR * C_L) / 100.0 }
/// Atomic unit of time ħ / E_h, seconds.
pub fn au_time_s() -> f64 { HBAR / hartree_j() }
/// Wavelength (nm, vacuum) of a photon carrying ΔE hartree: λ = h c / ΔE.
pub fn nm_from_hartree(de: f64) -> f64 { 2.0 * PI * HBAR * C_L / (de * hartree_j()) * 1e9 }
/// Masses in units of m_e (atomic units of mass).
pub fn proton_me() -> f64 { M_PROTON / M_E }
pub fn alpha_me() -> f64 { M_ALPHA / M_E }
/// Fine-structure constant α = e² / (4π ε0 ħ c) (derived, used for photon emission rates).
pub fn alpha_fs() -> f64 { E_CH * E_CH / (4.0 * PI * EPS0 * HBAR * C_L) }

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn derived_units_match_codata() {
        assert!((bohr_ang() - 0.529_177_210_9).abs() < 1e-9);
        assert!((hartree_ev() - 27.211_386_2).abs() < 1e-6);
        assert!((1.0 / alpha_fs() - 137.035_999).abs() < 1e-5);
    }
}

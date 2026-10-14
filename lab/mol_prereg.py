"""mol_prereg.py - pre-registration of study mol_ (docs/MOLECOLE.md), phase 1: LiH with four electrons, tier "esatto".
Writes lab/out/mol_prereg.json. Committed BEFORE the production runs. Later changes go in separate amendment files."""
import json, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
pre = {
    "study": "mol_ (Molecole piccole)",
    "date": "2026-10-13",
    "tiers": {"esatto": "few-electron exact: ECG with the rule's learned mass (mol_ecg.py), all electrons correlated",
              "approssimato": "grid Hartree-Fock on the agglomerated level + learned correlation correction (tavg_ code, "
                              "imported read-only); to be pre-registered in an amendment before its runs"},
    "LiH": {
        "tier": "esatto",
        "operator": "sum p^2/2m_kin - sum Z/r + sum 1/r_ij + nuclear repulsion; m_kin = 1.0000000045603115 (rule m0, "
                    "multi_e v5); nuclei clamped (Born-Oppenheimer, debt); singlet spin function "
                    "(12 singlet)x(34 singlet), all 4! permutations",
        "inputs": "m_e, alpha (in the rule), Z = 3 and 1, masses of the 7Li nucleus and of the proton (vibration only)",
        "config": {
            "R0": 3.015,
            "K_RE": 600, "nsweep_re": 3, "seed_re": 11,
            "K_LI": 300, "nsweep_li": 3, "seed_li": 7,
            "ntrial_grow": 60, "ntrial_ref": 30,
            "RGRID": [2.5, 2.7, 2.85, 3.015, 3.15, 3.3, 3.5, 3.8, 4.2, 5.0],
            "nsweep_curve": 2, "seed_curve": 100,
            "curve_start": "R0 basis with centres stretched z -> z R/R0, then nsweep_curve refinement sweeps",
            "fit_window": [2.7, 3.5], "fit": "polynomial degree 4 in R - R0 on the grid points in the window; R_e and "
                                            "omega_e from its minimum and curvature",
            "vibration": "sinc-DVR on a cubic spline of the curve, R in [2.5, 5.0], reduced mass of the nuclei",
            "dipole": "mu(R) = R (charge of H at R, Li at 0) - <sum z_i> from the ECG eigenvector; averaged over v = 0"
        },
        "development_observations_declared": [
            "Li atom 3e, K = 160, 2 sweeps: -7.477217 (0.84 mEh above -7.478060)",
            "LiH R = 3.015, 4e, K = 400, 2 sweeps at every 40/100: -8.068244 (2.31 mEh above -8.070553); 6.5 min",
            "K = 200: 5.95 mEh; K = 300: 3.59 mEh -> slow algebraic convergence, K = 600 expected ~1.5 mEh",
            "a sign error in the secular equation of the fast candidate test was found and fixed in development"
        ],
        "references_exam_only": {
            "E_BO_LiH_3.015": -8.070553, "E_Li": -7.478060,
            "R_e_a0": 3.015, "D_e_eV": 2.515, "omega_e_cm": 1405.6, "dipole_v0_D": 5.882,
            "DeltaG_half_cm_report_only": 1359.7
        },
        "exams": {
            "L0": {"what": "E(LiH, R0) ECG vs nonrelativistic BO -8.070553 Eh (machine gate)", "ref": -8.070553,
                   "tol_Eh": 0.0020},
            "L1": {"what": "bond length R_e from the curve minimum", "ref": 3.015, "tol_rel": 0.005},
            "L2": {"what": "D_e = E(Li) + E(H) - E(R_e) vs 2.515 eV", "ref": 2.515, "tol_eV": 0.0434},
            "L3": {"what": "harmonic omega_e vs 1405.6 cm-1", "ref": 1405.6, "tol_rel": 0.02},
            "L4": {"what": "dipole moment averaged over v = 0 vs 5.882 D", "ref": 5.882, "tol_rel": 0.02},
            "L5": {"what": "Li atom with the same machine vs -7.478060 Eh", "ref": -7.478060, "tol_Eh": 0.0005}
        },
        "honesty": [
            "This is VALIDATION (LiH is known), not discovery.",
            "Reference energies are nonrelativistic BO; experimental D_e, omega_e, dipole contain small effects not "
            "modelled (relativity, adiabatic corrections): below the thresholds.",
            "The ECG basis is chosen at random and refined: a solver, not the rule walking (debt, as aggmol phase 2).",
            "Scale-down if needed: if a run cannot finish in chunks, K is reduced and declared in an amendment."
        ]
    }
}
os.makedirs(OUT, exist_ok=True)
json.dump(pre, open(os.path.join(OUT, "mol_prereg.json"), "w"), indent=1, ensure_ascii=False)
print("written", os.path.join(OUT, "mol_prereg.json"))

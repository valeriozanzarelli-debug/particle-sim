"""big_prereg_b.py - content of the phase-B pre-registration (written AFTER R1/R2 results, BEFORE big_b.py runs).
Written to lab/out/big_prereg_B.json by `python big_prereg.py B`.

Status from phase R: R1 PASS (fp 0/25 in every null class, 3/3 true positives) -> the big step is allowed.
R2 FAIL at theta = 0.125 and 0.0625 (M = 1 relaxes as fast as M = 16) -> per big_prereg_R.json: theta = 0.125,
the Born-PINNING exam is NOT meaningful; B is run for SHAPE only. (Equivariance from |psi|^2 is still checked:
it does not depend on fake mixing.)
"""
B = {
    "status_from_R": "R1 pass; R2 fail at both thetas -> theta = 0.125, Born pinning not graded (reported only)",
    "rule": "pre_core.LatticeWalk, 24x24 BOX with mirror edges (chosen in advance; k-space Wilson needs translation "
            "invariance, which varying link weights break), theta = 0.125, surfers lam = 0, N = 20000",
    "reinforcement": {"eta": 0.3, "mu": 2, "wbar": 0.3, "w_min": 0.02, "K": 10,
                      "rule": "every K ticks: J = surfer crossings per interior link; w <- (1-eta) w + eta wbar J^mu/"
                              "mean(J^mu); clip [w_min, 1]; rescale mean to wbar (3 passes); final clip. Start w = wbar."},
    "runs_per_seed": {"reinf_uniform": "M = 16, uniform start, reinforcement on",
                      "plain_uniform": "M = 16, uniform start, eta = 0 (no reinforcement)",
                      "reinf_psi2": "M = 16, |psi|^2 start, reinforcement on"},
    "T": 5000, "seeds": {"1": "psi_seed 301, surf_seed 401", "2": "302, 402", "3": "303, 403"},
    "omega0": "1.8 theta = 0.225", "cells": "4x4 (F = 0.000875)",
    "exams": {
        "born_equivariance": "reinf_psi2: max_t H <= 3F in 3/3 seeds (graded)",
        "born_pinning": "NOT graded (R2 failed). Reported: mean H over the last 10% of ticks, reinf vs plain.",
        "highways": "reinf_uniform final w: Gini > 0.2 AND top-10% links carry > 30% of the weight, in 3/3 seeds (graded)",
        "shape": "repaired reader (big_core.read_shape) on the vacuum MI of the 24x24 box with the final reinforced "
                 "weights. 'Structure emerged' only if holes >= 1, OR pieces > 1, in 3/3 seeds with the SAME (pieces, "
                 "holes). Reference: uniform w = 0.3 box (already read in R1: 1 piece, 0 holes). Otherwise: 'no shape "
                 "emerged'.",
        "compact_hidden_direction": "claimed only if (D_small - D_large) of the reinforced geometry exceeds the "
                                    "reference by > 0.5 in 3/3 seeds (reference value seen in R1: 2.07 - 1.56 = 0.51).",
        "incastro": "localized region: the lowest-|omega| eigenmode of the final walk occupies < 10% of the sites "
                    "(1/(Ns sum rho^2) < 0.1) in 3/3 seeds while the reference does not.",
        "vortices": "phase-winding count (pre_core.windings) at start and end, reported",
        "mass_spectrum": "lowest 12 |omega| near omega = 0 (and near pi) of the final reinforced walk vs reference; "
                         "ratios omega_k/omega_1 reported; 'large ratio' = > 10 noted. No fit to any real mass.",
    },
}

"""nucleo1s_prereg.py - pre-registration of the study "1s pulito vicino al nucleo" (docs/NUCLEO_1S.md).
Written and committed BEFORE the production runs (lab/nucleo1s_run.py). Writes lab/out/nucleo1s_prereg.json.

Development already seen (declared, scratch scripts, not in the repo):
  * the graded-mesh walk reproduces RadialWalk and the Laboratorio Atomo numbers (1s +39.914 %, 2p3/2 +3.456 %,
    3d5/2 +1.828 %, 4f7/2 +1.270 % at Z = 60, a = 0.5);
  * uniform convergence at Z = 60 for 1s, 2p3/2, 3d5/2, 4f7/2, 2p1/2 down to a = 0.5/2^9; the analytic K and B for
    Z = 60 (computed after seeing the uniform data, but with no free number) and the check (e - K a^2g)/a^2 -> B;
  * nested zoom meshes at Z = 60 (a_out 0.5 ... 0.0625); a prototype of the automatic zoom (budgets 1 % and 0.3 %:
    true errors 0.67 % / 0.15 %, 139 / 213 cells); the agglomerated nucleus (exact reflection and the 2-number rule
    at R = 0.5, 1, 2);
  * Z = 1: uniform vs zoom at a = 0.2 ... 0.0125 (zoomed error / a^2 -> 0.332), crossover check at Z = 1, 2, 4, 8;
  * Z = 1 hybrid zig-zag + agglomerated rule, prototype runs (numbers in the JSON field dev_hybrid): coarse
    potential shell-averaged and point-sampled, J0 = 1, 2, 3, levels B = 256 and 512. HY1/HY2 below are therefore
    CONFIRMATION runs of a configuration chosen after seeing them (weak exams, declared); HY3 (Z = 2) is not seen.
NOT seen before this file: the uniform convergence at Z = 20, 40, 80, 92; any production zoom run with the frozen
algorithm below; the zoomed Z = 1 s-level extrapolation (H1z) and 1S-2S; the hybrid exams with the frozen settings.
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_core import ALPHA, gamma_of, dirac_binding
from nucleo1s_theory import inner_beta, contact_K, bulk_B

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
pred = {}
for Z in (20, 40, 60, 80, 92):
    Za = Z * ALPHA
    b = inner_beta(Za)
    pred[str(Z)] = dict(gamma=gamma_of(-1, Za), two_gamma=2 * gamma_of(-1, Za), beta=b, K=contact_K(Za, beta=b),
                        B=bulk_B(Za))

P = {
    "study": "nucleo1s", "doc": "docs/NUCLEO_1S.md",
    "inputs": "m_e (theta = a), alpha, Z (point nucleus); proton mass only in the reduced-mass factor (Z = 1)",
    "rule": "zig-zag radial rule (lab/zz_core.py RadialWalk) on graded meshes (lab/nucleo1s_core.py GradedWalk): "
            "rotation of each cell with dt = its width, links at c; equal cells = RadialWalk exactly",
    "predictions_law": {
        "statement": "relative 1s binding error e(a) = K a^(2 gamma) + B a^2 + higher; K from the scale-free inner "
                     "problem (beta), B from the symmetric BCH second-order Hamiltonian (analytic continuation)",
        "values": pred},
    "exams": {
        "L1": {"what": "uniform meshes a = 0.5/2^k, k = 0..11, Z in {20, 40, 80, 92} (not seen) and 60 (seen): "
                       "R = (e(a) - K a^2g)/(B a^2) at the finest a (0.5/2^11)",
               "threshold": "|R - 1| <= 0.02 for every Z"},
        "L2": {"what": "measured order of the contact term: p from (e - B a^2) between a = 0.5/2^10 and 0.5/2^11",
               "threshold": "|p - 2 gamma| <= 0.02 for every Z"},
        "L3": {"what": "p and d states (kappa = -2, -3) at the same Z: local order between the two finest a",
               "threshold": "|p - 2| <= 0.01"},
        "D1": {"what": "doubling: the one-tick walk at Z = 60, a = 0.5 has an eigenvalue at w_1s + pi/a whose "
                       "state is the staggered copy (-1)^i of the 1s with the OTHER boundary at r = 0 (refl +1); "
                       "check: w(refl -1 alias) - pi/a == w_1s(refl +1) and the 1s itself is not moved",
               "threshold": "agreement <= 1e-10 (relative)"},
        "Z1": {"what": "automatic zoom from the Laboratorio Atomo lattice (equal cells a = 0.5, Z = 60, rmax as "
                       "orbitali.rs), observable = 1s, budget (estimated error) 0.5 %",
               "threshold": "1s within 1 % of Dirac exact"},
        "Z2": {"what": "cost: zoom cells vs the equal-cell lattice that reaches the same 1s error "
                       "(found by bisection on a), for budgets 0.5 %, 0.1 %, 0.02 %",
               "threshold": "zoom cells <= 1/3 of the uniform cells, for each budget"},
        "Z3": {"what": "2p3/2, 3d5/2, 4f7/2 on the zoomed mesh (budget 0.5 %), extended to their own rmax with "
                       "cells 0.5: they must not get worse than on the equal-cell a = 0.5 lattice",
               "threshold": "err_zoom <= err_uniform + 0.05 percentage points"},
        "Z4": {"what": "honesty of the controller: estimated / true error of the 1s at the stop",
               "threshold": "0.5 <= ratio <= 2 for each budget"},
        "N1": {"what": "agglomerated nucleus: the zoom cells inside R = 0.5 (budget 0.1 % mesh) replaced by a "
                       "2-number reflection exp(i(d0 + tau w)) learned from the fine rule on w in [0.85, 1.0], "
                       "separately for kappa = -1 and kappa = -2, -3, -4; levels 1s..4s, 2p3/2, 3d5/2, 4f7/2",
               "threshold": "each level within 0.05 % (relative) of the same mesh with the cells"},
        "H1z": {"what": "Z = 1, s levels 1s, 2s, 3s with the zoom near the nucleus (outer cells theta = 0.05, "
                        "0.025, 0.0125; zoom 14 levels x 16 cells), extrapolated with a polynomial in theta^2 "
                        "through the three; vs Dirac (infinite mass). Same threshold as the failed H1 of ZIGZAG",
                "threshold": "|relative error| <= 1e-5 for each (report also vs 1e-7)"},
        "H2z": {"what": "Z = 1, 1S-2S frequency (reduced mass) from H1z vs the Dirac formula with reduced mass",
                "threshold": "|lattice error| <= 1 GHz (the QED Lamb shift ~8 GHz must then show up vs experiment)"},
        "HY1": {"what": "Z = 1 hybrid: zig-zag zoom for r < (J0 + 1/2) h, learned agglomerated rule (B = 256, "
                        "cell a0/16) outside, J0 = 3, point-sampled coarse Coulomb; 1s x reduced mass vs Dirac x "
                        "reduced mass (13.598468 eV). SEEN in development (-0.024 % without reduced mass): weak exam",
                "threshold": "within 0.1 %"},
        "HY2": {"what": "same with the a0/8 level (B = 512): 2s vs Dirac (re-judges A2b, +1.003 % with the 3-D "
                        "coarse rule). SEEN in development (-0.096 %): weak exam",
                "threshold": "within 1 %"},
        "HY3": {"what": "NOT seen: He+ (Z = 2) with the same rule and settings (B = 256, J0 = 3): 1s and 2s vs Dirac "
                        "at Z = 2 (the cell is now a_Z/8 for 1s)",
                "threshold": "within 0.5 % each"},
        "HYs": {"what": "systematics reported without threshold: J0 = 2, 3, 4; shell vs point potential; coarse "
                        "rule alone (no zig-zag zone)", "threshold": None},
    },
    "frozen_settings": {
        "uniform_law_runs": "a = 0.5/2^k, k = 0..11, rmax = 30 + 9/(Z alpha); kappa = -1, -2, -3",
        "autozoom": "start equal cells 0.5, split every cell -> eta (relative level change), eps = eta/(1-2^-p), "
                    "p = 2 (p = 2 gamma for the cell touching r = 0), estimate = |sum eps|; split cells with "
                    "eta >= 0.3 max eta; 2:1 balance; max 40 iterations",
        "rmax": "orbitali.rs: 30 + 9 n^2 / (Z alpha) lambda_C",
        "listening": "Newton on arg mu(w) with shift-invert at mu = 1, start 1.3 x Dirac (s) or 1.05 x Dirac",
        "hybrid": "fine zig-zag outer cell theta = a0/4096, 12 zoom levels x 16 cells; coarse s-wave stencil "
                  "|n| <= 48 from the direction average (400 Fibonacci directions) of the learned generator; "
                  "point-sampled Coulomb -alpha/r_j (consistent with point values u_j); J0 = 3; full 1-D kernel "
                  "(no truncation); energy from the generator (G/tau + V), not from the step",
    },
    "dev_hybrid": {"B256_shell": {"coarse_1s": -0.5984, "J0_2_1s": -0.2909, "J0_3_1s": -0.2835},
                   "B256_point": {"coarse_1s": 0.1619, "J0_2_1s": 0.1037, "J0_3_1s": -0.0243,
                                  "coarse_2s": 0.0994, "J0_3_2s": 0.0065},
                   "B512_point": {"coarse_1s": 0.5579, "J0_3_1s": -0.2448, "coarse_2s": 0.3044, "J0_3_2s": -0.0961},
                   "units": "percent vs Dirac, infinite nuclear mass"},
    "free_parameters": "0 physical; algorithm settings above (budget, 0.3, 2:1, J0, zoom levels) chosen by hand",
}
os.makedirs(OUT, exist_ok=True)
json.dump(P, open(os.path.join(OUT, "nucleo1s_prereg.json"), "w"), indent=1, default=float)
print(json.dumps(pred, indent=1, default=float))

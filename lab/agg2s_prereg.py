"""agg2s_prereg.py - PRE-REGISTRATION of the study "2s all'agglomeratore" (prefix agg2s_), written BEFORE the production
runs (lab/agg2s_run.py, lab/agg2s_contact.py). Writes lab/out/agg2s_prereg.json.

Development explorations already done (declared, not production): lab/agg2s_dev_box.py (box error of the 2s with the
EXACT Schroedinger kinetic, cells a0/2 and a0/4, boxes 16..48 a0, 4 twists), lab/agg2s_dev_box2p.py (same for the 2p),
lab/agg2s_dev_contact.py (contact density and <p^4> lost by the cell, 1s/2s/3s at a0/2, a0/4, 1s at a0/8; Born tail).
Their numbers fix the PREDICTIONS below; the production runs use the LEARNED rule (other kinetic, other states) and new
geometries (3s at L = 48 a0 with twists), so the exams are not the development numbers."""
import json, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
box = json.load(open(os.path.join(OUT, "agg2s_dev_box.json")))
b2p = json.load(open(os.path.join(OUT, "agg2s_dev_box2p.json")))

import agg2s_core as c
rel = lambda k: box[k]["2s"] / (-c.RY / 4) - 1
pred_box2s = rel("d4_N128_tw000") - 0.5 * (rel("d4_N128_tw000") + rel("d4_N128_tw111"))
pred_box2p = b2p["d4_N128_tw000_2p"]["E_rel"] - 0.5 * (b2p["d4_N128_tw000_2p"]["E_rel"] + b2p["d4_N128_tw111_2p"]["E_rel"])

P = dict(
    study="agg2s_ (2s all'agglomeratore)",
    date="2026-10-13",
    inputs="only -Z alpha/r (Z = 1) and the learned rules of AGGLOMERATORE COMPLETO (read-only): chain A, aggfull_relearn.pkl (K), "
           "aggfull_relearnT.pkl (T), aggfull_learn.pkl (coupling: gamma0 = 1, Darwin coefficient d = -gk). Dirac / Schroedinger "
           "/ hydrogen analytic numbers ONLY as exams.",
    new_pieces={
        "twisted boundary conditions": "psi(x + L e_a) = (+-1) psi(x): the box = a crystal of the atom's images; the single atom "
                                       "= average over the twists of the crystal band. Production average: (tw000 + tw111)/2 "
                                       "(cancels faces and corners of the image crystal exactly; edges survive). No parameter.",
        "Born tail of the nucleus cell": "beyond the coarse Brillouin cube the coarse rule's own equation (learned mass m, Coulomb "
                                         "input) gives psi(k) -> 8 pi m Z alpha psi(0) / k^4. Geometric constants of the cube: "
                                         "I = coef1 h (coef1 = <max|n_i|>_sphere 4pi/(pi (2pi)^3)), J = coef3 h^3 "
                                         "(coef3 = <max|n_i|^3> 4pi/(3 pi^3 (2pi)^3)). rho_B = rho_c/(1 - 8 pi m Z alpha I)^2; "
                                         "<k^4>_B = <k^4>_c + (8 pi m Z alpha)^2 rho_B I; cell energy dE_B = -2 m (4 pi Z alpha)^2 rho_c J. "
                                         "No parameter: everything from the cube geometry, the learned mass and -Z alpha/r.",
    },
    configurations=dict(
        main="cells a0/8 (chain A B = 512), N = 256, L = 32 a0 (+-16 a0, the aggfull box), alias-free Galerkin field (aggfull_core), "
             "gamma0 = 1; kinetic K (aggfull_relearn, B = 512) = MAIN, T (aggfull_relearnT) = SECOND, S = exact k^2/2 (m = 1) = cell reference; "
             "twists 000 and 111 (K also 100, 110: band-shape diagnostic); states 1s, 2s (deflated against the 1s of the same twist), 2p_z.",
        first_order="dV (learned Darwin + spin-orbit, first order) taken from the aggfull periodic runs qK8_1s / qK8_2s / qK8_2p "
                    "(1s: 6.3947e-4 eV; 2s: 7.9592e-5 eV; 2p: lower pair -3.01236e-5, upper quartet +1.51156e-5 eV); the box changes "
                    "them by < 1e-7 relative (contact quantities). Used also for T (difference of the states: < 1e-6 eV).",
        solver="lean single-vector LOBPCG (agg2s_core.lobpcg1), residual criterion 2e-8 mc^2 (as aggfull_run), warm start from the "
               "aggfull vectors (projected on the twisted sector)",
        contact_runs="lab/agg2s_contact.py: exact Schroedinger, cells a0/4, L = 48 a0 (N = 192), twists 000/111, 1s, 2s, 3s (3s is NEW: blind)",
        a16="1s at a0/16 (B = 256, N = 256, +-8 a0), kinetic T and K (vector warm start from aggfull eK16_1s), for the contact rule",
    ),
    predictions_from_development=dict(
        box_shift_2s_L32_rel=pred_box2s,
        box_shift_2p_L32_rel=pred_box2p,
        cell_error_S_2s_a0_8_rel="about -6e-5 (h^3 law from a0/2: -3.46e-3 and a0/4: -4.62e-4)",
        expected_2s_twisted_K_vs_Dirac="about +0.033 % (= mass of the fine rule theta^2/3 = +0.0373 % plus cell -0.006 % plus alpha^4 error +0.001 %)",
    ),
    exams={
        "BOX": "2s box shift with the learned rule K: (E(000) - E_avg)/|E_avg| within +-15 %% of the predicted %.5e" % pred_box2s,
        "BOX2p": "2p_z box shift with K: (E(000) - E_avg)/|E_avg| within +-20 %% of the predicted %.5e (sign: the periodic box "
                 "UNDER-binds the 2p, odd across the z faces)" % pred_box2p,
        "E2X": "four levels (1s1/2, 2s1/2, 2p1/2, 2p3/2) a0/8, kinetic K, twist average + dV, each within 0.1 % of Dirac (the "
               "failed E2 of aggfull, same threshold)",
        "E2Xs": "2s1/2 (K, twist average + dV) within 0.05 % of Dirac (stronger)",
        "E2XT": "same as E2X with the kinetic T",
        "BUD": "budget closure for the 2s (K): (E_K(000) + dV)/E_Dirac - 1 = mass (m_kin - 1) + cell (S runs, twist avg, vs -Ry/4) + "
               "alpha^4 ((kin.remainder_K + dV - (E_Dirac - E_Schr))/|E_Dirac|) + box (E_K(000) - E_K,avg)/|E| within 3e-5 (absolute, "
               "relative units); the same without the box term for the twisted average",
        "E2B": "with the Born cell correction dE_B (rho_c of the coarse state): 1s and 2s (K and T, twist avg + dV + dE_B) within "
               "0.005 % of m_kin x Dirac (the coarse level against the physics of its own mass)",
        "BORN_S": "exact Schroedinger at a0/8 (S runs, twist avg): Born-corrected rho_B and <k^4>_B within 0.5 % of the exact values "
                  "(1/(pi n^3), 8/n^3 - 3/n^4), and E_S + dE_B within 2e-5 relative of -Ry/n^2, for 1s and 2s",
        "BORN_3s": "BLIND 3s (contact runs, a0/4, L = 48, twist avg): rho_B and <k^4>_B within 2 % of the exact values; E + dE_B "
                   "within 2e-4 relative of -Ry/9; universality: the lost contact per unit coarse contact (rho - rho_c)/rho_c of 1s, 2s, 3s "
                   "equal within 5 %",
        "F3S": "alpha^4 shifts of the s states with the contact rule of the nucleus cell: A = c4 alpha^4 <k^4>_B + 4 pi d alpha^4 rho_B "
               "(c4 = isotropic k^4 coefficient of the learned band, d = -gk of the learned coupling, B = 512 or 256): 1s (a0/16), "
               "1s (a0/8), 2s (a0/8, twist avg) within 5 % of Dirac - Schroedinger (the F3 threshold), kinetic T (MAIN) and K",
        "F3Sstrong": "same, within 1.5 %",
    },
    declared_limits=[
        "the twist average leaves the edge terms of the image crystal (dev: ~2e-5 relative at L = 32 for the 2s)",
        "dV from the periodic aggfull runs",
        "the fine rule's own deviations from Dirac (mass theta^2/3 = 3.73e-4, contact theta^2 ln theta) are NOT a coarse error: "
        "E2X/E2Xs are against Dirac as in aggfull; E2B is against m_kin x Dirac",
        "the Born tail uses the asymptotic k^-4 form (corrections O((Z/a0)^2/k^2) at the cube surface, < 1e-3 at a0/8)",
    ],
    parameters="no new free parameter: twists fixed (0, pi), cube constants are geometry, Born tail from m, Z alpha",
)
json.dump(P, open(os.path.join(OUT, "agg2s_prereg.json"), "w"), indent=1)
print(json.dumps(P["predictions_from_development"], indent=1))

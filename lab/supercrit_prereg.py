"""supercrit_prereg.py - pre-registration of ATOMO SUPERCRITICO (docs/SUPERCRITICO.md).
Written and committed BEFORE the production runs (lab/supercrit_run_*.py). Writes lab/out/supercrit_prereg.json.

Development already SEEN before this file (scratch scripts, declared):
  * SWalk == RadialWalk to 0.0 for the point nucleus; Z = 60, a = 0.05: omega_1s = 0.89857 (Dirac 0.89905);
  * point nucleus: Z = 92 a = 0.01 / 0.001 relative binding error -2.7e-3 / -9.1e-5; Z = 137 a = 0.01 / 0.001:
    +0.44 / +0.25 (no convergence: the a^(2 gamma) law with 2 gamma = 0.046); Z = 150 a = 0.01 / 0.001:
    omega_1s = +0.162 / -0.187 (strong a dependence above Z alpha = 1);
  * shell nucleus R = 0.031 (Godunov et al.), Z = 170: omega = -0.6498, -0.6545, -0.6589, -0.6609 for
    a = 0.008, 0.004, 0.002, 0.001 (literature -0.664);
  * Z = 190 shell: continuum ODE phase scan E = -1.68487, Gamma = 6.12e-3 (literature -1.685, 6.06e-3); walk phase
    scan a = 0.004: (-1.6645, 5.58e-3), a = 0.002: (-1.6742, 5.81e-3);
  * continuum shell-model Z_cr(R): R = 1e-1 ... 1e-6 -> 205.22, 163.57, 150.42, 144.91, 142.17, 140.63
    (zeta_cr - 1 vs pi^2/(2 ln^2 R): 0.0575 vs 0.0582 at R = 1e-4; 0.0262 vs 0.0259 at 1e-6);
  * continuum shell Z_cr(R = 0.031) = 177.463.
NOT seen: everything else (all production a-series, Z_cr(a) of the point nucleus, uniform sphere, widths at
Z != 190, time-domain decay, spectra).
"""
import json, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

GODUNOV = {   # Godunov, Machet, Vysotsky, EPJC 77, 782 (2017), Table 1, shell R = 0.031/m, 'exact' columns:
    # Z: (Re eps, Gamma = 2 Im eps) in units of m c^2
    160: (-0.296, 0.0), 170: (-0.664, 0.0), 175: (-0.883, 0.0), 177: (-0.978, 0.0), 178: (-1.026, 0.0),
    180: (-1.127, 1.846e-6), 181: (-1.178, 1.895e-5), 182: (-1.231, 8.336e-5), 183: (-1.284, 2.366e-4),
    184: (-1.338, 5.182e-4), 185: (-1.394, 9.588e-4), 186: (-1.450, 1.581e-3), 187: (-1.507, 2.396e-3),
    188: (-1.565, 3.414e-3), 189: (-1.625, 4.636e-3), 190: (-1.685, 6.060e-3), 195: (-2.001, 1.606e-2),
    200: (-2.341, 3.034e-2), 210: (-3.094, 6.928e-2), 232: (-5.057, 1.928e-1)}

P = {
    "study": "supercrit", "doc": "docs/SUPERCRITICO.md",
    "rule": "zig-zag radial rule (lab/zz_core.py RadialWalk, frozen) with V(r) as a phase per tick; "
            "lab/supercrit_core.py SWalk only replaces V(r); kappa = -1 (s1/2)",
    "inputs": "m_e (theta = m a), alpha, Z. DEBTS: nuclear radius R (shell R = 0.031 lambda_C for the literature "
              "comparison; uniform sphere R = 1.2 (2.5 Z)^(1/3) fm, A = 2.5 Z because no nucleus with Z > 118 exists)",
    "free_parameters": 0,
    "error_models": {
        "point": "E(a) = E0 + K a^(2 gamma) + B a^2, gamma = sqrt(1 - (Z alpha)^2), fitted on the 3 finest a",
        "extended": "X(a) = X0 + c1 a + c2 a^2 on the finest 3 or 4 a (the kink of V at R gives order a)"},
    "configs": {
        "V1": {"Z": [20, 60, 92, 110, 120, 130, 137], "a": [0.5 / 2**k for k in range(4, 11)],
               "rmax": "max(30, 30/(Z alpha))"},
        "C1": {"a": [0.1, 0.03, 0.01, 0.003, 0.001, 3e-4, 1e-4, 3e-5], "rmax": 40,
               "method": "secant on omega_1s(Z) + 1 = 0 (shift-invert eigenvector of W, localisation checked)"},
        "C2": {"R": 0.031, "a": [0.004, 0.002, 0.001, 0.0005, 0.00025], "rmax": 40},
        "C3": {"R_fm": "1.2 (2.5 Z)^(1/3)", "a": [0.004, 0.002, 0.001, 0.0005, 0.00025], "rmax": 40},
        "G1": {"R": 0.031, "Z": [180, 182, 184, 186, 188, 190, 195, 200, 210, 232],
               "a": [0.004, 0.002, 0.001, 0.0005], "phase_site_r": 25.0,
               "method": "outward recursion of W chi = lambda chi, Bloch decomposition at r_m, Breit-Wigner fit "
                         "with quadratic background"},
        "G2": {"Z": [190, 200, 210, 232], "a": 0.002,
               "method": "initial state = eigenstate of the walk in a box closed at the middle of the Coulomb "
                         "barrier r_w = sqrt(r_in r_out) with r_in = zeta/(|E|+1), r_out = zeta/(|E|-1); evolve in "
                         "a large box with sponge; Gamma from the slope of ln P(r < r_w) between P = 0.7 and 0.05"},
        "G3": {"Z": [180, 181, 182, 183, 184, 185, 186, 187, 188], "a": "a -> 0 extrapolated values of G1 "
               "(extra Z 181..187 added at the finest two a only, with the same method)"},
        "D1": {"Z": 200, "a": 0.002, "detector_r": 30.0},
    },
    "exams": {
        "V1a": {"what": "point nucleus, Z <= 110: extrapolated 1s vs Dirac formula, relative binding error",
                "threshold": "<= 1e-4"},
        "V1b": {"what": "Z = 120", "threshold": "<= 1e-3"},
        "V1c": {"what": "Z = 130, 137: prediction = LIMIT (error at the finest a > 1e-2 for Z = 137; the "
                        "extrapolation is reported, threshold <= 1e-2 for Z = 130 and 137)",
                "threshold": "<= 1e-2 (expected to FAIL at 137)"},
        "C1a": {"what": "point nucleus beyond Z alpha = 1: the 1s reaches -mc^2 at a finite Z_cr(a) for every a "
                        "(the lattice is a cut-off); Z_cr(a) decreases monotonically as a decreases and stays "
                        "> 1/alpha", "threshold": "true for all 8 a"},
        "C1b": {"what": "effective nucleus: R_eff(a) such that the continuum shell model has Z_cr(R_eff) = "
                        "Z_cr(a); R_eff/a constant", "threshold": "max/min - 1 <= 0.10 over the 4 finest a; "
                        "prediction R_eff/a in [0.2, 2]"},
        "C2a": {"what": "shell R = 0.031: Z_cr walk (a -> 0) vs literature 177.46 (Godunov table, interpolated "
                        "between Z = 177 and 178; our continuum ODE 177.463)", "threshold": "|dZ| <= 0.3"},
        "C2b": {"what": "shell, Z = 170: E_1s (a -> 0) vs -0.664", "threshold": "relative <= 0.5 %"},
        "C3": {"what": "uniform sphere R = 1.2 (2.5 Z)^(1/3) fm: Z_cr (a -> 0)",
               "threshold": "in [169, 175] (literature ~170-173) AND |walk - continuum ODE| <= 0.3"},
        "G1a": {"what": "Gamma (a -> 0) vs Godunov exact, Z >= 182", "threshold": "|dGamma/Gamma| <= 10 %"},
        "G1b": {"what": "Gamma (a -> 0) vs Godunov exact, Z = 180", "threshold": "<= 30 %"},
        "G1c": {"what": "resonance energy (a -> 0) vs Godunov exact, all Z", "threshold": "relative <= 1 %"},
        "G2": {"what": "Gamma from the time decay vs Gamma from the phase scan, same a",
               "threshold": "<= 5 % (Z = 190, 200, 210), <= 10 % (Z = 232)"},
        "G3": {"what": "near-threshold Gamow law: ln Gamma = c0 - s 2 pi eta, eta = zeta |E_r| / sqrt(E_r^2 - 1); "
                       "fit over Z = 180..188", "threshold": "|s - 1| <= 0.20"},
        "D1": {"what": "Z = 200: spectrum of the outgoing wave at r = 30 (positron energy): peak vs |E_r| of the "
                       "phase scan, FWHM vs Gamma", "threshold": "peak within 1 %, FWHM within 15 %"},
    },
    "literature": {
        "Godunov_2017": "S.I. Godunov, B. Machet, M.I. Vysotsky, Eur. Phys. J. C 77, 782 (2017), arXiv:1707.07497",
        "table": {str(k): v for k, v in GODUNOV.items()},
        "Zcr_sphere": "~170 (Popov R_min = 1.2 (2.5Z)^(1/3) fm, e.g. arXiv:2209.00125), ~173 (Maltsev et al. PRL "
                      "123, 113401 (2019), Popov et al.)",
        "point_cutoff_law": "zeta_cr - 1 -> pi^2 / (2 ln^2 R) for R -> 0 (Zeldovich & Popov, Sov. Phys. Usp. 14, 673 "
                            "(1972)); checked here with our own continuum ODE (dev, declared)"},
}

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "supercrit_prereg.json"), "w") as f:
        json.dump(P, f, indent=1)
    print("written", os.path.join(OUT, "supercrit_prereg.json"))

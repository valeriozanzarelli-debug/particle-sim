"""box_prereg.py - pre-registration of the 'pairs in a box' study (A1-A4), written BEFORE the runs.
Writes lab/out/box_prereg.json."""
import json, os, subprocess, datetime

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
P = {
    "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head},
    "question": "Valerio: if the gamma photon passes between two plates very close together, does the e+e- pair "
                "fail to form because the electron's wave has no room?",
    "rule": "the R1 zig-zag walk of lab/redo_r1.py (x-hop R +1 / L -1 with Peierls phase, coin cos th, -i sin th; "
            "mass m = th), extended by ONE transverse direction y of d sites: per tick U = C(th) . W_y . S_x, with W_y "
            "a hop of one site in y conditioned on the sigma_y eigenvalue (+1 up, -1 down), so that in the bulk "
            "U ~ exp(-i(th sx + kx sz + ky sy)) (2+1 D Dirac, c = 1 in both directions). Hard walls: a ball that would "
            "leave the strip is reflected into the opposite y-mover on the same site with phase i (unitary). Field "
            "along x: Peierls phase A(t) on the x hop. Units hbar = c = a = tick = e = 1. th = 0.1.",
    "engine_fact_seen_before_prereg": "declared honestly: while building the engine I diagonalised U(kx=0) for the "
            "four reflection phases (1, i, -1, -i). Phases 1 and -1 leave a gapless mode (no confinement mass: rejected, "
            "not a hard wall for a massive particle); phases +-i give a half-gap m_eff(d) that matches "
            "sqrt(th^2 + (pi/(2d))^2) for d = 2..128, i.e. the lowest transverse momentum is pi/(2d) (Dirac 'bag' value), "
            "NOT the scalar-Dirichlet pi/d. So A1b (the formula as stated in the task, with pi/d) is expected to FAIL by "
            "a factor 2 in the transverse momentum; it is kept because the task asked for it.",
    "items": {
        "A1a_rate_vs_meff": "uniform field E along x; the whole sea of the strip (all lower states per kx) evolved "
            "exactly for kx from -1.2 to +1.2 (2.4/E ticks). Pairs = probability ending in states with quasi-energy in "
            "(0, pi/2) starting from states in (-pi/2, 0) (pi-gap doublers counted separately). Prediction: "
            "P_lat(d) = sum_j exp(-pi m_j^2/E), m_j = half-gaps at kx = 0 of the strip spectrum (lattice, measured). "
            "d in {4,6,8,12,16,24,32,64,128}; E = x m_eff(d)^2, x in {1,2,4}. PASS if |P_lat/P_pred - 1| <= 0.10 for "
            "every (d, x).",
        "A1b_meff_formula_as_stated": "m_eff(d) measured vs sqrt(th^2 + (pi/d)^2): PASS if within 10% for every d. "
            "(Expected FAIL, see engine fact.) Reported, not graded: vs sqrt(th^2 + (pi/(2d))^2).",
        "A1c_suppression_fixed_E": "fixed E = 2 th^2: P_lat(d) vs d, against the 1+1 D value exp(-pi th^2/E) = 0.208. "
            "PASS if P drops monotonically as d decreases AND at d <= 8 (d th <= 0.8) the suppression vs 1+1 D is at "
            "least a factor 10.",
        "A2_threshold": "photon = weak oscillating field A(t) = a0 sin(Omega t), a0 = 0.01, along x, on for T = 600 "
            "ticks; excitation summed over a kx grid. Absorption spectrum P(Omega), grid step 0.0025. Threshold Omega* = "
            "position of the absorption maximum (1+1 D van Hove edge just above the gap). PASS if Omega* = 2 m_eff(d) "
            "within 10% AND Omega*(d) rises above 2 th as d decreases; d in {8,16,32,128}.",
        "A3_crossover": "d* = width where m_eff(d*) = sqrt(2) th (transverse momentum = rest mass), interpolated. "
            "PASS if 0.5 <= d* th <= 5 (crossover at the Compton length 1/th).",
        "A4_formation_length": "1+1 D (R1 rule), Sauter step V(x) = (V0/2)(1+tanh(x/ell)) at FIXED peak field "
            "E = V0/(2 ell) = 2 th^2, ell from 2 to 160 sites. Yield Gamma(ell) = int dw/2pi N(w) (lattice transfer "
            "matrix, lattice_N of redo_r1) vs local-constant-field estimate Gamma_LCFA = int dx (E(x)/2pi) "
            "exp(-pi m^2/E(x)). Formation length ell_f = 2m/E = 10 sites. PASS if (a) Gamma/Gamma_LCFA within 15% of 1 "
            "for ell >= 8 ell_f AND (b) Gamma/Gamma_LCFA <= 0.5 for ell <= ell_f (exactly 0 when V0 < 2m). "
            "Reported, not graded: same with m -> m_eff of the strip d = 16.",
    },
    "no_fitted_parameters": "th, E, d, ell, a0 are inputs; nothing is adjusted.",
}
json.dump(P, open(os.path.join(OUT, "box_prereg.json"), "w"), indent=1)
print("written", os.path.join(OUT, "box_prereg.json"))

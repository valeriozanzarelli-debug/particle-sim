"""pre_prereg.py - PRE-REGISTRATION of the pre-experiments P1..P5 (written and committed BEFORE the runs).

Calibration done before this file (declared): one smoke test of pre_core.py (torus = geo_core rule to 4e-17,
unitarity with random link weights 1e-15, eigs timing n = 16..48, 200 ticks of |psi|^2-start surfers on a
24x24 box with M = 16: H-bar stayed at ~0.0009 = the noise floor). Nothing else was run.
Writes lab/out/pre_prereg.json with the git hash and the time.
"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

COMMON = {
    "rule": "zig-zag walk of zz_core/geo_core on an n x n lattice, fixed axis order x then y, flip theta = 0.5, "
            "4 components; open box = mirror edges (pre_core.LatticeWalk, link weights w; w = 1 = plain rule)",
    "psi": "sum of the M eigenmodes of the one-tick U closest to omega0 = 0.9 (positive band), random phases, seed 1",
    "surfers": "N = 20000, exactly equivariant local jump process per axis sub-step (pre_core.flows), "
               "lam = 'w' (default: lam = link weight; w = 1 -> ballistic); variant lam = 0 (Bell minimal flow)",
    "Hbar": "coarse-grained KL(surfers || |psi|^2) on cell x cell blocks; floor F = (cells - 1)/(2N)",
    "tau": "first time (H - F) <= (H0 - F)/e",
    "windings": "phase-winding plaquettes per spinor component on each sheet, mean of the 4 components "
                "(proxy: a 4-component spinor has no generic nodes of |psi|^2)",
}

PREREG = {
    "common": COMMON,
    "P1": {
        "setup": "box n = 24, cell 4 (36 cells, F = 0.000875), T = 400 ticks, M = 16 and M = 1 (control)",
        "P1a_floor": "start |psi|^2 (M = 16): max_t H <= 3F for lam = 'w' and lam = 0",
        "P1b_relax": "wrong starts (uniform; left half of the box) with M = 16: H(T)/H(0) < 0.10 (lam = 'w')",
        "P1c_control": "M = 1, uniform start: the continuum expectation (no vortices -> frozen) is "
                       "H(T)/H(0) > 0.9. OUR PREDICTION: this FAILS on the lattice, because the jump process is "
                       "stochastic: a Markov chain mixes even for a stationary psi (fake mixing).",
        "P1d_sheets": "start uniform on sheet (x+y) even only, M = 16. Sheet-resolved H (cells x sheet). PREDICTION: "
                      "with lam = 1 (ballistic, every surfer moves diagonally each tick) the sheet populations are "
                      "frozen -> sheet-resolved H(T) >= 0.5 H(0) (lattice artefact that BLOCKS mixing); "
                      "with lam = 0 surfers can change sheet -> relaxes (H(T)/H(0) < 0.1).",
    },
    "P2": {
        "pipeline": "geo_core blind reconstruction unchanged (components, -log(I/I0), Isomap, ball dimension, "
                    "Rips homology maxdim 2)",
        "positive_control": "vacuum MI of the rule on a 12x12 torus: 2 pieces, each b = (1, 2, *) ",
        "nulls": {
            "shuffle": "torus 12x12 vacuum MI, off-diagonal values randomly permuted (symmetric), diag kept; 20 draws",
            "random_graph": "random 4-regular graph, 144 sites, one fermion mode per site, H = -A, half filling, "
                            "MI of the ground state; 20 graphs",
            "noise": "symmetric matrix, log I off-diagonal Gaussian with the mean/std of log I (I > floor) of the "
                     "torus vacuum MI, same fraction below floor set to 0; diag = torus diag; 20 draws",
            "random_state": "Haar-random Gaussian state, 144 modes, 72 filled; 10 draws",
        },
        "false_positive": "on the largest piece: b1 + b2 >= 1 ('a hole or cavity'); reported separately: "
                          "D_small in [0.8, 3.5] ('looks low-dimensional')",
        "pass": "false-positive rate (holes) <= 10 % for every null class",
    },
    "P3": {
        "rule_fixed_in_advance": "every K ticks: J_l = number of surfer crossings of link l in the window; "
                                 "w <- (1 - eta) w + eta * wbar * J_l / mean(J); clip to [w_min, 1]; rescale "
                                 "mean(w) to wbar (3 clip/rescale passes). Box edge links stay 0.",
        "parameters": {"eta": 0.05, "wbar": 0.8, "w_min": 0.02, "K": 10, "count": 4},
        "setup": "box n = 12, cell 3 (16 cells), M = 16, N = 20000, T = 3000 ticks, lam = 'w', start w = wbar "
                 "everywhere, psi from the modes of U(w_start); same seeds with eta = 0 (no reinforcement)",
        "P3a_pinning": "'positive feedback pins non-equilibrium' if mean H over the last 10% of ticks with "
                       "reinforcement > 3 x the same without AND > 3F (uniform start)",
        "P3b_highways": "Gini(w) at the end > Gini(start) + 0.2 AND top 10% of links carry > 30% of the weight",
        "P3c_vortices": "winding count vs time reported (no pass/fail)",
        "P3d_shape": "geo_core reconstruction of the vacuum MI of U(w_end) vs U(w_start): nontrivial structure if "
                     "pieces or b1 or b2 change AND P2 false-positive rate <= 10 %",
        "P3e_equivariance": "|psi|^2 start WITH reinforcement: max H <= 3F (feedback must not create disorder)",
    },
    "P4": {
        "P4a": "one-tick Bloch U1(k) (fixed order): ||U1(k+(pi,pi)) - U1(k)|| < 1e-12 (sheet = exact copy) and "
               "U1(k+(pi,0)) = -U1(k) (pi-modes = temporal doubler)",
        "P4b": "alternating order (two ticks, zz_core.bloch_batch): still (pi,pi)-periodic -> alternation does NOT "
               "remove the sheets (prediction)",
        "P4c_fix": "Wilson-like coin in k-space: theta(k) = theta + r (2 - cos kx - cos ky), r = 0.5 (quasi-local, "
                   "Bessel tails). Pass if the vacuum MI on a 16x16 torus is ONE piece with b0, b1 = 1, 2",
        "P4d_labels": "cone windings at theta = 0 (R block): prediction all four Dirac points (0,0),(pi,pi),(pi,0),"
                      "(0,pi) have the SAME winding, L block opposite (no Nielsen-Ninomiya opposite-chirality "
                      "partner); the two sheets are related by a one-site translation: vacuum MI of sheet B "
                      "shifted by (1,0) equals sheet A to < 1e-10 -> no charge-like label",
    },
    "P5": {
        "setup": "uniform start, lam = 'w', tau vs n in {16, 24, 32} (M = 16) and vs M in {1, 4, 16, 36} (n = 24); "
                 "T = 400 (n <= 24), 600 (n = 32)",
        "report": "exponent a in tau ~ n^a (ballistic ~1, diffusive ~2); no pass/fail",
    },
}

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    try:
        git = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        git = "unknown"
    PREREG["_meta"] = {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head": git}
    with open(os.path.join(OUT, "pre_prereg.json"), "w", encoding="utf-8") as f:
        json.dump(PREREG, f, indent=1, ensure_ascii=False)
    print("pre-registration written", PREREG["_meta"])

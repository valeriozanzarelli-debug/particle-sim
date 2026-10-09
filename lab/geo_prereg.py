"""geo_prereg.py - PRE-REGISTRATION of every 'shape from the balls' test (written and run BEFORE the exams).

What was done before this file (calibration, declared): pilots on SMALLER shapes (ring 30, chain 30,
torus 10x10 and 14x14, line 20 x circle 9, line 30 x circle 12, stadium 20x12 and 24x14) to fix theta,
the maps and the thresholds below. The calibration found two things that the exams must now face:
  (1) DOUBLING. Every component of the rule moves exactly one site per axis per tick, so in 2-D the
      vacuum correlates only sites whose x AND y differences are both odd: the rule's space is two
      interleaved 'sheets' (x+y even / odd). On periodic directions the sheets stay separate; an open
      edge (the mirror) swaps sheet, so open edges GLUE the two sheets. In the quench (random product
      state) correlations need a common source: differences both even -> 2 sheets in 1-D, 4 in 2-D,
      glued in pairs at open edges.
  (2) The vacuum MI is exactly zero at even distances (1-D) - harmless after the geodesic step.
So each exam has two references: the shape we DREW (intended) and the shape the RULE connects
(the 'doubled' prediction, from mechanism (1)). Both are written here before the runs.
Writes lab/out/geo_prereg.json with the git hash and the time.
"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

PREREG = {
    "partA": {
        "rule": "zz_core zig-zag walk, 4 components (L0,L1,R0,R1), fixed axis order x then y (debt: no "
                "alternation), flip theta = 0.5, open edge = mirror (moving component -> opposite one, sign -1)",
        "theta": 0.5,
        "vacuum": "Dirac sea: fill single-ball modes with quasi-energy in (-pi, 0)",
        "quench": "random local product state (Haar 2-of-4 per site), seeds 2, MI averaged over seeds",
        "maps": {
            "I_floor": 1e-10,
            "vacuum_distance": "d = -log(I / I0), I0 = mean 2 S_i",
            "quench_distance": "d = first tick (log-interpolated) with I > 1e-3",
            "pieces": "connected components of the graph I > I_floor (quench: ever arrived)",
            "geodesic": "edges d <= 1.5 tau, tau = smallest connecting threshold; shortest paths",
            "dimension": "D_small = slope log N(r) vs log r on [1.5h, 4h]; D_large on [diam/4, diam/2]",
            "topology": "Rips (ripser) on geodesics; b0,b1 on all points, bar counts if persistence > 1.5 h; "
                        "b2 secondary on 220 max-min landmarks, persistence > 1.5 max(h, 2 eps_L)",
            "hidden_radius": "R = 3 * death(first-dying significant H1 bar) / (2 pi); visible extent = diam",
        },
        "exams": {
            "chain48":   {"intended": {"pieces": 1, "b": [1, 0], "D": 1},
                          "doubled_vacuum": {"pieces": 1, "b": [1, 0]},
                          "doubled_quench": {"pieces": 1, "b": [1, 1], "why": "even/odd chains glued at both ends = loop"},
                          "T": 34},
            "ring48":    {"intended": {"pieces": 1, "b": [1, 1], "D": 1, "R_over_extent": 0.3183},
                          "doubled_vacuum": {"pieces": 1, "b": [1, 1]},
                          "doubled_quench": {"pieces": 2, "b": [1, 1]},
                          "T": 20},
            "torus16":   {"intended": {"pieces": 1, "b": [1, 2, 1], "D": 2},
                          "doubled_vacuum": {"pieces": 2, "b": [1, 2, 1]},
                          "doubled_quench": {"pieces": 4, "b": [1, 2, 1]},
                          "T": 12},
            "cyl36x18":  {"intended": {"pieces": 1, "b": [1, 1, 0], "D_small": 2, "D_large": 1,
                                       "R_over_extent": (18 / (2 * 3.141592653589793)) / 35},
                          "doubled_vacuum": {"pieces": 1, "b": [1, 2], "why": "sheets glued at both ends: loop 2L x hidden circle"},
                          "doubled_quench": {"pieces": 2, "b": [1, 2], "note": "hidden circle = 9 sheet sites, marginal"},
                          "T": 26},
            "stadium30x16": {"intended": {"pieces": 1, "b": [1, 1, 0], "D": 2},
                             "doubled_vacuum": {"pieces": 1, "b": [1, 2], "why": "double of an annulus = torus"},
                             "doubled_quench": {"pieces": 1, "b": [1, 2]},
                             "T": 30},
        },
        "pass": {
            "fidelity": "Spearman(reconstructed geodesic, causal-graph distance of the rule) >= 0.9; "
                        "Spearman vs the drawn Euclidean distance reported",
            "betti": "b0, b1 exactly as predicted (both references reported)",
            "dimension": "round(D_small) = drawn dimension; cylinder: round(D_large) = 1",
            "radius": "|(R/diam)_est / (R/L)_true - 1| < 0.15 (ring: R/diam vs 1/pi; cylinder: vs (18/2pi)/35)",
        },
    },
    "partB_synthetic": {
        "true_shape": "flat oblique torus, sides a = 1, b = 1.37, angle 63 deg; base mass m0 = 0.35*(2 pi)",
        "masses": "E = sqrt(m0^2 + |k|^2), k in the dual lattice (2 pi units), sorted with multiplicity, 40 levels",
        "noise": "relative Gaussian 1e-3 on every mass, seed 20261009",
        "families": {"F1 square torus": 2, "F2 rectangle": 3, "F3 oblique torus (true family)": 4,
                     "F4 shape-free polynomial E^2 = sum_p c_p k^(p/2)": "P = 2..8"},
        "train_test": "train on levels 1..k (k = 3..20), predict k+1..40; score = median relative error",
        "expect": "F3 < 1e-2 for k >= 6; F1, F2 biased (fail); F4 with P >= k fails or degrades with distance",
    },
    "partB_isospectral": {
        "construction": "Gordon-Webb-Wolpert-type pair from 7 half-squares glued by three involutions of "
                        "GL(3,2) acting on points (drum 1) and on lines (drum 2) of the Fano plane "
                        "(Sunada / Buser-Conway-Doyle-Semmler); first triple found with both gluing graphs "
                        "trees and both tilings non-overlapping",
        "fem": "P1 elements, each half-square refined into n^2 congruent triangles (n = 16, 32)",
        "expect": "first 10 Dirichlet eigenvalues equal to < 1e-8 relative at the SAME n; couplings "
                  "g3 = int u1^3, g4 = int u1^4, |int u1 u2 u3| differ by more than 10x their n=16->32 change",
    },
    "partB_real": {
        "rule": "use measured masses for training only if the family has FEWER parameters than training "
                "numbers + held-out numbers; held-out numbers are the exam",
        "tests": {
            "KK2": "hidden 2-torus: charge circle (n1 = 1 for all charged leptons) x generation circle n2 = 0, 1, 2; "
                   "2 parameters (M, 1/R2) trained on m_e, m_mu; exam: m_tau and the absence of the n2 = -1 partner",
            "KOIDE_lep": "three points on a circle, sqrt m_k = M (1 + sqrt2 cos(delta + 2 pi k / 3)); 2 parameters "
                         "trained on m_e, m_mu; exam: m_tau (PDG 1776.86 +- 0.12 MeV), pass if within 3 sigma; "
                         "the larger root is taken",
            "KOIDE_up": "same family, new (M, delta): train on m_u, m_c, exam m_t; pass if within 10 %",
            "KOIDE_down": "same family: train on m_d, m_s, exam m_b; pass if within 10 %",
        },
        "warning": "the sqrt2 amplitude of the Koide family was chosen in 1981 knowing the leptons: a hidden "
                   "parameter (look-elsewhere); a pass on leptons is weak evidence, a fail on quarks counts",
    },
}

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    try:
        git = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        git = "unknown"
    PREREG["_meta"] = {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head": git}
    with open(os.path.join(OUT, "geo_prereg.json"), "w", encoding="utf-8") as f:
        json.dump(PREREG, f, indent=1, ensure_ascii=False)
    print("pre-registration written:", os.path.join(OUT, "geo_prereg.json"), PREREG["_meta"])

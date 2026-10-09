"""flow_prereg.py - PRE-REGISTRATION of phase L2 (the flow). Writes lab/out/flow_prereg.json. Committed together with
eye_prereg.json, BEFORE any L1 or L2 run.

DEV (declared, looked at): uniform hex lattice has Ollivier kappa = 0 exactly with alpha = 0.5 (flat = kappa 0 here);
Euler characteristic checks of the meshes (hex 0, torus Delaunay 0, sphere hull 2) on N = 60 / L = 6; 8 flow steps on a
L = 6 hex torus with 40 % perturbation (dev seed 5) with eps = 0.25: k_std 0.096 -> 0.047, w_std 0.24 -> 0.18 (slow),
2 ms per W1 linear program -> eps raised to 0.5 and MAXIT to 100 for the time budget.
"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
S0 = 20261009

L2 = {
    "flow": {
        "object": "triangulated surface; link LENGTHS w (Ni-Lin-Gao-Gu 2019 convention)",
        "transport": "lazy walk: stay alpha = 0.5, hop to neighbour z with prob (1-alpha) exp(-w_xz^2)/sum (p = 2)",
        "curvature": "Ollivier-Ricci kappa = 1 - W1(m_x, m_y)/w_xy, exact W1 by linear programming, ground cost = "
                     "shortest paths with the current lengths",
        "step": "w <- w - eps * kappa * w, eps = 0.5; then rescale to the initial total length (normalised flow); "
                "no surgery",
        "stop": "max_e |dw_e|/w_e < 1e-3, or 100 steps (whichever first)",
        "parameters": {"eps": 0.5, "alpha": 0.5, "tol": 1e-3, "maxit": 100, "p_exponent": 2},
        "eye_hopping": "balls for the eye hop with T_e = exp(1 - w_e^2) (same kernel), thermal beta 1 mu 0 (L1 eye, "
                       "THR2 from L1 calibration)",
        "holonomy": "nearest feasible check of the orb_core spinor test (that one needs a lattice symmetry; a weighted "
                    "graph has none): Regge geometry of the triangulation, deficit delta_v = 2pi - sum of corner angles "
                    "(law of cosines from the lengths). Spinor around v multiplied by exp(i sigma_z delta_v/2). "
                    "N_inv = 2 if every |delta_v| < 0.05 rad, else 0. Cycle holonomy of the torus not computed (for an "
                    "exactly flat closed orientable surface it is trivial, Bieberbach).",
    },
    "starts": {
        "a": "hex torus L = 15 (225 v, 675 e), w = 1 + 0.4 U(-1,1), seed default_rng(%d)" % (S0 + 41),
        "b": "random geometric triangulation of a flat 60-deg torus: 225 uniform points, periodic Delaunay, ALL lengths "
             "start at 1 (only the links are kept, the drawn positions are forgotten); seeds %d, %d. Deviation from "
             "'radius graph': a radius RGG is not a surface mesh, Regge holonomy needs triangles." % (S0 + 42, S0 + 43),
        "c": "sphere: 225 uniform points on S^2, convex hull triangulation, all lengths 1; seeds %d, %d" % (S0 + 44, S0 + 45),
        "z3": "hex torus L = 15, w_e = 1 + 0.2 f(midpoint), f(x) = sum_k cos(b_k . x), b_k the 3 shortest reciprocal "
              "vectors of the L-torus (120 deg related): extrema exactly at the 3 fixed points of the 120 deg rotation "
              "(0,0), (5,5), (10,10). No noise. Not a quotient.",
    },
    "exams": {
        "a_b_spread": "S = std of kappa over edges; PASS if S_final <= S_initial / 10",
        "a_b_eye": "final MI read by the eye: loops == 2 AND cones == 0 AND cavities == 1 (cavity graded 'if feasible': "
                   "reported separately)",
        "a_b_holonomy": "N_inv == 2 on the final lengths",
        "c": "PASS if mean kappa_final > 0 for both sphere seeds AND min over sphere seeds of mean kappa_final > max over "
             "torus (b) seeds of mean kappa_final (sphere rounder than the flat reference). Holonomy N_inv expected 0 "
             "(Gauss-Bonnet: sum delta = 4 pi always: a trivial consequence, not evidence). Eye: cavities 1, loops 0 "
             "reported.",
        "z3": "concentration C = sum_{v within 1 hop of the 3 fixed points} |delta_v| / sum_v |delta_v| (area baseline "
              "21/225 = 0.093); same for |vertex Ollivier curvature|. Verdict 'cone-like' if C_final >= 0.5 AND eye cones "
              "== 3 within 1 hop of the 3 fixed points; 'washed out' if max|delta|_final <= max|delta|_initial / 10; "
              "else 'neither'. Reported either way; no tuning.",
    },
    "go_no_go": "GO for a bigger emergence study only if L1 passes AND (a) passes all three exams AND z3 is 'cone-like'. "
                "Otherwise NO-GO with the reason.",
}


def meta():
    try:
        git = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        git = "unknown"
    return {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": git}


if __name__ == "__main__":
    with open(os.path.join(OUT, "flow_prereg.json"), "w") as f:
        json.dump({"meta": meta(), "L2": L2}, f, indent=1)
    print("written lab/out/flow_prereg.json")

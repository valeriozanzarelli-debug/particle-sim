"""eye_prereg.py - PRE-REGISTRATION of phase L1 (the eye). Writes lab/out/eye_prereg.json. Committed BEFORE eye_calib.py,
eye_nulls.py, eye_shapes.py run.

DEV SET (declared, outputs LOOKED AT while building the eye; none of these is an exam instance):
  * one null per class with seed default_rng(999) (timing + rel1/rel2/cone counts seen: shuffle gave 77 'cones' with the
    first naive one-scale detector -> led to the two-scale rule and the nearest-neighbour gate; max rel2 seen 0.18);
  * hex torus L = 12 and 15, T^2/Z3 from L = 12 (50 orbits), icosphere level 2 (162), T^4/Z2 n = 4 (136 orbits):
    rho profiles at 1.5..4 h and rel1/rel2 with 100 landmarks were seen (T^2/Z3-12: rho 0.43/0.37 at the 3 fixed points;
    icosphere 5-valent vertices: rho 0.86/0.84; T^4/Z2 n=4: 0.56 at 1.5 h but 0.92 at 2.5 h because the 2.5 h ball
    covers a quarter of that tiny space; flat torus rel2 0.38 (L=12) / 0.29 (L=15); sphere rel2 0.42).
  * ripser H2 timing: 150 landmarks took 15-80 s on dev nulls -> 100 landmarks (1.6-1.9 s).
Exam sizes are different from all dev sizes.
"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
S0 = 20261009

L1 = {
    "eye": {
        "input": "the MI matrix I only (one fermion mode per site); never the adjacency, coordinates or symmetry",
        "pipeline": "largest informational piece (I > 1e-10) -> d = -log(I/I0) -> geo_core.geodesic (factor 1.5) -> Gm, "
                    "diameter D = max Gm",
        "loops": "big_core H1 reader UNCHANGED: 220 max-min landmarks, ripser full filtration, loop = H1 bar with "
                 "persistence/D > 0.4",
        "cavities": "ripser H2 on 100 max-min landmarks (<= 150 as asked; 150 too slow, see dev), full filtration, "
                    "cavity = H2 bar with persistence/D > THR2. THR2 from calibration (below).",
        "cones": "h = median nearest-neighbour geodesic distance. For c in (1.5, 2.5): V_x(c) = #{y: Gm(x,y) <= c h} "
                 "incl. x; rho_x(c) = V_x(c)/median V(c); rho_bar = mean of the two. CONE POINT = rho_x(c) < 0.7 at BOTH "
                 "scales AND nearest-neighbour distance of x <= 1.2 h AND rho_bar(x) is the smallest in its 2.5 h ball "
                 "(ties: lowest index). Angle (2-D) = 2 pi rho_bar; solid-angle fraction = rho_bar. Expected: Z3 cone "
                 "(angle 2pi/3) rho ~ 1/3; 4-D Z2 point rho ~ 1/2; icosphere 5-valent vertex (angle 5/6 * 2pi) ~ 0.83.",
        "balls_on_drawn_shapes": "one fermion mode per site hopping on the drawn graph, H = -T, thermal Gaussian state "
                                 "C = 1/(1+exp(beta (H - mu))), beta = 1, mu = 0. Quotients: symmetric (untwisted) sector "
                                 "H_q = Q^T H Q, Q = normalised orbit indicator. Only MI(C) goes to the eye.",
    },
    "calibration": {
        "nulls": "4 classes x 10, seeds default_rng(%d) (fresh, never used): shuffle and noise of the R1 base (24x24 box "
                 "walker vacuum MI, theta 0.125, w 0.3), random_graph = random 6-regular graph on 576 nodes with the SAME "
                 "thermal balls as the drawn shapes (pipeline-matched; R1 used 4-regular ground state), random_state = "
                 "Haar Gaussian state 288/576" % (S0 + 31),
        "THR2_rule": "THR2 = max(0.2, 1.25 * largest rel2 over the 40 calibration nulls). Fixed before any exam null or "
                     "exam shape is read. Cone thresholds are NOT calibrated (fixed above).",
    },
    "exam_nulls": {
        "nulls": "same 4 classes x 25, seeds default_rng(%d) (fresh)" % (S0 + 32),
        "pass": "per class: fraction of instances with cavities >= 1 <= 5 % (<= 1/25) AND fraction with cones >= 1 "
                "<= 5 %. Loops (H1) false positives reported, not graded (graded in R1).",
    },
    "exam_shapes": {
        "T2/Z3": {"draw": "hex torus L = 27, 120 deg rotation, 245 orbits",
                  "pass": "cavities == 1 AND loops == 0 AND cones == 3, the 3 cones at hop distance <= 1 (orbit graph, "
                          "grading only) from 3 DISTINCT fixed orbits, each rho_bar in [0.2, 0.5] (angle 72..180 deg)"},
        "T4/Z2": {"draw": "hypercubic n = 6 torus, x -> -x, 656 orbits",
                  "pass": "cones == 16 at hop distance <= 1 from 16 distinct fixed orbits, each rho_bar in [0.35, 0.7]; "
                          "loops == 0",
                  "not_graded": "cavities (target b2 = 6). Declared expectation: H2 on 100 landmarks of a 656-site 4-D "
                                "space cannot resolve the six 2-cycles; H3, H4 out of budget (not computed)."},
        "T2_flat": {"draw": "hex torus L = 15 (225 sites)", "pass": "cones == 0",
                    "secondary": "loops == 2 and cavities == 1 (reported, graded separately as 'torus recognised')"},
        "S2": {"draw": "icosphere level 3 (642 sites; 12 five-valent vertices = mild cones of angle 300 deg)",
               "pass": "cones == 0", "secondary": "cavities == 1 and loops == 0"},
    },
    "L1_PASS": "exam nulls pass AND T2/Z3 pass AND T4/Z2 cone pass AND cones == 0 on both controls",
    "if_fail": "report plainly; no re-tuning of thresholds after any exam result",
}


def meta():
    try:
        git = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        git = "unknown"
    return {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": git}


if __name__ == "__main__":
    with open(os.path.join(OUT, "eye_prereg.json"), "w") as f:
        json.dump({"meta": meta(), "L1": L1}, f, indent=1)
    print("written lab/out/eye_prereg.json")

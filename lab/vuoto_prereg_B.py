"""vuoto_prereg_B.py - pre-registration of EXPLORATION B (written before its runs). Writes lab/out/vuoto_prereg_B.json.
Exploration, NOT the rule: flip angle th - mu on even ticks and th + mu on odd ticks (mean = th = the mass).
Development fact seen: this keeps the 4 massive species per pi-cell at mass th and gives the quasi-energy-pi cones
(massless in the rule) a mass mu (gap 2 mu per tick). Then the WHOLE sea (no species projector) is predicted to be
Pauli-Villars-regularized QED: physical species minus 'upside-down' regulators of mass mu/a."""
import json, os, datetime, subprocess
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = {"written": datetime.datetime.now().isoformat(timespec="seconds"),
     "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip(),
     "configs": {"theta": [0.05, 0.025], "mu": [0.3, 0.6], "directions": "x, y, z (average, per-direction subtraction)",
                 "q_over_m": [0.05, 0.025, 0.3, 1, 3], "quadrature": "chi_full_3d nseg = 16 (check 20)"},
     "exams": {
        "B1_uehling_from_whole_sea": "Pi_R,whole (xyz avg) / Pi_R,Uehling - 1 at q/m in {0.3, 1, 3}: PASS if within 2% at "
                                     "theta = 0.025 for both mu.",
        "B2_log_coefficient": "[pit0(theta=0.025) - pit0(theta=0.05)]/ln 2 at fixed mu: PASS if within 3% of 1/(6 pi^2) "
                              "for both mu.",
        "B3_regulator_mass_log": "[pit0(mu=0.6) - pit0(mu=0.3)]/ln 2 at fixed theta: predicted +1/(6 pi^2) "
                                 "(Pi(0) = (2 alpha/3pi) ln(M/m)). PASS if within 10% at both theta.",
        "B4_finite": "pit0 whole is finite and positive (screening) for every config.",
     },
     "status": "exploration only: a CHANGE of the rule would need Valerio's decision; not used for V1-V4."}
json.dump(P, open(os.path.join(OUT, "vuoto_prereg_B.json"), "w"), indent=1)

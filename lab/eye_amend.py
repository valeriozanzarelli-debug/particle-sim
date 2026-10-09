"""eye_amend.py - AMENDMENT to the L1 pre-registration, written AFTER the calibration nulls (eye_calib.py) and BEFORE
any exam null or exam shape was read. Writes lab/out/eye_amend.json.
Reason: on the 40 calibration nulls the pre-registered cone rule (v1) found 53-68 'cones' in EVERY shuffle null and
1-2 in 2 of 10 noise nulls. v1 would fail the exam nulls; it is still run and graded as registered.
v2 (eye_core.cones_v2) was then checked ONLY on the same calibration nulls: 0 cones in all 40 (shuffle refused,
f_low 0.15-0.18; noise cones removed by the third scale)."""
import json, os, datetime, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
A = {
    "THR2_from_calibration": "see eye_calib.json: max rel2 over 40 calibration nulls 0.190 -> THR2 = 0.237",
    "v1_status": "kept and graded exactly as pre-registered (expected to FAIL the shuffle class)",
    "cone_rule_v2": "scales (1.5, 2.5, 3.5) h: rho < 0.7 at ALL THREE, nn <= 1.2 h, local minimum of rho_bar in its "
                    "2.5 h ball; HOMOGENEITY GATE: if the fraction of sites meeting the low-volume condition exceeds 5 %, "
                    "the eye refuses to read cones (reports 0 and 'refused'). Angle = 2 pi rho_bar (mean of 3 scales).",
    "v2_grading": "same pass criteria as v1 (exam nulls <= 5 % per class; T2/Z3 3 cones rho_bar in [0.2,0.5] at distinct "
                  "fixed orbits; T4/Z2 16 cones rho_bar in [0.35,0.7]; controls 0). Reported side by side with v1. "
                  "L1 'PASS (v2)' is labelled as post-calibration, not as the original pre-registration.",
    "L2_uses": "the L2 eye readings use THR2 = 0.237 and report both cone rules",
}
if __name__ == "__main__":
    git = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    json.dump({"meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": git},
               "amendment": A}, open(os.path.join(HERE, "out", "eye_amend.json"), "w"), indent=1)

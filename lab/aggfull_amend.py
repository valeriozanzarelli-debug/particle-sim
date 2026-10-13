"""aggfull_amend.py - AMENDMENT to the pre-registration (written AFTER the first production results; declared).
Writes lab/out/aggfull_amend.json. The original exams and verdicts stay as they are; this adds new runs."""
import json, os, datetime, subprocess
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A = {
 "written": datetime.datetime.now().isoformat(timespec="seconds"),
 "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip(),
 "why": "F3 (alpha^4 shifts) is expected to FAIL for the s states: the post hoc diagnosis (aggfull_posthoc_kin.py) shows that "
        "the free generator of chain A carries a learning error at high k (~1e-3 relative at the top 1 % of the 1s momenta) "
        "that is 6x the relativistic -k^4/8 term at a0/8 (1s: <T - k^2/2m> = -5.97e-3 eV learned vs -0.99e-3 fine walk vs "
        "-0.98e-3 Dirac). Cause found in development tests: the 25 % of uniform samples over the Brillouin cube (kept for "
        "tameness) limits the fit to ~1e-4 relative, while the k^4 term is ~alpha^2/4 ~ 1e-5 relative.",
 "change": "Kinetic re-learned (debt #12, 'kinetic term not renormalised'): the same learner (agg_zz_core.learn_h), w = 4, "
           "nk = 8000, the uniform share kept (25 %) but with weight 1e-3 instead of 0.03 (a first try with no uniform share "
           "at all gave energy error 7e-7 but negative motion energies in the cube: rejected by the check; weights 3e-3..1e-4 "
           "tried on B = 256, 1e-3 chosen as the largest with energy error < 3e-6); level B = 256 from the fine walk, "
           "B = 512 from the re-learned B = 256 (repeat). Tameness is no longer in the fit but still CHECKED (no negative "
           "motion energy anywhere in the cube; agg_zz_core.kcheck state mode: packets, energy, bz_rel90, neg_gen): a rule "
           "failing the check is not used.",
 "new_runs": "atom jobs with suffix K (eK8_1s, eK8_2s, eK8_2p, eK16_1s; recoil pK8_*, pK16_1s) and their first-order jobs",
 "new_exams_same_thresholds": {
   "F3K": "alpha^4 shifts with the re-learned kinetic, within 5 % each (same threshold as F3)",
   "E1K/E2K": "levels vs Dirac within 0.1 % (same thresholds as E1/E2)",
   "R2K": "recoil shift ratio within [0.99, 1.01]",
 },
}
json.dump(A, open(os.path.join(OUT, "aggfull_amend.json"), "w"), indent=1)
print("written")

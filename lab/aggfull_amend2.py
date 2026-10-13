"""aggfull_amend2.py - SECOND amendment of AGGLOMERATORE COMPLETO (after aggfull_amend.json), written before the K atom runs.
Writes lab/out/aggfull_amend2.json"""
import json, os, subprocess, datetime
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A = {
 "written": datetime.datetime.now().isoformat(timespec="seconds"),
 "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip(),
 "why": "with the uniform weight 1e-3 of amendment 1 the re-learned B = 256 kinetic agrees with the fine walk on the 1s "
        "(a0/16) momenta (<T - k^2/2m> -7.44e-4 eV vs fine -7.65e-4, 3 %) and B = 512 on the 2s (-1.02e-4 vs -1.26e-4), but "
        "B = 512 on the 1s (a0/8) has the WRONG SIGN (+3.5e-3 eV vs fine -1.0e-3): the 1s at a0/8 lives up to the edge of "
        "the cell's zone, where 1e-3 of uniform weight still pulls the fit (aggfull_posthoc_kin.json).",
 "rule": "the uniform weight of each level is chosen by the LEARNER against its teacher (the fine walk), never against "
         "Dirac: the largest of 1e-3, 3e-4, 1e-4, 3e-5 whose state-distributed <T - k^2/2m_kin> agrees with the fine walk "
         "within 5 % (same samples), with the usual check passing (packets, energy, bz_rel90, no negative motion energy).",
 "dev_scan_B512": {"1e-3": "+9.2e-4 eV vs fine -1.29e-3 (30000 samples): rejected", "3e-4": "-1.316e-3 (2.3 %): chosen",
                   "1e-4": "-1.46e-3 (13 %)", "3e-5": "-1.320e-3 (2.6 %)"},
 "choice": {"B=256": 1e-3, "B=512": 3e-4},
 "exams": "unchanged: F3K, E1K, E2K, R2K of aggfull_amend.json with the same thresholds",
 "knobs": "one more hyper-parameter (the uniform weight) per level, chosen by a fixed rule; counted in the document",
}
json.dump(A, open(os.path.join(OUT, "aggfull_amend2.json"), "w"), indent=1)
print("written")

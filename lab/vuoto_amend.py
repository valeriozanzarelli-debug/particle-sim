"""vuoto_amend.py - declared amendment to lab/out/vuoto_prereg.json, written BEFORE any graded V3/V4 number was computed.
Writes lab/out/vuoto_amend.json."""
import json, os, datetime, subprocess
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A = {"written": datetime.datetime.now().isoformat(timespec="seconds"),
     "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip(),
     "A1_V3_V4e_smearing": "q_c = 0.3/theta -> 0.15/theta. Reason (found on the ANALYTIC Uehling alone, no lattice number "
        "seen): the lattice q-table stops at q_lat = 0.6 (= 0.6/theta in m units, the species cutoff); with q_c = 0.3/theta "
        "the damping there is only e^-4, and the abrupt end makes the sin(qr) integral ring: the analytic smeared Uehling "
        "at r = 3 lambda_C came out 0.015 x the true value instead of ~1. With q_c = qmax/4 (damping e^-16) the ringing is "
        "gone. The same smearing is applied to the reference, as pre-registered. Thresholds unchanged.",
     "A2_reference_values_seen": "references computed before the lattice runs (exam only): H 1S -215.10 MHz, H 2S "
        "-26.887 MHz (Schroedinger density + analytic Uehling, includes the (Z alpha)^5 finite-range term; leading formula "
        "-217.2 / -27.15), Nd59+ 1s -11.667 eV (Dirac 1s, leading formula -11.632)."}
json.dump(A, open(os.path.join(OUT, "vuoto_amend.json"), "w"), indent=1)
print(A)

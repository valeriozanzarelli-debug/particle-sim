"""vuoto_f2_amend_L2.py - declared AMENDMENT A-L2 to lab/out/vuoto_f2_prereg_L.json, written AFTER the pre-registered
Bethe-log run and after amendment A-L1 (both stay as they are; L1-L3 FAIL is the graded verdict). Writes
lab/out/vuoto_f2_amend_L2.json."""
import os, json
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A = {"amendment": "A-L2", "date": "2026-10-13", "status": "written after seeing the A-L1 results (declared)",
     "what_was_seen": "A-L1 (tighter projector quadrature) did NOT cure 2S at theta 0.05 (it got worse: ln k0(300 Ry) diff "
                      "+1.91 instead of +1.00), so the A-L1 diagnosis (projector leak) was wrong or incomplete. New diagnosis "
                      "(measured with eigs on the p-channel walk): the spectrum of W/lambda_n is mirrored (doubler partner at "
                      "phase +pi). The state degenerate with 2S (2P_1/2: |W/lambda - 1| = 1.6e-10 at theta 0.05, 1.4e-9 at "
                      "0.1) has a mirror at |W/lambda + 1| = 1.6e-10, i.e. at the pole of the Cayley T = i(W/l - 1)(W/l + 1)^-1: "
                      "any tiny weight of the source there is amplified ~1e10. For 1S the nearest is 1e-6 (harmless).",
     "change": "same method and numerics as A-L1, plus explicit removal (orthogonal projection, W unitary) of the 6 "
               "eigenvectors of W/lambda_n nearest to -1 from the projected source, before the Cayley transform. Same "
               "exams and thresholds; results in lab/out/vuoto_f2_bethe_A2.json, reported next to the pre-registered "
               "and A-L1 ones; NOT graded (the verdict stays L1, L2, L3 FAIL)."}
json.dump(A, open(os.path.join(OUT, "vuoto_f2_amend_L2.json"), "w"), indent=1)
print("written")

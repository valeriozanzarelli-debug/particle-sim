"""vuoto_f2_amend_L.py - declared AMENDMENT A-L1 to lab/out/vuoto_f2_prereg_L.json, written AFTER the pre-registered
Bethe-log run (which stays the graded result). Writes lab/out/vuoto_f2_amend_L.json."""
import os, json
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A = {"amendment": "A-L1", "date": "2026-10-13", "status": "written after seeing the pre-registered L results (declared)",
     "what_was_seen": "1S theta 0.05 agrees with NR at 30/100/300 Ry within 0.0063, but 2S theta 0.05 blows up (ln k0(300 Ry) "
                      "= 2.74 vs 1.73); D_K ratio 1S 0.968 (L3 fail). Diagnosis: the 2S dipole source has norm <r^2> ~ 8e5 "
                      "and the sign-function quadrature of the positive-energy projector (h = 0.5, tails e^-19 = 6e-9) leaves "
                      "a negative-energy leak whose (2m)^2-weighted contribution swamps the physical D ~ 1e-10 m.",
     "change": "same method, tighter numerics only: h = 0.35, tails 28 (quadrature error ~1e-12). Same exams and thresholds; "
               "results in lab/out/vuoto_f2_bethe_A.json, reported NEXT TO the pre-registered ones, which remain the "
               "graded verdict (L1, L2, L3 FAIL as pre-registered)."}
json.dump(A, open(os.path.join(OUT, "vuoto_f2_amend_L.json"), "w"), indent=1)
print("written")

"""nucleo1s_f7_prereg.py -> lab/out/nucleo1s_f7_prereg.json : PRE-REGISTRATION of phase 7 (before the runs).
(1) The phase-6 post-hoc setting becomes a real exam on charges never used: Z = 3 and Z = 5 (blind), outer cells
theta = 0.025 / 0.0125 / 0.00625, zoom 14 levels x 16 cells, rmax 45 a0/Z, quadratic extrapolation in theta^2.
Seen before (phase 6 post hoc, same setting): Z = 1: 1s +0.19 %, 2s +0.08 %; Z = 2: +0.10 %, +0.04 %; 3s -14 % / -3.5 %.
(2) Diagnosis of the 3s (explicit, each changes ONE thing, Z = 1, same theta set): box rmax 90 a0; zoom 18 levels;
16 -> 32 cells per level; a fourth theta point (0.003125) for the extrapolation."""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = dict(study="NUCLEO_1S phase 7", date="2026-10-13",
         config=dict(theta_out=[0.025, 0.0125, 0.00625], zoom_levels=14, n_per_level=16, rmax_a0=45, Z_exam=[3, 5],
                     diag=dict(box=dict(rmax_a0=90), levels=dict(zoom_levels=18), per_level=dict(n_per_level=32),
                               theta4=dict(theta_out=[0.025, 0.0125, 0.00625, 0.003125]))),
         exams=dict(
             B1=dict(what="Z = 3 and Z = 5 (blind): alpha^4 shift of 1s and 2s within 1 % of Dirac's", tol=0.01),
             B2=dict(what="Z = 3 and 5: the 1s and 2s shifts are better than 0.5 %", tol=0.005),
             B3=dict(what="3s diagnosis: at least one of the four single changes brings the Z = 1 3s shift within 2 %", tol=0.02)))
json.dump(P, open(os.path.join(OUT, "nucleo1s_f7_prereg.json"), "w"), indent=1)
print("written")

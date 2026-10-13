"""nucleo1s_f6_prereg.py -> lab/out/nucleo1s_f6_prereg.json : PRE-REGISTRATION of phase 6 of NUCLEO_1S (before the runs).
Goal (coordinator, 13/10): the alpha^4 shifts of the s states (Darwin + kinetic, = Dirac - Schroedinger) from the rule,
with fine zig-zag cells at the nucleus (cusp resolved) and coarser levels outside, for Z = 1 if it holds.
What is ALREADY SEEN (declared): with phase-1 data (zoom chain, outer cells theta = 0.05 / 0.025 / 0.0125, quadratic
extrapolation in theta^2) the Z = 1 shifts come out 1s +1.45 %, 2s +0.58 %, 3s -13.8 % from Dirac (computed today
from lab/out/nucleo1s_results.json, before this file).
Error budget that decides the architecture: the shift is 1.3e-5 of the 1s binding at Z = 1, so the rule must be right
to ~1e-7. The agglomerated coarse rules available outside (agg_zz_3d_A in the phase-1/2 hybrid: +0.039 %; aggfull
Taylor-fixed level 1: 3e-5 at Z = 10, docs/AGGLOMERATO_COMPLETO.md) are 100-1000x too coarse for this: gluing them
outside would bury the shift. So the outer levels here are the chain of levels of the SAME rule (graded zoom with delay
slots, phase 1: the agglomerator's 'block' applied to the rotations), and the coarse-rule hybrid is a design (doc §16).
Runs: radial zig-zag (kappa = -1), zoom 14 levels x 16 cells at the nucleus, outer theta in {0.1, 0.05, 0.025},
rmax = 45 a0/Z; Z = 1 and Z = 2 (He+, never run with the zoom); quadratic extrapolation in theta^2 through the three.
Shift_rule = B_rule - B_Schr, B_Schr = (Z alpha)^2 / (2 n^2) (exam reference only); Shift_Dirac = B_Dirac - B_Schr."""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = dict(study="NUCLEO_1S phase 6", date="2026-10-13",
         seen=dict(Z1_shift_rel_err={"1s": 0.0145, "2s": 0.0058, "3s": -0.138}, source="phase-1 z1 data"),
         config=dict(Z=[1, 2], states=[1, 2, 3], theta_out=[0.1, 0.05, 0.025], zoom_levels=14, n_per_level=16,
                     rmax_a0=45, extrapolation="quadratic in theta^2 through the 3 points"),
         exams=dict(
             A1=dict(what="Z = 1: alpha^4 shift of 1s and of 2s within 3 % of Dirac's (partly seen: other theta set)", tol=0.03),
             A2=dict(what="Z = 2 (blind): alpha^4 shift of 1s and 2s within 3 % of Dirac's", tol=0.03),
             A3=dict(what="the extrapolation is stable: dropping the coarsest point (linear in theta^2 through 0.05, 0.025) "
                          "changes each 1s/2s shift by <= 2 % (Z = 1 and 2)", tol=0.02),
             report=["3s shifts", "cells", "the part of the shift inside r < 30 lambda_C (fine zone of the old hybrid)"]))
json.dump(P, open(os.path.join(OUT, "nucleo1s_f6_prereg.json"), "w"), indent=1)
print("written")

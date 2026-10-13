"""nucleo1s_f8_prereg.py -> lab/out/nucleo1s_f8_prereg.json : PRE-REGISTRATION of phase 8 (before the production runs).
The project of section 16.2: zoom chain of the zig-zag rule near the nucleus (to R_f = (J0 + 1/2) h), cusp glue
(section 12.4, unchanged: NR Kato amplitude, tau = (K - T_NR) u_c, ghosts = large component G of the inner solution),
outside the TAYLOR-PINNED learned level of AGGLOMERATO_COMPLETO (lab/out/aggfull_relearnT.pkl, read only) reduced to
the s wave (lab/nucleo1s_f8_core.py). alpha^4 shift as a DIFFERENCE of two Z at the same geometry in a0/Z:
(Z = 1, B = 512, cell a0/8) and (Z = 2, B = 256, cell a0/16): delta = (B(2)/4 - B(1))/3, compared with the same formula
on the Dirac levels. Box 90 a0/Z. Control: the same hybrid with the exact band sqrt(1 + k^2) - 1 on the same cells.

DECLARED DEVELOPMENT RUNS (control only, never the learned level): exact band, J0 = 3: delta/Dirac - 1 = +0.96 (1s),
+0.31 (2s), +0.21 (3s); J0 = 5: +0.40 / +0.13 / +0.10. Variants tried and dropped: tau with the continuum Dirac band
(-2.6 at J0 = 3), ghosts in the Foldy-Wouthuysen form G (1 + (E - V)/4m) (+2.4 at J0 = 3, -0.55 with both at J0 = 5).
Reading of the dev numbers (to be tested, H3): the coarse rows (free band + point potential) have no Darwin coupling;
in the large-component form the Darwin energy is spread over the atom (density ~ e^{-2 Z r/a0} for the 1s), so the part
outside R_f is missing and the shift comes out too large; it must go to zero as R_f grows."""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = dict(study="NUCLEO_1S phase 8", date="2026-10-13",
         config=dict(pair=[dict(Z=1, B=512), dict(Z=2, B=256)], rmax_a0_over_Z=90, J0=[3, 4, 5, 6, 8, 12, 16, 24],
                     states=[1, 2, 3], bands=["learned", "exact"], zoom_levels=12, n_per_level=16, ndir=120,
                     glue="cusp (section 12.4), NR tau, ghosts = G", vpot="point"),
         dev=dict(exact_J0_3=[0.96, 0.31, 0.21], exact_J0_5=[0.40, 0.13, 0.10]),
         exams=dict(
             H1=dict(what="the project as asked: learned level, J0 = 3 (R_f = 3.5 cells = 0.44 a0/Z, 30 lambda_C at Z = 2): "
                          "1s and 2s alpha^4 shifts within 5 % of Dirac's (predicted FAIL from the control)", tol=0.05, J0=3),
             H2=dict(what="the learned level vs the exact band, same J0: |delta_L - delta_X| / delta_Dirac <= 5 % for 1s and 2s "
                          "at every J0 of the list (the agglomerated stencil's own error does cancel between the two Z)", tol=0.05),
             H3=dict(what="control (exact band): |error| of 1s and 2s does not grow with J0 (wiggle <= 0.01 allowed) and at J0 = 24 "
                          "(R_f = 3.06 a0/Z) |error 1s| <= 0.10, |error 2s| <= 0.20", tol1=0.10, tol2=0.20, wiggle=0.01),
             H4=dict(what="learned level at J0 = 24: |error 1s| <= 0.10 and |error 2s| <= 0.20", tol1=0.10, tol2=0.20),
             report="3s errors, absolute level errors, inner and outer cell counts, time (not judged)"))
json.dump(P, open(os.path.join(OUT, "nucleo1s_f8_prereg.json"), "w"), indent=1)
print("written")

"""nucleo1s_3d_amend.py -> lab/out/nucleo1s_3d_amend1.json : DECLARED AMENDMENT to nucleo1s_3d_prereg.json (commit 8e20ad2),
written after the pre-registration and before the runs it concerns. No threshold of the pre-registration is changed.

Why: (1) the coordinator asked (13/10) to measure the valley ("taste") split of the 1s separately in the error budget, to
say which part of the -11.6 % of debt #6 is valley split, and to try the nucleus at the CENTRE of a cell (on a site)
instead of on the corner of 8 cells (suggestion of the main thread, docs/BATTIMENTO_1S.md on master 83485bc, which found
the same 4-valley multiplet independently);  (2) the pre-registered tau glue (G1, G2) was found in development, AFTER
the pre-registration commit, to remove the J0-parity jumps but to leave a smooth drift (-0.36 % ... -0.24 % for
J0 = 2..6, Z = 1 1s): its window adds its own aliasing. G1/G2 stay as pre-registered (they will be judged and will
almost surely FAIL). A second glue, the CUSP glue, was then tried in development on the Z = 1 1s ONLY
(J0 = 2..6: +0.0378, +0.0384, +0.0391, +0.0390, +0.0391 %): it is pre-registered here for the states NOT seen
(Z = 1 2s, He+ 1s) and re-run for the Z = 1 1s (declared as seen).
"""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
am = dict(
    amends="lab/out/nucleo1s_3d_prereg.json (commit 8e20ad2)", date="2026-10-13",
    cusp_glue=dict(method="as the tau glue, but the defect comes from the Coulomb cusp only: tau = s * [K u_c - T_c u_c], "
                          "u_c = r e^{-lam r} (lam = 2 m_kin Z alpha), s fixed by the Kato law u''(0+) = -2 m_kin Z alpha u'(0) "
                          "with u'(0) from the inner zig-zag solution (linear fit of u/r on r in [0.5, 4] lambda_C); no window",
                   J0=[3, 4, 5, 6], states=[[1, 1], [1, 2], [2, 1]]),
    centre=dict(what="nucleus on a site (centre of a cell): origin = (a/2, a/2, a/2); the site at r = 0 gets the cube "
                     "average of -Z alpha/r, -(3 ln(2+sqrt3) - pi/2) Z alpha / a = -2.38008 Z alpha/a (hand-made: debt)",
                a=[0.5, 0.3], colours=4, nt_T=500.0),
    exams=dict(
        G1b=dict(what="cusp glue: Z = 1 1s moves by <= 0.02 % over J0 = 3..6 (SEEN in development)", tol=2e-4),
        G2b=dict(what="cusp glue: Z = 1 2s and He+ 1s move by <= 0.02 % over J0 = 3..6 and are within 0.1 % of Dirac; "
                      "Z = 1 1s within 0.1 %", tol_move=2e-4, tol_dirac=1e-3),
        V1=dict(what="valley budget (bare scans a = 0.4, 0.3, 0.24, 0.2): error of the 1s centroid = error of the "
                     "nucleus-blind pair U + valley part (centroid - U). Exam: at every a the valley part is >= 50 % of "
                     "|centroid error| (the multiplet, not the bulk, is the main 1s error)", min_frac=0.5),
        V2=dict(what="the -11.6 %: the ZIGZAG A1 extrapolation (quadratic in a^2 through a = 0.4, 0.3, 0.24) applied to "
                     "U alone is within 1 % of Dirac (the valley part is what spoiled it)", tol=0.01),
        V3=dict(what="valley split below threshold with the learned nucleus: max - min of the 1s lines (weight >= 0.05) "
                     "<= 0.002 mc^2 at every learned a (0.5 ... 0.2)", tol=0.002),
        CN1=dict(what="centre nucleus: the valley split (max - min of 1s lines, weight >= 0.05) at a = 0.3 is smaller than "
                      "with the corner nucleus, in every colour", ),
        CN2=dict(what="centre nucleus: the 4 colours are no longer equivalent (line positions differ by > 1e-4 between "
                      "the colour holding the nucleus and the others)", tol=1e-4),
    ),
)
json.dump(am, open(os.path.join(OUT, "nucleo1s_3d_amend1.json"), "w"), indent=1)
print("written")

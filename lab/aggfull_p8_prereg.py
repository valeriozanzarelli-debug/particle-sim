"""aggfull_p8_prereg.py - PRE-REGISTRATION of phase 8 (docs/AGGLOMERATO_COMPLETO.md s.17), committed BEFORE production.
Development (declared, numbers seen): pair (5,10) at J0 = 32 (R_f = 20 a0/Z): the 1s glue residual has NO root (the outer
region is a far tail, the matching is degenerate) -> J0 is mapped to the SAME glue radius in a0/Z for every pair:
J0(Z) = round(8 R / Z - 1/2) for R = 3.06, 4.06, 6.06, 8.06 a0/Z (= J0 24, 32, 48, 64 of the Z = 1 pair).
Dev with theta/4: pair 3: J10 1s +0.07 %, 2s +0.33 %; J21 +0.06 / +0.10 %; pair 5: J6 +1.79 / -3.50 %; J12 +0.07 / -0.16 %."""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
Rs = [3.0625, 4.0625, 6.0625, 8.0625]
J = {str(Z): [int(round(8 * R / Z - 0.5)) for R in Rs] for Z in (1, 3, 5)}
P = dict(
    study="AGGLOMERATO_COMPLETO phase 8: hybrid atom with the patches (chain cell theta/4 + T level + phase-6 coupling + glue)",
    date="2026-10-13",
    config=dict(variants=[dict(name="p8", theta_factor=0.25, J0=J),
                          dict(name="th1", theta_factor=1.0, J0={k: [v[-1]] for k, v in J.items()})],
                Z=[1, 3, 5], pairs="(Z, B=512) with (2Z, B=256): same geometry in a0/Z (cell Z/8 a0/Z)", states=[1, 2],
                rmax_a0_over_Z=90, zoom_levels=12, n_per_level=16, R_glue_a0_over_Z=Rs, J0_map=J),
    exams=dict(
        Q1=dict(what="pair (1,2): 1s and 2s alpha^4 within +-1 % of Dirac at J0 = 32, 48, 64", tol=0.01, prediction="PASS"),
        Q2=dict(what="pairs (3,6) and (5,10): 1s and 2s within +-1 % at R >= 4.06 a0/Z (J0 = 10/16/21 and 6/9/12)", tol=0.01,
                prediction="FAIL (dev: pair 5 at J0 = 6, only 6 coarse cells: 2s -3.5 %); PASS expected at R >= 6.06"),
        Q3=dict(what="pairs (3,6), (5,10) at R >= 6.06 a0/Z within +-1 %", tol=0.01, prediction="PASS"),
        Q4=dict(what="NEW prediction from phase 7: the unpatched chain (cell theta) at R = 8.06 a0/Z; the Z^3 artefact enters the pair "
                     "(Z, 2Z) as c3/(3 Z D4): 1s error = -11.87 % / Z (phase-7 c3) within 30 % for Z = 1, 3, 5", tol=0.30,
                prediction="PASS"),
        F3S=dict(what="agg2s F3S with aggfull_relearnT.pkl: NOT re-run (agg2s files are not mine); their committed F3S_T result is reported",
                 prediction="already PASS (agg2s commit 5d60a03)")),
)
json.dump(P, open(os.path.join(OUT, "aggfull_p8_prereg.json"), "w"), indent=1)
print(J)

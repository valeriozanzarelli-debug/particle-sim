"""aggfull_p7_prereg.py - PRE-REGISTRATION of phase 7 (docs/AGGLOMERATO_COMPLETO.md s.16), committed BEFORE production.
Development runs done before this file (declared, numbers seen): glue (aggfull_p7_glue.py dev) at J0 = 32, 64 with the
chain's outer cell theta, theta/2, theta/4 (same for both Z) -> 1s -10.5 %, -1.0 %, -0.03 %; 'geometric' scaling
theta/Z (-700 %: the walk's kinetic mass tan(theta)/theta is in lambda_C, it does not cancel); zoom levels 0/1/4/8/12/16
(-122 / -41 / -10 / -10.5 / -10.5 / -10.5 %); n_per_level 4 / 16 / 64 (-0.5 / -10.5 / -45 %); the chain ALONE
(aggfull_p7_chain.py, no coarse level) gives the same -10.48 %; four Z (1-4) at theta for the 1s: relative error after the
known mass is linear in Z (a Z^3 term): Z^3 + Z^4 fit gives Z^4 residual +1.0...+1.5 %."""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = dict(
    study="AGGLOMERATO_COMPLETO phase 7: the -10.5 % 1s plateau of the glue; the c4 of K",
    date="2026-10-13",
    config=dict(
        glue=dict(variants=[["th1", [1, 1], 12, 16], ["th2", [0.5, 0.5], 12, 16], ["th4", [0.25, 0.25], 12, 16]],
                  J0=[24, 32, 48, 64], states=[1, 2], lc=True, note="theta factor of the chain's outer cell, same for both Z"),
        chain=dict(variants=[["c1", 1.0, 12, 16], ["c2", 0.5, 12, 16]], Z=[1, 2, 3, 4], states=[1, 2], rmax_a0_over_Z=30,
                   fit="B/m_walk - B_Dirac = e3 Z^3 + e4 Z^4 (+ e2 Z^2 + e6 Z^6 reported); m_walk = tan(theta)/theta"),
    ),
    exams=dict(
        P1=dict(what="the chain ALONE (no coarse level, no glue) reproduces the 1s plateau of the glue at theta, J0 = 64, within 0.5 points",
                tol=0.005, prediction="PASS (seen in dev)"),
        P2=dict(what="glue with the chain's outer cell theta/4 (same for both Z): 1s and 2s alpha^4 within +-1 % of Dirac at every J0 in 24, 32, 48, 64",
                tol=0.01, prediction="PASS (seen in dev at 32, 64)"),
        P3=dict(what="term by term at theta (chain, four Z, 1s): -e3/3 (the Z^3 term read as alpha^4 by the two-Z formula) + e4 reproduce "
                     "the two-Z plateau within 1 point, and the Z^4 residual |e4 / (alpha^4/8)| <= 2 %",
                tol_point=0.01, tol_res=0.02, prediction="PASS (seen in dev)"),
        P4=dict(what="NEW prediction: the Z^3 term is a contact-like defect (weight psi_n(0)^2 ~ 1/n^3): e3(2s) / e3(1s) = 1/8 within 15 %",
                tol=0.15, prediction="PASS"),
        P5=dict(what="NEW prediction: the Z^3 term is a lattice effect of the chain: |e3(theta/2)| <= |e3(theta)| / 4", prediction="PASS"),
        K1=dict(what="point 2: the Taylor-pinned level T (aggfull_relearnT.pkl) is K re-learned with c4 pinned: its k^4 coefficient within 0.5 % "
                     "of Dirac (-k^4/8, m = 1) in all 13 directions at a0/16 and a0/8; parity-odd part <= 1e-10; (2p: A1 of phase 6 already PASS)",
                tol=0.005, prediction="PASS"),
    ),
    threshold_rule="plateau within +-1 % of Dirac (P2) OR the cause explained term by term (P3 + P4)",
)
json.dump(P, open(os.path.join(OUT, "aggfull_p7_prereg.json"), "w"), indent=1)
print("written")

"""nucleo1s_f5_prereg.py -> lab/out/nucleo1s_f5_prereg.json : PRE-REGISTRATION of phase 5 of NUCLEO_1S, written and
committed BEFORE the comparison runs. Thresholds of phases 1-4 unchanged.

(1) Can s0 (the factor on -Z alpha/r of the 2 sites nearest to the nucleus in each colour) be DERIVED instead of learned?
Predictions are computed HERE (lab/nucleo1s_f5_geo.py) and frozen in the json before any learned value at Z != 60 exists:
  P_A  colour-cell average (Voronoi cell of the colour's bcc lattice, truncated octahedron): s0 = 0.951, s1 = 1.0009,
       independent of Z and a. DECLARED: s0 of P_A was computed while the learned Z = 60 values (0.555 at a = 0.5,
       0.575 at a = 0.3) were already known: P_A for s0 is ALREADY SEEN TO FAIL at Z = 60. P_A for s1 is kept (it gives
       ~1.001, the learned s1 were 1.0035 / 1.0025).
  P_B1 (PRIMARY) exact scattering of the first cell, radial form of the same rule: s such that the radial 1s with uniform
       cells a and the Coulomb phase of the first site (r = a/2) multiplied by s equals the radial 1s of the same rule with
       the inner zone resolved (zoom of phase 1, 12 levels). Root closest to 1. Depends on Z and a.
  P_B2 (reported) same, first radial cell [0, sqrt3 a] (its site at the 3-D nearest-site radius sqrt3/2 a).
  The values of P_B at Z = 60 can be compared with the known learned ones only as "seen"; the test is Z = 20, 40, 80.
  Learned values at Z = 20, 40, 80: same internal criterion as phases 3-4 (one line AND at the bare nucleus-blind level U),
  a = 0.5, listening T = 2000; the bare U comes from a bare scan (T = 500) at the same Z and a. Box: Z = 20: 40 lambda_C
  (the 1s is 6.9 lambda_C wide), else 24. Listening window [1 - 2B, 1 - B/2] with B the Dirac binding (only where to
  listen; [0.80, 0.95] at Z = 60).
(2) the 9x9 junction in 3-D (lab/nucleo1s_junction3d.py), plane interface far from the nucleus, a = 0.25 <-> 0.125,
  junction coin m a/2 (from the 1-D prototype, phase 4).
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_f5_geo import P_A, P_B  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
PA = P_A()
pred = {"PB1": {}, "PB2": {}, "PA": {}}
raw = {}
for Z in (20, 40, 60, 80):
    for a in (0.5, 0.4, 0.3, 0.24, 0.2):
        for v in (1, 2):
            r = P_B(Z, a, v)
            raw[f"PB{v}_Z{Z}_a{a}"] = r
            pred[f"PB{v}"][f"Z{Z}_a{a}"] = [round(r["s"], 6) if r["s"] is not None else None, round(PA["s1"], 6)]
            print(v, Z, a, r["s"], r["roots"], flush=True)
        pred["PA"][f"Z{Z}_a{a}"] = [round(PA["s0"], 6), round(PA["s1"], 6)]
P = dict(
    study="NUCLEO_1S phase 5", date="2026-10-13",
    keeps="thresholds of phases 1-4 unchanged",
    P_A=PA, predictions=pred, P_B_raw=raw,
    config=dict(
        box={"20": 40.0, "40": 24.0, "60": 24.0, "80": 24.0}, bare_T=500.0, use_T=500.0,
        learn=dict(a=0.5, T=2000.0, start=[0.6, 1.0], simplex=[[0.6, 1.0], [0.8, 1.0], [0.6, 1.1]], maxfev=40),
        use=dict(Z60=[0.4, 0.24, 0.2], other=dict(Z=[20, 40, 80], a=0.5)),
        junction=dict(cases=[
            dict(a=0.25, k=[0.2, 0, 0], Ny=1, Nz=1, sig=10.0), dict(a=0.25, k=[0.4, 0, 0], Ny=1, Nz=1, sig=10.0),
            dict(a=0.25, k=[0.8, 0, 0], Ny=1, Nz=1, sig=10.0),
            dict(a=0.25, k=[0.4, 0, 0], Ny=1, Nz=1, sig=10.0, thj=0.0),
            dict(a=0.25, k=[0.4, 0.785398, 0], Ny=32, Nz=1, sig=10.0),
            dict(a=0.25, k=[0.8, 0.785398, 0], Ny=32, Nz=1, sig=10.0),
            dict(a=0.5, k=[0.4, 0, 0], Ny=1, Nz=1, sig=10.0)],
            note="R = probability left in the coarse side after the passage (+ what left through the far edges: "
                 "conservative). thj = 0 case and a = 0.5 case reported, not judged."),
    ),
    exams=dict(
        S1=dict(what="PRIMARY: P_B1 s0 within 0.03 of the learned s0 at Z = 20, 40, 80 (a = 0.5)", tol=0.03),
        S2=dict(what="P_A s1 within 0.01 of the learned s1 at Z = 20, 40, 80 (a = 0.5)", tol=0.01),
        S3=dict(what="use without learning, Z = 60: s = (P_B1 s0(a), P_A s1) at a = 0.4, 0.24, 0.2: 1s within 1 % of Dirac "
                     "and spread (weight >= 0.05, window) <= 0.002", tol=0.01, spread=0.002),
        S4=dict(what="use without learning, Z = 20, 40, 80 at a = 0.5 with P_B1: spread <= 0.002 (valleys closed)",
                spread=0.002),
        J1=dict(what="3-D junction exactly unitary: |norm + edge losses - 1| <= 1e-10 in every case", tol=1e-10),
        J2=dict(what="normal incidence, a = 0.25 <-> 0.125, k = 0.2, 0.4, 0.8: R_conservative <= 1e-3", tol=1e-3),
        J3=dict(what="oblique incidence (k_y = 0.785), k_x = 0.4, 0.8, a = 0.25: R_conservative <= 3e-3", tol=3e-3),
    ),
)
json.dump(P, open(os.path.join(OUT, "nucleo1s_f5_prereg.json"), "w"), indent=1, default=float)
print("written")

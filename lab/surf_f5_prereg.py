"""surf_f5_prereg.py -> lab/out/surf_f5_prereg.json : PHASE 5 of SURFISTI_ESATTI, written BEFORE the a = 0.2 run.
The law d(a) = d0 + c a^(2 gamma) was found AFTER phase 4 (amendment D1); here it is pre-registered and tested on a
lattice never used: a = 0.2 (120^3, measured peak memory 1.3 GB in a test of 6 ticks with 2e5 surfers)."""
import json
import os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R4 = json.load(open(os.path.join(OUT, "surf_f4_results.json")))
d4 = {float(k): v for k, v in R4["exams"]["G5"]["value"]["deficit"].items()}
ZA = 60 / 137.035999084
g2 = 2 * np.sqrt(1 - ZA**2)
x = np.array(sorted(d4)); y = np.array([d4[a] for a in x])
A = np.vstack([np.ones(len(x)), x**g2]).T
c, *_ = np.linalg.lstsq(A, y, rcond=None)
pred = float(c[0] + c[1] * 0.2**g2)
P = dict(
    study="SURFISTI_ESATTI fase 5", date="2026-10-13",
    law="d(a) = d0 + c a^(2 gamma), gamma = sqrt(1 - (Z alpha)^2), 2 gamma = %.4f (contact law of NUCLEO_1S)" % g2,
    fit_phase4=dict(points=d4, d0=float(c[0]), c=float(c[1])),
    prediction=dict(a=0.2, deficit=pred, v=1 - pred),
    setup=dict(a=0.2, n=120, box=24.0, nucleus="NUCLEO_1S sec. 14 transfer line: s(0.2) = (0.584934, 1.002032)",
               preparation="identical to phase 4 (spectrum pass + 2 filtering passes of 1024/a = 5120 ticks, projection on "
                           "valley k = 0)",
               run="shorter than phase 4 (time budget at 120^3, ~1 s per tick): physical time 600 lambda_C/c = 3000 ticks; "
                   "windows of 24 ticks (4.8 lambda_C/c) starting at t = 250, 500 and one closing at 600. v(0.2) = mean of the "
                   "windows at 250 and 500, exactly as in phase 4 (the phase-4 deficits use the same two windows)."),
    exams=dict(
        H1=dict(what="prediction: |d(0.2) measured - %.4f| <= 0.005" % pred, thr=0.005),
        H2=dict(what="continuum limit: 5-point fit d0 + c a^(2 gamma) (a = 0.5, 0.4, 0.3, 0.25, 0.2): |d0| <= 0.01", thr=0.01),
        H3=dict(what="the law fits: R^2 >= 0.999 on the 5 points", thr=0.999),
        H4=dict(what="a = 0.2 run is sound: prepared valley-0 weight >= 0.99, <Sigma_z> >= 0.9; fidelity >= 0.9 to the end; "
                     "surfer v_phi/(Za sin) in [0.75, 1.25] in all 3 windows; max |z| <= 4")),
)
json.dump(P, open(os.path.join(OUT, "surf_f5_prereg.json"), "w"), indent=1)
print("written lab/out/surf_f5_prereg.json; prediction d(0.2) =", pred, "v =", 1 - pred, "fit", c)

"""protetta_vestito.py - exams D1-D3 of lab/out/protetta_prereg.json: is the dressed limiting speed of the electron the
same at p >> m (what a flipless constituent of the composite photon carries) as at p << m?
D_j(pd) = dv_j(p = pd*theta) - dv_j(0.05*theta), exact discrete-time one loop (Debiti engine, read-only copy).
Checkpoints in lab/out/protetta_vestito.json.  Run: timeout 590 python3 lab/protetta_vestito.py  (exit 3 = rerun)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import protetta_vf3 as F

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "protetta_vestito.json")
T0 = time.time()
CUTS = {"K12": (0.6, 1.05), "K12/3": (0.2, 0.35)}
GRID = dict(nr=96, nct=20, nph=32)
THS = [0.0125, 0.00625, 0.003125, 0.0015625, 0.00078125]
TASKS = []
for cut in CUTS:
    for th in THS:
        for pd in (0.05, 1, 3, 10):
            if pd * th <= 0.04:
                TASKS.append(("R", cut, th, pd, (0, 1)))
for cut in CUTS:                                   # D3 massless limit
    for pd in (0.05, 25.6, 51.2, 76.8):
        TASKS.append(("R", cut, 3.90625e-4, pd, (0, 1)))
for th in (0.003125, 0.0015625, 0.00078125):       # D2 cyclic order
    for pd in (0.05, 10):
        TASKS.append(("C", "K12/3", th, pd, (0,)))


def key(t):
    return "%s|%s|%g|%g" % t[:4]


D = json.load(open(FN)) if os.path.exists(FN) else {}
for t in TASKS:
    k = key(t)
    if k in D:
        continue
    if time.time() - T0 > 450:
        sys.exit(3)
    rule, cut, th, pd, axes = t
    t1 = time.time()
    r = F.dv_exact(th, rule, CUTS[cut], pd=pd, axes=axes, **GRID)
    r["time_s"] = time.time() - t1
    D[k] = r
    json.dump(D, open(FN, "w"), indent=1)
    print(k, ["%.6f" % x if x is not None else None for x in r["dv_over_alpha"]], "%.0f s" % r["time_s"], flush=True)
print("ALL DONE")

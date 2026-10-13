"""vuoto_f2_se.py - phase 2 production (pre-registered in lab/out/vuoto_f2_prereg.json): self-energy loop of the rule
(electron walk x light walk in k): mass shift dm (S exams) and g - 2 (G exams). Checkpoints in lab/out/vuoto_f2_se.json.
Run: timeout 590 python lab/vuoto_f2_se.py   (exit 3 = rerun me; one task at most per ~9 min)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from vuoto_f2_core import self_energy, g2_all

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "vuoto_f2_se.json")
T0 = time.time()
TASKS = []
for lat in (False, True):
    for th in (0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125):
        TASKS.append(("dm", lat, th, 0.0, "std", ""))
GROUPS = {"A": ((2, 0), (2, 1), (0, 1)), "B": ((0, 2), (1, 0), (1, 2))}     # the 6 (Q axis, current axis) pairs
for lat in (False, True):
    for th in (0.05, 0.025, 0.0125):
        for lf in (1e-3, 3e-4):
            for gname in ("A", "B"):
                TASKS.append(("g2", lat, th, lf, "std", gname))
for gname in ("A", "B"):
    TASKS.append(("g2", True, 0.025, 1e-3, "fine", gname))   # angular/radial convergence check of the lattice loop
GRID = {"dm": {"std": dict(nr=96, nct=16, nph=32)},
        "g2": {"std": dict(nr=80, nct=16, nph=24), "fine": dict(nr=96, nct=20, nph=30)}}


def key(t):
    return "%s|%s|%g|%g|%s%s" % (t[0], "lat" if t[1] else "cont", t[2], t[3], t[4], ("|" + t[5]) if t[5] else "")


D = json.load(open(FN)) if os.path.exists(FN) else {}
for t in TASKS:
    k = key(t)
    if k in D:
        continue
    if time.time() - T0 > 150 and t[0] == "g2":
        sys.exit(3)
    if time.time() - T0 > 480:
        sys.exit(3)
    kind, lat, th, lf, grid, gname = t
    t1 = time.time()
    if kind == "dm":
        r = self_energy(th, lattice=lat, **GRID["dm"][grid])
        r = {kk: float(v) for kk, v in r.items()}
    else:
        r = g2_all(th, lattice=lat, lam=lf * th, Q=1e-2 * lf * th, pairs=GROUPS[gname], **GRID["g2"][grid])
    r["time_s"] = time.time() - t1
    D[k] = r
    json.dump(D, open(FN, "w"), indent=1)
    print(k, "%.0f s" % r["time_s"], flush=True)
print("ALL DONE")

"""vuoto_f3_g2.py - phase 3 production (pre-registered): g-2 orientation spread with the phase-2 H_eff machinery at the
4th lattice spacing theta = 0.00625 (slice 0), and at theta = 0.0125 on the other time slice (slice 1).
Same settings as lab/vuoto_f2_se.py (grid nr=80, nct=16, nph=24; lam/m = 1e-3, 3e-4; groups A, B).
Checkpoints in lab/out/vuoto_f3_g2.json.  Run: timeout 590 python lab/vuoto_f3_g2.py   (exit 3 = rerun me)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import vuoto_f2_core as G

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "vuoto_f3_g2.json")
T0 = time.time()
GROUPS = {"A": ((2, 0), (2, 1), (0, 1)), "B": ((0, 2), (1, 0), (1, 2))}
TASKS = [(0, 0.00625, lam, g) for lam in (1e-3, 3e-4) for g in "AB"] + [(1, 0.0125, lam, g) for lam in (1e-3, 3e-4) for g in "AB"]
D = json.load(open(FN)) if os.path.exists(FN) else {}
for t in TASKS:
    k = "g2|lat|slice%d|%g|%g|%s" % t
    if k in D:
        continue
    if time.time() - T0 > 300:
        sys.exit(3)
    sl, th, lam, g = t
    G.SLICE["s"] = sl
    t1 = time.time()
    r = G.g2_all(th, True, Q=1e-2 * lam * th, lam=lam * th, nr=80, nct=16, nph=24, pairs=GROUPS[g])
    r["time_s"] = time.time() - t1
    D[k] = r
    json.dump(D, open(FN, "w"), indent=1)
    print(k, "%.0f s" % r["time_s"], flush=True)
print("ALL DONE")

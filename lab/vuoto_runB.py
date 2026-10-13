"""vuoto_runB.py - exploration B (pre-registered in lab/out/vuoto_prereg_B.json). Checkpointed in lab/out/vuoto_B.json.
Run: timeout 590 python lab/vuoto_runB.py  (exit 3 = rerun me)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from vuoto_core import chi_full_3d
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "vuoto_B.json")
T0 = time.time()
T = [(th, mu, dn, qm, 16) for th in (0.05, 0.025) for mu in (0.3, 0.6) for dn in "xyz" for qm in (0.05, 0.025, 0.3, 1, 3)]
T += [(0.025, 0.3, "x", qm, 20) for qm in (0.05, 1)]
D = json.load(open(FN)) if os.path.exists(FN) else {}
for t in T:
    k = "%g|%g|%s|%g|%d" % t
    if k in D:
        continue
    if time.time() - T0 > 540:
        sys.exit(3)
    th, mu, dn, qm, ns = t
    q = np.zeros(3); q["xyz".index(dn)] = qm * th
    c = chi_full_3d(q, th, nseg=ns, mu=mu)
    D[k] = [c.real / (qm * th) ** 2, c.imag / (qm * th) ** 2]
    json.dump(D, open(FN, "w"), indent=0)
print("ALL DONE")

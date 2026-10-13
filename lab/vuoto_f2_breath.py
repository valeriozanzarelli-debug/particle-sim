"""vuoto_f2_breath.py - phase 2 production (pre-registered in lab/out/vuoto_f2_prereg.json, exams B5-B9): the
'breathing coin' (flip angle theta - mu on even ticks, theta + mu on odd ticks) with mu = sqrt(theta) chosen by a
principle (the regulator mass sits at the geometric mean of the electron mass and the lattice cutoff). EXPLORATION:
the frozen rule (lab/zz_core.py) is NOT changed. Checkpoints in lab/out/vuoto_f2_breath.json.
Run: timeout 590 python lab/vuoto_f2_breath.py   (exit 3 = rerun me)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from vuoto_core import chi_3d, chi_full_3d, species_weight
from vuoto_f2_core import band_gaps

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "vuoto_f2_breath.json")
T0 = time.time()
THS = (0.025, 0.0125)
QM = (0.05, 0.025, 1, 3)
TASKS = []
for th in THS:
    mu = float(np.sqrt(th))
    for m_ in (0.0, mu):
        for dn in "xyz":
            for qm in QM:
                TASKS.append(("phys", th, m_, dn, qm))
for th in THS:
    mu = float(np.sqrt(th))
    for dn in "xyz":
        for qm in QM:
            TASKS.append(("full", th, mu, dn, qm))

D = json.load(open(FN)) if os.path.exists(FN) else {}
if "bands" not in D:
    D["bands"] = {"%g" % th: band_gaps(th, float(np.sqrt(th))) for th in THS}
    json.dump(D, open(FN, "w"), indent=0)
w = species_weight()
for t in TASKS:
    k = "%s|%g|%.6g|%s|%g" % t
    if k in D:
        continue
    if time.time() - T0 > 480:
        sys.exit(3)
    kind, th, mu, dn, qm = t
    q = np.zeros(3); q["xyz".index(dn)] = qm * th
    if kind == "phys":
        c = chi_3d(q, th, weight=w, n=64, beta=6.0 if th < 0.02 else 5.0, mu=mu)
    else:
        c = chi_full_3d(q, th, nseg=20, mu=mu)
    D[k] = c.real / (qm * th) ** 2
    json.dump(D, open(FN, "w"), indent=0)
print("ALL DONE")

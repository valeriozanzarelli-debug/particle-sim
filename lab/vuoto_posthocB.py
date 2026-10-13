"""vuoto_posthocB.py - POST-HOC look (declared, NOT graded) after exploration B failed: does the whole sea approach
Uehling when the regulator mass mu is smaller (more 'continuum-like' regulators)? Writes lab/out/vuoto_posthocB.json."""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from vuoto_core import chi_full_3d
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "vuoto_posthocB.json")
T0 = time.time()
D = json.load(open(FN)) if os.path.exists(FN) else {}
for th, mu in ((0.0125, 0.1), (0.0125, 0.2), (0.025, 0.1)):
    for dn in "xyz":
        for qm in (0.05, 0.025, 1, 3):
            k = "%g|%g|%s|%g" % (th, mu, dn, qm)
            if k in D:
                continue
            if time.time() - T0 > 540:
                sys.exit(3)
            q = np.zeros(3); q["xyz".index(dn)] = qm * th
            c = chi_full_3d(q, th, nseg=20, mu=mu)
            D[k] = c.real / (qm * th) ** 2
            json.dump(D, open(FN, "w"), indent=0)
print("ALL DONE")

"""vuoto_f2_hproj.py - phase 2 production (pre-registered in lab/out/vuoto_f2_prereg.json, exams H): hydrogen 1S and 2S
Uehling shifts with the physical-species projector applied ALSO to the electron (the same w(|k a|), K1 = 0.6, K2 = 1.05
used for the vacuum; DST-II along the radius, lab/vuoto_levels.smooth_species). Vacuum table: theta_vac = 0.0125
(phase 1, lab/out/vuoto_3d*.json). Checkpoints in lab/out/vuoto_f2_hproj.json.
Run: timeout 590 python lab/vuoto_f2_hproj.py   (exit 3 = rerun me)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from zz_core import ALPHA
from vuoto_core import MHz
import vuoto_levels as VL

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "vuoto_f2_hproj.json")
T0 = time.time()
TASKS = [(name, th) for th in (0.05, 0.025, 0.0125) for name in ("1s", "2s")]
SPEC = {"1s": (-1, 1, 40 / ALPHA), "2s": (-1, 2, 80 / ALPHA)}

D = json.load(open(FN)) if os.path.exists(FN) else {}
qs, tab, p0 = VL.load_pitR(0.0125)
pit = VL.make_interp(qs, tab)
qmax = qs[-1]
if "ref_MHz" not in D:
    VL.ZREF[0] = 1
    lam = ALPHA
    D["ref_MHz"] = {"1s": MHz(VL.ref_shift(VL.P_1s(lam), 60 / lam)), "2s": MHz(VL.ref_shift(VL.P_2s(lam), 120 / lam))}
    json.dump(D, open(FN, "w"), indent=1)
for name, th in TASKS:
    k = "%s|%g" % (name, th)
    if k in D:
        continue
    if time.time() - T0 > 300:
        sys.exit(3)
    kappa, n, rmax = SPEC[name]
    t1 = time.time()
    om, r, wgt = VL.state(kappa, n, th, 1, rmax, smooth=True)
    v = VL.shift_q(pit, r, wgt, 1, qmax)
    D[k] = {"MHz": MHz(v), "omega": float(np.real(om)), "time_s": time.time() - t1}
    json.dump(D, open(FN, "w"), indent=1)
    print(k, D[k], flush=True)
print("ALL DONE")

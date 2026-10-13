"""vuoto_f2_bethe.py - phase 2 production (pre-registered in lab/out/vuoto_f2_prereg.json): Bethe logarithm of hydrogen
1S and 2S from the resolvent of the RADIAL rule (exams L, lab/out/vuoto_f2_prereg_L.json). Source = dipole r g of the
state, low-passed with the species projector, projected on positive energies and on the physical band (Omega below
Omega_n + 0.6). Checkpoints in lab/out/vuoto_f2_bethe.json.
Run: timeout 590 python lab/vuoto_f2_bethe.py   (exit 3 = rerun me)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from vuoto_f2_core import bethe_log

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
AMEND = os.environ.get("VUOTO_AMEND_L") == "1"      # amendment A-L1 (lab/out/vuoto_f2_amend_L.json)
FN = os.path.join(OUT, "vuoto_f2_bethe_A.json" if AMEND else "vuoto_f2_bethe.json")
T0 = time.time()
RMAX = {1: 40.0, 2: 30.0}                 # box in units of n a_B (n a_B x RMAX)
TASKS = [(n, kp, th) for th in (0.1, 0.05) for n in (1, 2) for kp in (1, -2)]
RY = 0.5 * 7.2973525693e-3 ** 2
KRY = [30, 100, 300, 1000, 3000]
KS = RY * np.array(KRY, float)


def key(t):
    return "%d|%d|%g" % t


D = json.load(open(FN)) if os.path.exists(FN) else {}
for t in TASKS:
    k = key(t)
    if k in D:
        continue
    if time.time() - T0 > 120:
        sys.exit(3)
    n, kp, th = t
    t1 = time.time()
    r = bethe_log(n, kp, th, Z=1, rmax_bohr=RMAX[n], Ks=KS, nk=200, h=0.35 if AMEND else 0.5, species=True, omega_cut=0.6,
                  tail=28.0 if AMEND else 19.0)
    r["time_s"] = time.time() - t1
    D[k] = r
    json.dump(D, open(FN, "w"), indent=1)
    print(k, "%.0f s" % r["time_s"], flush=True)
print("ALL DONE")

"""vuoto_f3_dv.py - phase 3 production (pre-registered in lab/out/vuoto_f3_prereg.json): exact dv_j of the dressed
electron for the frozen rule R (cutoffs K12, K12/3) and the proposal C (K12/3), 5 lattice spacings; continuum control.
Checkpoints in lab/out/vuoto_f3_dv.json.  Run: timeout 590 python lab/vuoto_f3_dv.py   (exit 3 = rerun me)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import vuoto_f3_core as F
import vuoto_f2_core as G
from vuoto_core import species_weight, PI
from zz_core import ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "vuoto_f3_dv.json")
T0 = time.time()
THS = [0.05, 0.025, 0.0125, 0.00625, 0.003125]
CUTS = {"K12": (0.6, 1.05), "K12/3": (0.2, 0.35)}
TASKS = ([("R", "K12", t) for t in THS] + [("R", "K12/3", t) for t in THS] + [("C", "K12/3", t) for t in THS]
         + [("cont", c, t) for c in ("K12", "K12/3") for t in THS] + [("check", "K12/3", 0.05)])
GRID = dict(nr=96, nct=20, nph=32)


def cont_dv(th, K12, pd=0.05):
    """continuum Dirac control (phase-2 loop, exact in the continuum): dv along z (isotropic)"""
    w = species_weight(*K12)
    K, WW = G.sphere_rule(GRID["nr"], GRID["nct"], GRID["nph"], 1e-4 * th, K12[1])
    WW = WW * w(K) / (2 * PI) ** 3

    def dsym(P):
        E, V = np.linalg.eigh(G.dirac_H(P[None], th)[0])
        S = [np.real(V[:, i].conj() @ G.sigma_matrix(E[i], P, th, False, 0.0, K, WW) @ V[:, i]) for i in range(4)]
        return 0.5 * (np.mean(S[2:]) - np.mean(S[:2])), E[2]
    d0, m = dsym(np.zeros(3))
    p = np.array([0, 0, pd * th])
    d1, E0 = dsym(p)
    return {"dv_over_alpha": float((E0 * d1 - m * d0) / (E0**2 - m**2) / ALPHA), "dm_over_m": float(d0 / th)}


def key(t):
    return "%s|%s|%g" % t


D = json.load(open(FN)) if os.path.exists(FN) else {}
for t in TASKS:
    k = key(t)
    if k in D:
        continue
    if time.time() - T0 > 400:
        sys.exit(3)
    rule, cut, th = t
    t1 = time.time()
    if rule == "cont":
        r = cont_dv(th, CUTS[cut])
    elif rule == "check":
        r = {}
        for rr in ("R", "C"):
            x = F.self_energy_exact(np.zeros(3), th, rr, CUTS[cut], nr=64, nct=12, nph=16, hole=False)
            r[rr] = {"imag_over_dT": float(np.max(np.abs(x["imag"])) / np.max(np.abs(x["dT"])))}
    else:
        r = F.dv_exact(th, rule, CUTS[cut], **GRID)
    r["time_s"] = time.time() - t1
    D[k] = r
    json.dump(D, open(FN, "w"), indent=1)
    print(k, "%.0f s" % r["time_s"], flush=True)
print("ALL DONE")

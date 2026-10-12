"""reticolo_r2a.py - R2a: speed of light on every level (sparse single level with dynamic regrid).
Packet lambda = 24 cells, sigma = 3 lambda, directions 0/30/45 deg, levels k = 0,1,5,10,20, 480 own ticks.
Also: the rule's prediction for the same packet (lattice group velocity averaged over the spectrum) and Maxwell's.
Writes lab/out/reticolo_r2a_<angles>.json. Usage: python reticolo_r2a.py [ang ...]."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reticolo_core import packet, single_level_world, kgrid, band_vectors, OUT

LAM, SIG, NT = 24, 48.0, 480
ANGS = [int(a) for a in sys.argv[1:]] or [0, 30, 45]
res = {"lambda_cells": LAM, "sigma_cells": SIG, "ticks": NT, "runs": []}
for ang in ANGS:
    n = 640
    F = packet((n, n), 2, LAM, ang, [n / 2, n / 2], SIG)
    # predictions from the spectrum
    Fk = np.fft.fftn(F, axes=(0, 1))
    wk = (np.abs(Fk) ** 2).sum(-1).reshape(-1)
    KA = kgrid((n, n), 2).reshape(-1, 3)
    sel = wk > wk.max() * 1e-8
    dirv = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
    vs, om = band_vectors(KA[sel])
    dk = 1e-5
    grad = []
    for j in range(2):
        Kp = KA[sel].copy(); Kp[:, j] += dk
        Km = KA[sel].copy(); Km[:, j] -= dk
        grad.append((band_vectors(Kp)[1] - band_vectors(Km)[1]) / (2 * dk))
    v_rule = float((wk[sel] * (grad[0] * dirv[0] + grad[1] * dirv[1])).sum() / wk[sel].sum())
    kn = np.linalg.norm(KA[sel, :2], axis=1)
    v_max = float((wk[sel] * (KA[sel, :2] @ dirv) / kn).sum() / wk[sel].sum())
    ref_arr = None
    for k in (0, 1, 10, 20):
        t0 = time.time()
        w = single_level_world(k, 2, 16, np.array([-n // 2, -n // 2]), F, tol=1e-6, every=8)
        L = w.levels[0]
        E0, M0 = w.moments()
        cells = [w.cells()]
        for t in range(NT):
            w.step_top(1)
            if t % 40 == 0:
                cells.append(w.cells())
        E1, M1 = w.moments()
        # physical units: positions in h0, time = NT ticks * dt_k = NT * 2^k h0/c
        v = float((M1 - M0) @ dirv / (NT * L.h))
        arr = L.get_dense(np.array([-n // 2, -n // 2]), (n + 2 * NT, n + 2 * NT)) if False else L.F[:L.nb].copy()
        same = None
        if ref_arr is None:
            ref_arr = (L.keys.copy(), arr)
        else:
            same = bool(np.array_equal(ref_arr[0], L.keys) and np.array_equal(ref_arr[1], arr))
        r = dict(ang=ang, k=k, h=L.h, v=v, err=v - 1, v_rule_pred=v_rule, v_maxwell_same_packet=v_max,
                 energy_drift=E1 / E0 - 1, cells_max=max(cells), cells_end=w.cells(),
                 bitwise_same_as_level0=same, wall=time.time() - t0)
        res["runs"].append(r)
        print(r, flush=True)
json.dump(res, open(os.path.join(OUT, "reticolo_r2a_%s.json" % "_".join(map(str, ANGS))), "w"), indent=1)

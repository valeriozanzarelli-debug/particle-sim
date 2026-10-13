"""nucleo1s_f5_pc.py - computes the P_C predictions of lab/out/nucleo1s_f5_amend2.json -> lab/out/nucleo1s_f5_pc.json"""
import os, sys, json
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA  # noqa: E402
from nucleo1s_f5_geo import _binding, target_binding  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def PC(Z, a, offset):
    Za = Z * ALPHA
    rmax = max(30.0, 10 * 137.036 / Z)
    tgt = target_binding(a, Za, rmax)
    if abs(offset - 0.5) < 1e-12:
        edges = np.arange(int(np.ceil(rmax / a)) + 1) * a
    else:   # sites at (i + offset) a: first cell [0, 2 offset a], then cells a
        R0 = 2 * offset * a
        edges = np.concatenate([[0.0], R0 + np.arange(int(np.ceil((rmax - R0) / a)) + 1) * a])
    f = lambda s: _binding(edges, Za, s, tgt) - tgt
    grid = np.linspace(0.02, 1.0, 50)
    v = [f(s) for s in grid]
    roots = [brentq(f, grid[i], grid[i + 1], xtol=1e-9) for i in range(len(grid) - 1) if v[i] * v[i + 1] < 0]
    return dict(s=max(roots) if roots else None, roots=roots, bare=f(1.0) + tgt, target=tgt)


out = {}
for name, off in (("PC1", 0.5), ("PC2", np.sqrt(3) / 2)):
    for Z in (20, 40, 60, 80):
        for a in (0.02, 0.01):
            r = PC(Z, a, off); out[f"{name}_Z{Z}_a{a}"] = r
            print(name, Z, a, r["s"], r["roots"], flush=True)
            json.dump(out, open(os.path.join(OUT, "nucleo1s_f5_pc.json"), "w"), indent=1)

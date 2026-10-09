"""zz_atom3d_radial.py - follow-up of zz_atom3d.py (written AFTER its run, information only).

The 3-D walk at Z = 60 missed the Dirac 1s binding by +29 / +19 / +10 % at a = 0.4 / 0.3 / 0.24 lambda_C and the
a^2 extrapolation failed (A1). Is that the 3-D rule, or the coarse lattice? Same rule in radial form, same Z,
same lattices, then finer ones. For Z alpha = 0.44 the Dirac 1s behaves like r^(gamma-1), gamma = 0.899, near
the nucleus: a lattice error ~ a^(2 gamma) = a^1.8 is expected there, not a^2.
Run: timeout 590 python lab/zz_atom3d_radial.py
"""
import json, os
import numpy as np
from zz_core import RadialWalk, dirac_formula, ALPHA

OUT = os.path.join(os.path.dirname(__file__), "out")
ZA = 60 * ALPHA
ref = 1 - dirac_formula(1, -1, ZA)
d3 = {}
try:
    for r in json.load(open(os.path.join(OUT, "zz_atom3d.json")))["runs"]:
        d3[r["a"]] = r["binding"]
except Exception:
    pass
res = {}
prev = None
for th in (0.4, 0.3, 0.24, 0.1, 0.05, 0.025):
    rw = RadialWalk(-1, th, 40.0, Zalpha=ZA)
    om, u, _ = rw.eigen(1 - 1.05 * ref)
    b = 1 - om
    res[th] = b
    s = f"theta = {th}: radial walk binding {b:.6f} mc^2 (Dirac {ref:.6f}, {100 * (b / ref - 1):+.2f} %)"
    if th in d3:
        s += f"; 3-D walk {d3[th]:.6f} ({100 * (d3[th] / b - 1):+.2f} % from the radial walk)"
    if prev is not None:
        s += f"; error ratio to previous {(prev[1] / ref - 1) / (b / ref - 1):.2f} for a ratio {prev[0] / th:.2f}"
    prev = (th, b)
    print(s, flush=True)
json.dump({str(k): v for k, v in res.items()}, open(os.path.join(OUT, "zz_atom3d_radial.json"), "w"), indent=1)

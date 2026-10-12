"""fp_g0.py - gate G0: the k-space rule equals zz_core.Walk('light').tick on a periodic box"""
import json, numpy as np
from zz_core import Walk
from fp_core import LightK, to_k, to_x
rng = np.random.default_rng(20261012)
n = 32
F = rng.normal(size=(3, n, n, n)) + 1j * rng.normal(size=(3, n, n, n))
w = Walk('light', n, 1.0)
lk = LightK((n, n, n))
a, Fk = F.copy(), to_k(F)
for t in range(4):
    a = w.tick(a)
    Fk = lk.tick(Fk)
d = np.abs(to_x(Fk) - a).max() / np.abs(a).max()
ok = bool(d <= 1e-10)
print("G0 max rel diff", d, "PASS" if ok else "FAIL")
json.dump({"G0_maxdiff": float(d), "pass": ok}, open("out/fp_g0.json", "w"), indent=1)

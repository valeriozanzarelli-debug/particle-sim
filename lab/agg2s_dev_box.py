"""agg2s_dev_box.py - DEVELOPMENT exploration (declared, not a production run): the box error of the 2s with the EXACT
Schroedinger kinetic k^2/2 and the alias-free Galerkin field (aggfull_core), at cells a0/4 (and a0/2), boxes L = 16..48 a0,
twists: n anti-periodic axes, n = 0..3 (Bloch band of the crystal of images). Writes lab/out/agg2s_dev_box.json.
Usage: timeout 590 python lab/agg2s_dev_box.py CELLDIV N [N ...]   (CELLDIV = 4 -> h = a0/4)"""
import os, sys, json, time
import numpy as np
import agg2s_core as c

JS = os.path.join(c.OUT, "agg2s_dev_box.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
div = int(sys.argv[1])
h = c.A0 / div
t00 = time.time()
for N in [int(a) for a in sys.argv[2:]]:
    f = c.galerkin_field(N, h, cache=False)
    for tw in [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)]:
        key = f"d{div}_N{N}_tw{''.join(map(str, tw))}"
        if key in res and "2s" in res[key]:
            continue
        T = c.kin_schroedinger(N, h, tw)
        at = c.TwistAtom(T, f, tw, sym="even")
        out = {}
        v1, E1, r1, it1, ok1 = c.solve_state(at, c.guess(at, "1s"), tol_res=1e-9, budget=200)
        v2, E2, r2, it2, ok2 = c.solve_state(at, c.guess(at, "2s"), Y=v1[:, None], tol_res=1e-9, budget=300)
        out = {"L_a0": N / div, "1s": E1 * c.MC2, "2s": E2 * c.MC2, "res": [r1, r2], "ok": [ok1, ok2]}
        res[key] = out
        json.dump(res, open(JS, "w"), indent=1)
        print(key, f"L={N/div:.0f} a0  1s {E1*c.MC2/(-c.RY)-1:+.4e}  2s {E2*c.MC2/(-c.RY/4)-1:+.4e}  res {r1:.1e} {r2:.1e} ({time.time()-t00:.0f}s)", flush=True)

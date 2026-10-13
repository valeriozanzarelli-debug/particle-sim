"""agg2s_dev_box2p.py - DEVELOPMENT exploration (declared): box error of the 2p_z (odd in z) with the exact Schroedinger
kinetic and the Galerkin field, cells a0/4, L = 32 and 48 a0, twists 000 / 111. Writes lab/out/agg2s_dev_box2p.json"""
import os, sys, json
import numpy as np
import agg2s_core as c
JS = os.path.join(c.OUT, "agg2s_dev_box2p.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
h = c.A0 / 4
for N in (128, 192):
    f = c.galerkin_field(N, h, cache=False)
    for tw in ((0, 0, 0), (1, 1, 1)):
        key = f"d4_N{N}_tw{''.join(map(str, tw))}_2p"
        if key in res:
            continue
        at = c.TwistAtom(c.kin_schroedinger(N, h, tw), f, tw, sym="oddz")
        v, E, r, it, ok = c.solve_state(at, c.guess(at, "2p"), tol_res=1e-9, budget=400)
        res[key] = dict(L_a0=N / 4, E_rel=E * c.MC2 / (-c.RY / 4) - 1, residual=r)
        json.dump(res, open(JS, "w"), indent=1)
        print(key, res[key], flush=True)

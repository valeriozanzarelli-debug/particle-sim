"""nucleo1s_f6_run.py - phase 6 runs (lab/out/nucleo1s_f6_prereg.json): s levels on the zoom chain, Z = 1 and 2.
Each level is saved as soon as it is computed (re-run until 'all done'). -> lab/out/nucleo1s_f6_results.json"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA  # noqa: E402
from nucleo1s_core import mesh_zoom, level, dirac_binding  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "nucleo1s_f6_prereg.json")))["config"]
RES = os.path.join(OUT, "nucleo1s_f6_results.json")
T0 = time.time()
res = json.load(open(RES)) if os.path.exists(RES) else {}
for Z in P["Z"]:
    Za = Z * ALPHA
    for th in P["theta_out"]:
        key = f"Z{Z}_th{th}"
        if all(f"{key}_n{n}" in res for n in P["states"]):
            continue
        m = mesh_zoom(th, P["rmax_a0"] / Za, P["zoom_levels"], P["n_per_level"])
        for n in P["states"]:
            k = f"{key}_n{n}"
            if k in res:
                continue
            if time.time() - T0 > 450:
                print("time budget: run again"); sys.exit(0)
            b, info = level(m, n, -1, Za, guess_rel=1 + th * th / 3)
            res[k] = dict(B=b, cells=m.N, resid=float(info["resid"]), rel=b / dirac_binding(n, -1, Za) - 1)
            json.dump(res, open(RES, "w"), indent=1)
            print(k, res[k], f"{time.time() - T0:.0f}s", flush=True)
print("all done")

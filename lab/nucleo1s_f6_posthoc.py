"""nucleo1s_f6_posthoc.py - DECLARED post hoc (after A1-A3 FAIL): finer outer cells theta = 0.0125, 0.00625 to see
whether the alpha^4 shifts converge when the coarsest point (0.1) is dropped. Not judged. -> out/nucleo1s_f6_posthoc.json"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA  # noqa: E402
from nucleo1s_core import mesh_zoom, level, dirac_binding  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R = json.load(open(os.path.join(OUT, "nucleo1s_f6_results.json")))
PH = os.path.join(OUT, "nucleo1s_f6_posthoc.json")
res = json.load(open(PH)) if os.path.exists(PH) else {}
T0 = time.time()
for Z in (1, 2):
    Za = Z * ALPHA
    for th in (0.0125, 0.00625):
        m = None
        for n in (1, 2, 3):
            k = f"Z{Z}_th{th}_n{n}"
            if k in res:
                continue
            if time.time() - T0 > 450:
                print("time budget"); sys.exit(0)
            m = m or mesh_zoom(th, 45 / Za, 14, 16)
            b, info = level(m, n, -1, Za, guess_rel=1 + th * th / 3)
            res[k] = dict(B=b, cells=m.N); json.dump(res, open(PH, "w"), indent=1); print(k, b, m.N, flush=True)
allr = {**R, **res}
out = {}
for Z in (1, 2):
    Za = Z * ALPHA
    for n in (1, 2, 3):
        S = Za**2 / (2 * n * n); D = dirac_binding(n, -1, Za)
        for name, ths in (("0.05-0.0125", (0.05, 0.025, 0.0125)), ("0.025-0.00625", (0.025, 0.0125, 0.00625))):
            t = np.array(ths); y = np.array([allr[f"Z{Z}_th{x}_n{n}"]["B"] for x in t])
            ext = np.polyfit(t**2, y, 2)[-1]
            out[f"Z{Z}_{n}s_{name}"] = (ext - S) / (D - S) - 1
res["shift_rel_err"] = out
json.dump(res, open(PH, "w"), indent=1)
for k, v in out.items():
    print(k, f"{100 * v:+.2f} %")

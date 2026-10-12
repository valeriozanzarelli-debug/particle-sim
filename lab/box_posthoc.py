"""box_posthoc.py - POST-HOC (not graded) check of A1a: the same whole-sea count with other kx windows
(K = 0.6, 0.8, 1.3 instead of the graded 1.0), to show how much the count depends on the window edges. Writes lab/out/box_posthoc.json."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from box_core import static_part, masses, floquet, sx_phase, schwinger_strip, PI

TH = 0.1
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


res = {}
for d in (4, 6, 8, 12, 16, 24, 32, 64):
    M = static_part(TH, d); m = masses(TH, d, M)
    for x in (1, 2, 4):
        E = x * m[0] ** 2
        r = {}
        for K in (0.6, 0.8, 1.3):
            r[f"whole_sea_ratio_K{K}"] = schwinger_strip(TH, d, E, K=K, M=M)[0] / float(np.exp(-PI * m**2 / E).sum())
        res[f"d{d}_x{x}"] = r
        print(d, x, {k: round(v, 3) for k, v in r.items()}, flush=True)
json.dump(res, open(os.path.join(OUT, "box_posthoc.json"), "w"), indent=1)

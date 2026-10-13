"""assist_v3diag.py - POST-HOC diagnostic (declared, not graded) for the failed exam V3a.
Is the small even-n splitting a true crossing (numerical floor of the minimiser) or a real gap?
Scan delta(omega) on finer and finer grids around the even-n minimum: a true crossing gives min ~ grid step.
Writes lab/out/assist_v3diag.json."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assist_core import floquet_split_walk, OUT, PI
v3 = json.load(open(os.path.join(OUT, "assist_val.json")))["V3"]
out = {}
for key in ("xi0.1_n2", "xi0.1_n4", "xi0.2_n2", "xi0.2_n4", "xi0.1_n5", "xi0.1_n3"):
    r = v3[key]; xi = r["xi"]; N = r["N_a0.025"]; w0 = r["omega_walk_a0.025"]
    rows = []
    c = w0
    for h in (1e-6, 1e-8, 1e-10, 1e-12):
        ws = c + h * np.arange(-200, 201)
        d = np.array([floquet_split_walk(w, xi, N) for w in ws])
        i = int(np.argmin(d)); c = ws[i]
        rows.append(dict(step=h, min_Delta=float(d[i] * c / (2 * PI)), omega=float(c)))
    out[key] = rows
    print(key, rows, flush=True)
json.dump(out, open(os.path.join(OUT, "assist_v3diag.json"), "w"), indent=1)

"""vuoto_f3_bands.py - phase 3, exam A7 (pre-registered): bare bands of the frozen rule R and of the proposal C.
Writes lab/out/vuoto_f3_bands.json.  Run: python lab/vuoto_f3_bands.py (~1 min)"""
import os, json
import numpy as np
import vuoto_f3_core as F

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
DIRS = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1], "111": [1, 1, 1], "110": [1, 1, 0]}
R = {}
for th in (0.05, 0.0125):
    for rule in ("R", "C"):
        nt = F.n_ticks(rule)
        d = {}
        E0 = np.sort(F.e_bands(np.zeros((1, 3)), th, rule)[0][0]) / nt
        d["rest_mass_err"] = float(np.max(np.abs(np.abs(E0) - th)))
        for name, v in DIRS.items():
            u = np.array(v, float) / np.linalg.norm(v)
            k = 0.01 * th
            E = np.sort(F.e_bands((k * u)[None], th, rule)[0][0])[-1] / nt
            d["v2_" + name] = float((E**2 - th**2) / k**2)
            E2 = np.sort(F.e_bands((0.2 * u)[None], th, rule)[0][0])[-1] / nt
            d["quartic_" + name] = float(E2**2 - th**2 - 0.04)
        for name in ("x", "111"):
            u = np.array(DIRS[name], float) / np.linalg.norm(DIRS[name])
            ks = np.linspace(0, 1.7, 1701)
            top = np.max(np.abs(F.e_bands(ks[:, None] * u[None], th, rule)[0]), 1)
            hit = np.nonzero(top > np.pi - 0.05)[0]
            d["kfold_" + name] = float(ks[hit[0]]) if len(hit) else None
            d["min_dist_to_pi_" + name] = float(np.pi - top.max())
        R["%s|%g" % (rule, th)] = d
X = {}
for th in (0.05, 0.0125):
    r, c = R["R|%g" % th], R["C|%g" % th]
    X["%g" % th] = {"i": c["rest_mass_err"] <= 1e-12,
                    "ii": all(abs(c["v2_" + n] - r["v2_" + n]) <= 1e-8 for n in DIRS),
                    "iii": all(abs(c["quartic_" + n]) <= 1.1 * abs(r["quartic_" + n]) for n in ("111", "110")),
                    "iv": all(c["kfold_" + n] is not None and r["kfold_" + n] is not None and c["kfold_" + n] >= 0.8 * r["kfold_" + n]
                              for n in ("x", "111"))}
R["A7"] = X
R["A7_PASS"] = bool(all(all(v.values()) for v in X.values()))
json.dump(R, open(os.path.join(OUT, "vuoto_f3_bands.json"), "w"), indent=1)
print(json.dumps(R, indent=1))

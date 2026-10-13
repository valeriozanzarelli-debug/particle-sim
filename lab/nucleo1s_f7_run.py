"""nucleo1s_f7_run.py - phase 7 runs (lab/out/nucleo1s_f7_prereg.json) -> lab/out/nucleo1s_f7_results.json (re-run until
'all done'), then the verdicts -> lab/out/nucleo1s_f7_judged.json"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, MC2_EV  # noqa: E402
from nucleo1s_core import mesh_zoom, level, dirac_binding  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "nucleo1s_f7_prereg.json")))
C = P["config"]
RES = os.path.join(OUT, "nucleo1s_f7_results.json")
res = json.load(open(RES)) if os.path.exists(RES) else {}
T0 = time.time()
runs = [(f"Z{Z}", Z, dict()) for Z in C["Z_exam"]] + [(f"diag_{k}", 1, v) for k, v in C["diag"].items()]
for tag, Z, ov in runs:
    cfg = {**{k: C[k] for k in ("theta_out", "zoom_levels", "n_per_level", "rmax_a0")}, **ov}
    Za = Z * ALPHA
    for th in cfg["theta_out"]:
        m = None
        for n in (1, 2, 3):
            k = f"{tag}_th{th}_n{n}"
            if k in res:
                continue
            if time.time() - T0 > 450:
                print("time budget: run again"); sys.exit(0)
            m = m or mesh_zoom(th, cfg["rmax_a0"] / Za, cfg["zoom_levels"], cfg["n_per_level"])
            b, info = level(m, n, -1, Za, guess_rel=1 + th * th / 3)
            res[k] = dict(B=b, cells=m.N, resid=float(info["resid"]))
            json.dump(res, open(RES, "w"), indent=1)
            print(k, b, m.N, f"{time.time() - T0:.0f}s", flush=True)
J = {}
for tag, Z, ov in runs:
    th = np.array(ov.get("theta_out", C["theta_out"]))
    Za = Z * ALPHA
    for n in (1, 2, 3):
        y = np.array([res[f"{tag}_th{t}_n{n}"]["B"] for t in th])
        ext = np.polyfit(th**2, y, 2)[-1]
        S = Za**2 / (2 * n * n); D = dirac_binding(n, -1, Za)
        J[f"{tag}_{n}s"] = dict(rel=(ext - S) / (D - S) - 1, dirac_shift_eV=(D - S) * MC2_EV,
                               rule_shift_eV=(ext - S) * MC2_EV, level_rel_err=ext / D - 1,
                               per_theta_rel=[(b - S) / (D - S) - 1 for b in y])
E = P["exams"]
J["B1"] = bool(all(abs(J[f"Z{Z}_{n}s"]["rel"]) <= E["B1"]["tol"] for Z in C["Z_exam"] for n in (1, 2)))
J["B2"] = bool(all(abs(J[f"Z{Z}_{n}s"]["rel"]) <= E["B2"]["tol"] for Z in C["Z_exam"] for n in (1, 2)))
J["B3"] = bool(any(abs(J[f"diag_{k}_3s"]["rel"]) <= E["B3"]["tol"] for k in C["diag"]))
json.dump(J, open(os.path.join(OUT, "nucleo1s_f7_judged.json"), "w"), indent=1)
for k, v in J.items():
    print(k, v if isinstance(v, bool) else f"{100 * v['rel']:+.3f} %  (level {v['level_rel_err']:+.2e})")
print("all done")

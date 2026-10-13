"""supercrit_run_v1.py - exam V1: point nucleus, Z <= 137, 1s of the rule vs the Dirac formula, a -> 0.
checkpointed in lab/out/supercrit_v1.json (rerun to continue)."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supercrit_core import SWalk, V_point, level, level_1s, ALPHA, dirac_formula
from supercrit_prereg import P

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "supercrit_v1.json")
cfg = P["configs"]["V1"]
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
t0 = time.time()
for Z in cfg["Z"]:
    Za = Z * ALPHA
    ED = dirac_formula(1, -1, Za)
    d = res.setdefault(str(Z), {"E_dirac": ED, "runs": {}})
    rmax = max(30.0, 30.0 / Za)
    guess = ED
    for a in cfg["a"]:
        if f"{a:.6g}" in d["runs"]:
            guess = d["runs"][f"{a:.6g}"]["omega"]
            continue
        if time.time() - t0 > 500:
            break
        w = SWalk(-1, a, rmax, lambda r: V_point(r, Za))
        om, nodes, resid, smooth = level_1s(w, max(min(guess, ED) - 0.1, -0.9), k=30)
        guess = om
        d["runs"][f"{a:.6g}"] = dict(a=a, omega=om, nodes=nodes, resid=resid, smooth=smooth, N=w.N,
                                     rel_err=(om - ED) / (1 - ED))
        print(Z, a, om, ED, (om - ED) / (1 - ED), resid, f"{time.time() - t0:.0f}s", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)

# extrapolation E(a) = E0 + K a^(2g) + B a^2 on the 3 finest a
for Z in cfg["Z"]:
    d = res[str(Z)]
    runs = sorted(d["runs"].values(), key=lambda x: x["a"])
    if len(runs) < len(cfg["a"]):
        continue
    g = np.sqrt(1 - (Z * ALPHA)**2)
    A = np.array([[1, x["a"]**(2 * g), x["a"]**2] for x in runs[:3]])
    E0, K, B = np.linalg.solve(A, [x["omega"] for x in runs[:3]])
    d["extrap"] = dict(E0=E0, K=K, B=B, rel_err=(E0 - d["E_dirac"]) / (1 - d["E_dirac"]),
                       finest_rel_err=runs[0]["rel_err"], two_gamma=2 * g)
    print("Z", Z, "extrap", E0, "Dirac", d["E_dirac"], "rel err", d["extrap"]["rel_err"], "finest", runs[0]["rel_err"])
json.dump(res, open(OUT, "w"), indent=1)

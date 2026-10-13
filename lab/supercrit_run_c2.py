"""supercrit_run_c2.py - exams C2 (shell nucleus R = 0.031) and C3 (uniform sphere R = 1.2 (2.5 Z)^(1/3) fm):
critical charge of the rule for each a, extrapolation a -> 0, and the continuum Dirac ODE reference.
Checkpointed in lab/out/supercrit_c2.json. Usage: python supercrit_run_c2.py [shell|sphere|cont]"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supercrit_core import SWalk, V_shell, V_sphere, R_sphere_fm, FM, level, ALPHA, dirac_critical_Z, dirac_bound
from supercrit_run_c1 import track_to_critical
from supercrit_prereg import P

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "supercrit_c2.json")
res = json.load(open(OUT)) if os.path.exists(OUT) else {}


def Vmodel(model, Z):
    Za = Z * ALPHA
    if model == "shell":
        return lambda r: V_shell(r, Za, 0.031)
    R = R_sphere_fm(Z) * FM
    return lambda r: V_sphere(r, Za, R)


def save():
    json.dump(res, open(OUT, "w"), indent=1)


t0 = time.time()
what = sys.argv[1] if len(sys.argv) > 1 else "shell"
if what in ("shell", "sphere"):
    cfg = P["configs"]["C2" if what == "shell" else "C3"]
    d = res.setdefault(what, {})
    for a in cfg["a"]:
        key = f"{a:.6g}"
        if key in d or time.time() - t0 > 420:
            continue
        mk = lambda Z: SWalk(-1, a, cfg["rmax"], Vmodel(what, Z))
        Zcr, track = track_to_critical(mk, 160.0, log=lambda s: None, sigma0=-0.3)
        e = {"a": a, "Zcr": Zcr, "track": track}
        if what == "shell":
            w = mk(170.0)
            e["E170"] = level(w, -0.66)[0]
        d[key] = e
        print(what, a, Zcr, e.get("E170"), f"{time.time() - t0:.0f}s", flush=True)
        save()
    runs = sorted(d.values(), key=lambda x: x["a"]) if all(k[0].isdigit() for k in d) else None
    if runs and len(runs) == len(cfg["a"]):
        A = np.array([x["a"] for x in runs])
        Zc = np.array([x["Zcr"] for x in runs])
        ex = {}
        for n in (3, 4):
            c = np.polyfit(A[:n], Zc[:n], 2)
            ex[f"Zcr_quad{n}"] = float(c[-1])
        c = np.polyfit(A[:2], Zc[:2], 1)
        ex["Zcr_lin2"] = float(c[-1])
        if what == "shell":
            Es = np.array([x["E170"] for x in runs])
            for n in (3, 4):
                ex[f"E170_quad{n}"] = float(np.polyval(np.polyfit(A[:n], Es[:n], 2), 0))
        res[what + "_extrap"] = ex
        print(what, ex)
        save()
elif what == "cont":
    c = res.setdefault("continuum", {})
    if "shell" not in c:
        c["shell"] = dirac_critical_Z(lambda Z: Vmodel("shell", Z), 170.0, 185.0)
        c["shell_E170"] = dirac_bound(Vmodel("shell", 170.0), -0.9, -0.4)
    if "sphere" not in c:
        c["sphere"] = dirac_critical_Z(lambda Z: Vmodel("sphere", Z), 160.0, 185.0)
        c["sphere_R_fm"] = R_sphere_fm(c["sphere"])
    print(c)
    save()

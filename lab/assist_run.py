"""assist_run.py - production runs of the assist_ study (pre-registered in lab/out/assist_prereg.json).
usage: python assist_run.py A | O | specA | specO | wkb
  A     : weak Sauter pulse, f(p=0) on the (gamma, E1, eps) grid, rule at 4 lattice spacings + continuum Dirac
  O     : weak oscillating field, same
  specA : spectra f(p) and pair numbers N at E1 = 0.1 (both fields, strong only, weak only)
  specO : same for the oscillating field
  wkb   : WKB spectra / N for the Sauter configurations (semiclassical comparison, report only)
Writes lab/out/assist_<part>.json (checkpointed after every configuration)."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assist_core import Field, walk_pairs, dirac_pairs, wkb_pairs, sauter_pulse_exact, OUT, PI

pre = json.load(open(os.path.join(OUT, "assist_prereg.json")))
part = sys.argv[1]
PATH = os.path.join(OUT, f"assist_{part}.json")
res = json.load(open(PATH)) if os.path.exists(PATH) else {}
t0 = time.time()
GS = pre["A_sauter"]["gs"]
THETAS = pre["A_sauter"]["thetas"]


def save():
    json.dump(res, open(PATH, "w"), indent=1)


def point(F):
    r = dict(dirac=float(dirac_pairs(F, [0.0], 0.01)[0]))
    for th in THETAS:
        r[f"th{th}"] = float(walk_pairs(F, [0.0], th)[0])
    r["rich"] = (4 * r["th0.025"] - r["th0.05"]) / 3
    return r


if part in ("A", "O"):
    shape = "sauter" if part == "A" else "cos"
    P = pre["A_sauter"] if part == "A" else pre["O_cos"]
    epss = [0.0, 0.01, -0.01, 0.005, -0.005, 0.1] if part == "A" else [0.0, 0.001, 0.01]
    for g in P["gammas"]:
        for E1 in P["E1_fit"]:
            key = f"g{g}_E{E1}"
            if key in res:
                continue
            r = {}
            for eps in epss:
                if eps == 0.0:
                    F = Field(E1, GS * E1)
                else:
                    F = Field(E1, GS * E1, eps * E1, g * E1, shape)
                r[f"eps{eps}"] = point(F)
            if part == "A" and E1 == 0.0625:
                F = Field(E1, GS * E1, 0.01 * E1, g * E1, shape)
                r["dt_check_eps0.01"] = float(dirac_pairs(F, [0.0], 0.005)[0])
            res[key] = r
            save()
            print(part, key, {k: (v["dirac"], v["rich"] / v["dirac"]) for k, v in r.items() if isinstance(v, dict)},
                  f"{time.time()-t0:.0f}s", flush=True)

if part in ("specA", "specO"):
    shape = "sauter" if part == "specA" else "cos"
    E1 = 0.1
    p = np.arange(-4, 4.00001, 0.01)
    gams = pre["A_sauter"]["gammas"] + [20.0] if shape == "sauter" else pre["O_cos"]["gammas"]
    epss = pre["A_sauter"]["eps_finite"] if shape == "sauter" else pre["O_cos"]["eps"]
    if "strong_only" not in res:
        F = Field(E1, GS * E1)
        r = dict(dirac=dirac_pairs(F, p, 0.01).tolist())
        for th in (0.1, 0.05, 0.025):
            r[f"th{th}"] = walk_pairs(F, p, th).tolist()
        res["strong_only"] = r; res["p"] = p.tolist(); save()
    for eps in epss:
        for g in gams:
            key = f"eps{eps}_g{g}"
            if key in res:
                continue
            r = {}
            for lab, F in (("both", Field(E1, GS * E1, eps * E1, g * E1, shape)),
                           ("weak_only", Field(0.0, GS * E1, eps * E1, g * E1, shape, E1_shape="none"))):
                d = dict(dirac=dirac_pairs(F, p, 0.01).tolist())
                for th in (0.1, 0.05, 0.025):
                    d[f"th{th}"] = walk_pairs(F, p, th).tolist()
                if lab == "weak_only" and shape == "sauter":
                    d["exact"] = sauter_pulse_exact(p, eps * E1, 1 / (g * E1)).tolist()
                r[lab] = d
            res[key] = r; save()
            print(part, key, f"{time.time()-t0:.0f}s", flush=True)

if part == "wkb":
    E1 = 0.1
    p = np.arange(-4, 4.00001, 0.1)
    for eps in pre["A_sauter"]["eps_finite"]:
        for g in pre["A_sauter"]["gammas"] + [20.0]:
            key = f"eps{eps}_g{g}"
            if key in res:
                continue
            F = Field(E1, GS * E1, eps * E1, g * E1, "sauter")
            f = np.array([wkb_pairs(F, pp) for pp in p])
            res[key] = dict(p=p.tolist(), f=f.tolist(), N=float(np.trapezoid(f, p) / (2 * PI)))
            save()
            print(part, key, res[key]["N"], f"{time.time()-t0:.0f}s", flush=True)
print("done", f"{time.time()-t0:.0f}s")

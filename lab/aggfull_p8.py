"""aggfull_p8.py - phase 8 (docs/AGGLOMERATO_COMPLETO.md s.17): the hybrid atom with the two proposed patches applied
(on MY copies, nucleo1s files untouched): zoom chain with outer cell a_out = theta_A/4 (aggfull_p7_glue.GlueT), T level
(aggfull_relearnT.pkl via nucleo1s_f8_core.CoarseT, read only) + the phase-6 field-motion coupling (large component),
cusp glue of nucleo1s_f8_core.Glue8. alpha^4 from two Z at the SAME geometry in a0/Z, using only the existing levels:
  (Z, B=512) and (2Z, B=256) for Z = 1, 3, 5  (cells (Z/8) a0/Z and the same for 2Z)
delta = (B2/Z2^2 - B1/Z1^2) / (Z2^2 - Z1^2), compared with the same formula on the Dirac levels (exam only).
Usage: timeout 590 python lab/aggfull_p8.py run   (repeat until 'all done') -> lab/out/aggfull_p8.json, aggfull_p8_err.json
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nucleo1s_f8_core as n8                      # read only
import aggfull_p6_glue as g6
import aggfull_p7_glue as g7
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
TH = n8.TH_A
DEV = os.environ.get("AGGFULL_DEV") == "1"


def run(cfg, resfile, budget=480):
    res = json.load(open(resfile)) if os.path.exists(resfile) else {}
    T0 = time.time()
    coarse, coefs = {}, {}
    for v in cfg["variants"]:
        for Zl in cfg["Z"]:
            for J0 in v["J0"][str(Zl)] if isinstance(v["J0"], dict) and str(Zl) in v["J0"] else (v["J0"][Zl] if isinstance(v["J0"], dict) else v["J0"]):
                for Z, B in ((Zl, 512), (2 * Zl, 256)):
                    for n in cfg["states"]:
                        key = f"{v['name']}_P{Zl}_J{J0}_Z{Z}_n{n}"
                        if key in res:
                            continue
                        if time.time() - T0 > budget:
                            print("time budget: run again", flush=True); return res, False
                        if B not in coarse:
                            coarse[B] = n8.CoarseT(B, "learned", ndir=120); coefs[B] = g6.lc_coeffs(B)
                        t = time.time()
                        g = g7.GlueT(coarse[B], Z, J0, cfg["rmax_a0_over_Z"] * n8.A0 / Z, coefs[B], TH * v["theta_factor"],
                                     cfg["zoom_levels"], cfg["n_per_level"])
                        r = g.solve(n)
                        res[key] = dict(B=None if r is None else r["B"], nodes=None if r is None else r["nodes"],
                                        cells_inner=g.mesh.N, sec=time.time() - t)
                        json.dump(res, open(resfile, "w"), indent=1, default=float)
                        print(key, res[key]["B"], f"{res[key]['sec']:.1f}s", flush=True)
    return res, True


def errors(res, cfg):
    out = {}
    for v in cfg["variants"]:
        for Zl in cfg["Z"]:
            for J0 in v["J0"][str(Zl)] if isinstance(v["J0"], dict) and str(Zl) in v["J0"] else (v["J0"][Zl] if isinstance(v["J0"], dict) else v["J0"]):
                for n in cfg["states"]:
                    b1 = res.get(f"{v['name']}_P{Zl}_J{J0}_Z{Zl}_n{n}", {}).get("B")
                    b2 = res.get(f"{v['name']}_P{Zl}_J{J0}_Z{2 * Zl}_n{n}", {}).get("B")
                    if b1 is None or b2 is None:
                        continue
                    d = n8.delta_pair(b1, b2, Zl, 2 * Zl)
                    out[f"{v['name']}_P{Zl}_J{J0}_{n}s"] = d / n8.delta_dirac(n, Zl, 2 * Zl) - 1
    return out


if __name__ == "__main__":
    pre = os.path.join(OUT, "aggfull_p8_prereg.json")
    cfg = json.load(open(pre))["config"] if os.path.exists(pre) and not DEV else dict(
        variants=[dict(name="dev", theta_factor=0.25, J0={1: [32], 3: [10, 21], 5: [6, 12]})], Z=[3, 5], states=[1, 2], rmax_a0_over_Z=90, zoom_levels=12, n_per_level=16)
    res, done = run(cfg, os.path.join(OUT, ("dev_" if DEV else "") + "aggfull_p8.json"))
    E = errors(res, cfg)
    for k, v in E.items():
        print(k, f"{v:+.4f}")
    if done and not DEV:
        json.dump(E, open(os.path.join(OUT, "aggfull_p8_err.json"), "w"), indent=1)
        print("all done")

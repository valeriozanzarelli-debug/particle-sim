"""aggfull_p7_glue.py - phase 7 (docs/AGGLOMERATO_COMPLETO.md s.16): where does the -10.5 % plateau of the 1s come from?
The zoom chain of nucleo1s_f8_core.Glue8 (imported READ ONLY) uses outer chain cells a_out = theta_A (fixed in lambda_C)
for both Z: at Z = 2 the chain is twice as coarse in units a0/Z, so its lattice error (~ (Z alpha theta)^2 relative, i.e.
~ Z^4) does NOT cancel in the two-Z difference and mimics an alpha^4 term. Here the chain's cell is a free parameter:
  GlueT(c, Z, J0, rmax, coef, a_out, zoom_levels)  = aggfull_p6_glue.GlueLC with the mesh rebuilt with cell a_out.
Usage: timeout 590 python lab/aggfull_p7_glue.py run   (repeat until 'all done') -> lab/out/aggfull_p7_glue.json
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nucleo1s_f8_core as n8                     # read only
from nucleo1s_core import mesh_zoom               # read only
import aggfull_p6_glue as g6
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
TH = n8.TH_A


class GlueT(g6.GlueLC):
    def __init__(self, c, Z, J0, rmax, coef, a_out=TH, zoom_levels=12, npl=16):
        super().__init__(c, Z, J0, rmax, coef, zoom_levels, npl)
        self.mesh = mesh_zoom(a_out, (J0 + 0.5) * c.h, zoom_levels, npl)


def solve_one(coarse, coef, Z, J0, n, a_out, zl, rmax_a0=90, npl=16):
    g = GlueT(coarse, Z, J0, rmax_a0 * n8.A0 / Z, coef, a_out, zl, npl)
    r = g.solve(n)
    return None if r is None else r["B"]


def run(cfg, resfile, budget=480):
    res = json.load(open(resfile)) if os.path.exists(resfile) else {}
    T0 = time.time()
    coarse, coefs = {}, {}
    for v in cfg["variants"]:
        for J0 in cfg["J0"]:
            for Z, B in ((1, 512), (2, 256)):
                for n in cfg["states"]:
                    key = f"{v['name']}_J{J0}_Z{Z}_n{n}"
                    if key in res:
                        continue
                    if time.time() - T0 > budget:
                        print("time budget: run again", flush=True); return res, False
                    if B not in coarse:
                        coarse[B] = n8.CoarseT(B, "learned", ndir=120); coefs[B] = g6.lc_coeffs(B)
                    t = time.time()
                    a_out = TH * v["theta_factor"][Z - 1]
                    res[key] = dict(B=solve_one(coarse[B], coefs[B] if v["lc"] else None, Z, J0, n, a_out, v["zoom"], npl=v.get("npl", 16)),
                                    a_out=a_out, sec=time.time() - t)
                    json.dump(res, open(resfile, "w"), indent=1, default=float)
                    print(key, res[key]["B"], f"{res[key]['sec']:.1f}s", flush=True)
    return res, True


def errors(res, cfg):
    out = {}
    for v in cfg["variants"]:
        for J0 in cfg["J0"]:
            for n in cfg["states"]:
                b1, b2 = (res.get(f"{v['name']}_J{J0}_Z{Z}_n{n}", {}).get("B") for Z in (1, 2))
                if b1 is not None and b2 is not None:
                    out[f"{v['name']}_J{J0}_{n}s"] = n8.delta_pair(b1, b2) / n8.delta_dirac(n) - 1
    return out


if __name__ == "__main__":
    if sys.argv[1] == "run":
        P = json.load(open(os.path.join(OUT, "aggfull_p7_prereg.json")))["config"]["glue"]
        cfg = dict(variants=[dict(name=v[0], theta_factor=v[1], zoom=v[2], npl=v[3], lc=P["lc"]) for v in P["variants"]],
                   J0=P["J0"], states=P["states"])
        res, done = run(cfg, os.path.join(OUT, "aggfull_p7_glue.json"))
        if done:
            json.dump(errors(res, cfg), open(os.path.join(OUT, "aggfull_p7_glue_err.json"), "w"), indent=1)
            print("all done")
    if sys.argv[1] == "dev":
        cfg = dict(variants=[dict(name=s.split(":")[0], theta_factor=[float(x) for x in s.split(":")[1].split("/")],
                                  zoom=int(s.split(":")[2]), npl=int(s.split(":")[3]) if s.count(":") > 2 else 16, lc=True) for s in sys.argv[3].split(",")],
                   J0=[int(x) for x in sys.argv[2].split(",")], states=[1, 2])
        res, done = run(cfg, os.path.join(OUT, "dev_aggfull_p7_glue.json"))
        for k, v in errors(res, cfg).items():
            print(k, f"{v:+.4f}")

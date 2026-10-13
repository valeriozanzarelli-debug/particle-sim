"""nucleo1s_f8_run.py - phase 8 runs (lab/out/nucleo1s_f8_prereg.json) -> lab/out/nucleo1s_f8_results.json
(re-run until 'all done'), then the verdicts -> lab/out/nucleo1s_f8_judged.json"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_f8_core import CoarseT, Glue8, delta_pair, delta_dirac, A0  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "nucleo1s_f8_prereg.json")))
C = P["config"]
RES = os.path.join(OUT, "nucleo1s_f8_results.json")
res = json.load(open(RES)) if os.path.exists(RES) else {}
T0 = time.time()
coarse = {}
for band in C["bands"]:
    for J0 in C["J0"]:
        for pz in C["pair"]:
            Z, B = pz["Z"], pz["B"]
            for n in C["states"]:
                k = f"{band}_J{J0}_Z{Z}_n{n}"
                if k in res:
                    continue
                if time.time() - T0 > 480:
                    print("time budget: run again"); sys.exit(0)
                if (band, B) not in coarse:
                    coarse[(band, B)] = CoarseT(B, band, ndir=C["ndir"])
                t = time.time()
                g = Glue8(coarse[(band, B)], Z, J0, C["rmax_a0_over_Z"] * A0 / Z, C["zoom_levels"], C["n_per_level"])
                r = g.solve(n)
                r.update(cells_inner=g.mesh.N, cells_outer=g.N - J0, sec=time.time() - t, imag=g.imag_ratio)
                res[k] = r
                json.dump(res, open(RES, "w"), indent=1, default=float)
                print(k, f"{r['rel']:+.3e}", r["nodes"], f"{r['sec']:.1f}s", flush=True)
J = {}
for band in C["bands"]:
    for J0 in C["J0"]:
        for n in C["states"]:
            d = delta_pair(res[f"{band}_J{J0}_Z1_n{n}"]["B"], res[f"{band}_J{J0}_Z2_n{n}"]["B"])
            J[f"{band}_J{J0}_{n}s"] = d / delta_dirac(n) - 1
E = P["exams"]
err = lambda b, j, n: J[f"{b}_J{j}_{n}s"]
J["H1"] = bool(all(abs(err("learned", 3, n)) <= E["H1"]["tol"] for n in (1, 2)))
J["H2"] = bool(all(abs(err("learned", j, n) - err("exact", j, n)) <= E["H2"]["tol"] for j in C["J0"] for n in (1, 2)))
mono = all(abs(err("exact", C["J0"][i + 1], n)) <= abs(err("exact", C["J0"][i], n)) + E["H3"]["wiggle"]
           for i in range(len(C["J0"]) - 1) for n in (1, 2))
J["H3"] = bool(mono and abs(err("exact", 24, 1)) <= E["H3"]["tol1"] and abs(err("exact", 24, 2)) <= E["H3"]["tol2"])
J["H4"] = bool(abs(err("learned", 24, 1)) <= E["H4"]["tol1"] and abs(err("learned", 24, 2)) <= E["H4"]["tol2"])
json.dump(J, open(os.path.join(OUT, "nucleo1s_f8_judged.json"), "w"), indent=1)
for k, v in J.items():
    print(k, v if isinstance(v, bool) else f"{v:+.4f}")
print("all done")

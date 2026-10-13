"""aggfull_p7_chain.py - phase 7: the zoom chain ALONE (no coarse level, no glue): the zig-zag radial walk of
nucleo1s (zz_regular, read only) on mesh_zoom(theta, rmax, levels, npl) out to rmax = 30 a0/Z; eigenvalue by shooting
(the regular solution must vanish at rmax). alpha^4 by the same two-Z difference as the glue test.
Usage: timeout 590 python lab/aggfull_p7_chain.py run|dev ...  -> lab/out/aggfull_p7_chain.json"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_hybrid import zz_regular, TH_A, A0          # read only
from nucleo1s_core import mesh_zoom, dirac_binding, ALPHA  # read only
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def chain_level(Z, n, theta, levels, npl, rmax_a0=30, span=0.02, npts=13):
    Za = Z * ALPHA
    mesh = mesh_zoom(theta, rmax_a0 * A0 / Z, levels, npl)
    iref = int(np.searchsorted(mesh.r, 0.3 * A0 / Z))

    def f(E):
        r, g, _ = zz_regular(mesh, -1, Za, 1 + E)
        x = g[-1] / g[iref]                               # zz_regular phases g at its maximum: use a ratio for the sign
        return float(np.real(x) / (1 + abs(x)))
    Ed = -dirac_binding(n, -1, Za)
    m = np.tan(theta) / theta                            # the walk's kinetic mass (property of the fine rule)
    Es = np.linspace(Ed * m * (1 + span), Ed * m * (1 - span), npts)
    fs = [f(E) for E in Es]
    br = [i for i in range(npts - 1) if np.sign(fs[i]) != np.sign(fs[i + 1])]
    if not br:
        return None, mesh.N
    i = min(br, key=lambda i: abs(0.5 * (Es[i] + Es[i + 1]) - Ed * m))
    E = brentq(f, Es[i], Es[i + 1], xtol=1e-22, rtol=1e-15, maxiter=200)
    return -E, mesh.N


def delta_err(B1, B2, n):
    d = (B2 / 4 - B1) / 3
    dd = (dirac_binding(n, -1, 2 * ALPHA) / 4 - dirac_binding(n, -1, ALPHA)) / 3
    return d / dd - 1


def run(variants, states, resfile, budget=500):
    res = json.load(open(resfile)) if os.path.exists(resfile) else {}
    T0 = time.time()
    for name, tf, lv, npl in variants:
        for n in states:
            for Z in (1, 2):
                key = f"{name}_Z{Z}_n{n}"
                if key in res:
                    continue
                if time.time() - T0 > budget:
                    print("time budget: run again", flush=True); return res, False
                t = time.time()
                B, N = chain_level(Z, n, TH_A * tf, lv, npl)
                res[key] = dict(B=B, cells=N, theta=TH_A * tf, levels=lv, npl=npl,
                                rel=None if B is None else B / dirac_binding(n, -1, Z * ALPHA) - 1, sec=time.time() - t)
                json.dump(res, open(resfile, "w"), indent=1, default=float)
                print(key, res[key]["rel"], N, f"{res[key]['sec']:.0f}s", flush=True)
    return res, True


def errors(res, variants, states):
    out = {}
    for name, *_ in variants:
        for n in states:
            a, b = res.get(f"{name}_Z1_n{n}", {}).get("B"), res.get(f"{name}_Z2_n{n}", {}).get("B")
            if a is not None and b is not None:
                out[f"{name}_{n}s"] = delta_err(a, b, n)
    return out


def run_multi(cfg, resfile, budget=500):
    res = json.load(open(resfile)) if os.path.exists(resfile) else {}
    T0 = time.time()
    for name, tf, lv, npl in cfg["variants"]:
        for n in cfg["states"]:
            for Z in cfg["Z"]:
                key = f"{name}_Z{Z}_n{n}"
                if key in res:
                    continue
                if time.time() - T0 > budget:
                    print("time budget: run again", flush=True); return res, False
                t = time.time()
                B, N = chain_level(Z, n, TH_A * tf, lv, npl, rmax_a0=cfg["rmax_a0_over_Z"])
                res[key] = dict(B=B, cells=N, theta=TH_A * tf, rel=None if B is None else B / dirac_binding(n, -1, Z * ALPHA) - 1,
                                sec=time.time() - t)
                json.dump(res, open(resfile, "w"), indent=1, default=float)
                print(key, res[key]["rel"], f"{res[key]['sec']:.0f}s", flush=True)
    return res, True


if __name__ == "__main__":
    if sys.argv[1] == "run":
        cfg = json.load(open(os.path.join(OUT, "aggfull_p7_prereg.json")))["config"]["chain"]
        res, done = run_multi(cfg, os.path.join(OUT, "aggfull_p7_chain.json"))
        if done:
            print("all done")
    if sys.argv[1] == "dev":
        V = [(s.split(":")[0], float(s.split(":")[1]), int(s.split(":")[2]), int(s.split(":")[3])) for s in sys.argv[2].split(",")]
        res, _ = run(V, [1], os.path.join(OUT, "dev_aggfull_p7_chain.json"))
        for k, v in errors(res, V, [1]).items():
            print(k, f"{v:+.4f}")

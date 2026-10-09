"""big_b.py - phase B, the big step: surfers (lam = 0) reinforce the links they use; does a geometry build itself?
Pre-registered in lab/out/big_prereg_B.json.
usage: python big_b.py SEED   (1, 2, 3; ~3 min each)  -> lab/out/big_b_seed{SEED}.json/.npz
       python big_b.py ref    (uniform w = 0.3 reference geometry)  -> lab/out/big_b_ref.json/.npz
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import big_core as B
from pre_core import LatticeWalk

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
n, TH, T = 24, 0.125, 5000
RE = dict(eta=0.3, mu=2, wbar=0.3, wmin=0.02, K=10)


def geometry(wk):
    I, gap = B.vacuum_mi(wk)
    r = B.read_shape(I)
    sp = B.low_spectrum(wk)
    om = sp["zero"]
    out = dict(pieces=r["pieces"], sizes=r["sizes"], holes=r["holes"], rel=r["rel"], D_small=r["D_small"],
               D_large=r["D_large"], diam=r["diam"], gap=float(gap),
               omega_zero=om.tolist(), omega_pi=sp["pi"].tolist(), ratios=(om / om[0]).tolist(),
               ipr_frac=sp["ipr_frac"].tolist(), ipr_frac_lowest=float(sp["ipr_frac"][0]))
    arr = dict(X=r["X"], sites=r["sites"], rho_low=sp["rho_low"])
    return out, arr


if __name__ == "__main__":
    t0 = time.time()
    arg = sys.argv[1]
    if arg == "ref":
        wk = LatticeWalk(n, TH, periodic=False, w=[np.full((n, n), RE["wbar"])] * 2)
        g, a = geometry(wk)
        print("ref", g, f"{time.time()-t0:.0f}s", flush=True)
        with open(os.path.join(OUT, "big_b_ref.json"), "w") as f:
            json.dump(g, f, indent=1, default=float)
        np.savez(os.path.join(OUT, "big_b_ref.npz"), **a)
        sys.exit()
    s = int(arg)
    ps, ss = 300 + s, 400 + s
    res, arr = {"seed": s}, {}
    for name, start, rf in (("reinf_uniform", "uniform", RE), ("plain_uniform", "uniform", None),
                            ("reinf_psi2", "psi2", RE)):
        r = B.run(n, TH, 16, start, T, psi_seed=ps, surf_seed=ss, lam=0.0, every=10, reinf=rf, w0=RE["wbar"],
                  vort_every=T, snap=(0, 1000, T))
        k10 = max(1, len(r["H"]) // 10)
        e = dict(H0=r["H0"], H_last10=float(r["H"][-k10:].mean()), maxH=float(r["H"].max()), F=r["F"],
                 tau=r["tau"], windings_start_end=r["V"].tolist())
        if rf:
            wi = B.interior(r["walk"].W, n)
            e.update(gini_end=B.gini(wi), top10_end=B.top10(wi), w_max=float(wi.max()),
                     frac_at_wmin=float(np.mean(wi <= RE["wmin"] + 1e-9)), frac_at_1=float(np.mean(wi >= 1 - 1e-9)))
            arr[name + "_G"] = r["Gini"]
            arr[name + "_W"] = np.array(r["walk"].W)
            arr[name + "_W1000"] = np.array(r["snaps"][1000])
        arr[name + "_H"] = r["H"]; arr["t"] = r["t"]
        res[name] = e
        print(name, e, f"{time.time()-t0:.0f}s", flush=True)
        if name == "reinf_uniform":
            g, a = geometry(r["walk"])
            res["shape_reinf_uniform"] = g
            arr.update({"shape_" + k: v for k, v in a.items()})
            print("shape", g, f"{time.time()-t0:.0f}s", flush=True)
    res["time_s"] = time.time() - t0
    with open(os.path.join(OUT, f"big_b_seed{s}.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    np.savez(os.path.join(OUT, f"big_b_seed{s}.npz"), **arr)
    print("done", res["time_s"])

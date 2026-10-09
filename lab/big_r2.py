"""big_r2.py - R2: is the fake (jump-noise) mixing gone with lam = 0 surfers at small theta?
M = 1 control vs M = 16, 24x24 box, uniform start, 3 surfer seeds (curves averaged). Pre-registered in
lab/out/big_prereg_R.json.  usage: python big_r2.py THETA   -> lab/out/big_r2_theta{THETA}.json, .npz
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import big_core as B
from pre_core import fit_tau

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

if __name__ == "__main__":
    th = float(sys.argv[1])
    T = int(round(200 / th)); every = int(round(2.5 / th))
    t0 = time.time()
    res, arr = {"theta": th, "T": T}, {}
    for M in (1, 16):
        Hs = []
        for s in (201, 202, 203):
            r = B.run(24, th, M, "uniform", T, psi_seed=101, surf_seed=s, lam=0.0, every=every,
                      vort_every=T // 4)
            Hs.append(r["H"])
            print(th, M, s, "H0", r["H0"], "ratio", r["ratio"], "tau", r["tau"], "wind", r["V"].tolist(),
                  f"{time.time()-t0:.0f}s", flush=True)
        H = np.mean(Hs, 0)
        F = r["F"]
        res[f"M{M}"] = dict(H0=float(H[0]), HT=float(H[-1]), ratio=float(H[-1] / H[0]), tau=fit_tau(r["t"], H, F),
                            tau_seeds=[fit_tau(r["t"], h, F) for h in Hs], ratio_seeds=[float(h[-1] / h[0]) for h in Hs],
                            F=F, omega=r["omega"].tolist(), windings=r["V"].tolist())
        arr[f"M{M}"] = H; arr["t"] = r["t"]
    a, b = res["M1"], res["M16"]
    c1 = a["tau"] >= 3 * b["tau"]
    c2 = a["ratio"] > 0.10 and b["ratio"] < 0.03
    res["verdict"] = {"tau_M1>=3tau_M16": bool(c1), "M1>10%_and_M16<3%": bool(c2), "R2_pass": bool(c1 or c2),
                      "tau_ratio": float(a["tau"] / b["tau"]) if np.isfinite(b["tau"]) and b["tau"] > 0 else None}
    res["time_s"] = time.time() - t0
    print(json.dumps(res, indent=1, default=float))
    with open(os.path.join(OUT, f"big_r2_theta{th}.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    np.savez(os.path.join(OUT, f"big_r2_theta{th}.npz"), **arr)

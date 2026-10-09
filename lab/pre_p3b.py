"""pre_p3b.py - POST-HOC follow-up of P3 (declared). With wbar = 0.8 and the cap w <= 1 the top 10% of links can
carry at most 1/(10*0.8) = 12.5% of the weight, so the pre-registered highway test (> 30%) was impossible by
construction. Here: low budget wbar = 0.3, strong superlinear feedback (eta = 0.3, mu = 2), uniform and |psi|^2
starts, plus the no-feedback reference at the same wbar. Writes lab/out/pre_p3b.json, pre_p3b.npz (~1 min)."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pre_p3 as P

if __name__ == "__main__":
    t0 = time.time()
    P.WBAR = 0.3
    res, arr = {}, {}
    for name, eta, start, shp in (("low_strong_uniform", 0.3, "uniform", True), ("low_plain_uniform", 0.0, "uniform", True),
                                  ("low_strong_psi2", 0.3, "psi2", False)):
        r = P.run(eta, start, do_shape=shp, mu=2.0)
        k10 = max(1, len(r["H"]) // 10)
        res[name] = dict(wbar=0.3, eta=eta, mu=2.0, H0=float(r["H"][0]), H_last10=float(r["H"][-k10:].mean()),
                         maxH=float(r["H"].max()), F=r["F"], gini_end=r["gini_end"], top10_end=r["top10_end"],
                         vort_first=float(r["V"][:3].mean()), vort_last=float(r["V"][-3:].mean()),
                         shapes={str(k): v for k, v in r["shapes"].items()})
        arr[name + "_H"] = r["H"]; arr[name + "_G"] = r["Gini"]; arr["t"] = r["t"]
        if shp:
            for k, W in r["Wsnap"].items():
                arr[f"{name}_W{k}"] = W
        print(name, {k: v for k, v in res[name].items() if k != "shapes"}, f"{time.time()-t0:.0f}s", flush=True)
    S, Pl = res["low_strong_uniform"], res["low_plain_uniform"]
    res["checks"] = {"pinning": bool(S["H_last10"] > 3 * Pl["H_last10"] and S["H_last10"] > 3 * S["F"]),
                     "highways(gini>+0.2 & top10>0.3)": bool(S["gini_end"] > 0.2 and S["top10_end"] > 0.3),
                     "equivariance_psi2(maxH<=3F)": bool(res["low_strong_psi2"]["maxH"] <= 3 * res["low_strong_psi2"]["F"])}
    print(res["checks"], time.time() - t0)
    with open(os.path.join(P.OUT, "pre_p3b.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    np.savez(os.path.join(P.OUT, "pre_p3b.npz"), **arr)

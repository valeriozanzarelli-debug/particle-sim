"""pre_p3.py - P3: reinforcement vs mixing (toy). Surfers reinforce the links they cross (Physarum-like),
rule and its 4 parameters fixed in lab/out/pre_prereg.json BEFORE the run.
  every K ticks: J_l = surfer crossings of link l; w <- (1-eta) w + eta wbar J_l/mean(J); clip [w_min,1];
  rescale mean(w) = wbar (3 passes). Box edge links stay 0.
Measures H(t) with/without reinforcement, highways (Gini), windings, and the geo_core reconstruction of the
vacuum MI of U(w) at t = 0, 1000, 2000, 3000. Writes lab/out/pre_p3.json, pre_p3.npz (~1-2 min).
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import geo_core as G
from pre_core import LatticeWalk, modes_state, surfer_tick, sample_sites, hbar, windings, SEED

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
ETA, WBAR, WMIN, K = 0.05, 0.8, 0.02, 10
n, M, N, T, CELL = 12, 16, 20000, 3000, 3


def interior(W):
    return np.concatenate([W[0][:-1, :].ravel(), W[1][:, :-1].ravel()])


def set_interior(vals):
    W0 = np.zeros((n, n)); W1 = np.zeros((n, n))
    m = (n - 1) * n
    W0[:-1, :] = vals[:m].reshape(n - 1, n)
    W1[:, :-1] = vals[m:].reshape(n, n - 1)
    return [W0, W1]


def gini(x):
    x = np.sort(x); k = len(x)
    return float((2 * np.arange(1, k + 1) - k - 1).dot(x) / (k * x.sum()))


def top10(x):
    x = np.sort(x)[::-1]
    return float(x[:max(1, len(x) // 10)].sum() / x.sum())


def shape_of(wk):
    C, gap, nocc = G.vacuum(wk.matrix())
    I = G.mutual_information(C, n * n)
    nc, rs = G.reconstruct_from_I(I, maxdim=2)
    sizes = sorted([len(r["sites"]) for r in rs], reverse=True)
    big = max(rs, key=lambda r: len(r["sites"]))
    return dict(pieces=int(nc), sizes=sizes[:6], betti=[int(b) for b in big["betti"]], D=float(big["dim"]),
                gap=float(gap), all_betti=[[int(b) for b in r["betti"]] for r in rs if len(r["sites"]) >= 8])


def run(eta, start, seed=SEED, do_shape=False, mu=1.0):
    wk = LatticeWalk(n, 0.5, periodic=False, w=[np.full((n, n), WBAR)] * 2)
    psi, om = modes_state(wk, M)
    rng = np.random.default_rng(seed)
    rho = np.sum(np.abs(psi)**2, 0)
    X = sample_sites(rng, N, rho if start == "psi2" else np.ones((n, n)))
    counts = [np.zeros((n, n)), np.zeros((n, n))]
    ts, H, V, Gs, shapes, Wsnap = [], [], [], [], {}, {}
    for t in range(T + 1):
        if t % 10 == 0:
            rho = np.sum(np.abs(psi)**2, 0)
            h, nc = hbar(X, rho, CELL)
            ts.append(t); H.append(h); Gs.append(gini(interior(wk.W)))
            if t % 100 == 0:
                V.append(windings(psi))
        if do_shape and t % 1000 == 0:
            shapes[t] = shape_of(wk)
            Wsnap[t] = np.array(wk.W)
            print("  shape t", t, shapes[t], flush=True)
        if t == T:
            break
        psi, X = surfer_tick(wk, psi, X, rng, "w", counts)
        if eta > 0 and (t + 1) % K == 0:
            J = interior(counts)
            w = interior(wk.W)
            Jm = J**mu
            w = (1 - eta) * w + eta * WBAR * Jm / max(Jm.mean(), 1e-12)
            for _ in range(3):
                w = np.clip(w, WMIN, 1.0)
                w *= WBAR / w.mean()
            w = np.clip(w, WMIN, 1.0)
            wk.set_w(set_interior(w))
            counts = [np.zeros((n, n)), np.zeros((n, n))]
        elif (t + 1) % K == 0:
            counts = [np.zeros((n, n)), np.zeros((n, n))]
    F = (nc - 1) / (2 * N)
    wi = interior(wk.W)
    return dict(t=np.array(ts), H=np.array(H), V=np.array(V), Gini=np.array(Gs), F=F, gini_end=gini(wi),
                top10_end=top10(wi), shapes=shapes, Wsnap=Wsnap, w_end=wk.W)


if __name__ == "__main__":
    t0 = time.time()
    res, arr = {}, {}
    for name, eta, start, shp in (("reinf_uniform", ETA, "uniform", True), ("plain_uniform", 0.0, "uniform", False),
                                  ("reinf_psi2", ETA, "psi2", False), ("plain_psi2", 0.0, "psi2", False)):
        r = run(eta, start, do_shape=shp)
        k10 = max(1, len(r["H"]) // 10)
        res[name] = dict(H0=float(r["H"][0]), H_last10=float(r["H"][-k10:].mean()), maxH=float(r["H"].max()),
                         F=r["F"], gini_end=r["gini_end"], top10_end=r["top10_end"],
                         vort_first=float(r["V"][:3].mean()), vort_last=float(r["V"][-3:].mean()),
                         shapes={str(k): v for k, v in r["shapes"].items()})
        arr[name + "_H"] = r["H"]; arr[name + "_V"] = r["V"]; arr[name + "_G"] = r["Gini"]
        arr["t"] = r["t"]
        if shp:
            for k, W in r["Wsnap"].items():
                arr[f"{name}_W{k}"] = W
        print(name, {k: v for k, v in res[name].items() if k != "shapes"}, f"{time.time()-t0:.0f}s", flush=True)
    R, P = res["reinf_uniform"], res["plain_uniform"]
    g0 = 0.0
    s0, s1 = R["shapes"]["0"], R["shapes"][str(T)]
    shape_changed = (s0["pieces"] != s1["pieces"]) or (s0["betti"][1:] != s1["betti"][1:])
    res["verdict"] = {
        "P3a_pinning": bool(R["H_last10"] > 3 * P["H_last10"] and R["H_last10"] > 3 * R["F"]),
        "P3b_highways": bool(R["gini_end"] > g0 + 0.2 and R["top10_end"] > 0.30),
        "P3d_shape_changed(before P2 check)": bool(shape_changed),
        "P3e_equivariance(reinf psi2 max H <= 3F)": bool(res["reinf_psi2"]["maxH"] <= 3 * res["reinf_psi2"]["F"]),
    }
    # POST-HOC (declared): strong, superlinear Physarum-like feedback (Tero et al.: J^mu, mu > 1)
    for name, start in (("posthoc_strong_uniform", "uniform"), ("posthoc_strong_psi2", "psi2")):
        r = run(0.3, start, do_shape=(start == "uniform"), mu=2.0)
        k10 = max(1, len(r["H"]) // 10)
        res[name] = dict(eta=0.3, mu=2.0, H0=float(r["H"][0]), H_last10=float(r["H"][-k10:].mean()),
                         maxH=float(r["H"].max()), F=r["F"], gini_end=r["gini_end"], top10_end=r["top10_end"],
                         vort_first=float(r["V"][:3].mean()), vort_last=float(r["V"][-3:].mean()),
                         shapes={str(k): v for k, v in r["shapes"].items()})
        arr[name + "_H"] = r["H"]; arr[name + "_V"] = r["V"]; arr[name + "_G"] = r["Gini"]
        if start == "uniform":
            for k, W in r["Wsnap"].items():
                arr[f"{name}_W{k}"] = W
        print(name, {k: v for k, v in res[name].items() if k != "shapes"}, flush=True)
    res["time_s"] = time.time() - t0
    print(json.dumps(res["verdict"], indent=1), res["time_s"])
    with open(os.path.join(OUT, "pre_p3.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    np.savez(os.path.join(OUT, "pre_p3.npz"), **arr)

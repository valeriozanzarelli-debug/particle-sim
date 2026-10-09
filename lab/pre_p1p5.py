"""pre_p1p5.py - P1 (relaxation on the lattice rule itself) and P5 (tau vs size and number of modes).
Pre-registered in lab/out/pre_prereg.json. Writes lab/out/pre_p1p5.json and pre_p1p5.npz.  (~1 min)
"""
import json, os, time, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from pre_core import LatticeWalk, modes_state, surfer_tick, sample_sites, hbar, windings, fit_tau, SEED

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
N = 20000


def run(n, M, start, lam, T, cell=4, seed=SEED, sheets=False, periodic=False, theta=0.5, omega0=0.9, every=5):
    wk = LatticeWalk(n, theta, periodic=periodic)
    psi, om = modes_state(wk, M, omega0=omega0)
    rng = np.random.default_rng(seed)
    rho = np.sum(np.abs(psi)**2, 0)
    xx, yy = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    if start == "psi2":
        prob = rho
    elif start == "uniform":
        prob = np.ones((n, n))
    elif start == "half":
        prob = (xx < n // 2).astype(float)
    elif start == "sheet0":
        prob = ((xx + yy) % 2 == 0).astype(float)
    X = sample_sites(rng, N, prob)
    ts, H, Hs, V = [], [], [], []
    for t in range(T + 1):
        if t % every == 0:
            rho = np.sum(np.abs(psi)**2, 0)
            h, ncell = hbar(X, rho, cell)
            ts.append(t); H.append(h)
            Hs.append(hbar(X, rho, cell, sheets=True)[0] if sheets else np.nan)
            if t % (8 * every) == 0:
                V.append(windings(psi))
        if t < T:
            psi, X = surfer_tick(wk, psi, X, rng, lam)
    F = (ncell - 1) / (2 * N)
    Fs = (2 * ncell - 1) / (2 * N)
    return dict(t=np.array(ts), H=np.array(H), Hs=np.array(Hs), F=F, Fs=Fs, vort=float(np.mean(V)),
                tau=fit_tau(ts, H, F), ratio=float((H[-1]) / H[0]) if H[0] > 0 else float("nan"))


if __name__ == "__main__":
    t0 = time.time()
    res, curves = {}, {}
    verdict = {}
    # ---------------- P1
    n, T = 24, 400
    for lam in ("w", 0.0):
        r = run(n, 16, "psi2", lam, T)
        key = f"P1a_lam{lam}"
        res[key] = dict(maxH=float(r["H"].max()), F=r["F"], vort=r["vort"])
        curves[key] = r["H"]
        verdict[key] = bool(r["H"].max() <= 3 * r["F"])
        print(key, res[key], verdict[key], flush=True)
    for start in ("uniform", "half"):
        for M in (16, 1):
            for lam in ("w", 0.0):
                r = run(n, M, start, lam, T)
                key = f"P1_{start}_M{M}_lam{lam}"
                res[key] = dict(H0=float(r["H"][0]), HT=float(r["H"][-1]), ratio=r["ratio"], tau=r["tau"],
                                F=r["F"], vort=r["vort"])
                curves[key] = r["H"]
                print(key, res[key], flush=True)
    verdict["P1b_uniform"] = res["P1_uniform_M16_lamw"]["ratio"] < 0.1
    verdict["P1b_half"] = res["P1_half_M16_lamw"]["ratio"] < 0.1
    verdict["P1c_control_frozen(M=1 uniform ratio>0.9)"] = res["P1_uniform_M1_lamw"]["ratio"] > 0.9
    for lam in (1.0, 0.0):
        r = run(n, 16, "sheet0", lam, T, sheets=True)
        key = f"P1d_sheet0_lam{lam}"
        res[key] = dict(Hs0=float(r["Hs"][0]), HsT=float(r["Hs"][-1]), ratio=float(r["Hs"][-1] / r["Hs"][0]),
                        Fs=r["Fs"], H_cells_T=float(r["H"][-1]))
        curves[key] = r["Hs"]
        print(key, res[key], flush=True)
    verdict["P1d_lam1_frozen(ratio>=0.5)"] = res["P1d_sheet0_lam1.0"]["ratio"] >= 0.5
    verdict["P1d_lam0_relaxes(ratio<0.1)"] = res["P1d_sheet0_lam0.0"]["ratio"] < 0.1
    # ---------------- P5
    p5 = {}
    for nn, TT in ((16, 400), (24, 400), (32, 600)):
        r = run(nn, 16, "uniform", "w", TT)
        p5[f"n{nn}_M16"] = dict(tau=r["tau"], ratio=r["ratio"], vort=r["vort"], H0=float(r["H"][0]))
        curves[f"P5_n{nn}_M16"] = r["H"]; curves[f"P5_n{nn}_M16_t"] = r["t"]
        print("P5", nn, p5[f"n{nn}_M16"], flush=True)
    for M in (1, 4, 36):
        r = run(24, M, "uniform", "w", 400)
        p5[f"n24_M{M}"] = dict(tau=r["tau"], ratio=r["ratio"], vort=r["vort"], H0=float(r["H"][0]))
        print("P5 M", M, p5[f"n24_M{M}"], flush=True)
    p5["n24_M16"] = p5.get("n24_M16")
    taus = np.array([p5[f"n{k}_M16"]["tau"] for k in (16, 24, 32)])
    if np.all(np.isfinite(taus)):
        p5["exponent_a"] = float(np.polyfit(np.log([16, 24, 32]), np.log(taus), 1)[0])
    res["P5"] = p5
    res["verdict"] = {k: bool(v) for k, v in verdict.items()}
    res["time_s"] = time.time() - t0
    print(json.dumps(res["verdict"], indent=1), "time", res["time_s"])
    with open(os.path.join(OUT, "pre_p1p5.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    np.savez(os.path.join(OUT, "pre_p1p5.npz"), t=np.arange(0, 401, 5), **{k: v for k, v in curves.items()})

"""redo_r3.py [main|fine] - R3: deterministic lattice surfer (monotone CDF map along each line, = integrated j/rho
of the minimal net current) replacing the random jumps of pre_core. Same setup as pre_p1p5 P1.
main: n = 24, th = 0.5, cells 4x4 (pre-registered exams).  fine: n = 96, th = 0.125, cells 16x16 (secondary).
Writes lab/out/redo_r3_{main,fine}.json / .npz.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from pre_core import LatticeWalk, modes_state, sample_sites, hbar, fit_tau, SEED

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
NS = 20000


def cdf_move(X, j, rho0, rho1):
    """move continuous positions X[:, j] along axis j so that a density uniform inside each site, rho0 -> rho1,
    is mapped exactly (per line, monotone). Lines are indexed by the other coordinate (integer part)."""
    n = rho0.shape[0]
    r0 = rho0 if j == 0 else rho0.T          # r[line_coord_j, other] -> make lines rows: R[other, along]
    r1 = rho1 if j == 0 else rho1.T
    R0, R1 = r0.T, r1.T                      # R[line, i_along]
    tot0 = R0.sum(1, keepdims=True); tot1 = R1.sum(1, keepdims=True)
    C0 = np.concatenate([np.zeros((n, 1)), np.cumsum(R0, 1)], 1) / np.maximum(tot0, 1e-300)
    C1 = np.concatenate([np.zeros((n, 1)), np.cumsum(R1, 1)], 1) / np.maximum(tot1, 1e-300)
    line = np.floor(X[:, 1 - j]).astype(int)
    x = X[:, j]
    i = np.clip(np.floor(x).astype(int), 0, n - 1)
    fr = x - i
    c = C0[line, i] + fr * (C0[line, i + 1] - C0[line, i])           # in [0, 1]
    # inverse of C1 on each line: global search on offset CDFs (line l occupies [2l, 2l+1])
    G1 = (C1 + 2 * np.arange(n)[:, None]).ravel()                     # increasing (non-strict inside a line)
    q = c + 2 * line
    pos = np.searchsorted(G1, q, side="right") - 1                    # index into flattened (n, n+1)
    pos = np.clip(pos, line * (n + 1), line * (n + 1) + n - 1)
    k = pos - line * (n + 1)
    lo, hi = G1[pos], G1[pos + 1]
    frac = np.where(hi > lo, (q - lo) / np.maximum(hi - lo, 1e-300), 0.5)
    xn = k + np.clip(frac, 0, 1 - 1e-12)
    dead = (tot0[line, 0] <= 1e-300) | (tot1[line, 0] <= 1e-300)
    xn = np.where(dead, x, xn)
    X = X.copy(); X[:, j] = xn
    return X


def det_tick(walk, psi, X):
    for j in (0, 1):
        rho0 = np.sum(np.abs(psi) ** 2, 0)
        psi = walk.axis_step(psi, j)
        rho1 = np.sum(np.abs(psi) ** 2, 0)
        X = cdf_move(X, j, rho0, rho1)
    return walk.coin(psi), X


def run(n, M, start, T, theta, cell, seed=SEED, omega0=0.9, every=5):
    wk = LatticeWalk(n, theta, periodic=False)
    psi, om = modes_state(wk, M, omega0=omega0)
    rng = np.random.default_rng(seed)
    rho = np.sum(np.abs(psi) ** 2, 0)
    xx, yy = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    prob = {"psi2": rho, "uniform": np.ones((n, n)), "half": (xx < n // 2).astype(float)}[start]
    X = sample_sites(rng, NS, prob).astype(float) + rng.random((NS, 2))
    ts, H = [], []
    for t in range(T + 1):
        if t % every == 0:
            rho = np.sum(np.abs(psi) ** 2, 0)
            h, ncell = hbar(np.floor(X).astype(int), rho, cell)
            ts.append(t); H.append(h)
        if t < T:
            psi, X = det_tick(wk, psi, X)
    F = (ncell - 1) / (2 * NS)
    H = np.array(H)
    return dict(t=np.array(ts), H=H, F=F, tau=fit_tau(ts, H, F), ratio=float(H[-1] / H[0]), H0=float(H[0]),
                HT=float(H[-1]), maxH=float(H.max()))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "main"
    n, th, cell, T = (24, 0.5, 4, 400) if mode == "main" else (96, 0.125, 16, 1600)
    t0 = time.time()
    res, curves = {}, {}
    for M, start in ((16, "psi2"), (16, "uniform"), (1, "uniform"), (16, "half"), (1, "half")):
        r = run(n, M, start, T, th, cell)
        key = f"M{M}_{start}"
        res[key] = {k: r[k] for k in ("F", "tau", "ratio", "H0", "HT", "maxH")}
        curves[key] = r["H"]; curves["t"] = r["t"]
        print(key, res[key], f"{time.time()-t0:.0f}s", flush=True)
    v = {}
    v["R3a_equivariance"] = res["M16_psi2"]["maxH"] <= 3 * res["M16_psi2"]["F"]
    v["R3b_relax_uniform"] = res["M16_uniform"]["ratio"] < 0.10
    v["R3b_relax_half"] = res["M16_half"]["ratio"] < 0.10
    for s in ("uniform", "half"):
        a, b = res[f"M1_{s}"], res[f"M16_{s}"]
        tr = a["tau"] / b["tau"] if b["tau"] > 0 else float("nan")
        res[f"tau_ratio_{s}"] = tr
        v[f"R3c_control_{s}"] = bool((tr >= 3) or (a["ratio"] > 0.10 and b["ratio"] < 0.03))
    v["PASS(primary: a, b, c_uniform)"] = bool(v["R3a_equivariance"] and v["R3b_relax_uniform"] and
                                               v["R3b_relax_half"] and v["R3c_control_uniform"])
    res["verdict"] = {k: bool(x) for k, x in v.items()}
    res["time_s"] = time.time() - t0
    print(json.dumps(res["verdict"], indent=1), res["time_s"])
    json.dump(res, open(os.path.join(OUT, f"redo_r3_{mode}.json"), "w"), indent=1, default=float)
    np.savez(os.path.join(OUT, f"redo_r3_{mode}.npz"), **curves)

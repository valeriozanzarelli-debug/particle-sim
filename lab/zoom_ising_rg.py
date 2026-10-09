"""ZOOM test b, part 2: learning the block-spin rule at every level, fixed point, exponents, checks.

Reads lab/out/zoom_ising_mc.npz written by lab/zoom_ising.py (Swendsen-Wang Monte Carlo, L = 256,
blocked with zoom_core.block at every sweep).  Exact: T_c = 2/ln(1+sqrt2), y_t = 1 (nu = 1),
y_h = 15/8 (eta = 1/4).
Run: timeout 590 python lab/zoom_ising_rg.py
"""
import os, sys, time, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zoom_core as zc
from zoom_ising import TC, KC, L0, NLEV, OUT, lattice_state, Ops  # noqa

T0 = time.time()
NAMES = {"below": "T = 0.96 Tc", "Tc": "T = Tc", "above": "T = 1.04 Tc"}


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


def mcrg(ops, n, jk=20):
    """Swendsen MCRG between levels n-1 -> n: T = A^-1 B, A = <S^n S^n>_c, B = <S^n S^(n-1)>_c.
    Returns y_t, err, y_h, err (jackknife, 20 blocks)."""
    ev, od = [0, 1, 2, 3], [4, 5]

    def y(sub, sel):
        X = sub[:, n, :][:, sel]; Y = sub[:, n - 1, :][:, sel]
        X = X - X.mean(0); Y = Y - Y.mean(0)
        w = np.linalg.eigvals(np.linalg.solve(X.T @ X, X.T @ Y))
        return np.log2(np.max(w.real))
    yt, yh = y(ops, ev), y(ops, od)
    nb = len(ops) // jk
    jt, jh = [], []
    for b in range(jk):
        sub = ops[np.r_[0:b * nb, (b + 1) * nb:len(ops)]]
        jt.append(y(sub, ev)); jh.append(y(sub, od))
    return yt, np.sqrt((jk - 1) * np.var(jt)), yh, np.sqrt((jk - 1) * np.var(jh))


def learn_level(cfgs, L, sites, rng):
    g = zc.Geometry(lattice_state(L, cfgs[0]))
    states = [lattice_state(L, c) for c in cfgs]
    m = zc.LocalModel()
    info = m.learn_stat(states, geom=g, alpha=1e-6, sites=sites, rng=rng)
    return m, info, g


def coupl(m):
    c = m.couplings()
    return dict(K1=float(m.phi(1.0)), K2=float(m.phi(np.sqrt(2))), K3=float(m.phi(2.0)),
                K5=float(m.phi(np.sqrt(5))), C4=float(c["C4"]), C3=float(c["C3"]), H=float(c["H"]))


def main():
    D = np.load(os.path.join(OUT, "zoom_ising_mc.npz"))
    rng = np.random.default_rng(7)
    res = dict(Tc=TC, flows={}, mcrg=[], checks=[])
    # ---- level-0 sanity: the learner on the fine lattice must find K1 = 1/T_c, everything else 0
    cf = D["Tc_small_store0"][-150:]
    m0, _, _ = learn_level(cf, 64, None, rng)
    c0 = coupl(m0)
    halves = [coupl(learn_level(h, 64, None, rng)[0])["K1"] for h in (cf[:75], cf[75:])]
    log(f"level-0 sanity (L=64 at T_c, 150 configs): K1={c0['K1']:.4f} (exact {KC:.4f}; two halves {halves[0]:.4f}, {halves[1]:.4f}) "
        f"K_sqrt2={c0['K2']:+.4f} K_2={c0['K3']:+.4f} C4={c0['C4']:+.4f} C3={c0['C3']:+.4f} H={c0['H']:+.4f} (exact 0)")
    res["level0"] = c0; res["level0_halves"] = halves
    # ---- flows of the learned couplings
    models = {}
    for name in ("Tc", "below", "above"):
        flow = []
        for n in range(1, NLEV + 1):
            Ln = L0 >> n
            cfg = D[f"{name}_store{n}"]
            cfg = cfg[-min(len(cfg), 120 if Ln > 64 else 200):]
            m, info, g = learn_level(cfg, Ln, 1500 if Ln > 64 else None, rng)
            c = coupl(m); c["level"] = n; c["L"] = Ln
            flow.append(c)
            models[(name, n)] = (m, g)
            log(f"  {name:5s} level {n} ({Ln:3d}^2): K1={c['K1']:+.4f} K_sqrt2={c['K2']:+.4f} K_2={c['K3']:+.4f} K_sqrt5={c['K5']:+.4f} C4={c['C4']:+.4f} C3={c['C3']:+.4f} H={c['H']:+.4f}")
        res["flows"][name] = flow

    def vec(c):
        return np.array([c["K1"], c["K2"], c["K3"], c["C4"]])
    for name in ("Tc", "below", "above"):
        f = res["flows"][name]
        d = [float(np.linalg.norm(vec(f[i + 1]) - vec(f[i]))) for i in range(len(f) - 1)]
        log(f"  |K^(n+1) - K^(n)| ({name}): " + ", ".join(f"{x:.4f}" for x in d))
        res["flows"][name + "_steps"] = d
    # ---- MCRG exponents at T_c
    ops = D["Tc_ops"]
    for n in range(1, NLEV + 1):
        yt, et, yh, eh = mcrg(ops, n)
        res["mcrg"].append(dict(level=n, yt=yt, et=et, yh=yh, eh=eh))
        log(f"  MCRG {n - 1}->{n}: y_t = {yt:.3f} +- {et:.3f} (exact 1; nu = {1 / yt:.3f})   y_h = {yh:.4f} +- {eh:.4f} (exact 1.875; eta = 2(2-y_h) = {2 * (2 - yh):.3f})")
    fb, fa = res["flows"]["below"], res["flows"]["above"]
    split = [fb[i]["K1"] - fa[i]["K1"] for i in range(NLEV)]
    yt_fl = [float(np.log2(split[i + 1] / split[i])) for i in range(NLEV - 1)]
    log(f"  learned-coupling flow: K1(0.96Tc)-K1(1.04Tc) per level {np.round(split, 4)} -> y_t from ratios {np.round(yt_fl, 3)} (crude: |t| = 4% is not infinitesimal)")
    res["yt_flow"] = yt_fl; res["split"] = split
    # ---- spot-checks: heat bath of the learned rule from blocked fine configurations
    for n in (1, 2, 3):
        Ln = L0 >> n
        m, g = models[("Tc", n)]
        store = D[f"Tc_store{n}"]
        st = lattice_state(Ln, store[0])
        hb = zc.HeatBath(m, st)
        mn = zc.LocalModel(); mn.mode = "stat"; mn.geom = g          # naive: the fine rule, unchanged
        th = np.zeros(g.nb + 3); th[int(round((1.0 - zc.DR) / zc.DR))] = KC; mn.theta = th
        hbn = zc.HeatBath(mn, st)
        opsn = Ops(st)
        starts = store[::max(1, len(store) // 12)][:12]
        ref = np.array([opsn(c.astype(float)) for c in store]) / Ln ** 2
        refm, refe = ref.mean(0), ref.std(0) / np.sqrt(len(ref) / 2)
        out = {"learned": [], "naive": []}
        nsw = 30 if n == 1 else 60
        t0 = time.time()
        for c in starts:
            for key, h in (("learned", hb), ("naive", hbn)):
                q = h.sweep(c.astype(float), rng, nsw)
                out[key].append(opsn(q) / Ln ** 2)
        tchk = time.time() - t0
        row = dict(level=n, sweeps=nsw, ms_per_sweep=1000 * tchk / (2 * nsw * len(starts)), ref=refm[:4].tolist(), ref_err=refe[:4].tolist())
        for key in out:
            v = np.array(out[key])
            err = np.abs(v.mean(0)[:4] - refm[:4]) / np.abs(refm[:4])
            row[key] = v.mean(0)[:4].tolist(); row[key + "_err"] = err.tolist()
            row[key + "_sd"] = (v.std(0)[:4] / np.sqrt(len(v)) / np.abs(refm[:4])).tolist()
            row[key + "_flag"] = bool(err[0] > 0.02)
        res["checks"].append(row)
        log(f"  check level {n} ({Ln}^2, {len(starts)} starts x {nsw} heat-bath sweeps, {row['ms_per_sweep']:.1f} ms/sweep): <s s>_nn per site: fine {refm[0]:.4f}+-{refe[0]:.4f}, "
            f"learned rule {row['learned'][0]:.4f} (err {100 * row['learned_err'][0]:.1f}% +- {100 * row['learned_sd'][0]:.1f}%{' FLAG' if row['learned_flag'] else ''}), "
            f"naive fine rule {row['naive'][0]:.4f} (err {100 * row['naive_err'][0]:.1f}%{' FLAG' if row['naive_flag'] else ''}); "
            f"4-ball: fine {refm[3]:.4f} learned {row['learned'][3]:.4f} naive {row['naive'][3]:.4f}")
    meta = D["Tc_meta"]
    res["cost"] = dict(sw_ms=float(D["sw_ms"][0]), mc_ms_with_blocking=float(1000 * meta[2] / meta[3]))
    log(f"  cost: bare SW sweep at L=256 {D['sw_ms'][0]:.1f} ms; heat bath of the learned rule: " +
        ", ".join(f"level {r['level']} {r['ms_per_sweep']:.1f} ms" for r in res["checks"]))
    json.dump(res, open(os.path.join(OUT, "zoom_ising.json"), "w"), indent=1)
    # ---- figures
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3, NLEV + 1, figsize=(1.6 * (NLEV + 1), 5.2), dpi=80)
        for row, name in enumerate(("below", "Tc", "above")):
            for n in range(NLEV + 1):
                Ln = L0 >> n
                a = ax[row, n]
                a.imshow(D[f"{name}_store{n}"][-1].reshape(Ln, Ln), cmap="gray", vmin=-1, vmax=1, interpolation="nearest")
                a.set_xticks([]); a.set_yticks([])
                if row == 0:
                    a.set_title(f"livello {n}\n{Ln}x{Ln}", fontsize=8)
                if n == 0:
                    a.set_ylabel(NAMES[name], fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "zoom_ising_snapshots.png")); plt.close(fig)
        fig, ax = plt.subplots(1, 1, figsize=(4.2, 3), dpi=90)
        tt = {"below": 0.96, "Tc": 1.0, "above": 1.04}
        for name, col in (("below", "C3"), ("Tc", "k"), ("above", "C0")):
            f = res["flows"][name]
            ax.plot([0] + [c["level"] for c in f], [1 / (tt[name] * TC)] + [c["K1"] for c in f], "o-", color=col, label=NAMES[name])
        ax.set_xlabel("livello di zoom n"); ax.set_ylabel("K1 imparato (primi vicini)"); ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "zoom_ising_flow.png")); plt.close(fig)
    except Exception as e:
        log("figure failed:", e)
    log("done")


if __name__ == "__main__":
    main()

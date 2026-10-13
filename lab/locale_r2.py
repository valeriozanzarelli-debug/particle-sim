"""locale_r2.py - REGOLA LOCALE, production run 2: the fermion exams D1-D8 (pre-registration lab/out/locale_prereg.json).
Writes lab/out/locale_r2.json.  Run: OMP_NUM_THREADS=1 timeout 590 python lab/locale_r2.py"""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
from locale_core import *
from locale_cands import fermion

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
t0 = time.time()
rng = np.random.default_rng(SEED + 2)
res, verdict = {}, []


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def E_low(w, K):
    """smallest |E| (per unit time) of the period operator at each k"""
    eps, _ = w.phases(np.atleast_2d(K))
    return np.min(np.abs(eps), axis=1) / w.T


def dirs(dim):
    d = [np.array([1, 0, 0.]), np.array([1, 1, 0]) / 2 ** .5]
    if dim == 3:
        d.append(np.array([1, 1, 1]) / 3 ** .5)
    return d


def d_exams(w, th, dim):
    out = {}
    out["E0_rel"] = float(abs(E_low(w, np.zeros((1, 3)))[0] / th - 1))
    errs = []
    for d in dirs(dim):
        P = np.linspace(0.1, 1, 10)[:, None] * th * d[None]
        E = E_low(w, P)
        errs.append(np.max(np.abs(E / np.hypot(np.linalg.norm(P, axis=1), th) - 1)))
    out["disp"] = float(max(errs))
    mk = []
    for d in dirs(dim):
        ps = np.linspace(0.02, 0.1, 9) * th
        E = E_low(w, ps[:, None] * d[None])
        E0 = E_low(w, np.zeros((1, 3)))[0]
        A = np.stack([ps ** 2, ps ** 4], 1)
        c = np.linalg.lstsq(A, E ** 2 - E0 ** 2, rcond=None)[0]
        mk.append(1.0 / c[0])                 # c[0] = m/m_kin  ->  m_kin/m = 1/c[0]
    out["mkin_rel"] = [float(m - 1) for m in mk]          # m_kin/m - 1 per direction
    # drift
    p = 1e-5
    g = []
    for j in range(dim):
        e = np.zeros(3); e[j] = p
        eps_p = np.sort(w.phases(e[None])[0][0]); eps_m = np.sort(w.phases(-e[None])[0][0])
        g.append(np.max(np.abs(eps_p - eps_m)) / w.T / (2 * p))
    out["drift"] = float(max(g))
    return out


def doublers(w, dim, kc=0.4, n=40):
    g = (np.arange(n) + 0.5) / n * PI - PI / 2
    grids = np.meshgrid(*([g] * dim), indexing="ij")
    K = np.zeros((n ** dim, 3))
    for a in range(dim):
        K[:, a] = grids[a].ravel()
    kn = np.linalg.norm(K, axis=1)
    K = K[kn >= kc]
    E = E_low(w, K)
    order = np.argsort(E)[:24]
    best = (np.inf, None)

    def f(x):
        k = np.zeros(3); k[:dim] = x
        k[:dim] = (k[:dim] + PI / 2) % PI - PI / 2
        pen = max(0.0, kc - np.linalg.norm(k)) * 100
        return float(E_low(w, k[None])[0] + pen)
    for i in order:
        r = minimize(f, K[i, :dim], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxiter=600))
        if r.fun < best[0]:
            best = (r.fun, r.x)
    kb = (np.array(best[1]) + PI / 2) % PI - PI / 2
    return float(best[0]), kb.tolist(), float(E[order[0]])


def pi_gap(w, dim, n=24):
    K = rng.uniform(-PI, PI, (n ** dim, 3)); K[:, dim:] = 0
    eps, _ = w.phases(K)
    return float(np.min(PI - np.abs(eps)))


res["D"] = {}
for name, dim in (("F0", 3), ("F2", 3), ("F2b", 3), ("F3", 3), ("F0", 2), ("F2", 2)):
    key = f"{name}_{dim}D"
    th = 0.1
    w = fermion(name, dim, th)
    r = d_exams(w, th, dim)
    check(f"D1 {key}", r["E0_rel"] <= 1e-12, f"|E(0)/m - 1| = {r['E0_rel']:.1e}")
    check(f"D2 {key}", r["disp"] <= 0.01, f"max |E/sqrt(p^2+m^2) - 1| (p <= m) = {r['disp']:.2e}")
    mx = max(abs(x) for x in r["mkin_rel"])
    check(f"D3 {key}", mx <= 0.01, f"m_kin/m - 1 = {', '.join(f'{x:+.3e}' for x in r['mkin_rel'])} (axis, face{', body' if dim == 3 else ''})")
    check(f"D4 {key}", r["drift"] <= 1e-9, f"|grad E(0)| = {r['drift']:.1e} c")
    if name != "F3":
        for thd in (0.1, 0.02):
            wd = fermion(name, dim, thd)
            Emin, kb, Egrid = doublers(wd, dim)
            r[f"doubler_theta{thd}"] = dict(Emin=Emin, k=kb, Egrid=Egrid)
            check(f"D5 {key} theta={thd}", Emin >= 0.3,
                  f"lowest |E| outside |k|<0.4 = {Emin:.4f}/a at k/pi = {np.round(np.array(kb) / PI, 3).tolist()} (m = {thd})")
        r["pi_gap"] = pi_gap(w, dim)
        print(f"     (info) {key}: smallest distance of a quasi-energy from pi (per period) = {r['pi_gap']:.4f} rad", flush=True)
    res["D"][key] = r
    print(f"     ({time.time() - t0:.0f} s)", flush=True)

# D6: independence of the Wilson strength
print("=== D6 Wilson strength h: low-energy numbers vs the frozen rule ===", flush=True)
ref = res["D"]["F0_3D"]
res["D6"] = {}
for name in ("F2", "F2b"):
    for h in (0.05, 0.1, 0.2):
        w = fermion(name, 3, 0.1, h=h)
        r = d_exams(w, 0.1, 3)
        dd = abs(r["disp"] - ref["disp"])
        dm = max(abs(a - b) for a, b in zip(r["mkin_rel"], ref["mkin_rel"]))
        res["D6"][f"{name}_h{h}"] = dict(disp=r["disp"], mkin_rel=r["mkin_rel"], d_disp=dd, d_mkin=dm)
        check(f"D6 {name} h={h}", max(dd, dm) <= 1e-6, f"|d disp| = {dd:.1e}, |d (m_kin/m)| = {dm:.1e} vs frozen")

# D8: order in theta
print("=== D8 kinetic-mass error vs theta ===", flush=True)
res["D8"] = {}
for name, lo, hi in (("F0", 3, 5), ("F3", 12, 20), ("F2", 3, 5)):
    e = {}
    for th in (0.2, 0.1, 0.05, 0.025):
        r = d_exams(fermion(name, 3, th), th, 3)
        e[th] = max(abs(x) for x in r["mkin_rel"])
    ratio = e[0.1] / e[0.05]
    res["D8"][name] = dict(err=e, ratio=ratio)
    check(f"D8 {name}", lo <= ratio <= hi, f"|m_kin/m - 1|: " + ", ".join(f"th={t}: {v:.2e}" for t, v in e.items()) + f"; ratio 0.1/0.05 = {ratio:.2f}")

# D7: lattice light cone
print("=== D7 largest group velocity over the zone ===", flush=True)
res["D7"] = {}
K = rng.uniform(-PI, PI, (4000, 3))
for name, dim in (("F0", 3), ("F2", 3), ("F2b", 3), ("F3", 3), ("F0", 2), ("F2", 2)):
    w = fermion(name, dim, 0.1)
    KK = K.copy(); KK[:, dim:] = 0
    v, gap = group_velocities(w, KK)
    ok = gap > 1e-3
    kk = np.repeat(np.linalg.norm(KK, axis=1)[:, None], v.shape[1], 1)
    sup = (v > 1 + 1e-6) & ok
    res["D7"][f"{name}_{dim}D"] = dict(vmax=float(v[ok].max()), frac_superluminal=float(sup.mean()),
                                        kmin_superluminal=float(kk[sup].min()) if sup.any() else None)
    check(f"D7 {name}_{dim}D", v[ok].max() <= 1 + 1e-6,
          f"max v_g = {v[ok].max():.4f} c; fraction of (k, band) above c = {sup.mean():.2e}"
          + (f"; smallest |k| with v > c: {kk[sup].min():.3f}" if sup.any() else ""))
    print(f"     ({time.time() - t0:.0f} s)", flush=True)

res["verdict"] = verdict
res["seconds"] = time.time() - t0
with open(os.path.join(OUT, "locale_r2.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
npass = sum(v for _, v in verdict)
print(f"=== {npass} PASS / {len(verdict) - npass} FAIL ({time.time() - t0:.0f} s)", flush=True)

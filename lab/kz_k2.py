"""kz_k2.py - K2: equilibrium 3D XY at beta_c, does the simulation know its universality class (nu, Tc)?
usage: python kz_k2.py L nmeas      -> lab/out/kz_k2_L{L}.npz   (time series of E and |m|)
       python kz_k2.py fit          -> lab/out/kz_k2.json
Swendsen-Wang (embedded Ising) + 1 over-relaxation sweep per step, numpy (single thread).
"""
import sys, os, json, time, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from kz_core import sw_update, overrelax_np, TC3, NU_XY

BC = 1.0 / TC3
OUT = os.path.join(HERE, "out")


def run(L, nmeas, ntherm=1000):
    rng = np.random.default_rng(12345 + L)
    th = rng.uniform(-np.pi, np.pi, (L,) * 3)
    par = np.indices(th.shape).sum(0) % 2
    N = th.size
    E = np.empty(nmeas); M = np.empty(nmeas)
    t0 = time.time()
    for i in range(ntherm + nmeas):
        th = sw_update(th, BC, rng)
        th = overrelax_np(th, par)
        if i >= ntherm:
            k = i - ntherm
            E[k] = -sum(np.cos(np.roll(th, -1, d) - th).sum() for d in range(3))
            M[k] = math.hypot(np.cos(th).sum(), np.sin(th).sum()) / N
    np.savez(os.path.join(OUT, "kz_k2_L%d.npz" % L), E=E, M=M, L=L, beta=BC)
    print(L, nmeas, "%.0fs" % (time.time() - t0), "U4=%.4f" % (np.mean(M ** 4) / np.mean(M ** 2) ** 2), flush=True)


def obs(E, M, beta=BC, beta0=BC):
    """reweighted observables at beta from data at beta0."""
    x = -(beta - beta0) * (E - E.mean()); w = np.exp(x - x.max()); w /= w.sum()
    a = lambda v: (w * v).sum()
    m1, m2, m4, e = a(M), a(M ** 2), a(M ** 4), a(E)
    U4 = m4 / m2 ** 2
    dlnm = (m1 * e - a(M * E)) / m1          # d ln<|m|>/d beta
    dlnm2 = (m2 * e - a(M ** 2 * E)) / m2    # d ln<m^2>/d beta
    return U4, dlnm, dlnm2


def jack(E, M, f, nb=20):
    n = len(E) // nb
    full = np.array(f(E, M))
    js = np.array([f(np.delete(E, slice(i * n, (i + 1) * n)), np.delete(M, slice(i * n, (i + 1) * n))) for i in range(nb)])
    err = np.sqrt((nb - 1) / nb * ((js - js.mean(0)) ** 2).sum(0))
    return full, err, js


def fit():
    Ls = [8, 12, 16, 24, 32]
    data = {L: np.load(os.path.join(OUT, "kz_k2_L%d.npz" % L)) for L in Ls}
    res = {"L": Ls}
    vals = {}; jks = {}
    for L in Ls:
        E, M = data[L]["E"], data[L]["M"]
        v, e, js = jack(E, M, lambda E, M: obs(E, M))
        vals[L] = (v, e); jks[L] = js
        res["L%d" % L] = {"n": len(E), "U4": [v[0], e[0]], "dlnm_dbeta": [v[1], e[1]], "dlnm2_dbeta": [v[2], e[2]]}
    # power-law fits: ln |deriv| = c + (1/nu) ln L
    def nu_fit(idx, Lsel, jk=None):
        y = np.array([np.log(abs(vals[L][0][idx] if jk is None else jks[L][jk][idx])) for L in Lsel])
        p = np.polyfit(np.log(Lsel), y, 1)
        return 1.0 / p[0]
    nb = 20
    est = {}
    for name, idx in (("dlnm", 1), ("dlnm2", 2)):
        for win in ((8, 12, 16, 24, 32), (12, 16, 24, 32)):
            nu = nu_fit(idx, win)
            js = np.array([nu_fit(idx, win, j) for j in range(nb)])
            err = math.sqrt((nb - 1) / nb * ((js - js.mean()) ** 2).sum())
            est["%s_L%d-32" % (name, win[0])] = [nu, err]
    main = 0.5 * (est["dlnm_L8-32"][0] + est["dlnm2_L8-32"][0])
    stat = 0.5 * math.hypot(est["dlnm_L8-32"][1], est["dlnm2_L8-32"][1])
    sys_est = 0.5 * abs(est["dlnm_L8-32"][0] - est["dlnm2_L8-32"][0])
    sys_win = 0.5 * abs(0.5 * (est["dlnm_L12-32"][0] + est["dlnm2_L12-32"][0]) - main)
    err = math.sqrt(stat ** 2 + sys_est ** 2 + sys_win ** 2)
    res["nu_estimates"] = est
    res["nu"] = [main, err]
    # Binder crossings by reweighting
    def U_at(L, T):
        return obs(data[L]["E"], data[L]["M"], beta=1.0 / T)[0]
    cross = {}
    for (L1, L2) in ((8, 16), (12, 24), (16, 32)):
        Ts = np.linspace(2.18, 2.225, 451)
        d = np.array([U_at(L1, T) - U_at(L2, T) for T in Ts])
        s = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
        if len(s):
            i = s[0]; Tx = Ts[i] - d[i] * (Ts[i + 1] - Ts[i]) / (d[i + 1] - d[i])
            cross["%d-%d" % (L1, L2)] = [Tx, U_at(L2, Tx)]
        else:
            cross["%d-%d" % (L1, L2)] = None
    res["binder_crossings_T_U4"] = cross
    Tx = cross.get("16-32")
    res["PASS"] = bool(abs(main - NU_XY) <= 0.05 and Tx is not None and abs(Tx[0] / TC3 - 1) <= 0.005)
    res["criteria"] = "|nu - 0.6717| <= 0.05 and T_cross(16,32) within 0.5% of 2.20184"
    with open(os.path.join(OUT, "kz_k2.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "fit":
        fit()
    else:
        run(int(sys.argv[1]), int(sys.argv[2]))

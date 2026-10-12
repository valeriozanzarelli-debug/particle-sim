"""glue_analysis.py - correlators, variational (GEVP) analysis, fixed plateau rule, jackknife.

usage: python glue_analysis.py [prefix]      (prefix '' = production chunks b{beta}_c*.npz, 'pilot_' = pilot)
Writes lab/out/glue_{prefix}results.json. The rules used here are the ones fixed in lab/out/glue_prereg.json.
"""
import os, sys, json, glob
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.linalg import eigh

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out"); RAW = os.path.join(OUT, "glue_raw")
RULE = dict(n_bins=40, prune=1e-3, t0=0, td=1, sig_min=3.0, rel_err_max=0.20)


def load(tag):
    # main chain first, then a second independent chain '{tag}x' if present (declared deviation for beta 12)
    files = sorted(glob.glob(os.path.join(RAW, f"{tag}_c*.npz"))) + sorted(glob.glob(os.path.join(RAW, f"{tag}x_c*.npz")))
    d = [np.load(f) for f in files]
    names = [str(x) for x in d[0]["names"]]; pnames = [str(x) for x in d[0]["pnames"]]
    ops = np.concatenate([x["ops"] for x in d]); pops = np.concatenate([x["pops"] for x in d])
    W = np.concatenate([x["W"] for x in d]); plaq = np.concatenate([x["plaq"] for x in d])
    act = np.concatenate([x["act"] for x in d])
    return names, pnames, ops, pops, W, plaq, act


def binned_corr(O, nb, complex_=False):
    """O: (N, nops, T) per-config operator time series. Returns per-bin sums S (nb, nops, nops, T),
    per-bin operator means M (nb, nops), counts n (nb). C_ij(t) = <O_i(t0+t) O_j(t0)>_t0."""
    N, K, T = O.shape
    edges = np.linspace(0, N, nb + 1).astype(int)
    S = np.zeros((nb, K, K, T)); M = np.zeros((nb, K)); n = np.zeros(nb)
    for b in range(nb):
        X = O[edges[b]:edges[b + 1]]
        F = np.fft.fft(X, axis=-1)                                  # (n, K, T)
        # sum_t0 O_i(t0+t) O_j(t0) = IFFT( F_i * conj(F_j) )
        P = np.einsum("nkw,nlw->klw", F, F.conj())
        c = np.fft.ifft(P, axis=-1) / T
        S[b] = c.real
        M[b] = X.mean(axis=-1).real.sum(axis=0) if not complex_ else 0.0
        n[b] = len(X)
    return S, M, n


def jk_corrs(S, M, n, subtract=True):
    """jackknife correlator samples (nb+1, K, K, T): index 0 = full sample."""
    nb = len(n)
    out = []
    for b in [-1] + list(range(nb)):
        sel = np.ones(nb, bool)
        if b >= 0:
            sel[b] = False
        N = n[sel].sum()
        C = S[sel].sum(0) / N
        if subtract:
            m = M[sel].sum(0) / N
            C = C - np.outer(m, m)[:, :, None]
        C = 0.5 * (C + np.roll(C[:, :, ::-1], 1, axis=-1))         # C(t) <- (C(t)+C(T-t))/2
        C = 0.5 * (C + C.transpose(1, 0, 2))
        out.append(C)
    return np.array(out)


def jk_err(samples):
    """samples[0] = full, samples[1:] = jackknife. returns (value, error)."""
    x = np.asarray(samples)
    nb = len(x) - 1
    jk = x[1:]
    return x[0], np.sqrt((nb - 1) / nb * np.sum((jk - jk.mean(0)) ** 2, axis=0))


def gevp(Cfull):
    """normalise, prune (keep eigen-directions of C(0) with lambda >= prune*lambda_max),
    solve C(td) v = lambda C(t0) v; return projection vectors (K, nstates) ordered by decreasing lambda."""
    d = 1.0 / np.sqrt(np.diag(Cfull[:, :, 0]))
    Cn = Cfull * d[:, None, None] * d[None, :, None]
    w, V = eigh(Cn[:, :, RULE["t0"]])
    keep = w >= RULE["prune"] * w.max()
    R = V[:, keep] / np.sqrt(w[keep])                            # whitening of C(t0)
    A = R.T @ Cn[:, :, RULE["td"]] @ R
    lam, X = eigh(0.5 * (A + A.T))
    order = np.argsort(lam)[::-1]
    vec = (R @ X[:, order]) * d[:, None]
    return vec, lam[order], int(keep.sum())


def meff_cosh(c, t, T):
    if not (np.isfinite(c[t]) and np.isfinite(c[t + 1])) or c[t] <= 0 or c[t + 1] <= 0 or c[t] <= c[t + 1]:
        return np.nan
    r = c[t] / c[t + 1]
    f = lambda m: np.cosh(m * (t - T / 2)) / np.cosh(m * (t + 1 - T / 2)) - r
    try:
        return brentq(f, 1e-6, 20.0)
    except ValueError:
        return np.nan


def analyse_state(Cs, T):
    """Cs: (nb+1, T) jackknife samples of one projected correlator. Applies the fixed plateau rule."""
    th = T // 2
    c0, dc = jk_err(Cs[:, :th + 1])
    snr = c0 / np.where(dc > 0, dc, np.inf)
    t_sig = 0
    while t_sig + 1 <= th and snr[t_sig + 1] >= RULE["sig_min"]:
        t_sig += 1
    me = np.array([[meff_cosh(Cs[k], t, T) for t in range(th)] for k in range(len(Cs))])
    m0, dm = jk_err(me)
    t_min = None
    for t in range(1, t_sig):
        ok = (np.isfinite(dm[t]) and np.isfinite(dm[t + 1]) and np.isfinite(m0[t]) and np.isfinite(m0[t + 1])
              and dm[t] <= RULE["rel_err_max"] * m0[t] and dm[t + 1] <= RULE["rel_err_max"] * m0[t + 1])
        if ok and abs(m0[t] - m0[t + 1]) <= np.hypot(dm[t], dm[t + 1]):
            t_min = t; break
    flag = "plateau"
    if t_min is None:
        t_min = 1; flag = "no_plateau_fallback_t1"
    t_max = max(t_sig, t_min + 1)
    ts = np.arange(t_min, t_max + 1)
    sig = np.where(dc[ts] > 0, dc[ts], 1.0)
    f = lambda t, A, m: A * np.cosh(m * (t - T / 2))
    res = []
    for k in range(len(Cs)):
        m_init = m0[t_min] if np.isfinite(m0[t_min]) else 1.0
        try:
            p, _ = curve_fit(f, ts, Cs[k, ts], p0=[Cs[k, t_min] / np.cosh(m_init * (t_min - T / 2)), m_init],
                             sigma=sig, maxfev=20000)
            res.append(p[1])
        except Exception:
            res.append(np.nan)
    m, em = jk_err(np.array(res))
    return dict(m=float(m), err=float(em), t_min=int(t_min), t_max=int(t_max), flag=flag, t_sig=int(t_sig),
                meff=[None if not np.isfinite(x) else float(x) for x in m0],
                meff_err=[None if not np.isfinite(x) else float(x) for x in dm],
                c_over_c0=[float(x) for x in c0 / c0[0]])


def channel(Cjk, idx, T, nstates):
    C = Cjk[:, idx][:, :, idx]
    vec, lam, nkeep = gevp(C[0])
    out = []
    for s in range(nstates):
        v = vec[:, s]
        Cs = np.einsum("i,kijt,j->kt", v, C, v)
        r = analyse_state(Cs, T)
        r["overlap_t0"] = float(lam[s]); r["n_ops"] = len(idx); r["n_kept"] = nkeep
        out.append(r)
    return out


def sigma_from_torelon(mP, L):
    """Nambu-Goto ground state of a flux tube of length L (D = 2+1): m = sigma L sqrt(1 - pi/(3 sigma L^2))."""
    return (np.pi / 3 + np.sqrt(np.pi ** 2 / 9 + 4 * L ** 2 * mP ** 2)) / (2 * L ** 2)


def creutz(W, nb=20):
    """Creutz ratios chi(R,R) = -ln[W(R,R) W(R-1,R-1) / W(R,R-1)^2], jackknife over bins of measurements."""
    Ws = 0.5 * (W + W.transpose(0, 2, 1))
    N = len(Ws); edges = np.linspace(0, N, nb + 1).astype(int)
    bins = np.array([Ws[edges[b]:edges[b + 1]].mean(0) for b in range(nb)])
    samp = [bins.mean(0)] + [np.delete(bins, b, 0).mean(0) for b in range(nb)]
    R = Ws.shape[1] - 1
    res = {}
    for r in range(2, R + 1):
        vals = []
        for Wm in samp:
            x = Wm[r, r] * Wm[r - 1, r - 1] / Wm[r, r - 1] ** 2
            vals.append(-np.log(x) if x > 0 else np.nan)
        v, e = jk_err(np.array(vals))
        res[r] = (float(v), float(e))
    return res


def analyse_beta(tag, plan):
    names, pnames, ops, pops, W, plaq, act = load(tag)
    T, L = plan["T"], plan["L"]
    nb = RULE["n_bins"]
    S, M, n = binned_corr(ops, nb)
    Cjk = jk_corrs(S, M, n)
    A1 = [i for i, q in enumerate(names) if q.startswith("A1_")]
    E = [i for i, q in enumerate(names) if q.startswith("E_")]
    P = [i for i, q in enumerate(names) if q.startswith("P_")]
    res = dict(n_meas=int(len(ops)), plaquette=float(plaq.mean()), T=T, L=L)
    res["A1"] = channel(Cjk, A1, T, 2)
    res["E"] = channel(Cjk, E, T, 1)
    res["P"] = channel(Cjk, P, T, 1)
    # momentum p = 2 pi / L (A1 operators), both directions
    Sp = 0
    for d in range(2):
        Sd, _, n_ = binned_corr(pops[:, :, :, d].astype(np.complex128), nb, complex_=True)
        Sp = Sp + Sd / 2
    Cp = jk_corrs(Sp, np.zeros((nb, len(pnames))), n_, subtract=False)
    res["A1p"] = channel(Cp, list(range(len(pnames))), T, 1)
    # string tension from the torelon (jackknife through the fit is approximated by error propagation)
    mP, emP = res["P"][0]["m"], res["P"][0]["err"]
    s = sigma_from_torelon(mP, L); ds = abs(sigma_from_torelon(mP + emP, L) - s)
    res["sigma_a2"] = float(s); res["sigma_a2_err"] = float(ds)
    res["a_sqrt_sigma"] = float(np.sqrt(s)); res["a_sqrt_sigma_err"] = float(ds / (2 * np.sqrt(s)))
    res["creutz"] = {str(k): v for k, v in creutz(W).items()} if len(W) >= 20 else {}
    res["n_creutz_configs"] = int(len(W))
    # energy check data: mean action per slice (vacuum), for the record
    res["action_per_site"] = [float(x) for x in act.mean(axis=(0, 1)) / (L * L)]
    return res


def main():
    prefix = sys.argv[1] if len(sys.argv) > 1 else ""
    planfile = "glue_pilot_plan.json" if prefix == "pilot_" else "glue_prereg.json"
    plans = json.load(open(os.path.join(OUT, planfile)))["runs"]
    allres = {}
    for b, plan in plans.items():
        tag = f"{prefix}b{b}"
        if not glob.glob(os.path.join(RAW, f"{tag}_c*.npz")):
            continue
        r = analyse_beta(tag, plan)
        allres[b] = r
        a0, a1, e0 = r["A1"][0], r["A1"][1], r["E"][0]
        print(f"beta {b}: N={r['n_meas']} a sqrt(sigma)={r['a_sqrt_sigma']:.4f}({r['a_sqrt_sigma_err']:.4f}) "
              f"am0={a0['m']:.4f}({a0['err']:.4f}) [{a0['t_min']},{a0['t_max']}] "
              f"am0*={a1['m']:.4f}({a1['err']:.4f}) [{a1['t_min']},{a1['t_max']}] "
              f"am2={e0['m']:.4f}({e0['err']:.4f}) [{e0['t_min']},{e0['t_max']}] "
              f"amP={r['P'][0]['m']:.4f}({r['P'][0]['err']:.4f}) E(p)={r['A1p'][0]['m']:.4f}({r['A1p'][0]['err']:.4f})")
    json.dump(dict(rule=RULE, results=allres), open(os.path.join(OUT, f"glue_{prefix}results.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

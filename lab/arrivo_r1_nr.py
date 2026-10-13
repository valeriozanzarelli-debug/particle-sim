"""arrivo_r1_nr.py CFG - NR reference runs (exam), configs NR_A, NR_B, NR_C of lab/out/arrivo_prereg.json.
Writes lab/out/arrivo_r1_<CFG>.json (summary) and the samples + arrival times to the scratch npz (not in git)."""
import json
import os
import sys
import time
import numpy as np
import arrivo_nr as A

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCR = os.environ.get("ARRIVO_SCRATCH", "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad")
CFG = dict(
    NR_A=dict(L=100.0, omega=1000.0, alphas=[np.pi / 2, 0.0], N=20000, seed=20261013, tau_end=60.0),
    NR_B=dict(L=10.0, omega=1000.0, alphas=[np.pi / 2, 0.0], N=4000, seed=20261014, tau_end=40.0),
    NR_C=dict(L=10.0, omega=20.0, alphas=[np.pi / 2, 0.0, np.pi / 4, 3 * np.pi / 4, np.pi], N=2000, seed=20261015,
              tau_end=40.0),
)


def stats(T):
    f = np.isfinite(T)
    Tf = T[f]
    return dict(arrived=float(f.mean()), mean=float(Tf.mean()), q50=float(np.median(Tf)),
                q90=float(np.percentile(Tf, 90)), q99=float(np.percentile(Tf, 99)), max=float(Tf.max()))


def flux_cdf(L, tau):
    tt = np.linspace(1e-4, tau.max(), 200001)
    P = A.flux(L, tt)
    C = np.concatenate([[0], np.cumsum(0.5 * (P[1:] + P[:-1]) * np.diff(tt))])
    return np.interp(tau, tt, C), tt, P


def ks_vs_cdf(T, cdf_fn, tau_max):
    Ts = np.sort(T[np.isfinite(T) & (T <= tau_max)])
    n = np.isfinite(T).size
    emp_hi = np.arange(1, Ts.size + 1) / n
    emp_lo = np.arange(0, Ts.size) / n
    c = cdf_fn(Ts)
    return float(max(np.max(np.abs(emp_hi - c)), np.max(np.abs(emp_lo - c))))


if __name__ == "__main__":
    name = sys.argv[1]
    c = CFG[name]
    t0 = time.time()
    rng = np.random.default_rng(c["seed"])
    X0, Y0, Z0 = A.sample_initial(c["N"], c["omega"], rng)
    res = dict(config=name, L=c["L"], omega=c["omega"], N=c["N"], seed=c["seed"], runs={})
    npz = {"X0": X0, "Y0": Y0, "Z0": Z0}
    dtf = (lambda t: 2e-4 if t < 0.05 else (2e-3 if t < 2 else 1e-2))
    if c["omega"] < 100:
        dtf = (lambda t: 2e-4 if t < 0.05 else (2e-3 if t < 2 else 5e-3))
    for al in c["alphas"]:
        T = A.bohm_arrivals(c["N"], c["L"], c["omega"], al, c["tau_end"], None, dt_fn=dtf, XYZ0=(X0, Y0, Z0))
        key = "alpha=%.4f" % al
        res["runs"][key] = stats(T)
        npz["T_" + key] = T
        print(name, key, res["runs"][key], "%.0f s" % (time.time() - t0), flush=True)
    # dt/2 check on the first 1000 samples, alpha = pi/2
    T1 = A.bohm_arrivals(1000, c["L"], c["omega"], np.pi / 2, c["tau_end"], None, dt_fn=lambda t: 0.5 * dtf(t),
                         XYZ0=(X0[:1000], Y0[:1000], Z0[:1000]))
    Tref = npz["T_alpha=%.4f" % (np.pi / 2)][:1000]
    res["dt_half_check"] = dict(stats_half=stats(T1), stats_ref=stats(Tref),
                                median_abs_diff=float(np.nanmedian(np.abs(T1 - Tref))))
    # flux CDF and KS for spin along the guide
    Tz = npz["T_alpha=0.0000"]
    tau_g = np.linspace(0, c["tau_end"], 4001)
    C, tt, P = flux_cdf(c["L"], tau_g)
    res["flux_cdf_at_tau_end"] = float(C[-1])
    res["KS_spin_up_vs_flux"] = ks_vs_cdf(Tz, lambda x: np.interp(x, tau_g, C), c["tau_end"])
    tmx = res["runs"]["alpha=%.4f" % (np.pi / 2)]["max"]
    res["P_up_after_taumax_x"] = float(np.mean(Tz[np.isfinite(Tz)] > tmx) * np.isfinite(Tz).mean()
                                       + (1 - np.isfinite(Tz).mean()))
    res["flux_mean_trunc"] = float(np.trapezoid(tt * P, tt) / np.trapezoid(P, tt))
    res["secs"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "arrivo_r1_%s.json" % name), "w"), indent=1)
    np.savez(os.path.join(SCR, "arrivo_r1_%s.npz" % name), **npz)
    print(json.dumps({k: v for k, v in res.items() if k != "runs"}, indent=1))

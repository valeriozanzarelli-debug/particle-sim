"""emiss6_run.py - production runs of pre-registration v6 (lab/out/emiss6_prereg.json).  Each job writes
lab/out/emiss6_res_<job>.json (small numbers; series stay in lab/out/emiss_ck).  Exit 3 = checkpoint, run again.
Usage: python lab/emiss6_run.py JOB   (jobs: q60_<a>, q1_<a>, rabi60_<a>, scan60, pk60_<a>, pk1_<a>)"""
import sys, os, json
import numpy as np
from emiss6_core import *

job = sys.argv[1]
budget = float(sys.argv[2]) if len(sys.argv) > 2 else 540


def save(obj):
    json.dump(obj, open(os.path.join(OUT, f"emiss6_res_{job}.json"), "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in obj.items() if not isinstance(v, (list, dict))}, default=float))


def slope_fit(t, y):
    return np.polyfit(t, y, 1)


kind, sa = job.split("_") if "_" in job else (job, None)
a = float(sa) if sa else None
if kind == "q60":
    Tover = 10.0 if a == 0.025 else 6.0
    res, meta = qe_run(job, 60, a, 40.0, Tover, budget=budget)
    if res is None:
        sys.exit(3)
    t, n, D = res["t"], res["norm"], res["D"]
    ArW = meta["A_rule"] * meta["W2"]
    i = t >= 0.5 / ArW
    A_fit = -slope_fit(t[i], np.log(n[i]))[0]
    dev = float(np.max(np.abs(n - np.exp(-A_fit * t))))
    j = (t >= 0.5 / ArW) & (t <= 3 / ArW)
    ph = np.unwrap(np.angle(D[j]))
    w_emit = -slope_fit(t[j], ph)[0]
    win = []
    for k in range(6):
        jj = (t >= k / ArW) & (t < (k + 1) / ArW) if k else (t >= 100) & (t < 1 / ArW)
        if np.sum(jj) > 10:
            win.append(float(-slope_fit(t[jj], np.log(n[jj]))[0] / ArW))
    out = dict(meta, A_fit=A_fit, A_fit_over_rule=A_fit / ArW, A_over_W2=A_fit / meta["W2"], max_dev_exp=dev,
               max_p2p_minus_norm=float(np.max(np.abs(res["p2p"] - n))), omega_emit=w_emit,
               omega_emit_over_rule=w_emit / meta["omega"], windows_over_rule=win, norm_end=float(n[-1]))
    if a == 0.025:
        J = res["J"] @ np.conj(VEC)
        nfft = 8 * len(J)
        S = np.abs(np.fft.fft(J, nfft))**2
        fr = -2 * np.pi * np.fft.fftfreq(nfft, d=meta["b"])     # J ~ e^{-i w t}
        o = np.argsort(fr); fr, S = fr[o], S[o]
        from scipy.optimize import curve_fit
        k0 = np.argmax(S)
        lor = lambda x, h, c, g: h / (1 + ((x - c) / (g / 2))**2)
        sel = np.abs(fr - fr[k0]) < 10 * A_fit
        p, _ = curve_fit(lor, fr[sel], S[sel], p0=[S[k0], fr[k0], A_fit])
        sel2 = np.abs(fr - p[1]) < 5 * abs(p[2])
        out.update(fwhm=abs(p[2]), fwhm_over_A=abs(p[2]) / A_fit, line_centre=p[1], line_centre_over_rule=p[1] / meta["omega"],
                   lor_resid=float(np.max(np.abs(lor(fr[sel2], *p) - S[sel2])) / p[0]))
    save(out)
elif kind == "q1":
    lev = levels(1, a, 2500.0)
    P = 2 * np.pi / lev["w"]
    res, meta = qe_run(job, 1, a, 2500.0, None, T_abs=6 * P, budget=budget)
    if res is None:
        sys.exit(3)
    t, n, D = res["t"], res["norm"], res["D"]
    j = t >= 2 * P
    A_fit = -slope_fit(t[j], np.log(n[j]))[0]
    w_emit = -slope_fit(t[j], np.unwrap(np.angle(D[j])))[0]
    ArW = meta["A_rule"] * meta["W2"]
    win = [float(-slope_fit(t[(t >= k * P) & (t < (k + 1) * P)], np.log(n[(t >= k * P) & (t < (k + 1) * P)]))[0] / ArW) for k in range(6)]
    save(dict(meta, A_fit=A_fit, A_fit_over_rule=A_fit / ArW, A_over_W2=A_fit / meta["W2"], omega_emit=w_emit,
              omega_emit_over_rule=w_emit / meta["omega"], windows_over_rule=win, loss=float(1 - n[-1])))
elif kind == "rabi60":
    res, meta = rabi_run(job, 60, a, 40.0, omega_frac=0.005, budget=budget)
    if res is None:
        sys.exit(3)
    Pm, W, r = fit_rabi(res["t"], res["pe"])
    Om = W * np.sqrt(Pm)
    Rd = Om * np.sqrt(6) / (2 * meta["eE0"])
    save(dict(meta, Pmax=Pm, W=W, Omega_meas=Om, Omega_over_rule=Om / meta["Om_rule"], fit_resid=r, R_dyn=Rd,
              f32=2 / 3 * meta["omega"] * 4 / 6 * Rd**2))
elif kind == "scan60":
    out = {}
    for d in (-2.0, -1.0, -0.5, 0.5, 1.0, 2.0):
        res, meta = rabi_run(f"{job}_{d:+g}", 60, 0.05, 40.0, omega_frac=0.005, detune=d, budget=budget)
        if res is None:
            sys.exit(3)
        out[f"{d:+g}"] = dict(Pmax=float(res["pe"].max()), theory=1 / (1 + d * d), diff=float(res["pe"].max() - 1 / (1 + d * d)))
    save(dict(scan=out, max_abs_diff=max(abs(v["diff"]) for v in out.values())))
elif kind in ("pk60", "pk1"):
    Z = 60 if kind == "pk60" else 1
    r = packet_absorb(Z, a, 40.0 if Z == 60 else 2500.0)
    save(dict(r, R_over_rule=r["R_dyn"] / r["R_rule"]))
else:
    raise SystemExit("unknown job " + job)

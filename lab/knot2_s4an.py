"""knot2_s4an.py - S4 spectral exam (pre-registered): omega_est(k) of the outgoing waves in the shell 4.5-6.3,
averaged over t = 16, 20, 24; fit omega^2 = c^2 k^2 + m^2 over the bins with 90 % of the spectral energy.
Adds 'spectral' to lab/out/knot2_s4_v<v>.json."""
import os, sys, json, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
for v in sys.argv[1:]:
    fn = os.path.join(OUT, f"knot2_s4_v{v}.json"); d = json.load(open(fn))
    kb = np.array(d["kbins"]); kc = np.concatenate([[0], 0.5 * (kb[1:] + kb[:-1]), [kb[-1]]])   # bin index i -> [kb[i-1], kb[i])
    ts = [t for t in ("16.0", "20.0", "24.0") if t in d["specs"]]
    P = sum(np.array(d["specs"][t]["P"]) for t in ts); Pd = sum(np.array(d["specs"][t]["Pd"]) for t in ts)
    WL = sum(np.array(d["specs"][t]["Wlat_w"]) for t in ts)
    ok = (P > 0) & (Pd > 0)
    w_est = np.where(ok, np.sqrt(Pd / np.where(P > 0, P, 1)), 0); w_lat = np.where(Pd > 0, WL / np.where(Pd > 0, Pd, 1), 0)
    order = np.argsort(-Pd); cum = np.cumsum(Pd[order]) / Pd.sum()
    sel = np.zeros(len(Pd), bool); sel[order[: int(np.searchsorted(cum, 0.9)) + 1]] = True
    sel &= ok & (kc > 0)
    A = np.c_[kc[sel] ** 2, np.ones(sel.sum())]; W = Pd[sel]
    coef = np.linalg.lstsq(A * np.sqrt(W)[:, None], w_est[sel] ** 2 * np.sqrt(W), rcond=None)[0]
    c = math.sqrt(max(coef[0], 0)); m2 = coef[1]
    ratio = float((Pd[sel] * w_est[sel] / w_lat[sel]).sum() / Pd[sel].sum())
    kmean = float((Pd * kc).sum() / Pd.sum())
    sp = dict(times=ts, c_fit=c, m2_fit=m2, ratio_to_lattice=ratio, k_mean=kmean, wavelength=2 * math.pi / kmean,
              k_bins_used=kc[sel].tolist(), omega_est=w_est.tolist(), omega_lat=w_lat.tolist(), k_center=kc.tolist(), Pd=Pd.tolist(),
              pass_=bool(abs(c - 1) < 0.1 and m2 < 0.25**2 and abs(ratio - 1) < 0.1))
    d["spectral"] = sp
    # POST-HOC (not pre-registered, added after the spectral exam failed): speed of the outgoing energy front =
    # slope of r vs the time of the peak of the shell energy e_r(t), shells 3.0 <= r <= 6.4 (inside the sponge-free zone)
    rb = np.array(d["rbins"]); E = np.array([r["e_r"] for r in d["records"]]); t = np.array([r["t"] for r in d["records"]])
    js = [j for j in range(1, len(rb)) if 3.0 <= rb[j - 1] and rb[j] <= 6.4 + 1e-9]
    rr = np.array([0.5 * (rb[j - 1] + rb[j]) for j in js]); tp = np.array([t[np.argmax(E[:, j] * (t > 3))] for j in js])
    cf = np.polyfit(tp, rr, 1)
    d["front_posthoc"] = dict(r=rr.tolist(), t_peak=tp.tolist(), speed=float(cf[0]),
                              note="post-hoc: energy-front speed from peak times, not pre-registered")
    json.dump(d, open(fn, "w"))
    print(v, {k: val for k, val in sp.items() if not isinstance(val, list)}, "front speed", d["front_posthoc"]["speed"])

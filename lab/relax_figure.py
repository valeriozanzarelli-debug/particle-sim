"""relax_figure.py - scores the pre-registered criteria (relax_prereg.json) and draws lab/out/relax_figure.png.
Reads lab/out/relax_box_*.npz and relax_slit_{1d,2d}.json; writes lab/out/relax_results.json.
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from relax_core import Box, PI

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def load(name):
    p = os.path.join(OUT, f"relax_box_{name}.npz")
    return dict(np.load(p)) if os.path.exists(p) else None


def tau(t, H, F):
    h = H - F
    target = h[0] / np.e
    i = np.nonzero(h <= target)[0]
    if len(i) == 0:
        return float("inf")
    i = i[0]
    return float(t[i - 1] + (t[i] - t[i - 1]) * (h[i - 1] - target) / (h[i - 1] - h[i]))


def r2fit(t, H, F, tmax):
    h = H - F
    m = (t <= tmax) & (h > 0)
    if m.sum() < 4:
        return float("nan"), float("nan")
    p = np.polyfit(t[m], np.log(h[m]), 1)
    pred = np.polyval(p, t[m]); y = np.log(h[m])
    return float(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()), float(-p[0])


def main():
    R = {"runs": {}, "criteria": {}}
    runs = {}
    for M in (1, 4, 16, 25):
        for s in ("psi4", "unif"):
            d = load(f"M{M}_{s}")
            if d is not None:
                runs[f"M{M}_{s}"] = d
    for extra in ("M16_eq", "M16_cos2", "M16_cos12"):
        d = load(extra)
        if d is not None:
            runs[extra] = d
    for k, d in runs.items():
        n = int(d["n"]); t = d["t"]
        row = {"vortices_mean": float(d["vort"].mean()), "secs": float(d["secs"]), "t_end": float(d["t_end"])}
        for K in (4, 8, 16, 32):
            F = (K * K - 1) / (2 * n); H = d[f"H{K}"]
            row[f"H{K}_0"] = float(H[0]); row[f"H{K}_end"] = float(H[-1]); row[f"F{K}"] = F
            row[f"tau{K}"] = tau(t, H, F)
        F = (16 * 16 - 1) / (2 * n)
        row["R2_16"], row["rate_fit16"] = r2fit(t, d["H16"], F, 3 * row["tau16"] if np.isfinite(row["tau16"]) else t[-1])
        R["runs"][k] = row
    rr = R["runs"]; C = R["criteria"]
    if "M1_psi4" in rr and "M1_unif" in rr:
        v = [rr[f"M1_{s}"]["H16_end"] / rr[f"M1_{s}"]["H16_0"] for s in ("psi4", "unif")]
        C["C1_control_M1"] = {"ratio_end_over_0": v, "pass": all(x > 0.9 for x in v)}
    if "M16_eq" in rr:
        d = runs["M16_eq"]; F = (255) / (2 * int(d["n"]))
        C["C2_equivariance"] = {"max_H16": float(d["H16"].max()), "3F": 3 * F, "pass": bool(d["H16"].max() <= 3 * F)}
    vals = {}
    for M in (16, 25):
        for s in ("psi4", "unif"):
            k = f"M{M}_{s}"
            if k in rr:
                r = rr[k]; F = r["F16"]
                vals[k] = (r["H16_end"] - F) / (r["H16_0"] - F)
    C["C3_relax"] = {"remaining_fraction": vals, "pass": bool(vals) and all(v < 0.1 for v in vals.values())}
    rate = {}; vort = {}
    for M in (4, 16, 25):
        k = f"M{M}_psi4"
        if k in rr:
            rate[M] = 1 / rr[k]["tau16"]; vort[M] = rr[k]["vortices_mean"]
    if len(rate) == 3:
        order_ok = rate[4] < rate[16] and rate[16] <= 1.2 * rate[25]
        rank_r = np.argsort(np.argsort([rate[M] for M in (4, 16, 25)]))
        rank_v = np.argsort(np.argsort([vort[M] for M in (4, 16, 25)]))
        C["C4_order"] = {"rate_psi4": rate, "vortices": vort,
                         "rate_unif": {M: 1 / rr[f"M{M}_unif"]["tau16"] for M in (4, 16, 25) if f"M{M}_unif" in rr},
                         "vortices_increase": bool(vort[4] < vort[16] < vort[25]),
                         "same_rank": bool((rank_r == rank_v).all()),
                         "pass": bool(order_ok and vort[4] < vort[16] < vort[25] and (rank_r == rank_v).all())}
    r2 = {k: rr[k]["R2_16"] for k in rr if k.startswith(("M16_psi4", "M16_unif", "M25_psi4", "M25_unif"))}
    C["C5_exponential"] = {"R2": r2, "pass": bool(r2) and all(v >= 0.9 for v in r2.values())}
    # test 2
    lag = {}
    for k in ("M16_psi4", "M25_psi4"):
        if k in rr:
            lag[k] = {"tau4": rr[k]["tau4"], "tau16": rr[k]["tau16"], "ratio": rr[k]["tau4"] / rr[k]["tau16"]}
    C["T2a_scale_lag"] = {"values": lag, "pass": bool(lag) and all(v["ratio"] > 1.2 for v in lag.values())}
    if "M16_cos2" in rr and "M16_cos12" in rr:
        a, b = rr["M16_cos2"]["tau32"], rr["M16_cos12"]["tau32"]
        C["T2b_wavelength_lag"] = {"tau_long_k2": a, "tau_short_k12": b, "ratio": a / b, "pass": bool(a > 1.2 * b)}
    # test 3
    s1 = os.path.join(OUT, "relax_slit_1d.json"); s2 = os.path.join(OUT, "relax_slit_2d.json")
    d1 = json.load(open(s1)) if os.path.exists(s1) else None
    d2 = json.load(open(s2)) if os.path.exists(s2) else None
    if d1:
        below = [w < m - 3 * s for w, (m, s) in zip(d1["wrong"], d1["band"])]
        C["T3a_1D"] = {"wrong_corr": d1["wrong"], "band": d1["band"], "below_band_each_z": below,
                       "pass": all(below),
                       "posthoc_lr_corr": d1["wrong_lr"], "posthoc_KS": d1["ks"]}
    if d2:
        g = {k: (v["wrong"], v["eq"], v["band"]) for k, v in d2.items() if k.startswith("M")}
        def inband(k):
            w, e, (m, s) = g[k]; return w >= m - 3 * s
        keys = list(g)
        m16 = [k for k in keys if k.startswith("M16")]; m1 = [k for k in keys if k.startswith("M1_")]
        p = (inband(m16[-1]) and not inband(m16[0]) and not any(inband(k) for k in m1))
        C["T3b_2D"] = {"runs": {k: {"wrong": g[k][0], "eq": g[k][1], "ideal": g[k][2]} for k in keys}, "pass": bool(p)}
    json.dump(R, open(os.path.join(OUT, "relax_results.json"), "w"), indent=1, default=float)
    for k, v in C.items():
        print(k, "PASS" if v.get("pass") else "FAIL", {a: b for a, b in v.items() if a != "pass"})

    # ---------------- figure
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
    col = {1: "#888888", 4: "#1f77b4", 16: "#d62728", 25: "#2ca02c"}
    a = ax[0, 0]
    for M in (1, 4, 16, 25):
        for s, ls in (("psi4", "-"), ("unif", "--")):
            k = f"M{M}_{s}"
            if k in runs:
                d = runs[k]; F = 255 / (2 * int(d["n"]))
                a.semilogy(d["t"], np.maximum(d["H16"] - F, 1e-4), ls, color=col[M],
                           label=f"M={M} {'|psi|^4' if s=='psi4' else 'uniforme'} ({rr[k]['vortices_mean']:.0f} vortici)")
    if "M16_eq" in runs:
        d = runs["M16_eq"]; F = 255 / (2 * int(d["n"]))
        a.semilogy(d["t"], np.maximum(d["H16"] - F, 1e-4), ":", color="k", label="M=16 partenza |psi|^2 (controllo)")
    a.set_xlabel("t"); a.set_ylabel("H̄(t) - pavimento (celle 16x16)"); a.set_title("1. Rilassamento: H̄(t) vs numero di modi M")
    a.legend(fontsize=7); a.set_ylim(1e-4, 3)
    # snapshot
    a = ax[0, 1]
    if "M16_psi4" in runs:
        d = runs["M16_psi4"]; box = Box(16, seed=1); te = float(d["t_end"])
        ng = 200; psi = box.dens_grid(te, ng)
        a.imshow((np.abs(psi) ** 2).T, origin="lower", extent=[0, PI, 0, PI], cmap="Greys", alpha=0.8)
        ph = np.angle(psi); w = lambda p, q: (q - p + PI) % (2 * PI) - PI
        sw = (w(ph[:-1, :-1], ph[1:, :-1]) + w(ph[1:, :-1], ph[1:, 1:]) + w(ph[1:, 1:], ph[:-1, 1:]) + w(ph[:-1, 1:], ph[:-1, :-1]))
        I, J = np.nonzero(np.abs(sw) > PI); g = (np.arange(ng) + 1.0) * PI / ng
        a.scatter(d["x"][:4000], d["y"][:4000], s=0.6, color="#d62728", alpha=0.5, label="surfisti (partiti |psi|^4)")
        a.scatter(g[I], g[J], marker="o", s=60, facecolors="none", edgecolors=np.where(sw[I, J] > 0, "b", "orange"),
                  linewidths=1.8, label="vortici (blu +, arancio -)")
        a.set_title(f"Istantanea M=16, t={te:.1f}: surfisti mescolati su |psi|^2"); a.legend(fontsize=7, loc="upper right")
    # scale lag
    a = ax[0, 2]
    if "M16_psi4" in runs:
        d = runs["M16_psi4"]
        for K, c in zip((4, 8, 16, 32), ("#08306b", "#2171b5", "#6baed6", "#c6dbef")):
            F = (K * K - 1) / (2 * int(d["n"]))
            h = (d[f"H{K}"] - F) / (d[f"H{K}"][0] - F)
            a.semilogy(d["t"], np.maximum(h, 1e-3), color=c, label=f"celle {K}x{K} (tau={rr['M16_psi4'][f'tau{K}']:.2f})")
        for k, ls in (("M16_cos2", "-."), ("M16_cos12", ":")):
            if k in runs:
                d = runs[k]; F = 1023 / (2 * int(d["n"]))
                h = (d["H32"] - F) / (d["H32"][0] - F)
                a.semilogy(d["t"], np.maximum(h, 1e-3), ls, color="k", label=f"{k[4:]} su 32x32 (tau={rr[k]['tau32']:.2f})")
        a.set_title("2. Le scale grandi restano indietro?"); a.set_xlabel("t"); a.set_ylabel("H̄/H̄(0) (sopra il pavimento)")
        a.legend(fontsize=7); a.set_ylim(1e-3, 1.5)
    # double slit 1D
    a = ax[1, 0]
    if d1:
        z = d1["screens"]
        for k, c, lab in (("eq", "k", "partenza |psi|^2"), ("wrong", "#1f77b4", "uniforme nelle fenditure"),
                          ("wrong4", "#2ca02c", "|psi|^4"), ("wrong_lr", "#d62728", "2/3 a sinistra")):
            a.plot(z, d1["ks"][k], "o-", color=c, label=lab)
        a.set_xscale("log"); a.set_xlabel("distanza dello schermo z"); a.set_ylabel("distanza KS clic vs |psi|^2")
        a.set_title("3a. Doppia fenditura 1D (nessun vortice): l'errore resta congelato"); a.legend(fontsize=7)
    a = ax[1, 1]
    if d1:
        x = np.array(d1["screen_x"]); r = np.array(d1["screen_ref"])
        e = np.linspace(x[0], x[-1], 121); cx = 0.5 * (e[1:] + e[:-1])
        for k, c, lab in (("screen_eq", "k", "partenza |psi|^2"), ("screen_lr", "#d62728", "2/3 a sinistra")):
            h, _ = np.histogram(d1[k], e, density=True); a.step(cx, h, where="mid", color=c, label=lab)
        a.plot(x, r / (r.sum() * (x[1] - x[0])), color="0.6", lw=0.8, label="|psi|^2 esatto")
        a.set_title("3a. Schermo a z=20: frange storte per sempre"); a.legend(fontsize=7); a.set_xlabel("x")
    a = ax[1, 2]
    if d2:
        keys = [k for k in d2 if k.startswith("M")]
        xs = np.arange(len(keys))
        a.bar(xs - 0.2, [d2[k]["wrong"] for k in keys], 0.4, color="#d62728", label="partenza sbagliata (uniforme)")
        a.bar(xs + 0.2, [d2[k]["eq"] for k in keys], 0.4, color="0.5", label="partenza |psi|^2 alle fenditure")
        for i, k in enumerate(keys):
            m, s = d2[k]["band"]; a.plot([i - 0.45, i + 0.45], [m - 3 * s] * 2, "k--", lw=0.8)
        a.set_xticks(xs); a.set_xticklabels([k.replace("_Z", "\nZ0=") for k in keys], fontsize=8)
        lo = min(min(d2[k]["wrong"] for k in keys), 0.9); a.set_ylim(lo - 0.01, 1.0)
        a.set_ylabel("correlazione frange"); a.set_title("3b. Guida multimodo prima delle fenditure (2D)"); a.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "relax_figure.png"), dpi=110)
    print("figure written")


if __name__ == "__main__":
    main()

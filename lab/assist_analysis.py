"""assist_analysis.py - grades the pre-registered exams A1-A3, O1-O3 and builds the report numbers (A4, O4).
Reads lab/out/assist_A.json, assist_O.json, assist_specA.json, assist_specO.json, assist_wkb.json (if present).
Writes lab/out/assist_analysis.json."""
import json, os, sys
import numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assist_core import S_wi, OUT, PI

pre = json.load(open(os.path.join(OUT, "assist_prereg.json")))
L = lambda n: json.load(open(os.path.join(OUT, f"assist_{n}.json"))) if os.path.exists(os.path.join(OUT, f"assist_{n}.json")) else None
A, O, SA, SO, W = L("A"), L("O"), L("specA"), L("specO"), L("wkb")
GS = pre["A_sauter"]["gs"]
E1s = np.array(pre["A_sauter"]["E1_fit"])
x = 1 / E1s
out = {}


def fit3(lq):
    M = np.vstack([-x, np.log(x), np.ones_like(x)]).T
    c, *_ = np.linalg.lstsq(M, lq, rcond=None)
    return dict(S=float(c[0]), nu=float(c[1]), c=float(c[2]), maxres=float(np.max(np.abs(lq - M @ c))))


def Q(r, src, e=0.01):
    return (r[f"eps{e}"][src] + r[f"eps{-e}"][src] - 2 * r["eps0.0"][src]) / (2 * e * e)


# ------------------------------------------------------------------ A (weak Sauter pulse)
if A:
    a1 = dict(maxdev_f_th0025=0.0, maxdev_f_rich=0.0, maxdev_Q=0.0, order_ratios=[])
    for key, r in A.items():
        for e, v in r.items():
            if not isinstance(v, dict):
                continue
            a1["maxdev_f_th0025"] = max(a1["maxdev_f_th0025"], abs(v["th0.025"] / v["dirac"] - 1))
            a1["maxdev_f_rich"] = max(a1["maxdev_f_rich"], abs(v["rich"] / v["dirac"] - 1))
            a1["order_ratios"].append((v["th0.1"] / v["dirac"] - 1) / (v["th0.05"] / v["dirac"] - 1))
        a1["maxdev_Q"] = max(a1["maxdev_Q"], abs(Q(r, "rich") / Q(r, "dirac") - 1))
    a1["order_ratio_range"] = [float(np.min(a1["order_ratios"])), float(np.max(a1["order_ratios"]))]
    del a1["order_ratios"]
    a1["dt_check_maxdev"] = max(abs(r["dt_check_eps0.01"] / r["eps0.01"]["dirac"] - 1) for k, r in A.items()
                                if "dt_check_eps0.01" in r)
    a1["Q_eps_check_maxdev"] = max(abs(Q(r, "rich", 0.005) / Q(r, "rich", 0.01) - 1) for r in A.values())
    a1["PASS"] = bool(a1["maxdev_f_th0025"] <= 0.03 and a1["maxdev_f_rich"] <= 0.01 and a1["maxdev_Q"] <= 0.05)
    out["A1"] = a1
    # POST-HOC (declared, not graded): same exam without the points where f < 1e-20 (destructive interference of
    # a0 and eps*a1 at eps < 0 near threshold, E1 = 0.0625, f down to 1e-25, close to the round-off floor)
    ph = dict(maxdev_f_th0025=0.0, maxdev_f_rich=0.0, n_excluded=0)
    for key, r in A.items():
        for e, v in r.items():
            if not isinstance(v, dict):
                continue
            if v["dirac"] < 1e-20:
                ph["n_excluded"] += 1
                continue
            ph["maxdev_f_th0025"] = max(ph["maxdev_f_th0025"], abs(v["th0.025"] / v["dirac"] - 1))
            ph["maxdev_f_rich"] = max(ph["maxdev_f_rich"], abs(v["rich"] / v["dirac"] - 1))
    ph["n_total"] = sum(1 for r in A.values() for v in r.values() if isinstance(v, dict))
    out["A1_posthoc"] = ph
    gams = pre["A_sauter"]["gammas"]
    a2 = {}
    S0 = fit3(np.log([A[f"g{gams[0]}_E{E}"]["eps0.0"]["rich"] for E in E1s]))
    for g in gams:
        lq = np.log([Q(A[f"g{g}_E{E}"], "rich") for E in E1s])
        lqd = np.log([Q(A[f"g{g}_E{E}"], "dirac") for E in E1s])
        fr, fd = fit3(lq), fit3(lqd)
        a2[str(g)] = dict(S_fit=fr["S"], nu=fr["nu"], c=fr["c"], maxres=fr["maxres"], S_fit_dirac=fd["S"],
                          S_wi=S_wi(g, GS), ratio=fr["S"] / S_wi(g, GS),
                          enh_exponent_E01=float((S_wi(0.1, GS) - S_wi(g, GS)) / 0.1))
    a2["S0_fit_eps0"] = S0
    a2["PASS"] = bool(all(abs(a2[str(g)]["ratio"] - 1) <= 0.03 for g in gams if g <= 6))
    out["A2"] = a2
    gc_pred = pre["A_sauter"]["predictions"]["gamma_crit_gs"]
    gg = [g for g in gams if g <= 6]
    Sf = np.array([a2[str(g)]["S_fit"] for g in gg])
    obj = lambda gc: float(np.sum((Sf - np.array([S_wi(g * gc_pred / gc, GS) for g in gg])) ** 2))
    r = minimize_scalar(obj, bounds=(1.0, 2.5), method="bounded", options=dict(xatol=1e-6))
    out["A3"] = dict(gc_fit=float(r.x), gc_pred=gc_pred, ratio=float(r.x / gc_pred), rms=float(np.sqrt(r.fun / len(gg))),
                     PASS=bool(abs(r.x / gc_pred - 1) <= 0.05))
    # zeroth + first order in eps (Torgrimsson et al. 2017) vs exact at eps = 0.1, E1 = 0.1
    t = {}
    for g in gams:
        rr = A[f"g{g}_E0.1"]
        f0 = rr["eps0.0"]["rich"]
        lin = (rr["eps0.01"]["rich"] - rr["eps-0.01"]["rich"]) / 0.02
        q = Q(rr, "rich")
        f1 = f0 + 0.1 * lin + 0.01 * q
        t[str(g)] = dict(f_exact=rr["eps0.1"]["rich"], f_first_order=f1, ratio=rr["eps0.1"]["rich"] / f1,
                         enh_eps01=rr["eps0.1"]["rich"] / f0, enh_eps001=rr["eps0.01"]["rich"] / f0)
    out["A4_first_order_E01"] = t

# ------------------------------------------------------------------ O (weak oscillating field)
if O:
    o1 = dict(maxdev_f_th0025=0.0, maxdev_f_rich=0.0)
    for key, r in O.items():
        for e, v in r.items():
            o1["maxdev_f_th0025"] = max(o1["maxdev_f_th0025"], abs(v["th0.025"] / v["dirac"] - 1))
            o1["maxdev_f_rich"] = max(o1["maxdev_f_rich"], abs(v["rich"] / v["dirac"] - 1))
    o1["PASS"] = bool(o1["maxdev_f_th0025"] <= 0.03 and o1["maxdev_f_rich"] <= 0.01)
    out["O1"] = o1
    o2, ok2, ok3 = {}, True, True
    gams = pre["O_cos"]["gammas"]
    S0 = fit3(np.log([O[f"g{gams[0]}_E{E}"]["eps0.0"]["rich"] for E in E1s]))
    for eps in pre["O_cos"]["eps"]:
        for g in gams:
            pr = pre["O_cos"]["predictions"][f"eps{eps}_g{g}"]
            fr = fit3(np.log([O[f"g{g}_E{E}"][f"eps{eps}"]["rich"] for E in E1s]))
            d = dict(S_fit=fr["S"], nu=fr["nu"], maxres=fr["maxres"], S_WKB=pr["S_WKB"], ratio_WKB=fr["S"] / pr["S_WKB"],
                     s_fit=fr["S"] / S0["S"], s_Linder=pr["s_Linder"],
                     ratio_Linder=(fr["S"] / S0["S"]) / pr["s_Linder"] if pr["s_Linder"] > 0 else None,
                     graded_O2=pr["graded_O2"], graded_O3=pr["graded_O3"],
                     f0_E01=O[f"g{g}_E0.1"][f"eps{eps}"]["rich"], enh_E01=O[f"g{g}_E0.1"][f"eps{eps}"]["rich"] / O[f"g{g}_E0.1"]["eps0.0"]["rich"])
            if pr["graded_O2"]:
                ok2 &= abs(d["ratio_WKB"] - 1) <= 0.03
            if pr["graded_O3"]:
                ok3 &= abs(d["ratio_Linder"] - 1) <= 0.03
            o2[f"eps{eps}_g{g}"] = d
    o2["S0_fit"] = S0
    out["O2_O3"] = o2
    out["O2_PASS"] = bool(ok2); out["O3_PASS"] = bool(ok3)

# ------------------------------------------------------------------ spectra: pair numbers and enhancement at E1 = 0.1
for name, S in (("specA", SA), ("specO", SO)):
    if not S:
        continue
    p = np.array(S["p"])
    N = lambda f: float(np.trapezoid(np.array(f), p) / (2 * PI))
    rich = lambda d: (4 * np.array(d["th0.025"]) - np.array(d["th0.05"])) / 3
    N1 = N(rich(S["strong_only"])); N1d = N(S["strong_only"]["dirac"])
    tab = dict(N_strong_rule=N1, N_strong_dirac=N1d)
    for key, r in S.items():
        if not key.startswith("eps"):
            continue
        Nb, Nw = N(rich(r["both"])), N(rich(r["weak_only"]))
        Nbd, Nwd = N(r["both"]["dirac"]), N(r["weak_only"]["dirac"])
        d = dict(N_both=Nb, N_weak=Nw, N_both_dirac=Nbd, N_weak_dirac=Nwd, enhancement=Nb / (N1 + Nw),
                 conv_th01=N(r["both"]["th0.1"]) / Nbd, conv_th005=N(r["both"]["th0.05"]) / Nbd,
                 conv_th0025=N(r["both"]["th0.025"]) / Nbd, rich_over_dirac=Nb / Nbd)
        if "exact" in r["weak_only"]:
            d["N_weak_exact"] = N(r["weak_only"]["exact"])
        if W and name == "specA" and key in W:
            d["N_WKB"] = W[key]["N"]; d["N_WKB_over_rule"] = W[key]["N"] / Nb
        tab[key] = d
    out[name] = tab

json.dump(out, open(os.path.join(OUT, "assist_analysis.json"), "w"), indent=1)
for k, v in out.items():
    if isinstance(v, dict) and "PASS" in v:
        print(k, "PASS" if v["PASS"] else "FAIL", {a: b for a, b in v.items() if not isinstance(b, dict)})
for k in ("O2_PASS", "O3_PASS"):
    if k in out:
        print(k, out[k])

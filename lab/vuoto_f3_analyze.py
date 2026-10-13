"""vuoto_f3_analyze.py - grades phase 3 (lab/out/vuoto_f3_prereg.json). Reads vuoto_f3_dv.json, vuoto_f3_g2.json,
vuoto_f3_bands.json and (for the first 3 g-2 lattice spacings) vuoto_f2_se.json. Writes lab/out/vuoto_f3_exams.json."""
import os, json
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
L = lambda n: json.load(open(os.path.join(OUT, n)))
DV, G3, BANDS, SE = L("vuoto_f3_dv.json"), L("vuoto_f3_g2.json"), L("vuoto_f3_bands.json"), L("vuoto_f2_se.json")
THS = [0.05, 0.025, 0.0125, 0.00625, 0.003125]
R = {}


def fit(y):
    A = np.stack([np.ones(len(THS)), THS], 1)
    c, res, *_ = np.linalg.lstsq(A, np.array(y), rcond=None)
    r = np.array(y) - A @ c
    s = np.sqrt(np.sum(r**2) / (len(y) - 2) * np.linalg.inv(A.T @ A)[0, 0])
    return float(c[0]), float(c[1]), float(s)


tab = {}
for rule, cut in (("R", "K12"), ("R", "K12/3"), ("C", "K12/3")):
    rows = []
    for th in THS:
        d = DV["%s|%s|%g" % (rule, cut, th)]
        v = d["dv_over_alpha"]
        rows.append({"theta": th, "dv": v, "dvT": d["dvT_over_alpha"], "dvC": d["dvC_over_alpha"],
                     "aniso": v[1] - 0.5 * (v[0] + v[2]), "xz": v[0] - v[2], "iso": float(np.mean(v)),
                     "spread": max(v) - min(v), "dm_over_m": d["dmT_over_m"] + d["dmC_over_m"]})
    t = {"rows": rows}
    t["aniso_fit"] = fit([r["aniso"] for r in rows])
    t["iso_fit"] = fit([r["iso"] for r in rows])
    tab["%s|%s" % (rule, cut)] = t
for cut in ("K12", "K12/3"):
    y = [DV["cont|%s|%g" % (cut, th)]["dv_over_alpha"] for th in THS]
    tab["cont|%s" % cut] = {"dv": y, "iso_fit": fit(y)}
R["table"] = tab
# A1
R["A1_PASS"] = bool(all(abs(r["xz"]) <= 1e-3 for k in ("R|K12", "R|K12/3") for r in tab[k]["rows"]) and
                    all(abs(r["dv"][1] - r["dv"][0]) >= 1e-2 for r in tab["R|K12"]["rows"]))
# A2, A3, A6 decisions
for k in ("R|K12", "R|K12/3"):
    R["A2_aniso_a0|" + k] = tab[k]["aniso_fit"][0]; R["A2_survives|" + k] = bool(abs(tab[k]["aniso_fit"][0]) >= 0.01)
for k in ("R|K12", "R|K12/3", "C|K12/3", "cont|K12", "cont|K12/3"):
    R["A3_iso_a0|" + k] = tab[k]["iso_fit"][0]; R["A3_survives|" + k] = bool(abs(tab[k]["iso_fit"][0]) >= 0.01)
# A4, A5
R["A4_spreads"] = [r["spread"] for r in tab["C|K12/3"]["rows"]]
R["A4_PASS"] = bool(max(R["A4_spreads"]) <= 1e-3)
R["A5_ratios"] = [c["dm_over_m"] / r["dm_over_m"] for c, r in zip(tab["C|K12/3"]["rows"], tab["R|K12/3"]["rows"])]
R["A5_PASS"] = bool(all(abs(x - 1) <= 0.01 for x in R["A5_ratios"]))
# A7
R["A7"] = BANDS["A7"]; R["A7_PASS"] = BANDS["A7_PASS"]
# A8: g-2 spread, lambda -> 0 (linear from 1e-3 and 3e-4), slice 0 at 4 spacings (+ slice 1 at 0.0125)


def pairs(get):
    out = {}
    for lam in (1e-3, 3e-4):
        out[lam] = {}
        for g in "AB":
            out[lam].update(get(lam, g)["pairs"])
    return {p: out[3e-4][p]["F2_over_alpha2pi"] + (out[3e-4][p]["F2_over_alpha2pi"] - out[1e-3][p]["F2_over_alpha2pi"]) * 3e-4 / 7e-4
            for p in out[1e-3]}


g2 = {}
for th in (0.05, 0.025, 0.0125):
    g2["slice0|%g" % th] = pairs(lambda lam, g: SE["g2|lat|%g|%g|std|%s" % (th, lam, g)])
g2["slice0|0.00625"] = pairs(lambda lam, g: G3["g2|lat|slice0|0.00625|%g|%s" % (lam, g)])
g2["slice1|0.0125"] = pairs(lambda lam, g: G3["g2|lat|slice1|0.0125|%g|%s" % (lam, g)])
R["g2_pairs"] = g2
sp = {k: max(v.values()) - min(v.values()) for k, v in g2.items()}
mean = {k: float(np.mean(list(v.values()))) for k, v in g2.items()}
R["g2_spread"] = sp; R["g2_mean"] = mean
ths4 = [0.05, 0.025, 0.0125, 0.00625]
A = np.stack([np.ones(4), ths4], 1)
c, *_ = np.linalg.lstsq(A, np.array([sp["slice0|%g" % t] for t in ths4]), rcond=None)
R["A8_s0"] = float(c[0]); R["A8_survives"] = bool(c[0] >= 0.02)
R["A8_slice_dependence_mean"] = mean["slice1|0.0125"] - mean["slice0|0.0125"]
# A9
R["A9"] = DV["check|K12/3|0.05"]
R["A9_PASS"] = bool(all(v["imag_over_dT"] <= 1e-12 for k, v in DV["check|K12/3|0.05"].items() if k in ("R", "C")))
json.dump(R, open(os.path.join(OUT, "vuoto_f3_exams.json"), "w"), indent=1)
for k, v in R.items():
    if k not in ("table", "g2_pairs"):
        print(k, v)
for k, t in tab.items():
    if "rows" in t:
        print(k)
        for r in t["rows"]:
            print("  th %-9g dv %s  aniso %+.5f iso %+.5f xz %.1e dm/m %.6f" % (r["theta"], np.round(r["dv"], 5), r["aniso"], r["iso"], r["xz"], r["dm_over_m"]))
    else:
        print(k, np.round(t["dv"], 5))
for k, v in g2.items():
    print(k, {p: round(x, 4) for p, x in v.items()})

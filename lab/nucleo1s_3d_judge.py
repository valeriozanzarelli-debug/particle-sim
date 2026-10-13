"""nucleo1s_3d_judge.py - judges phase 2 (3-D) against lab/out/nucleo1s_3d_prereg.json and nucleo1s_3d_amend1.json,
writes the verdicts into lab/out/nucleo1s_3d_results.json ('judge') and the figure lab/out/nucleo1s_3d_figure.png."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, dirac_formula  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R = json.load(open(os.path.join(OUT, "nucleo1s_3d_results.json")))
P = json.load(open(os.path.join(OUT, "nucleo1s_3d_prereg.json")))
AM = json.load(open(os.path.join(OUT, "nucleo1s_3d_amend1.json")))
Za = 60 * ALPHA
D1 = dirac_formula(1, -1, Za)
B1D = 1 - D1
J = {}


def rel(om):
    return (1 - om) / B1D - 1


def U_of(L):
    """the nucleus-blind upper pair: the highest-omega lines whose summed weight reaches ~1/2"""
    s = sorted(L, key=lambda x: -x["omega"])
    w, o = 0.0, 0.0
    for x in s:
        w += x["weight"]; o += x["weight"] * x["omega"]
        if w >= 0.4:
            break
    return o / w


# ---- B1, B2 (lab), B3 (a = 0.16)
J["B1"] = R["lab"]["B1"]; J["B2"] = R["lab"]["B2"]
L16 = sorted(R["scan_bare_0.16"]["lines"], key=lambda x: x["omega"])
U16 = U_of(L16)
d1, d2 = U16 - L16[0]["omega"], U16 - L16[1]["omega"]
e = P["exams"]["B3"]
J["B3"] = dict(U_L1=d1, pred=e["pred_U_L1"], rel1=d1 / e["pred_U_L1"] - 1, U_L2=d2, pred2=e["pred_U_L2"],
               rel2=d2 / e["pred_U_L2"] - 1,
               PASS=bool(abs(d1 / e["pred_U_L1"] - 1) <= e["tol_U_L1"] and abs(d2 / e["pred_U_L2"] - 1) <= e["tol_U_L2"]))
# ---- learning
J["N1"] = R["learn"].get("N1")
J["learned_s"] = R["learn"]["s"]
A_L = [0.5, 0.4, 0.3, 0.24, 0.2]
tab = {}
for a in [0.5, 0.4, 0.3, 0.24, 0.2, 0.16]:
    b = R.get(f"scan_bare_{a}"); l = R.get(f"scan_learned_{a}")
    row = {}
    if b:
        row.update(bare_centroid=b["quartet"]["centroid"], bare_rel=rel(b["quartet"]["centroid"]), bare_std=b["quartet"]["std"],
                   bare_spread=b["quartet"]["spread"], U=U_of(b["lines"]), U_rel=rel(U_of(b["lines"])),
                   L=[x["omega"] for x in sorted(b["lines"], key=lambda x: x["omega"])])
        row["valley_part"] = row["bare_rel"] - row["U_rel"]
    if l:
        strong = [x for x in l["lines"] if x["weight"] >= 0.05]
        row.update(learned_centroid=l["quartet"]["centroid"], learned_rel=rel(l["quartet"]["centroid"]),
                   learned_std=l["quartet"]["std"],
                   learned_spread=float(max(x["omega"] for x in strong) - min(x["omega"] for x in strong)))
    tab[str(a)] = row
J["table"] = tab
J["N2"] = dict(ratios={a: tab[str(a)]["learned_std"] / tab[str(a)]["bare_std"] for a in A_L[1:] if "learned_std" in tab[str(a)]})
J["N2"]["PASS"] = bool(len(J["N2"]["ratios"]) == 4 and all(v <= P["exams"]["N2"]["max_ratio"] for v in J["N2"]["ratios"].values()))
J["N3"] = dict(rels={a: tab[str(a)].get("learned_rel") for a in A_L})
J["N3"]["PASS"] = bool(all(v is not None and abs(v) <= 0.01 for v in J["N3"]["rels"].values()))


def a1_extrap(vals):
    aa = np.array([0.4, 0.3, 0.24])
    return float(np.polyfit(aa**2, np.array(vals), 2)[-1])


bind_l = [1 - tab[str(a)]["learned_centroid"] for a in (0.4, 0.3, 0.24)]
ex = a1_extrap(bind_l)
J["N4"] = dict(extrap_binding=ex, rel=ex / B1D - 1, PASS=bool(abs(ex / B1D - 1) <= 0.01))
lab = R.get("scan_learned_0.5")
if lab:
    strong = [x for x in lab["lines"] if x["weight"] >= 0.05]
    sp = max(x["omega"] for x in strong) - min(x["omega"] for x in strong)
    J["N5"] = dict(spread=sp, period_ticks=(2 * np.pi / sp / 0.5) if sp > 0 else np.inf, PASS=bool(sp <= P["exams"]["N5"]["tol"]))
# ---- valley budget (amendment)
fr = {a: tab[str(a)]["valley_part"] / abs(tab[str(a)]["bare_rel"]) for a in (0.4, 0.3, 0.24, 0.2) if "valley_part" in tab[str(a)]}
J["V1"] = dict(frac=fr, PASS=bool(len(fr) == 4 and all(v >= 0.5 for v in fr.values())))
Ub = [1 - tab[str(a)]["U"] for a in (0.4, 0.3, 0.24)]
Cb = [1 - tab[str(a)]["bare_centroid"] for a in (0.4, 0.3, 0.24)]
J["V2"] = dict(U_extrap_rel=a1_extrap(Ub) / B1D - 1, centroid_extrap_rel=a1_extrap(Cb) / B1D - 1,
               old_A1=-0.116, PASS=bool(abs(a1_extrap(Ub) / B1D - 1) <= 0.01))
J["V3"] = dict(spreads={a: tab[str(a)].get("learned_spread") for a in A_L})
J["V3"]["PASS"] = bool(all(v is not None and v <= 0.002 for v in J["V3"]["spreads"].values()))
# ---- cost (C1): bare law fit K a^{2g} + B a^2 on the bare centroid, a_1%
g = np.sqrt(1 - Za**2)
aa = np.array([a for a in (0.3, 0.24, 0.2, 0.16)])
yy = np.array([tab[str(a)]["bare_rel"] for a in aa])
M = np.stack([aa**(2 * g), aa**2], 1)
KB = np.linalg.lstsq(M, yy, rcond=None)[0]
from scipy.optimize import brentq
a1 = brentq(lambda a: KB[0] * a**(2 * g) + KB[1] * a**2 - 0.01, 1e-3, 0.3)
passing = [a for a in A_L if tab[str(a)].get("learned_rel") is not None and abs(tab[str(a)]["learned_rel"]) <= 0.01]
J["C1"] = dict(bare_fit_K=KB[0], bare_fit_B=KB[1], a_1pct_bare=a1, cells_bare=(24 / a1)**3,
               learned_coarsest_pass=max(passing) if passing else None,
               cells_learned=(24 / max(passing))**3 if passing else None)
J["C1"]["ratio"] = J["C1"]["cells_bare"] / J["C1"]["cells_learned"] if passing else None
J["C1"]["PASS"] = bool(passing and J["C1"]["ratio"] >= 10)
# ---- excited (E1), isotropy (I1)
DIR = {"2s": dirac_formula(2, -1, Za), "2p1/2": dirac_formula(2, 1, Za), "2p3/2,3/2": dirac_formula(2, -2, Za),
       "2p3/2,1/2": dirac_formula(2, -2, Za)}
ex1, iso = {}, {}
ok = True
for a in (0.5, 0.4):
    b, l = R.get(f"exc_bare_{a}"), R.get(f"exc_learned_{a}")
    if not (b and l):
        ok = False; continue
    for k in DIR:
        sb = max(b[k], key=lambda x: x["weight"]) if b[k] else None
        sl = max(l[k], key=lambda x: x["weight"]) if l[k] else None
        if sb is None or sl is None:
            ok = False; continue
        bind = 1 - DIR[k]
        eb, el = abs(sb["omega"] - DIR[k]) / bind, abs(sl["omega"] - DIR[k]) / bind
        ex1[f"{k}@{a}"] = dict(bare=sb["omega"], learned=sl["omega"], dirac=DIR[k], err_bare=eb, err_learned=el,
                               ok=bool(el <= eb + P["exams"]["E1"]["tol"]))
        ok = ok and ex1[f"{k}@{a}"]["ok"]
    for cfg, rr in (("bare", b), ("learned", l)):
        s3 = max(rr["2p3/2,3/2"], key=lambda x: x["weight"])["omega"]
        s1 = max(rr["2p3/2,1/2"], key=lambda x: x["weight"])["omega"]
        iso[f"{cfg}@{a}"] = abs(s3 - s1) / (1 - DIR["2p3/2,3/2"])
J["E1"] = dict(items=ex1, PASS=bool(ok))
if "learned@0.4" in iso:
    J["I1"] = dict(aniso=iso, PASS=bool(iso["learned@0.4"] <= 0.02 and iso["learned@0.4"] <= 1.1 * iso["bare@0.4"] + 1e-12))
# ---- glue
if "glue" in R:
    gl = R["glue"]
    def spread(rec, suf):
        v = [rec[f"J0={j}{suf}"] for j in (3, 4, 5, 6)]
        return (max(v) - min(v)) if all(x is not None for x in v) else None, v
    s1, v1 = spread(gl["Z1_1s"], "")
    J["G1"] = dict(spread=s1, vals=v1, PASS=bool(s1 is not None and s1 <= 2e-4))
    J["G2"] = dict(vals=v1, PASS=bool(s1 is not None and all(abs(x) <= 1e-3 for x in v1)
                                    and all((spread(gl[k], "")[0] or 1) <= 5e-4 for k in ("Z1_2s", "Z2_1s"))))
    c1, cv1 = spread(gl["Z1_1s"], "_cusp")
    c2, cv2 = spread(gl["Z1_2s"], "_cusp")
    c3, cv3 = spread(gl["Z2_1s"], "_cusp")
    J["G1b"] = dict(spread=c1, vals=cv1, PASS=bool(c1 is not None and c1 <= 2e-4))
    J["G2b"] = dict(spread_2s=c2, vals_2s=cv2, spread_He=c3, vals_He=cv3,
                    PASS=bool(None not in (c1, c2, c3) and c2 <= 2e-4 and c3 <= 2e-4 and
                              all(abs(x) <= 1e-3 for x in cv1 + cv2 + cv3)))
    J["glue_phase1_row"] = {k: [gl[k][f"J0={j}_phase1"] for j in (3, 4, 5, 6)] for k in ("Z1_1s", "Z1_2s", "Z2_1s")}
# ---- centre nucleus
for a in (0.5, 0.3):
    c = R.get(f"centre_{a}")
    if c and len(c) == 4:
        J[f"centre_{a}"] = {q: dict(split=v["split"], centroid_rel=rel(v["quartet"]["centroid"]) if v["quartet"] else None,
                                    lines=[(round(x["omega"], 6), round(x["weight"], 3)) for x in v["lines"]]) for q, v in c.items()}
if "centre_0.3" in J:
    corner = tab["0.3"]["bare_spread"]
    J["CN1"] = dict(corner_split=corner, centre_splits={q: v["split"] for q, v in J["centre_0.3"].items()},
                    PASS=bool(all(v["split"] is not None and v["split"] < corner for v in J["centre_0.3"].values())))
    # (a colour whose 1s window [0.80, 0.95] holds no line of weight >= 0.05 has split None: counted as FAIL)
    c0 = sorted(x[0] for x in J["centre_0.3"]["(0, 0, 0)"]["lines"] if x[1] >= 0.05)
    c1_ = sorted(x[0] for x in J["centre_0.3"]["(1, 0, 0)"]["lines"] if x[1] >= 0.05)
    dmax = max(min(abs(x - y) for y in c1_) for x in c0) if c0 and c1_ else None
    J["CN2"] = dict(max_diff=dmax, PASS=bool(dmax is not None and dmax > 1e-4))
R["judge"] = J
json.dump(R, open(os.path.join(OUT, "nucleo1s_3d_results.json"), "w"), indent=1, default=float)
for k, v in J.items():
    if isinstance(v, dict) and "PASS" in v:
        print(k, "PASS" if v["PASS"] else "FAIL", {kk: vv for kk, vv in v.items() if kk not in ("items",)})
print("table:")
for a, r in tab.items():
    print(a, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items() if k != "L"})

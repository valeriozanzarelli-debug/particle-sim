"""nucleo1s_f4_judge.py - judges phase 4 (exams Q1..Q4 of lab/out/nucleo1s_f4_prereg.json); verdicts into
lab/out/nucleo1s_3d_results.json ('judge_f4') and lab/out/nucleo1s_f4_results.json."""
import os, sys, json, fcntl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, dirac_formula  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
RP = os.path.join(OUT, "nucleo1s_3d_results.json")
R = json.load(open(RP))
E = json.load(open(os.path.join(OUT, "nucleo1s_f4_prereg.json")))["exams"]
Za = 60 * ALPHA
D1 = dirac_formula(1, -1, Za)
rel = lambda om: (1 - om) / (1 - D1) - 1


def spread(L, wmin=0.05):
    s = [x["omega"] for x in L if x["weight"] >= wmin and 0.80 <= x["omega"] <= 0.95]
    return (max(s) - min(s)) if s else None


def U_of(L):
    s = sorted(L, key=lambda x: -x["omega"]); w = o = 0.0
    for x in s:
        w += x["weight"]; o += x["weight"] * x["omega"]
        if w >= 0.4:
            break
    return o / w


J = {}
l3 = R["learn3"]
best = min(l3["evals"], key=lambda e: e["obj"])
sp = spread(best["lines"], 0.03)
J["Q1"] = dict(s=l3["s"], U=l3["U_bare"], centroid=best["centroid"], dU=best["centroid"] - l3["U_bare"], spread=sp,
               n_evals=len(l3["evals"]), PASS=bool(sp <= E["Q1"]["spread"] and abs(best["centroid"] - l3["U_bare"]) <= E["Q1"]["dU"]))
tab = {}
for cfg in ("transfer", "transferP"):
    for a in (0.4, 0.24, 0.2, 0.16):
        r = R.get(f"scan_{cfg}_{a}")
        if r:
            tab[f"{cfg}@{a}"] = dict(s=r["s"], rel=r["rel"], spread=spread(r["lines"]), centroid=r["quartet"]["centroid"],
                                    rel_U_bare=rel(U_of(R[f"scan_bare_{a}"]["lines"])) if f"scan_bare_{a}" in R else None,
                                    lines=[(round(x["omega"], 6), round(x["weight"], 3)) for x in r["lines"] if x["weight"] >= 0.01])
test = [f"transfer@{a}" for a in (0.4, 0.24, 0.2)]
ok = all(k in tab and abs(tab[k]["rel"]) <= E["Q2"]["tol"] and tab[k]["spread"] is not None
         and tab[k]["spread"] <= E["Q2"]["spread"] for k in test)
J["Q2"] = dict(PASS=bool(ok), table={k: v for k, v in tab.items()})
b, l = R.get("exc_bare_0.4"), R.get("exc_transfer_0.4")
DIR = {"2s": dirac_formula(2, -1, Za), "2p1/2": dirac_formula(2, 1, Za), "2p3/2,3/2": dirac_formula(2, -2, Za),
       "2p3/2,1/2": dirac_formula(2, -2, Za)}
if b and l:
    eb = {k: abs(max(b[k], key=lambda x: x["weight"])["omega"] - DIR[k]) / (1 - DIR[k]) for k in DIR}
    el = {k: abs(max(l[k], key=lambda x: x["weight"])["omega"] - DIR[k]) / (1 - DIR[k]) for k in DIR}
    J["Q3"] = dict(err_bare=eb, err_transfer=el, mean_bare=sum(eb.values()) / 4, mean_transfer=sum(el.values()) / 4,
                   PASS=bool(sum(el.values()) <= sum(eb.values())))
if ok:
    cb = R["judge"]["C1"]["cells_bare"]
    cl = (24 / 0.4)**3
    J["Q4"] = dict(cells_bare=cb, cells_transfer=cl, ratio=cb / cl, PASS=bool(cb / cl >= E["Q4"]["min_ratio"]))
else:
    J["Q4"] = dict(PASS=False, note="not judgeable: Q2 fails")
print(json.dumps({k: (v if k != "Q2" else {"PASS": v["PASS"]}) for k, v in J.items()}, indent=1))
for k, v in tab.items():
    print(k, v["s"], f"{100 * v['rel']:+.3f} %", v["spread"], "U_bare", v["rel_U_bare"])
with open(RP + ".lock", "w") as lk:
    fcntl.flock(lk, fcntl.LOCK_EX)
    R = json.load(open(RP)); R["judge_f4"] = J
    json.dump(R, open(RP + ".tmp", "w")); os.replace(RP + ".tmp", RP)
json.dump(J, open(os.path.join(OUT, "nucleo1s_f4_results.json"), "w"), indent=1)

"""nucleo1s_f3_judge.py - judges phase 3 (exams M1..M5 of lab/out/nucleo1s_f3_prereg.json); verdicts into
lab/out/nucleo1s_3d_results.json ('judge_f3') and a compact copy lab/out/nucleo1s_f3_results.json."""
import os, sys, json, fcntl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, dirac_formula  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
RP = os.path.join(OUT, "nucleo1s_3d_results.json")
R = json.load(open(RP))
P = json.load(open(os.path.join(OUT, "nucleo1s_f3_prereg.json")))
E = P["exams"]
Za = 60 * ALPHA
D1 = dirac_formula(1, -1, Za)
rel = lambda om: (1 - om) / (1 - D1) - 1
J = {}


def spread(L, wmin=0.05):
    s = [x["omega"] for x in L if x["weight"] >= wmin and 0.80 <= x["omega"] <= 0.95]
    return (max(s) - min(s)) if s else None


l2 = R.get("learn2", {})
if l2.get("done"):
    best = min(l2["evals"], key=lambda e: e["obj"])
    sp = spread(best["lines"], 0.03)
    J["M1"] = dict(s=l2["s"], U=l2["U_bare"], centroid=best["centroid"], dU=best["centroid"] - l2["U_bare"], spread=sp,
                   n_evals=len(l2["evals"]),
                   PASS=bool(sp <= E["M1"]["spread"] and abs(best["centroid"] - l2["U_bare"]) <= E["M1"]["dU"]))
sc = {}
for a in (0.5, 0.4, 0.3, 0.24, 0.2):
    r = R.get(f"scan_learned2_{a}")
    if r:
        sc[str(a)] = dict(rel=r["rel"], centroid=r["quartet"]["centroid"], spread=spread(r["lines"]),
                          rel_learned1=R[f"scan_learned_{a}"]["rel"], rel_bare=R[f"scan_bare_{a}"]["rel"])
judged = [str(a) for a in (0.4, 0.3, 0.24, 0.2)]
if all(k in sc for k in judged):
    J["M2"] = dict(scan=sc, PASS=bool(all(abs(sc[k]["rel"]) <= E["M2"]["tol"] and sc[k]["spread"] is not None
                                          and sc[k]["spread"] <= E["M2"]["spread"] for k in judged)))
    passing = [float(k) for k in judged if abs(sc[k]["rel"]) <= E["M2"]["tol"] and (sc[k]["spread"] or 1) <= E["M2"]["spread"]]
    if J["M2"]["PASS"]:
        cb = R["judge"]["C1"]["cells_bare"]; cl = (24 / max(passing))**3
        J["M5"] = dict(cells_bare=cb, cells_learned2=cl, ratio=cb / cl, PASS=bool(cb / cl >= E["M5"]["min_ratio"]))
    else:
        J["M5"] = dict(PASS=False, note="not judgeable as pre-registered: M2 fails")
else:
    J["M2"] = dict(scan=sc, PASS=None, note="incomplete")
if "0.5" in sc:
    J["M3"] = dict(spread=sc["0.5"]["spread"], PASS=bool(sc["0.5"]["spread"] is not None and
                                                          sc["0.5"]["spread"] <= E["M3"]["max_spread"]))
DIR = {"2s": dirac_formula(2, -1, Za), "2p1/2": dirac_formula(2, 1, Za), "2p3/2,3/2": dirac_formula(2, -2, Za),
       "2p3/2,1/2": dirac_formula(2, -2, Za)}
b, l = R.get("exc_bare_0.4"), R.get("exc_learned2_0.4")
if b and l:
    eb, el = {}, {}
    for k in DIR:
        bind = 1 - DIR[k]
        eb[k] = abs(max(b[k], key=lambda x: x["weight"])["omega"] - DIR[k]) / bind
        el[k] = abs(max(l[k], key=lambda x: x["weight"])["omega"] - DIR[k]) / bind
    mb, ml = sum(eb.values()) / 4, sum(el.values()) / 4
    J["M4"] = dict(err_bare=eb, err_learned2=el, mean_bare=mb, mean_learned2=ml, PASS=bool(ml <= mb))
print(json.dumps(J, indent=1))
with open(RP + ".lock", "w") as lk:
    fcntl.flock(lk, fcntl.LOCK_EX)
    R = json.load(open(RP)); R["judge_f3"] = J
    json.dump(R, open(RP + ".tmp", "w")); os.replace(RP + ".tmp", RP)
json.dump(J, open(os.path.join(OUT, "nucleo1s_f3_results.json"), "w"), indent=1)

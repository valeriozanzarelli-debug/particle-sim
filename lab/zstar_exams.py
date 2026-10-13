"""zstar_exams.py - judge the ZSTAR runs against lab/out/zstar_prereg.json. Writes lab/out/zstar_exams.json.
Base-rule walk values (kappa = -1, s = 1/2, refl0 = -1) are the giro-1 runs lab/out/supercrit_c1.json (same tracker)."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zstar_core import predict_Zcr
from zz_core import ALPHA

HERE = os.path.dirname(os.path.abspath(__file__))
PRE = json.load(open(os.path.join(HERE, "out", "zstar_prereg.json")))
RUN = json.load(open(os.path.join(HERE, "out", "zstar_runs.json")))            # pre-registered tracker
SGN = json.load(open(os.path.join(HERE, "out", "zstar_runs_sign.json")))       # amendment A3 (sign method)
C1 = json.load(open(os.path.join(HERE, "out", "supercrit_c1.json")))
PR = PRE["predictions"]


def walk(v, a, method="sign"):
    """Z_cr of the walk. method 'sign' (A3, used for the exams), 'tracker' (pre-registered; base = giro-1 C1)"""
    if method == "sign":
        x = SGN.get(f"{v}|{a:.6g}")
        if x is not None:
            return x["Zcr"]
    if v == "base":
        k = f"{a:.6g}"
        return C1[k]["Zcr"] if k in C1 else None
    k = f"{v}|{a:.6g}"
    return RUN[k]["Zcr"] if k in RUN else None


def tracker(v, a):
    return walk(v, a, method="tracker")


def match(v, a):
    return PR[v]["Zcr_matching"][f"{a:.6g}"]


def zinf(v, alist=(1e-3, 3e-4, 1e-4)):
    """walk-only extrapolation Z = Z_inf + c1 a^(2 gamma*) + c2 a (3 points, exact solve)"""
    zs = PR[v]["zeta_roots"]
    kappa = PR[v]["kappa"]
    g2 = 2 * np.sqrt(kappa**2 - zs[-1]**2) if zs else None
    Z = np.array([walk(v, a) for a in alist])
    if g2 is None or any(z is None for z in Z):
        return None
    A = np.array([[1.0, a**g2, a] for a in alist])
    return float(np.linalg.solve(A, Z)[0])


ex = {}
# F1 retrospective
zm, _ = predict_Zcr(3e-5, -1, 0.5, -1.0, zlo=0.85, zhi=1.25, n=60)
ex["F1"] = dict(Z_walk=walk("base", 3e-5), Z_match=zm / ALPHA, diff=walk("base", 3e-5) - zm / ALPHA,
                PASS=abs(walk("base", 3e-5) - zm / ALPHA) <= 0.01, note="retrospective (seen)")
ex["F1"]["table"] = {f"{a:.6g}": dict(walk=walk("base", a), match=match("base", a) if f"{a:.6g}" in PR["base"]["Zcr_matching"] else None)
                     for a in (3e-3, 1e-3, 3e-4, 1e-4)}
# F2
f2 = {}
for v in ("R+1", "R+i", "R-i", "S0.25", "K+1", "K+1R+1"):
    d4 = walk(v, 1e-4) - match(v, 1e-4) if walk(v, 1e-4) is not None else None
    d3 = walk(v, 1e-3) - match(v, 1e-3)
    f2[v] = dict(Z_star=PR[v]["Z_roots"], walk={f"{a:.6g}": walk(v, a) for a in PRE["configs"]["a"]},
                 match={f"{a:.6g}": match(v, a) for a in PRE["configs"]["a"]}, d_1em4=d4, d_1em3=d3,
                 PASS=(d4 is not None and abs(d4) <= 0.05 and abs(d3) <= 0.2))
ex["F2"] = dict(per_variant=f2, PASS=all(x["PASS"] for x in f2.values()))
# F3
f3 = {}
for v in ("R+1", "S0.25", "K+1"):
    zi = zinf(v)
    f3[v] = dict(Z_inf=zi, Z_star=PR[v]["Z_roots"][-1], diff=None if zi is None else zi - PR[v]["Z_roots"][-1],
                 PASS=zi is not None and abs(zi - PR[v]["Z_roots"][-1]) <= 0.1)
f3["base (giro 1 data, same fit)"] = dict(Z_inf=zinf("base"), Z_star=PR["base"]["Z_roots"][-1])
ex["F3"] = dict(per_variant=f3, PASS=all(x["PASS"] for k, x in f3.items() if k in ("R+1", "S0.25", "K+1")))
# F4 duality
zb, zd = zinf("base"), zinf("K+1R+1")
zr, zk = zinf("R+1"), zinf("K+1")
ex["F4"] = dict(base=zb, K1R1=zd, R1=zr, K1=zk,
                PASS=None not in (zb, zd, zr, zk) and abs(zb - zd) <= 0.1 and abs(zr - zk) <= 0.1)
# U1
zst = {v: PR[v]["Z_roots"][-1] for v in ("base", "R+1", "R+i", "R-i", "S0.25")}
spread = max(zst.values()) - min(zst.values())
ex["U1"] = dict(Z_star=zst, spread=spread, universal=spread < 1, PASS_prediction_not_universal=spread >= 1)
# R1 remedy
r1 = {}
for v in ("S1", "S0.866"):
    Zw = [walk(v, a) for a in PRE["configs"]["a"]]
    ok = all(z is not None and z >= 137.036 for z in Zw) and all(Zw[i] > Zw[i + 1] for i in range(len(Zw) - 1))
    d4 = Zw[-1] - match(v, 1e-4) if Zw[-1] is not None else None
    r1[v] = dict(walk=Zw, match=[match(v, a) for a in PRE["configs"]["a"]], match_1em6=match(v, 1e-6),
                 match_1em9=match(v, 1e-9), d_1em4=d4, PASS=ok and d4 is not None and abs(d4) <= 0.2)
r1["S0.75 (not an exam)"] = dict(walk=[walk("S0.75", a) for a in PRE["configs"]["a"]],
                                 match=[match("S0.75", a) for a in PRE["configs"]["a"]])
ex["R1"] = dict(per_variant=r1, PASS=all(r1[v]["PASS"] for v in ("S1", "S0.866")))
# K2
zk2 = [walk("K-2", a) for a in PRE["configs"]["a"]]
ex["K2"] = dict(walk=zk2, Z_star=PR["K-2"]["Z_roots"], PASS=zk2[-1] is not None and abs(zk2[-1] - 206.0207) <= 0.3)
# K3 (A3: crossing of omega = 0 and omega = -1 by the contact state, sign method)
AK = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
zx0 = {f"{a:.6g}": SGN.get(f"K+2|{a:.6g}|omega=0", {}).get("Zcr") for a in AK}
zx1 = {f"{a:.6g}": SGN.get(f"K+2|{a:.6g}", {}).get("Zcr") for a in AK}
zx_pre = {f"{a:.6g}": RUN.get(f"K+2|{a:.6g}", {}).get("Zx") for a in (1e-3, 3e-4, 1e-4)}
ex["K3"] = dict(Z_x_omega0=zx0, Z_x_omega_minus1=zx1, tracker_grid_prereg=zx_pre, Z_star=PR["K+2"]["Z_roots"],
                PASS=zx0["0.0001"] is not None and abs(zx0["0.0001"] - 42.0531) <= 0.1)
# O1-O3
ex["O1"] = dict(**RUN["order"], PASS=RUN["order"]["max_dlambda"] < 1e-10)
ds = walk("split", 1e-4) - walk("base", 1e-4) if walk("split", 1e-4) else None
ex["O2"] = dict(split={f"{a:.6g}": walk("split", a) for a in PRE["configs"]["a"]},
                base={f"{a:.6g}": walk("base", a) for a in PRE["configs"]["a"]}, d_1em4=ds,
                PASS=ds is not None and abs(ds) <= 0.05)
dr = walk("refl1", 1e-3) - walk("base", 1e-3) if walk("refl1", 1e-3) else None
ex["O3"] = dict(refl1=walk("refl1", 1e-3), base=walk("base", 1e-3), diff=dr, PASS=dr is not None and abs(dr) <= 1e-3)
ex["A3_validation"] = {v: {f"{a:.6g}": dict(sign=walk(v, a), tracker=tracker(v, a)) for a in PRE["configs"]["a"]}
                       for v in ("base", "R+1", "K+1", "K+1R+1", "R-i", "R+i", "S0.25", "S1", "S0.866", "S0.75", "split", "K-2")}
ex["fine_a"] = {v: {f"{a:.6g}": walk(v, a) for a in (1e-4, 3e-5, 1e-5)} for v in
                ("base", "R+1", "K+1", "K+1R+1", "R-i", "R+i", "S0.25", "S1", "S0.866", "S0.75", "split", "K-2")}
ex["fine_a_Zinf"] = {v: zinf(v, (1e-4, 3e-5, 1e-5)) for v in ("base", "R+1", "K+1", "K+1R+1", "R-i", "R+i", "S0.25")}
json.dump(ex, open(os.path.join(HERE, "out", "zstar_exams.json"), "w"), indent=1, default=float)
for k, v in ex.items():
    print(k, "PASS" if v.get("PASS", v.get("PASS_prediction_not_universal")) else "FAIL",
          {kk: vv for kk, vv in v.items() if kk not in ("per_variant", "table")} if k not in ("F2", "R1", "F3") else "")
print("A3 validation", json.dumps(ex["A3_validation"], indent=0)[:3000])
print("fine a", ex["fine_a"], ex["fine_a_Zinf"])
for k in ("F2", "F3", "R1"):
    for v, x in ex[k]["per_variant"].items():
        print("  ", k, v, {kk: vv for kk, vv in x.items() if kk not in ("walk", "match")})

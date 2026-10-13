"""nucleo1s_f5_judge.py - judges phase 5 (exams S1-S4, J1-J3 of lab/out/nucleo1s_f5_prereg.json, amendment
nucleo1s_f5_amend1.json) -> lab/out/nucleo1s_f5_judged.json"""
import os, json
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = json.load(open(os.path.join(OUT, "nucleo1s_f5_prereg.json")))
R = json.load(open(os.path.join(OUT, "nucleo1s_f5_results.json")))
R3 = json.load(open(os.path.join(OUT, "nucleo1s_3d_results.json")))
JN = json.load(open(os.path.join(OUT, "nucleo1s_junction3d.json")))
E = P["exams"]; pred = P["predictions"]
J = {}


def spread(L, lo, hi, wmin=0.05):
    s = [x["omega"] for x in L if x["weight"] >= wmin and lo <= x["omega"] <= hi]
    return (max(s) - min(s)) if s else None


learned = {"60": dict(s=R3["learn2"]["s"], seen=True)}
for Z in (20, 40, 80):
    d = R.get(f"learn_Z{Z}", {})
    if d.get("done"):
        learned[str(Z)] = dict(s=d["s"], std=d["std"], dU=d["dU"], U=d["U_bare"], n_evals=len(d["evals"]))
tab = {}
for Z in ("20", "40", "60", "80"):
    if Z in learned:
        tab[Z] = dict(learned=learned[Z]["s"], PA=pred["PA"][f"Z{Z}_a0.5"], PB1=pred["PB1"][f"Z{Z}_a0.5"],
                      PB2=pred["PB2"][f"Z{Z}_a0.5"])
test = [z for z in ("20", "40", "80")]
ok1 = all(z in tab and tab[z]["PB1"][0] is not None and abs(tab[z]["PB1"][0] - tab[z]["learned"][0]) <= E["S1"]["tol"] for z in test)
J["S1"] = dict(table=tab, PASS=bool(ok1))
ok2 = all(z in tab and abs(tab[z]["PA"][1] - tab[z]["learned"][1]) <= E["S2"]["tol"] for z in test)
J["S2"] = dict(PASS=bool(ok2), s1={z: (tab[z]["learned"][1] if z in tab else None) for z in test}, PA_s1=P["P_A"]["s1"])
u3 = {}
for a in (0.4, 0.24, 0.2):
    r = R.get(f"use_PB1_Z60_{a}")
    if r:
        u3[str(a)] = dict(s=r["s"], rel=r["rel"], spread=spread(r["lines"], 0.798, 0.950))
J["S3"] = dict(use=u3, PASS=bool(len(u3) == 3 and all(abs(v["rel"]) <= E["S3"]["tol"] and v["spread"] is not None
                                                   and v["spread"] <= E["S3"]["spread"] for v in u3.values())))
u4 = {}
from zz_core import ALPHA, dirac_formula  # noqa: E402
for Z in (20, 40, 80):
    r = R.get(f"use_PB1_Z{Z}_0.5")
    Bd = 1 - dirac_formula(1, -1, Z * ALPHA)
    u4[str(Z)] = (dict(s=r["s"], rel=r["rel"], spread=spread(r["lines"], 1 - 2 * Bd, 1 - 0.5 * Bd)) if r else
                  dict(note="no P_B1 prediction exists (no root)"))
J["S4"] = dict(use=u4, PASS=bool(all("spread" in v and v["spread"] is not None and v["spread"] <= E["S4"]["spread"]
                                     for v in u4.values())))
pre = [json.dumps(c) for c in P["config"]["junction"]["cases"]]
rows = {json.dumps(r["case"]): r for r in JN}
J["J1"] = dict(max_drift=max(r["norm_drift"] for r in JN), PASS=bool(all(r["norm_drift"] <= E["J1"]["tol"] for r in JN)))
nrm = [c for c in P["config"]["junction"]["cases"] if c["a"] == 0.25 and c["k"][1] == 0 and "thj" not in c]
obl = [c for c in P["config"]["junction"]["cases"] if c["k"][1] != 0]
J["J2"] = dict(R={str(c["k"][0]): rows[json.dumps(c)]["R_conservative"] for c in nrm if json.dumps(c) in rows},
               PASS=bool(all(json.dumps(c) in rows and rows[json.dumps(c)]["R_conservative"] <= E["J2"]["tol"] for c in nrm)))
J["J3"] = dict(R={str(c["k"]): rows[json.dumps(c)]["R_conservative"] for c in obl if json.dumps(c) in rows},
               PASS=bool(all(json.dumps(c) in rows and rows[json.dumps(c)]["R_conservative"] <= E["J3"]["tol"] for c in obl)))
J["junction_all"] = [dict(case=r["case"], R=r["R_conservative"], drift=r["norm_drift"]) for r in JN]
json.dump(J, open(os.path.join(OUT, "nucleo1s_f5_judged.json"), "w"), indent=1, default=float)
for k, v in J.items():
    if isinstance(v, dict):
        print(k, v.get("PASS"), {kk: vv for kk, vv in v.items() if kk != "PASS"})

# ---- amendment 2 (nucleo1s_f5_amend2.json): P_C, the scale-free first-cell prediction
A2 = json.load(open(os.path.join(OUT, "nucleo1s_f5_amend2.json")))
PC = json.load(open(os.path.join(OUT, "nucleo1s_f5_pc.json")))
s5 = {}
for Z in ("20", "40", "60", "80"):
    p1 = PC[f"PC1_Z{Z}_a0.01"]["s"]; p2 = PC[f"PC2_Z{Z}_a0.01"]["s"]
    ln = ls_ = learned[Z]["s"][0]
    s5[Z] = dict(PC1=p1, PC1_a002=PC[f"PC1_Z{Z}_a0.02"]["s"], PC2=p2, learned=ln,
                 diff=(p1 - ln) if p1 is not None else None)
J["S5"] = dict(table=s5, note="comparison with values already seen (declared)",
               PASS=bool(all(v["diff"] is not None and abs(v["diff"]) <= A2["exams"]["S5"]["tol"] for v in s5.values())))
s6 = {}
for Z, a in A2["exams"]["S6"]["points"]:
    r = R.get(f"pc_Z{Z}_{a}")
    Bd = 1 - dirac_formula(1, -1, Z * ALPHA)
    if r:
        s6[f"Z{Z}_a{a}"] = dict(s=r["s"], rel=r["rel"], spread=spread(r["lines"], 1 - 2 * Bd, 1 - 0.5 * Bd))
ok6 = len(s6) == 4 and all(v["spread"] is not None and v["spread"] <= A2["exams"]["S6"]["spread"] for v in s6.values()) \
    and all(abs(v["rel"]) <= A2["exams"]["S6"]["tol_Z60"] for k, v in s6.items() if k.startswith("Z60"))
J["S6"] = dict(use=s6, PASS=bool(ok6))
json.dump(J, open(os.path.join(OUT, "nucleo1s_f5_judged.json"), "w"), indent=1, default=float)
for k in ("S5", "S6"):
    print(k, J[k]["PASS"], {kk: vv for kk, vv in J[k].items() if kk != "PASS"})

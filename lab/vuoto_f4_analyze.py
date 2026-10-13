"""phase 4: exams against lab/out/vuoto_f4_prereg.json -> lab/out/vuoto_f4_exams.json (+ Ward rerun G1)"""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vuoto_f4_core as F4
from zz_core import ALPHA

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
g2 = json.load(open(os.path.join(D, "vuoto_f4_g2.json")))
li = json.load(open(os.path.join(D, "vuoto_f4_light.json")))
CASES = [("R", "K12"), ("R", "K12/3"), ("C", "K12/3")]
TH = [0.05, 0.025, 0.0125, 0.00625]
KM = [1e-3, 3e-4]
PAIRS = ["20", "21", "01", "02", "10", "12"]


def lin0(x, y):
    c = np.polyfit(x, y, 1)
    return float(c[1])


E = {"G": {}, "L": {}, "K": {}}
# ------------------------------------------------ g - 2
for rule, cn in CASES:
    tab = {}
    for th in TH:
        runs = [g2["runs"]["%s|%s|%g|%g" % (rule, cn, th, k)] for k in KM]
        a = np.array([[r[p]["a_over_a2pi"] for p in PAIRS] for r in runs])          # (kmin, pair)
        a0 = a[0] - (a[0] - a[1]) / (KM[0] - KM[1]) * KM[0]                          # kmin -> 0 per pair
        tab[th] = {"pairs_kmin0": dict(zip(PAIRS, a0.tolist())), "mean": float(a0.mean()),
                   "spread": float(a0.max() - a0.min()), "mean_kmin": a.mean(1).tolist(),
                   "max_imag": max(abs(r[p]["imag_r"]) for r in runs for p in PAIRS),
                   "max_resid": max(r[p]["resid"] for r in runs for p in PAIRS)}
    th = np.array(TH)
    mean0 = lin0(th, [tab[t]["mean"] for t in TH])
    pairs0 = {p: lin0(th, [tab[t]["pairs_kmin0"][p] for t in TH]) for p in PAIRS}
    spread0 = max(pairs0.values()) - min(pairs0.values())
    sp = np.array([tab[t]["spread"] for t in TH])
    pw = float(np.polyfit(np.log(th), np.log(np.maximum(sp, 1e-300)), 1)[0])
    key = "%s %s" % (rule, cn)
    E["G"][key] = {"per_theta": tab, "mean_a0": mean0, "pairs_a0": pairs0, "spread_a0": spread0,
                   "spread_power_of_theta": pw,
                   "G2": "PASS" if abs(mean0 - 1) <= 0.03 else "FAIL",
                   "G3": "survives" if spread0 >= 0.01 else "vanishes",
                   "G5": "PASS" if max(max(tab[t]["max_imag"], tab[t]["max_resid"]) for t in TH) <= 1e-8 else "FAIL"}
    if rule == "C":
        E["G"][key]["G4"] = "PASS" if max(sp) <= 1e-6 else "FAIL"
    # electron speed (kmin 1e-3 run)
    dv = np.array([g2["runs"]["%s|%s|%g|%g" % (rule, cn, t, KM[0])]["dmv"]["dv"] for t in TH]) / ALPHA
    E["K"].setdefault(key, {})["electron_dv_over_alpha"] = dv.tolist()
    E["K"][key]["electron_iso_a0"] = lin0(th, dv.mean(1))
    E["K"][key]["electron_aniso_a0"] = lin0(th, dv[:, 1] - 0.5 * (dv[:, 0] + dv[:, 2]))
# ------------------------------------------------ light
for rule, cn in CASES:
    key = "%s %s" % (rule, cn)
    per = {}
    for t in TH:
        row = {}
        for an in ("x", "y", "z", "d"):
            r = li["runs"]["%s|%s|%g|%s" % (rule, cn, t, an)]
            row[an] = {"dc": r["dc_over_alpha"], "mean": r["dc_mean"], "biref": r["biref"],
                       "imag": max(r["imag_rel_over_alpha"]), "mass2": r["mass2_over_alpha"]}
        per[t] = row
    th = np.array(TH)
    iso = [np.mean([per[t][a]["mean"] for a in "xyz"]) for t in TH]
    ani = [per[t]["y"]["mean"] - 0.5 * (per[t]["x"]["mean"] + per[t]["z"]["mean"]) for t in TH]
    bire0 = {a: lin0(th, [abs(per[t][a]["biref"]) for t in TH]) for a in "xyzd"}
    L = {"per_theta": per, "iso_a0": lin0(th, iso), "aniso_a0": lin0(th, ani), "biref_a0": bire0,
         "L4": "PASS" if max(per[t][a]["imag"] for t in TH for a in "xyzd") <= 0.01 else "FAIL",
         "L3": {a: ("seen" if abs(bire0[a]) >= 0.002 else "not seen") for a in "xyzd"}}
    if rule == "R":
        dxz = max(abs(np.sort(per[t]["x"]["dc"]) - np.sort(per[t]["z"]["dc"])).max() for t in TH)
        L["L1"] = {"max_|x-z|": float(dxz), "result": "PASS" if dxz <= 0.002 else "FAIL"}
    else:
        spr = max(max(per[t][a]["mean"] for a in "xyzd") - min(per[t][a]["mean"] for a in "xyzd") for t in TH)
        L["L1"] = {"max_spread": float(spr), "result": "PASS" if spr <= 0.002 else "FAIL"}
    E["L"][key] = L
    K = E["K"][key]
    K["light_iso_a0"] = L["iso_a0"]; K["light_aniso_a0"] = L["aniso_a0"]
    K["Delta_iso_a0"] = K["electron_iso_a0"] - L["iso_a0"]
    K["Delta_aniso_a0"] = K["electron_aniso_a0"] - L["aniso_a0"]
    K["Delta_iso_per_theta"] = [float(np.mean(K["electron_dv_over_alpha"][i]) - iso[i]) for i in range(len(TH))]
    K["verdict"] = "physical difference" if abs(K["Delta_iso_a0"]) >= 0.01 else "below 0.01 alpha"
# systematics L2
base = li["runs"]
syst = {}
for tag in ("grid2", "q2", "wsym"):
    d = []
    for an in "xyz":
        b = np.sort(base["R|K12|0.05|%s" % an]["dc_over_alpha"]); s = np.sort(base["R|K12|0.05|%s|%s" % (an, tag)]["dc_over_alpha"])
        d.append(float(np.abs(b - s).max()))
    syst[tag] = max(d)
E["L"]["L2"] = {"max_change": syst, "result": "PASS" if max(syst.values()) <= 0.005 else "FAIL"}
# ------------------------------------------------ G1 Ward rerun at the production grid
W = {}
for rule, K12, j in (("R", (0.6, 1.05), 0), ("R", (0.6, 1.05), 1), ("C", (0.2, 0.35), 0)):
    th = 0.05; g = dict(nr=64, nct=12, nph=16, kmin_fac=1e-3, K12=K12)
    p0 = np.zeros(3); p0[j] = 0.3 * th
    r = F4.vertex_all(th, rule, p0, p0, j, grid=g)
    h = 1e-4 * th; dd = []
    for s in (1, -1):
        pp = p0.copy(); pp[j] += s * h
        q = F4.dmu_onshell(th, rule, pp, g); dd.append((q["Om"][2:], q["dmu"][2:]))
    dmu0 = (np.exp(-1j * dd[0][0]) - np.exp(-1j * dd[1][0])).sum() / (2 * h)
    ward = 1 + (dd[0][1] - dd[1][1]).sum() / (2 * h) / dmu0
    v = np.trace(r["V"]) / np.trace(r["F1"])
    W["%s axis %d" % (rule, j)] = {"V/F1": [v.real, v.imag], "Ward": [ward.real, ward.imag], "diff": float(abs(v - ward))}
E["G"]["G1"] = {"runs": W, "result": "PASS" if max(w["diff"] for w in W.values()) <= 1e-6 else "FAIL"}
json.dump(E, open(os.path.join(D, "vuoto_f4_exams.json"), "w"), indent=1)
for k in E["K"]:
    print(k, {x: (round(v, 4) if isinstance(v, float) else v) for x, v in E["K"][k].items() if "a0" in x or x == "verdict"})
for k in E["G"]:
    if k != "G1":
        print(k, round(E["G"][k]["mean_a0"], 4), round(E["G"][k]["spread_a0"], 4), E["G"][k]["G2"], E["G"][k]["G3"], E["G"][k].get("G4"))
print("G1", E["G"]["G1"]["result"], "L2", E["L"]["L2"])
for k in E["L"]:
    if k != "L2":
        print(k, E["L"][k]["L1"], E["L"][k]["L3"], E["L"][k]["L4"], {a: round(v, 5) for a, v in E["L"][k]["biref_a0"].items()})

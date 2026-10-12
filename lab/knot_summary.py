"""knot_summary.py - K1 summary (Richardson, box trend, ball-ODE curve) and pass/fail of all exams -> knot_results.json"""
import os, sys, json, glob, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
J = lambda n: json.load(open(os.path.join(OUT, n)))
ex = lambda n: os.path.exists(os.path.join(OUT, n))
R = {}
k1 = {t: J(f"knot_k1_{t}.json") for t in ["h0.2_L32", "h0.2_L48", "h0.2_L64", "h0.15_L32", "h0.15_L48", "h0.15_L64",
                                          "h0.2_L48_crude", "h0.25_L38"] if ex(f"knot_k1_{t}.json")}
e = {t: d["E_over_12pi2"] for t, d in k1.items()}
cont = e["h0.15_L64"] + (e["h0.15_L64"] - e["h0.2_L48"]) * 0.15**2 / (0.2**2 - 0.15**2)
ballR = [2.4, 3.2, 3.6, 4.8, 6.4, 8.0, 10.0]
ballE = [kc.ball_energy(r) / (12 * math.pi**2) for r in ballR]
S = dict(continuum=cont, ball_R=ballR, ball_E=ballE, E=e, B={t: d["B"] for t, d in k1.items()},
         Lambda={t: d["Lambda"] for t, d in k1.items()}, iterations={t: d["iterations"] for t, d in k1.items()},
         Fmax={t: d["Fmax"] for t, d in k1.items()})
d1 = e["h0.2_L32"] - e["h0.2_L48"]; d2 = e["h0.2_L48"] - e["h0.2_L64"]
S["diff_ratio_h0.2"] = d1 / d2
if "h0.15_L32" in e:
    S["diff_ratio_h0.15"] = (e["h0.15_L32"] - e["h0.15_L48"]) / (e["h0.15_L48"] - e["h0.15_L64"])
json.dump(S, open(os.path.join(OUT, "knot_k1_summary.json"), "w"), indent=1)
R["K1a_continuum"] = dict(value=cont, dev=cont / 1.2315 - 1, pass_=abs(cont / 1.2315 - 1) < 0.01)
R["K1b_box"] = dict(value=abs(e["h0.2_L64"] - e["h0.2_L48"]) / e["h0.2_L64"], pass_=abs(e["h0.2_L64"] - e["h0.2_L48"]) / e["h0.2_L64"] < 0.005)
R["K1c_power"] = dict(value=d1 / d2, pass_=d1 / d2 >= 2.5)
prod = [t for t in k1 if t != "h0.25_L38"]
R["K1d_topology"] = dict(Bmaxdev=max(abs(k1[t]["B"] - 1) for t in prod),
                         crude_dev=abs(e["h0.2_L48_crude"] / e["h0.2_L48"] - 1) if "h0.2_L48_crude" in e else None)
R["K1d_topology"]["pass_"] = R["K1d_topology"]["Bmaxdev"] < 0.01 and (R["K1d_topology"]["crude_dev"] or 1) < 0.002
# K2
for f in sorted(glob.glob(os.path.join(OUT, "knot_k2_*.json"))):
    d = json.load(open(f)); v = d["meta"]["v"]; key = os.path.basename(f)[8:-5]
    rest = os.path.join(OUT, f"knot_k2_h{d['meta']['h']:g}_v0.json")
    r0 = json.load(open(rest))["ratio_b_late"] if os.path.exists(rest) else 1.0
    g = d["meta"]["gamma"]
    item = dict(E_over_gE0=d["E_over_gE0"], contraction=d["ratio_b_late"] / r0 * g, rest_ratio=r0,
                shape_kept=d["e_sx_end_over_t0"], v_fit=d["v_fit"], drift=d["drift"], B_range=[d["B_min"], d["B_max"]],
                r_half_contracted_over_h=d["meta"]["r_half_contracted_over_h"], e_ratio_late_times_gamma=d["ratio_e_late"] * g)
    if v > 0:
        item["pass_"] = dict(K2a=abs(item["E_over_gE0"] - 1) < 0.03, K2b=abs(item["contraction"] - 1) < 0.03,
                             K2c=abs(item["shape_kept"] - 1) < 0.03, K2d=abs(item["v_fit"] / v - 1) < 0.03,
                             valid=item["drift"] < 0.005)
    else:
        item["pass_"] = dict(rest_ratio_1=abs(r0 - 1) < 0.01, valid=item["drift"] < 0.005)
    R["K2_" + key] = item
# K3
for f in sorted(glob.glob(os.path.join(OUT, "knot_k3_*.json"))):
    d = json.load(open(f)); key = os.path.basename(f)[8:-5]
    r = d["records"]
    item = dict(frac_out=d["frac_out"], absB_min=d["absB_min"], absB_final=d["absB_final"], t_absB_below_0p2=d["t_absB_below_0p2"],
                B_final=d["B_final"], E0=d["E0"], E_end=r[-1]["E"], I_end=[r[-1]["Ixx"], r[-1]["Iyy"], r[-1]["Izz"]],
                pos_end=[r[-1]["bpos_x"], r[-1]["bpos_y"], r[-1]["bpos_z"]])
    if d["kind"] == "SA":
        item["pass_"] = dict(K3a=d["absB_min"] < 0.2, K3b=d["frac_out"] >= 0.7)
    elif d["kind"] == "SS":
        item["pass_"] = dict(B2=abs(d["B_final"] / 2 - 1) < 0.02)
    else:
        item["pass_"] = dict(control=abs(d["B_final"] - 1) < 0.01)
    R["K3_" + key] = item
if ex("knot_k4.json"):
    R["K4"] = J("knot_k4.json")["pass"]
if ex("knot_k5.json"):
    R["K5"] = J("knot_k5.json")["pass"]
print(json.dumps(R, indent=1))
json.dump(R, open(os.path.join(OUT, "knot_results.json"), "w"), indent=1)

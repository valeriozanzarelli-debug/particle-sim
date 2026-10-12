"""knot2_summary.py - pass/fail of S1-S4 against knot2_prereg.json -> lab/out/knot2_results.json"""
import os, sys, json, glob, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
J = lambda n: json.load(open(os.path.join(OUT, n)))
ex = lambda n: os.path.exists(os.path.join(OUT, n))
R = {}
if ex("knot2_s1.json"):
    s = J("knot2_s1.json")
    R["S1"] = {k: {a: b for a, b in v.items() if a.startswith("max") or a == "code_check_pass"} for k, v in s.items() if k.startswith("S1a")}
    R["S1"]["S1b"] = s["S1b_old_rest_knot"]; R["S1"]["S1c"] = s["S1c_cherenkov_kh"]; R["S1"]["S1d"] = s["S1d_knot_spectrum_fraction"]
# S2
rest = {}
for f in glob.glob(os.path.join(OUT, "knot2_boost_*_v0.json")):
    d = json.load(open(f)); rest[(d["meta"]["scheme"], d["meta"]["h"])] = d["ratio_b_late"]
for f in sorted(glob.glob(os.path.join(OUT, "knot2_boost_*.json"))):
    d = json.load(open(f)); m = d["meta"]; v = m["v"]; key = os.path.basename(f)[12:-5]
    r0 = rest.get((m["scheme"], m["h"]), rest.get((m["scheme"], 0.2), 1.0))
    g = m["gamma"]
    it = dict(E_over_gE0=d["E_over_gE0"], contraction=d["ratio_b_late"] / r0 * g, rest_ratio=r0, shape_kept=d["e_sx_end_over_t0"],
              v_fit=d["v_fit"], drift=d["drift"], B_range=[d["B_min"], d["B_max"]], balls=m["r_half_contracted_over_h"])
    if v > 0:
        it["pass_"] = dict(E=abs(it["E_over_gE0"] - 1) < 0.03, contraction=abs(it["contraction"] - 1) < 0.03,
                           speed=abs(it["v_fit"] / v - 1) < 0.03, shape=abs(it["shape_kept"] - 1) < 0.05, valid=it["drift"] < 0.005)
        it["all"] = all(it["pass_"].values())
    else:
        it["pass_"] = dict(rest_ratio_1=abs(r0 - 1) < 0.01, valid=it["drift"] < 0.005)
    R["S2_" + key] = it
ok = [R[k]["balls"] for k in R if k.startswith("S2_c4") and R[k].get("all")]
R["S2_min_balls_passing_c4"] = min(ok) if ok else None
if ex("knot2_s3.json"):
    s = J("knot2_s3.json")
    R["S3"] = dict(primary=s["primary"], old_richardson=s["secondary_old_richardson"], fit=s["secondary_fit_c4_h0.2"],
                   mu_ratio=s["mu"]["lattice_primary"]["mu_ratio"], mu_pass=s["mu"]["pass_"], mu_dev_exp=s["mu"]["dev_exp"],
                   mu_dev_ode=s["mu"]["dev_ode"])
for f in sorted(glob.glob(os.path.join(OUT, "knot2_s4_v*.json"))):
    d = json.load(open(f)); key = os.path.basename(f)[9:-5]
    R["S4_" + key] = {k: d[k] for k in ("frac_out_24", "frac_out_T", "t_cross_70", "absB_min", "t_absB_below_0p2", "B_final", "E0", "T")}
    if "spectral" in d:
        R["S4_" + key]["spectral"] = {a: b for a, b in d["spectral"].items() if not isinstance(b, list) or a == "times"}
    if "front_posthoc" in d:
        R["S4_" + key]["front_speed_posthoc"] = d["front_posthoc"]["speed"]
    R["S4_" + key]["pass_"] = dict(radiation=(d["frac_out_24"] or 0) >= 0.7, unknot=d["absB_min"] < 0.2,
                                   massless=bool(d.get("spectral", {}).get("pass_", False)))
print(json.dumps(R, indent=1))
R = json.loads(json.dumps(R, default=float))
json.dump(R, open(os.path.join(OUT, "knot2_results.json"), "w"), indent=1)

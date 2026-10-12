"""zoomknot_summary.py - verdicts Z0-Z4 from the result files; writes lab/out/zoomknot_results.json and prints the
error-vs-cost / error-budget table."""
import os, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
PRE = json.load(open(os.path.join(OUT, "zoomknot_prereg.json")))
ld = lambda t: json.load(open(os.path.join(OUT, f"zoomknot_{t}.json"))) if os.path.exists(os.path.join(OUT, f"zoomknot_{t}.json")) else None
R = {t: ld(t) for t in ["U0", "Z1_h0.05_R0.6", "Z1_h0.05_R1.0", "Z1_h0.05_R1.4", "Z1_h0.0667_R1.0", "Z1_h0.0667_R1.4",
                        "Z2", "Z2_learned"]}
L = json.load(open(os.path.join(OUT, "zoomknot_learn.json")))
BUD = 64 ** 3
rows = []
for t, m in R.items():
    if m is None:
        continue
    lv = m["levels"]
    fine = (2 * lv[0]["N"] + 1) ** 3 if len(lv) > 1 else 0
    rows.append(dict(tag=t, sites=m["sites"], fine_sites=fine, agglo_sites=m["sites"] - fine,
                     E12=m["E_over_12pi2"], err_E=m["err_E"], p0=m["p0_central"], err_p0=m["err_p0"],
                     err_p0_fb=m["err_p0_fb"], core_dev=m["core_max_dev_central"], core_dev_fb=m["core_max_dev_fb"],
                     Lam=m["Lam"], err_Lam=m["err_Lam"], virial=m["virial"], vonLaue_fb=m["vonLaue_fb"],
                     vonLaue_c=m["vonLaue_central"], B=m["B"], r0=m["r0_lattice"], stop=m["stop"],
                     steps=m["outer_steps"], wall=m["wall_s"]))
V = {}
u = R["U0"]
if u:
    V["Z0"] = dict(E12=u["E_over_12pi2"], rel_vs_P4=u["E_over_12pi2"] / 1.2057818 - 1, core_dev=u["core_max_dev_central"],
                   pass_=bool(abs(u["E_over_12pi2"] / 1.2057818 - 1) < 2e-3 and u["core_max_dev_central"] > 0.30))
z1 = {t: m for t, m in R.items() if t.startswith("Z1") and m}
if z1:
    ok = [t for t, m in z1.items() if m["sites"] <= BUD and abs(m["err_p0"]) < 0.05 and m["core_max_dev_central"] < 0.05]
    V["Z1a"] = dict(passing=ok, pass_=bool(ok))
    seq = [R[f"Z1_h0.05_R{x}"] for x in ("0.6", "1.0", "1.4")]
    if all(seq):
        d = [m["core_max_dev_central"] for m in seq]
        V["Z1b"] = dict(core_dev=d, pass_=bool(d[0] > d[1] > d[2]))
    inb = {t: m for t, m in z1.items() if m["sites"] <= BUD}
    best = min(inb, key=lambda t: inb[t]["core_max_dev_central"])
    V["Z1c"] = dict(best=best, err_E=inb[best]["err_E"], pass_=bool(abs(inb[best]["err_E"]) < 3e-3))
z = R["Z2_learned"]
if z:
    V["Z2a"] = dict(err_p0=z["err_p0"], core_dev=z["core_max_dev_central"],
                    pass_=bool(abs(z["err_p0"]) < 0.05 and z["core_max_dev_central"] < 0.05))
    V["Z2b"] = dict(err_E=z["err_E"], pass_2pct=bool(abs(z["err_E"]) < 0.02), pass_0p3pct=bool(abs(z["err_E"]) < 3e-3))
    V["Z2c"] = dict(err_Lam=z["err_Lam"], pass_=bool(abs(z["err_Lam"]) < 0.03))
    zb = R["Z2"]
    V["Z2d"] = dict(levels_learned_better=L["n_levels_learned_better"], i_pass=L["Z2d_i_pass"],
                    err_E_learned=z["err_E"], err_E_bare=zb["err_E"] if zb else None,
                    ii_pass=bool(zb and abs(z["err_E"]) < abs(zb["err_E"])))
    V["Z2d"]["pass_"] = bool(V["Z2d"]["i_pass"] and V["Z2d"]["ii_pass"])
V["Z2e"] = dict(sum_static_err_learned=L["sum_static_err_learned_over_E"], sum_static_err_bare=L["sum_static_err_bare_over_E"],
                pass_=L["Z2e_pass"])
V["Z4"] = dict(status="SKIPPED (budget: 4 shared cores, 1 thread)")
out = dict(rows=rows, verdicts=V, learned=dict(k2=L["k2"], k4=L["k4"]))
json.dump(out, open(os.path.join(OUT, "zoomknot_results.json"), "w"), indent=1)
print(f"{'config':18s} {'sites':>7s} {'fine':>7s} {'agglo':>7s} {'E/12pi2':>9s} {'dE':>7s} {'p0':>6s} {'dp0':>7s} "
      f"{'core':>6s} {'Lam':>6s} {'dLam':>7s} {'vir':>7s} {'vL_fb':>7s} {'r0':>6s} stop")
for r in rows:
    print(f"{r['tag']:18s} {r['sites']:7d} {r['fine_sites']:7d} {r['agglo_sites']:7d} {r['E12']:9.6f} {r['err_E']:+7.2%} "
          f"{r['p0']:6.2f} {r['err_p0']:+7.2%} {r['core_dev']:6.1%} {r['Lam']:6.1f} {r['err_Lam']:+7.2%} "
          f"{r['virial']:+7.2%} {r['vonLaue_fb']:+7.3f} {r['r0']:6.3f} {r['stop']}({r['steps']}, {r['wall']:.0f}s)")
print(json.dumps(V, indent=1))

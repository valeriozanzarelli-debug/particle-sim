"""surf_figure.py -> lab/out/surf_figure.png and lab/out/surf_results.json (summary PASS/FAIL against the prereg)"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
L = lambda n: json.load(open(os.path.join(OUT, f"surf_{n}.json")))
pre = L("prereg")["exams"]
f6, cont, h1s, beat, orb = L("f6"), L("cont"), L("h1s"), L("beat"), L("orb")
a1, a12, afr = L("a3d_1s"), L("a3d_1s2p"), L("a3d_free")
# AMENDMENT (surf_amend.json): TV ratio against an ideal sampler of the FINAL density (surf_tvnull.py)
if os.path.exists(os.path.join(OUT, "surf_tvnull.json")):
    tvn = L("tvnull")
    for nm, d in (("1s", a1), ("1s2p", a12), ("free", afr)):
        d["tv_ratio_final_t0ref"] = d["tv_ratio_final"]
        d["tv_ratio_final"] = tvn[nm]["tv_ratio_final"]

R = {}
def ok(name, val, cond):
    R[name] = dict(value=val, PASS=bool(cond))

# E0
for nm, d in (("1s", a1), ("1s2p", a12), ("free", afr)):
    ok(f"E0b Liouville {nm}", [d["liouville_median"], d["liouville_q99"]],
       d["liouville_median"] <= pre["E0b"]["thr"]["median"] and d["liouville_q99"] <= pre["E0b"]["thr"]["q99"])
    ok(f"E0c reversibility {nm}", [d["reverse_median"], d["reverse_q99"]],
       d["reverse_median"] <= pre["E0c"]["thr"]["median"] and d["reverse_q99"] <= pre["E0c"]["thr"]["q99"])
ok("E1 F6 flow quantile", f6["quantile_err_flow"], f6["quantile_err_flow"] <= pre["E1"]["thr"]["flow"])
R["E1 old reader (report)"] = dict(value=f6["quantile_err_old_reader"])
D = cont["D"]
ok("E2 continuum: monotone", D, D[0] > D[1] > D[2])
ok("E2 continuum: order", cont["order"], cont["order"] >= pre["E2"]["thr"]["order_min"])
ok("E2 Bohm quantile identity", cont["bohm_quantile_identity_err"], cont["bohm_quantile_identity_err"] <= pre["E2"]["thr"]["bohm_identity"])
ok("E3 h1s flow max dq", h1s["flow"]["max_dq"], h1s["flow"]["max_dq"] <= pre["E3"]["thr"]["flow_max_dq"])
for k in ("bell", "bellflux"):
    ok(f"E3 h1s {k} jumps", [h1s[k]["jumps_per_surfer_per_tick"], h1s[k]["fraction_ever_jumped"]],
       h1s[k]["jumps_per_surfer_per_tick"] <= 1e-6 and h1s[k]["fraction_ever_jumped"] <= 1e-2)
ok("E4 beat flow quantile", beat["flow_quantile_err_max"], beat["flow_quantile_err_max"] <= pre["E4"]["thr"]["flow_quantile_err"])
for k in pre["E4"]["thr"]["schemes"]:
    ok(f"E4 beat z {k}", beat["max_abs_z"][k], beat["max_abs_z"][k] <= 4.0)
R["E4 beat old (report)"] = dict(value=beat["max_abs_z"]["old"])
for E, d in (("E5", a1), ("E6", a12), ("E7", afr)):
    for k in pre[E]["thr"]["schemes"]:
        ok(f"{E} z {k}", d["max_abs_z"][k], d["max_abs_z"][k] <= 4.0)
        ok(f"{E} TV ratio {k}", d["tv_ratio_final"][k], d["tv_ratio_final"][k] <= 1.2)
    for k in d["max_abs_z"]:
        if k not in pre[E]["thr"]["schemes"]:
            R[f"{E} z {k} (report)"] = dict(value=d["max_abs_z"][k])
th = pre["E5"]["thr"]
ok("E5 flow2 drift growth", a1["drift_growth_1000_over_100"]["flow2"], a1["drift_growth_1000_over_100"]["flow2"] <= th["flow2_drift_growth_1000_over_100"])
ok("E5 flow2 neq H ratio", a1["H_ratio"]["neq_flow2"], a1["H_ratio"]["neq_flow2"] >= th["flow2_neq_H_ratio_min"])
v = a1["vphi"]["flow2"]
ok("E5 flow2 v_phi/(Za sin)", v["vphi_over_Za_sin"], th["flow2_vphi_over_Za_sin"][0] <= v["vphi_over_Za_sin"] <= th["flow2_vphi_over_Za_sin"][1])
ok("E5 flow2 v_r rms/Za", v["vr_rms_over_Za"], v["vr_rms_over_Za"] <= th["flow2_vr_rms_over_Za_max"])
ok("E8 radial flow neq H ratio", beat["neq_H_ratio"]["flow"], beat["neq_H_ratio"]["flow"] >= 0.999)
ok("E9 orbitals proiettato z", orb["max_abs_z_flow"], orb["max_abs_z_flow"] <= 4.0)
R["E9 old orbital reader (report)"] = dict(value=[orb["max_abs_z_old"], orb["max_abs_z_old_rebirth"]])
npass = sum(1 for r in R.values() if r.get("PASS") is True)
nfail = sum(1 for r in R.values() if r.get("PASS") is False)
json.dump(dict(results=R, n_pass=npass, n_fail=nfail), open(os.path.join(OUT, "surf_results.json"), "w"), indent=1, default=float)
print(f"PASS {npass}  FAIL {nfail}")
for k, r in R.items():
    print(f"  {k:38s} {'PASS' if r.get('PASS') else ('FAIL' if r.get('PASS') is False else '----')}  {r['value']}")

# ------------------------------------------------------------------ figure
COL = dict(flow2="#1f77b4", flow="#17becf", bellflux="#2ca02c", bell="#bcbd22", ballistic="#ff7f0e", cdfasse="#9467bd",
           old="#d62728", old_rebirth="#e377c2")
NAME = dict(flow2="flusso (2 tick)", flow="flusso (1 tick)", bellflux="Bell in tempo di flusso", bell="Bell minimo per asse",
            ballistic="Bell balistico", cdfasse="CDF per asse (R3)", old="lettore vecchio", old_rebirth="vecchio + rinascite")
fig, ax = plt.subplots(2, 3, figsize=(17, 10.5))
# (a) 3-D 1s+2p: |z| vs time
x = a12["rec"]["t"]
for k in ("flow2", "flow", "bellflux", "ballistic", "old", "old_rebirth"):
    ax[0, 0].semilogy(x, np.maximum(np.abs(a12["rec"]["z"][k]), 1e-2), "-o", ms=3, color=COL[k], label=NAME[k])
ax[0, 0].axhline(4, color="k", ls="--", lw=1)
ax[0, 0].set_title("Atomo zig-zag 3-D (Z=60, 48³), 1s+2p che batte\n|z| dell'istogramma dei surfisti contro |ψ|² (soglia 4)")
ax[0, 0].set_xlabel("tick"); ax[0, 0].legend(fontsize=8)
# (b) trajectories of 1s (xy projection)
for k, c in (("flow2", COL["flow2"]), ("old", COL["old"]), ("ballistic", COL["ballistic"])):
    tr = np.array(a1["traj"][k])[:300]              # first 600 ticks
    for i in range(3):
        ax[0, 1].plot(tr[:, i, 0], tr[:, i, 1], "-", color=c, lw=0.8, alpha=0.9, label=NAME[k] if i == 0 else None)
ax[0, 1].plot(0, 0, "r*", ms=12)
ax[0, 1].set_aspect("equal"); ax[0, 1].legend(fontsize=8)
ax[0, 1].set_title("1s 3-D: tre surfisti per regola (vista dall'alto, 600 tick)\nBohm-Dirac: cerchi attorno all'asse dello spin")
ax[0, 1].set_xlabel("x (λ_C)"); ax[0, 1].set_ylabel("y (λ_C)")
# (c) fake mixing: H-bar of a |psi|^4 start, relative
tb = np.linspace(0, 20, len(beat["neq_H"]["flow"]))
for k in ("flow", "bellflux", "bell", "ballistic"):
    h = np.array(beat["neq_H"][k]); ax[0, 2].plot(tb, h / h[0], "-", color=COL[k], label="radiale: " + NAME[k])
t1 = np.array(a1["rec"]["t"])
for k in ("flow2", "flow", "bellflux", "ballistic", "cdfasse"):
    h = np.array(a1["H_neq"]["neq_" + k]); ax[0, 2].plot(np.concatenate([[0], t1]) / 50, h / h[0], "--", color=COL[k], label="1s 3-D: " + NAME[k])
ax[0, 2].set_ylim(0, 1.1)
ax[0, 2].set_title("Mescolamento finto: partenza da |ψ|⁴\n(stati dove Bohm non mescola; 1 = niente)")
ax[0, 2].set_xlabel("periodi di battimento (radiale) / tick/50 (1s 3-D)"); ax[0, 2].legend(fontsize=7)
# (d) randomness per surfer per tick
ks = ["flow2", "flow", "cdfasse", "bellflux", "bell", "ballistic"]
jm = [a1["stats"][k]["jumps"] for k in ks]
ax[1, 0].bar(range(len(ks)), jm, color=[COL[k] for k in ks])
ax[1, 0].set_xticks(range(len(ks))); ax[1, 0].set_xticklabels([NAME[k] for k in ks], rotation=25, fontsize=8)
ax[1, 0].set_title("1s 3-D: salti casuali per surfista per tick\n(le regole a flusso: zero numeri casuali)")
# (e) continuum limit
aa = np.array([0.2, 0.1, 0.05])
ax[1, 1].loglog(aa, D, "o-", color=COL["flow"], label="surfisti a flusso vs Bohm-Dirac continuo (max distanza)")
sp = [cont["per_a"][str(a)]["bell_rms_vs_flow_T40"] for a in (0.2, 0.1, 0.05)]
ax[1, 1].loglog(aa, sp, "s-", color=COL["bellflux"], label="Bell: dispersione attorno al flusso (t=40)")
ax[1, 1].loglog(aa, D[0] * aa / 0.2, "k:", lw=1, label="∝ a")
ax[1, 1].loglog(aa, sp[0] * np.sqrt(aa / 0.2), "k--", lw=1, label="∝ √a")
ax[1, 1].set_xlabel("passo del reticolo a (λ_C)"); ax[1, 1].set_ylabel("λ_C"); ax[1, 1].legend(fontsize=8)
ax[1, 1].set_title("Limite continuo in 1-D (pacchetto misto, zitterbewegung)")
# (f) radial beat and orbitals z
tt = np.linspace(0, 20, len(beat["chi2_z"]["flow"]))
for k in ("flow", "bellflux", "ballistic", "old"):
    ax[1, 2].semilogy(tt, np.maximum(np.abs(beat["chi2_z"][k]), 1e-2), "-", color=COL[k], label="radiale 1s+2s: " + NAME[k])
to = np.linspace(0, 10, len(orb["z_flow"]))
ax[1, 2].semilogy(to, np.maximum(np.abs(orb["z_flow"]), 1e-2), "--", color=COL["flow2"], label="orbitali analitici: proiettato")
ax[1, 2].semilogy(to, np.maximum(np.abs(orb["z_old"]), 1e-2), "--", color=COL["old"], label="orbitali: lettore vecchio")
ax[1, 2].semilogy(to, np.maximum(np.abs(orb["z_old_rebirth"]), 1e-2), "--", color=COL["old_rebirth"], label="orbitali: vecchio + rinascite")
ax[1, 2].axhline(4, color="k", ls="--", lw=1)
ax[1, 2].set_xlabel("periodi di battimento"); ax[1, 2].legend(fontsize=7)
ax[1, 2].set_title("Equivarianza su molti periodi: |z|")
fig.suptitle(f"Surfisti esatti — PASS {npass}, FAIL {nfail} (pre-registrazione surf_prereg.json)", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "surf_figure.png"), dpi=110)
print("written surf_figure.png")

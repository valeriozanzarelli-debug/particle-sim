"""tunnel_r2_verdict.py - giudizio dei pacchetti (R2) contro la pre-registrazione. Uscita lab/out/tunnel_r2_verdict.json"""
import json, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
PRE = json.load(open(os.path.join(OUT, "tunnel_prereg.json")))
AS = [0.05, 0.025, 0.0125]
R = {c: [json.load(open(os.path.join(OUT, f"tunnel_r2_{c}_{a}.json"))) for a in AS] for c in ("PK1", "PK2")}


def extrap(vals):
    a = np.array(AS)
    M = np.stack([np.ones(3), a, a ** 2], 1)
    return float(np.linalg.solve(M, np.array(vals))[0])


V = {}
V["E1_cono_di_luce"] = dict(exact_all=all(r["cone_exact"] for c in R for r in R[c]),
                            front_speed=[r["front_speed"] for c in R for r in R[c]])
V["E1_cono_di_luce"]["PASS"] = bool(V["E1_cono_di_luce"]["exact_all"] and
                                    max(abs(s - 1) for s in V["E1_cono_di_luce"]["front_speed"]) < 1e-10)
e5 = {}
ok = True
for c in R:
    P = PRE["pacchetti"][c]
    adv = extrap([r["peak"]["advance"] for r in R[c]])
    mar = extrap([r["peak"]["margin_vs_light"] for r in R[c]])
    e5[c] = dict(advance=adv, advance_pred=P["peak_advance_pred"], err_over_L=abs(adv - P["peak_advance_pred"]) / P["L"],
                 margin=mar, margin_pred=P["margin_vs_light_pred"], superluminal_peak=bool(mar > 0),
                 superluminal_pred=P["superluminal_peak_pred"])
    ok &= e5[c]["err_over_L"] < 0.02 and e5[c]["superluminal_peak"] == e5[c]["superluminal_pred"]
V["E5_picco_pacchetto"] = dict(data=e5, PASS=bool(ok))
e6 = max(abs(r["dwell_identity"] / r["dwell"] - 1) for c in R for r in R[c])
V["E6_larmor_dwell_pacchetto"] = dict(max_rel_err=e6, soglia=1e-4, PASS=bool(e6 < 1e-4))
ey = max(abs(r["tau_y"] / r["stationary"]["tau_y"] - 1) for c in R for r in R[c])
ez = max(abs(r["tau_z"] / r["stationary"]["tau_z"] - 1) for c in R for r in R[c])
emz = max(abs(r["tau_mz"] / r["stationary"]["tau_mz"] - 1) for c in R for r in R[c])
emy = max(abs(r["tau_my"] - r["stationary"]["tau_my"]) / r["stationary"]["tau_y"] for c in R for r in R[c])
V["E7_larmor_pacchetto"] = dict(max_rel_err_tau_y=ey, max_rel_err_tau_z=ez, PASS=bool(ey < 1e-3 and ez < 1e-2),
                                extra_mass_clock=dict(max_rel_err_tau_mz=emz, max_abs_err_tau_my_over_tau_y=emy))
e12 = {}
ok = True
for c in R:
    L = PRE["pacchetti"][c]["L"]
    for r, a in zip(R[c], AS):
        b = r["bohm"]
        idn = abs(r["PT"] * b["tauB_T_traj"] + r["PR"] * b["tauB_R"] - r["dwell"]) / r["dwell"]
        e12[f"{c}_{a}"] = dict(tauB_T=b["tauB_T_traj"], tauB_T_min=b["tauB_T_min"], L=L, tauB_R=b["tauB_R"], identity_err=idn,
                               traj_vs_density=abs(b["tauB_T_traj"] / b["tauB_T_density"] - 1))
        ok &= b["tauB_T_min"] >= L - 2 * a and idn < 1e-3
V["E12_bohm"] = dict(data=e12, PASS=bool(ok))
e13 = {}
for c in R:
    r = R[c][-1]
    A = r["arrival"]
    e13[c] = dict(int_absorber=A["int_F"], passed=A["passed"], ratio=A["int_F"] / A["passed"],
                  mean_t_absorber=A["mean_t_absorber"], mean_t_flux=A["mean_t_flux"], mean_t_free=A["mean_t_free"],
                  min_J_over_max=A["min_J_over_max"])
V["E13_arrivo"] = dict(data=e13, soglia=1e-6,
                       PASS=bool(all(abs(e13[c]["ratio"] - 1) < 1e-6 for c in e13)))
json.dump(V, open(os.path.join(OUT, "tunnel_r2_verdict.json"), "w"), indent=1, default=float)
for k, v in V.items():
    print(k, "PASS" if v["PASS"] else "FAIL", json.dumps({kk: vv for kk, vv in v.items() if kk != "PASS"}, default=float)[:600])

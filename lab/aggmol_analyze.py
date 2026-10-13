"""aggmol_analyze.py - verdicts of study aggmol_ against lab/out/aggmol_prereg.json. Writes lab/out/aggmol_verdict.json."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_core as mc
import aggmol_exact as ex

O = mc.OUT
L = json.load(open(os.path.join(O, "aggmol_learn.json")))
A = json.load(open(os.path.join(O, "aggmol_atom.json")))
H = json.load(open(os.path.join(O, "aggmol_h2p.json")))
FR = json.load(open(os.path.join(O, "aggfull_fineref.json")))["levels"]
CHEM = 1.594e-3
EHeV = mc.ALPHA ** 2 * mc.MC2
th = mc.TH
mk = np.tan(th) / th
V = {}


def verdict(name, ok, **kw):
    V[name] = dict(PASS=bool(ok), **kw)
    print(f"{name}: {'PASS' if ok else 'FAIL'}  " + ", ".join(f"{k}={v}" for k, v in kw.items()))


# ---- L
L1 = (L.get("chosen_two_atom_level_B") == 1024 and not L["B2048_from512"]["free"]["tried"][-1]["ok"]
      and not L["B2048_from1024"]["free"]["tried"][-1]["ok"])
verdict("L1", L1, chosen_B=L.get("chosen_two_atom_level_B"),
        a0_2_energy_rel=[t["check"]["energy_rel"] for t in L["B2048_from512"]["free"]["tried"] + L["B2048_from1024"]["free"]["tried"]],
        a0_4_energy_rel=L["B1024_from512"]["free"]["tried"][-1]["check"]["energy_rel"])
md = {k: L[k]["free"]["markov_w4"]["max_c4_dev"] for k in ("B1024_from512", "B2048_from512", "B2048_from1024")}
verdict("L2", all(v > 1 for v in md.values()), markov_c4_dev=md)
ph = L["B1024_from512"]["coupling"]["physics"]
dev = dict(gk=ph["gk"][0] / -0.125 - 1, gx=ph["gx"][0] / 0.25 - 1, gso=ph["gso"][0] / 0.25 - 1)
verdict("L3", all(abs(v) <= 0.02 for v in dev.values()), stencil=L["B1024_from512"]["coupling"]["chosen"], rel_dev=dev)

# ---- A
s, p = A["s_a0_4"]["states"], A["p_a0_4"]["states"]
ps, pp = A["pert_s_a0_4"], A["pert_p_a0_4"]["2p"]
ev = np.array(pp["dV_eV"])
E = {"1s1/2": s["1s"]["E0_eV"] + np.mean(ps["1s"]["dV_eV"]), "2s1/2": s["2s"]["E0_eV"] + np.mean(ps["2s"]["dV_eV"]),
     "2p1/2": p["2p"]["E0_eV"] + ev[:2].mean(), "2p3/2": p["2p"]["E0_eV"] + ev[2:].mean()}
rows = {}
okA1 = okA2 = True
for k, e in E.items():
    fine, dirac = -FR[k]["E_bind_eV"], -FR[k]["dirac_bind_eV"]
    r1, r2 = e / fine - 1, e / dirac - 1
    rows[k] = dict(E_eV=e, vs_fine_rel=r1, vs_fine_mEh=1e3 * (e - fine) / EHeV, vs_dirac_rel=r2, vs_dirac_mEh=1e3 * (e - dirac) / EHeV)
    okA1 &= abs(r1) <= 2.5e-3 and abs(e - fine) / EHeV <= CHEM
    okA2 &= abs(r2) <= 2.5e-3 and abs(e - dirac) / EHeV <= CHEM
verdict("A1", okA1, rows=rows)
verdict("A2", okA2)
fsD = FR["2p3/2"]["dirac_bind_eV"] * -1 - FR["2p1/2"]["dirac_bind_eV"] * -1
verdict("A3", abs(pp["splitting_eV"] / fsD - 1) <= 0.02 and abs(pp["interval_ratio"] / 0.5 - 1) <= 0.01,
        FS_ueV=1e6 * pp["splitting_eV"], dirac_ueV=1e6 * fsD, rel=pp["splitting_eV"] / fsD - 1, ratio=pp["interval_ratio"])

# ---- C
c = {R: H[f"curve_4_R{R}"] for R in (1.4, 2.0)}
verdict("C1", all(abs(v["err_mEh"]) <= 1e3 * CHEM and abs(v["rel"]) <= 1e-3 for v in c.values()),
        err_mEh={R: v["err_mEh"] for R, v in c.items()}, rel={R: v["rel"] for R, v in c.items()})


def De(f):
    return H[f"atom_{f}"]["E_Eh"] - H[f"curve_{f}_R2.0"]["E_tot_Eh"]


DeX = -0.5 - (ex.energy(2.0)[0] + 0.5)
verdict("C2", abs(De(4) - DeX) <= CHEM, De_mEh=1e3 * De(4), exact_mEh=1e3 * DeX, diff_mEh=1e3 * (De(4) - DeX))
Rs = np.array([1.6, 1.8, 2.0, 2.2, 2.4])


def fitmin(E):
    cf = np.polyfit(Rs, E, 4)
    rr = np.roots(np.polyder(cf)); rr = rr[np.abs(rr.imag) < 1e-9].real; rr = rr[(rr > 1.6) & (rr < 2.4)][0]
    return rr, np.polyval(np.polyder(cf, 2), rr)


Rr, kr = fitmin(np.array([H[f"curve_4_R{R}"]["E_tot_Eh"] for R in Rs]))
Rx, kx = fitmin(np.array([ex.energy(R)[0] + 1 / R for R in Rs]))
verdict("C3", abs(Rr / Rx - 1) <= 5e-3 and abs(kr / kx - 1) <= 0.02, Re_rule=Rr, Re_exact_samefit=Rx, Re_rel=Rr / Rx - 1, k_rel=kr / kx - 1)
exc = {st: H[f"exc_4_{st}"] for st in ("1su", "1pu", "2sg")}
verdict("C4", all(abs(v["rel"]) <= 1e-3 for v in exc.values()), rel={k: v["rel"] for k, v in exc.items()})
E0 = H["curve_4_R2.0"]["E_el_Eh"]
sh = [H[f"shift_4_{i}"]["E_el_Eh"] / E0 - 1 for i in range(3)]
verdict("C5", max(abs(x) for x in sh) <= 1e-4, rel=sh)
verdict("C6", abs(H["rot_4"]["E_el_Eh"] / E0 - 1) <= 1e-4, rel=H["rot_4"]["E_el_Eh"] / E0 - 1)
n = exc["2sg"]
verdict("C7", n["xi_spread_rel"] <= 5e-3 and n["xi_spread_rel"] <= n["r_spread_rel"] / 5 and abs(n["xi0_err_rel"]) <= 0.01,
        xi_mean=n["xi_mean"], xi0_exact=n["xi0_exact"], xi_spread=n["xi_spread_rel"], r_spread=n["r_spread_rel"], xi0_err_rel=n["xi0_err_rel"])
key = mk * ex.energy(mk * 2.0)[0]
lad = {f: H[f"curve_{f}_R2.0"]["E_el_Eh"] / key - 1 for f in (4, 8, 16)}
des = {f: 1e3 * (De(f) - DeX) for f in (4, 8, 16)}
verdict("C8", abs(lad[4]) > abs(lad[8]) > abs(lad[16]) and abs(lad[16]) <= 2e-4 and all(abs(d) <= 1e3 * CHEM for d in des.values()),
        rel_vs_masskey=lad, rel_vs_exact={f: H[f"curve_{f}_R2.0"]["rel"] for f in (4, 8, 16)}, De_diff_mEh=des, m_kin=mk)
g = {R: H[f"gap_4_R{R}"] for R in (4.0, 6.0)}
verdict("C9", all(abs(v["rel"]) <= 0.01 for v in g.values()), rel={R: v["rel"] for R, v in g.items()},
        period_au={R: v["period_au"] for R, v in g.items()})
npass = sum(v["PASS"] for v in V.values())
V["summary"] = f"{npass} PASS / {len(V) - npass} FAIL"
print(V["summary"])
json.dump(V, open(os.path.join(O, "aggmol_verdict.json"), "w"), indent=1, default=float)

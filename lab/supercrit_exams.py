"""supercrit_exams.py - collects all runs and judges them against lab/out/supercrit_prereg.json.
Writes lab/out/supercrit_exams.json."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda n: json.load(open(os.path.join(HERE, "out", f"supercrit_{n}.json")))
from supercrit_core import ALPHA, dirac_critical_Z, V_shell
from supercrit_prereg import GODUNOV
v1, c1, c2, g1, g2 = L("v1"), L("c1"), L("c2"), L("g1"), L("g2")
pw, pd = L("posthoc_wide"), L("posthoc_decay")
E = {}
pf = lambda ok: "PASS" if ok else "FAIL"
# V1
errs = {Z: v1[Z]["extrap"]["rel_err"] for Z in v1}
E["V1a"] = dict(values={Z: errs[Z] for Z in ("20", "60", "92", "110")}, result=pf(all(abs(errs[Z]) <= 1e-4 for Z in ("20", "60", "92", "110"))))
E["V1b"] = dict(values=errs["120"], result=pf(abs(errs["120"]) <= 1e-3))
fin137 = sorted(v1["137"]["runs"].values(), key=lambda x: x["a"])[0]
E["V1c"] = dict(values={"130": errs["130"], "137_finest_level_err": fin137["rel_err"],
                        "137_note": "at a = 0.03125 the lowest level is at -0.866; for a <= 0.0156 it has dived below -1 "
                                    "(see C1) and the level that tends to Dirac's 1s (0.0229) is at 0.231 at a = 4.9e-4"},
                result="Z=130 " + pf(abs(errs["130"]) <= 1e-2) + "; Z=137 FAIL (as predicted: limit)")
# C1
rows = sorted([v for k, v in c1.items() if k[0].isdigit()], key=lambda x: x["a"])
A = np.array([r["a"] for r in rows]); Zc = np.array([r["Zcr"] for r in rows])
from scipy.optimize import curve_fit
p, cov = curve_fit(lambda a, Zs, c, q: Zs + c * a**q, A[:5], Zc[:5], p0=[131.3, 50, 0.6])
gstar = np.sqrt(1 - (p[0] * ALPHA)**2)
mono = bool(np.all(np.diff(Zc) > 0))
E["C1a"] = dict(values=dict(a=A.tolist(), Zcr=Zc.tolist(), monotonic=mono, all_above_137=bool(np.all(Zc > 1 / ALPHA))),
                result="FAIL (Z_cr(a) < 137 for a <= 0.01: monotonic yes, above 1/alpha no)")
E["C1_fit"] = dict(Z_star=p[0], zeta_star=p[0] * ALPHA, c=p[1], power=p[2], two_gamma_star=2 * gstar,
                   note="post-hoc fit Z_cr(a) = Z* + c a^q on the 5 finest a; theory of a pole of the lattice inner "
                        "boundary number at zeta*: q = 2 gamma(zeta*)")
E["C1b"] = dict(values="R_eff exists only for a = 0.1, 0.03 (Z_cr > 137); for a <= 0.01 Z_cr < 137 has no cut-off "
                       "equivalent", result="FAIL")
# C2, C3
ex = c2["shell_extrap"]
E["C2a"] = dict(values=dict(walk=ex["Zcr_quad3"], walk_quad4=ex["Zcr_quad4"], lit=177.46,
                            continuum=c2["continuum"]["shell"]), result=pf(abs(ex["Zcr_quad3"] - 177.46) <= 0.3))
E["C2b"] = dict(values=dict(walk=ex["E170_quad3"], lit=-0.664), result=pf(abs(ex["E170_quad3"] / -0.664 - 1) <= 0.005))
es = c2["sphere_extrap"]
E["C3"] = dict(values=dict(walk=es["Zcr_quad3"], continuum=c2["continuum"]["sphere"], R_fm=c2["continuum"]["sphere_R_fm"]),
               result=pf(169 <= es["Zcr_quad3"] <= 175 and abs(es["Zcr_quad3"] - c2["continuum"]["sphere"]) <= 0.3))
# G1
gt = {}
for Z in g1:
    x = g1[Z]["extrap"]
    if int(Z) in GODUNOV:
        gt[Z] = dict(E=x["E_res_quad3"], G=x["Gamma_quad3"], E_lit=GODUNOV[int(Z)][0], G_lit=GODUNOV[int(Z)][1],
                     dE=x["dE_rel"], dG=x["dG_rel"])
big = [Z for Z in gt if int(Z) >= 182]
E["G1a"] = dict(values={Z: gt[Z]["dG"] for Z in big},
                result=pf(all(abs(gt[Z]["dG"]) <= 0.10 for Z in big)) + " (pre-registered fit, window +-12 Gamma: fails at "
                       "Z = 210, 232); post-hoc window +-3 Gamma: " +
                       str({Z: pw[Z]["walk_win3_extrap"]["dG_rel"] for Z in pw if "walk_win3_extrap" in pw[Z]}))
E["G1b"] = dict(values=gt["180"]["dG"], result=pf(abs(gt["180"]["dG"]) <= 0.3))
E["G1c"] = dict(values={Z: gt[Z]["dE"] for Z in gt}, result=pf(all(abs(gt[Z]["dE"]) <= 0.01 for Z in gt)))
# G2
r = {Z: g2[Z]["ratio_time_phase"] for Z in g2}
ok = all(abs(r[Z] - 1) <= (0.10 if Z == "232" else 0.05) for Z in r)
E["G2"] = dict(values=r, result=pf(ok) + "; post-hoc (floor fit, narrow-window phase): " +
               str({Z: round(pd[Z]["ratio"], 4) for Z in pd}))
# G3 Gamow
Zs = [Z for Z in sorted(g1, key=int) if 180 <= int(Z) <= 188]
eta, lg = [], []
for Z in Zs:
    Er, G = g1[Z]["extrap"]["E_res_quad3"], g1[Z]["extrap"]["Gamma_quad3"]
    eta.append(int(Z) * ALPHA * abs(Er) / np.sqrt(Er**2 - 1)); lg.append(np.log(G))
s = -np.polyfit(2 * np.pi * np.array(eta), lg, 1)[0]
E["G3"] = dict(values=dict(slope=s, Z=Zs, eta=eta), result=pf(abs(s - 1) <= 0.2))
# D1
d = g2["200"]
E["D1"] = dict(values=dict(peak=d["spec_peak"], E_phase=d["E_phase"], fwhm=d["spec_fwhm"], Gamma=d["Gamma_phase"]),
               result=pf(abs(d["spec_peak"] / d["E_phase"] - 1) <= 0.01 and abs(d["spec_fwhm"] / d["Gamma_phase"] - 1) <= 0.15))
for k, v in E.items():
    print(k, v["result"] if "result" in v else "", "|", str(v.get("values", v))[:300])
json.dump(E, open(os.path.join(HERE, "out", "supercrit_exams.json"), "w"), indent=1, default=float)

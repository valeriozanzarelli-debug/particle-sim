"""glue_verdict.py - grades lab/out/glue_results.json against lab/out/glue_prereg.json.
Writes lab/out/glue_verdict.json and prints a summary."""
import os, json
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
pre = json.load(open(os.path.join(OUT, "glue_prereg.json")))
res = json.load(open(os.path.join(OUT, "glue_results.json")))["results"]
AT = pre["literature_checked_before_prereg"]["AT2017"]
LITC = AT["continuum_m_over_sqrt_sigma"]
R2L, R0L = 7.762 / 4.7367, 6.861 / 4.7367
betas = sorted(res.keys(), key=float)


def ratio(a, ea, b, eb):
    r = a / b
    return r, r * np.hypot(ea / a, eb / b)


def grade(val, err, lit, tol):
    if not np.isfinite(val) or not np.isfinite(err):
        return "INCONCLUSIVE"
    if err > 0.08 * abs(val):
        return "INCONCLUSIVE"
    return "PASS" if abs(val - lit) <= 2 * err + tol * lit else "FAIL"


V = {"per_beta": {}, "grades": {}}
for b in betas:
    r = res[b]
    s, es = r["a_sqrt_sigma"], r["a_sqrt_sigma_err"]
    m0, e0 = r["A1"][0]["m"], r["A1"][0]["err"]
    m1, e1 = r["A1"][1]["m"], r["A1"][1]["err"]
    m2, e2 = r["E"][0]["m"], r["E"][0]["err"]
    d = dict(a_sqrt_sigma=(s, es), am0=(m0, e0), am0s=(m1, e1), am2=(m2, e2),
             R2=ratio(m2, e2, m0, e0), R0s=ratio(m1, e1, m0, e0),
             m0_s=ratio(m0, e0, s, es), m0s_s=ratio(m1, e1, s, es), m2_s=ratio(m2, e2, s, es),
             flags={"0++": r["A1"][0]["flag"], "0++*": r["A1"][1]["flag"], "2++": r["E"][0]["flag"],
                    "torelon": r["P"][0]["flag"]},
             windows={"0++": [r["A1"][0]["t_min"], r["A1"][0]["t_max"]], "0++*": [r["A1"][1]["t_min"], r["A1"][1]["t_max"]],
                      "2++": [r["E"][0]["t_min"], r["E"][0]["t_max"]], "torelon": [r["P"][0]["t_min"], r["P"][0]["t_max"]]},
             n_meas=r["n_meas"])
    # G4 dispersion
    E, eE = r["A1p"][0]["m"], r["A1p"][0]["err"]
    p2 = (2 * np.pi / r["L"]) ** 2
    c2 = (E ** 2 - m0 ** 2) / p2
    ec2 = np.hypot(2 * E * eE, 2 * m0 * e0) / p2
    d["E_p"] = (E, eE); d["c2"] = (c2, ec2)
    # G1b Creutz
    cr = r["creutz"]; sig = s * s
    best = None
    for R in sorted(cr, key=int):
        v, e = cr[R]
        if np.isfinite(v) and e <= 0.10 * abs(v):
            best = int(R)
    if best is not None:
        v, e = cr[str(best)]
        corr = v - np.pi / (24 * (best - 0.5) ** 2)
        d["creutz_best"] = dict(R=best, chi=v, err=e, chi_minus_luscher=corr, ratio_to_sigma=corr / sig)
    V["per_beta"][b] = d

G = V["grades"]
G["G1a_sigma_vs_literature"] = {b: ("PASS" if abs(V["per_beta"][b]["a_sqrt_sigma"][0] / AT["lattice"][b]["a_sqrt_sigma"] - 1) <= 0.03
                                    else "FAIL") + f" ({V['per_beta'][b]['a_sqrt_sigma'][0]:.4f} vs {AT['lattice'][b]['a_sqrt_sigma']})"
                                for b in betas}
G["G1b_creutz"] = {}
for b in betas:
    cb = V["per_beta"][b].get("creutz_best")
    G["G1b_creutz"][b] = "NO DATA" if cb is None else (
        ("PASS" if abs(cb["ratio_to_sigma"] - 1) <= 0.15 else "FAIL") + f" (R={cb['R']}, ratio {cb['ratio_to_sigma']:.3f})")
G["G2a_mass_gap_vs_literature"] = {b: ("PASS" if abs(V["per_beta"][b]["am0"][0] / AT["lattice"][b]["am0"] - 1) <= 0.05
                                       else "FAIL") + f" ({V['per_beta'][b]['am0'][0]:.4f}({V['per_beta'][b]['am0'][1]:.4f}) vs {AT['lattice'][b]['am0']})"
                                   for b in betas}
fin = "12" if "12" in betas else betas[-1]
pf = V["per_beta"][fin]
G["G2b_R2_finest"] = grade(*pf["R2"], R2L, 0.05)
G["G2c_R0s_finest"] = grade(*pf["R0s"], R0L, 0.07)
G["G2d_m_over_sqrt_sigma_finest"] = {"0++": grade(*pf["m0_s"], LITC["0++"][0], 0.04),
                                     "2++": grade(*pf["m2_s"], LITC["2++"][0], 0.06),
                                     "0++*": grade(*pf["m0s_s"], LITC["0++*"][0], 0.08)}

# continuum: c + d x, x = (a sqrt sigma)^2, weighted
def cfit(key):
    x = np.array([V["per_beta"][b]["a_sqrt_sigma"][0] ** 2 for b in betas])
    y = np.array([V["per_beta"][b][key][0] for b in betas]); e = np.array([V["per_beta"][b][key][1] for b in betas])
    A = np.vstack([np.ones_like(x), x]).T / e[:, None]
    cov = np.linalg.inv(A.T @ A)
    p = cov @ (A.T @ (y / e))
    chi2 = float(np.sum(((y - (p[0] + p[1] * x)) / e) ** 2))
    return float(p[0]), float(np.sqrt(cov[0, 0])), float(p[1]), chi2

if len(betas) >= 3:
    C = {k: cfit(k) for k in ("R2", "R0s", "m0_s", "m2_s", "m0s_s")}
    V["continuum"] = C
    G["G2e_continuum"] = {"R2": grade(C["R2"][0], C["R2"][1], R2L, 0.05),
                          "R0s": grade(C["R0s"][0], C["R0s"][1], R0L, 0.07),
                          "m0/sqrt(sigma)": grade(C["m0_s"][0], C["m0_s"][1], LITC["0++"][0], 0.04),
                          "m2/sqrt(sigma)": grade(C["m2_s"][0], C["m2_s"][1], LITC["2++"][0], 0.06)}
    # G3: beta * a sqrt(sigma) vs 1/beta, linear
    bb = np.array([float(b) for b in betas])
    y = np.array([float(b) * V["per_beta"][b]["a_sqrt_sigma"][0] for b in betas])
    e = np.array([float(b) * V["per_beta"][b]["a_sqrt_sigma"][1] for b in betas])
    A = np.vstack([np.ones_like(bb), 1 / bb]).T / e[:, None]
    cov = np.linalg.inv(A.T @ A); p = cov @ (A.T @ (y / e))
    sg = p[0] / 4; esg = np.sqrt(cov[0, 0]) / 4
    V["G3"] = dict(beta_a_sqrt_sigma=y.tolist(), err=e.tolist(), sqrt_sigma_over_g2=float(sg), err_sg=float(esg),
                   slope=float(p[1]))
    G["G3_no_mass_parameter"] = ("PASS" if (abs(sg / 0.3353 - 1) <= 0.05 and np.all(np.diff(y) < 0)) else "FAIL") + \
        f" (sqrt(sigma)/g^2 = {sg:.4f}({esg:.4f}) vs 0.3353; beta a sqrt(sigma) = {np.round(y, 4).tolist()})"
g4 = {}
for b in ("9", "12"):
    if b in V["per_beta"]:
        c2, ec2 = V["per_beta"][b]["c2"]
        g4[b] = ("INCONCLUSIVE" if ec2 > 0.3 else ("PASS" if abs(c2 - 1) <= 2 * ec2 else "FAIL")) + f" (c^2 = {c2:.3f}({ec2:.3f}))"
G["G4_mass_is_energy"] = g4

json.dump(V, open(os.path.join(OUT, "glue_verdict.json"), "w"), indent=1)
for b in betas:
    d = V["per_beta"][b]
    f = lambda t: f"{t[0]:.4f}({t[1]:.4f})"
    print(f"beta {b} N={d['n_meas']}: a sqrt(s)={f(d['a_sqrt_sigma'])} am0={f(d['am0'])} am0*={f(d['am0s'])} am2={f(d['am2'])} "
          f"R2={f(d['R2'])} R0s={f(d['R0s'])} m0/s={f(d['m0_s'])} m0*/s={f(d['m0s_s'])} m2/s={f(d['m2_s'])} c2={f(d['c2'])}")
    print("   windows", d["windows"], d["flags"])
if "continuum" in V:
    for k, v in V["continuum"].items():
        print(f"continuum {k}: {v[0]:.4f}({v[1]:.4f}) slope {v[2]:.2f} chi2 {v[3]:.2f}")
print(json.dumps(G, indent=1))

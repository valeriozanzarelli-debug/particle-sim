"""assist_limits.py - the honest limits of the assist_ study (report only, not graded).
1. 1+1 -> 3+1: for a homogeneous linearly polarised E(t) the 3+1 spinor QED kinetic equation is the 1+1 Dirac problem with
   m -> m_perp = sqrt(m^2 + p_perp^2), times 2 (spin):  n_3D = 2 int d^2p_perp/(2pi)^2 N_1D(m_perp)
   = (1/(2 pi)) int_1^inf d(m_perp^2) N_1D(m_perp).  N_1D(m_perp) from the rule (zig-zag, scaled units) and from Dirac.
   Scaling: f(p; m_perp, E, w, Om) = f(p/m_perp; 1, E/m_perp^2, w/m_perp, Om/m_perp).
2. real lasers: E/E_S, omega/m, combined gamma and the worldline exponent for laser-class numbers.
Writes lab/out/assist_limits.json."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assist_core import Field, walk_pairs, dirac_pairs, S_wi, s_wi_sauter, OUT, PI

t0 = time.time()
res = {}
E1, GS, eps, g = 0.1, 0.2, 0.1, 3.0
m2 = np.arange(1.0, 1.40001, 0.02)
rows = []
for M2 in m2:
    M = np.sqrt(M2)
    p = np.arange(-4, 4.00001, 0.01)
    row = dict(m_perp2=float(M2))
    for lab, F in (("both", Field(E1 / M2, GS * E1 / M, eps * E1 / M2, g * E1 / M, "sauter")),
                   ("strong", Field(E1 / M2, GS * E1 / M)),
                   ("weak", Field(0.0, GS * E1 / M, eps * E1 / M2, g * E1 / M, "sauter", E1_shape="none"))):
        fd = dirac_pairs(F, p, 0.01)
        fl = walk_pairs(F, p, 0.025)
        # N per unit length in physical units: dp_phys = M dp_scaled
        row[lab + "_dirac"] = float(M * np.trapezoid(fd, p) / (2 * PI))
        row[lab + "_rule"] = float(M * np.trapezoid(fl, p) / (2 * PI))
    rows.append(row)
    print(row, f"{time.time()-t0:.0f}s", flush=True)
res["mperp_rows"] = rows
d = {}
for k in ("both", "strong", "weak"):
    for s in ("dirac", "rule"):
        d[f"n3D_{k}_{s}"] = float(np.trapezoid([r[f"{k}_{s}"] for r in rows], m2) / (2 * PI))
    d[f"N1D_{k}"] = rows[0][f"{k}_rule"]
d["enh_1D"] = rows[0]["both_rule"] / (rows[0]["strong_rule"] + rows[0]["weak_rule"])
d["enh_3D"] = d["n3D_both_rule"] / (d["n3D_strong_rule"] + d["n3D_weak_rule"])
d["tail_check_last_over_first_both"] = rows[-1]["both_rule"] / rows[0]["both_rule"]
d["config"] = dict(E1=E1, gs=GS, eps=eps, gamma=g)
res["3D"] = d
print(d)

# ---- real lasers (SI from zz_core constants)
c, eps0, ES = 299792458.0, 8.8541878128e-12, 1.32328e18
lasers = []
for name, I in (("10 PW class today (ELI-NP / Apollon), 1e23 W/cm^2", 1e23),
                ("next generation (~100 PW), 1e24 W/cm^2", 1e24), ("1e25 W/cm^2 (beyond plans)", 1e25)):
    E = np.sqrt(2 * I * 1e4 / (c * eps0))
    lasers.append(dict(name=name, I_Wcm2=I, E_Vm=float(E), E_over_ES=float(E / ES),
                       schwinger_exponent=float(PI * ES / E)))
fast = [("XFEL 10 keV", 10e3), ("XFEL 25 keV", 25e3), ("gamma 100 keV", 100e3), ("gamma 300 keV", 300e3)]
tab = []
for L in lasers:
    for fn, w in fast:
        wm = w / 510998.95
        gam = wm / L["E_over_ES"]
        S = PI * s_wi_sauter(gam) / L["E_over_ES"]
        tab.append(dict(strong=L["name"], fast=fn, omega_over_m=wm, gamma=float(gam), exponent_WI=float(S),
                        large_gamma_limit_2pi_m_over_omega=float(2 * PI / wm)))
res["lasers"] = lasers
res["assisted_table"] = tab
for r in tab:
    print(r)
json.dump(res, open(os.path.join(OUT, "assist_limits.json"), "w"), indent=1)
print("done", f"{time.time()-t0:.0f}s")

"""supercrit_posthoc_decay.py - POST-HOC (after G2 failed): decay rate of P_in(t) fitted as A exp(-G t) + C
(C = weight of the initial box state on the true bound states 2s, 3s, ... that never decays) and compared with
the narrow-window (+-3 Gamma) phase width at the same a = 0.002. Writes lab/out/supercrit_posthoc_decay.json"""
import json, os
import numpy as np
from scipy.optimize import curve_fit
HERE = os.path.dirname(os.path.abspath(__file__))
g2 = json.load(open(os.path.join(HERE, "out", "supercrit_g2.json")))
pw = json.load(open(os.path.join(HERE, "out", "supercrit_posthoc_wide.json")))
res = {}
for Z in ("190", "200", "210", "232"):
    d = np.load(os.path.join(HERE, "out", f"supercrit_g2_{Z}.npz"))
    t, P = d["t"], d["Pin"]
    G0 = g2[Z]["Gamma_time"]
    s = P < 0.9
    f = lambda t, A, G, C: A * np.exp(-G * t) + C
    p, _ = curve_fit(f, t[s], P[s], p0=[1.0, G0, 0.005])
    Gref = pw[Z]["walk_win3_extrap"]["Gammas"][2] if "walk_win3_extrap" in pw.get(Z, {}) else g2[Z]["Gamma_phase"]
    res[Z] = dict(Gamma_fit_floor=p[1], floor=p[2], Gamma_phase_win3_a0002=Gref, ratio=p[1] / Gref,
                  spec_fwhm_ratio=g2[Z]["spec_fwhm"] / Gref, raw_ratio=G0 / Gref)
    print(Z, res[Z])
json.dump(res, open(os.path.join(HERE, "out", "supercrit_posthoc_decay.json"), "w"), indent=1)

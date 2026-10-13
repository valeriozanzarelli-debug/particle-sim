"""supercrit_posthoc_wide.py - POST-HOC (written after G1 showed +11 % / +22 % at Z = 210 / 232): is the gap a
property of the rule or of the Breit-Wigner fit? Same fit (window +-12 Gamma, quadratic background) on the
CONTINUUM Dirac phase of the same shell nucleus, and on narrower windows (+-3 Gamma) for walk and continuum.
Writes lab/out/supercrit_posthoc_wide.json"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supercrit_core import (SWalk, V_shell, ALPHA, dirac_phase, fit_breit_wigner, bloch_amplitudes)
from supercrit_prereg import GODUNOV

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "supercrit_posthoc_wide.json")
g1 = json.load(open(os.path.join(os.path.dirname(OUT), "supercrit_g1.json")))
res = {}
for Z in (190, 200, 210, 232):
    Za = Z * ALPHA
    Vf = lambda r: V_shell(r, Za, 0.031)
    ex = g1[str(Z)]["extrap"]
    Er, G = ex["E_res_quad3"], ex["Gamma_quad3"]
    out = {}
    for win in (12, 3):
        E = np.linspace(Er - win * G, Er + win * G, 121)
        ph = np.unwrap([dirac_phase(Vf, e, 25.0) for e in E])
        Ec, Gc, rms = fit_breit_wigner(E, ph)
        a = 0.0005
        w = SWalk(-1, a, 27.0, Vf)
        phw = np.unwrap([bloch_amplitudes(w, e, int(25.0 / a), 10)[0] for e in E])
        Ew, Gw, rmsw = fit_breit_wigner(E, phw)
        out[f"win{win}"] = dict(cont_E=Ec, cont_G=Gc, cont_rms=rms, walk_a0005_E=Ew, walk_a0005_G=Gw)
        print(Z, win, "continuum", Ec, Gc, " walk a=5e-4", Ew, Gw, " Godunov", GODUNOV[Z], flush=True)
    res[str(Z)] = out
json.dump(res, open(OUT, "w"), indent=1)

# walk a -> 0 with the narrow window (+-3 Gamma), same extrapolation as G1 (quadratic on 3 finest a)
for Z in (195, 200, 210, 232):
    Za = Z * ALPHA
    Vf = lambda r: V_shell(r, Za, 0.031)
    d = g1[str(Z)]
    Gs, Es, As = [], [], []
    for a in (0.0005, 0.001, 0.002):
        r0 = d[f"{a:.6g}"]
        Er, G = r0["E_res"], r0["Gamma"]
        w = SWalk(-1, a, 27.0, Vf)
        for _ in range(2):
            E = np.linspace(Er - 3 * G, Er + 3 * G, 121)
            phw = np.unwrap([bloch_amplitudes(w, e, int(25.0 / a), 10)[0] for e in E])
            Er, G, rms = fit_breit_wigner(E, phw)
        Gs.append(G); Es.append(Er); As.append(a)
    G0 = float(np.polyval(np.polyfit(As, Gs, 2), 0)); E0 = float(np.polyval(np.polyfit(As, Es, 2), 0))
    res.setdefault(str(Z), {})["walk_win3_extrap"] = dict(Gamma=G0, E_res=E0, Gammas=Gs, Es=Es, a=As,
                                                         dG_rel=(G0 - GODUNOV[Z][1]) / GODUNOV[Z][1])
    print(Z, "walk win3 a->0", E0, G0, "vs", GODUNOV[Z], (G0 - GODUNOV[Z][1]) / GODUNOV[Z][1], flush=True)
json.dump(res, open(OUT, "w"), indent=1)

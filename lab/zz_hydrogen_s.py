"""zz_hydrogen_s.py - follow-up of zz_hydrogen.py (written AFTER its first run, information only):
why the s levels extrapolate worse (H1 FAIL, 1s 7e-5) than p and d levels (1e-8).

s states are the only ones whose density does not vanish at the nucleus. Near r = 0 the walk's rotation
angle per tick (kappa/r) dt is large (2 rad on the first site), so the lattice error there is not a pure
power series in a^2: it contains a contact term ~ |psi(0)|^2 (expected ~ 1/n^3: 1s : 2s : 3s = 1 : 1/8 : 1/27).
Here: (1) the s levels on finer lattices (theta down to 0.0125), (2) fits with and without a
theta^2 ln(theta) term, (3) the other boundary choice at r = 0 (f(0) = 0 instead of g(0) = 0).
Run: timeout 590 python lab/zz_hydrogen_s.py
"""
import json, os, time
import numpy as np
import scipy.sparse.linalg as spl
from zz_core import RadialWalk, dirac_formula, ALPHA, MC2_EV

OUT = os.path.join(os.path.dirname(__file__), "out")
A0 = 1 / ALPHA
t0 = time.time()
res = {}
ref = {n: 1 - dirac_formula(n, -1) for n in (1, 2)}
E = {1: {}, 2: {}}
for th in (0.2, 0.1, 0.05, 0.025, 0.0125):
    rw = RadialWalk(-1, th, 45 * A0)
    w, v = spl.eigs(rw.matrix(), k=2, sigma=np.exp(-1j * (1 - 1.02 * ref[1]) * th), which="LM", tol=0)
    b = np.sort(1 - (-np.angle(w) / th))[::-1]
    E[1][th], E[2][th] = b[0], b[1]
    print(f"theta={th}: 1s {b[0] * MC2_EV:.7f} eV (err {(b[0] / ref[1] - 1):+.2e}), 2s {b[1] * MC2_EV:.7f} eV (err {(b[1] / ref[2] - 1):+.2e})  ({time.time() - t0:.0f} s)", flush=True)
res["levels"] = {f"{n}s": {str(t): E[n][t] * MC2_EV for t in E[n]} for n in (1, 2)}
for n in (1, 2):
    th = np.array(sorted(E[n])); y = np.array([E[n][t] for t in th]) / ref[n] - 1
    # model A: c2 th^2 + c4 th^4 (+ intercept) on the 3 finest; model B: + th^2 ln th on the 4 finest
    A = np.polyfit(th[:3]**2, y[:3], 2)[-1]
    M = np.column_stack([np.ones(4), th[:4]**2, th[:4]**2 * np.log(th[:4]), th[:4]**4])
    B = np.linalg.solve(M, y[:4])[0]
    res[f"{n}s_extrapolation_rel_err"] = dict(poly_theta2=A, with_log=B)
    print(f"{n}s: relative error of the a->0 extrapolation: polynomial in theta^2 {A:+.1e}; with theta^2 ln theta {B:+.1e}", flush=True)
for refl in (+1.0,):
    rw = RadialWalk(-1, 0.05, 45 * A0, refl0=refl)
    w, v = spl.eigs(rw.matrix(), k=1, sigma=np.exp(-1j * (1 - 1.02 * ref[1]) * 0.05), which="LM", tol=0)
    b = 1 - (-np.angle(w[0]) / 0.05)
    res["boundary_f0"] = b * MC2_EV
    print(f"boundary f(0)=0 (reflection +1), theta=0.05: 1s {b * MC2_EV:.7f} eV vs g(0)=0: {E[1][0.05] * MC2_EV:.7f} eV", flush=True)
# lines with the refined s levels (polynomial extrapolation of the three finest lattices) and the p levels
# of zz_hydrogen.json (which extrapolate to 1e-8)
try:
    from zz_core import HC_EVNM, H_EVS, MP_ME
    MU = 1 / (1 + 1 / MP_ME)
    hj = json.load(open(os.path.join(OUT, "zz_hydrogen.json")))
    e1 = (1 + res["1s_extrapolation_rel_err"]["poly_theta2"]) * ref[1] * MC2_EV
    e2 = (1 + res["2s_extrapolation_rel_err"]["poly_theta2"]) * ref[2] * MC2_EV
    p12, p32 = hj["levels"]["2p1/2"]["extrapolated_eV"], hj["levels"]["2p3/2"]["extrapolated_eV"]
    lya = HC_EVNM / ((e1 - (2 * p32 + p12) / 3) * MU)
    s12 = (e1 - e2) * MU / H_EVS
    res["lines_refined"] = dict(E1s_eV=e1, E2s_eV=e2, lyman_alpha_nm=lya, f1S2S_Hz=s12)
    print(f"refined: 1s {e1:.7f} eV, 2s {e2:.7f} eV -> Lyman-alpha {lya:.5f} nm ({(lya / 121.567 - 1) * 1e6:+.1f} ppm vs 121.567),"
          f" 1S-2S {s12:.9e} Hz ({(s12 - 2466061413187035) / 1e9:+.2f} GHz vs measured: QED Lamb shifts)", flush=True)
except Exception as ex:
    print("refined lines skipped:", ex)
json.dump(res, open(os.path.join(OUT, "zz_hydrogen_s.json"), "w"), indent=1, default=float)
print(f"total {time.time() - t0:.0f} s")

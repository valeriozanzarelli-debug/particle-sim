"""zoomknot_ballcheck.py - POST-HOC diagnostic (not pre-registered): exact hedgehog in a ball of radius R with the
Dirichlet wall f(R) = 0 (continuum, no lattice). Gives the part of the Z1 errors due only to the finite box.
(The Z1 box is a cube of half-side 4.8: the ball of radius 4.8 is inscribed, so this is an upper bound on the effect.)"""
import json, os, numpy as np
from scipy.integrate import solve_bvp
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
res = {}
for R in (4.8, 6.0, 89.6):
    def rhs(r, y):
        f, fp = y; s2 = np.sin(f) ** 2; s2f = np.sin(2 * f)
        return np.vstack([fp, -(2 * r * fp + s2f * fp**2 - s2f - s2 * s2f / r**2) / (r**2 + 2 * s2)])
    r0 = 1e-4
    rg = np.concatenate([np.linspace(r0, 3, 2000), np.geomspace(3.001, R, 2000)]) if R > 3.01 else np.linspace(r0, R, 3000)
    y = np.vstack([np.pi * np.exp(-rg / 0.9) * (1 - rg / R), -np.pi / 0.9 * np.exp(-rg / 0.9)])
    s = solve_bvp(rhs, lambda a, b: np.array([a[0] - np.pi, b[0]]), rg, y, tol=1e-9, max_nodes=300000)
    r = np.linspace(r0, R, 400001) if R < 10 else np.concatenate([np.linspace(r0, 10, 400001), np.geomspace(10.0001, R, 100000)])
    f, fp = s.sol(r); S = np.sin(f) ** 2
    E = 4 * np.pi * np.trapezoid(r**2 * fp**2 + 2 * S + 2 * S * fp**2 + S**2 / r**2, r)
    Lam = (16 * np.pi / 3) * np.trapezoid(r**2 * S * (1 + fp**2 + S / r**2), r)
    a = -float(np.polyfit(r[5:400], fp[5:400], 0)[0]) if False else -float(fp[np.searchsorted(r, 1e-3)])
    p0 = a**2 * (a**2 - 1)
    res[str(R)] = dict(E_over_12pi2=E / (12 * np.pi**2), a=a, p0=p0, err_p0=p0 / 12.236 - 1, Lam=Lam,
                       err_E=E / (12 * np.pi**2) / 1.231446 - 1, err_Lam=Lam / 106.8 - 1, status=int(s.status))
    print(R, res[str(R)])
json.dump(res, open(os.path.join(OUT, "zoomknot_ballcheck.json"), "w"), indent=1)

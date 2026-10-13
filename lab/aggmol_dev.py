"""aggmol_dev.py - DEVELOPMENT exploration (declared as such in the pre-registration, not an exam):
exact parabola kinetic k^2/2 (not the learned rule) + Galerkin field, to size the cell error of the 1s and of H2+ at
a0/2, a0/4, a0/8, and the time per Hamiltonian application. Writes lab/out/aggmol_dev.json."""
import json, os, time, sys
import numpy as np
import aggmol_core as mc
import aggmol_exact as ex

res = {}
JS = os.path.join(mc.OUT, "aggmol_dev.json")
for frac, N in ((2, 32), (4, 64), (8, 128)):
    h = mc.A0 / frac
    T = mc.kinetic_exact(N, h)
    t0 = time.time()
    f = mc.MultiField(N, h, [(0, 0, 0)])
    tf = time.time() - t0
    m = mc.Mol(T, f, sym=(1, 1, 1))
    X, E, r, conv, its = mc.solve_until(m, mc.guess_mol(N, h, [(0, 0, 0)], "1s")[:, None], 400, tag=f"1s a0/{frac}")
    t1 = time.time(); [m.H(X[:, 0].reshape((N,) * 3)) for _ in range(3)]; tH = (time.time() - t1) / 3
    res[f"1s_a0/{frac}"] = dict(E_Eh=E[0] / mc.ALPHA ** 2, rel=E[0] / mc.ALPHA ** 2 / -0.5 - 1, its=its, tH=tH, tfield=tf)
    print(res[f"1s_a0/{frac}"], flush=True)
    for R in (2.0,):
        zc = R * mc.A0 / 2 / h
        f2 = mc.MultiField(N, h, [(0, 0, zc), (0, 0, -zc)])
        m2 = mc.Mol(T, f2, sym=(1, 1, 1))
        X, E, r, conv, its = mc.solve_until(m2, mc.guess_mol(N, h, [(0, 0, zc), (0, 0, -zc)], "1sg")[:, None], 400, tag=f"H2+ a0/{frac}")
        Eex = ex.energy(R)[0]
        res[f"h2p_R{R}_a0/{frac}"] = dict(E_Eh=E[0] / mc.ALPHA ** 2, exact=Eex, rel=E[0] / mc.ALPHA ** 2 / Eex - 1, its=its)
        De = -0.5 - E[0] / mc.ALPHA ** 2 * 0 - 0  # placeholder
        res[f"h2p_R{R}_a0/{frac}"]["De_rule_Eh"] = res[f"1s_a0/{frac}"]["E_Eh"] - (E[0] / mc.ALPHA ** 2 + 1 / R)
        res[f"h2p_R{R}_a0/{frac}"]["De_exact_Eh"] = -0.5 - (Eex + 1 / R)
        print(res[f"h2p_R{R}_a0/{frac}"], flush=True)
    json.dump(res, open(JS, "w"), indent=1)

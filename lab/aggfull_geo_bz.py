"""aggfull_geo_bz.py - geometric idea test (declared post hoc, not pre-registered): how much of the cubic Brillouin zone
does the coarse atom need? Exact k^2/2 kinetic + alias-free Galerkin field, 1s and 2p... on a0/8 cells, box +-8 a0;
modes kept: the full cube, or only the ball |k| < k_c (fraction of modes (4 pi/3)(k_c h/2pi)^3 ...). If the ball inscribed in
the cube (52 % of the modes) gives almost the same energy, a lattice whose zone is rounder than a cube (body-centred cubic:
zone = rhombic dodecahedron, or face-centred: truncated octahedron) needs ~1/0.52-0.68 fewer cells for the same accuracy.
Writes lab/out/aggfull_geo_bz.json"""
import json, os, sys, time
import numpy as np
import aggfull_core as ac, aggfull_atom as aa
N, B = 128, 512
h = B * ac.TH_A
f = ac.GalerkinField(N, h)
k2 = aa.kgrid2(N, h)
kmax = np.pi / h
out = {}
for name, kc in (("cube", None), ("ball_inscribed", kmax), ("ball_0.9", 0.9 * kmax), ("ball_0.8", 0.8 * kmax), ("ball_0.7", 0.7 * kmax)):
    T = k2 / 2
    frac = 1.0
    if kc is not None:
        T = np.where(k2 < kc ** 2, T, 50.0)          # modes outside the ball removed (penalty 50 mc^2)
        frac = float(np.mean(k2 < kc ** 2))
    at = aa.Atom(T, f, gamma0=1.0, real=True, sym="even")
    v, r, it = aa.solve(at, aa.guess(f.x, "1s"), maxiter=40)
    E, res = at.energy(v)
    out[name] = dict(mode_fraction=frac, E_eV=E * ac.MC2, rel=E * ac.MC2 / (-ac.RY) - 1, residual=res)
    print(name, f"modes {frac:.3f}  1s rel {E * ac.MC2 / (-ac.RY) - 1:+.3e}", flush=True)
json.dump(out, open(os.path.join(ac.OUT, "aggfull_geo_bz.json"), "w"), indent=1)

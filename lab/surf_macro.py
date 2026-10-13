"""surf_macro.py - EXPLORATION (not pre-registered), idea "surfisti agglomerati": the exact flow works on ANY block
level, because summing the fine face fluxes over the faces of a 2^k block gives the exact flux of the block (discrete
Stokes). Surfers at level k (cells of 2^k x 2^k x 2^k) are exactly equivariant for the block measure, cost 8^k less
grid work, and see the two sheets of the zig-zag merged. Test on the listened 1s (two-tick flux), 400 ticks:
Stokes check, chi2 on the blocks, relaxation of a |psi|^4 start in equal-probability radial shells, r drift.
Writes lab/out/surf_macro.json."""
import json
import numpy as np
from surf_core import *
from surf_run import setup_3d, A3


def coarse(rho, b):
    m = rho.shape[0] // b
    return rho.reshape(m, b, m, b, m, b).sum((1, 3, 5))


def coarse_flux(F, j, b):
    sl = [slice(None)] * 3
    sl[j] = slice(b - 1, None, b)
    G = F[tuple(sl)]                                   # fine links on the +j faces of the blocks
    shp = []
    for a in range(3):
        shp += [G.shape[a]] if a == j else [G.shape[a] // b, b]
    G = G.reshape(shp)
    inblock = []
    pos = 0
    for a in range(3):
        if a == j:
            pos += 1
        else:
            inblock.append(pos + 1); pos += 2
    return G.sum(axis=tuple(inblock))


rng = np.random.default_rng(SEED + 11)
W, ps = setup_3d("1s")
rho = density(ps)
T = 400
L = {}
for b in (1, 2):
    rc = coarse(rho, b)
    m = rc.shape[0]
    h = A3 * b
    xc = (np.arange(m) - m / 2 + 0.5) * h
    Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing="ij")
    Rr = np.sqrt(Xc**2 + Yc**2 + Zc**2).ravel()
    o = np.argsort(Rr); C = np.cumsum(rc.ravel()[o]); C /= C[-1]
    edges = np.concatenate([[0], np.interp(np.linspace(0, 1, 21)[1:-1], C, Rr[o]), [1e9]])
    c, x = born_sample(rng, rc, 100000)
    cn, xn = born_sample(rng, rc**2, 50000)
    L[b] = dict(c=c, x=x, cn=cn, xn=xn, m=m, h=h, edges=edges, z=[], stokes=0.0)
def rad(c, x, m, h):
    p = (c + x - m / 2) * h
    return np.sqrt((p**2).sum(1))
def H(r, e):
    cnt = np.histogram(r, bins=e)[0] / len(r); mm = cnt > 0
    return float(np.sum(cnt[mm] * np.log(cnt[mm] * (len(e) - 1))))
for b, R in L.items():
    R["r0"] = rad(R["c"], R["x"], R["m"], R["h"]); R["H0"] = H(rad(R["cn"], R["xn"], R["m"], R["h"]), R["edges"])
for t in range(0, T, 2):
    r0 = density(ps); F2 = [0, 0, 0]
    for _ in range(2):
        ps, subs = W.tick(ps)
        for (j, A, B) in subs:
            F2[j] = F2[j] + net_flux(j, A, B)
    r1 = density(ps)
    for b, R in L.items():
        a0, a1 = coarse(r0, b), coarse(r1, b)
        Fk = [coarse_flux(F2[j], j, b) for j in range(3)]
        div = sum(Fk[j] - np.roll(Fk[j], 1, axis=j) for j in range(3))
        R["stokes"] = max(R["stokes"], float(np.max(np.abs(a0 - div - a1))))
        flow_tick(R["c"], R["x"], a0, a1, Fk)
        flow_tick(R["cn"], R["xn"], a0, a1, Fk)
        if t % 40 == 0:
            R["z"].append(chi2_z(R["c"], a1, 1)[0])
out = {}
for b, R in L.items():
    dr = rad(R["c"], R["x"], R["m"], R["h"]) - R["r0"]
    out[f"block_{b}"] = dict(stokes_err=R["stokes"], max_abs_z=float(np.max(np.abs(R["z"]))),
                             neq_H_ratio=H(rad(R["cn"], R["xn"], R["m"], R["h"]), R["edges"]) / R["H0"],
                             rms_dr_over_r=float(np.sqrt(np.mean(dr**2)) / np.mean(R["r0"])))
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(OUT, "surf_macro.json"), "w"), indent=1)

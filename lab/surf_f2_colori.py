"""surf_f2_colori.py - POST-HOC (declared B5 in surf_f2_amend.json), after the SOTTORETICOLO warning (master 178743d):
did E5 measure one colour (sublattice) or the sum of the 4? For the E5 state and the exact levels: weight, <Sigma_z> and
the two-tick-flux v_phi of each colour c = ((i+j)%2, (j+k)%2) evolved ALONE (colours never exchange amplitude, and the
link flux A - B is a sum of |amplitude|^2 at single cells, so the full flux is exactly the sum of the colour fluxes).
-> lab/out/surf_f2_colori.json"""
import json
import os
import numpy as np
from surf_core import OUT
from surf_f2_valli import sigma_z, walk, CK, N3, X, Y, Z, R, ZA
from surf_core import net_flux


def flux_vphi(W, psi):
    """as surf_f2_valli.flux_vphi but density-weighted as a ratio of sums (sum J_phi R/rc / sum rho over the same
    region), so it is defined on one colour (where rho vanishes on 3/4 of the cells) and linear in the flux."""
    r0 = np.sum(np.abs(psi)**2, 0); F2 = [0, 0, 0]; W.tick_count = 0
    for _ in range(2):
        psi, subs = W.tick(psi)
        for (j, A, B) in subs:
            F2[j] = F2[j] + net_flux(j, A, B)
    W.tick_count = 0
    Fc = [0.5 * (F2[j] + np.roll(F2[j], 1, axis=j)) / 2 for j in range(3)]
    rc = np.sqrt(X**2 + Y**2)
    m = (R > 1) & (R < 6) & (np.abs(Z) < 0.5 * R)
    Jphi = (-Y * Fc[0] + X * Fc[1]) / rc
    return float(np.sum((Jphi * R / rc)[m]) / np.sum(r0[m]) / ZA)

i, j, k = np.meshgrid(*(np.arange(N3),) * 3, indexing="ij")
COL = {f"{a}{b}": ((i + j) % 2 == a) & ((j + k) % 2 == b) for a in (0, 1) for b in (0, 1)}
names = {"E5_listened": "state_1s.npy", "0.89930": "f2_level_0.89930.npy", "0.89960": "f2_level_0.89960.npy"}
W = walk()
REG = (R > 1) & (R < 6) & (np.abs(Z) < 0.5 * R)
res = {}
for nm, fn in names.items():
    p = np.load(os.path.join(CK, fn)); p = p / np.linalg.norm(p)
    d = dict(total=dict(sigma_z=sigma_z(p), vphi=flux_vphi(W, p)), colours={})
    for c, m in COL.items():
        pc = p * m[None]
        w = float(np.sum(np.abs(pc)**2))
        pcn = pc / np.sqrt(w)
        # the colour evolved alone: check it stays on its colour after two ticks
        q = pcn.copy(); W.tick_count = 0
        for _ in range(2):
            q, _s = W.tick(q)
        W.tick_count = 0
        d["colours"][c] = dict(weight=w, sigma_z=sigma_z(pcn), vphi=flux_vphi(W, pcn),
                               leak_after_2_ticks=float(np.sum(np.abs(q * (~m)[None])**2)),
                               rho_region=float(np.sum(np.sum(np.abs(pcn)**2, 0)[REG])))
    d["weighted_vphi_sum"] = sum(v["weight"] * v["vphi"] * v["rho_region"] for v in d["colours"].values()) / float(np.sum(np.sum(np.abs(p)**2, 0)[REG]))
    res[nm] = d
    print(nm, {kk: round(v, 4) for kk, v in d["total"].items()}, "sum_w_vphi", round(d["weighted_vphi_sum"], 4),
          {c: {q: round(v, 4) for q, v in x.items()} for c, x in d["colours"].items()}, flush=True)
json.dump(res, open(os.path.join(OUT, "surf_f2_colori.json"), "w"), indent=1, default=float)

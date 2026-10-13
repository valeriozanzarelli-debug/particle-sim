"""surf_f2_posthoc.py - POST-HOC diagnostic (declared, not pre-registered) of part 1: is the extra mixing (M = 1 relaxes
to ~55 %, tau gets SHORTER on the finer lattice) made by circulation at the scale of one cell? Same states and start
(|psi|^4) as surf_f2_valentini.py, but the surfers move on blocks of b x b cells with the exact block flux (sum of the
fine face fluxes, exact by discrete Stokes; idea 3 of SURFISTI_ESATTI). Block side in physical units: 0.5 and 1.0.
usage: python surf_f2_posthoc.py <case 0..5>. Writes lab/out/surf_f2_posthoc.json."""
import json
import numpy as np
from surf_core import *
from pre_core import LatticeWalk, modes_state, hbar, fit_tau
from surf_f2_valentini import tick_flux, CFG, NS


def block(r, b):
    m = r.shape[0] // b
    return r.reshape(m, b, m, b).sum((1, 3))


def block_flux(F, j, b):
    sl = [slice(None)] * 2
    sl[j] = slice(b - 1, None, b)
    G = F[tuple(sl)]
    if j == 0:
        return G.reshape(G.shape[0], G.shape[1] // b, b).sum(2)
    return G.reshape(G.shape[0] // b, b, G.shape[1]).sum(1)


import sys
CASES = ((48, 1, 2), (48, 16, 2), (96, 1, 2), (96, 1, 4), (96, 16, 2), (96, 16, 4))
jf = os.path.join(OUT, "surf_f2_posthoc.json")
out = json.load(open(jf)) if os.path.exists(jf) else {}
for n, M, b in [CASES[int(sys.argv[1])]]:
    c = CFG[n]
    wk = LatticeWalk(n, c["theta"], periodic=False)
    psi, om = modes_state(wk, M, omega0=2 * c["theta"])
    rng = np.random.default_rng(SEED + 100 * M + n)
    rho = density(psi)
    cl, xi = born_sample(rng, block(rho**2, b) if False else block(rho, b)**2, NS)
    cellH = c["cell"] // b
    K = (n // c["cell"])**2
    floor = (K - 1) / (2 * NS)
    ts, H = [], []
    for t in range(c["ticks"] + 1):
        if t % c["every"] == 0:
            ts.append(t); H.append(hbar(cl, block(density(psi), b), cellH)[0])
        if t == c["ticks"]:
            break
        r0 = density(psi)
        psi, F = tick_flux(wk, psi)
        flow_tick(cl, xi, block(r0, b), block(density(psi), b), [block_flux(F[j], j, b) for j in range(2)])
    tau = fit_tau(ts, H, floor)
    key = f"n{n}_M{M}_b{b}"
    out[key] = dict(block_phys=b * c["theta"], ratio=(H[-1] - floor) / (H[0] - floor), tau_ticks=tau,
                    tau_phys=tau * c["theta"], t=ts, H=H, floor=floor)
    print(key, round(out[key]["ratio"], 3), round(out[key]["tau_phys"], 1), flush=True)
json.dump(out, open(jf, "w"), indent=1, default=float)

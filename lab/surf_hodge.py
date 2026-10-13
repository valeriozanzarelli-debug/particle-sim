"""surf_hodge.py - EXPLORATION (not pre-registered), for the section "Idee geometriche e visive".
The two-tick net flux of the walk is a discrete 2-form on the faces of the cubic lattice. Hodge/Helmholtz split on the
lattice:  F = -grad(phi) + F_rot,  L phi = -div F  (FFT, periodic):  the gradient part carries ALL the change of the
density (breathing), the rotational part is divergence-free (circulation: for the 1s the spin current of Bohm-Dirac,
Gordon decomposition). Question: which fraction of the surfer motion is circulation, in the stationary 1s and along the
beat of 1s + 2p? (colour idea: surfers coloured by the local share of rotation vs gradient).
Writes lab/out/surf_hodge.json."""
import json
import numpy as np
from surf_core import *
from surf_run import setup_3d, N3, A3, ZA, grid3


def hodge(F):
    n = F[0].shape[0]
    div = sum(F[j] - np.roll(F[j], 1, axis=j) for j in range(3))
    k = 2 * np.pi * np.fft.fftfreq(n)
    lam = sum((2 * np.cos(kk) - 2) for kk in np.meshgrid(k, k, k, indexing="ij"))
    lam[0, 0, 0] = 1.0
    phi = np.fft.ifftn(np.fft.fftn(div) / lam).real              # L phi = div F
    G = [np.roll(phi, -1, axis=j) - phi for j in range(3)]       # gradient part, div G = div F
    Rr = [F[j] - G[j] for j in range(3)]
    return G, Rr, div


out = {}
for state, nticks in (("1s", 8), ("1s2p", 168)):
    W, psi = setup_3d(state)
    X, Y, Z = grid3()
    rows = []
    for t in range(0, nticks, 2):
        F2 = [0, 0, 0]
        r0 = density(psi)
        for _ in range(2):
            psi, subs = W.tick(psi)
            for (j, A, B) in subs:
                F2[j] = F2[j] + net_flux(j, A, B)
        G, Rr, div = hodge(F2)
        nG = sum(np.sum(g**2) for g in G); nR = sum(np.sum(r**2) for r in Rr)
        # azimuthal share of the rotational part (circulation around the spin axis z), cell-centred
        Rc = [0.5 * (Rr[j] + np.roll(Rr[j], 1, axis=j)) for j in range(3)]
        rc = np.sqrt(X**2 + Y**2) + 1e-9
        rphi = (-Y * Rc[0] + X * Rc[1]) / rc
        az = float(np.sum(rphi**2) / max(sum(np.sum(c**2) for c in Rc), 1e-300))
        dz = float(np.sum(density(psi) * Z) / np.sum(density(psi)))
        rows.append(dict(t=t + 2, grad_share=float(nG / (nG + nR)), azimuthal_share_of_rot=az, z_mean=dz))
    out[state] = rows
    g = [r["grad_share"] for r in rows]
    print(state, "gradient share of |F|^2: min %.3e max %.3e ; azimuthal share of rotation: %.3f" %
          (min(g), max(g), np.mean([r["azimuthal_share_of_rot"] for r in rows])), flush=True)
json.dump(out, open(os.path.join(OUT, "surf_hodge.json"), "w"), indent=1)

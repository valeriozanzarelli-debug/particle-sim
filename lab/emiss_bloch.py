"""emiss_bloch.py - the light rule's OWN Larmor law, from its Bloch matrices (no time stepping).
The two-tick light walk W(k) is a real rotation by angle Phi(k) about an axis n(k): its transverse modes have
frequency +-Phi/(2b) and the "sin^2 theta" of a dipole along z becomes 1 - (n.z)^2.  Golden rule:
  P_lattice / (P_Larmor * W2) = Int dOmega [k*^2 / (dPhi/dk / 2b)] (1 - (n.z)^2) |W(k*)|^2 / Int dOmega w^2 sin^2 |W(w n)|^2
with k* on each ray solving Phi(k* n) = 2 b w.  W = discrete Fourier factor of the Gaussian average (sigma cells)."""
import numpy as np
from zz_core import bloch_batch
from emiss_core import Light3D


def ratio(cells_per_lambda, sigma=1.0, nth=48, nph=96):
    b = 1.0
    w = 2 * np.pi / cells_per_lambda
    L = Light3D(9, b, sponge=0, sigma=sigma)
    o, wt = L.offs, L.w
    from scipy.special import roots_legendre
    u, wu = roots_legendre(nth)
    ph = (np.arange(nph) + 0.5) * 2 * np.pi / nph
    num = den = 0.0
    ks = np.linspace(0.5 * w, 1.6 * w, 400)
    for iu, cu in enumerate(u):
        s = np.sqrt(1 - cu * cu)
        for p in ph:
            nhat = np.array([s * np.cos(p), s * np.sin(p), cu])
            K = ks[:, None] * nhat[None]
            R = bloch_batch(K, spin=1).real
            Phi = np.arccos(np.clip((np.trace(R, axis1=1, axis2=2) - 1) / 2, -1, 1))
            i = np.searchsorted(Phi, 2 * b * w)
            f = (2 * b * w - Phi[i - 1]) / (Phi[i] - Phi[i - 1])
            ks_ = ks[i - 1] + f * (ks[i] - ks[i - 1])
            dPhi = (Phi[i] - Phi[i - 1]) / (ks[i] - ks[i - 1])
            Rk = bloch_batch(ks_ * nhat[None], spin=1).real[0]
            ax = np.array([Rk[2, 1] - Rk[1, 2], Rk[0, 2] - Rk[2, 0], Rk[1, 0] - Rk[0, 1]]); ax /= np.linalg.norm(ax)
            Wl = np.sum(wt * np.exp(1j * ks_ * (o @ nhat)))
            Wc = np.sum(wt * np.exp(1j * w * (o @ nhat)))
            num += wu[iu] * ks_**2 / (dPhi / (2 * b)) * (1 - ax[2]**2) * abs(Wl)**2
            den += wu[iu] * w**2 * (1 - cu * cu) * abs(Wc)**2
    return num / den


if __name__ == "__main__":
    for c in (20, 30, 40):
        print(f"{c} cells per wavelength: P_lattice / (Larmor W2) = {ratio(c):.5f}")

"""kk_speed_core.py - the KK zig-zag walk of kk_core.py with ONE change: the hop along the hidden
circle w has its own speed c_w (c = 1 along x, y).

Hop along w at speed c_w: in one tick the L/R parts slide around the hidden circle by -/+ c_w cells
(exact band-limited translation: hidden mode n gets the phase exp(-+ i c_w n / R)). c_w = 1 is the
rule of kk_core.py exactly. Nothing else changes (same flip, same alternated axis order).

Dynamic "time twist" (task 2): at the end of every tick every ball slides around the hidden circle by
phi(x) (a g_tw component of the geometry = the electric potential A_0 = phi / R written as a slide).
Mode n then picks up the phase exp(-i n phi / R): potential energy (n / R) phi. The field phi is NOT
a background: it is computed every tick from the charge of the OTHER packet with its own stiffness K
(Gauss law in 1+1 D, exact: K phi'' = -rho, rho = sum_n (n / R) |psi_n|^2).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from zz_core import BASIS_HALF, order_of_tick
import kk_core as K


class KKWalkCw(K.KKWalk):
    def __init__(self, nx, ny, N, theta, cw=1.0, **kw):
        super().__init__(nx, ny, N, theta, **kw)
        self.cw = float(cw)
        self.phi = None                       # time-twist slide per tick, shape (nx, ny) or None

    def _axis(self, blk, j, h):
        if j != 2 or self.cw == 1.0:
            return super()._axis(blk, j, h)
        lam, P = BASIS_HALF[j]
        c = np.tensordot(P.conj().T, blk, axes=(1, 0))
        one = np.ones(self.shape[:2])
        for q in range(2):
            c[q] = self._hshift(c[q], h * lam[q] * self.cw * one)
        return np.tensordot(P, c, axes=(1, 0))

    def tick(self, psi):
        out = super().tick(psi)
        if self.phi is not None:
            for q in range(4):
                out[q] = self._hshift(out[q], self.phi)
        return out


def charge_density(psi, N):
    """rho(x) = sum_n (n/R) |psi_n(x)|^2  (coupling units: the phase per unit slide is n/R)"""
    R = N / (2 * np.pi)
    F = np.fft.fft(psi, axis=-1)
    w = np.sum(np.abs(F)**2, axis=0) / N            # (nx, ny, N) weight of each hidden mode
    return np.sum(w * (K.nvals(N) / R)[None, None, :], axis=-1)


def gauss_1d(rho, Kst):
    """phi(x) solving K phi'' = -rho on an OPEN line (no periodic images): phi = -(1/2K) sum rho |x-x'|"""
    L = rho.size
    G = np.abs(np.arange(-L + 1, L))                  # |x - x'|
    full = np.convolve(rho, G)                         # length 3L-2, centre aligned at index L-1
    return -full[L - 1:2 * L - 1] / (2 * Kst)


def mean_x_k(psi):
    """<x> and <k_x> of a packet psi (4, nx, 1, N) on the x ring"""
    nx = psi.shape[1]
    rx = np.sum(np.abs(psi)**2, axis=(0, 2, 3))
    x = np.arange(nx)
    F = np.fft.fft(psi, axis=1)
    rk = np.sum(np.abs(F)**2, axis=(0, 2, 3))
    k = 2 * np.pi * np.fft.fftfreq(nx)
    return float(np.sum(x * rx) / rx.sum()), float(np.sum(k * rk) / rk.sum())


def packet_cw(nx, ny, N, theta, n, cw=1.0, k0=(0.0, 0.0), sigma=5.0, center=None, spinor=(1, 0, 1, 0)):
    """positive-energy packet with hidden mode n for the walk with hop speed cw along w.
    Same construction as kk_core.packet; the Bloch matrix sees the hidden phase cw * 2 pi n / N per hop."""
    from zz_core import bloch_batch
    kx = 2 * np.pi * np.fft.fftfreq(nx)
    ky = 2 * np.pi * np.fft.fftfreq(ny)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    kw = 2 * np.pi * n / N
    cx, cy = center if center is not None else (nx // 2, ny // 2)
    dkx = (KX - k0[0] + np.pi) % (2 * np.pi) - np.pi
    dky = (KY - k0[1] + np.pi) % (2 * np.pi) - np.pi
    g = np.exp(-sigma**2 * (dkx**2 + (dky**2 if ny > 1 else 0 * dky)) - 1j * (KX * cx + KY * cy))
    KA = np.stack([KX.ravel(), KY.ravel(), np.full(KX.size, cw * kw)], axis=1)
    W = bloch_batch(KA, theta)
    lam, V = np.linalg.eig(W)
    om = -np.angle(lam) / 2
    s0 = np.asarray(spinor, complex) / np.linalg.norm(spinor)
    coef = np.empty((KA.shape[0], 4), complex)
    for i in range(KA.shape[0]):
        Vp = V[i][:, om[i] > 0]
        Pp = Vp @ np.linalg.solve(Vp.conj().T @ Vp, Vp.conj().T)
        coef[i] = Pp @ s0
    coef = coef.reshape(nx, ny, 4).transpose(2, 0, 1) * g[None]
    psi2 = np.fft.ifft2(coef, axes=(1, 2))
    w = np.arange(N)
    psi = psi2[..., None] * np.exp(1j * kw * w)[None, None, None, :]
    return psi / np.sqrt(np.sum(np.abs(psi)**2))


def rest_energy(N, theta, n, cw, ticks=600):
    """rest energy of hidden mode n measured on the walk (phase of <psi0|psi(t)> per tick, 2-tick period)"""
    psi = packet_cw(1, 1, N, theta, n, cw=cw, sigma=0.0)
    W = KKWalkCw(1, 1, N, theta, cw=cw)
    p0 = psi.copy(); ov = np.empty(ticks, complex)
    for t in range(ticks):
        ov[t] = np.vdot(p0, psi); psi = W.tick(psi)
    return float(-np.angle(np.mean(ov[2::2] / ov[:-2:2])) / 2)

"""senzadoppi_real.py - real-space runner for the rules of senzadoppi_ (the same factors as the Bloch matrices of
senzadoppi_core / senzadoppi_r1): axis steps for L and R in a given order, coin (mass flip), twisted Wilson link
rotations, Coulomb (or any) potential as a phase applied together with the coin.  Periodic box n^3, psi (4, n, n, n)."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import numpy as np
import zz_core as zz

BH = zz.BASIS_HALF


def axis_shift(blk, j, h):
    """component with sigma_j eigenvalue lam moves by h*lam along axis j (blk shape (2, n, n, n))"""
    lam, P = BH[j]
    c = np.tensordot(P.conj().T, blk, axes=(1, 0))
    for q in range(2):
        c[q] = np.roll(c[q], h * lam[q], axis=j)
    return np.tensordot(P, c, axes=(1, 0))


def rot4(psi, M):
    return np.tensordot(M, psi, axes=(1, 0))


def link(psi, v, a, G):
    """S_v^dag e^{-i a G} S_v e^{+i a G}; S_v moves the R part by +v (R(x) <- R(x - v))"""
    ea = np.cos(a) * np.eye(4) - 1j * np.sin(a) * G
    eb = np.cos(a) * np.eye(4) + 1j * np.sin(a) * G
    psi = rot4(psi, eb)
    v = tuple(int(x) for x in v)
    psi[2:] = np.roll(psi[2:], v, axis=(1, 2, 3))
    psi = rot4(psi, ea)
    psi[2:] = np.roll(psi[2:], tuple(-x for x in v), axis=(1, 2, 3))
    return psi


class RealRule:
    """ticks: list of (orderR, orderL); mode 'each' (coin theta each tick, Wilson around each coin, mirrored on odd
    ticks) or 'once' (one coin T*theta per period). h, vecs: Wilson. V: potential array (n,n,n) or None, applied as
    exp(-i V dt) with dt = duration between coins, at the coin."""

    def __init__(self, ticks, theta, mode="each", h=0.0, vecs=(), V=None):
        self.ticks, self.theta, self.mode, self.h, self.vecs = ticks, theta, mode, h, list(vecs)
        self.T = len(ticks)
        self.V = V
        YCHI = np.kron(zz.SY, zz.I2)
        self.G = YCHI
        XCHI = np.kron(zz.SX, zz.I2)
        ang = theta if mode == "each" else self.T * theta
        self.C = np.cos(ang) * np.eye(4) - 1j * np.sin(ang) * XCHI
        dt = 1 if mode == "each" else self.T
        self.ph = None if V is None else np.exp(-1j * V * dt)

    def _coin_block(self, psi, n):
        if self.h:
            v = self.vecs if n % 2 == 0 or self.mode == "once" else self.vecs[::-1]
            for w in v:
                psi = link(psi, w, self.h / 4, self.G)
        psi = rot4(psi, self.C)
        if self.ph is not None:
            psi = psi * self.ph[None]
        if self.h:
            v = self.vecs if n % 2 == 0 or self.mode == "once" else self.vecs[::-1]
            for w in v[::-1]:
                psi = link(psi, w, self.h / 4, self.G)
        return psi

    def period(self, psi):
        psi = psi.copy()
        for n, (oR, oL) in enumerate(self.ticks):
            L, R = psi[:2], psi[2:]
            for j in oR:
                R = axis_shift(R, j, +1)
            for j in oL:
                L = axis_shift(L, j, -1)
            psi = np.concatenate([L, R])
            if self.mode == "each":
                psi = self._coin_block(psi, n)
        if self.mode == "once":
            psi = self._coin_block(psi, 0)
        return psi

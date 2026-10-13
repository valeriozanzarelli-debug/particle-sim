"""qg_pacchetti.py - PART 2 initial conditions: gluon wave packets of a fixed colour direction.

A packet = a solution of the LINEAR lattice wave equation (transverse, Gauss law exact in the abelian sense),
embedded along one su(3) generator T^a:  U_i = exp(i A_i T^a),  E_i = e_i T^a.
Inside one packet everything commutes, so the non-abelian Gauss law is exactly satisfied too. Two packets of different
colour are put far apart (Gaussian tails overlap ~ e^{-(sep/2w)^2 ...}: the Gauss residual this leaves is measured).
"""
import numpy as np

from qg_su3 import TA, expi, mm


def packet(shape, centre, k0, pol, sig_par, sig_perp, amp, direction):
    """vector potential A_i (3,*shape) and its time derivative (= our E_i) for a packet moving along +x (direction=+1)
    or -x (direction=-1), carrier wave number k0 along x, polarisation unit vector pol (transverse)."""
    X, Y, Z = shape
    x = np.arange(X)[:, None, None] - centre[0]
    y = np.arange(Y)[None, :, None] - centre[1]
    z = np.arange(Z)[None, None, :] - centre[2]
    env = np.exp(-x ** 2 / (2 * sig_par ** 2) - (y ** 2 + z ** 2) / (2 * sig_perp ** 2))
    A0 = np.stack([amp * p * env * np.cos(k0 * x) for p in pol])
    kx = 2 * np.pi * np.fft.fftfreq(X)[:, None, None]
    ky = 2 * np.pi * np.fft.fftfreq(Y)[None, :, None]
    kz = 2 * np.pi * np.fft.fftfreq(Z)[None, None, :]
    K = [np.exp(1j * kx) - 1, np.exp(1j * ky) - 1, np.exp(1j * kz) - 1]       # forward difference (links)
    K = [np.broadcast_to(k, shape) for k in K]
    K2 = sum(np.abs(k) ** 2 for k in K)
    a = np.stack([np.fft.fftn(A0[i]) for i in range(3)])
    # transverse projection: conj(K) . a = 0  (lattice Coulomb gauge; keeps only the physical, radiating part)
    div = sum(np.conj(K[i]) * a[i] for i in range(3))
    safe = np.where(K2 > 0, K2, 1.0)
    a = np.stack([a[i] - K[i] * div / safe for i in range(3)])
    a[:, 0, 0, 0] = 0
    w = np.sqrt(K2)
    sg = np.sign(np.broadcast_to(np.fft.fftfreq(X)[:, None, None], shape)) * direction
    adot = -1j * w * sg * a                      # right mover for direction=+1: A(x - t)
    A = np.stack([np.fft.ifftn(a[i]).real for i in range(3)])
    Ad = np.stack([np.fft.ifftn(adot[i]).real for i in range(3)])
    return A, Ad


def build(ym, packets):
    """packets: list of dicts (centre, k0, pol, sig_par, sig_perp, amp, direction, gen). Sets ym.U, ym.E."""
    U = ym.U
    E = ym.E
    E[:] = 0
    for p in packets:
        A, Ad = packet(ym.shape, p["centre"], p["k0"], p["pol"], p["sig_par"], p["sig_perp"], p["amp"], p["direction"])
        T = TA[p["gen"]]
        for i in range(3):
            H = np.einsum("...,jk->jk...", A[i], T)
            U[i] = mm(expi(H, order=12), U[i])
            # our E is dA/dt (dU/dt = i E U): for one packet U = exp(iA T) -> E = (dA/dt) T exactly
            E[i] += np.einsum("...,jk->jk...", Ad[i], T)
    return ym

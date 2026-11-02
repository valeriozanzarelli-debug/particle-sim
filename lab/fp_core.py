"""fp_core.py - FOTONE A PALLINE: the zig-zag light rule in k-space + energy balls that ride the field.

The rule is zz_core.Walk('light') (frozen 2026-10-09). On a periodic box a roll by +-1 cell along axis j is
exactly a multiplication by exp(-+i k_j) in Fourier space, so one axis of the rule is the 3x3 matrix
    M_j(k_j) = P_j diag(exp(-i lam_q k_j a)) P_j^dag      (lam = +1, 0, -1 the eigenvalues of S_j)
and a tick is M_z M_y M_x (even ticks) or M_x M_y M_z (odd ticks). Gate G0 checks it against Walk.tick.

Balls (the drawing): seeded ONCE with density ~ u = (E^2 + B^2)/2, then moved by the energy velocity
    v = 2 E x B / (E^2 + B^2)            (c = 1, |v| <= 1)
Because du/dt + div(u v) = 0 (Poynting), a cloud seeded ~ u stays ~ u: no re-birth is needed (exam F1).
Units: lattice cell a = 1, tick = 1, c = 1. F = E + iB.
"""
import numpy as np
import scipy.fft as sfft
from scipy.ndimage import map_coordinates
from zz_core import BASIS_ONE, order_of_tick

WORKERS = 4


def axis_mats(n, j):
    """M_j for every k_j on an axis of n sites: shape (n, 3, 3)"""
    lam, P = BASIS_ONE[j]
    k = 2 * np.pi * sfft.fftfreq(n)
    ph = np.exp(-1j * np.outer(k, lam))                      # (n, 3): roll by +lam
    return np.einsum('iq,kq,jq->kij', P, ph, P.conj())


class LightK:
    """the zig-zag light rule, applied in Fourier space (identical to the lattice roll on a periodic box)"""

    def __init__(self, shape):
        self.shape = tuple(shape)
        self.M = [axis_mats(shape[j], j).astype(np.complex128) for j in range(3)]
        self.tick_count = 0

    def _apply(self, Fk, j):
        M = self.M[j]
        sl = [None, None, None]
        sl[j] = slice(None)
        out = np.empty_like(Fk)
        for i in range(3):
            acc = 0
            for l in range(3):
                acc = acc + M[(slice(None), i, l)][tuple(sl)] * Fk[l]
            out[i] = acc
        return out

    def tick(self, Fk):
        for j in order_of_tick(self.tick_count):
            Fk = self._apply(Fk, j)
        self.tick_count += 1
        return Fk


def to_k(F):
    return sfft.fftn(F, axes=(1, 2, 3), workers=WORKERS)


def to_x(Fk):
    return sfft.ifftn(Fk, axes=(1, 2, 3), workers=WORKERS)


def pulse(shape, lam, w0, sig, pol, center=None):
    """Gaussian pulse flying along +z, built wave by wave from Maxwell plane waves: E perp q, B = q^ x E.
    pol: 'circ+', 'circ-', 'lin_x'. Returns F = E + iB on the lattice (complex, shape (3,)+shape)."""
    nx, ny, nz = shape
    if center is None:
        center = (nx / 2, ny / 2, nz / 2)
    qx = 2 * np.pi * sfft.fftfreq(nx)[:, None, None]
    qy = 2 * np.pi * sfft.fftfreq(ny)[None, :, None]
    qz = 2 * np.pi * sfft.fftfreq(nz)[None, None, :]
    k0 = 2 * np.pi / lam
    A = np.exp(-(qx**2 + qy**2) * w0**2 / 4 - (qz - k0)**2 * sig**2 / 2) * (qz > 0)
    A = A * np.exp(-1j * (qx * center[0] + qy * center[1] + qz * center[2]))   # place the focus
    qn = np.sqrt(qx**2 + qy**2 + qz**2)
    qn[qn == 0] = 1
    qh = [np.broadcast_to(qx / qn, A.shape), np.broadcast_to(qy / qn, A.shape), np.broadcast_to(qz / qn, A.shape)]
    s = {'circ+': 1, 'circ-': -1, 'lin_x': 0}[pol]
    e0 = np.array([1, 1j * s, 0]) / (np.sqrt(2) if s else 1)
    dot = qh[0] * e0[0] + qh[1] * e0[1] + qh[2] * e0[2]
    Ek = np.array([A * (e0[i] - qh[i] * dot) for i in range(3)])
    Bk = np.array([qh[1] * Ek[2] - qh[2] * Ek[1], qh[2] * Ek[0] - qh[0] * Ek[2], qh[0] * Ek[1] - qh[1] * Ek[0]])
    E = sfft.ifftn(Ek, axes=(1, 2, 3), workers=WORKERS).real
    B = sfft.ifftn(Bk, axes=(1, 2, 3), workers=WORKERS).real
    nrm = np.sqrt(((E**2 + B**2).sum()) / 2)
    return (E + 1j * B) / nrm


def energy(F):
    return 0.5 * (np.abs(F)**2).sum(0)          # (E^2 + B^2)/2


def poynting(F):
    E, B = F.real, F.imag
    return np.array([E[1] * B[2] - E[2] * B[1], E[2] * B[0] - E[0] * B[2], E[0] * B[1] - E[1] * B[0]])


def seed(u, n, rng):
    """n balls with density ~ u (cell chosen ~ u, uniform jitter inside the cell centred on the site)"""
    p = u.ravel() / u.sum()
    idx = rng.choice(p.size, size=n, p=p)
    ijk = np.array(np.unravel_index(idx, u.shape), float)
    return ijk + rng.random((3, n)) - 0.5


def sample(F, pos):
    """trilinear E, B at ball positions (periodic)"""
    out = np.empty((6, pos.shape[1]))
    E, B = F.real, F.imag
    for i in range(3):
        out[i] = map_coordinates(E[i], pos, order=1, mode='grid-wrap')
        out[3 + i] = map_coordinates(B[i], pos, order=1, mode='grid-wrap')
    return out[:3], out[3:]


def vel(E, B):
    S = np.cross(E.T, B.T).T
    u2 = (E**2).sum(0) + (B**2).sum(0)
    return 2 * S / np.maximum(u2, 1e-300)


def advance(pos, F0, F1):
    """Heun step over one tick with the fields at the start (F0) and end (F1) of the tick"""
    E, B = sample(F0, pos)
    v0 = vel(E, B)
    p1 = pos + v0
    E, B = sample(F1, p1)
    v1 = vel(E, B)
    return pos + 0.5 * (v0 + v1), 0.5 * (v0 + v1)

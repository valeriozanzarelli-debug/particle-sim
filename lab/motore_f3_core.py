"""motore_f3_core.py - phase 3 of the "motore ottimizzato" study (docs/MOTORE_OTTIMIZZATO.md, section 15).

(1) Bell family for the local (lapped) level change: the kinetic-energy error of phase 2 (G3c) comes from the spectral
    tails of the bell at the block edges. A smoother bell is still orthogonal (conservative) and still local (same
    support): it is the 'correction at the block edge'. Family:
      'sin<k>'  : beta = sine iterated k times (phase 2 used k = 2), contact order 2^k at the overlap edges;
      'cinf'    : C-infinity bell, beta(s) = smooth step built from exp(-1/x) (Coifman-Meyer), super-algebraic tails.
    Blocks L, half-overlap e (fine); coarse L/2, e/2.
(2)/(3) helpers for the octant atom: more states (3s, 3d0 in the even sector, 3p in the odd one), dipole matrix
    elements between the sectors, exact Schroedinger kick probabilities."""
import os, sys
import numpy as np
import scipy.fft as sfft
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_f2_core as f2

OUT = f2.OUT
A0 = f2.A0


# =============================== (1) bells ===============================
def _smoothstep_cinf(s):
    """C-infinity step on [-1, 1] -> [-1, 1], odd, beta(+-1) = +-1 with all derivatives 0 there"""
    x = (s + 1) / 2                                    # 0..1
    f = lambda t: np.where(t > 0, np.exp(-1 / np.maximum(t, 1e-300)), 0.0)
    st = f(x) / (f(x) + f(1 - x))                      # 0..1, C-inf
    return 2 * st - 1


def bell(e, kind="sin2"):
    """rho(s), rho(-s) sampled at s = (t + 1/2)/e, t = 0..e-1; rho^2(s) + rho^2(-s) = 1"""
    s = (np.arange(e) + 0.5) / e
    if kind.startswith("sin"):
        b = s
        for _ in range(int(kind[3:])):
            b = np.sin(np.pi / 2 * b)
    elif kind == "cinf":
        b = _smoothstep_cinf(s)
    else:
        raise ValueError(kind)
    return np.sin(np.pi / 4 * (1 + b)), np.cos(np.pi / 4 * (1 + b))


def _fold(x, L, e, kind, sign, inverse):
    n = x.shape[0]
    rp, rm = bell(e, kind)
    sh = (e,) + (1,) * (x.ndim - 1)
    rp, rm = rp.reshape(sh), rm.reshape(sh)
    sg = -sign if inverse else sign
    for a in range(0, n, L):
        ip = (a + np.arange(e)) % n
        im = (a - 1 - np.arange(e)) % n
        u, v = x[ip].copy(), x[im].copy()
        x[ip] = rp * u + sg * rm * v
        x[im] = rp * v - sg * rm * u
    return x


def lot_fwd(x, L, e, axis=0, kind="sin2", sign=1):
    x = np.moveaxis(np.asarray(x), axis, 0).copy()
    x = _fold(x, L, e, kind, sign, False)
    n = x.shape[0]
    y = sfft.dct(x.reshape((n // L, L) + x.shape[1:]), type=4, axis=1, norm="ortho")
    return np.moveaxis(y.reshape(x.shape), 0, axis)


def lot_bwd(c, L, e, axis=0, kind="sin2", sign=1):
    c = np.moveaxis(np.asarray(c), axis, 0)
    n = c.shape[0]
    x = sfft.idct(c.reshape((n // L, L) + c.shape[1:]), type=4, axis=1, norm="ortho").reshape(c.shape).copy()
    x = _fold(x, L, e, kind, sign, True)
    return np.moveaxis(x, 0, axis)


def lot_restrict(F, L, e, axes, kind="sin2"):
    C = F
    for ax in axes:
        c = np.moveaxis(lot_fwd(C, L, e, ax, kind), ax, 0)
        n = c.shape[0]
        cb = c.reshape((n // L, L) + c.shape[1:])[:, : L // 2]
        C = lot_bwd(np.moveaxis(cb.reshape((n // 2,) + c.shape[1:]), 0, ax), L // 2, e // 2, ax, kind) / np.sqrt(2)
    return C


def lot_prolong(Fc, L, e, axes, kind="sin2"):
    C = Fc
    for ax in axes:
        c = np.moveaxis(lot_fwd(C * np.sqrt(2), L // 2, e // 2, ax, kind), ax, 0)
        nc = c.shape[0]
        cb = c.reshape((nc // (L // 2), L // 2) + c.shape[1:])
        full = np.zeros((cb.shape[0], L) + c.shape[1:], dtype=c.dtype)
        full[:, : L // 2] = cb
        C = lot_bwd(np.moveaxis(full.reshape((nc * 2,) + c.shape[1:]), 0, ax), L, e, ax, kind)
    return C


def electron_roundtrip(k0a, L, e, kind, halvings=4, nf=196608, sigma_a0=0.5):
    """G3c measurement (as motore_f2_run.part_lot): fine 1-D walk theta = a0/4096, positive packet, `halvings` local
    restrictions then prolongations; relative norm and fine-rule kinetic-energy changes"""
    import motore_core as mc
    th, _ = mc.load_levels_1d()
    xf = (np.arange(nf) - nf / 2 + 0.5) * th
    k = 2 * np.pi * np.fft.fftfreq(nf)
    Ekin = np.arccos(np.cos(th) * np.cos(k)) / th - 1.0
    psi = mc.positive_spinor_packet(xf, k0a / mc.A0, sigma_a0 * mc.A0, th)
    N0 = np.sum(np.abs(psi) ** 2)
    c = psi
    for _ in range(halvings):
        c = lot_restrict(c, L, e, (1,), kind)
    b = c
    for _ in range(halvings):
        b = lot_prolong(b, L, e, (1,), kind)
    fb, f0 = np.fft.fft(b, axis=-1), np.fft.fft(psi, axis=-1)
    Eb = np.sum(np.abs(fb) ** 2 * Ekin) / np.sum(np.abs(fb) ** 2)
    E0 = np.sum(np.abs(f0) ** 2 * Ekin) / np.sum(np.abs(f0) ** 2)
    return dict(k0_a0=k0a, L=L, e=e, bell=kind, dN_coarse=float(2 ** halvings * np.sum(np.abs(c) ** 2) / N0 - 1),
                dN_roundtrip=float(np.sum(np.abs(b) ** 2) / N0 - 1), dEkin_roundtrip=float(Eb / E0 - 1),
                max_err=float(np.abs(b - psi).max() / np.abs(psi).max()))


# =============================== (2) exact hydrogen kick (Schroedinger, a0 units) ===============================
def _radial(n, l, r):
    from scipy.special import genlaguerre, factorial
    rho = 2 * r / n
    norm = np.sqrt((2 / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho ** l * genlaguerre(n - l - 1, 2 * l + 1)(rho)


def exact_kick_probs(q, states=((1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)), nr=6000, rmax=120.0):
    """P(n l m=0) = |<n l 0| e^{iqz} |1s>|^2 = (2l+1) |int R_nl j_l(q r) R_10 r^2 dr|^2 (Rayleigh expansion)"""
    from scipy.special import spherical_jn
    r = np.linspace(0, rmax, nr + 1)[1:]
    dr = r[1] - r[0]
    R10 = _radial(1, 0, r)
    out = {}
    for n, l in states:
        I = np.sum(_radial(n, l, r) * spherical_jn(l, q * r) * R10 * r ** 2) * dr
        out["%d%s" % (n, "spd"[l])] = float((2 * l + 1) * I ** 2)
    return out


def exact_dipole_strengths(states=((2, 1), (3, 1)), nr=6000, rmax=120.0):
    """|<n p0| z |1s>|^2 in a0^2 (sum over all n and continuum = <1s|z^2|1s> = 1)"""
    r = np.linspace(0, rmax, nr + 1)[1:]
    dr = r[1] - r[0]
    R10 = _radial(1, 0, r)
    return {"%dp" % n: float((np.sum(_radial(n, 1, r) * R10 * r ** 3) * dr) ** 2 / 3) for n, l in states}

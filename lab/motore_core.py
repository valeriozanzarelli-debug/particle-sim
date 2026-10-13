"""motore_core.py - the PHYSICS side of the multilevel engine ("motore ottimizzato", docs/MOTORE_OTTIMIZZATO.md).

What lives here (all Python, numpy; the zig-zag rule is imported from zz_core / agg_zz_core, never changed):
  * LIGHT 2-D: the spin-1 zig-zag light walk evolved EXACTLY in Fourier space (two-tick Bloch matrix W(k) to the
    power m, by repeated squaring) and the CONTINUUM Maxwell reference exp(-i t k.S) (real light, the exam).
    The same initial data go to both. A common co-moving shift (exact Fourier translation) keeps the packet in the
    box. Measures: energy centroid, overlap phase, fidelity, helicity, width.
  * ERROR MODEL: phase and group-velocity errors of a level (cells per wavelength N, direction) from the Bloch
    dispersion of the rule itself -> predicted envelope drift and fringe drift after a time T -> error in PIXELS
    for a given footprint (wavelengths per pixel). Level choice = coarsest level whose projected error stays < 1 px
    (Nanite-like) and that is above the DYNAMICS FLOOR.
  * CONE: real-space ticks with a FROZEN outside (only the backward light cone of the viewed region is updated):
    exact by construction of the rule (one cell per axis per tick).
  * ELECTRON 1-D: the learned agglomerated levels of agg_zz_1d (lab/out/agg_zz_1d.npz) vs the exact fine walk.
Units: cell of the level = 1, tick = 1 (c = 1). For the electron: hbar = c = m_e = 1, fine cell = theta."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zz_core as zz

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SPIN1 = [np.asarray(S) for S in zz.SPIN1]


# =============================== light, 2-D, Fourier-exact ===============================
def kgrid2(n):
    k1 = 2 * np.pi * np.fft.fftfreq(n)
    KX, KY = np.meshgrid(k1, k1, indexing="ij")
    KA = np.zeros((n * n, 3))
    KA[:, 0], KA[:, 1] = KX.ravel(), KY.ravel()
    return KA


def lattice_W(n):
    """two-tick light walk matrix W(k) on an n x n periodic grid, shape (n*n, 3, 3)"""
    return zz.bloch_batch(kgrid2(n), 0.0, spin=1)


def matpow(W, m):
    """W^m for a batch of matrices, by repeated squaring"""
    R = np.broadcast_to(np.eye(W.shape[-1], dtype=complex), W.shape).copy()
    P = W.copy()
    while m > 0:
        if m & 1:
            R = P @ R
        P = P @ P
        m >>= 1
    return R


def continuum_U(n, t, scale=1.0):
    """continuum Maxwell propagator exp(-i t k.S) on the same grid (k in 1/cell; scale = cell size if needed)"""
    KA = kgrid2(n) / scale
    k = np.linalg.norm(KA, axis=1)
    kh = KA / np.where(k > 0, k, 1.0)[:, None]
    M = np.einsum("nj,jab->nab", kh, np.array(SPIN1))
    M2 = M @ M
    th = (k * t)[:, None, None]
    I = np.eye(3)[None]
    return I - 1j * np.sin(th) * M + (np.cos(th) - 1) * M2


def fft2F(F):
    return np.fft.fft2(F, axes=(0, 1)).reshape(-1, 3)


def ifft2F(Fk, n):
    return np.fft.ifft2(Fk.reshape(n, n, 3), axes=(0, 1))


def apply_k(U, F):
    n = F.shape[0]
    return ifft2F(np.einsum("nab,nb->na", U, fft2F(F)), n)


def shift_k(F, d):
    """exact Fourier translation by the real vector d (cells): F(x) -> F(x - d)"""
    n = F.shape[0]
    KA = kgrid2(n)
    ph = np.exp(-1j * (KA[:, 0] * d[0] + KA[:, 1] * d[1]))
    return ifft2F(fft2F(F) * ph[:, None], n)


def light_packet(n, lam, ang_deg, sigma_par, sigma_perp, center=None):
    """helicity +1 continuum packet sampled on the grid (NOT projected on the lattice band: same data for both)"""
    if center is None:
        center = (n / 2, n / 2)
    a = np.radians(ang_deg)
    kh = np.array([np.cos(a), np.sin(a), 0.0])
    ph = np.array([-np.sin(a), np.cos(a)])
    k0 = 2 * np.pi / lam
    x = np.arange(n) + 0.5
    X, Y = np.meshgrid(x - center[0], x - center[1], indexing="ij")
    u = X * kh[0] + Y * kh[1]
    v = X * ph[0] + Y * ph[1]
    env = np.exp(-u ** 2 / (2 * sigma_par ** 2) - v ** 2 / (2 * sigma_perp ** 2))
    e = zz_pol_plus(kh)
    return (env * np.exp(1j * k0 * u))[..., None] * e


def zz_pol_plus(kh):
    M = sum(kh[j] * SPIN1[j] for j in range(3))
    w, v = np.linalg.eigh(M)
    return v[:, np.argmax(w)]


def energy(F):
    return float(np.sum(np.abs(F) ** 2))


def centroid(F):
    n = F.shape[0]
    rho = np.sum(np.abs(F) ** 2, axis=2)
    x = np.arange(n) + 0.5
    E = rho.sum()
    return np.array([np.sum(rho.sum(1) * x), np.sum(rho.sum(0) * x)]) / E


def width_perp(F, ang_deg):
    n = F.shape[0]
    rho = np.sum(np.abs(F) ** 2, axis=2)
    c = centroid(F)
    a = np.radians(ang_deg)
    x = np.arange(n) + 0.5
    X, Y = np.meshgrid(x - c[0], x - c[1], indexing="ij")
    v = -X * np.sin(a) + Y * np.cos(a)
    return float(np.sqrt(np.sum(rho * v ** 2) / rho.sum()))


def overlap(Fa, Fb):
    """<Fa|Fb> / (|Fa||Fb|): phase and fidelity"""
    o = np.vdot(Fa.ravel(), Fb.ravel()) / np.sqrt(energy(Fa) * energy(Fb))
    return complex(o)


def helicity_plus_fraction(F):
    """fraction of energy on the CONTINUUM helicity +1 band"""
    n = F.shape[0]
    KA = kgrid2(n)
    k = np.linalg.norm(KA, axis=1)
    kh = KA / np.where(k > 0, k, 1.0)[:, None]
    M = np.einsum("nj,jab->nab", kh, np.array(SPIN1))
    Fk = fft2F(F)
    # projector on eigenvalue +1 of M (spin-1): (M^2 + M)/2
    P = 0.5 * (M @ M + M)
    Fp = np.einsum("nab,nb->na", P, Fk)
    return float(np.sum(np.abs(Fp) ** 2) / np.sum(np.abs(Fk) ** 2))


# =============================== error model (Bloch) ===============================
def omega_plus(ka):
    w, _ = zz.quasi_energies(np.asarray(ka, float), 0.0, 1.0)
    return float(w.max())


def level_errors(N, ang_deg, h=1e-6):
    """for a carrier with N cells per wavelength along ang: (domega/omega, v_g error vector (2,), |v_g| - 1)
    per tick of the level (c = 1). Positive band (helicity +1)."""
    a = np.radians(ang_deg)
    k = 2 * np.pi / N
    kv = k * np.array([np.cos(a), np.sin(a), 0.0])
    w = omega_plus(kv)
    vg = np.array([(omega_plus(kv + h * e) - omega_plus(kv - h * e)) / (2 * h) for e in np.eye(3)])[:2]
    vc = np.array([np.cos(a), np.sin(a)])
    return (w - k) / k, vg - vc, float(np.linalg.norm(vg) - 1)


def predicted_drift(N, ang_deg, T_lambda):
    """after T wavelengths of flight (T_lambda * N ticks): envelope drift (in wavelengths, vector) and
    overlap phase (rad) of the lattice packet with respect to real light"""
    dw, dvg, _ = level_errors(N, ang_deg)
    t = T_lambda * N
    k = 2 * np.pi / N
    dx = dvg * t / N                       # cells -> wavelengths
    dphi = -dw * k * t                     # overlap phase = -(omega_lat - omega) t
    return dx, dphi


def screen_error_px(dx_lambda, dphi, footprint_lambda_per_px):
    """visible displacement in pixels: envelope drift always; fringe drift (dphi/2pi wavelengths) only if the
    fringes are resolved on screen (lambda >= 2 px, Nyquist), otherwise the pixel averages them away"""
    p = footprint_lambda_per_px
    env = float(np.linalg.norm(dx_lambda)) / p
    fringes_resolved = (1.0 / p) >= 2.0
    fr = abs(dphi) / (2 * np.pi) / p if fringes_resolved else 0.0
    return max(env, fr), env, fr, fringes_resolved


# =============================== causal cone (real space) ===============================
def light_tick_real(F, ntick):
    """one tick of the 2-D light rule in real space (zz_core Walk 'light', axes x, y; palindromic order)"""
    order = (0, 1) if ntick % 2 == 0 else (1, 0)
    for j in order:
        lam, P = zz.BASIS_ONE[j]
        c = np.einsum("ba,xyb->xya", P.conj(), F)
        for q in range(3):
            s = int(lam[q])
            if s != 0:
                c[..., q] = np.roll(c[..., q], s, axis=j)
        F = np.einsum("ab,xyb->xya", P, c)
    return F


def dirac_tick_real(psi, ntick, theta, phase=None):
    """one tick of the 2-D massive zig-zag (4 comps L0 L1 R0 R1 on the last axis), axes x, y"""
    order = (0, 1) if ntick % 2 == 0 else (1, 0)
    L, R = psi[..., 0:2], psi[..., 2:4]
    for j in order:
        lam, P = zz.BASIS_HALF[j]
        for blk, hsign in ((0, -1), (1, +1)):
            X = L if blk == 0 else R
            c = np.einsum("ba,xyb->xya", P.conj(), X)
            for q in range(2):
                c[..., q] = np.roll(c[..., q], hsign * int(lam[q]), axis=j)
            X = np.einsum("ab,xyb->xya", P, c)
            if blk == 0:
                L = X
            else:
                R = X
    ct, st = np.cos(theta), np.sin(theta)
    out = np.concatenate([ct * L - 1j * st * R, ct * R - 1j * st * L], axis=-1)
    if phase is not None:
        out = out * phase[..., None]
    return out


def chebyshev_dist_to_box(n, box):
    """L-infinity distance (cells) of every site to the rectangle box = (x0, x1, y0, y1) (inclusive-exclusive)"""
    x = np.arange(n)
    dx = np.maximum(0, np.maximum(box[0] - x, x - (box[1] - 1)))
    dy = np.maximum(0, np.maximum(box[2] - x, x - (box[3] - 1)))
    return np.maximum(dx[:, None], dy[None, :])


# =============================== electron 1-D (agglomerated levels) ===============================
def load_levels_1d():
    """the automatic 1-D hierarchy of agg_zz_1d (theta, {B: Rule})"""
    import json
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from agg_zz_core import Rule
    d = np.load(os.path.join(OUT, "agg_zz_1d.npz"))
    js = json.load(open(os.path.join(OUT, "agg_zz_1d.json")))
    th = js["theta"]
    lv = {}
    for B in (4, 16, 64, 192, 384):
        lv[B] = Rule(d[f"auto_B{B}_K"], d[f"auto_B{B}_offs"], 1, 2, unitary=True)
    return th, lv


def fine1_power(k, th, m):
    """exact m-tick fine 1-D walk multiplier (n, 2, 2) at fine wave numbers k"""
    c, s = np.cos(th), np.sin(th)
    T = np.zeros(k.shape + (2, 2), complex)
    T[..., 0, 0] = np.exp(-1j * k)
    T[..., 1, 1] = np.exp(1j * k)
    C = np.array([[c, -1j * s], [-1j * s, c]])
    return matpow(C @ T, m)


def block1(psi, b):
    """agglomerator block (agg_zz_core.block, d = 1): mean x sqrt(b), psi shape (2, n)"""
    n = psi.shape[-1]
    return psi.reshape(2, n // b, b).mean(-1) * np.sqrt(b)


def omega_plus_batch(KA):
    """positive (helicity +1) band of the light walk, per tick, for many k (N, 3)"""
    W = zz.bloch_batch(KA, 0.0, spin=1)
    lam = np.linalg.eigvals(W)
    return np.max(-np.angle(lam) / 2.0, axis=1)


def predicted_drift_spectral(F, t, h=1e-5):
    """error model weighted by the packet's own spectrum (the region's FFT, which the engine has):
    envelope drift = t (<grad omega_lat> - <grad |k|>), overlap phase = -t <omega_lat - |k|>,
    weights = energy of the continuum helicity +1 part at each k. Returns (dx in cells (2,), dphi)."""
    n = F.shape[0]
    KA = kgrid2(n)
    k = np.linalg.norm(KA, axis=1)
    kh = KA / np.where(k > 0, k, 1.0)[:, None]
    M = np.einsum("nj,jab->nab", kh, np.array(SPIN1))
    P = 0.5 * (M @ M + M)
    Fk = fft2F(F)
    wgt = np.sum(np.abs(np.einsum("nab,nb->na", P, Fk)) ** 2, axis=1)
    sel = wgt > 1e-12 * wgt.max()
    Ks, ws = KA[sel], wgt[sel] / wgt[sel].sum()
    om = omega_plus_batch(Ks)
    g = np.zeros((len(Ks), 2))
    for j in range(2):
        e = np.zeros(3); e[j] = h
        g[:, j] = (omega_plus_batch(Ks + e) - omega_plus_batch(Ks - e)) / (2 * h)
    ks = np.linalg.norm(Ks, axis=1)
    gc = Ks[:, :2] / np.where(ks > 0, ks, 1.0)[:, None]
    dx = t * np.sum(ws[:, None] * (g - gc), axis=0)
    dphi = -t * np.sum(ws * (om - ks))
    return dx, float(dphi)


# ----- electron scene helpers -----
A0 = 1.0 / zz.ALPHA


def positive_spinor_packet(xf, k0, sigma, th, x0=0.0):
    """Gaussian packet on the FINE grid (positions xf, fine cells of size th), mean momentum k0 (1/lambda_C),
    projected on the positive-energy band of the fine walk (so no Zitterbewegung unless asked)."""
    n = len(xf)
    g = np.exp(-(xf - x0) ** 2 / (4 * sigma ** 2) + 1j * k0 * (xf - x0))
    psi = np.stack([g, g]) / np.sqrt(2)
    return project_band_fine(psi, th, +1)


def fine_eig(n, th):
    k = 2 * np.pi * np.fft.fftfreq(n)
    c, s = np.cos(th), np.sin(th)
    T = np.zeros((n, 2, 2), complex)
    T[:, 0, 0] = np.exp(-1j * k)
    T[:, 1, 1] = np.exp(1j * k)
    C = np.array([[c, -1j * s], [-1j * s, c]])
    W = C @ T
    lam, V = np.linalg.eig(W)
    om = -np.angle(lam)              # per tick
    return k, W, lam, V, om


def project_band_fine(psi, th, sign):
    n = psi.shape[-1]
    k, W, lam, V, om = fine_eig(n, th)
    # positive band = the one continuously connected to +omega (|omega| larger root has sign of energy): choose om*sign>0
    sel = (om * sign > 0).astype(float)
    Vi = np.linalg.inv(V)
    P = np.einsum("nab,nb,nbc->nac", V, sel, Vi)
    f = np.fft.fft(psi, axis=-1).T            # (n, 2)
    return np.fft.ifft(np.einsum("nab,nb->na", P, f).T, axis=-1)


def fine_evolve(psi, th, m):
    n = psi.shape[-1]
    k, W, lam, V, om = fine_eig(n, th)
    Vi = np.linalg.inv(V)
    Wm = np.einsum("nab,nb,nbc->nac", V, lam ** m, Vi)
    f = np.fft.fft(psi, axis=-1).T
    return np.fft.ifft(np.einsum("nab,nb->na", Wm, f).T, axis=-1)


def x_mean_width(psi, x):
    rho = np.sum(np.abs(psi) ** 2, axis=0)
    N = rho.sum()
    m = np.sum(rho * x) / N
    return float(m), float(np.sqrt(np.sum(rho * (x - m) ** 2) / N)), float(N)

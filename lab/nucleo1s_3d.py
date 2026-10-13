"""nucleo1s_3d.py - the FULL 3-D zig-zag rule (cubic lattice, 4-component balls, lab/zz_core.py Walk('dirac') /
zz_fast.FastDirac3D / src/atomo/zz.rs) on the Coulomb problem, solved cleanly.  Study NUCLEO_1S, phase 2.

Two exact facts of the rule are used (nothing is changed in the rule):
 1. FOUR COLOURS. In one tick a ball moves one cell along x, one along y and one along z (split step), so it always
    goes from a site with index parities (p, p, p) to (1-p, 1-p, 1-p) ... more precisely parity vector q -> q + (1,1,1).
    The 8 parity classes form 4 pairs {q, q+(1,1,1)} that NEVER talk to each other: the cubic lattice is 4 independent
    body-centred-cubic (bcc) lattices.  The local rotation (coin + Coulomb phase) is per site, so this holds with any
    potential.  Here we keep one colour: q = (0,0,0) ("E", all indices even) + (1,1,1) ("O", all odd).
    The nucleus (centre of the box, between 8 sites) has, in each colour, its two nearest sites on ONE body diagonal:
    each colour sees an axial (D3d) nucleus, the 4 colours are copies rotated onto the 4 body diagonals.
 2. FUSED HOP. The three axis shifts of one tick are one gather from the 8 body-diagonal neighbours:
        out(s) = sum_{e in {+-1}^3} P_z^{ez} P_y^{ey} P_x^{ex} psi_L(s + e)   (left-handed: sigma_j=e moves by -e)
        out(s) = sum_e        P_z^{ez} P_y^{ey} P_x^{ex} psi_R(s - e)   (right-handed: moves by +e)
    (order xyz on even ticks; on odd ticks P_x P_y P_z), then the per-site rotation
        L' = ph (cos th L - i sin th R), R' = ph (cos th R - i sin th L), ph = e^{-i V a}, th = m a.
    Checked against zz_fast.FastDirac3D in the self-test.

Levels: the overlap C(t) = <psi0|psi(t)> recorded while the rule runs ("listening"), then harmonic inversion
(filter diagonalisation, Mandelshtam-Taylor 1997) instead of a Fourier window: it separates close levels with a
finite listening time.  Units hbar = c = m_e = 1.  Inputs: m_e, alpha, Z.  Dirac numbers = CHECK ONLY.
"""
import os, sys, time
import numpy as np
import numba as nb
import scipy.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import BASIS_HALF, ALPHA, RadialWalk, dirac_formula  # noqa: E402


def _proj():
    P = np.zeros((3, 2, 2, 2), complex)        # P[j, e_index] (e_index 0: +1, 1: -1)
    for j in range(3):
        lam, V = BASIS_HALF[j]
        for q in range(2):
            P[j, q] = np.outer(V[:, q], V[:, q].conj())
    M = np.zeros((2, 8, 2, 2), complex)          # [order, combo]
    E = np.zeros((8, 3), np.int64)
    c = 0
    for qx in range(2):
        for qy in range(2):
            for qz in range(2):
                E[c] = (1 - 2 * qx, 1 - 2 * qy, 1 - 2 * qz)
                M[0, c] = P[2, qz] @ P[1, qy] @ P[0, qx]
                M[1, c] = P[0, qx] @ P[1, qy] @ P[2, qz]
                c += 1
    return M, E


MATS, EVEC = _proj()


@nb.njit(cache=True)
def _hop(src, dst, M, E, to_odd, ph, ct, st):
    """dst = coin(hop(src)); src, dst: (m, m, m, 4) = (L0, L1, R0, R1). to_odd: dst is the O class.
    Full index of E-site I is 2I, of O-site I is 2I+1. A displacement d in {-1,+1} from O-site full index 2I+1 to source
    full 2I+1-d: for dst O, src E index I + (1-d)//2; for dst E (full 2I), source 2I-d -> O index I - (1+d)//2."""
    m = src.shape[0]
    for i in range(m):
        for j in range(m):
            for k in range(m):
                l0 = 0j; l1 = 0j; r0 = 0j; r1 = 0j
                for c in range(8):
                    ex, ey, ez = E[c, 0], E[c, 1], E[c, 2]
                    # left-handed: component e moves by -e  -> came from s + e  (d = -e)
                    # right-handed: moves by +e -> came from s - e (d = +e)
                    if to_odd:
                        il = (i + (1 + ex) // 2) % m; jl = (j + (1 + ey) // 2) % m; kl = (k + (1 + ez) // 2) % m
                        ir = (i + (1 - ex) // 2) % m; jr = (j + (1 - ey) // 2) % m; kr = (k + (1 - ez) // 2) % m
                    else:
                        il = (i - (1 - ex) // 2) % m; jl = (j - (1 - ey) // 2) % m; kl = (k - (1 - ez) // 2) % m
                        ir = (i - (1 + ex) // 2) % m; jr = (j - (1 + ey) // 2) % m; kr = (k - (1 + ez) // 2) % m
                    a0 = src[il, jl, kl, 0]; a1 = src[il, jl, kl, 1]
                    b0 = src[ir, jr, kr, 2]; b1 = src[ir, jr, kr, 3]
                    l0 += M[c, 0, 0] * a0 + M[c, 0, 1] * a1
                    l1 += M[c, 1, 0] * a0 + M[c, 1, 1] * a1
                    r0 += M[c, 0, 0] * b0 + M[c, 0, 1] * b1
                    r1 += M[c, 1, 0] * b0 + M[c, 1, 1] * b1
                p = ph[i, j, k]
                dst[i, j, k, 0] = p * (ct * l0 - 1j * st * r0)
                dst[i, j, k, 1] = p * (ct * l1 - 1j * st * r1)
                dst[i, j, k, 2] = p * (ct * r0 - 1j * st * l0)
                dst[i, j, k, 3] = p * (ct * r1 - 1j * st * l1)


@nb.njit(cache=True)
def _dot(a, b):
    s = 0j
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            for k in range(a.shape[2]):
                for q in range(4):
                    s += np.conj(a[i, j, k, q]) * b[i, j, k, q]
    return s


class Colour3D:
    """one colour (E + O classes) of the cubic 3-D rule, n sites per axis (n even), spacing a, periodic box.
    V: function r -> potential (default point Coulomb -Za/r at the sites, as in the Laboratorio)."""

    def __init__(self, n, a, Za, mass=1.0, Vfun=None, q=(0, 0, 0)):
        assert n % 2 == 0
        self.n, self.a, self.Za, self.m, self.q = n, a, Za, n // 2, tuple(q)
        x = (np.arange(n) - n / 2 + 0.5) * a
        I = np.arange(n // 2)
        # colour q: class E' = full indices 2I + q_j, class O' = 2I + q_j + 1 (periodic)
        self.xE3 = [x[(2 * I + q[j]) % n] for j in range(3)]
        self.xO3 = [x[(2 * I + q[j] + 1) % n] for j in range(3)]
        self.xE, self.xO = x[0::2], x[1::2]
        self.ct, self.st = np.cos(mass * a), np.sin(mass * a)
        Vf = Vfun if Vfun is not None else (lambda r: -Za / r)
        self.rE = self._r(self.xE)
        self.rO = self._r(self.xO)
        self.VE, self.VO = Vf(self.rE), Vf(self.rO)
        self.phE = np.exp(-1j * self.VE * a)
        self.phO = np.exp(-1j * self.VO * a)
        self.ticks = 0

    def _r(self, x):
        X, Y, Z = np.meshgrid(*(self.xE3 if x is self.xE else self.xO3), indexing="ij")
        return np.sqrt(X**2 + Y**2 + Z**2)

    def coords(self, cls):
        return np.meshgrid(*(self.xE3 if cls == 0 else self.xO3), indexing="ij")

    def tick(self, E, O, E2, O2):
        """one tick: (E, O) -> (E2, O2). returns the new pair (buffers swapped by the caller)."""
        o = self.ticks % 2
        _hop(E, O2, MATS[o], EVEC, True, self.phO, self.ct, self.st)
        _hop(O, E2, MATS[o], EVEC, False, self.phE, self.ct, self.st)
        self.ticks += 1
        return E2, O2

    def from_full(self, psi):
        """restrict a full-lattice field (4, n, n, n) to this colour, as (E, O) arrays (m, m, m, 4)"""
        E = np.ascontiguousarray(np.moveaxis(psi[:, 0::2, 0::2, 0::2], 0, -1))
        O = np.ascontiguousarray(np.moveaxis(psi[:, 1::2, 1::2, 1::2], 0, -1))
        return E, O

    def field(self, fun):
        """fill (E, O) with fun(X, Y, Z) -> (..., 4) spinor field (chiral L0 L1 R0 R1)"""
        out = []
        for c in (0, 1):
            X, Y, Z = self.coords(c)
            out.append(np.ascontiguousarray(fun(X, Y, Z)).astype(complex))
        nr = np.sqrt(np.sum(np.abs(out[0])**2) + np.sum(np.abs(out[1])**2))
        return out[0] / nr, out[1] / nr

    def listen(self, E, O, nt, every=1, filt=None, log=None):
        """run nt ticks from (E, O); record C(t) = <psi0|psi(t)> each tick. filt: list of omegas for which the
        windowed sum sum_t w(t) e^{i omega t a} psi(t) is accumulated (returns the filtered states)."""
        E0, O0 = E.copy(), O.copy()
        Ea, Oa = E.copy(), O.copy()
        Eb, Ob = np.empty_like(E), np.empty_like(O)
        C = np.zeros(nt + 1, complex)
        acc = None
        if filt is not None:
            acc = [(np.zeros_like(E), np.zeros_like(O)) for _ in filt]
        t0 = time.time()
        for t in range(nt + 1):
            C[t] = _dot(E0, Ea) + _dot(O0, Oa)
            if acc is not None and t % 2 == 0:      # even ticks only: eigenstates of the two-tick rule
                w = np.sin(np.pi * (t + 0.5) / (nt + 1))**2
                for (ae, ao), om in zip(acc, filt):
                    f = w * np.exp(1j * om * t * self.a)
                    ae += f * Ea
                    ao += f * Oa
            if t == nt:
                break
            Eb, Ob = self.tick(Ea, Oa, Eb, Ob)
            Ea, Eb, Oa, Ob = Eb, Ea, Ob, Oa
            if log is not None and t in (9, nt // 4, nt // 2):
                log(f"    tick {t + 1}/{nt}: {(time.time() - t0) / (t + 1) * 1e3:.1f} ms/tick")
        return C, (Ea, Oa), acc


# ------------------------------------------------------------------ harmonic inversion
def fdm(c, a, wlo, whi, K=None, Mh=None, pad=1.5):
    """filter diagonalisation of c_n = sum_k d_k lambda_k^n (lambda = e^{-i omega a}), window [wlo, whi] in omega.
    Returns list of (omega, |d|, d, residual-ish) sorted by omega, for roots inside the window."""
    N = len(c) - 2
    Mh = Mh if Mh is not None else N // 2
    T = Mh * a
    if K is None:
        K = max(4, int(np.ceil(pad * (whi - wlo) * T / (2 * np.pi))) + 2)
    phis = np.linspace(wlo, whi, K) * a                 # basis frequencies
    z = np.exp(-1j * phis)                              # |z| = 1, z_j = lambda at omega_j
    n = np.arange(Mh + 1)
    F = z[None, :] ** (-n[:, None])                     # (Mh+1, K) : z_j^{-n}
    H0 = sla.hankel(c[0:Mh + 1], c[Mh:2 * Mh + 1])
    H1 = sla.hankel(c[1:Mh + 2], c[Mh + 1:2 * Mh + 2])
    U0 = F.T @ H0 @ F
    U1 = F.T @ H1 @ F
    # regularise: drop the near-null space of U0
    s, Q = np.linalg.eigh(0.5 * (U0 + U0.conj().T)) if False else (None, None)
    w, B = sla.eig(U1, U0)
    out = []
    Cv = F.T @ c[:Mh + 1]
    for k in range(len(w)):
        lam = w[k]
        if not np.isfinite(lam) or abs(lam) < 1e-6:
            continue
        b = B[:, k]
        nrm = b @ U0 @ b
        if abs(nrm) < 1e-300:
            continue
        d = (b @ Cv)**2 / nrm
        om = -np.angle(lam) / a
        if wlo <= om <= whi:
            out.append((om, abs(d), d, abs(abs(lam) - 1)))
    out.sort()
    return out


def pencil(c, a, wlo, whi, rank_tol=1e-7, L_frac=0.4):
    """harmonic inversion by band-pass + decimation + matrix pencil (SVD-truncated, Hua-Sarkar). c_n = sum d_k lam_k^n,
    lam = e^{-i omega a}.  Returns [(omega, |d|, d, |1-|lam||)] for roots in [wlo, whi]."""
    c = np.asarray(c, complex)
    N = len(c)
    n = np.arange(N)
    wc, bw = 0.5 * (wlo + whi), 0.5 * (whi - wlo)
    x = c * np.exp(1j * wc * a * n)                       # band centre -> 0
    # FIR low-pass (Hann-windowed sinc), cutoff 1.6 bw, applied without phase shift (valid part only)
    D = max(1, int(np.floor(np.pi / (a * 3.2 * bw))))     # decimation keeps |omega - wc| < pi/(D a)
    half = 6 * D
    k = np.arange(-half, half + 1)
    h = np.sinc(2 * 1.6 * bw * a / (2 * np.pi) * k) * np.hanning(2 * half + 1)
    h /= h.sum()
    y = np.convolve(x, h, mode="valid")[::D]            # sample j corresponds to n = half + j D
    M = len(y)
    L = max(2, int(L_frac * M))
    Y = np.array([y[i:i + L + 1] for i in range(M - L)])
    U, s, Vh = np.linalg.svd(Y, full_matrices=False)
    r = int(np.sum(s > rank_tol * s[0]))
    V = Vh[:r].T                                       # rows of Y span (z^0..z^L)
    V1, V2 = V[:-1], V[1:]
    z = np.linalg.eigvals(np.linalg.pinv(V1) @ V2)       # z = lam^D * e^{i wc a D}
    om = wc - np.angle(z) / (a * D)
    # amplitudes by least squares on y
    jj = np.arange(M)
    Z = z[None, :] ** jj[:, None]
    dd = np.linalg.lstsq(Z, y, rcond=None)[0]
    # undo the shift of the start (n0 = half): y_j ~ sum d lam^{half}... and the filter gain at that frequency
    out = []
    for q in range(len(z)):
        if not (wlo <= om[q] <= whi):
            continue
        gain = np.sum(h * np.exp(-1j * (om[q] - wc) * a * (-k)))   # filter response
        lam_half = np.exp(-1j * (om[q] - wc) * a * half)
        d = dd[q] / (gain * lam_half) if abs(gain) > 1e-12 else np.nan
        out.append((float(om[q]), float(abs(d)), d, float(abs(1 - abs(z[q])))))
    out.sort()
    return out


def fdm_stable(c, a, wlo, whi, dmin=1e-4, tol=None):
    """roots that are stable between two signal lengths (full, 0.8 x) and have weight > dmin"""
    r1 = pencil(c, a, wlo, whi)
    r2 = pencil(c[:int(0.8 * len(c))], a, wlo, whi)
    tol = tol if tol is not None else 2e-4
    good = []
    for (om, ad, d, g) in r1:
        if ad < dmin:
            continue
        if r2:
            k = np.argmin([abs(om - q[0]) for q in r2])
            dd = abs(om - r2[k][0])
        else:
            dd = np.inf
        good.append(dict(omega=om, weight=ad, decay=g, stab=dd, ok=bool(dd < tol)))
    return good


# ------------------------------------------------------------------ start states (as in the Laboratorio)
def radial_state(kappa, nlev, a, Za, rmax=40.0, guess=None):
    """level (kappa, n) of the radial form of the same rule with the same a (as src/atomo/zz.rs livello()):
    returns omega, and the sites r_i and u-functions (g, f) (r G, r F)"""
    w = RadialWalk(kappa, a, rmax, Za)
    g0 = guess if guess is not None else 1.0 - 0.5 * Za * Za / nlev**2 * 1.2
    om, u, res = w.eigen(g0, iters=10)
    g, f = w.gf(u)
    return om, w.r, g.real, f.real


def spinor_field(kappa, r_sites, g, f, a, mj=0.5, axis=None):
    """function (X, Y, Z) -> chiral spinor field of the radial level (j = 1/2: kappa = -1 s, +1 p1/2), m = +1/2 along
    z, linear interpolation in r as in zz.rs spinore().  Dirac rep (phi, chi) -> chiral L = (phi - chi)/sqrt2, R = (phi+chi)/sqrt2."""
    y00 = 0.5 / np.sqrt(np.pi)

    def fun(X, Y, Z):
        r = np.maximum(np.sqrt(X**2 + Y**2 + Z**2), 1e-9)
        x = r / a - 0.5
        i = np.clip(np.floor(np.maximum(x, 0)).astype(int), 0, len(g) - 2)
        t = np.clip(x - i, 0, 1)
        G = (g[i] * (1 - t) + g[i + 1] * t) / r
        F = (f[i] * (1 - t) + f[i + 1] * t) / r
        cz, cx, cy = Z / r, X / r, Y / r
        om_s = np.stack([y00 * np.ones_like(r), 0 * r], -1).astype(complex)
        om_p = np.stack([-y00 * cz, -y00 * (cx + 1j * cy)], -1)
        up, lo = (om_s, om_p) if kappa < 0 else (om_p, om_s)
        phi = up * G[..., None]
        chi = 1j * lo * F[..., None]
        s = 1 / np.sqrt(2)
        return np.concatenate([(phi - chi) * s, (phi + chi) * s], -1)
    return fun


def dirac_omega(n, kappa, Za):
    return dirac_formula(n, kappa, Za)


# ------------------------------------------------------------------ self-test
if __name__ == "__main__":
    from zz_fast import FastDirac3D
    rng = np.random.default_rng(1)
    n, a, Za = 12, 0.5, 60 * ALPHA
    col = Colour3D(n, a, Za)
    x = (np.arange(n) - n / 2 + 0.5) * a
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    V = -Za / np.sqrt(X**2 + Y**2 + Z**2)
    W = FastDirac3D(n, a, a, V)
    psi = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n))
    mask = ((np.arange(n)[:, None, None] % 2 == np.arange(n)[None, :, None] % 2) &
            (np.arange(n)[None, :, None] % 2 == np.arange(n)[None, None, :] % 2))
    psi *= mask[None]
    E, O = col.from_full(psi)
    E2, O2 = np.empty_like(E), np.empty_like(O)
    for t in range(5):
        psi = W.tick(psi)
        E2, O2 = col.tick(E, O, E2, O2)
        E, E2, O, O2 = E2, E, O2, O
    Ef, Of = col.from_full(psi)
    print("colour walk vs full FastDirac3D after 5 ticks:", np.max(np.abs(Ef - E)), np.max(np.abs(Of - O)))
    print("leak into other colours:", np.max(np.abs(psi * (~mask)[None])))
    print("norm:", np.sum(np.abs(E)**2) + np.sum(np.abs(O)**2), np.sum(np.abs(psi)**2))
    # FDM self test
    om = np.array([0.9, 0.91, 0.95]); d = np.array([0.6, 0.3, 0.1])
    tt = np.arange(1201)
    c = (d[None, :] * np.exp(-1j * om[None, :] * tt[:, None] * a)).sum(1) + 1e-3 * np.exp(-1j * 0.3 * tt * a)
    print("fdm:", [(round(o, 8), round(w, 6)) for o, w, _, _ in fdm(c, a, 0.85, 1.0) if w > 1e-3])
    print("pencil:", [(round(o, 8), round(w, 6)) for o, w, _, _ in pencil(c, a, 0.85, 1.0) if w > 1e-3])

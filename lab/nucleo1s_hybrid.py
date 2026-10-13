"""nucleo1s_hybrid.py - true hydrogen (Z = 1) in radial form with TWO rules glued together:
  * near the nucleus (r < R_f): the zig-zag radial rule itself (Dirac, 2 walkers per site), on the zoom mesh;
  * far away: the AGGLOMERATED rule learned in AGGLOMERATO_ZIGZAG.md (lab/out/agg_zz_3d_A.pkl, level B = 256: cell
    a0/16 = 8.56 lambda_C, or B = 512: a0/8), reduced to the s wave: its learned generator G(k) (spin-blind part),
    averaged over directions, T(k) = G_iso(k)/tau, acts on u = r psi as a 1-D stencil with u(-r) = -u(r).
Glue (stationary states, energy E): the coarse cells j < J0 are "ghosts" filled with the zig-zag solution regular at
the nucleus (large component G = r g at r_j, same E); the coarse equations j >= J0 are solved with them as boundary
data; E is the energy at which the coarse equation of the last ghost cell (j = J0 - 1) is also satisfied.
Nothing is fitted to the answer: the inner rule is the zig-zag, the outer rule is the learned one, the coarse
potential is the shell average of -alpha/r (written by hand: debt 9 of CONTESTO_CLOUD, unchanged).
"""
import os, sys, pickle
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_core import Mesh, mesh_zoom, rot_chi, gamma_of, ALPHA, MC2_EV, dirac_binding  # noqa
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / ALPHA
TH_A = A0 / 4096                     # fine cell of chain A (AGGLOMERATO_ZIGZAG.md section 4)


def fib_dirs(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5**0.5) * i
    d = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], -1)
    return d


class CoarseRadial:
    """the learned level as an s-wave radial operator on cells h (coarse units: k in 1/cell)"""

    def __init__(self, B=256, ndir=120, nk=1024, nfine=32768, vmode="point"):
        self.vmode = vmode
        st = pickle.load(open(os.path.join(OUT, "agg_zz_3d_A.pkl"), "rb"))
        lv = [l for l in st["levels"] if l["B"] == B][0]
        self.B, self.h, self.tau = B, B * TH_A, B * TH_A
        offs = lv["offs"].astype(float)
        K = lv["K"]
        Gs = 0.5 * (K[0, 0] + K[1, 1])                 # spin-blind part of the generator stencil
        self.spin_part = float(np.max(np.abs(K[0, 1])) / np.max(np.abs(Gs)))
        dirs = fib_dirs(ndir)
        ks = np.linspace(0, np.pi, nk)
        Giso = np.zeros(nk)
        for c0 in range(0, nk, 128):
            kv = ks[c0:c0 + 128, None, None] * dirs[None, :, :]
            ph = np.exp(1j * np.tensordot(kv, offs.T, axes=(-1, 0)))
            Giso[c0:c0 + 128] = np.real(ph @ Gs).mean(axis=1)
        self.G0 = Giso[0]
        self.ks = ks
        self.Tk = (Giso - self.G0) / self.tau          # kinetic energy (m c^2) vs k (1/cell)
        kk = (self.h / A0) * dirs
        g1 = np.real(np.exp(1j * kk @ offs.T) @ Gs)
        self.aniso_1a0 = float((g1.max() - g1.min()) / max(abs(g1.mean() - self.G0), 1e-300))
        from scipy.interpolate import CubicSpline
        self._spl = CubicSpline(ks, self.Tk)
        self._kf = np.linspace(0, np.pi, nfine)
        self._Tf = self._spl(self._kf)
        kphys = ks[1:6] / self.h
        self.m_kin = float(np.mean(kphys**2 / (2 * self.Tk[1:6])))
        self._tcache = {}

    def t(self, nmax):
        """the full 1-D kernel t_n = (1/pi) int_0^pi T(k) cos(nk) dk, n = 0..nmax (no truncation of the tail)"""
        if nmax not in self._tcache:
            kf, Tf = self._kf, self._Tf
            w = np.full(len(kf), kf[1] - kf[0]); w[0] *= 0.5; w[-1] *= 0.5
            n = np.arange(nmax + 1)
            out = np.empty(nmax + 1)
            for c0 in range(0, nmax + 1, 64):
                out[c0:c0 + 64] = (np.cos(np.outer(n[c0:c0 + 64], kf)) @ (w * Tf)) / np.pi
            self._tcache[nmax] = out
        return self._tcache[nmax]

    def V(self, N, Za=ALPHA):
        j = np.arange(N)
        if self.vmode == "shell":
            inv = 1.5 * ((j + 1)**2 - j**2) / ((j + 1)**3 - j**3) / self.h      # shell average of 1/r
        else:
            inv = 1.0 / ((j + 0.5) * self.h)                                   # value at the point r_j
        return -Za * inv

    def H(self, N, Za=ALPHA):
        """dense coarse Hamiltonian with the odd reflection at r = 0 (u_{-j-1} = -u_j), u = 0 beyond N"""
        t = self.t(2 * N + 1)
        j = np.arange(N)
        Hm = t[np.abs(j[:, None] - j[None, :])] - t[j[:, None] + j[None, :] + 1]
        return Hm + np.diag(self.V(N, Za))


def zz_regular(mesh, kappa, Za, w, refl0=-1.0):
    """the zig-zag stationary solution regular at the nucleus at frequency w, run outward on the mesh;
    returns r_i and G_i = r g (large component) at the palindromic point, normalised to G real."""
    M = rot_chi(mesh.h, kappa, mesh.r, Za)
    Mh = rot_chi(0.5 * mesh.h, kappa, mesh.r, Za)
    N = mesh.N
    p = np.zeros(N, complex); q = np.zeros(N + 1, complex)
    q[0] = 1.0
    inp = refl0 * np.exp(1j * w * mesh.l_in) * q[0]
    lk = np.concatenate([mesh.l, [mesh.l_out]])
    for i in range(N):
        inm = (q[i] - M[i, 1, 0] * inp) / M[i, 1, 1]
        p[i] = M[i, 0, 0] * inp + M[i, 0, 1] * inm
        q[i + 1] = inm * np.exp(-1j * w * lk[i])
        inp = np.exp(1j * w * lk[i]) * p[i]
    q = q[:N]
    u0 = np.conj(Mh[:, 0, 0]) * p + np.conj(Mh[:, 1, 0]) * q
    u1 = np.conj(Mh[:, 0, 1]) * p + np.conj(Mh[:, 1, 1]) * q
    g = (u0 + u1) / np.sqrt(2)
    f = 1j * (u0 - u1) / np.sqrt(2)
    k = np.argmax(np.abs(g))
    ph = np.exp(-1j * np.angle(g[k]))
    return mesh.r, g * ph, f * ph


class Hybrid:
    def __init__(self, coarse, rmax, J0=2, a_fine=TH_A, zoom_levels=12, npl=16, Za=ALPHA, kappa=-1):
        self.c, self.J0, self.Za, self.kappa = coarse, J0, Za, kappa
        h = coarse.h
        self.N = int(np.ceil(rmax / h))
        self.r = (np.arange(self.N) + 0.5) * h
        self.Hm = coarse.H(self.N, Za)
        R_f = (J0 + 0.5) * h
        self.mesh = mesh_zoom(a_fine, R_f, zoom_levels, npl)
        self.cells_fine = self.mesh.N

    def ghosts(self, E):
        r, g, f = zz_regular(self.mesh, self.kappa, self.Za, 1.0 + E)
        self.imag_ratio = float(np.max(np.abs(g.imag)) / np.max(np.abs(g.real)))
        return np.interp(self.r[:self.J0], r, g.real)

    def residual(self, E, return_u=False):
        J0, N = self.J0, self.N
        ug = self.ghosts(E)
        A = self.Hm - E * np.eye(N)
        ui = np.linalg.solve(A[J0:, J0:], -A[J0:, :J0] @ ug)
        u = np.concatenate([ug, ui])
        res = A[J0 - 1] @ u
        if return_u:
            return res, u
        return res / np.abs(ug).max()

    def solve(self, E_lo, E_hi):
        E = brentq(self.residual, E_lo, E_hi, xtol=1e-18, rtol=1e-14, maxiter=200)
        res, u = self.residual(E, return_u=True)
        return E, u


def coarse_only(coarse, rmax, n_state=1, Za=ALPHA):
    N = int(np.ceil(rmax / coarse.h))
    vals = np.linalg.eigvalsh(coarse.H(N, Za))
    return vals[n_state - 1]

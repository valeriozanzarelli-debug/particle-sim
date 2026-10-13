"""aggfull_atom.py - the coarse hydrogen atom with the FULLY LEARNED rule (no hand-written pieces in the rule).

Coarse generator (one coarse step = e^{-i H tau}, tau = h/c):
    H = T_e(k) [+ T_p(-k)] + V_c
  * T_e = G_e(k)/tau: learned free generator of chain A (lab/out/agg_zz_3d_A.pkl, spin-blind part; its spin part is
    1e-10), symmetrised k <-> -k (the fine walk is exactly parity symmetric, checked to 2e-16; the learned stencils have an
    odd part made of round-off, reported).
  * T_p = G_p(-k)/tau: the PROTON's learned generator (lab/out/aggfull_proton.pkl), for the recoil runs: the pair state in
    the sector of zero total momentum (exact for translation-invariant rules), k = relative momentum.
  * V_c: the learned coupling (lab/out/aggfull_learn.pkl): <k'|V_c|k> = Vt(k'-k) Gamma(k',k), Gamma = sum C_jl e^{-ik'j} e^{ikl}.
    The eigenproblem is solved with Gamma -> gamma0 = Gamma(0,0) (spin-blind, k-independent leading term); the rest
    dGamma = Gamma - gamma0 (Darwin-like, spin-orbit, ...; size alpha^2 relative) is added in FIRST ORDER (error
    O(alpha^4) relative, below everything else here). Vt = field of the nucleus (input -Z alpha / r) restricted to the
    coarse Fourier modes without aliasing (aggfull_core.GalerkinField).
Eigenstates by LOBPCG (real arithmetic for parity-symmetric problems: the 1s, 2s and 2p_z are real).
Usage: python lab/aggfull_atom.py JOB [JOB ...]   (jobs defined in JOBS; results in lab/out/aggfull_atom.json; vectors .npy)
"""
import os, sys, time, json, pickle
import numpy as np
import scipy.sparse.linalg as spl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_core as ac
from aggfull_core import fftn, ifftn

DEV = os.environ.get("AGGFULL_DEV") == "1"
TAG = "dev_aggfull_atom" if DEV else "aggfull_atom"
JS = os.path.join(ac.OUT, TAG + ".json")
ALPHA = ac.zz.ALPHA


# ------------------------------------------------------------------ learned pieces on a grid
def parity_even(G):
    """(G(k) + G(-k))/2 on an FFT grid"""
    Gm = np.roll(np.flip(G, axis=(0, 1, 2)), 1, axis=(0, 1, 2))
    return 0.5 * (G + Gm), float(np.max(np.abs(G - Gm)) / 2)


def kinetic(rec, N, h):
    """learned spin-blind generator per unit time on the N^3 grid (cell h), parity-even part; returns T, odd, spin"""
    r = ac.genrule(rec)
    Gk = r.gen_grid(N)
    spin = float(np.max(np.abs(Gk[..., 0, 1])) + np.max(np.abs(Gk[..., 0, 0] - Gk[..., 1, 1])))
    g = np.real(0.5 * (Gk[..., 0, 0] + Gk[..., 1, 1]))
    del Gk
    g0 = float(np.real(0.5 * np.trace(r.gen_at(np.zeros((1, 3)))[0])))
    g, odd = parity_even(g - g0)
    return g / h, odd / h, spin


def kgrid2(N, h):
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij", sparse=True)
    return KX ** 2 + KY ** 2 + KZ ** 2


def kinetic_mass(rec, h, m_ref=1.0):
    """kinetic mass from the parity-even learned generator at k = 1/64..1/16 per cell (x, diagonal, (1,2,0)); mean"""
    r = ac.genrule(rec)
    g0 = float(np.real(0.5 * np.trace(r.gen_at(np.zeros((1, 3)))[0])))
    out = []
    for u in (np.array([1.0, 0, 0]), np.ones(3) / np.sqrt(3), np.array([1.0, 2.0, 0]) / np.sqrt(5)):
        kc = np.array([1 / 64, 1 / 32, 1 / 16])[:, None] * u[None]
        gp = np.real(0.5 * np.trace(r.gen_at(kc), axis1=1, axis2=2))
        gm = np.real(0.5 * np.trace(r.gen_at(-kc), axis1=1, axis2=2))
        e = (0.5 * (gp + gm) - g0) / h
        k = np.linalg.norm(kc, axis=1) / h
        out.append(k ** 2 / (2 * e) / m_ref)
    out = np.array(out)
    return float(out[:, 0].mean()), out.tolist()


# ------------------------------------------------------------------ the atom
class Atom:
    def __init__(self, T, field, gamma0=1.0, real=True, sym=None):
        """sym: None, 'even' (even in x, y, z about the nucleus) or 'oddz' (odd in z, even in x, y): the parity sector
        is kept exactly (the operator conserves it; projecting stops round-off from sliding into the ground state)"""
        self.T, self.f, self.g0, self.real, self.sym = T, field, gamma0, real, sym
        self.N = T.shape[0]
        self.pre = 1.0 / (T - T.min() + 0.5 * ALPHA ** 2)
        self.napply = 0

    def project(self, v):
        if self.sym is None:
            return v
        fl = lambda f, ax: np.roll(np.flip(f, axis=ax), 1, axis=ax)
        sz = -1.0 if self.sym == "oddz" else 1.0
        v = 0.5 * (v + fl(v, 0))
        v = 0.5 * (v + fl(v, 1))
        return 0.5 * (v + sz * fl(v, 2))

    def H(self, v):
        v = self.project(v)
        F = fftn(v)
        out = ifftn(self.T * F + self.g0 * self.f.apply_k(F))
        return self.project(out.real if self.real else out)

    def _apply(self, X, fn):
        X = np.asarray(X)
        one = X.ndim == 1
        X = X.reshape(self.N ** 3, -1)
        out = np.empty_like(X)
        for c in range(X.shape[1]):
            out[:, c] = fn(X[:, c].reshape((self.N,) * 3)).ravel()
        return out[:, 0] if one else out

    def A(self, X):
        self.napply += 1
        return self._apply(X, self.H)

    def P(self, X):
        def pre(v):
            o = ifftn(self.pre * fftn(self.project(v)))
            return self.project(o.real if self.real else o)
        return self._apply(X, pre)

    def energy(self, v):
        f = v.reshape((self.N,) * 3)
        Hf = self.H(f)
        nrm = np.vdot(f, f).real
        E = np.vdot(f, Hf).real / nrm
        return float(E), float(np.linalg.norm(Hf - E * f) / np.sqrt(nrm))


def guess(x, kind, a0=ac.A0):
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    g = {"1s": np.exp(-r / a0), "2s": (1 - r / (2 * a0)) * np.exp(-r / (2 * a0)), "2p": Z * np.exp(-r / (2 * a0))}[kind]
    return (g / np.linalg.norm(g)).ravel()


def solve(atom, x0, Y=None, tol=1e-10, maxiter=40):
    N3 = atom.N ** 3
    dt = float if atom.real else complex
    Aop = spl.LinearOperator((N3, N3), matvec=atom.A, matmat=atom.A, dtype=dt)
    Pop = spl.LinearOperator((N3, N3), matvec=atom.P, matmat=atom.P, dtype=dt)
    x0 = x0.astype(dt)[:, None]
    if Y is not None:
        x0 = x0 - Y @ (Y.conj().T @ x0)
    w, v, hist = spl.lobpcg(Aop, x0, M=Pop, Y=Y, tol=tol, maxiter=maxiter, largest=False, retResidualNormsHistory=True)
    return v[:, 0], (float(np.max(np.abs(hist[-1]))) if len(hist) else float("nan")), len(hist)


# ------------------------------------------------------------------ first-order pieces
def shifted(phi, l):
    """(T_l phi)(X) = phi(X + l)"""
    return np.roll(phi, shift=(-l[0], -l[1], -l[2]), axis=(0, 1, 2))


def coupling_blocks(field, phis, offs):
    """W[p][l] = V_G (T_l phi_p) for all sites l of the stencil (Galerkin field), and the overlaps
    S[p', j, p, l] = <T_j phi_p' | V_G | T_l phi_p>"""
    W = [[field.apply(shifted(ph, o)) for o in offs] for ph in phis]
    n, P = len(offs), len(phis)
    S = np.zeros((P, n, P, n), complex)
    for p2 in range(P):
        for j, oj in enumerate(offs):
            a = shifted(phis[p2], oj)
            for p in range(P):
                for l in range(n):
                    S[p2, j, p, l] = np.vdot(a, W[p][l])
    return S


def first_order_matrix(C, offs, S, gamma0):
    """matrix of dV = sum_jl (C_jl - gamma0 delta_j0 delta_l0) T_j^dag V T_l on the basis phi_p x spin_a:
    M[(p',a),(p,b)] = sum_jl dC[a,b,j,l] S[p',j,p,l]"""
    dC = C.copy()
    i0 = [i for i, o in enumerate(offs) if not np.any(o)][0]
    dC[:, :, i0, i0] -= gamma0 * np.eye(2)
    P = S.shape[0]
    M = np.einsum("abjl,qjpl->qapb", dC, S).reshape(2 * P, 2 * P)
    return 0.5 * (M + M.conj().T)


def kinetic_remainder(T, h, phi, m_kin):
    """<phi| T(k) - k^2/(2 m_kin) |phi>: the part of the learned kinetic energy beyond the parabola (contains -k^4/8m^3)"""
    k2 = kgrid2(T.shape[0], h)
    F = fftn(phi)
    w = np.abs(F) ** 2
    return float(np.sum(w * (T - k2 / (2 * m_kin))) / np.sum(w))

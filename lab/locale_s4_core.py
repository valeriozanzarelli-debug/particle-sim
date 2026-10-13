"""locale_s4_core.py - phase 2 of REGOLA LOCALE: tools for S4 (the 'passo del gambero' composition) in 1-D,
in the radial (hydrogen) form, and for causality runs. Uses lab/locale_core.py and lab/zz_core.py (frozen, read only).

S4 = palindromic composition of Strang blocks B(m) = C(m/2) . Shift^m . C(m/2) with m in SEQ (sum m^3 = 0).
  1-D Dirac: Shift moves R by +1, L by -1 (one tick of the frozen dirac1 walk); C(m/2) = exp(-i (m/2) theta X)
             (and the potential/field phase for time m/2).  B(1)B(1)... = the frozen 1-D walk up to a conjugation.
  radial:    Shift = zz_core.RadialWalk.shift (chi+ out, chi- in, reflection at r = 0); C = the local rotation
             exp(-i dt (m sigma_z + kappa/r sigma_x + V)) for time m dt / 2.
SEQ_FROZEN = (1,) gives exactly (a conjugate of) the frozen rule.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from zz_core import RadialWalk, U_CHI, ALPHA

PI = np.pi
SEQ_S4 = (1, 1, 1, 1, -2, 1, 1, 1, 1)
SEQ_FROZEN = (1,)


def seq_c(c):
    """palindrome (1 x c^3/2, -c, 1 x c^3/2): sum m^3 = 0, lattice cone (c^2+1)/(c^2-1); c even"""
    h = c ** 3 // 2
    return tuple([1] * h + [-c] + [1] * h)


def cone(seq):
    return sum(abs(m) for m in seq) / sum(seq)


# ------------------------------------------------------------------ 1-D Dirac, k-space
def dirac1_bloch(k, theta, seq, A=None, t0=0.0):
    """period matrix (2x2, basis (R, L)) at momentum k; A(t) = vector potential (minimal coupling k -> k - A(t)),
    evaluated at the middle of every shift (time runs backwards inside a negative block). returns W, T."""
    W = np.eye(2, dtype=complex)
    t = t0
    X = np.array([[0, 1], [1, 0]], complex)

    def C(a):
        return np.cos(a) * np.eye(2) - 1j * np.sin(a) * X
    for m in seq:
        W = C(m * theta / 2) @ W
        ke = k - (A(t + m / 2) if A is not None else 0.0)
        W = np.diag([np.exp(-1j * m * ke), np.exp(1j * m * ke)]) @ W
        W = C(m * theta / 2) @ W
        t += m
    return W, float(sum(seq))


# ------------------------------------------------------------------ radial (hydrogen)
class RadialS4:
    """the zig-zag radial walk for one kappa, composed with the integer sequence seq (period time sum(seq) dt)."""

    def __init__(self, kappa, a, rmax, seq, Zalpha=ALPHA, mass=1.0):
        self.rw = RadialWalk(kappa, a, rmax, Zalpha=Zalpha, mass=mass)
        self.a, self.seq, self.N = a, seq, self.rw.N
        r = self.rw.r
        V = -Zalpha / r
        hx, hz = kappa / r, mass * np.ones_like(r)
        nrm = np.sqrt(hx ** 2 + hz ** 2)
        self._loc = (V, hx, hz, nrm)
        N = self.N
        # shift matrix S (interleaved index 2 i + s) from RadialWalk.matrix with M = identity
        i = np.arange(N)
        rows, cols, vals = [], [], []
        src_p = np.where(i >= 1, 2 * (i - 1), 1); fac_p = np.where(i >= 1, 1.0, self.rw.refl0)
        src_q = np.where(i <= N - 2, 2 * (i + 1) + 1, 2 * (N - 1)); fac_q = np.where(i <= N - 2, 1.0, self.rw.refl1)
        rows = np.concatenate([2 * i, 2 * i + 1]); cols = np.concatenate([src_p, src_q]); vals = np.concatenate([fac_p, fac_q])
        self.S = sp.csc_matrix((vals.astype(complex), (rows, cols)), shape=(2 * N, 2 * N))
        self.Sinv = self.S.conj().T.tocsc()

    def local(self, tau):
        """block-diagonal local rotation exp(-i tau (m sz + kappa/r sx + V)) in the chi basis, tau = time"""
        V, hx, hz, nrm = self._loc
        ph = tau * nrm
        c, s = np.cos(ph), np.sin(ph)
        e = np.exp(-1j * V * tau)
        M = np.zeros((self.N, 2, 2), complex)
        M[:, 0, 0] = c - 1j * s * hz / nrm; M[:, 1, 1] = c + 1j * s * hz / nrm
        M[:, 0, 1] = -1j * s * hx / nrm; M[:, 1, 0] = -1j * s * hx / nrm
        M *= e[:, None, None]
        M = np.einsum("ab,nbc,cd->nad", U_CHI.conj().T, M, U_CHI)
        return self._bd(M)

    def _bd(self, M):
        N = self.N
        i = np.arange(N)
        rows = np.concatenate([2 * i, 2 * i, 2 * i + 1, 2 * i + 1])
        cols = np.concatenate([2 * i, 2 * i + 1, 2 * i, 2 * i + 1])
        vals = np.concatenate([M[:, 0, 0], M[:, 0, 1], M[:, 1, 0], M[:, 1, 1]])
        return sp.csc_matrix((vals, (rows, cols)), shape=(2 * N, 2 * N))

    def matrix(self):
        W = sp.identity(2 * self.N, dtype=complex, format="csc")
        for m in self.seq:
            Ch = self.local(m * self.a / 2)
            Sm = self.S if m > 0 else self.Sinv
            B = Ch
            for _ in range(abs(m)):
                B = Sm @ B
            B = Ch @ B
            W = (B @ W).tocsc()
        return W

    def eigen(self, omega_guess, iters=10, W=None):
        W = self.matrix() if W is None else W
        T = sum(self.seq) * self.a
        lam0 = np.exp(-1j * omega_guess * T)
        lu = spl.splu((W - lam0 * sp.identity(2 * self.N, format="csc")).tocsc())
        x = np.ones(2 * self.N, complex); x /= np.linalg.norm(x)
        for _ in range(iters):
            y = lu.solve(x); x = y / np.linalg.norm(y)
        lam = np.vdot(x, W @ x)
        resid = float(np.linalg.norm(W @ x - lam * x))
        om = -np.angle(lam) / T
        # unfold: the bound states sit near omega_guess
        per = 2 * PI / T
        om = omega_guess + ((om - omega_guess) + per / 2) % per - per / 2
        return float(om), x, resid

"""aggmol_ecg.py - H2 with TWO electrons (phase 2 of study aggmol_, after Valerio lifted the molecule gate for H2):
the EXACT few-body level of the many-electron rule of the main thread (docs/MOLTI_ELETTRONI.md: explicitly correlated
Gaussians with the rule's operator), extended to TWO CENTRES (Gaussians with centres on the molecular axis).

The operator (the main thread's rule as it is, its debts named there):
  H = sum_i p_i^2 / (2 m_kin) - sum_i sum_nuclei 1/|r_i - X_n| + 1/|r_1 - r_2|     (Hartree units of m_e)
  m_kin = m0 of the v3 level (multi_e_v3.json on master: R3v3 m0 = 1.0000000010922): the learned band, quadratic part
  (DEBT-4 there: the alpha^2 rest of the band is not included; ~1e-5 Eh, declared); e-e as instantaneous Coulomb (DEBT-1),
  antisymmetry put in (DEBT-2: singlet spin -> symmetric space), nuclei clamped (Born-Oppenheimer) for the curve; the
  nuclear motion is then the 1-D vibration with the proton mass (an allowed input).
Basis: phi(r) = exp(-sum_ij A_ij (r_i - s_i).(r_j - s_j)), s_i on the z axis (keeps the Sigma symmetry exactly),
projected on (P12 symmetric) x (inversion gerade). Elements analytic (Gaussian moments, Coulomb = erf / |mu|).
Nonlinear parameters by the stochastic variational method (random trials, keep the best; refinement sweeps).
Usage: python lab/aggmol_ecg.py gate | curve | ..."""
import os, sys, json, time
import numpy as np
import scipy.linalg as sla
from scipy.special import erf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
M_KIN = 1.0000000010921748            # multi_e_v3.json (origin/master), R3v3 m0
SQ2 = np.sqrt(2.0)


def coul(mu, s2):
    """E[1/|x|], x ~ N(mu, s2 I_3); mu (...,3), s2 (...)"""
    m = np.linalg.norm(mu, axis=-1)
    s = np.sqrt(s2)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(m > 1e-10 * s, erf(m / (SQ2 * s)) / m, np.sqrt(2 / np.pi) / s)
    return v


def elems(A1, s1, A2, s2, nuclei, m=M_KIN):
    """batched <phi1|1, T, V|phi2>; A (..., N, N), s (..., N, 3); returns S, T, V (...)"""
    N = A1.shape[-1]
    C = A1 + A2
    Ci = np.linalg.inv(C)
    det = np.linalg.det(C)
    v = A1 @ s1 + A2 @ s2
    R0 = Ci @ v
    c0 = np.einsum("...ix,...ij,...jx->...", s1, A1, s1) + np.einsum("...ix,...ij,...jx->...", s2, A2, s2)
    ex = np.einsum("...ix,...ij,...jx->...", v, Ci, v) - c0
    S = (np.pi ** N / det) ** 1.5 * np.exp(ex)
    Sig = 0.5 * Ci
    M = A1 @ A2
    d1, d2 = R0 - s1, R0 - s2
    Tk = 4 * (3 * np.einsum("...ij,...ji->...", M, Sig) + np.einsum("...ix,...ij,...jx->...", d1, M, d2)) / (2 * m)
    V = 0.0
    for i in range(N):
        for X in nuclei:
            V = V - coul(R0[..., i, :] - X, Sig[..., i, i])
    if N == 2:
        V = V + coul(R0[..., 0, :] - R0[..., 1, :], Sig[..., 0, 0] + Sig[..., 1, 1] - 2 * Sig[..., 0, 1])
    return S, Tk * S, V * S


def sym_ops(N):
    """(permutation matrix, inversion sign) for the projector: singlet (space symmetric) x gerade"""
    ops = []
    perms = [np.eye(N)] + ([np.array([[0., 1.], [1., 0.]])] if N == 2 else [])
    for P in perms:
        for inv in (1.0, -1.0):
            ops.append((P, inv))
    return ops


class Basis:
    def __init__(self, N, nuclei):
        self.N, self.nuc = N, [np.array(x, float) for x in nuclei]
        self.A = np.zeros((0, N, N)); self.s = np.zeros((0, N, 3))
        self.ops = sym_ops(N)
        self.H = np.zeros((0, 0)); self.S = np.zeros((0, 0))

    def row(self, A, s, A2=None, s2=None):
        """symmetrised elements of candidates (A, s: (n, N, N), (n, N, 3)) against functions (A2, s2) (K, ...)"""
        A2 = self.A if A2 is None else A2
        s2 = self.s if s2 is None else s2
        n, K = len(A), len(A2)
        Sr = np.zeros((n, K)); Hr = np.zeros((n, K))
        if K == 0:
            return Sr, Hr
        for P, inv in self.ops:
            Ap = np.einsum("ji,kjl,lm->kim", P, A2, P)
            sp = inv * np.einsum("ij,kjx->kix", P, s2)
            S, T, V = elems(A[:, None], s[:, None], Ap[None], sp[None], self.nuc)
            Sr += S; Hr += T + V
        return Sr, Hr

    def diag(self, A, s):
        S, H = self.row(A, s, A, s)
        return np.diag(S).copy() if S.ndim == 2 else S, np.diag(H).copy()

    def energy(self, H, S, nstate=0):
        try:
            L = np.linalg.cholesky(S)
        except np.linalg.LinAlgError:
            return np.inf
        Li = sla.solve_triangular(L, np.eye(len(S)), lower=True)
        Hm = Li @ H @ Li.T
        w = np.linalg.eigvalsh(0.5 * (Hm + Hm.T))
        # reject near-singular overlap (condition) as a numerical safety
        if np.min(np.linalg.eigvalsh(S)) < 1e-13 * np.max(np.diag(S)):
            return np.inf
        return w[nstate]

    def normed(self, A, s):
        """normalisation factors of the symmetrised candidates"""
        d = np.array([self.row(A[i:i + 1], s[i:i + 1], A[i:i + 1], s[i:i + 1])[0][0, 0] for i in range(len(A))])
        return 1 / np.sqrt(d)

    def add_best(self, gen, ntrial, rng):
        A, s = gen(rng, ntrial)
        nr = self.normed(A, s)
        Sr, Hr = self.row(A, s)
        if len(self.A):
            Sr *= self.nrm[None]; Hr *= self.nrm[None]
        Sd = np.array([self.row(A[i:i + 1], s[i:i + 1], A[i:i + 1], s[i:i + 1]) for i in range(ntrial)])
        K = len(self.A)
        best, bi = np.inf, -1
        for i in range(ntrial):
            if not np.isfinite(nr[i]):
                continue
            S = np.zeros((K + 1, K + 1)); H = np.zeros((K + 1, K + 1))
            S[:K, :K] = self.S; H[:K, :K] = self.H
            S[K, :K] = S[:K, K] = Sr[i] * nr[i]; H[K, :K] = H[:K, K] = Hr[i] * nr[i]
            S[K, K] = Sd[i, 0][0, 0] * nr[i] ** 2; H[K, K] = Sd[i, 1][0, 0] * nr[i] ** 2
            if K and np.max(np.abs(S[K, :K])) > 0.995:
                continue
            e = self.energy(H, S)
            if e < best:
                best, bi, Sb, Hb = e, i, S, H
        if bi < 0:
            return np.inf
        self.A = np.concatenate([self.A, A[bi:bi + 1]]); self.s = np.concatenate([self.s, s[bi:bi + 1]])
        self.nrm = np.concatenate([getattr(self, "nrm", np.zeros(0)), [nr[bi]]])
        self.S, self.H = Sb, Hb
        return best

    def refine(self, gen, ntrial, rng):
        """replace each function in turn by the best of ntrial candidates (if it lowers E)"""
        K = len(self.A)
        e0 = self.energy(self.H, self.S)
        for k in range(K):
            keep = np.arange(K) != k
            A, s = gen(rng, ntrial)
            h = ntrial // 2                       # half the trials: small random moves of the function itself
            L = np.linalg.cholesky(self.A[k])
            for i in range(h):
                Lp = L * np.exp(rng.normal(0, 0.15, L.shape)) * np.tril(np.ones_like(L))
                A[i] = Lp @ Lp.T
                s[i] = self.s[k] + np.concatenate([np.zeros((self.N, 2)), rng.normal(0, 0.1, (self.N, 1))], 1) / np.sqrt(np.diag(self.A[k]))[:, None] * 0.5
            nr = self.normed(A, s)
            Sr, Hr = self.row(A, s, self.A[keep], self.s[keep])
            Sr *= self.nrm[keep][None] * nr[:, None]; Hr *= self.nrm[keep][None] * nr[:, None]
            Ssub = self.S[np.ix_(keep, keep)]; Hsub = self.H[np.ix_(keep, keep)]
            for i in range(ntrial):
                if not np.isfinite(nr[i]) or np.max(np.abs(Sr[i])) > 0.995:
                    continue
                Sd, Hd = self.row(A[i:i + 1], s[i:i + 1], A[i:i + 1], s[i:i + 1])
                S = np.zeros((K, K)); H = np.zeros((K, K))
                S[:K - 1, :K - 1] = Ssub; H[:K - 1, :K - 1] = Hsub
                S[K - 1, :K - 1] = S[:K - 1, K - 1] = Sr[i]; H[K - 1, :K - 1] = H[:K - 1, K - 1] = Hr[i]
                S[K - 1, K - 1] = Sd[0, 0] * nr[i] ** 2; H[K - 1, K - 1] = Hd[0, 0] * nr[i] ** 2
                e = self.energy(H, S)
                if e < e0 - 1e-12:
                    e0 = e
                    order = np.concatenate([np.where(keep)[0], [k]])
                    self.A[k], self.s[k], self.nrm[k] = A[i], s[i], nr[i]
                    # rebuild matrices in the original order
                    inv = np.argsort(order)
                    self.S, self.H = S[np.ix_(inv, inv)], H[np.ix_(inv, inv)]
                    Ssub = self.S[np.ix_(keep, keep)]; Hsub = self.H[np.ix_(keep, keep)]
        return e0

    def rebuild(self):
        S, H = self.row(self.A, self.s)
        n = self.normed(self.A, self.s)
        self.nrm = n
        self.S = S * n[:, None] * n[None]; self.H = H * n[:, None] * n[None]
        return self.energy(self.H, self.S)


def make_gen(N, R, scale=1.0):
    """random trial functions: widths log-uniform, e-e correlation term, centres near the nuclei or between"""
    def gen(rng, n):
        a = np.exp(rng.uniform(np.log(0.03), np.log(30.0), (n, N))) * scale ** -2
        A = np.zeros((n, N, N))
        A[:, range(N), range(N)] = a
        if N == 2:
            b = np.exp(rng.uniform(np.log(1e-3), np.log(3.0), n)) * rng.choice([1, 1, 1, -0.3], n)
            b = np.maximum(b, -0.45 * np.min(a, axis=1))
            A[:, 0, 0] += b; A[:, 1, 1] += b; A[:, 0, 1] = A[:, 1, 0] = -b
        s = np.zeros((n, N, 3))
        cen = rng.choice([-0.5, 0.5, 0.0, -0.5, 0.5], (n, N)) * R + rng.normal(0, 0.25, (n, N)) * scale
        s[..., 2] = cen
        return A, s
    return gen


def build(N, R, K, ntrial=24, nref=2, seed=1, log=print, budget=500, state=None):
    """grow a basis of K functions for N electrons, protons at +-R/2 on z; returns energy and the basis"""
    t0 = time.time()
    rng = np.random.default_rng(seed)
    nuc = [(0, 0, R / 2), (0, 0, -R / 2)] if R > 0 else [(0, 0, 0)]
    B = Basis(N, nuc)
    gen = make_gen(N, max(R, 1.0))
    if state is not None:
        B.A, B.s = state["A"], state["s"]
        e = B.rebuild()
    e = np.inf
    while len(B.A) < K and time.time() - t0 < budget:
        e = B.add_best(gen, ntrial, rng)
    for _ in range(nref):
        if time.time() - t0 > budget:
            break
        e = B.refine(gen, ntrial, rng)
    return e, B, time.time() - t0

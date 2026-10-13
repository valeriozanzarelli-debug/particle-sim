"""mol_ecg.py - study mol_ (docs/MOLECOLE.md), tier "esatto": explicitly correlated Gaussians (ECG) for a FEW electrons
around SEVERAL nuclei on the z axis, any spin state given as a spin function. Extends the two-centre, two-electron machine
of aggmol_ (lab/aggmol_ecg.py, docs/AGGLOMERATO_DUE_ATOMI.md §12) to N electrons and nuclei of any charge (LiH: 4 e).

Operator (the main thread's many-electron rule as it is, same as aggmol phase 2; debts named there):
  H = sum_i p_i^2/(2 m_kin) - sum_i sum_n Z_n/|r_i - X_n| + sum_{i<j} 1/|r_i - r_j| + sum_{n<m} Z_n Z_m/|X_n - X_m|
  m_kin = aggmol_ecg.M_KIN (learned band, v3 m0); nuclei clamped (Born-Oppenheimer, debt).
Basis: phi(r) = exp(-sum_ij A_ij (r_i - s_i).(r_j - s_j)), centres s_i on the z axis (Sigma symmetry exact).
Antisymmetry: Psi = A[phi chi], chi a fixed spin function; matrix elements <phi_k|O| sum_P c_P P phi_l>,
c_P = sign(P) <chi|P chi>. Imported read-only: aggmol_ecg.coul, aggmol_ecg.M_KIN."""
import os, sys, itertools, time
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_ecg as AE

OUT = AE.OUT
M_KIN = 1.0000000045603115      # learned rule mass m0, multi_e v5 (R3v5), as tav_ecg.M0 and the helium export; aggmol used v3 (1.0000000010922): 3.5e-9 relative, ~3e-8 Eh at LiH
coul = AE.coul
UP, DN = np.array([1., 0.]), np.array([0., 1.])


def spin_singlet_pairs(npair, extra_up=0):
    """chi = prod singlet(2k, 2k+1) x alpha^extra_up"""
    th = (np.kron(UP, DN) - np.kron(DN, UP)) / np.sqrt(2)
    v = np.ones(1)
    for _ in range(npair):
        v = np.kron(v, th)
    for _ in range(extra_up):
        v = np.kron(v, UP)
    return v


def spin_terms(chi, N):
    chi = chi.reshape((2,) * N)
    out = []
    for p in itertools.permutations(range(N)):
        Pm = np.eye(N)[list(p)]
        c = np.linalg.det(Pm) * np.sum(chi * np.transpose(chi, p))
        if abs(c) > 1e-12:
            out.append((Pm, float(c)))
    return out


class System:
    def __init__(self, N, nuclei, chi, m=M_KIN):
        """nuclei: list of (Z, z position)"""
        self.N, self.m = N, m
        self.Z = np.array([z for z, _ in nuclei], float)
        self.X = [np.array((0., 0., x)) for _, x in nuclei]
        self.ops = spin_terms(chi, N)
        self.Enn = sum(self.Z[a] * self.Z[b] / abs(nuclei[a][1] - nuclei[b][1])
                       for a in range(len(nuclei)) for b in range(a + 1, len(nuclei)))

    def raw(self, A1, s1, A2, s2, dip=False):
        """batched <phi1|1, H|phi2> (no symmetriser); optionally <phi1| sum z_i |phi2>"""
        N = self.N
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
        T = 4 * (3 * np.einsum("...ij,...ji->...", M, Sig) + np.einsum("...ix,...ij,...jx->...", d1, M, d2)) / (2 * self.m)
        V = self.Enn
        for i in range(N):
            for Z, X in zip(self.Z, self.X):
                V = V - Z * coul(R0[..., i, :] - X, Sig[..., i, i])
            for j in range(i + 1, N):
                V = V + coul(R0[..., i, :] - R0[..., j, :], Sig[..., i, i] + Sig[..., j, j] - 2 * Sig[..., i, j])
        if dip:
            return S, (T + V) * S, R0[..., :, 2].sum(-1) * S
        return S, (T + V) * S

    def sym(self, A1, s1, A2, s2, dip=False):
        """<phi1| O P_chi |phi2>, broadcasting (A1 (...,N,N) vs A2)"""
        Ss = Hs = Ds = 0.0
        for Pm, c in self.ops:
            Ap = np.einsum("ji,...jl,lm->...im", Pm, A2, Pm)
            sp = np.einsum("ji,...jx->...ix", Pm, s2)
            r = self.raw(A1, s1, Ap, sp, dip)
            Ss = Ss + c * r[0]; Hs = Hs + c * r[1]
            if dip:
                Ds = Ds + c * r[2]
        return (Ss, Hs, Ds) if dip else (Ss, Hs)

    def fast_setup(self):
        import mol_fast
        self.F = mol_fast
        self.Pm = np.array([P for P, _ in self.ops]); self.cops = np.array([c for _, c in self.ops])
        self.Zz = self.Z.copy(); self.Xz = np.array([x[2] for x in self.X])

    def perm_basis(self, A, s):
        """(nops, K, N, N), (nops, K, N): P^T A P, P^T s (z only)"""
        Ap = np.einsum("oji,kjl,olm->okim", self.Pm, A, self.Pm)
        sp = np.einsum("oji,kj->oki", self.Pm, s[..., 2])
        return np.ascontiguousarray(Ap), np.ascontiguousarray(sp)

    def frows(self, A, s, Ap, sp):
        return self.F.rows(np.ascontiguousarray(A), np.ascontiguousarray(s[..., 2]), Ap, sp, self.cops, self.Zz,
                           self.Xz, self.m, self.Enn)

    def fdiag(self, A, s):
        return self.F.diag(np.ascontiguousarray(A), np.ascontiguousarray(s[..., 2]), self.Pm, self.cops, self.Zz,
                           self.Xz, self.m, self.Enn)

    def mats(self, A, s, dip=False):
        r = self.sym(A[:, None], s[:, None], A[None], s[None], dip)
        return tuple(0.5 * (x + x.T) for x in r)


def lowest(H, S, k=0):
    try:
        w = sla.eigh(H, S, eigvals_only=True, subset_by_index=[0, k])
    except (np.linalg.LinAlgError, ValueError):
        return np.inf
    return w[k]


def bordered_lowest(Ss, Hs, Sr, Hr, Hd):
    """lowest eigenvalue after adding each candidate (rows Sr, Hr; normalised self overlap 1, self energy Hd) to the
    basis (Ss, Hs): one eigensolve of the basis, then a secular equation per candidate (vectorised bisection+Newton)"""
    nc = len(Hd)
    if len(Ss) == 0:
        return Hd.copy()
    w, U = sla.eigh(Hs, Ss)
    p = Sr @ U; q = Hr @ U
    n2 = 1 - np.sum(p * p, 1)
    out = np.full(nc, np.inf)
    ok = (n2 > 1e-9) & (np.max(np.abs(Sr), 1) < 0.99)
    if not ok.any():
        return out
    p, q, n2, hd = p[ok], q[ok], n2[ok], Hd[ok]
    n = np.sqrt(n2)
    g = (q - w[None] * p) / n[:, None]
    d = (hd - 2 * np.sum(p * q, 1) + np.sum(w[None] * p * p, 1)) / n2
    g2 = g * g
    hi = np.full(len(d), w[0]) - 1e-14 * max(1.0, abs(w[0]))
    lo = np.minimum(d, w[0]) - np.sqrt(np.sum(g2, 1)) - 1e-6
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = mid - d + np.sum(g2 / (w[None] - mid[:, None]), 1)
        pos = f > 0
        hi = np.where(pos, mid, hi); lo = np.where(pos, lo, mid)
    out[ok] = 0.5 * (lo + hi)
    return out


class SVM:
    """stochastic variational method: grow by best-of-n trials, refine by replacing functions"""

    def __init__(self, sysm, gen, rng):
        self.sy, self.gen, self.rng = sysm, gen, rng
        N = sysm.N
        self.A = np.zeros((0, N, N)); self.s = np.zeros((0, N, 3))
        self.S = np.zeros((0, 0)); self.H = np.zeros((0, 0)); self.E = np.inf
        self.nrm = np.zeros(0)
        sysm.fast_setup()
        self.Ap, self.sp = sysm.perm_basis(self.A, self.s)

    def set(self, A, s):
        self.A, self.s = A.copy(), s.copy()
        self.Ap, self.sp = self.sy.perm_basis(self.A, self.s)
        S, H = self.sy.frows(self.A, self.s, self.Ap, self.sp)
        S = 0.5 * (S + S.T); H = 0.5 * (H + H.T)
        n = 1 / np.sqrt(np.diag(S))
        self.nrm = n
        self.S = S * n[:, None] * n[None]; self.H = H * n[:, None] * n[None]
        self.E = lowest(self.H, self.S)
        return self.E

    def _cands(self, A, s, keep):
        Sd, Hd = self.sy.fdiag(A, s)
        ok = np.isfinite(Sd) & (Sd > 0)
        nr = np.where(ok, 1 / np.sqrt(np.where(ok, Sd, 1)), 0)
        Sr, Hr = self.sy.frows(A, s, np.ascontiguousarray(self.Ap[:, keep]), np.ascontiguousarray(self.sp[:, keep]))
        Sr = Sr * nr[:, None] * self.nrm[keep][None]; Hr = Hr * nr[:, None] * self.nrm[keep][None]
        return ok, nr, Sr, Hr, Hd * nr ** 2

    def _try(self, Ss, Hs, sr, hr, hd):
        K = len(Ss)
        S = np.empty((K + 1, K + 1)); H = np.empty((K + 1, K + 1))
        S[:K, :K] = Ss; H[:K, :K] = Hs
        S[K, :K] = S[:K, K] = sr; H[K, :K] = H[:K, K] = hr
        S[K, K] = 1.0; H[K, K] = hd
        if K and np.max(np.abs(sr)) > 0.99:
            return np.inf, None, None
        return lowest(H, S), S, H

    def grow(self, ntrial):
        A, s = self.gen(self.rng, ntrial)
        keep = np.ones(len(self.A), bool)
        ok, nr, Sr, Hr, Hd = self._cands(A, s, keep)
        el = bordered_lowest(self.S, self.H, Sr, Hr, np.where(ok, Hd, 0.0))
        el[~ok] = np.inf
        i = int(np.argmin(el))
        if not np.isfinite(el[i]):
            return self.E
        e, S, H = self._try(self.S, self.H, Sr[i], Hr[i], Hd[i])
        if not np.isfinite(e):
            return self.E
        self.A = np.concatenate([self.A, A[i:i + 1]]); self.s = np.concatenate([self.s, s[i:i + 1]])
        self.nrm = np.concatenate([self.nrm, [nr[i]]])
        a, b = self.sy.perm_basis(A[i:i + 1], s[i:i + 1])
        self.Ap = np.ascontiguousarray(np.concatenate([self.Ap, a], 1)); self.sp = np.ascontiguousarray(np.concatenate([self.sp, b], 1))
        self.S, self.H, self.E = S, H, e
        return e

    def refine_one(self, k, ntrial):
        K = len(self.A)
        keep = np.arange(K) != k
        A, s = self.gen(self.rng, ntrial)
        h = ntrial // 2
        Lk = np.linalg.cholesky(self.A[k])
        for i in range(h):                          # small moves of the function itself
            Lp = Lk * np.exp(self.rng.normal(0, 0.12, Lk.shape)) * np.tril(np.ones_like(Lk))
            A[i] = Lp @ Lp.T
            s[i] = self.s[k].copy()
            s[i][:, 2] += self.rng.normal(0, 0.15, self.sy.N) / np.sqrt(np.diag(self.A[k]))
        ok, nr, Sr, Hr, Hd = self._cands(A, s, keep)
        Ss = self.S[np.ix_(keep, keep)]; Hs = self.H[np.ix_(keep, keep)]
        el = bordered_lowest(Ss, Hs, Sr, Hr, np.where(ok, Hd, 0.0))
        el[~ok] = np.inf
        i = int(np.argmin(el))
        if not el[i] < self.E - 1e-12:
            return self.E
        e, S, H = self._try(Ss, Hs, Sr[i], Hr[i], Hd[i])
        if not e < self.E - 1e-12:
            return self.E
        order = np.concatenate([np.where(keep)[0], [k]]); inv = np.argsort(order)
        self.A[k], self.s[k], self.nrm[k] = A[i], s[i], nr[i]
        a, b = self.sy.perm_basis(A[i:i + 1], s[i:i + 1])
        self.Ap[:, k] = a[:, 0]; self.sp[:, k] = b[:, 0]
        self.S, self.H, self.E = S[np.ix_(inv, inv)], H[np.ix_(inv, inv)], e
        return e

    def state(self):
        w, U = sla.eigh(self.H, self.S)
        c = U[:, 0] / np.sqrt(U[:, 0] @ self.S @ U[:, 0])
        return w[0], c * self.nrm              # coefficients for the unnormalised exp(-...)


def make_gen(N, groups):
    """groups: list of (prob, electron-kinds) where each electron has (centre z mean, centre spread, log10 width range).
    kinds: dict electron index -> (zc, dz, lo, hi). Correlation couplings between all pairs, log-uniform."""
    probs = np.array([g[0] for g in groups], float); probs /= probs.sum()

    def gen(rng, n):
        A = np.zeros((n, N, N)); s = np.zeros((n, N, 3))
        for k in range(n):
            kinds = groups[rng.choice(len(groups), p=probs)][1]
            a = np.zeros(N)
            for i in range(N):
                zc, dz, lo, hi = kinds[i]
                a[i] = 10 ** rng.uniform(lo, hi)
                s[k, i, 2] = zc + rng.normal(0, dz)
            Ak = np.diag(a)
            for i in range(N):
                for j in range(i + 1, N):
                    if rng.random() < 0.6:
                        b = 10 ** rng.uniform(-3, np.log10(0.5 * min(a[i], a[j]) + 1e-3) + 0.3)
                        if rng.random() < 0.15:
                            b = -0.2 * min(a[i], a[j]) * rng.random()
                        Ak[i, i] += b; Ak[j, j] += b; Ak[i, j] -= b; Ak[j, i] -= b
            if np.linalg.eigvalsh(Ak).min() <= 1e-6:
                Ak = np.diag(a)
            A[k] = Ak
        return A, s
    return gen

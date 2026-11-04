"""acqua_chol.py - integrali fattorizzati costruiti da noi (decomposizione di Cholesky a pivot della matrice
(mu nu|la si), soglia TAU_CH, nessuna base ausiliaria esterna), RHF fattorizzato, orbitali naturali MP2 congelati
(FNO, soglia TAU_NO sulle occupazioni virtuali, correzione Delta-MP2) e CCSD(T) a guscio chiuso sopra (acqua_qc).
Serve per il dimero (H2O)2 (studio ACQUA): memoria ~ N_ch * M^2 invece di M^4.
calc(...) lavora a pezzi con checkpoint (npz NON tracciati in lab/out/acqua_amp/)."""
import os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numba import njit
import acqua_int as I, acqua_qc as Q

TAU_CH = 1e-6


@njit(cache=True)
def eri_cols(sL, sa, sC, sN, off, CART, SPH, C, D, out):
    """out[mu, nu, k] = (mu nu|c d) for the components k of shell pair (C, D)"""
    ns = len(sL)
    Eab = np.zeros((15, 15, 9, 9, 9)); Ecd = np.zeros((15, 15, 9, 9, 9))
    R = np.zeros((17, 17, 17, 17)); F = np.zeros(17); Gk = np.zeros((9, 9, 9))
    blk = np.zeros((15, 15, 15, 15)); t1 = np.zeros((9, 15, 15, 15)); t2 = np.zeros((9, 9, 15, 15))
    t3 = np.zeros((9, 9, 9, 15)); o = np.zeros((9, 9, 9, 9))
    I._pairE(sL[C], sL[D], sa[C], sa[D], sC[C], sC[D], CART, Ecd)
    mc, md = 2 * sL[C] + 1, 2 * sL[D] + 1
    for A in range(ns):
        for B in range(A + 1):
            I._pairE(sL[A], sL[B], sa[A], sa[B], sC[A], sC[B], CART, Eab)
            I._quartet(A, B, C, D, sL, sa, sC, CART, SPH, Eab, Ecd, R, F, Gk, blk, t1, t2, t3, o)
            nrm = sN[A] * sN[B] * sN[C] * sN[D]
            for m in range(2 * sL[A] + 1):
                for n in range(2 * sL[B] + 1):
                    for p in range(mc):
                        for q in range(md):
                            x = o[m, n, p, q] * nrm
                            out[off[A] + m, off[B] + n, p * md + q] = x
                            out[off[B] + n, off[A] + m, p * md + q] = x


@njit(cache=True)
def eri_diag(sL, sa, sC, sN, off, CART, SPH, Dg):
    ns = len(sL)
    Eab = np.zeros((15, 15, 9, 9, 9))
    R = np.zeros((17, 17, 17, 17)); F = np.zeros(17); Gk = np.zeros((9, 9, 9))
    blk = np.zeros((15, 15, 15, 15)); t1 = np.zeros((9, 15, 15, 15)); t2 = np.zeros((9, 9, 15, 15))
    t3 = np.zeros((9, 9, 9, 15)); o = np.zeros((9, 9, 9, 9))
    for A in range(ns):
        for B in range(A + 1):
            I._pairE(sL[A], sL[B], sa[A], sa[B], sC[A], sC[B], CART, Eab)
            I._quartet(A, B, A, B, sL, sa, sC, CART, SPH, Eab, Eab, R, F, Gk, blk, t1, t2, t3, o)
            nrm = (sN[A] * sN[B]) ** 2
            for m in range(2 * sL[A] + 1):
                for n in range(2 * sL[B] + 1):
                    Dg[off[A] + m, off[B] + n] = o[m, n, m, n] * nrm
                    Dg[off[B] + n, off[A] + m] = o[m, n, m, n] * nrm


def cholesky(B, tau=TAU_CH, nmax=None):
    """returns L (Nch, M, M) with (mu nu|la si) ~ sum_P L[P,mu,nu] L[P,la,si]"""
    M = B.M; nmax = nmax or 12 * M
    shell_of = np.concatenate([[k] * (2 * B.sL[k] + 1) for k in range(len(B.sL))])
    Dg = np.zeros((M, M)); eri_diag(B.sL, B.sa, B.sC, B.sN, B.off, I.CART, I.SPH, Dg)
    D = Dg.ravel().copy()
    L = np.zeros((nmax, M * M)); n = 0
    while True:
        piv = int(np.argmax(D)); dmax = D[piv]
        if dmax < tau or n >= nmax:
            break
        mu, nu = divmod(piv, M); Cs, Ds = shell_of[mu], shell_of[nu]
        cols = np.zeros((M, M, (2 * B.sL[Cs] + 1) * (2 * B.sL[Ds] + 1)))
        eri_cols(B.sL, B.sa, B.sC, B.sN, B.off, I.CART, I.SPH, Cs, Ds, cols)
        cols = cols.reshape(M * M, -1)
        idx = [(B.off[Cs] + p) * M + B.off[Ds] + q for p in range(2 * B.sL[Cs] + 1) for q in range(2 * B.sL[Ds] + 1)]
        while True:          # take the components of this shell pair while they are large
            k = int(np.argmax(D[idx])); pv = idx[k]; dv = D[pv]
            if dv < max(tau, 1e-2 * dmax) or n >= nmax:
                break
            v = cols[:, k] - L[:n].T @ L[:n, pv]
            v /= math.sqrt(dv); L[n] = v; n += 1
            D -= v * v; D[pv] = 0.0
    assert n < nmax, "Cholesky: troppi vettori"
    out = L[:n].reshape(n, M, M).copy(); del L
    return out


def scf_df(B, ints, L, nocc, tol=1e-9, maxit=100):
    S, h = ints["S"], ints["T"] + ints["V"]; Enuc = B.enuc()
    w, U = np.linalg.eigh(S); keep = w > 1e-8; X = U[:, keep] / np.sqrt(w[keep])
    e, V = np.linalg.eigh(X.T @ h @ X); C = X @ V; hist = []; Eold = 0.0
    Lf = L.reshape(len(L), -1)
    for it in range(maxit):
        Co = C[:, :nocc]; Dm = Co @ Co.T
        J = (Lf.T @ (Lf @ Dm.ravel())).reshape(B.M, B.M)
        Y = np.einsum("Pmn,ni->Pmi", L, Co, optimize=True)
        K = np.einsum("Pmi,Pni->mn", Y, Y, optimize=True)
        Fk = h + 2 * J - K
        E = float(np.sum(Dm * (h + Fk))) + Enuc
        err = Fk @ Dm @ S - S @ Dm @ Fk
        hist.append((Fk, err)); hist = hist[-8:]
        if len(hist) > 1:
            n = len(hist); Bm = -np.ones((n + 1, n + 1)); Bm[-1, -1] = 0
            for i in range(n):
                for j in range(n):
                    Bm[i, j] = np.sum(hist[i][1] * hist[j][1])
            r = np.zeros(n + 1); r[-1] = -1
            Fk = sum(c * hh[0] for c, hh in zip(np.linalg.solve(Bm, r)[:n], hist))
        if np.abs(err).max() < tol and abs(E - Eold) < 1e-10:
            break
        Eold = E
        e, V = np.linalg.eigh(X.T @ Fk @ X); C = X @ V
    Fmo = C.T @ (h + 2 * J - K) @ C
    return dict(E=E, C=C, eps=np.diag(Fmo).copy(), conv=bool(np.abs(err).max() < tol), it=it)


class DFG:
    """lazy (pq|rs) from B[P, p, q] (MO); supports G[s1, s2, s3, s4] with slices or one leading int"""

    def __init__(self, Bm):
        self.B = Bm; self.shape = (Bm.shape[1],) * 4

    def __getitem__(self, k):
        a, b, c, d = k
        if isinstance(a, (int, np.integer)):
            return np.tensordot(self.B[:, a, b], self.B[:, c, d], axes=([0], [0]))
        X = self.B[:, a, b]; Y = self.B[:, c, d]
        P = X.shape[0]
        return (X.reshape(P, -1).T @ Y.reshape(P, -1)).reshape(X.shape[1:] + Y.shape[1:])


def mp2_fno(Bm, eps, nocc, tau_no):
    """canonical MP2 (full), MP2 virtual natural orbitals with occupation > tau_no, semicanonical; returns
    (B_new, eps_new, E_mp2_full, E_mp2_trunc, n_keep)"""
    nv = Bm.shape[1] - nocc; o = slice(0, nocc); v = slice(nocc, None)
    G = DFG(Bm); ovov = G[o, v, o, v]
    eo, ev = eps[o], eps[v]
    D = eo[:, None, None, None] - ev[None, :, None, None] + eo[None, None, :, None] - ev[None, None, None, :]
    t = ovov / D                                  # t[i,a,j,b]
    E2 = float(np.sum((2 * t - t.transpose(0, 3, 2, 1)) * ovov))
    T = t.transpose(0, 2, 1, 3)                   # [i,j,a,b]
    Dv = 2 * np.einsum("ijac,ijbc->ab", 2 * T - T.transpose(0, 1, 3, 2), T, optimize=True)
    occ, U = np.linalg.eigh(Dv); order = np.argsort(occ)[::-1]; occ, U = occ[order], U[:, order]
    k = int(np.sum(occ > tau_no)); U = U[:, :k]
    Fv = U.T @ np.diag(ev) @ U; e2, W = np.linalg.eigh(Fv); U = U @ W
    Cn = np.zeros((Bm.shape[1], nocc + k)); Cn[:nocc, :nocc] = np.eye(nocc); Cn[nocc:, nocc:] = U
    Bn = np.einsum("Ppq,pr,qs->Prs", Bm, Cn, Cn, optimize=True)
    en = np.concatenate([eo, e2])
    o2 = DFG(Bn)[o, v, o, v]
    D2 = eo[:, None, None, None] - e2[None, :, None, None] + eo[None, None, :, None] - e2[None, None, None, :]
    t2 = o2 / D2
    E2t = float(np.sum((2 * t2 - t2.transpose(0, 3, 2, 1)) * o2))
    return Bn, en, E2, E2t, k, occ


def mo_B(L, C):
    """B[P] = C^T L[P] C, in chunks (returns a new array (Nch, n, n))"""
    out = np.empty((len(L), C.shape[1], C.shape[1]))
    for p0 in range(0, len(L), 64):
        out[p0:p0 + 64] = np.einsum("Pmn,mp,nq->Ppq", L[p0:p0 + 64], C, C, optimize=True)
    return out

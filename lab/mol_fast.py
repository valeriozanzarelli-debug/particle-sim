"""mol_fast.py - numba kernel for the ECG elements of mol_ecg.System when all centres lie on the z axis
(same formulas as mol_ecg.System.raw/sym; checked against them in mol_ecg.py job 'check')."""
import numpy as np
import numba as nb
from math import erf, sqrt, exp, pi

SQ2 = sqrt(2.0)
SQ2PI = sqrt(2.0 / pi)


@nb.njit(cache=True, fastmath=False)
def _coul(m, s2):
    s = sqrt(s2)
    if m > 1e-10 * s:
        return erf(m / (SQ2 * s)) / m
    return SQ2PI / s


@nb.njit(cache=True)
def _pair(A1, s1, A2, s2, Z, Xz, m, Enn, C, L, w, R0, Sig):
    N = A1.shape[0]
    for i in range(N):
        for j in range(N):
            C[i, j] = A1[i, j] + A2[i, j]
    # Cholesky C = L L^T
    for i in range(N):
        for j in range(i + 1):
            acc = C[i, j]
            for k in range(j):
                acc -= L[i, k] * L[j, k]
            if i == j:
                if acc <= 0:
                    return 0.0, 0.0
                L[i, i] = sqrt(acc)
            else:
                L[i, j] = acc / L[j, j]
    det = 1.0
    for i in range(N):
        det *= L[i, i] * L[i, i]
    # Ci = C^-1 (via L^-1)
    Li = np.zeros((N, N))
    for i in range(N):
        Li[i, i] = 1.0 / L[i, i]
        for j in range(i):
            acc = 0.0
            for k in range(j, i):
                acc -= L[i, k] * Li[k, j]
            Li[i, j] = acc / L[i, i]
    for i in range(N):
        for j in range(N):
            acc = 0.0
            for k in range(max(i, j), N):
                acc += Li[k, i] * Li[k, j]
            Sig[i, j] = 0.5 * acc            # Sig = Ci / 2
    # v = A1 s1 + A2 s2 ; R0 = Ci v
    for i in range(N):
        acc = 0.0
        for j in range(N):
            acc += A1[i, j] * s1[j] + A2[i, j] * s2[j]
        w[i] = acc
    ex = 0.0
    for i in range(N):
        acc = 0.0
        for j in range(N):
            acc += 2.0 * Sig[i, j] * w[j]
        R0[i] = acc
        ex += w[i] * acc
    for i in range(N):
        for j in range(N):
            ex -= s1[i] * A1[i, j] * s1[j] + s2[i] * A2[i, j] * s2[j]
    S = (pi ** N / det) ** 1.5 * exp(ex)
    # kinetic: 4 (3 tr(A1 A2 Sig) + d1^T A1 A2 d2) / 2m
    tr = 0.0
    dd = 0.0
    for i in range(N):
        for j in range(N):
            mij = 0.0
            for k in range(N):
                mij += A1[i, k] * A2[k, j]
            tr += mij * Sig[j, i]
            dd += (R0[i] - s1[i]) * mij * (R0[j] - s2[j])
    T = 4.0 * (3.0 * tr + dd) / (2.0 * m)
    V = Enn
    for i in range(N):
        for n in range(Z.shape[0]):
            V -= Z[n] * _coul(abs(R0[i] - Xz[n]), Sig[i, i])
        for j in range(i + 1, N):
            V += _coul(abs(R0[i] - R0[j]), Sig[i, i] + Sig[j, j] - 2 * Sig[i, j])
    return S, (T + V) * S


@nb.njit(cache=True)
def rows(A1, s1, A2p, s2p, cops, Z, Xz, m, Enn):
    """A1 (n,N,N), s1 (n,N) z-centres; A2p (nops,K,N,N), s2p (nops,K,N) permuted basis; returns S, H (n, K)"""
    n, N = A1.shape[0], A1.shape[1]
    nops, K = A2p.shape[0], A2p.shape[1]
    So = np.zeros((n, K)); Ho = np.zeros((n, K))
    C = np.zeros((N, N)); L = np.zeros((N, N)); w = np.zeros(N); R0 = np.zeros(N); Sig = np.zeros((N, N))
    for a in range(n):
        for b in range(K):
            ss = 0.0; hh = 0.0
            for o in range(nops):
                s_, h_ = _pair(A1[a], s1[a], A2p[o, b], s2p[o, b], Z, Xz, m, Enn, C, L, w, R0, Sig)
                ss += cops[o] * s_; hh += cops[o] * h_
            So[a, b] = ss; Ho[a, b] = hh
    return So, Ho


@nb.njit(cache=True)
def diag(A1, s1, Pm, cops, Z, Xz, m, Enn):
    """<phi_a| P |phi_a> for each candidate (self elements)"""
    n, N = A1.shape[0], A1.shape[1]
    nops = Pm.shape[0]
    So = np.zeros(n); Ho = np.zeros(n)
    C = np.zeros((N, N)); L = np.zeros((N, N)); w = np.zeros(N); R0 = np.zeros(N); Sig = np.zeros((N, N))
    Ap = np.zeros((N, N)); sp = np.zeros(N)
    for a in range(n):
        for o in range(nops):
            P = Pm[o]
            # Ap = P^T A P ; sp = P^T s
            for i in range(N):
                acc = 0.0
                for j in range(N):
                    acc += P[j, i] * s1[a, j]
                sp[i] = acc
                for k in range(N):
                    acc2 = 0.0
                    for j in range(N):
                        for l in range(N):
                            acc2 += P[j, i] * A1[a, j, l] * P[l, k]
                    Ap[i, k] = acc2
            s_, h_ = _pair(A1[a], s1[a], Ap, sp, Z, Xz, m, Enn, C, L, w, R0, Sig)
            So[a] += cops[o] * s_; Ho[a] += cops[o] * h_
    return So, Ho

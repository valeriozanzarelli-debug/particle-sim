"""acqua_int.py - motore di integrali gaussiani a piu' centri (studio ACQUA).

Primitive gaussiane NON contratte, sferiche reali: phi = N r^l Y_lm(r^) exp(-a r^2), con r^l Y_lm scritto come
combinazione di monomi cartesiani (coefficienti ricavati qui, per minimi quadrati esatti su punti a caso: niente
tabelle). Integrali con McMurchie-Davidson (coefficienti di Hermite E, integrali R con la funzione di Boys),
compilati con numba: sovrapposizione, cinetica, attrazione nucleare, dipolo, repulsione a due elettroni (simmetria
8 volte, screening di Schwarz 1e-13). Unita' atomiche (bohr, Hartree).
"""
import os, math
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from numba import njit
from scipy.special import sph_harm_y, gamma

LMAX = 4


def cart_list(L):
    return [(a, b, L - a - b) for a in range(L, -1, -1) for b in range(L - a, -1, -1)]


def _sph_coef(L):
    """C[m, k]: r^L Y_Lm(real, orthonormal on the sphere) = sum_k C[m,k] x^ax y^ay z^az"""
    cl = cart_list(L); rng = np.random.default_rng(7 + L)
    n = 4 * len(cl) + 20
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    th = np.arccos(np.clip(v[:, 2], -1, 1)); ph = np.arctan2(v[:, 1], v[:, 0])
    A = np.array([[x ** a * y ** b * z ** c for (a, b, c) in cl] for x, y, z in v])
    C = np.zeros((2 * L + 1, len(cl)))
    for i, m in enumerate(range(-L, L + 1)):
        Y = sph_harm_y(L, abs(m), th, ph)
        y = Y.real if m == 0 else (math.sqrt(2) * (-1) ** m * (Y.real if m > 0 else Y.imag))
        C[i] = np.linalg.lstsq(A, y, rcond=None)[0]
        assert np.abs(A @ C[i] - y).max() < 1e-10
    return C


CART = np.zeros((LMAX + 1, (LMAX + 1) * (LMAX + 2) // 2, 3), dtype=np.int64)
SPH = np.zeros((LMAX + 1, 2 * LMAX + 1, (LMAX + 1) * (LMAX + 2) // 2))
for _L in range(LMAX + 1):
    CART[_L, :len(cart_list(_L))] = np.array(cart_list(_L))
    SPH[_L, :2 * _L + 1, :len(cart_list(_L))] = _sph_coef(_L)


def norm_prim(L, a):
    return 1.0 / math.sqrt(gamma(L + 1.5) / (2.0 * (2.0 * a) ** (L + 1.5)))


class Basis:
    """shells: list of (center index, L, exponent); centers: (n,3) bohr; Z: charges"""

    def __init__(self, centers, Z, shells):
        self.R = np.asarray(centers, float); self.Z = np.asarray(Z, float)
        self.shells = list(shells)
        self.sL = np.array([s[1] for s in shells], dtype=np.int64)
        self.sa = np.array([s[2] for s in shells], float)
        self.sC = np.array([self.R[s[0]] for s in shells], float)
        self.sN = np.array([norm_prim(s[1], s[2]) for s in shells])
        self.off = np.concatenate([[0], np.cumsum(2 * self.sL + 1)]).astype(np.int64)
        self.M = int(self.off[-1])
        self.atom_of = np.concatenate([[s[0]] * (2 * s[1] + 1) for s in shells]).astype(int)

    def enuc(self):
        E = 0.0
        for i in range(len(self.Z)):
            for j in range(i):
                E += self.Z[i] * self.Z[j] / np.linalg.norm(self.R[i] - self.R[j])
        return E


# ------------------------------------------------------------------ numba kernels
@njit(cache=True)
def boys(nmax, T, F):
    if T < 1e-13:
        for n in range(nmax + 1):
            F[n] = 1.0 / (2 * n + 1) - T / (2 * n + 3)
        return
    if T > 45.0:
        eT = math.exp(-T)
        F[0] = 0.5 * math.sqrt(math.pi / T) * math.erf(math.sqrt(T))
        for n in range(nmax):
            F[n + 1] = ((2 * n + 1) * F[n] - eT) / (2 * T)
        return
    eT = math.exp(-T)
    term = 1.0 / (2 * nmax + 1); s = term; k = 1
    while True:
        term *= 2 * T / (2 * nmax + 2 * k + 1); s += term; k += 1
        if term < 1e-17 * s:
            break
    F[nmax] = eT * s
    for n in range(nmax - 1, -1, -1):
        F[n] = (2 * T * F[n + 1] + eT) / (2 * n + 1)


@njit(cache=True)
def hermE(la, lb, a, b, XAB, E):
    p = a + b; mu = a * b / p; XPA = -b / p * XAB; XPB = a / p * XAB; h = 0.5 / p
    for i in range(la + 1):
        for j in range(lb + 1):
            for t in range(la + lb + 1):
                E[i, j, t] = 0.0
    E[0, 0, 0] = math.exp(-mu * XAB * XAB)
    for i in range(la + 1):
        for j in range(lb + 1):
            if i == 0 and j == 0:
                continue
            if i > 0:
                i0, j0, X = i - 1, j, XPA
            else:
                i0, j0, X = i, j - 1, XPB
            for t in range(i + j + 1):
                v = X * E[i0, j0, t]
                if t > 0:
                    v += h * E[i0, j0, t - 1]
                if t + 1 <= i0 + j0:
                    v += (t + 1) * E[i0, j0, t + 1]
                E[i, j, t] = v


@njit(cache=True)
def hermR(L, alpha, X, Y, Z, R, F):
    """R[0, t, u, v] for t+u+v <= L (aux index n in first axis)"""
    boys(L, alpha * (X * X + Y * Y + Z * Z), F)
    f = 1.0
    for n in range(L + 1):
        R[n, 0, 0, 0] = f * F[n]; f *= -2 * alpha
    for N in range(1, L + 1):
        for n in range(L - N + 1):
            for t in range(N + 1):
                for u in range(N - t + 1):
                    v = N - t - u
                    if t > 0:
                        val = X * R[n + 1, t - 1, u, v]
                        if t > 1:
                            val += (t - 1) * R[n + 1, t - 2, u, v]
                    elif u > 0:
                        val = Y * R[n + 1, t, u - 1, v]
                        if u > 1:
                            val += (u - 1) * R[n + 1, t, u - 2, v]
                    else:
                        val = Z * R[n + 1, t, u, v - 1]
                        if v > 1:
                            val += (v - 1) * R[n + 1, t, u, v - 2]
                    R[n, t, u, v] = val


@njit(cache=True)
def _sph2(blk, la, lb, SPH, out):
    """cartesian (na,nb) -> spherical (2la+1, 2lb+1)"""
    na = (la + 1) * (la + 2) // 2; nb = (lb + 1) * (lb + 2) // 2
    for m in range(2 * la + 1):
        for n in range(2 * lb + 1):
            s = 0.0
            for i in range(na):
                ci = SPH[la, m, i]
                if ci == 0.0:
                    continue
                for j in range(nb):
                    s += ci * SPH[lb, n, j] * blk[i, j]
            out[m, n] = s


@njit(cache=True)
def _kin(E, i_, j_, b, sq):
    v_ = -2 * b * b * E[i_, j_ + 2, 0] + b * (2 * j_ + 1) * E[i_, j_, 0]
    if j_ >= 2:
        v_ += -0.5 * j_ * (j_ - 1) * E[i_, j_ - 2, 0]
    return v_ * sq


@njit(cache=True)
def one_e(sL, sa, sC, sN, off, Zn, Rn, CART, SPH, S, T, V, D):
    ns = len(sL); L4 = 2 * 4 + 3
    Ex = np.zeros((7, 7, 14)); Ey = np.zeros((7, 7, 14)); Ez = np.zeros((7, 7, 14))
    R = np.zeros((L4, L4, L4, L4)); F = np.zeros(L4)
    bs = np.zeros((15, 15)); bt = np.zeros((15, 15)); bv = np.zeros((15, 15))
    bd = np.zeros((3, 15, 15)); o = np.zeros((9, 9))
    for A in range(ns):
        for B in range(A + 1):
            la, lb = sL[A], sL[B]; a, b = sa[A], sa[B]; p = a + b
            hermE(la, lb + 2, a, b, sC[A, 0] - sC[B, 0], Ex)
            hermE(la, lb + 2, a, b, sC[A, 1] - sC[B, 1], Ey)
            hermE(la, lb + 2, a, b, sC[A, 2] - sC[B, 2], Ez)
            P0 = (a * sC[A, 0] + b * sC[B, 0]) / p; P1 = (a * sC[A, 1] + b * sC[B, 1]) / p
            P2 = (a * sC[A, 2] + b * sC[B, 2]) / p
            sq = math.sqrt(math.pi / p)
            na = (la + 1) * (la + 2) // 2; nb = (lb + 1) * (lb + 2) // 2
            for i in range(na):
                for j in range(nb):
                    bv[i, j] = 0.0
            for c in range(len(Zn)):
                hermR(la + lb, p, P0 - Rn[c, 0], P1 - Rn[c, 1], P2 - Rn[c, 2], R, F)
                for i in range(na):
                    ax, ay, az = CART[la, i, 0], CART[la, i, 1], CART[la, i, 2]
                    for j in range(nb):
                        bx, by, bz = CART[lb, j, 0], CART[lb, j, 1], CART[lb, j, 2]
                        s = 0.0
                        for t in range(ax + bx + 1):
                            for u in range(ay + by + 1):
                                for v in range(az + bz + 1):
                                    s += Ex[ax, bx, t] * Ey[ay, by, u] * Ez[az, bz, v] * R[0, t, u, v]
                        bv[i, j] += -Zn[c] * 2 * math.pi / p * s
            for i in range(na):
                ax, ay, az = CART[la, i, 0], CART[la, i, 1], CART[la, i, 2]
                for j in range(nb):
                    bx, by, bz = CART[lb, j, 0], CART[lb, j, 1], CART[lb, j, 2]
                    sx = Ex[ax, bx, 0] * sq; sy = Ey[ay, by, 0] * sq; sz = Ez[az, bz, 0] * sq
                    bs[i, j] = sx * sy * sz

                    bt[i, j] = _kin(Ex, ax, bx, b, sq) * sy * sz + sx * _kin(Ey, ay, by, b, sq) * sz \
                        + sx * sy * _kin(Ez, az, bz, b, sq)
                    # dipole about the origin: x = x_P + P_x
                    bd[0, i, j] = (Ex[ax, bx, 1] + P0 * Ex[ax, bx, 0]) * sq * sy * sz
                    bd[1, i, j] = sx * (Ey[ay, by, 1] + P1 * Ey[ay, by, 0]) * sq * sz
                    bd[2, i, j] = sx * sy * (Ez[az, bz, 1] + P2 * Ez[az, bz, 0]) * sq
            nrm = sN[A] * sN[B]
            for q in range(6):
                if q == 0:
                    _sph2(bs, la, lb, SPH, o)
                elif q == 1:
                    _sph2(bt, la, lb, SPH, o)
                elif q == 2:
                    _sph2(bv, la, lb, SPH, o)
                else:
                    _sph2(bd[q - 3], la, lb, SPH, o)
                for m in range(2 * la + 1):
                    for n in range(2 * lb + 1):
                        x = o[m, n] * nrm; I = off[A] + m; J = off[B] + n
                        if q == 0:
                            S[I, J] = x; S[J, I] = x
                        elif q == 1:
                            T[I, J] = x; T[J, I] = x
                        elif q == 2:
                            V[I, J] = x; V[J, I] = x
                        else:
                            D[q - 3, I, J] = x; D[q - 3, J, I] = x


@njit(cache=True)
def _pairE(la, lb, a, b, CA, CB, CART, Eab):
    """Eab[i, j, t, u, v] for cart comps i of A, j of B; returns P"""
    Ex = np.zeros((5, 5, 9)); Ey = np.zeros((5, 5, 9)); Ez = np.zeros((5, 5, 9))
    hermE(la, lb, a, b, CA[0] - CB[0], Ex); hermE(la, lb, a, b, CA[1] - CB[1], Ey)
    hermE(la, lb, a, b, CA[2] - CB[2], Ez)
    na = (la + 1) * (la + 2) // 2; nb = (lb + 1) * (lb + 2) // 2; L = la + lb
    for i in range(na):
        for j in range(nb):
            for t in range(L + 1):
                for u in range(L + 1 - t):
                    for v in range(L + 1 - t - u):
                        Eab[i, j, t, u, v] = Ex[CART[la, i, 0], CART[lb, j, 0], t] * Ey[CART[la, i, 1], CART[lb, j, 1], u] \
                            * Ez[CART[la, i, 2], CART[lb, j, 2], v]


@njit(cache=True)
def _quartet(A, B, C, D, sL, sa, sC, CART, SPH, Eab, Ecd, R, F, Gk, blk, t1, t2, t3, out):
    la, lb, lc, ld = sL[A], sL[B], sL[C], sL[D]
    a, b, c, d = sa[A], sa[B], sa[C], sa[D]
    p = a + b; q = c + d; alpha = p * q / (p + q)
    P = (a * sC[A] + b * sC[B]) / p; Q = (c * sC[C] + d * sC[D]) / q
    Lab = la + lb; Lcd = lc + ld; L = Lab + Lcd
    hermR(L, alpha, P[0] - Q[0], P[1] - Q[1], P[2] - Q[2], R, F)
    na = (la + 1) * (la + 2) // 2; nb = (lb + 1) * (lb + 2) // 2
    nc = (lc + 1) * (lc + 2) // 2; nd = (ld + 1) * (ld + 2) // 2
    pref = 2 * math.pi ** 2.5 / (p * q * math.sqrt(p + q))
    for k in range(nc):
        for l in range(nd):
            # Gk[t,u,v] = sum_{tau nu phi} (-1)^(tau+nu+phi) Ecd R[t+tau, ...]
            for t in range(Lab + 1):
                for u in range(Lab + 1 - t):
                    for v in range(Lab + 1 - t - u):
                        s = 0.0
                        for t2_ in range(Lcd + 1):
                            for u2 in range(Lcd + 1 - t2_):
                                for v2 in range(Lcd + 1 - t2_ - u2):
                                    e = Ecd[k, l, t2_, u2, v2]
                                    if e == 0.0:
                                        continue
                                    sg = 1.0 if (t2_ + u2 + v2) % 2 == 0 else -1.0
                                    s += sg * e * R[0, t + t2_, u + u2, v + v2]
                        Gk[t, u, v] = s
            for i in range(na):
                for j in range(nb):
                    s = 0.0
                    for t in range(Lab + 1):
                        for u in range(Lab + 1 - t):
                            for v in range(Lab + 1 - t - u):
                                s += Eab[i, j, t, u, v] * Gk[t, u, v]
                    blk[i, j, k, l] = pref * s
    # cartesian -> spherical, one index at a time
    ma, mb, mc, md = 2 * la + 1, 2 * lb + 1, 2 * lc + 1, 2 * ld + 1
    for m in range(ma):
        for j in range(nb):
            for k in range(nc):
                for l in range(nd):
                    s = 0.0
                    for i in range(na):
                        s += SPH[la, m, i] * blk[i, j, k, l]
                    t1[m, j, k, l] = s
    for m in range(ma):
        for n in range(mb):
            for k in range(nc):
                for l in range(nd):
                    s = 0.0
                    for j in range(nb):
                        s += SPH[lb, n, j] * t1[m, j, k, l]
                    t2[m, n, k, l] = s
    for m in range(ma):
        for n in range(mb):
            for o in range(mc):
                for l in range(nd):
                    s = 0.0
                    for k in range(nc):
                        s += SPH[lc, o, k] * t2[m, n, k, l]
                    t3[m, n, o, l] = s
    for m in range(ma):
        for n in range(mb):
            for o in range(mc):
                for r in range(md):
                    s = 0.0
                    for l in range(nd):
                        s += SPH[ld, r, l] * t3[m, n, o, l]
                    out[m, n, o, r] = s


@njit(cache=True)
def eri_full(sL, sa, sC, sN, off, CART, SPH, G, thr):
    ns = len(sL)
    Eab = np.zeros((15, 15, 9, 9, 9)); Ecd = np.zeros((15, 15, 9, 9, 9))
    R = np.zeros((17, 17, 17, 17)); F = np.zeros(17); Gk = np.zeros((9, 9, 9))
    blk = np.zeros((15, 15, 15, 15)); t1 = np.zeros((9, 15, 15, 15)); t2 = np.zeros((9, 9, 15, 15))
    t3 = np.zeros((9, 9, 9, 15)); out = np.zeros((9, 9, 9, 9))
    # Schwarz: Q[A,B] = max |(ab|ab)|^0.5
    Qs = np.zeros((ns, ns))
    for A in range(ns):
        for B in range(A + 1):
            _pairE(sL[A], sL[B], sa[A], sa[B], sC[A], sC[B], CART, Eab)
            _quartet(A, B, A, B, sL, sa, sC, CART, SPH, Eab, Eab, R, F, Gk, blk, t1, t2, t3, out)
            mx = 0.0
            for m in range(2 * sL[A] + 1):
                for n in range(2 * sL[B] + 1):
                    mx = max(mx, abs(out[m, n, m, n]))
            Qs[A, B] = Qs[B, A] = math.sqrt(mx) * sN[A] * sN[B]
    for A in range(ns):
        for B in range(A + 1):
            _pairE(sL[A], sL[B], sa[A], sa[B], sC[A], sC[B], CART, Eab)
            for C in range(A + 1):
                Dmax = B if C == A else C
                for D in range(Dmax + 1):
                    if Qs[A, B] * Qs[C, D] < thr:
                        continue
                    _pairE(sL[C], sL[D], sa[C], sa[D], sC[C], sC[D], CART, Ecd)
                    _quartet(A, B, C, D, sL, sa, sC, CART, SPH, Eab, Ecd, R, F, Gk, blk, t1, t2, t3, out)
                    nrm = sN[A] * sN[B] * sN[C] * sN[D]
                    for m in range(2 * sL[A] + 1):
                        I = off[A] + m
                        for n in range(2 * sL[B] + 1):
                            J = off[B] + n
                            for o in range(2 * sL[C] + 1):
                                K = off[C] + o
                                for r in range(2 * sL[D] + 1):
                                    L_ = off[D] + r
                                    x = out[m, n, o, r] * nrm
                                    G[I, J, K, L_] = x; G[J, I, K, L_] = x; G[I, J, L_, K] = x; G[J, I, L_, K] = x
                                    G[K, L_, I, J] = x; G[L_, K, I, J] = x; G[K, L_, J, I] = x; G[L_, K, J, I] = x


def integrals(B, eri=True):
    M = B.M
    S = np.zeros((M, M)); T = np.zeros((M, M)); V = np.zeros((M, M)); D = np.zeros((3, M, M))
    one_e(B.sL, B.sa, B.sC, B.sN, B.off, B.Z, B.R, CART, SPH, S, T, V, D)
    out = dict(S=S, T=T, V=V, D=D)
    if eri:
        G = np.zeros((M, M, M, M))
        eri_full(B.sL, B.sa, B.sC, B.sN, B.off, CART, SPH, G, 1e-13)
        out["G"] = G
    return out

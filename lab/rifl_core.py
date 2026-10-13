"""rifl_core.py - RIFLUSSO QUANTISTICO (quantum backflow) sulla regola zig-zag 1+1D.

Il problema (Bracken & Melloy 1994): tra gli stati con solo impulsi positivi (e, per Dirac, energia
positiva) qual e' la massima probabilita' che passa da destra a sinistra attraverso un punto in un tempo T?

Regola zig-zag 1-D (zz_core.Walk kind='dirac1'): per ogni tick
    R'(x) = cos th R(x-1) - i sin th L(x+1),   L'(x) = cos th L(x+1) - i sin th R(x-1)
(prima corri a c: R va a +1, L a -1; poi ruota). th = m c a / hbar, tick dt = a / c.
Flusso esatto attraverso il legame (0|1) nel tick n: J_n = |R_n(0)|^2 - |L_n(1)|^2 (prima dello spostamento:
R in 0 passa in 1, L in 1 passa in 0; la rotazione e' locale e non sposta probabilita').
Riflusso in N tick: Delta_N = -sum_{n<N} J_n = P_sinistra(N) - P_sinistra(0)   (sinistra = siti <= 0).

Onde piane: psi(x) = u(k) e^{ikx}, U(k) = [[c e^{-ik}, -i s e^{ik}], [-i s e^{-ik}, c e^{ik}]],
autovalori e^{-+i w}, cos w = cos th cos k. Banda positiva: e^{-i w}, w in (0, pi).
Stato psi_0 = int_W dk f(k) u(k) e^{ikx} / sqrt(2 pi), int |f|^2 = 1. Allora Delta_N = <f|K|f> con
    K(k,k') = (1/2pi) S_N(w_k - w_k') [ conj(uL_k) uL_k' e^{-i(k-k')} - conj(uR_k) uR_k' ],
    S_N(d) = sum_{n=0}^{N-1} e^{i d n}.
Il massimo riflusso e' l'autovalore piu' grande di K ristretto alla finestra W (quadratura di Gauss-Legendre).

Riferimenti di continuo (SOLO esami, mai ingressi):
  NR (Bracken-Melloy): K(u,v) = -(1/pi) sin((u^2-v^2)/2)/(u-v), u = p sqrt(T/2m) in (0,U);
     c_bm = lim_{U->inf} = 0.0384517 (Penz et al. 2006).
  Dirac 1+1 (Melloy-Bracken 1998): H = [[p, m],[m, -p]] (stessa rappresentazione chirale della regola),
     K(p,p') = (1/pi) sin(d T/2)/d (L L' - R R'), d = E - E'.  eps = sqrt(4 hbar/(m c^2 T)).
Unita': hbar = c = 1. Nella regola: m = th/a, T = N a.
"""
import numpy as np
import scipy.linalg as sl
import scipy.sparse.linalg as spl

C_BM = 0.0384517          # Penz, Gruebl, Kreidl, Wagner, J. Phys. A 39, 423 (2006): SOLO esame


def gl_composite(lo, hi, nint, npt):
    """Gauss-Legendre composta: nint sottointervalli uguali, npt punti ciascuno"""
    x, w = np.polynomial.legendre.leggauss(npt)
    edges = np.linspace(lo, hi, nint + 1)
    a, b = edges[:-1, None], edges[1:, None]
    X = (a + b) / 2 + (b - a) / 2 * x[None, :]
    W = (b - a) / 2 * w[None, :]
    return X.ravel(), W.ravel()


def top_eig(A, vec=False):
    """autovalore piu' grande (e autovettore) di una matrice hermitiana densa"""
    M = A.shape[0]
    if M <= 3000:
        w, v = sl.eigh(A, subset_by_index=[M - 1, M - 1])
        return (w[0], v[:, 0]) if vec else w[0]
    w, v = spl.eigsh(A, k=1, which="LA", tol=1e-13, maxiter=20000)
    return (w[0], v[:, 0]) if vec else w[0]


# ---------------------------------------------------------------- regola zig-zag
def walk_modes(k, theta):
    """banda positiva della regola: w(k), uR(k), uL(k) (autovettore di U(k) con autovalore e^{-iw})"""
    c, s = np.cos(theta), np.sin(theta)
    w = np.arccos(np.clip(c * np.cos(k), -1, 1))
    lam = np.exp(-1j * w)
    # (U - lam) u = 0. Due scelte di riga; uso quella con norma maggiore (stabile)
    a1 = 1j * s * np.exp(1j * k)                     # uR da riga 1
    b1 = c * np.exp(-1j * k) - lam                   # uL da riga 1
    a2 = c * np.exp(1j * k) - lam                    # uR da riga 2: (-i s e^{-ik}) uR + (c e^{ik}-lam) uL = 0
    b2 = 1j * s * np.exp(-1j * k)
    n1 = np.abs(a1)**2 + np.abs(b1)**2
    n2 = np.abs(a2)**2 + np.abs(b2)**2
    use1 = n1 >= n2
    uR = np.where(use1, a1, a2)
    uL = np.where(use1, b1, b2)
    nrm = np.sqrt(np.abs(uR)**2 + np.abs(uL)**2)
    return w, uR / nrm, uL / nrm


def S_N(d, N):
    """sum_{n=0}^{N-1} e^{i d n}, stabile per d -> 0"""
    d = np.asarray(d, float)
    half = d / 2
    small = np.abs(np.sin(half)) < 1e-12
    with np.errstate(all="ignore"):
        r = np.exp(1j * half * (N - 1)) * np.sin(N * half) / np.sin(half)
    return np.where(small, N + 0j, r)


def walk_kernel(k, wq, theta, N):
    w, uR, uL = walk_modes(k, theta)
    d = w[:, None] - w[None, :]
    g = (np.conj(uL)[:, None] * uL[None, :] * np.exp(-1j * (k[:, None] - k[None, :]))
         - np.conj(uR)[:, None] * uR[None, :])
    K = S_N(d, N) * g / (2 * np.pi)
    sw = np.sqrt(wq)
    A = sw[:, None] * K * sw[None, :]
    return (A + A.conj().T) / 2


def walk_backflow(theta, N, kmax, nint, npt=24, kmin=0.0, vec=False):
    k, wq = gl_composite(kmin, kmax, nint, npt)
    A = walk_kernel(k, wq, theta, N)
    r = top_eig(A, vec)
    if vec:
        lam, v = r
        return lam, k, wq, v / np.sqrt(wq)
    return r


# ---------------------------------------------------------------- continuo (solo esami)
def nr_backflow(U, nint, npt=24):
    u, w = gl_composite(0, U, nint, npt)
    d = u[:, None] - u[None, :]
    s = (u[:, None]**2 - u[None, :]**2) / 2
    with np.errstate(all="ignore"):
        K = -np.sin(s) / d / np.pi
    K[np.arange(len(u)), np.arange(len(u))] = -u / np.pi
    sw = np.sqrt(w)
    return top_eig(sw[:, None] * K * sw[None, :])


def dirac_modes(p, m=1.0):
    E = np.sqrt(p**2 + m**2)
    R = np.sqrt((1 + p / E) / 2)
    L = np.sqrt((1 - p / E) / 2)
    return E, R, L


def dirac_backflow(T, pmax, nint, npt=24, m=1.0, pmin=0.0):
    p, w = gl_composite(pmin, pmax, nint, npt)
    E, R, L = dirac_modes(p, m)
    d = E[:, None] - E[None, :]
    with np.errstate(all="ignore"):
        F = np.sin(d * T / 2) / d
    F[np.abs(d) < 1e-14] = T / 2
    K = F * (L[:, None] * L[None, :] - R[:, None] * R[None, :]) / np.pi
    sw = np.sqrt(w)
    return top_eig(sw[:, None] * K * sw[None, :])


def eps_of(theta, N):
    """eps = sqrt(4 hbar / (m c^2 T)) con m = theta/a, T = N a"""
    return np.sqrt(4.0 / (theta * N))


# ---------------------------------------------------------------- controllo con la regola vera
def lattice_check(theta, N, f_k, k, Lring):
    """costruisce psi_0 dalla f(k) su un anello di Lring siti (k sulla griglia 2 pi j / Lring),
    lo fa evolvere con zz_core.Walk('dirac1') per N tick e somma il flusso verso sinistra al taglio (0|1)."""
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from zz_core import Walk
    x = np.arange(Lring) - Lring // 2          # sito 0 al centro; sinistra = x <= 0
    _, uR, uL = walk_modes(k, theta)
    ph = np.exp(1j * np.outer(x, k))
    psi = np.array([ph @ (f_k * uR), ph @ (f_k * uL)])
    psi /= np.sqrt(np.sum(np.abs(psi)**2))
    # Walk dirac1 usa np.roll sull'array: R(x) <- R(x-1): coerente con la convenzione sopra
    W = Walk("dirac1", Lring, 1.0, dim=1, theta=theta)
    i0 = int(np.where(x == 0)[0][0])
    back = 0.0
    for _ in range(N):
        # flusso del tick: L in 1 va in 0 (verso sinistra), R in 0 va in 1 (verso destra).
        # Si misura solo il taglio (0|1): l'anello ha un secondo taglio lontano che non conta.
        back += abs(psi[1, i0 + 1])**2 - abs(psi[0, i0])**2
        psi = W.tick(psi)
    return back

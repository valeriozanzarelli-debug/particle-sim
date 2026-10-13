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


def gl_adaptive(lo, hi, rate, cosc=6.0, npt=24, hmax=None):
    """Gauss-Legendre a pannelli: lunghezza del pannello h(x) = cosc / rate(x) (rate = fase per unita' di x
    del nucleo oscillante), cosi' ogni pannello contiene ~cosc/(2 pi) oscillazioni. Raddoppiare la
    risoluzione = dimezzare cosc."""
    edges = [lo]
    x = lo
    while x < hi:
        h = cosc / max(rate(x), 1e-300)
        if hmax is not None:
            h = min(h, hmax)
        # correzione: la fase cresce dentro il pannello; uso il ritmo al bordo destro (difetto trovato
        # in E2: il primo pannello era troppo lungo vicino a k = 0)
        for _ in range(4):
            h = min(h, cosc / max(rate(min(hi, x + h)), 1e-300))
        x = min(hi, x + h)
        edges.append(x)
    edges = np.array(edges)
    xg, wg = np.polynomial.legendre.leggauss(npt)
    a, b = edges[:-1, None], edges[1:, None]
    X = (a + b) / 2 + (b - a) / 2 * xg[None, :]
    W = (b - a) / 2 * wg[None, :]
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
    # w = arccos(c cos k) in forma stabile per theta, k piccoli: 1 - c cos k = 2 sin^2(th/2) + 2 c sin^2(k/2)
    w = 2 * np.arcsin(np.sqrt(np.clip(np.sin(theta / 2)**2 + c * np.sin(k / 2)**2, 0, 1)))
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


def walk_kernel(k, wq, theta, N, chunk=512):
    """matrice A = sqrt(w) K sqrt(w), costruita a blocchi di righe (memoria: una sola M x M complessa)"""
    w, uR, uL = walk_modes(k, theta)
    M = len(k)
    sw = np.sqrt(wq)
    aL = uL * np.exp(1j * k) * sw          # conj(uL_k) e^{-ik} ... : g = conj(aL_k) aL_k' - conj(aR_k) aR_k'
    aR = uR * sw
    A = np.empty((M, M), complex)
    for i0 in range(0, M, chunk):
        i1 = min(M, i0 + chunk)
        d = w[i0:i1, None] - w[None, :]
        A[i0:i1] = S_N(d, N) * (np.conj(aL[i0:i1, None]) * aL[None, :] - np.conj(aR[i0:i1, None]) * aR[None, :])
    A /= 2 * np.pi
    A += A.conj().T
    A /= 2
    return A


def walk_grid(theta, N, kmax, cosc=6.0, kmin=0.0, npt=24):
    """griglia adattiva in k: fase del nucleo = v_g(k) N / 2 per unita' di k (v_g = dw/dk)"""
    c = np.cos(theta)
    def rate(k):
        w = np.arccos(np.clip(c * np.cos(k), -1, 1))
        vg = c * np.sin(k) / max(np.sin(w), 1e-300)
        return max(vg * N / 2, 1.0, np.sqrt(N / (2 * np.tan(max(theta, 1e-12)))) * 0.5)
    return gl_adaptive(kmin, kmax, rate, cosc, npt, hmax=(kmax - kmin) / 4)


def walk_backflow(theta, N, kmax, nint=None, npt=24, kmin=0.0, vec=False, cosc=6.0):
    """nint = None: griglia adattiva (cosc); altrimenti Gauss-Legendre composta uniforme"""
    if nint is None:
        k, wq = walk_grid(theta, N, kmax, cosc, kmin, npt)
    else:
        k, wq = gl_composite(kmin, kmax, nint, npt)
    A = walk_kernel(k, wq, theta, N)
    r = top_eig(A, vec)
    if vec:
        lam, v = r
        return lam, k, wq, v / np.sqrt(wq)
    return r


# ---------------------------------------------------------------- continuo (solo esami)
def nr_backflow(U, nint=None, npt=24, cosc=6.0):
    if nint is None:
        u, w = gl_adaptive(0, U, lambda x: max(x, 1.0), cosc, npt)
    else:
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


def dirac_grid(T, pmax, cosc=6.0, m=1.0, pmin=0.0, npt=24, geo=False):
    """geo=True (per eps grandi, T << 1/m): sotto p ~ 1/T il nucleo non oscilla; i pannelli crescono
    in modo geometrico (ritmo 2/(p+m)) invece di avere lunghezza fissa."""
    def rate(p):
        v = p / np.sqrt(p * p + m * m)
        if geo:
            return max(v * T / 2, 2.0 / (p + m))
        return max(v * T / 2, 1.0, 0.5 * np.sqrt(T / (2 * m)))
    return gl_adaptive(pmin, pmax, rate, cosc, npt, hmax=(pmax - pmin) / 4)


def dirac_backflow(T, pmax, nint=None, npt=24, m=1.0, pmin=0.0, cosc=6.0, geo=False):
    if nint is None:
        p, w = dirac_grid(T, pmax, cosc, m, pmin, npt, geo)
    else:
        p, w = gl_composite(pmin, pmax, nint, npt)
    E, R, L = dirac_modes(p, m)
    # memoria: una sola matrice M x M costruita sul posto
    K = E[:, None] - E[None, :]                      # d
    small = np.abs(K) < 1e-14
    K[small] = 1.0
    F = np.sin(K * (T / 2))
    F /= K
    del K
    F[small] = T / 2
    del small
    F *= L[:, None] * L[None, :] - R[:, None] * R[None, :]
    sw = np.sqrt(w) / np.sqrt(np.pi)
    F *= sw[:, None]
    F *= sw[None, :]
    return top_eig(F)


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


# ---------------------------------------------------------------- forma "di posizione" (Gram), per N grandi
def walk_backflow_gram(theta, N, kmax, cosc=6.0, kmin=0.0, vec=False):
    """Stesso massimo di walk_backflow, ma come problema 2N x 2N.
    B = sum_n U^-n J U^n con J = |L,1><L,1| - |R,0><R,0|  =>  P B P = X D X^dag, X = [P U^-n |a>].
    Autovalori non nulli = quelli di G^1/2 D G^1/2, G = X^dag X (Toeplitz nel ritardo n - n'):
    G_(n a),(n' b) = int_W dk/2pi e^{-i w (n-n')} c_a(k) conj(c_b(k)), c_R = uR, c_L = uL e^{ik}.
    La dimensione non dipende dalla finestra: utile per la banda intera e N fino a ~2000."""
    k, wq = walk_grid(theta, max(N, 2), kmax, cosc, kmin)
    w, uR, uL = walk_modes(k, theta)
    cR, cL = uR, uL * np.exp(1j * k)
    lags = np.arange(-(N - 1), N)
    E = np.exp(-1j * np.outer(lags, w)) * (wq / (2 * np.pi))[None, :]       # (2N-1, M)
    I = {}
    for na, ca in (("R", cR), ("L", cL)):
        for nb, cb in (("R", cR), ("L", cL)):
            I[na + nb] = E @ (ca * np.conj(cb))
    n = np.arange(N)
    idx = (n[:, None] - n[None, :]) + (N - 1)
    G = np.empty((2 * N, 2 * N), complex)
    G[:N, :N] = I["LL"][idx]
    G[:N, N:] = I["LR"][idx]
    G[N:, :N] = I["RL"][idx]
    G[N:, N:] = I["RR"][idx]
    G = (G + G.conj().T) / 2
    s, V = sl.eigh(G)
    s = np.clip(s, 0, None)
    Gh = (V * np.sqrt(s)[None, :]) @ V.conj().T
    D = np.concatenate([np.ones(N), -np.ones(N)])
    H = Gh @ (D[:, None] * Gh)
    H = (H + H.conj().T) / 2
    lam, v = sl.eigh(H, subset_by_index=[2 * N - 1, 2 * N - 1])
    return (lam[0], v[:, 0]) if vec else lam[0]

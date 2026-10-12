"""box_core.py - the zig-zag walk of lab/redo_r1.py with ONE transverse direction y of d sites and hard walls.

One tick (for one plane wave e^{i kx x} along x, field along x as a Peierls phase A):
    S_x : R <- e^{-i(kx+A)} R,  L <- e^{+i(kx+A)} L          (the R1 x-hop, sigma_z in the (R, L) basis)
    W_y : the ball hops one site in y, up if its sigma_y = +1, down if sigma_y = -1  (= e^{-i ky sigma_y} in the bulk)
          walls (AMENDED, lab/out/box_amend.json): the y-ring has d strip sites + W = 12 wall sites where the coin
          angle is th_w = 1.5 ("mass wall", same sign as th: no edge states); the first version (reflection into the
          opposite y-mover with phase i, kept as static_part_reflect) carries chiral edge states that close the gap
    C   : the R1 coin, cos th, -i sin th  (e^{-i th sigma_x}, mass m = th)
In the bulk U ~ exp(-i (th sx + kx sz + ky sy)): a 2+1 D Dirac particle with c = 1 in x and y.
Units hbar = c = a = tick = e = 1.
"""
import numpy as np

PI = np.pi


def coin(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -1j * s], [-1j * s, c]])


TH_W, W_SITES = 1.5, 12
W_CUT = 1.4          # quasi-energies |w| < W_CUT: strip states (the wall bands start at |w| = TH_W)


def static_part(th, d, thw=TH_W, w=W_SITES):
    """C . W_y on the ring of L = d + w sites (strip y = 0..d-1 with coin th, wall with coin thw); basis (R y, L y)"""
    L = d + w
    sy = np.array([[0, -1j], [1j, 0]])
    _, V = np.linalg.eigh(sy)
    N = 2 * L
    W = np.zeros((N, N), complex)
    for y in range(L):
        W[L + (y + 1) % L, L + y] = 1
        W[(y - 1) % L, y] = 1
    T = np.kron(V, np.eye(L))
    W = T @ W @ T.conj().T
    ths = np.array([th] * d + [thw] * w)
    c, s = np.cos(ths), np.sin(ths)
    C = np.block([[np.diag(c), np.diag(-1j * s)], [np.diag(-1j * s), np.diag(c)]])
    return C @ W


def static_part_reflect(th, d, rho=1j):
    """C . W_y as a (2d x 2d) matrix in the basis (R y=0..d-1, L y=0..d-1)"""
    sy = np.array([[0, -1j], [1j, 0]])
    w, V = np.linalg.eigh(sy)                      # w = -1 (down-mover), +1 (up-mover)
    N = 2 * d
    W = np.zeros((N, N), complex)                  # sigma_y basis: block 0 = down, block 1 = up
    for y in range(d):
        if y + 1 < d:
            W[d + y + 1, d + y] = 1
        else:
            W[y, d + y] = rho
        if y - 1 >= 0:
            W[y - 1, y] = 1
        else:
            W[d + y, y] = rho
    T = np.kron(V, np.eye(d))
    W = T @ W @ T.conj().T
    return np.kron(coin(th), np.eye(d)) @ W


def sx_phase(kx, d):
    """diagonal of S_x for wave number(s) kx: shape (2d,) or (2d, n)"""
    kx = np.asarray(kx, float)
    return np.concatenate([np.exp(-1j * np.multiply.outer(np.ones(d), kx)),
                           np.exp(1j * np.multiply.outer(np.ones(d), kx))], 0)


def floquet(M, kx, d=None):
    """eigen-decomposition of U(kx) = M . S_x(kx): quasi-energies w = -angle(lambda), sorted, and eigenvectors"""
    d = M.shape[0] // 2
    U = M * sx_phase(kx, d)[None, :]
    lam, v = np.linalg.eig(U)
    w = -np.angle(lam)
    o = np.argsort(w)
    # U is unitary and (for distinct w) its eigenvectors are orthogonal; re-orthonormalise degenerate blocks
    v = v[:, o]
    q, _ = np.linalg.qr(v)
    # keep the eigen-spaces: QR mixes only within-tolerance degenerate vectors if v is already near-unitary
    return w[o], v if np.allclose(v.conj().T @ v, np.eye(len(w)), atol=1e-8) else q


def masses(th, d, M=None):
    """half-gaps m_j at kx = 0: the positive quasi-energies below W_CUT (one per transverse strip mode), sorted"""
    M = static_part(th, d) if M is None else M
    w, _ = floquet(M, 0.0)
    return np.sort(w[(w > 0) & (w < W_CUT)])


def min_gap(th, d, M=None, kmax=0.3, n=121):
    """smallest half-gap over kx in [-kmax, kmax] (checks that the gap minimum is at kx = 0: no edge states)"""
    M = static_part(th, d) if M is None else M
    ks = np.linspace(-kmax, kmax, n)
    g = [np.min(np.abs(floquet(M, k)[0])) for k in ks]
    return float(min(g)), float(ks[int(np.argmin(g))])


def schwinger_strip(th, d, E, K=1.0, M=None):
    """the whole sea of the strip (states with w in (-W_CUT, 0) at kx = -K) evolved exactly in the uniform field
    A(t) = E t until kx = +K; returns (pairs = probability in (0, W_CUT), probability in |w| > W_CUT, ticks)"""
    M = static_part(th, d) if M is None else M
    d = M.shape[0] // 2
    w0, v0 = floquet(M, -K)
    sel = (w0 > -W_CUT) & (w0 < 0)
    psi = v0[:, sel]
    n = int(round(2 * K / E))
    kx = -K
    for t in range(n):
        psi = M @ (sx_phase(kx, d)[:, None] * psi)
        kx += E
    w1, v1 = floquet(M, kx)
    up = (w1 > 0) & (w1 < W_CUT)
    far = np.abs(w1) > W_CUT
    amp = v1.conj().T @ psi
    P0 = float((np.abs(amp[up]) ** 2).sum())
    Pf = float((np.abs(amp[far]) ** 2).sum())
    return P0, Pf, n

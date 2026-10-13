"""motore_f2_core.py - phase 2 of the "motore ottimizzato" study (docs/MOTORE_OTTIMIZZATO.md, section 14).

(1) ATOM in an octant. The learned coarse rule of chain A (lab/out/agg_zz_3d_A.pkl, levels B = 256 -> a0/16 'Ultra' and
    B = 512 -> a0/8 'Intel'), spin-blind part (the spin part is 1e-10), stencil SYMMETRIZED over the 48 cube
    operations (declared: the learned stencil carries a 0.1-0.3 % anisotropy at atomic k, learning noise + the xyz/zyx
    order), nucleus on a cell corner = origin. States with m = 0 symmetry are even in x, y and even or odd in z: they
    live in ONE OCTANT with real-to-real transforms (DCT-II, or DST-II along z for the odd sector). This makes the box
    +-12 a0 possible at a0/16 (octant 192^3 instead of 384^3).
    H = G(k)/tau + V (V = cell average of -alpha/r), states by LOBPCG; energies read from the RULE
    W = e^{-iV tau/2} e^{-iG} e^{-iV tau/2} as E = -arg<x|W|x>/tau (as in agg_zz_atom).
(2) LOCAL conservative level change: a lapped orthogonal transform (local cosine basis, Coifman-Meyer / Malvar folding
    with a sine bell) per block; keep the coefficients below the coarse Nyquist, synthesize on the coarse grid with the
    same bell. Orthonormal and local (a block talks only to its neighbours within the overlap).
(3) positive-energy descriptor for the 1-D learned levels (band projector of the coarse rule)."""
import os, sys, pickle, itertools
import numpy as np
import scipy.fft as sfft
import scipy.sparse.linalg as spl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zz_core as zz

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
MC2 = zz.MC2_EV
TH3 = A0 / 16 / 256                      # fine cell of chain A (lambda_C)


# =============================== (1) atom in an octant ===============================
def load_level3(B):
    st = pickle.load(open(os.path.join(OUT, "agg_zz_3d_A.pkl"), "rb"))
    L = [l for l in st["levels"] if l["B"] == B][0]
    offs = np.asarray(L["offs"])
    w = int(offs.max())
    n = 2 * w + 1
    G3 = (0.5 * (L["K"][0, 0] + L["K"][1, 1])).reshape(n, n, n)
    acc = 0
    for perm in itertools.permutations(range(3)):
        for fl in itertools.product([1, -1], repeat=3):
            acc = acc + np.transpose(G3, perm)[tuple(slice(None, None, f) for f in fl)]
    gs = (acc / 48).real
    asym = float(np.abs(G3 - gs).max() / np.abs(gs).max())
    return dict(B=B, w=w, g=gs, phi0=float(L["phi0"]), asym=asym, h=B * TH3, tau=B * TH3)


def cell_avg_inv_r_octant(n, h, nsub=16, near=4):
    x = (np.arange(n) + 0.5) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij", sparse=True)
    inv = 1.0 / np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    u = (np.arange(nsub) + 0.5) / nsub - 0.5
    U = np.stack(np.meshgrid(u, u, u, indexing="ij"), -1).reshape(-1, 3)
    for i in range(near):
        for j in range(near):
            for k in range(near):
                p = np.array([x[i], x[j], x[k]]) + U * h
                inv[i, j, k] = np.mean(1 / np.linalg.norm(p, axis=1))
    return inv, x


class OctAtom:
    def __init__(self, lev, n, sector="even"):
        self.n, self.lev, self.sector = n, lev, sector
        self.h, self.tau = lev["h"], lev["tau"]
        w = lev["w"]
        j = np.arange(-w, w + 1)
        kc = np.pi * np.arange(n) / n                 # DCT-II wave numbers (1/cell)
        ks = np.pi * (np.arange(n) + 1) / n           # DST-II wave numbers
        Cx = np.cos(np.outer(kc, j))
        Cz = np.cos(np.outer(ks if sector == "odd" else kc, j))
        self.G = np.einsum("abc,pa,qb,rc->pqr", lev["g"], Cx, Cx, Cz, optimize=True)
        inv, self.x = cell_avg_inv_r_octant(n, self.h)
        self.V = -zz.ALPHA * inv
        self.pre = 1.0 / (self.G / self.tau + 0.5 * zz.ALPHA ** 2)
        self.D = np.exp(-0.5j * self.V * self.tau)
        self.EG = np.exp(-1j * self.G)
        self.napply = 0

    # transforms: even in x, y always; z even (DCT) or odd (DST)
    def fwd(self, f):
        g = sfft.dctn(f, type=2, axes=(0, 1), norm="ortho")
        return (sfft.dst if self.sector == "odd" else sfft.dct)(g, type=2, axis=2, norm="ortho")

    def bwd(self, F):
        g = (sfft.idst if self.sector == "odd" else sfft.idct)(F, type=2, axis=2, norm="ortho")
        return sfft.idctn(g, type=2, axes=(0, 1), norm="ortho")

    def H(self, f):
        return self.bwd(self.G / self.tau * self.fwd(f)) + self.V * f

    def W(self, f):
        return self.D * self.bwd(self.EG * self.fwd(self.D * f))

    def _op(self, fun, X):
        X = np.asarray(X)
        one = X.ndim == 1
        X = X.reshape(self.n ** 3, -1)
        out = np.empty_like(X)
        for c in range(X.shape[1]):
            out[:, c] = fun(X[:, c].reshape((self.n,) * 3)).ravel()
            self.napply += 1
        return out[:, 0] if one else out

    def A(self, X):
        return self._op(self.H, X)

    def P(self, X):
        return self._op(lambda f: self.bwd(self.pre * self.fwd(f)), X)

    def grid(self):
        return np.meshgrid(self.x, self.x, self.x, indexing="ij")

    def guess(self, kind):
        X, Y, Z = self.grid()
        r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
        g = {"1s": np.exp(-r / A0), "2s": (1 - r / (2 * A0)) * np.exp(-r / (2 * A0)),
             "2p": Z * np.exp(-r / (2 * A0))}[kind]
        return (g / np.linalg.norm(g)).ravel()

    def energy(self, v):
        f = v.reshape((self.n,) * 3).astype(complex)
        Wf = self.W(f)
        lam = np.vdot(f, Wf) / np.vdot(f, f)
        resid = float(np.linalg.norm(Wf - lam * f) / np.linalg.norm(f))
        return float(-np.angle(lam) / self.tau), resid, complex(lam)


def solve(atom, kind, Y=None, x0=None, maxiter=60, tol=1e-10):
    N3 = atom.n ** 3
    Aop = spl.LinearOperator((N3, N3), matvec=atom.A, matmat=atom.A, dtype=float)
    Pop = spl.LinearOperator((N3, N3), matvec=atom.P, matmat=atom.P, dtype=float)
    x = (atom.guess(kind) if x0 is None else x0)[:, None]
    if Y is not None:
        x = x - Y @ (Y.T @ x)
    w, v, hist = spl.lobpcg(Aop, x, M=Pop, Y=Y, tol=tol, maxiter=maxiter, largest=False, retResidualNormsHistory=True)
    v = v[:, 0] / np.linalg.norm(v[:, 0])
    return v, float(np.max(hist[-1])) if len(hist) else float("nan"), len(hist)


def exact_kick_2p(q_a0, nr=4000, nmu=400):
    """|<2p0| e^{i q z} |1s>|^2 for hydrogen (Schroedinger, a0 units), by quadrature"""
    r = np.linspace(0, 60, nr + 1)[1:]
    dr = r[1] - r[0]
    mu, wmu = np.polynomial.legendre.leggauss(nmu)
    R, M = np.meshgrid(r, mu, indexing="ij")
    p1s = np.exp(-R) / np.sqrt(np.pi)
    p2p = R * np.exp(-R / 2) * M / np.sqrt(32 * np.pi)
    integrand = p2p * p1s * np.sin(q_a0 * R * M) * R ** 2 * 2 * np.pi
    amp = np.sum(integrand * wmu[None, :]) * dr
    return float(amp ** 2)


# =============================== (2) local cosine (lapped orthogonal) transform ===============================
BELL_ITER = 2


def _bell(e, k=None):
    """rho(s), rho(-s) on the half-overlap, rho^2(s) + rho^2(-s) = 1; beta iterated sine (k times): smoother bell,
    faster decay of the coefficients of smooth content"""
    k = BELL_ITER if k is None else k
    s = (np.arange(e) + 0.5) / e
    b = s
    for _ in range(k):
        b = np.sin(np.pi / 2 * b)
    return np.sin(np.pi / 4 * (1 + b)), np.cos(np.pi / 4 * (1 + b))


def lot_fwd(x, L, e, axis=0, sign=1):
    """lapped orthogonal transform along `axis` of a periodic array: blocks of L samples, folding half-width e at every
    block edge, DCT-IV in each block. Returns coefficients (same shape): block b, index m at position b*L + m."""
    x = np.moveaxis(np.asarray(x), axis, 0).copy()
    n = x.shape[0]
    rp, rm = _bell(e)
    sh = (e,) + (1,) * (x.ndim - 1)
    rp, rm = rp.reshape(sh), rm.reshape(sh)
    for a in range(0, n, L):
        ip = (a + np.arange(e)) % n                   # a + t
        im = (a - 1 - np.arange(e)) % n               # a - 1 - t
        u, v = x[ip].copy(), x[im].copy()
        x[ip] = rp * u + sign * rm * v
        x[im] = rp * v - sign * rm * u
    y = x.reshape((n // L, L) + x.shape[1:])
    y = sfft.dct(y, type=4, axis=1, norm="ortho")
    return np.moveaxis(y.reshape(x.shape), 0, axis)


def lot_bwd(c, L, e, axis=0, sign=1):
    c = np.moveaxis(np.asarray(c), axis, 0)
    n = c.shape[0]
    y = c.reshape((n // L, L) + c.shape[1:])
    x = sfft.idct(y, type=4, axis=1, norm="ortho").reshape(c.shape).copy()
    rp, rm = _bell(e)
    sh = (e,) + (1,) * (x.ndim - 1)
    rp, rm = rp.reshape(sh), rm.reshape(sh)
    for a in range(0, n, L):
        ip = (a + np.arange(e)) % n
        im = (a - 1 - np.arange(e)) % n
        u, v = x[ip].copy(), x[im].copy()
        x[ip] = rp * u - sign * rm * v
        x[im] = rp * v + sign * rm * u
    return np.moveaxis(x, 0, axis)


def lot_restrict(F, L, e, axes, sign=1):
    """fine -> coarse (factor 2 per axis): keep the low half of each block's coefficients, synthesize on the coarse
    grid with blocks of L/2 and overlap e/2 (the same bell, sampled at half resolution); amplitude / sqrt(2) per axis
    so that sum|F|^2 h^d is kept exactly for the kept coefficients."""
    C = F
    for ax in axes:
        c = lot_fwd(C, L, e, ax, sign)
        c = np.moveaxis(c, ax, 0)
        n = c.shape[0]
        cb = c.reshape((n // L, L) + c.shape[1:])[:, : L // 2]
        cc = cb.reshape((n // 2,) + c.shape[1:]) / np.sqrt(2) * np.sqrt(2)   # coefficient energy kept; samples halve
        C = lot_bwd(np.moveaxis(cc, 0, ax), L // 2, e // 2, ax, sign) / np.sqrt(2)
    return C


def lot_prolong(Fc, L, e, axes, sign=1):
    C = Fc
    for ax in axes:
        c = lot_fwd(C * np.sqrt(2), L // 2, e // 2, ax, sign)
        c = np.moveaxis(c, ax, 0)
        nc = c.shape[0]
        cb = c.reshape((nc // (L // 2), L // 2) + c.shape[1:])
        full = np.zeros((cb.shape[0], L) + c.shape[1:], dtype=c.dtype)
        full[:, : L // 2] = cb
        C = lot_bwd(np.moveaxis(full.reshape((nc * 2,) + c.shape[1:]), 0, ax), L, e, ax, sign)
    return C

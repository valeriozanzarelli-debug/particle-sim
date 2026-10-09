"""pre_core.py - shared tools for the pre-experiments (P1..P5) before the "self-built geometry" step.

1. LatticeWalk: the zig-zag rule of zz_core / geo_core on an n x n lattice (fixed axis order x then y,
   flip theta), with a WEIGHT w in [0,1] on every link. On a link the two components that want to cross it
   (a: at i moving +, b: at i+1 moving -) go through a 2x2 unitary beam splitter
        a'(i+1) = t a + i r b,   b'(i) = t b + i r a,    t = sqrt(w), r = sqrt(1 - w).
   w = 1: the plain rule (pure translation). w = 0: a mirror (the component turns around in place).
   Open box = periodic arrays with the wrap-around links set to w = 0 (edge reflection phase i instead of
   geo_core's -1: a boundary-condition choice, checked against geo_core on the torus where both agree).
2. Surfers on the lattice (no continuum!): per axis sub-step, an exactly equivariant LOCAL jump process.
   On each link: A = mass at i heading +, B = mass at i+1 heading -, after the step A' = mass arriving at
   i+1 through the link, B' = mass left at i. Any flows F_fwd (i->i+1), F_bwd (i+1->i) with
   F_fwd - F_bwd = A' - B, 0 <= F_fwd <= A, 0 <= F_bwd <= B keep the surfer histogram = |psi|^2 exactly.
   PRE-REGISTERED family: F_bwd = max(0, B - A') + lam * (min(B, B') - max(0, B - A')), lam in [0, 1]
       lam = 0  : Bell minimal net flow (surfers move only when the net current forces them)
       lam = 1  : maximal crossing (w = 1 -> fully ballistic, every surfer moves one site, like the walk)
       lam = w  : the link weight (default: walls never let surfers through, open links are ballistic)
   A surfer at i jumps + with prob F_fwd(i)/rho_i, - with prob F_bwd(i-1)/rho_i, else stays.
3. H-bar on cells, phase-winding count, MI helpers.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from zz_core import BASIS_HALF

SEED = 20261009
PI = np.pi


class LatticeWalk:
    def __init__(self, n, theta=0.5, periodic=False, w=None):
        self.n, self.theta, self.periodic = n, theta, periodic
        # W[j][x, y] = weight of the link from site (x,y) to (x,y)+e_j (wrapping)
        if w is None:
            w = [np.ones((n, n)), np.ones((n, n))]
        self.W = [np.array(w[0], float), np.array(w[1], float)]
        if not periodic:
            self.W[0][-1, :] = 0.0
            self.W[1][:, -1] = 0.0
        self.P = [BASIS_HALF[j][1] for j in (0, 1)]

    def set_w(self, w):
        self.W = [np.array(w[0], float), np.array(w[1], float)]
        if not self.periodic:
            self.W[0][-1, :] = 0.0
            self.W[1][:, -1] = 0.0

    # psi: (4, n, n, ...) ; comps 0,1 = L (moves along -sigma_j), 2,3 = R (along +sigma_j)
    def axis_step(self, psi, j, need_flux=False):
        P = self.P[j]
        Pd = P.conj().T
        L = np.tensordot(Pd, psi[0:2], axes=(1, 0))
        R = np.tensordot(Pd, psi[2:4], axes=(1, 0))
        # eigenvalue lam = (+1, -1): L block q moves -lam[q] -> q=1 moves +, q=0 moves -;
        #                            R block q moves +lam[q] -> q=0 moves +, q=1 moves -
        w = self.W[j]
        extra = (None,) * (psi.ndim - 3)
        t = np.sqrt(w)[(slice(None), slice(None)) + extra]
        r = np.sqrt(np.maximum(1 - w, 0))[(slice(None), slice(None)) + extra]
        ax = j
        flux = None
        if need_flux:
            A = np.abs(L[1])**2 + np.abs(R[0])**2                      # at i heading +
            B = np.roll(np.abs(L[0])**2 + np.abs(R[1])**2, -1, axis=ax)  # at i+1 heading - (aligned to link i)
        outs = []
        for blk, (pp, mm) in ((L, (1, 0)), (R, (0, 1))):
            a = blk[pp]
            bn = np.roll(blk[mm], -1, axis=ax)
            ap = t * a + 1j * r * bn          # arrives at i+1 (aligned to link i)
            bp = t * bn + 1j * r * a          # stays at i heading -
            nb = np.empty_like(blk)
            nb[pp] = np.roll(ap, 1, axis=ax)
            nb[mm] = bp
            outs.append(nb)
            if need_flux:
                if flux is None:
                    flux = [np.abs(ap)**2, np.abs(bp)**2]
                else:
                    flux[0] = flux[0] + np.abs(ap)**2
                    flux[1] = flux[1] + np.abs(bp)**2
        L2 = np.tensordot(P, outs[0], axes=(1, 0))
        R2 = np.tensordot(P, outs[1], axes=(1, 0))
        out = np.concatenate([L2, R2], axis=0)
        if need_flux:
            return out, (A, B, flux[0], flux[1])
        return out

    def coin(self, psi):
        ct, st = np.cos(self.theta), np.sin(self.theta)
        L, R = psi[0:2], psi[2:4]
        return np.concatenate([ct * L - 1j * st * R, ct * R - 1j * st * L], axis=0)

    def tick(self, psi):
        for j in (0, 1):
            psi = self.axis_step(psi, j)
        return self.coin(psi)

    def matrix(self):
        """one-tick unitary as sparse matrix, index a*Ns + (x*n + y) (same ordering as geo_core)"""
        n = self.n
        D = 4 * n * n
        E = np.eye(D, dtype=complex).reshape(4, n, n, D)
        U = self.tick(E).reshape(D, D)
        U[np.abs(U) < 1e-14] = 0
        return sp.csr_matrix(U)


def modes_state(walk, M, omega0=0.9, seed=1):
    """psi = sum of the M eigenmodes of the one-tick U closest to e^{-i omega0}, random phases (M = 1: one mode)"""
    U = walk.matrix().tocsc()
    D = U.shape[0]
    k = max(M, 2)
    vals, vecs = spl.eigs(U, k=k, sigma=np.exp(-1j * omega0), which="LM", tol=1e-10)
    om = -np.angle(vals)
    o = np.argsort(np.abs(om - omega0))[:M]
    vecs = vecs[:, o]
    vecs /= np.linalg.norm(vecs, axis=0)
    rng = np.random.default_rng(seed)
    ph = np.exp(1j * rng.uniform(0, 2 * PI, M)) if M > 1 else np.ones(1)
    psi = (vecs * ph[None, :]).sum(1)
    psi /= np.linalg.norm(psi)
    return psi.reshape(4, walk.n, walk.n), om[o]


def flows(A, B, Ap, Bp, lam):
    """equivariant link flows (see module doc). A, B, Ap (= mass arriving at i+1), Bp (= mass left at i)"""
    lo = np.maximum(0.0, B - Ap)
    Fb = lo + lam * (np.minimum(B, Bp) - lo)
    Ff = Ap - B + Fb
    return np.clip(Ff, 0, None), np.clip(Fb, 0, None)


def surfer_axis(walk, psi, X, j, rng, lam, count=None):
    """advance psi by one axis sub-step and move surfers X (N,2 int) consistently. lam: scalar, 'w', or array.
    count: optional array (n,n) accumulating crossings of links along axis j."""
    n = walk.n
    psi2, (A, B, Ap, Bp) = walk.axis_step(psi, j, need_flux=True)
    lamv = walk.W[j] if (isinstance(lam, str) and lam == "w") else lam
    Ff, Fb = flows(A, B, Ap, Bp, lamv)
    rho = np.sum(np.abs(psi)**2, axis=0)
    rho = np.maximum(rho, 1e-300)
    pf = Ff / rho                                  # at site i, jump + along j
    pb = np.roll(Fb, 1, axis=j) / rho              # at site i, jump - along j (link i-1)
    x, y = X[:, 0], X[:, 1]
    u = rng.random(len(X))
    p1 = pf[x, y]; p2 = pb[x, y]
    fwd = u < p1
    bwd = (~fwd) & (u < p1 + p2)
    if count is not None:
        if j == 0:
            np.add.at(count, (x[fwd], y[fwd]), 1)
            np.add.at(count, ((x[bwd] - 1) % n, y[bwd]), 1)
        else:
            np.add.at(count, (x[fwd], y[fwd]), 1)
            np.add.at(count, (x[bwd], (y[bwd] - 1) % n), 1)
    step = fwd.astype(int) - bwd.astype(int)
    X = X.copy()
    X[:, j] = (X[:, j] + step) % n
    return psi2, X


def surfer_tick(walk, psi, X, rng, lam="w", counts=None):
    for j in (0, 1):
        psi, X = surfer_axis(walk, psi, X, j, rng, lam, None if counts is None else counts[j])
    return walk.coin(psi), X


def sample_sites(rng, N, prob):
    p = prob.ravel() / prob.sum()
    idx = rng.choice(len(p), size=N, p=p)
    n = prob.shape[0]
    return np.stack([idx // n, idx % n], axis=1)


def hbar(X, rho, cell=4, sheets=False):
    """coarse-grained H = sum p ln(p/q) on cell x cell blocks (optionally split by sheet (x+y) mod 2)"""
    n = rho.shape[0]
    K = n // cell
    cx, cy = X[:, 0] // cell, X[:, 1] // cell
    ok = (cx < K) & (cy < K)
    if sheets:
        s = (X[:, 0] + X[:, 1]) % 2
        cnt = np.zeros((K, K, 2))
        np.add.at(cnt, (cx[ok], cy[ok], s[ok]), 1)
        xx, yy = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        ref = np.zeros((K, K, 2))
        m = (xx // cell < K) & (yy // cell < K)
        np.add.at(ref, ((xx // cell)[m], (yy // cell)[m], ((xx + yy) % 2)[m]), rho[m])
    else:
        cnt = np.zeros((K, K))
        np.add.at(cnt, (cx[ok], cy[ok]), 1)
        ref = rho[:K * cell, :K * cell].reshape(K, cell, K, cell).sum((1, 3))
    p = cnt / cnt.sum(); q = ref / ref.sum()
    m = p > 0
    return float((p[m] * np.log(p[m] / np.maximum(q[m], 1e-300))).sum()), int(cnt.size)


def windings(psi):
    """phase-winding plaquettes of each spinor component, counted on each sheet separately
    (sheet neighbours of (x,y): (x+1,y+1), (x+1,y-1)); plaquette (x,y),(x+1,y+1),(x+2,y),(x+1,y-1).
    Returns the mean count over the 4 components (both sheets summed). Interior only."""
    def w(a, b):
        d = np.angle(b) - np.angle(a)
        return (d + PI) % (2 * PI) - PI
    tot = 0
    for c in psi:
        p0 = c[:-2, 1:-1]; p1 = c[1:-1, 2:]; p2 = c[2:, 1:-1]; p3 = c[1:-1, :-2]
        s = w(p0, p1) + w(p1, p2) + w(p2, p3) + w(p3, p0)
        tot += int((np.abs(s) > PI).sum())
    return tot / 4.0


def fit_tau(t, h, floor):
    """tau = first time at which (H - floor) falls below (H0 - floor)/e (linear interp); inf if never"""
    y = np.asarray(h) - floor
    target = y[0] / np.e
    idx = np.where(y <= target)[0]
    if len(idx) == 0:
        return float("inf")
    k = idx[0]
    if k == 0:
        return 0.0
    f = (y[k - 1] - target) / max(y[k - 1] - y[k], 1e-12)
    return float(t[k - 1] + f * (t[k] - t[k - 1]))

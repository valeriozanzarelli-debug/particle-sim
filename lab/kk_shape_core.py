"""kk_shape_core.py - the zig-zag rule with a hidden SPACE THAT HAS A SHAPE (Valerio's "forma nascosta").

Same rule as lab/zz_core.py / lab/kk_core.py: every ball translates at c along each axis in turn (split on
the eigenvectors of sigma_j; L moves opposite to R; order 0,1,2 on even ticks, 2,1,0 on odd ticks) and
then rotates L<->R by the flip angle theta = m a (mass). Nothing new is added to the rule:

  * axis 0 = the hidden CIRCLE of N0 sites (charge n, as in kk_core). For a ball at rest in real space we
    only need the hidden axes; axis 0 is handled exactly in Fourier space (charge n -> phase e^{-i k0 s}).
  * axes 1, 2 = a hidden PLANE (grid n1 x n2). The SHAPE is drawn in it with the flip angle itself:
    theta = theta_in inside the shape, theta = pi/2 outside ("mass wall": an infinitely heavy outside, the
    MIT-bag / Berry-Mondragon 'neutrino billiard' boundary). With theta = pi/2 the flip swaps L and R
    completely, so a ball cannot go more than ~1 site into the wall: the shape confines by itself.
  * 1-D hidden LOOPS (a closed curve of N sites, parametrized by arclength) use the 1-D zig-zag rule
    (R moves +1, L moves -1, then flip theta(s)).

Lattice units a = dt = c = hbar = 1. Seeded (SEED below). Uses <= 2 threads.

Doublers (lattice artefact, handled exactly): with two hidden plane axes the sign maps (-1)^{w1}, (-1)^{w2}
commute with the two-tick operator U2 and U2 maps each parity class (w1 mod 2, w2 mod 2) to itself, so the
spectrum of U2 is 4 identical copies (the physical one + 3 fermion doublers). We keep ONE copy by
restricting U2 to the (even, even) sites. At charge n = 0 the operator M = tau_x (x) sigma_x (L/R x spin)
commutes with the rule and splits each copy in two blocks; we keep the M = +1 block.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np
import scipy.sparse as sp
from scipy.special import jv
from scipy.optimize import brentq
from zz_core import BASIS_HALF, order_of_tick

SEED = 20261009
WALL = np.pi / 2


# ======================================== shapes ========================================
def shape_mask(kind, area, par=None, margin=4):
    """boolean mask (n1, n2) of a hidden 2-D shape of given AREA (in sites) centred in an even-sized box.
    kind: 'disk' | 'ellipse' (par = aspect a/b) | 'triangle' (equilateral) | 'stadium' (par = a/r, straight
    half-length over cap radius; Bunimovich 'standard' a = r) | 'qstadium' / 'qdisk' / 'qellipse' (quarter,
    desymmetrized: the shape restricted to x > 0, y > 0, walls on the cuts).
    returns mask, centre (c1, c2) and the shape parameters actually used."""
    q = kind.startswith("q")
    base = kind[1:] if q else kind
    A = area * (4 if q else 1)                       # quarter shapes: area of the quarter = area
    if base == "disk":
        R = np.sqrt(A / np.pi); h1 = h2 = R
        f = lambda x, y: x**2 + y**2 < R**2
        info = dict(R=R)
    elif base == "ellipse":
        asp = par; a = np.sqrt(A * asp / np.pi); b = a / asp; h1, h2 = a, b
        f = lambda x, y: (x / a)**2 + (y / b)**2 < 1
        info = dict(a=a, b=b)
    elif base == "stadium":
        al = 1.0 if par is None else par
        r = np.sqrt(A / (np.pi + 4 * al)); a = al * r; h1, h2 = a + r, r
        f = lambda x, y: ((np.abs(x) <= a) & (np.abs(y) < r)) | ((np.abs(x) - a)**2 + y**2 < r**2)
        info = dict(a=a, r=r)
    elif base == "triangle":
        s = np.sqrt(4 * A / np.sqrt(3)); hgt = s * np.sqrt(3) / 2
        # centroid at origin, base at y = -hgt/3
        y0 = -hgt / 3
        f = lambda x, y: (y > y0) & (y < y0 + hgt - np.sqrt(3) * np.abs(x))
        h1, h2 = s / 2, 2 * hgt / 3
        info = dict(side=s)
    else:
        raise ValueError(kind)
    if q:
        n1 = 2 * int(np.ceil((h1 + 2 * margin) / 2)) + 2; n2 = 2 * int(np.ceil((h2 + 2 * margin) / 2)) + 2
        c1, c2 = margin - 0.5, margin - 0.5          # cut lines between sites
        X, Y = np.meshgrid(np.arange(n1) - c1, np.arange(n2) - c2, indexing="ij")
        mask = f(X, Y) & (X > 0) & (Y > 0)
    else:
        n1 = 2 * int(np.ceil(h1 + margin)); n2 = 2 * int(np.ceil(h2 + margin))
        c1, c2 = (n1 - 1) / 2, (n2 - 1) / 2
        X, Y = np.meshgrid(np.arange(n1) - c1, np.arange(n2) - c2, indexing="ij")
        mask = f(X, Y)
    info.update(kind=kind, area_target=area, sites=int(mask.sum()), n1=n1, n2=n2, c1=c1, c2=c2)
    return mask, info


def soft_wall(mask, w):
    """flip-angle map with a SMOOTH mass wall: theta rises from 0 at the shape boundary to pi/2 over w sites
    (smoothstep of the distance outside the shape). w = 0 -> sharp wall."""
    from scipy.ndimage import distance_transform_edt
    if w <= 0:
        return np.where(mask, 0.0, WALL)
    d = distance_transform_edt(~mask)
    u = np.clip(d / w, 0, 1)
    return WALL * u * u * (3 - 2 * u)


# ============================ the rule as a sparse operator ============================
def _shift(n1, n2, g, s):
    """periodic roll of a (n1, n2) field by s along grid axis g, as a sparse matrix (np.roll semantics)"""
    idx = np.arange(n1 * n2).reshape(n1, n2)
    src = np.roll(idx, s, axis=g).ravel()
    N = n1 * n2
    return sp.csr_matrix((np.ones(N), (np.arange(N), src)), shape=(N, N))


def _axis_block(j, h, n1, n2, k0):
    """2x2 block (list of lists of sparse) of the translation along walk axis j for one handedness h"""
    lam, P = BASIS_HALF[j]
    N = n1 * n2
    sh = []
    for q in range(2):
        s = int(h * lam[q])
        if j == 0:
            sh.append(np.exp(-1j * k0 * s) * sp.identity(N, format="csr", dtype=complex))
        else:
            sh.append(_shift(n1, n2, j - 1, s).astype(complex))
    return [[sum(P[a, q] * np.conj(P[b, q]) * sh[q] for q in range(2)) for b in range(2)] for a in range(2)]


def tick_operator(theta, t, n=0, N0=None, V=None):
    """one tick of the zig-zag rule on the hidden plane (axes 1,2) x hidden circle in Fourier mode n (axis 0).
    theta: (n1, n2) flip-angle map (pi/2 = wall). Vector layout: comp-major, comp in (L0, L1, R0, R1)."""
    n1, n2 = theta.shape
    N = n1 * n2
    k0 = 0.0 if n == 0 else 2 * np.pi * n / N0
    I = sp.identity(N, format="csr", dtype=complex)
    TL = [[I, None], [None, I]]; TR = [[I, None], [None, I]]
    TL = sp.bmat(TL, format="csr"); TR = sp.bmat(TR, format="csr")
    for j in order_of_tick(t):
        TL = sp.bmat(_axis_block(j, -1, n1, n2, k0), format="csr") @ TL
        TR = sp.bmat(_axis_block(j, +1, n1, n2, k0), format="csr") @ TR
    T = sp.block_diag([TL, TR], format="csr")
    c = np.cos(theta).ravel(); s = np.sin(theta).ravel()
    C = sp.diags(np.tile(c, 2)); S = sp.diags(-1j * np.tile(s, 2))
    F = sp.bmat([[C, S], [S, C]], format="csr")
    U = F @ T
    if V is not None:
        U = sp.diags(np.tile(np.exp(-1j * V.ravel()), 4)) @ U
    return U.tocsr()


def two_tick(theta, n=0, N0=None, V=None):
    return (tick_operator(theta, 1, n, N0, V) @ tick_operator(theta, 0, n, N0, V)).tocsr()


def ee_index(n1, n2):
    """indices of the 4 components on the (even, even) sublattice (comp-major layout)"""
    I1, I2 = np.meshgrid(np.arange(n1), np.arange(n2), indexing="ij")
    sites = np.flatnonzero(((I1 % 2) == 0).ravel() & ((I2 % 2) == 0).ravel())
    N = n1 * n2
    return np.concatenate([c * N + sites for c in range(4)]), sites


def M_basis(nsites):
    """sparse isometry onto the M = tau_x (x) sigma_x = +1 block (2 of the 4 components per site)"""
    Mloc = np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, 1], [1, 0]])).astype(complex)
    w, v = np.linalg.eigh(Mloc)
    B = v[:, w > 0]                                   # 4 x 2
    return sp.kron(sp.csr_matrix(B), sp.identity(nsites, format="csr")).tocsr()


def reduced_U2(theta, n=0, N0=None, block=True):
    """U2 restricted to one doubler copy (ee sublattice) and, at n = 0 and block=True, to the M = +1 block"""
    n1, n2 = theta.shape
    U2 = two_tick(theta, n, N0)
    idx, sites = ee_index(n1, n2)
    Ue = U2[idx][:, idx]
    if n == 0 and block:
        P = M_basis(len(sites))
        Ue = (P.conj().T @ Ue @ P).tocsr()
        Ue.eliminate_zeros()
    return Ue, sites


def rest_masses(theta, n=0, N0=None, Emax=0.6, block=True, vectors=False):
    """rest energies (positive quasi-energies E = -arg(lambda)/2 of the two-tick operator) below Emax, sorted.
    Dense diagonalisation of the reduced operator. The wall flat band sits at E = pi/2 (excluded)."""
    Ue, sites = reduced_U2(theta, n, N0, block)
    lam, V = np.linalg.eig(Ue.toarray())
    E = -np.angle(lam) / 2
    keep = (E > 1e-9) & (E < Emax)
    o = np.argsort(E[keep])
    if vectors:
        return E[keep][o], V[:, keep][:, o], sites
    return E[keep][o]


def rest_masses_sparse(theta, k=24, n=0, N0=None, block=True, vectors=False):
    """the k lowest positive rest energies via shift-invert around lambda = 1 (E = 0)"""
    from scipy.sparse.linalg import eigs
    Ue, sites = reduced_U2(theta, n, N0, block)
    lam, V = eigs(Ue.tocsc(), k=2 * k + 4, sigma=1.0, which="LM")
    E = -np.angle(lam) / 2
    keep = E > 1e-9
    o = np.argsort(E[keep])
    if vectors:
        return E[keep][o][:k], V[:, keep][:, o][:, :k], sites
    return E[keep][o][:k]


def pair_up(E):
    """RESIDUAL DOUBLER: with the palindromic two-tick order the walk has a second Dirac cone at
    k = (pi/2, pi/2); the mass wall mixes it with the physical one, so every level comes as a near-degenerate
    pair (bonding/antibonding). Returns (midpoints, splittings, worst splitting / next gap)."""
    E = np.sort(E)
    m = len(E) // 2
    a, b = E[0:2 * m:2], E[1:2 * m:2]
    mid, spl = (a + b) / 2, b - a
    gaps = a[1:] - b[:-1]
    worst = float(np.max(spl[:-1] / np.maximum(gaps, 1e-300))) if m > 1 else np.nan
    return mid, spl, worst


def mode_density(vec, sites, n1, n2, block=True):
    """|psi|^2 on the hidden plane of a reduced eigenvector (shown on the ee sublattice, upsampled x2)"""
    ns = len(sites)
    if block:
        vec = (M_basis(ns) @ vec)
    rho = np.sum(np.abs(vec.reshape(4, ns))**2, axis=0)
    out = np.zeros(n1 * n2); out[sites] = rho
    out = out.reshape(n1, n2)
    # fill the 3 other sublattices by nearest ee value (display only)
    out[1::2, :] = out[0:n1 - 1:2, :][: out[1::2, :].shape[0]]
    out[:, 1::2] = out[:, 0:n2 - 1:2][:, : out[:, 1::2].shape[1]]
    return out


# ============================ direct time stepping (no matrices) ============================
class ShapeWalk:
    """the zig-zag rule on a real grid of hidden axes: axis 0 = circle (N0 sites, periodic, may be 1),
    axes 1, 2 = hidden plane with flip map theta (n1, n2). psi shape (4, N0, n1, n2)."""

    def __init__(self, theta, N0=1, V=None):
        self.theta = theta[None]                       # broadcast over axis 0
        self.N0 = N0
        self.phase = None if V is None else np.exp(-1j * V)
        self.t = 0

    @staticmethod
    def _axis(blk, j, h):
        lam, P = BASIS_HALF[j]
        c = np.tensordot(P.conj().T, blk, axes=(1, 0))
        for q in range(2):
            c[q] = np.roll(c[q], int(h * lam[q]), axis=j)
        return np.tensordot(P, c, axes=(1, 0))

    def tick(self, psi):
        L, R = psi[0:2], psi[2:4]
        for j in order_of_tick(self.t):
            L = self._axis(L, j, -1); R = self._axis(R, j, +1)
        ct, st = np.cos(self.theta), np.sin(self.theta)
        out = np.concatenate([ct * L - 1j * st * R, ct * R - 1j * st * L])
        if self.phase is not None:
            out *= self.phase
        self.t += 1
        return out


def circle_charge(psi):
    """<n> on the hidden circle (axis 1 of psi)"""
    N0 = psi.shape[1]
    F = np.fft.fft(psi, axis=1)
    w = np.sum(np.abs(F)**2, axis=(0, 2, 3)); w /= w.sum()
    nv = np.rint(np.fft.fftfreq(N0) * N0)
    return float(np.sum(nv * w))


def Jz(psi, c1, c2):
    """<J_z> = <L_z + sigma_x/2> in the hidden plane (rotation generator for axes 1 -> sigma_y, 2 -> sigma_z:
    [sigma_x/2, sigma_y] = i sigma_z, same for L and R). L_z = w1 p2 - w2 p1 with spectral derivatives
    (the field vanishes near the box edge, inside the walls)."""
    _, N0, n1, n2 = psi.shape
    k1 = 2 * np.pi * np.fft.fftfreq(n1); k2 = 2 * np.pi * np.fft.fftfreq(n2)
    F = np.fft.fft2(psi, axes=(2, 3))
    d1 = np.fft.ifft2(1j * k1[None, None, :, None] * F, axes=(2, 3))
    d2 = np.fft.ifft2(1j * k2[None, None, None, :] * F, axes=(2, 3))
    X1 = (np.arange(n1) - c1)[None, None, :, None]; X2 = (np.arange(n2) - c2)[None, None, None, :]
    Lz = -1j * (X1 * d2 - X2 * d1)
    sx = np.empty_like(psi)
    sx[0], sx[1], sx[2], sx[3] = psi[1], psi[0], psi[3], psi[2]
    num = np.sum(np.conj(psi) * Lz) + 0.5 * np.sum(np.conj(psi) * sx)
    return float(np.real(num / np.sum(np.abs(psi)**2)))


# ================================ hidden 1-D loops ================================
def loop_curve(kind, par=None, M=4096):
    """closed curve with total length 1, sampled densely. returns arclength s (M,), curvature kappa(s) (M,),
    points xy (M, 2). kinds: 'circle', 'ellipse' (par = aspect), 'trefoil' (r = 1 + par cos 3phi),
    'egg' (r = 1 + par[0] cos 2phi + par[1] cos 3phi), 'stadium' (par = a/r)"""
    if kind == "stadium":
        al = 1.0 if par is None else par
        r = 1.0; a = al
        P = 2 * np.pi * r + 4 * a
        s = np.arange(M) / M * P
        xy = np.zeros((M, 2)); kap = np.zeros(M)
        for i, si in enumerate(s):
            if si < 2 * a:                     # bottom straight, left to right
                xy[i] = (-a + si, -r); kap[i] = 0
            elif si < 2 * a + np.pi * r:      # right cap
                ph = -np.pi / 2 + (si - 2 * a) / r; xy[i] = (a + r * np.cos(ph), r * np.sin(ph)); kap[i] = 1 / r
            elif si < 4 * a + np.pi * r:      # top straight
                xy[i] = (a - (si - 2 * a - np.pi * r), r); kap[i] = 0
            else:
                ph = np.pi / 2 + (si - 4 * a - np.pi * r) / r; xy[i] = (-a + r * np.cos(ph), r * np.sin(ph)); kap[i] = 1 / r
        return s / P, kap * P, xy / P
    phi = np.linspace(0, 2 * np.pi, 20 * M, endpoint=False)
    if kind == "circle":
        x, y = np.cos(phi), np.sin(phi)
    elif kind == "ellipse":
        x, y = np.sqrt(par) * np.cos(phi), np.sin(phi) / np.sqrt(par)
    elif kind in ("trefoil", "egg"):
        rr = 1 + (par * np.cos(3 * phi) if kind == "trefoil" else par[0] * np.cos(2 * phi) + par[1] * np.cos(3 * phi))
        x, y = rr * np.cos(phi), rr * np.sin(phi)
    else:
        raise ValueError(kind)
    dphi = phi[1] - phi[0]
    dd = lambda f: (np.roll(f, -1) - np.roll(f, 1)) / (2 * dphi)       # periodic derivative
    xp, yp = dd(x), dd(y)
    xpp, ypp = dd(xp), dd(yp)
    kap = (xp * ypp - yp * xpp) / (xp**2 + yp**2)**1.5
    ds = np.sqrt(xp**2 + yp**2) * dphi
    s = np.concatenate([[0], np.cumsum(ds)[:-1]]); P = s[-1] + ds[-1]
    su = np.arange(M) / M * P
    kap_u = np.interp(su, s, kap, period=P)
    xy = np.stack([np.interp(su, s, x, period=P), np.interp(su, s, y, period=P)], 1)
    return su / P, kap_u * P, xy / P            # unit length; kappa in units of 1/length (mean 2 pi)


def loop_sample(kind, par, N):
    """curvature (units: per unit loop length, mean 2pi for convex curves) at N equal-arclength sites"""
    s, kap, xy = loop_curve(kind, par, M=max(4096, 8 * N))
    idx = (np.arange(N) * len(s)) // N
    return kap[idx], xy[idx]


def loop_U(theta, V=None):
    """one tick of the 1-D zig-zag rule on a loop of N sites (2 comps R, L): R -> +1, L -> -1, then flip"""
    N = len(theta)
    Sp = _shift(N, 1, 0, +1).astype(complex); Sm = _shift(N, 1, 0, -1).astype(complex)
    T = sp.block_diag([Sp, Sm])
    C = sp.diags(np.cos(theta)); S = sp.diags(-1j * np.sin(theta))
    U = sp.bmat([[C, S], [S, C]]) @ T
    if V is not None:
        U = sp.diags(np.tile(np.exp(-1j * V), 2)) @ U
    return U.tocsr()


def loop_masses(theta, V=None, Emax=None):
    """positive quasi-energies (rest masses) of the loop walk, sorted (doublers sit near E = pi, cut by Emax)"""
    lam = np.linalg.eigvals(loop_U(theta, V).toarray())
    E = -np.angle(lam)
    Emax = np.pi / 2 if Emax is None else Emax
    E = np.sort(E[(E >= -1e-12) & (E < Emax)])
    return np.abs(E)


# ================================ continuum references ================================
def dirac_disk_levels(R, kmax, sign=(+1, -1)):
    """2-D massless Dirac in a disk with infinite-mass wall (Berry & Mondragon 1987):
    J_l(kR) = s J_{l+1}(kR), l = 0, 1, 2, ..., s = +-1 (the two signs = the two 2-component blocks / both
    orientations). returns sorted k values below kmax with (l, s) labels."""
    out = []
    for s in sign:
        for l in range(0, 200):
            f = lambda x: jv(l, x) - s * jv(l + 1, x)
            xs = np.linspace(1e-6, kmax * R, 4000)
            fv = f(xs)
            roots = [brentq(f, xs[i], xs[i + 1]) for i in range(len(xs) - 1) if fv[i] * fv[i + 1] < 0]
            if not roots:
                if l * 1.0 > kmax * R + 5:
                    break
                continue
            out += [(x / R, l, s) for x in roots]
    out.sort()
    return out


def laplacian(mask):
    """5-point Dirichlet Laplacian (-Delta) on the sites of mask, h = 1"""
    n1, n2 = mask.shape
    idx = -np.ones(mask.shape, int); pts = np.argwhere(mask); idx[mask] = np.arange(len(pts))
    rows, cols, vals = [], [], []
    for p, (i, j) in enumerate(pts):
        rows.append(p); cols.append(p); vals.append(4.0)
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            a, b = i + di, j + dj
            if 0 <= a < n1 and 0 <= b < n2 and mask[a, b]:
                rows.append(p); cols.append(idx[a, b]); vals.append(-1.0)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(pts), len(pts))), pts


# ================================ statistics, Koide ================================
def ratio_stat(E):
    """mean of r~_n = min(s_n, s_{n+1}) / max(...) (Oganesyan-Huse; Atas et al. PRL 110, 084101 (2013)).
    Poisson 0.3863, GOE 0.5307, GUE 0.5996, GSE 0.6744. No unfolding needed."""
    s = np.diff(np.sort(E))
    s = s[:-1], s[1:]
    r = np.minimum(*s) / np.maximum(np.maximum(*s), 1e-300)
    return float(np.mean(r)), r


def unfold(E, deg=6):
    """unfold levels with a smooth polynomial fit of the staircase N(E) -> spacings of mean 1"""
    E = np.sort(E)
    Nst = np.arange(1, len(E) + 1)
    c = np.polyfit(E, Nst, deg)
    x = np.polyval(c, E)
    return np.diff(x)


def koide(m):
    m = np.asarray(m, float)
    return float(np.sum(m) / np.sum(np.sqrt(m))**2)


LEPTON = dict(me=0.51099895, mmu=105.6583755, mtau=1776.86)
R2_OBS = LEPTON["mmu"] / LEPTON["me"]          # 206.768
R3_OBS = LEPTON["mtau"] / LEPTON["me"]         # 3477.2
Q_OBS = koide([LEPTON["me"], LEPTON["mmu"], LEPTON["mtau"]])

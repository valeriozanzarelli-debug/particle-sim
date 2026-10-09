"""eye_core.py - L1 "the eye": the repaired shape reader (big_core, H1 only) plus CAVITIES (H2) and CONE POINTS
(local ball-volume deficit), all read blind from the mutual-information matrix I of the balls.

Pipeline (input = I only, never the adjacency or coordinates):
  largest informational piece -> d = -log(I/I0) -> Isomap geodesic Gm (geo_core, factor 1.5) -> diameter D
  LOOPS    : big_core H1 reader unchanged (220 max-min landmarks, full filtration, rel = pers/D > 0.4)
  CAVITIES : ripser H2 on <= 150 max-min landmarks, full filtration, rel2 = pers/D > THR2 (THR2 from the prereg /
             calibration on nulls, see eye_prereg.py)
  CONES    : h = median nearest-neighbour geodesic distance; for c in SCALES = (1.5, 2.5): V_x(c) = #{y : Gm(x,y) <= c h}
             (incl. x), rho_x(c) = V_x(c) / median_y V_y(c); rho_bar = mean over the two scales.
             A CONE POINT = site with rho_x(c) < RHO_THR at BOTH scales, nearest-neighbour distance <= NN_GATE h (not an
             isolated outlier), and the smallest rho_bar within its own 2.5 h ball (ties: lowest index).
             Angle estimate (2-D) = 2 pi rho_bar; solid-angle fraction (any D) = rho_bar.

Drawn shapes for the exams: the balls are one fermion mode per site hopping on the drawn graph (H = -T, T = hopping
matrix), in the thermal Gaussian state C = 1/(1 + exp(beta (H - mu))). For a quotient (orbifold) the balls live in the
symmetric (untwisted) sector: H_q = Q^T H Q with Q the normalised orbit indicator. Only I = MI(C) is handed to the eye.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np
import scipy.sparse as sp
import geo_core as G
import big_core as B
import orb_core as O
from pre_p2 import mi_one_mode, random_regular

BETA, MU = 1.0, 0.0           # thermal state of the balls (pre-registered)
NLAND2 = 100                  # landmarks for H2 (150 timed at 15-80 s per null in dev: too slow for 140 nulls)
SCALES = (1.5, 2.5)           # cone ball radii in units of h (half-integers: avoid ties with lattice hop multiples)
RHO_THR = 0.7                 # cone threshold on the volume ratio, at BOTH scales
NN_GATE = 1.2                 # a cone point must be part of the fabric: nearest-neighbour distance <= 1.2 h
THR2_DEFAULT = 0.4            # cavity threshold before calibration


# ------------------------------------------------------------------ balls -> MI
def thermal_C(H, beta=BETA, mu=MU):
    w, V = np.linalg.eigh(H)
    f = 0.5 * (1 - np.tanh(0.5 * beta * (w - mu)))
    return (V * f) @ V.conj().T


def mi_of_hopping(T, beta=BETA, mu=MU):
    H = -np.asarray(T.toarray() if sp.issparse(T) else T, float)
    return mi_one_mode(thermal_C(H, beta, mu))


# ------------------------------------------------------------------ drawn graphs (hopping matrices)
def hex_torus(L):
    A = sum(O.hex_shift(L, e) for e in O.HEX_E)
    return (A + A.T).toarray()


def hex_xy(L):
    i, j = np.divmod(np.arange(L * L), L)
    return np.stack([i + 0.5 * j, np.sqrt(3) / 2 * j], 1)


def orbit_matrix(perms, N):
    """perms = list of index arrays g(x) (group elements incl. identity). returns Q (N x norb) normalised, orbit id"""
    oid = -np.ones(N, int); orbits = []
    for x in range(N):
        if oid[x] >= 0:
            continue
        orb = sorted({int(p[x]) for p in perms})
        for y in orb:
            oid[y] = len(orbits)
        orbits.append(orb)
    Q = np.zeros((N, len(orbits)))
    for k, orb in enumerate(orbits):
        Q[orb, k] = 1 / np.sqrt(len(orb))
    return Q, oid, orbits


def t2z3(L):
    """T^2/Z3 from the hex L x L torus (L % 3 == 0): returns hopping on orbits, orbit sizes, fixed-orbit ids"""
    A = hex_torus(L)
    P = O.hex_rot_perm(L).tocoo()
    g = np.empty(L * L, int); g[P.col] = P.row
    perms = [np.arange(L * L), g, g[g]]
    Q, oid, orbits = orbit_matrix(perms, L * L)
    fixed = [k for k, o in enumerate(orbits) if len(o) == 1]
    return Q.T @ A @ Q, np.array([len(o) for o in orbits]), fixed, orbits


def hypercubic(n, d):
    N = n**d
    idx = np.arange(N).reshape((n,) * d)
    A = np.zeros((N, N))
    for ax in range(d):
        nb = np.roll(idx, -1, axis=ax).ravel()
        A[np.arange(N), nb] += 1; A[nb, np.arange(N)] += 1
    return A, idx


def t4z2(n):
    A, idx = hypercubic(n, 4)
    coords = np.array(np.unravel_index(np.arange(n**4), (n,) * 4)).T
    g = np.ravel_multi_index(tuple(((-coords) % n).T), (n,) * 4)
    Q, oid, orbits = orbit_matrix([np.arange(n**4), g], n**4)
    fixed = [k for k, o in enumerate(orbits) if len(o) == 1]
    return Q.T @ A @ Q, np.array([len(o) for o in orbits]), fixed, orbits


def icosphere(levels):
    t = (1 + 5**0.5) / 2
    V = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
         (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    F = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6),
         (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10),
         (8, 6, 7), (9, 8, 1)]
    V = [np.array(v, float) / np.linalg.norm(v) for v in V]
    for _ in range(levels):
        cache = {}
        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in cache:
                m = V[a] + V[b]; V.append(m / np.linalg.norm(m)); cache[k] = len(V) - 1
            return cache[k]
        NF = []
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            NF += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        F = NF
    X = np.array(V); N = len(X)
    A = np.zeros((N, N))
    for a, b, c in F:
        for u, v in ((a, b), (b, c), (c, a)):
            A[u, v] = A[v, u] = 1
    return A, X, np.array(F)


# ------------------------------------------------------------------ the eye
def _h2_rel(Gm, diam, nland=NLAND2):
    from ripser import ripser
    L = G.landmarks(Gm, nland)
    dg = ripser(Gm[np.ix_(L, L)], distance_matrix=True, maxdim=2)["dgms"]
    D2 = dg[2]
    if len(D2) == 0:
        return np.zeros(0), D2
    assert np.all(np.isfinite(D2[:, 1]))
    pers = D2[:, 1] - D2[:, 0]
    o = np.argsort(pers)[::-1]
    return pers[o] / diam, D2[o]


def cones(Gm, scales=SCALES, rho_thr=RHO_THR, nn_gate=NN_GATE):
    n = len(Gm)
    Gd = Gm + np.diag(np.full(n, np.inf))
    nn = Gd.min(axis=1)
    h = float(np.median(nn))
    R = []
    for c in scales:
        V = (Gm <= c * h + 1e-12).sum(1).astype(float)
        R.append(V / np.median(V))
    R = np.array(R)
    rbar = R.mean(0)
    inb = Gm <= scales[-1] * h + 1e-12
    cand = []
    for x in np.argsort(rbar, kind="stable"):
        if not (np.all(R[:, x] < rho_thr) and nn[x] <= nn_gate * h):
            continue
        nb = np.where(inb[x])[0]
        if all((rbar[x] < rbar[y]) or (rbar[x] == rbar[y] and x <= y) for y in nb):
            cand.append(int(x))
    return dict(h=h, rho=rbar, rho_scales=R, cone_sites=cand, cones=len(cand), rho_min=float(rbar.min()),
                cone_rho=[float(rbar[c]) for c in cand])


def eye(I, thr2=THR2_DEFAULT, do_h2=True, keep=False):
    nc, lab = G.components(I)
    sizes = np.bincount(lab)
    c = int(np.argmax(sizes))
    idx = np.where(lab == c)[0]
    d = G.mi_distance(I[np.ix_(idx, idx)])
    Gm, tau = G.geodesic(d)
    rel1, diam, _ = B._h1_rel(Gm)
    out = dict(pieces=int(nc), size=int(len(idx)), diam=diam, rel1=rel1[:5].tolist(), loops=int((rel1 > B.THR).sum()))
    if do_h2:
        rel2, _ = _h2_rel(Gm, diam)
        out.update(rel2=rel2[:5].tolist(), cavities=int((rel2 > thr2).sum()))
    cn = cones(Gm)
    out.update(cones=cn["cones"], cone_rho=cn["cone_rho"], rho_min=cn["rho_min"], h=cn["h"],
               cone_sites=idx[cn["cone_sites"]].tolist())
    if keep:
        out.update(Gm=Gm, sites=idx, rho=cn["rho"])
    return out


def strip(r):
    return {k: v for k, v in r.items() if k not in ("Gm", "sites", "rho")}


# ------------------------------------------------------------------ nulls (R1 classes, fresh seeds)
def null_base():
    from pre_core import LatticeWalk
    n = 24
    wk = LatticeWalk(n, 0.125, periodic=False, w=[np.full((n, n), 0.3)] * 2)
    Ib, gap = B.vacuum_mi(wk)
    return Ib


def make_null(cls, Ib, rng):
    Ns = len(Ib)
    iu = np.triu_indices(Ns, 1)
    off = Ib[iu]
    pos = off > G.I_FLOOR
    mu, sd, frac = np.log(off[pos]).mean(), np.log(off[pos]).std(), pos.mean()
    if cls == "shuffle":
        v = off.copy(); rng.shuffle(v)
    elif cls == "noise":
        v = np.where(rng.random(len(off)) < frac, np.exp(rng.normal(mu, sd, len(off))), 0.0)
    elif cls == "random_graph":
        return mi_of_hopping(random_regular(Ns, 6, rng))
    elif cls == "random_state":
        Z = rng.normal(size=(Ns, Ns)) + 1j * rng.normal(size=(Ns, Ns))
        Q, _ = np.linalg.qr(Z)
        Vo = Q[:, :Ns // 2]
        return mi_one_mode(Vo @ Vo.conj().T)
    I = np.zeros_like(Ib); I[iu] = v; I = I + I.T; np.fill_diagonal(I, np.diag(Ib))
    return I


CLASSES = ("shuffle", "noise", "random_graph", "random_state")

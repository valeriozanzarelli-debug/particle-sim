"""flow_core.py - L2 "the flow": discrete Ricci flow of link LENGTHS driven by the balls' own transport.

Graph = triangulated surface (vertices, edges, triangles); w_e = edge length (Ni-Lin-Gao-Gu 2019 convention).
Balls' transport (lazy random walk): from x stay with prob ALPHA, else hop to neighbour z with prob
    (1 - ALPHA) exp(-w_xz^2) / sum_z' exp(-w_xz'^2)                     (Ni et al. 2019, exponent p = 2)
Ollivier-Ricci curvature  kappa(x,y) = 1 - W1(m_x, m_y) / w_xy, ground cost = shortest-path distance with lengths w.
Flow step (fixed): w <- w - EPS * kappa * w, then rescale so that sum w = initial sum (normalised flow). No surgery.
Stop: max_e |dw_e| / w_e < TOL, or MAXIT steps.
Holonomy (nearest feasible check of orb_core's spinor test on a weighted graph): Regge geometry of the triangulation
(angles from the lengths by the law of cosines). Deficit at a vertex delta_v = 2 pi - sum of corner angles. A 2-spinor
carried once around v is multiplied by exp(i sigma_z delta_v / 2). It comes back equal iff delta_v = 0 (mod 4 pi).
N_inv (invariant spinors of the whole surface) = 2 if every |delta_v| < HOL_TOL, else 0.
Hopping of the balls for the eye: T_xy = exp(1 - w_xy^2) on the links (= 1 at w = 1), same kernel as the walk.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import shortest_path
from scipy.optimize import linprog
from scipy.spatial import Delaunay, ConvexHull

ALPHA, EPS, TOL, MAXIT = 0.5, 0.5, 1e-3, 100
HOL_TOL = 0.05


# ------------------------------------------------------------------ triangulations
class Mesh:
    def __init__(self, N, tris, pos=None, kind=""):
        tris = np.array(sorted({tuple(sorted(t)) for t in np.asarray(tris).tolist()}))
        self.N, self.T, self.pos, self.kind = N, tris, pos, kind
        E = set()
        for a, b, c in tris:
            E |= {(a, b), (b, c), (a, c)}
        self.E = np.array(sorted(E))
        self.eid = {tuple(e): k for k, e in enumerate(self.E)}
        self.nbr = [[] for _ in range(N)]
        for k, (a, b) in enumerate(self.E):
            self.nbr[a].append((b, k)); self.nbr[b].append((a, k))
        self.euler = N - len(self.E) + len(tris)

    def tri_edges(self):
        e = self.eid
        return np.array([[e[(b, c)], e[(a, c)], e[(a, b)]] for a, b, c in self.T])   # edge opposite a, b, c


def hex_mesh(L):
    idx = lambda i, j: (i % L) * L + (j % L)
    tris = []
    for i in range(L):
        for j in range(L):
            tris.append((idx(i, j), idx(i + 1, j), idx(i, j + 1)))
            tris.append((idx(i, j), idx(i + 1, j), idx(i + 1, j - 1)))
    i, j = np.divmod(np.arange(L * L), L)
    pos = np.stack([i + 0.5 * j, np.sqrt(3) / 2 * j], 1)
    return Mesh(L * L, tris, pos, "hex")


def torus_delaunay(N, rng, L=1.0):
    """N uniform points on a rhombic (60 deg) torus of side L, periodic Delaunay triangulation"""
    a1, a2 = np.array([L, 0.0]), np.array([L / 2, L * np.sqrt(3) / 2])
    s = rng.random((N, 2))
    P = s[:, :1] * a1 + s[:, 1:] * a2
    big, ids = [], []
    for u in (-1, 0, 1):
        for v in (-1, 0, 1):
            big.append(P + u * a1 + v * a2); ids.append(np.arange(N))
    big = np.vstack(big); ids = np.concatenate(ids)
    D = Delaunay(big)
    centre = slice(4 * N, 5 * N)
    tris = []
    for t in D.simplices:
        if np.any((t >= 4 * N) & (t < 5 * N)):
            tris.append(tuple(ids[t]))
    m = Mesh(N, tris, P, "torus_delaunay")
    return m


def sphere_hull(N, rng):
    X = rng.normal(size=(N, 3)); X /= np.linalg.norm(X, axis=1)[:, None]
    return Mesh(N, ConvexHull(X).simplices, X, "sphere_hull")


# ------------------------------------------------------------------ geometry
def apsp(m, w):
    A = sp.csr_matrix((np.r_[w, w], (np.r_[m.E[:, 0], m.E[:, 1]], np.r_[m.E[:, 1], m.E[:, 0]])), shape=(m.N, m.N))
    return shortest_path(A, directed=False)


def measure(m, w, x):
    nb = m.nbr[x]
    z = np.array([b for b, _ in nb]); wz = np.array([w[k] for _, k in nb])
    p = np.exp(-(wz**2 - np.min(wz**2))); p = (1 - ALPHA) * p / p.sum()   # shifted: same measure, no underflow
    return np.r_[x, z], np.r_[ALPHA, p]


def w1(cost, a, b):
    """exact W1 (transport LP, HiGHS). If HiGHS presolve wrongly declares the problem infeasible (seen when some
    masses are ~1e-7), retry once without presolve (same LP)."""
    na, nb = len(a), len(b)
    Aeq = np.zeros((na + nb, na * nb))
    for i in range(na):
        Aeq[i, i * nb:(i + 1) * nb] = 1
    for j in range(nb):
        Aeq[na + j, j::nb] = 1
    r = linprog(cost.ravel(), A_eq=Aeq, b_eq=np.r_[a, b], bounds=(0, None), method="highs")
    if r.status != 0:
        r = linprog(cost.ravel(), A_eq=Aeq, b_eq=np.r_[a, b], bounds=(0, None), method="highs",
                    options={"presolve": False})
    if r.status != 0:
        raise RuntimeError(f"W1 linear program failed: {r.message}")
    return r.fun


def ollivier(m, w, D=None):
    if D is None:
        D = apsp(m, w)
    meas = [measure(m, w, x) for x in range(m.N)]
    k = np.empty(len(m.E))
    for e, (x, y) in enumerate(m.E):
        sx, px = meas[x]; sy, py = meas[y]
        k[e] = 1 - w1(D[np.ix_(sx, sy)], px, py) / w[e]
    return k


def vertex_curv(m, k):
    out = np.zeros(m.N); cnt = np.zeros(m.N)
    np.add.at(out, m.E[:, 0], k); np.add.at(out, m.E[:, 1], k)
    np.add.at(cnt, m.E[:, 0], 1); np.add.at(cnt, m.E[:, 1], 1)
    return out / cnt


def deficits(m, w):
    """Regge deficit per vertex; NaN where a triangle violates the triangle inequality"""
    te = m.tri_edges()
    a, b, c = w[te[:, 0]], w[te[:, 1]], w[te[:, 2]]          # opposite vertex 0, 1, 2
    def ang(o, p, q):
        x = (p**2 + q**2 - o**2) / (2 * p * q)
        return np.where(np.abs(x) <= 1, np.arccos(np.clip(x, -1, 1)), np.nan)
    A = np.stack([ang(a, b, c), ang(b, a, c), ang(c, a, b)], 1)
    s = np.zeros(m.N)
    for col in range(3):
        np.add.at(s, m.T[:, col], A[:, col])
    return 2 * np.pi - s, int(np.isnan(A).any(axis=1).sum())


def holonomy(m, w):
    d, bad = deficits(m, w)
    ok = np.isfinite(d)
    passing = ok & (np.abs(d) < HOL_TOL)
    return dict(N_inv=int(2 if passing.all() else 0), frac_vertices_pass=float(passing.mean()),
                max_abs_deficit=float(np.nanmax(np.abs(d))) if ok.any() else float("nan"),
                sum_deficit_over_2pi=float(np.nansum(d) / (2 * np.pi)), bad_triangles=bad), d


def hopping(m, w):
    T = np.zeros((m.N, m.N))
    t = np.exp(1 - w**2)
    T[m.E[:, 0], m.E[:, 1]] = t; T[m.E[:, 1], m.E[:, 0]] = t
    return T


# ------------------------------------------------------------------ the flow
def flow(m, w0, log=print, maxit=MAXIT):
    import time
    t0 = time.time()
    w = np.array(w0, float); total = w.sum()
    hist = []
    k = ollivier(m, w)
    k0 = k.copy()
    stop = "maxit"
    for it in range(maxit):
        hist.append(dict(it=it, k_mean=float(k.mean()), k_std=float(k.std()), w_std=float(w.std() / w.mean()),
                         w_min=float(w.min()), w_max=float(w.max())))
        wn = w - EPS * k * w
        wn *= total / wn.sum()
        rel = float(np.max(np.abs(wn - w) / w))
        w = wn
        k = ollivier(m, w)
        log(f"  it {it:3d} k_mean {hist[-1]['k_mean']:+.4f} k_std {hist[-1]['k_std']:.4f} "
            f"w_std {hist[-1]['w_std']:.4f} maxrel {rel:.2e} {time.time()-t0:.0f}s")
        if rel < TOL:
            stop = "tol"
            break
    hist.append(dict(it=len(hist), k_mean=float(k.mean()), k_std=float(k.std()), w_std=float(w.std() / w.mean()),
                     w_min=float(w.min()), w_max=float(w.max())))
    return w, k0, k, hist, stop

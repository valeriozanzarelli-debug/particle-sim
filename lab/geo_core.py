"""geo_core.py - "space from the balls": read the geometry OFF the zig-zag rule, do not feed it in.

Valerio's principle: we do not use the shape to move the balls; we move the balls and ask which shape
they follow. Here:
  1. a HIDDEN shape = a set of sites on a square grid (mask) with, per axis, open or periodic ends.
     The shape is used ONLY to build the one-tick unitary of the ball rule (zz_core: translate at c
     along each axis on the eigenvectors of sigma_j, L opposite to R, then flip L<->R by theta).
     At an open edge the walker that would leave turns around (moving component -> opposite moving
     component at the same site, sign -1, as RadialWalk at r = 0): a permutation, so U stays unitary.
  2. the many-ball state is a free-fermion (Gaussian) state: either the VACUUM (Dirac sea: every
     single-ball eigenmode of U with quasi-energy omega in (-pi, 0) is filled) or a QUENCH (random
     local product state, then evolved by the rule). Everything is exact through the correlation
     matrix C (projector on the filled single-ball states).
  3. the RECONSTRUCTION sees only the mutual-information matrix I(i, j) between sites (4 fermion
     modes per site). It is NOT told the adjacency, the coordinates, the dimension or the topology.
     distance d -> threshold graph -> geodesics (Isomap) -> classical MDS embedding, ball-growth
     dimension, Vietoris-Rips persistent homology (ripser) -> Betti numbers and the hidden radius.
Debt (declared): the alternating axis order of zz_core (x,y on even ticks, y,x on odd) is replaced by a
fixed order x then y. Reason: with alternation the period is two ticks and the two-tick phases of the
positive and negative bands coincide (omega and omega - pi alias), so "the vacuum" is not defined.
With a fixed order the one-tick U has gaps at omega = 0 and omega = pi and the sea is unambiguous.
Units: lattice a = tick = c = hbar = 1.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import shortest_path, connected_components
from zz_core import BASIS_HALF

SEED = 20261009


# ------------------------------------------------------------------ hidden shapes
class Shape:
    """sites = True entries of mask (nx, ny); periodic = (px, py); axes = lattice axes that exist.
    true_xy(i) and true_dist() are used ONLY for grading, never by the reconstruction."""

    def __init__(self, name, mask, periodic, axes):
        self.name, self.mask = name, np.asarray(mask, bool)
        self.periodic, self.axes = periodic, axes
        self.nx, self.ny = self.mask.shape
        self.xy = np.argwhere(self.mask)                     # (Ns, 2) integer coordinates
        self.Ns = len(self.xy)
        self.index = -np.ones(self.mask.shape, int)
        self.index[self.mask] = np.arange(self.Ns)

    def neighbour(self, i, j, s):
        """site reached from site i by s steps along axis j, or -1 if it leaves the shape"""
        p = self.xy[i].copy()
        p[j] += s
        n = (self.nx, self.ny)[j]
        if self.periodic[j]:
            p[j] %= n
        elif not (0 <= p[j] < n):
            return -1
        return self.index[p[0], p[1]]

    def true_dist(self):
        """flat Euclidean geodesic distance inside the shape (8-neighbour grid graph, weights 1 and
        sqrt 2; minimum image on periodic axes). Grading only."""
        rows, cols, w = [], [], []
        for i in range(self.Ns):
            for (dx, dy) in ((1, 0), (0, 1), (1, 1), (1, -1)):
                p = self.xy[i] + (dx, dy)
                ok = True
                for j, n in enumerate((self.nx, self.ny)):
                    if self.periodic[j]:
                        p[j] %= n
                    elif not (0 <= p[j] < n):
                        ok = False
                if not ok or not self.mask[p[0], p[1]]:
                    continue
                if dx and dy:   # a diagonal step must not cut a corner of the shape
                    a = self.xy[i] + (dx, 0); b = self.xy[i] + (0, dy)
                    a %= (self.nx, self.ny); b %= (self.nx, self.ny)
                    if not (self.mask[a[0], a[1]] and self.mask[b[0], b[1]]):
                        continue
                k = self.index[p[0], p[1]]
                rows.append(i); cols.append(k); w.append(np.hypot(dx, dy))
        G = sp.coo_matrix((w, (rows, cols)), shape=(self.Ns, self.Ns)).tocsr()
        return shortest_path(G, directed=False)


def chain(n):
    return Shape(f"catena {n}", np.ones((n, 1)), (False, True), (0,))


def ring(n):
    return Shape(f"anello {n}", np.ones((n, 1)), (True, True), (0,))


def torus(n1, n2):
    return Shape(f"toro {n1}x{n2}", np.ones((n1, n2)), (True, True), (0, 1))


def cylinder(L, N):
    """visible open line of L sites x hidden circle of N sites"""
    return Shape(f"linea {L} x cerchio {N}", np.ones((L, N)), (False, True), (0, 1))


def stadium_with_hole(W=30, H=16, hole_c=(19.0, 9.5), hole_r=3.6):
    """irregular hidden shape: a stadium (two half discs + rectangle) with an off-centre round hole"""
    r = H / 2.0
    x = np.arange(W)[:, None] + 0.5; y = np.arange(H)[None, :] + 0.5
    xc = np.clip(x, r, W - r)
    inside = (x - xc)**2 + (y - r)**2 <= r**2
    hole = (x - hole_c[0])**2 + (y - hole_c[1])**2 <= hole_r**2
    return Shape(f"stadio con buco {W}x{H}", inside & ~hole, (False, False), (0, 1))


# ------------------------------------------------------------------ the ball rule on a shape
def axis_operator(shape, j):
    """one-axis translation of the zig-zag rule on the shape, as a sparse unitary on (4 comps x Ns).
    comps 0,1 = L (moves along -S_j), 2,3 = R (moves along +S_j); S_j = sigma_j (BASIS_HALF)."""
    Ns = shape.Ns
    lam, P = BASIS_HALF[j]                   # lam = (+1, -1); P columns = eigenvectors of sigma_j
    blocks = []
    for h in (-1, +1):                       # L block, R block
        # in the eigen-basis: component q moves by s = h*lam[q]; at an open edge it turns into the
        # other component at the same site with sign -1 (reflection)
        rows, cols, vals = [], [], []
        for q in range(2):
            s = h * lam[q]
            for i in range(Ns):
                k = shape.neighbour(i, j, s)
                if k >= 0:
                    rows.append(q * Ns + k); cols.append(q * Ns + i); vals.append(1.0)
                else:
                    rows.append((1 - q) * Ns + i); cols.append(q * Ns + i); vals.append(-1.0)
        Perm = sp.csr_matrix((vals, (rows, cols)), shape=(2 * Ns, 2 * Ns))
        Pk = sp.kron(sp.csr_matrix(P), sp.identity(Ns), format="csr")
        blocks.append(Pk @ Perm @ Pk.conj().T)
    return sp.block_diag(blocks, format="csr")


def rule_unitary(shape, theta):
    """one tick: translate along each existing axis (fixed order), then flip L<->R by theta"""
    Ns = shape.Ns
    U = sp.identity(4 * Ns, dtype=complex, format="csr")
    for j in shape.axes:
        U = axis_operator(shape, j) @ U
    ct, st = np.cos(theta), np.sin(theta)
    I2 = sp.identity(2 * Ns, format="csr")
    Cflip = sp.bmat([[ct * I2, -1j * st * I2], [-1j * st * I2, ct * I2]], format="csr")
    return (Cflip @ U).tocsr()


def vacuum(U):
    """Dirac sea: projector on eigenmodes with omega in (-pi, 0), omega = -arg(lambda).
    U is normal and U v = e^{-i omega} v gives (U - U^dag)/(2i) v = -sin(omega) v, so the sea is the
    span of the POSITIVE eigenvalues of that Hermitian matrix (orthonormal eigenvectors via eigh).
    Returns C (projector) and the smallest |sin omega| (gap check)."""
    Ud = U.toarray()
    Hs = (Ud - Ud.conj().T) / 2j
    w, V = np.linalg.eigh(Hs)
    occ = w > 0
    Vo = V[:, occ]
    return Vo @ Vo.conj().T, float(np.min(np.abs(w))), int(occ.sum())


def random_product_state(Ns, rng):
    """each site: a random (Haar) 2-dim subspace of its 4 local modes is filled; no correlations"""
    C = np.zeros((4 * Ns, 4 * Ns), complex)
    for i in range(Ns):
        Z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        Q, R = np.linalg.qr(Z)
        Q = Q * (np.diag(R) / np.abs(np.diag(R)))[None, :]
        V = Q[:, :2]
        idx = np.arange(4) * Ns + i
        C[np.ix_(idx, idx)] = V @ V.conj().T
    return C


# ------------------------------------------------------------------ mutual information
def _h(nu):
    nu = np.clip(nu.real, 1e-300, 1.0)
    mu = np.clip(1.0 - nu, 1e-300, 1.0)
    return -(nu * np.log(nu) + mu * np.log(mu))


def mutual_information(C, Ns, pairs=None, chunk=20000):
    """I(i,j) = S_i + S_j - S_ij for a Gaussian fermion state with correlation matrix C.
    Each site owns modes {a*Ns + i, a = 0..3}. Returns full symmetric matrix (diag = 2 S_i) or values
    on the requested pairs."""
    modes = np.arange(4)[None, :] * Ns + np.arange(Ns)[:, None]          # (Ns, 4)
    Cs = C[modes[:, :, None], modes[:, None, :]]
    S1 = _h(np.linalg.eigvalsh(Cs)).sum(axis=1)
    if pairs is None:
        ii, jj = np.triu_indices(Ns, 1)
    else:
        ii, jj = pairs
    out = np.empty(len(ii))
    for a in range(0, len(ii), chunk):
        pm = np.concatenate([modes[ii[a:a + chunk]], modes[jj[a:a + chunk]]], axis=1)   # (P, 8)
        Cp = C[pm[:, :, None], pm[:, None, :]]
        out[a:a + chunk] = S1[ii[a:a + chunk]] + S1[jj[a:a + chunk]] - _h(np.linalg.eigvalsh(Cp)).sum(axis=1)
    if pairs is not None:
        return out, S1
    M = np.zeros((Ns, Ns))
    M[ii, jj] = out; M[jj, ii] = out
    M[np.arange(Ns), np.arange(Ns)] = 2 * S1
    return M


# ------------------------------------------------------------------ blind reconstruction
I_FLOOR = 1e-10   # PRE-REGISTERED: below this the mutual information is "zero" (numerical noise ~1e-15)


def components(I, floor=I_FLOOR):
    """PRE-REGISTERED: the pieces of space that share ANY information (graph I_ij > floor)"""
    A = sp.csr_matrix(I > floor)
    return connected_components(A, directed=False)


def mi_distance(I, floor=I_FLOOR):
    """PRE-REGISTERED map: d_ij = -log(I_ij / I_0), I_0 = mean self-information 2 S_i (so d_ii = 0).
    Pairs with I below the floor get d = -log(floor / I_0) (far)."""
    I0 = np.mean(np.diag(I))
    d = -np.log(np.maximum(I, floor) / I0)
    np.fill_diagonal(d, 0.0)
    return np.maximum(d, 0.0)


def geodesic(d, factor=1.5):
    """PRE-REGISTERED Isomap graph: tau = smallest threshold that makes the graph {d_ij <= tau}
    connected; keep the edges with d_ij <= factor * tau; geodesic = shortest paths."""
    vals = np.unique(np.round(d[np.triu_indices(len(d), 1)], 12))
    lo, hi = 0, len(vals) - 1
    while lo < hi:                           # bisection on the connectivity threshold
        mid = (lo + hi) // 2
        A = sp.csr_matrix(d <= vals[mid])
        if connected_components(A, directed=False)[0] == 1:
            hi = mid
        else:
            lo = mid + 1
    tau = vals[lo]
    W = np.where((d <= factor * tau) & (d > 0), d, 0.0)
    G = shortest_path(sp.csr_matrix(W), directed=False)
    return G, float(tau)


def mds(G, k=3):
    """classical MDS; returns coordinates (n, k) and all eigenvalues of the centred Gram matrix"""
    n = len(G)
    J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ (G**2) @ J
    w, V = np.linalg.eigh(B)
    o = np.argsort(w)[::-1]
    w, V = w[o], V[:, o]
    X = V[:, :k] * np.sqrt(np.maximum(w[:k], 0))
    return X, w


def _slope(Gd, lo, hi):
    rs = np.geomspace(lo, hi, 10)
    Nr = np.array([np.mean(np.sum(Gd <= r, axis=1)) for r in rs])
    ok = Nr > 0
    if ok.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log(rs[ok]), np.log(Nr[ok]), 1)[0])


def ball_dimension(G):
    """PRE-REGISTERED intrinsic dimension: slope of log <N(r)> vs log r, N(r) = number of OTHER points
    within geodesic distance r. h = median nearest-neighbour distance, diam = largest geodesic.
      D_small: fit on r in [1.5 h, 4 h]        (short distances)
      D_large: fit on r in [diam/4, diam/2]    (long distances)
    A compact hidden direction should show up as D_small > D_large."""
    n = len(G)
    Gd = G + np.diag(np.full(n, np.inf))
    h = float(np.median(Gd.min(axis=1)))
    diam = float(np.max(G[np.isfinite(G)]))
    return _slope(Gd, 1.5 * h, 4 * h), _slope(Gd, diam / 4, diam / 2), h, diam


def dimension_curve(G, nr=30):
    """local slope d log N / d log r (for the figure)"""
    n = len(G)
    Gd = G + np.diag(np.full(n, np.inf))
    h = np.median(Gd.min(axis=1)); diam = np.max(G[np.isfinite(G)])
    rs = np.geomspace(h, diam, nr)
    Nr = np.array([np.mean(np.sum(Gd <= r, axis=1)) for r in rs])
    return rs, Nr


N_LANDMARKS = 220   # PRE-REGISTERED: Rips homology on at most this many max-min landmarks


def landmarks(G, m=N_LANDMARKS):
    """greedy max-min (farthest point) subsample, deterministic start at point 0"""
    n = len(G)
    if n <= m:
        return np.arange(n)
    Gf = np.where(np.isfinite(G), G, 1e9)
    L = [0]
    dmin = Gf[0].copy()
    for _ in range(m - 1):
        k = int(np.argmax(dmin))
        L.append(k)
        dmin = np.minimum(dmin, Gf[k])
    return np.array(sorted(L))


def homology(G, h, maxdim=2, min_pers=1.5):
    """PRE-REGISTERED topology: Vietoris-Rips persistent homology of the geodesic matrix (ripser).
    b0, b1 (primary): all points; a bar counts if persistence (death - birth) > min_pers * h
    (h = median nearest-neighbour distance); the H0 bar that never dies counts.
    b2 (secondary, cost): on max-min landmarks (<= N_LANDMARKS points), counted if persistence >
    min_pers * 2 eps_L (eps_L = covering radius of the landmarks: the Rips stability bound).
    Returns [b0, b1, b2], diagrams [H0, H1, H2] and the scale h."""
    from ripser import ripser
    dg = ripser(G, distance_matrix=True, maxdim=1)["dgms"]
    betti, bars = [], []
    for D in dg:
        pers = D[:, 1] - D[:, 0]
        betti.append(int(((pers > min_pers * h) | ~np.isfinite(D[:, 1])).sum()))
        bars.append(D)
    if maxdim >= 2:
        L = landmarks(G)
        eps = float(np.max(np.min(np.where(np.isfinite(G[:, L]), G[:, L], 1e9), axis=1)))
        D2 = ripser(G[np.ix_(L, L)], distance_matrix=True, maxdim=2)["dgms"][2]
        pers = D2[:, 1] - D2[:, 0]
        betti.append(int(((pers > min_pers * max(h, 2 * eps)) | ~np.isfinite(D2[:, 1])).sum()))
        bars.append(D2)
    return betti, bars, h


def hidden_radius(bars, h, min_pers=1.5):
    """PRE-REGISTERED hidden radius: a geodesic circle of length Lc has its Rips H1 bar dying at Lc/3
    (Adamaszek-Adams). The HIDDEN circle = the significant H1 bar that dies FIRST (the smallest hole).
    R_est = 3 * death / (2 pi), same units as G. Also returns the death of the largest significant hole."""
    D = bars[1]
    if len(D) == 0:
        return float("nan"), float("nan")
    pers = D[:, 1] - D[:, 0]
    sig = pers > min_pers * h
    if not sig.any():
        return float("nan"), float("nan")
    deaths = D[sig, 1]
    return 3.0 * float(np.min(deaths)) / (2 * np.pi), float(np.max(deaths))


def causal_distance(U, Ns):
    """the shape the rule ACTUALLY connects: site graph with an edge i - j when one tick of U carries
    amplitude from any mode of i to any mode of j; distance = number of ticks (BFS). Reference only."""
    A = abs(U).tocoo()
    si, sj = A.row % Ns, A.col % Ns
    M = sp.csr_matrix((np.ones(len(si)), (si, sj)), shape=(Ns, Ns))
    M = ((M + M.T) > 0).astype(float)
    M.setdiag(0); M.eliminate_zeros()
    return shortest_path(M, directed=False, unweighted=True)


def fidelity(G, Dtrue):
    """Pearson and Spearman correlation between reconstructed and true distances (all pairs)"""
    from scipy.stats import pearsonr, spearmanr
    iu = np.triu_indices(len(G), 1)
    a, b = G[iu], Dtrue[iu]
    ok = np.isfinite(a)
    return float(pearsonr(a[ok], b[ok])[0]), float(spearmanr(a[ok], b[ok])[0])


def reconstruct_from_I(I, maxdim=2):
    """blind pipeline from the mutual-information matrix: split into informational components, then
    reconstruct each component (sites of the component in the 'sites' entry, for grading only)"""
    nc, lab = components(I)
    out = []
    for c in range(nc):
        idx = np.where(lab == c)[0]
        r = reconstruct(mi_distance(I[np.ix_(idx, idx)]), maxdim=maxdim)
        r["sites"] = idx
        out.append(r)
    return nc, out


def reconstruct(d, maxdim=2):
    """the whole blind pipeline from a distance matrix"""
    G, tau = geodesic(d)
    dsmall, dlarge, h, diam = ball_dimension(G)
    betti, bars, scale = homology(G, h, maxdim=maxdim)
    R, big = hidden_radius(bars, scale)
    X, ev = mds(G, 3)
    return dict(G=G, tau=tau, dim=dsmall, dim_large=dlarge, h=h, scale=scale, diam=diam, betti=betti, bars=bars,
                R=R, big_death=big, X=X, ev=ev)


# ------------------------------------------------------------------ dynamic version (quench)
EPS_ARRIVAL = 1e-3   # PRE-REGISTERED: correlation has "arrived" when I_ij first exceeds this (nats)


def arrival_times(U, Ns, T, seeds, eps=EPS_ARRIVAL):
    """start from random local product states (no correlations at all), let the rule run, and record
    for every pair the (log-interpolated) tick at which I_ij first exceeds eps; I is averaged over the
    seeds first. Never arrived within T ticks -> inf. Returns the arrival matrix (diag 0)."""
    rng = np.random.default_rng(SEED)
    Cs = [random_product_state(Ns, rng) for _ in seeds]
    ii, jj = np.triu_indices(Ns, 1)
    arr = np.full(len(ii), np.inf)
    prev = np.zeros(len(ii))
    Uh = U.conj().T.tocsr()
    for t in range(1, T + 1):
        Cs = [np.asarray(U @ C @ Uh) for C in Cs]
        todo = np.where(~np.isfinite(arr))[0]
        if len(todo) == 0:
            break
        Iv = np.mean([mutual_information(C, Ns, pairs=(ii[todo], jj[todo]))[0] for C in Cs], axis=0)
        hit = Iv > eps
        k = todo[hit]
        a = np.log(np.maximum(prev[k], 1e-30)); b = np.log(Iv[hit])
        frac = np.clip((np.log(eps) - a) / np.maximum(b - a, 1e-12), 0, 1)
        arr[k] = (t - 1) + frac
        prev[todo] = Iv
    A = np.zeros((Ns, Ns))
    A[ii, jj] = arr; A[jj, ii] = arr
    return A


def reconstruct_from_arrival(A, maxdim=2):
    """blind pipeline for the dynamic version: pieces = pairs that ever correlate; d = arrival time"""
    nc, lab = connected_components(sp.csr_matrix(np.isfinite(A) & (A > 0)), directed=False)
    out = []
    for c in range(nc):
        idx = np.where(lab == c)[0]
        d = A[np.ix_(idx, idx)].copy()
        d[~np.isfinite(d)] = 1e6
        r = reconstruct(d, maxdim=maxdim)
        r["sites"] = idx
        out.append(r)
    return nc, out

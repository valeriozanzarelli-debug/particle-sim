"""geo_drums.py - Part B: different shapes, same masses (Gordon-Webb-Wolpert 1992), and what breaks the tie.

Construction (Sunada 1985; Buser, Conway, Doyle, Semmler 1994): take 7 copies of a half-square
(right isosceles triangle) with edges a, b (legs) and c (hypotenuse). Pick three involutions A, B, C of
the group GL(3,2) (the symmetries of the Fano plane, 7 points and 7 lines). Drum 1: tile t is glued
along its edge a to tile A(t) (A acting on the POINTS), etc.; a tile fixed by A keeps edge a as boundary.
Drum 2: the same with A, B, C acting on the LINES. The two actions are 'almost conjugate' (Gassmann
triple), so the Dirichlet spectra are equal (transplantation). We search the first triple whose two
gluing graphs are trees and whose two tilings do not overlap in the plane.
Then: P1 finite elements, each tile refined into n^2 congruent triangles (the same reference mesh in
every tile, so even the discrete spectra must coincide), first 12 eigenvalues for n = 16 and 32, and
three 'coupling' numbers that are NOT spectral: g3 = int u1^3, g4 = int u1^4, |int u1 u2 u3|
(in a Kaluza-Klein reduction these are the cubic / quartic couplings between the modes = 'particles').
Writes lab/out/geo_drums.json, .log, .npz (polygons and first modes for the figure).
"""
import os, sys, json, time, itertools
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# ---------------------------------------------------------------- GL(3,2) and its two actions
VEC = [np.array(v) for v in itertools.product((0, 1), repeat=3) if any(v)]
KEY = {tuple(v): i for i, v in enumerate(VEC)}


def gl32():
    els = []
    for bits in itertools.product((0, 1), repeat=9):
        M = np.array(bits).reshape(3, 3)
        if round(np.linalg.det(M)) % 2 == 1:
            els.append(M)
    return els


def inv2(M):
    for N in GL:
        if np.array_equal((M @ N) % 2, np.eye(3, dtype=int)):
            return N


GL = gl32()
INV = [M for M in GL if np.array_equal((M @ M) % 2, np.eye(3, dtype=int)) and not np.array_equal(M, np.eye(3))]


def perm_points(M):
    return [KEY[tuple((M @ v) % 2)] for v in VEC]


def perm_lines(M):
    Mt = inv2(M).T
    return [KEY[tuple((Mt @ v) % 2)] for v in VEC]


# ---------------------------------------------------------------- tiling by reflections
REFL = {  # reflections of the reference triangle O(0,0) A(1,0) B(0,1) across its edges, as (L, c)
    "a": (np.array([[1, 0], [0, -1]]), np.array([0, 0])),     # edge OA (y = 0)
    "b": (np.array([[-1, 0], [0, 1]]), np.array([0, 0])),     # edge OB (x = 0)
    "c": (np.array([[0, -1], [-1, 0]]), np.array([1, 1])),    # edge AB (x + y = 1)
}
REF = np.array([[0, 0], [1, 0], [0, 1]], float)


def tile_maps(perms):
    """perms = {'a': sigma_a, 'b': ..., 'c': ...}; returns affine maps of the 7 tiles or None"""
    edges = set()
    for e, s in perms.items():
        for t in range(7):
            if s[t] != t:
                edges.add((min(t, s[t]), max(t, s[t]), e))
    if len(edges) != 6:
        return None
    maps = {0: (np.eye(2), np.zeros(2))}
    todo = [0]
    while todo:
        t = todo.pop()
        Mt, ct = maps[t]
        for e, s in perms.items():
            u = s[t]
            if u == t or u in maps:
                continue
            L, c = REFL[e]
            maps[u] = (Mt @ L, Mt @ c + ct)           # g_u = g_t o r_e
            todo.append(u)
    if len(maps) != 7:
        return None
    return maps


def overlaps(maps):
    g = (np.arange(-60, 60) + 0.5) / 8.0
    X, Y = np.meshgrid(g, g)
    P = np.stack([X.ravel(), Y.ravel()], 1)
    cover = np.zeros(len(P), int)
    for t, (M, c) in maps.items():
        q = (P - c) @ np.linalg.inv(M).T            # back to reference coordinates
        cover += (q[:, 0] > 0) & (q[:, 1] > 0) & (q[:, 0] + q[:, 1] < 1)
    return cover.max() > 1


def _fingerprint(m):
    P = np.array([M @ REF.mean(0) + c for (M, c) in m.values()])
    return np.sort(np.round(np.linalg.norm(P[:, None] - P[None], axis=2).ravel(), 9))


def congruent(m1, m2):
    """isometry-invariant test: sorted pairwise distances of the tile centroids"""
    return np.allclose(_fingerprint(m1), _fingerprint(m2))


def find_pair():
    for A, B, C in itertools.product(INV, repeat=3):
        pp = {"a": perm_points(A), "b": perm_points(B), "c": perm_points(C)}
        pl = {"a": perm_lines(A), "b": perm_lines(B), "c": perm_lines(C)}
        m1, m2 = tile_maps(pp), tile_maps(pl)
        if m1 is None or m2 is None or overlaps(m1) or overlaps(m2):
            continue
        if congruent(m1, m2):      # added after run 1 (geo_drums_run1_congruent.log): the first triple
            continue               # gave two congruent copies of one drum - a trivial 'pair'
        return pp, pl, m1, m2
    return None


# ---------------------------------------------------------------- P1 finite elements
def fem(perms, maps, n):
    ref_nodes = [(i, j) for j in range(n + 1) for i in range(n + 1 - j)]
    rid = {p: k for k, p in enumerate(ref_nodes)}
    tris = []
    for j in range(n):
        for i in range(n - j):
            tris.append((rid[(i, j)], rid[(i + 1, j)], rid[(i, j + 1)]))
            if i + j < n - 1:
                tris.append((rid[(i + 1, j)], rid[(i + 1, j + 1)], rid[(i, j + 1)]))
    gid, coords, bnd, elems = {}, [], set(), []
    for t in range(7):
        M, c = maps[t]
        local = []
        for (i, j) in ref_nodes:
            x = M @ np.array([i / n, j / n]) + c
            key = tuple(np.round(x * 2 * n).astype(int))
            if key not in gid:
                gid[key] = len(coords); coords.append(x)
            g = gid[key]
            local.append(g)
            on = {"a": j == 0, "b": i == 0, "c": i + j == n}
            if any(on[e] and perms[e][t] == t for e in "abc"):
                bnd.add(g)
        elems += [(local[a], local[b], local[d]) for (a, b, d) in tris]
    X = np.array(coords); E = np.array(elems)
    N = len(X)
    rows, cols, kv, mv = [], [], [], []
    for tri in E:
        p = X[tri]
        Bm = np.array([p[1] - p[0], p[2] - p[0]]).T
        area = abs(np.linalg.det(Bm)) / 2
        G = np.linalg.inv(Bm).T @ np.array([[-1, 1, 0], [-1, 0, 1]])
        Ke = area * G.T @ G
        Me = area / 12 * (np.ones((3, 3)) + np.eye(3))
        for a in range(3):
            for b in range(3):
                rows.append(tri[a]); cols.append(tri[b]); kv.append(Ke[a, b]); mv.append(Me[a, b])
    K = sp.csr_matrix((kv, (rows, cols)), shape=(N, N))
    Mm = sp.csr_matrix((mv, (rows, cols)), shape=(N, N))
    free = np.array(sorted(set(range(N)) - bnd))
    Kf, Mf = K[free][:, free], Mm[free][:, free]
    w, V = eigsh(Kf, k=12, M=Mf, sigma=0, which="LM")
    o = np.argsort(w); w, V = w[o], V[:, o]
    U = np.zeros((N, 12)); U[free] = V
    lump = np.asarray(Mm.sum(axis=1)).ravel()
    for k in range(12):                               # normalize int u^2 = 1, sign: int u > 0
        U[:, k] /= np.sqrt(U[:, k] @ (Mm @ U[:, k]))
        if lump @ U[:, k] < 0:
            U[:, k] *= -1
    u1, u2, u3 = U[:, 0], U[:, 1], U[:, 2]
    coup = dict(g3=float(lump @ u1**3), g4=float(lump @ u1**4), g123=float(abs(lump @ (u1 * u2 * u3))),
                umax=float(np.abs(u1).max()))
    return w, coup, X, E, U[:, :3]


def outline(maps, perms):
    segs = []
    vid = {"a": (0, 1), "b": (0, 2), "c": (1, 2)}
    for t, (M, c) in maps.items():
        for e in "abc":
            if perms[e][t] == t:
                i, j = vid[e]
                segs.append([(M @ REF[i] + c).tolist(), (M @ REF[j] + c).tolist()])
    return segs


def main():
    t0 = time.time()
    log = []
    say = lambda s: (print(s), log.append(s))
    found = find_pair()
    if found is None:
        say("no admissible triple found"); return
    pp, pl, m1, m2 = found
    say(f"gluing (points): {pp}\ngluing (lines):  {pl}")
    # are the two drums congruent? compare the sorted pairwise distances of the tile centroids
    cen = lambda m: np.array([M @ REF.mean(0) + c for (M, c) in m.values()])
    dd = lambda P: np.sort(np.linalg.norm(P[:, None] - P[None], axis=2).ravel())
    say(f"centroid-distance fingerprints differ: {not np.allclose(dd(cen(m1)), dd(cen(m2)))}")
    res = {"gluing_points": pp, "gluing_lines": pl, "n": {}}
    npz = {}
    for n in (16, 32):
        w1, c1, X1, E1, U1 = fem(pp, m1, n)
        w2, c2, X2, E2, U2 = fem(pl, m2, n)
        rel = np.abs(w1 / w2 - 1)
        say(f"n={n}: eigenvalues drum 1 {np.round(w1[:6], 5).tolist()}")
        say(f"      eigenvalues drum 2 {np.round(w2[:6], 5).tolist()}")
        say(f"      max relative difference (12 modes) {rel.max():.1e}")
        say(f"      couplings drum 1 {c1}")
        say(f"      couplings drum 2 {c2}")
        res["n"][n] = dict(eig1=w1.tolist(), eig2=w2.tolist(), maxrel=float(rel.max()), coup1=c1, coup2=c2)
        if n == 32:
            npz = dict(X1=X1, E1=E1, U1=U1, X2=X2, E2=E2, U2=U2)
    # pre-registered verdicts
    a, b = res["n"][16], res["n"][32]
    res["prereg_isospectral"] = bool(a["maxrel"] < 1e-8 and b["maxrel"] < 1e-8)
    tie = {}
    for key in ("g3", "g4", "g123"):
        diff = abs(b["coup1"][key] - b["coup2"][key])
        disc = max(abs(b["coup1"][key] - a["coup1"][key]), abs(b["coup2"][key] - a["coup2"][key]))
        tie[key] = dict(drum_difference=diff, discretization_change=disc, broken=bool(diff > 10 * disc))
    res["tie_breakers"] = tie
    say(f"PRE-REG isospectral (< 1e-8 at same n): {res['prereg_isospectral']}")
    say(f"PRE-REG tie broken by couplings: {tie}")
    res["outline1"], res["outline2"] = outline(m1, pp), outline(m2, pl)
    with open(os.path.join(OUT, "geo_drums.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "geo_drums.log"), "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")
    np.savez_compressed(os.path.join(OUT, "geo_drums.npz"), **npz)
    say(f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()

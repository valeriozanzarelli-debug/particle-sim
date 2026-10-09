"""orb_o3.py - O3: hidden-circle twists (Wilson lines / AB phases) on T2/Z3 and how they change light-mode counts.
(a) flux-free T2/Z3 lattice: which AB phases (phi1, phi2) on the two cycles admit a Z3 operator commuting with the hopping;
(b) magnetised T2/Z3 (uniform AB flux, M quanta): lowest Landau level (= zero modes of the chiral Dirac operator)
    split by the Z3 operator into (n0, n1, n2), for every allowed Wilson line;
(c) T6/Z3 fixed-point classes under 0/1/2/3 Z3 Wilson lines (bookkeeping).
Writes lab/out/orb_o3.json."""
import json, time, itertools
import numpy as np
from orb_core import *

t0 = time.time()
R = {}
DIRS = [(1, 0), (0, 1), (-1, 1)]


def hopping(L, M, phi1, phi2, sgn=1):
    """triangular-lattice hopping H = -sum e^{i theta} |x+e><x| + h.c., Landau gauge A = (0, alpha i) in skew coords,
    alpha = 2 pi M / L^2, magnetic-translation gauge on the i boundary; AB phases phi1 (i-cycle), phi2 (j-cycle)."""
    N = L * L
    al = 2 * np.pi * M / L ** 2
    H = np.zeros((N, N), complex)
    links = []
    for i in range(L):
        for j in range(L):
            for di, dj in DIRS:
                i2, j2 = i + di, j + dj
                th = al * (i + di / 2) * dj
                wi, wj = (i2 // L), (j2 // L)
                th += -sgn * wi * al * L * (j2 % L)           # magnetic translation gauge
                th += wi * phi1 + wj * phi2
                u, v = i * L + j, (i2 % L) * L + (j2 % L)
                H[v, u] += -np.exp(1j * th)
                links.append((u, v))
    H = H + H.conj().T
    return H, links


def flux_check(H, L):
    """phase of every elementary triangle"""
    idx = hex_index(L)
    ph = []
    for i in range(L):
        for j in range(L):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i, j + 1), idx(i + 1, j + 1)
            for tri in ((a, b, c), (b, d, c)):
                p = -H[tri[1], tri[0]] * -H[tri[2], tri[1]] * -H[tri[0], tri[2]]
                ph.append(np.angle(p))
    return np.array(ph)


def z3_operator(H, L, links):
    """U = G* P with P the 120-degree permutation and G diagonal such that U H U^dag = H; None if impossible."""
    N = L * L
    perm = np.zeros(N, int)
    for i in range(L):
        for j in range(L):
            perm[i * L + j] = ((-i - j) % L) * L + (i % L)
    Hp = np.zeros_like(H); Hp[np.ix_(perm, perm)] = H        # Hp[Rx, Ry] = H[x, y]
    g = np.full(N, np.nan + 0j); g[0] = 1
    adj = [[] for _ in range(N)]
    for u, v in links:
        adj[u].append(v); adj[v].append(u)
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if np.isnan(g[v]):          # Hp[u,v] = g_u H[u,v] conj(g_v)
                g[v] = np.conj(Hp[u, v] / (g[u] * H[u, v]))
                stack.append(v)
    resid = np.abs(Hp - (g[:, None] * H * g.conj()[None, :])).max()
    if resid > 1e-8:
        return None, resid
    P = np.zeros((N, N)); P[perm, np.arange(N)] = 1
    U = np.diag(g.conj()) @ P
    U3 = np.linalg.matrix_power(U, 3)
    c = U3[0, 0]
    assert np.allclose(U3, c * np.eye(N))
    U = U * np.exp(-1j * np.angle(c) / 3)
    return U, resid


# ---------------------------------------------------------------- (a) Wilson-line compatibility (no flux)
L = 12
comp = {}
for f1, f2 in itertools.product([0, 1 / 3, 2 / 3, 1 / 2], repeat=2):
    H, links = hopping(L, 0, 2 * np.pi * f1, 2 * np.pi * f2)
    U, res = z3_operator(H, L, links)
    comp[f"{f1:.3f},{f2:.3f}"] = U is not None
R["a_wilson_compatible_L12"] = comp
allowed = [k for k, v in comp.items() if v]
print("allowed Wilson lines (phi1, phi2)/2pi:", allowed, flush=True)

# ---------------------------------------------------------------- (b) magnetised T2/Z3
L = 36
sgn_ok = None
for sgn in (1, -1):
    H, links = hopping(L, 3, 0, 0, sgn)
    ph = flux_check(H, L)
    if np.allclose(ph, np.pi * 3 / L ** 2, atol=1e-9):
        sgn_ok = sgn
print("gauge sign", sgn_ok, flush=True)
fixed = hex_fixed_points(L)
res_b = {}
three = []
for M in range(1, 13):
    for wl in [(0, 0), (1 / 3, 2 / 3), (2 / 3, 1 / 3), (1 / 3, 1 / 3), (1 / 2, 1 / 2), (1 / 3, 0)]:
        H, links = hopping(L, M, 2 * np.pi * wl[0], 2 * np.pi * wl[1], sgn_ok)
        U, rs = z3_operator(H, L, links)
        key = f"M={M} W=({wl[0]:.3f},{wl[1]:.3f})"
        if U is None:
            res_b[key] = "no Z3 operator (incompatible)"
            continue
        w, v = np.linalg.eigh(H)
        V = v[:, :M]
        gap = float(w[M] - w[M - 1]); spread = float(w[M - 1] - w[0])
        ev = np.linalg.eigvals(V.conj().T @ U @ V)
        lab = np.round(np.angle(ev) / (2 * np.pi / 3)).astype(int) % 3
        n = [int(np.sum(lab == q)) for q in range(3)]
        # weight of each sector's modes at the 3 fixed points (relative to uniform density)
        wu, vu = np.linalg.eig(V.conj().T @ U @ V)
        labs = np.round(np.angle(wu) / (2 * np.pi / 3)).astype(int) % 3
        modes = V @ vu
        loc = {}
        for q in range(3):
            sel = modes[:, labs == q]
            if sel.shape[1]:
                sel = np.linalg.qr(sel)[0]
            loc[q] = [float(np.sum(np.abs(sel[fp[0] * L + fp[1]]) ** 2) * L * L / max(1, sel.shape[1]))
                      if sel.shape[1] else 0.0 for fp in fixed]
        res_b[key] = dict(counts=n, LLL_gap=gap, LLL_spread=spread, weight_at_fixed_points=loc,
                          U_residual=float(np.abs(U @ H - H @ U).max()))
        if 3 in n:
            three.append(key + f" sector {n.index(3)}")
    print(M, {k.split(' ')[1]: (v['counts'] if isinstance(v, dict) else 'X') for k, v in res_b.items() if k.startswith(f"M={M} ")},
          f"{time.time()-t0:.1f}s", flush=True)
R["b_magnetised"] = res_b
R["b_configs_with_exactly_3"] = three

# ---------------------------------------------------------------- (c) T6/Z3 fixed-point classes
fp = list(itertools.product(range(3), repeat=3))
cls = {}
for nwl in range(4):
    labels = [tuple(f[t] for t in range(nwl)) for f in fp]
    sizes = sorted(np.unique([str(l) for l in labels], return_counts=True)[1].tolist())
    cls[f"{nwl} independent Wilson lines"] = dict(classes=len(set(labels)), sizes=sizes)
R["c_T6_classes"] = cls
R["c_note"] = ("local twist at fixed point (n1,n2,n3) = V + sum n_t a_t; with 2 independent lines the 27 points split "
               "into 9 classes of 3 (the 3 = fixed points of the untouched T2). Which classes host massless families "
               "needs the E8xE8 shift/Wilson vectors and level matching: not computed.")
R["grading"] = dict(a_allowed=allowed,
                    a_as_expected=sorted(allowed) == sorted(["0.000,0.000", "0.333,0.667", "0.667,0.333"]),
                    b_all_LLL_equal_M=all(sum(v["counts"]) == int(k.split()[0][2:]) and v["LLL_gap"] > 5 * v["LLL_spread"]
                                          for k, v in res_b.items() if isinstance(v, dict)),
                    b_toy_reaches_3=len(three) > 0,
                    c_classes=[cls[k]["sizes"] for k in cls])
R["seconds"] = time.time() - t0
print(json.dumps(R["grading"], indent=1), three[:20])
json.dump(R, open(f"{OUT}/orb_o3.json", "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))

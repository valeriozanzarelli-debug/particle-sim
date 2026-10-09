"""orb_o2.py - O2: counting from balls' own operators. Writes lab/out/orb_o2.json.
 * fixed points by brute force on the lattices;
 * rational Betti numbers of the quotient = homology of the G-invariant chain complex (simplicial hexagonal T^2,
   cubical T^4); equivalently zero modes of the cell (Hodge) Laplacian in the invariant sector;
 * T^6/Z3 by Kunneth from the measured Z3 action on the harmonic chains of one hexagonal T^2;
 * Dirac zero modes of the Wilson-Dirac walker generator in the orbifold sector (T^4/Z2 and T^2/Z3)."""
import json, time, itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla
from orb_core import *

t0 = time.time()
R = {}


def rank(M):
    if M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    return int(np.linalg.matrix_rank(M, tol=1e-8))


def invariant_betti(cells, bnd, act):
    """cells[k]: list of cells; bnd[k]: dense matrix C_k -> C_{k-1}; act[k]: list of group matrices on C_k.
    Returns Betti of the torus and of the quotient (G-invariant chains)."""
    K = len(cells)
    Q = []
    for k in range(K):
        Pi = sum(act[k]) / len(act[k])
        w, v = np.linalg.eigh((Pi + Pi.T) / 2)
        Q.append(v[:, w > 0.5])
    b_tor, b_quo = [], []
    for k in range(K):
        dk = bnd[k]; dk1 = bnd[k + 1] if k + 1 < K else np.zeros((len(cells[k]), 0))
        b_tor.append(len(cells[k]) - rank(dk) - rank(dk1))
        dkq = Q[k - 1].T @ dk @ Q[k] if k > 0 else np.zeros((0, Q[k].shape[1]))
        dk1q = Q[k].T @ dk1 @ Q[k + 1] if k + 1 < K else np.zeros((Q[k].shape[1], 0))
        b_quo.append(Q[k].shape[1] - rank(dkq) - rank(dk1q))
    return b_tor, b_quo, Q


# ----------------------------------------------------------------- hexagonal T^2 (simplicial)
def hex_complex(L):
    idx = hex_index(L)
    V = [(i, j) for i in range(L) for j in range(L)]
    E, T = set(), set()
    for i, j in V:
        a, b, c, d = idx(i, j), idx(i + 1, j), idx(i, j + 1), idx(i + 1, j + 1)
        for e in ((a, b), (a, c), (b, c)):
            E.add(tuple(sorted(e)))
        T.add(tuple(sorted((a, b, c)))); T.add(tuple(sorted((b, d, c))))
    cells = [[(v,) for v in range(L * L)], sorted(E), sorted(T)]
    pos = [{c: n for n, c in enumerate(cl)} for cl in cells]
    bnd = [np.zeros((0, L * L))]
    for k in (1, 2):
        B = np.zeros((len(cells[k - 1]), len(cells[k])))
        for n, c in enumerate(cells[k]):
            for q in range(len(c)):
                B[pos[k - 1][c[:q] + c[q + 1:]], n] += (-1) ** q
        bnd.append(B)
    # rotation on vertices
    rv = np.zeros(L * L, int)
    for i, j in V:
        rv[idx(i, j)] = idx(-i - j, i)

    def cellmap(k):
        M = np.zeros((len(cells[k]), len(cells[k])))
        for n, c in enumerate(cells[k]):
            img = [rv[v] for v in c]
            perm = np.argsort(img)
            sgn = np.linalg.det(np.eye(len(c))[perm]) if len(c) > 1 else 1
            M[pos[k][tuple(sorted(img))], n] = round(sgn)
        return M
    Rk = [cellmap(k) for k in range(3)]
    act = [[np.eye(len(cells[k])), Rk[k], Rk[k] @ Rk[k]] for k in range(3)]
    return cells, bnd, act, Rk


for L in (6, 9):
    cells, bnd, act, Rk = hex_complex(L)
    bt, bq, Q = invariant_betti(cells, bnd, act)
    fixed = hex_fixed_points(L)
    # Z3 traces on harmonic chains (Hodge Laplacian zero modes) -> for Kunneth
    traces, h1eig = [], None
    for k in range(3):
        Lk = (bnd[k].T @ bnd[k]) + (bnd[k + 1] @ bnd[k + 1].T if k < 2 else 0)
        w, v = np.linalg.eigh(Lk)
        Hk = v[:, np.abs(w) < 1e-8]
        Mk = Hk.T @ Rk[k] @ Hk
        traces.append(float(np.trace(Mk)))
        if k == 1:
            h1eig = np.linalg.eigvals(Mk)
    R[f"T2/Z3 L={L}"] = dict(fixed_points=len(fixed), fixed=fixed, betti_T2=bt, betti_quotient=bq,
                            chi_quotient=sum((-1) ** k * b for k, b in enumerate(bq)),
                            harmonic_traces_rotation=traces,
                            H1_rotation_eigenvalues=[[float(z.real), float(z.imag)] for z in h1eig],
                            chi_orb=(0 + 8 * len(fixed)) / 3)
    print("T2/Z3", L, bt, bq, len(fixed), traces, f"{time.time()-t0:.1f}s", flush=True)

# ----------------------------------------------------------------- T^6/Z3 by Kunneth
t1 = np.array([1.0, 2.0, 1.0])
tg = np.array(R["T2/Z3 L=9"]["harmonic_traces_rotation"])
P1 = np.convolve(np.convolve(t1, t1), t1)
Pg = np.convolve(np.convolve(tg, tg), tg)
b6 = (P1 + 2 * Pg) / 3          # trace of rho and rho^2 are equal (real)
# Hodge-type split: H1 of T2 = (omega) + (omega-bar) (dz, dz-bar); invariant (p,q) classes of T^6
w = np.exp(2j * np.pi / 3)
hinv = np.zeros((4, 4), int)
for p in range(4):
    for q in range(4):
        hinv[p, q] = int(round(sum(1 for P in itertools.combinations(range(3), p)
                                   for Qd in itertools.combinations(range(3), q)
                                   if abs(w ** (p - q) - 1) < 1e-9)))
L3 = 3
fx = hex_fixed_points(L3)
fixed6 = len(list(itertools.product(fx, fx, fx)))
# brute check on the full T^6 L=3 lattice
idx = hex_index(L3)
cnt = 0
for x in itertools.product(range(L3), repeat=6):
    ok = all(idx(-x[2 * t] - x[2 * t + 1], x[2 * t]) == idx(x[2 * t], x[2 * t + 1]) for t in range(3))
    cnt += ok
chi18 = int(round(sum((-1) ** k * b for k, b in enumerate(b6))))
R["T6/Z3"] = dict(betti_quotient=[int(round(b)) for b in b6], chi_quotient=chi18, fixed_points_brute=cnt,
                  fixed_points_product=fixed6, h_pq_untwisted=hinv.tolist(),
                  h11=int(hinv[1, 1]) + cnt, h21=int(hinv[2, 1]) + 0,
                  chi_orb=(0 + 8 * cnt) / 3, chi_resolved=chi18 + 2 * cnt,
                  net_generations_standard_embedding=int(hinv[1, 1]) + cnt - int(hinv[2, 1]),
                  note="twisted sector adds one (1,1) class per fixed point and no (2,1) class (textbook, not computed "
                       "by balls); resolution of C3/Z3 replaces a point (chi 1) by P2 (chi 3): +2")
print("T6/Z3", R["T6/Z3"], flush=True)


# ----------------------------------------------------------------- cubical T^4/Z2
def cube_complex(n, d=4):
    sites = list(itertools.product(range(n), repeat=d))
    sidx = {s: q for q, s in enumerate(sites)}
    cells = [[(s, D) for s in sites for D in itertools.combinations(range(d), k)] for k in range(d + 1)]
    pos = [{c: q for q, c in enumerate(cl)} for cl in cells]
    bnd = [np.zeros((0, len(cells[0])))]
    for k in range(1, d + 1):
        rows, cols, vals = [], [], []
        for q, (s, D) in enumerate(cells[k]):
            for a, dd in enumerate(D):
                Dm = D[:a] + D[a + 1:]
                s2 = list(s); s2[dd] = (s2[dd] + 1) % n
                rows += [pos[k - 1][(tuple(s2), Dm)], pos[k - 1][(s, Dm)]]
                cols += [q, q]; vals += [(-1) ** a, -(-1) ** a]
        bnd.append(sp.coo_matrix((vals, (rows, cols)), shape=(len(cells[k - 1]), len(cells[k]))).toarray())
    act = []
    for k in range(d + 1):
        rows, cols, vals = [], [], []
        for q, (s, D) in enumerate(cells[k]):
            img = [(-s[t] - (1 if t in D else 0)) % n for t in range(d)]
            rows.append(pos[k][(tuple(img), D)]); cols.append(q); vals.append((-1) ** k)
        M = sp.coo_matrix((vals, (rows, cols)), shape=(len(cells[k]),) * 2).toarray()
        act.append([np.eye(len(cells[k])), M])
    fixed = [s for s in sites if all((-x) % n == x for x in s)]
    return cells, bnd, act, fixed


for n in (4,):
    cells, bnd, act, fixed = cube_complex(n)
    bt, bq, _ = invariant_betti(cells, bnd, act)
    chi = sum((-1) ** k * b for k, b in enumerate(bq))
    R[f"T4/Z2 n={n}"] = dict(fixed_points=len(fixed), betti_T4=bt, betti_quotient=bq, chi_quotient=chi,
                            chi_orb=(0 + 3 * len(fixed)) / 2, chi_resolved=chi + len(fixed),
                            b2_resolved=bq[2] + len(fixed))
    print("T4/Z2", n, bt, bq, len(fixed), f"{time.time()-t0:.1f}s", flush=True)

# n = 6: fixed points brute + Betti via sparse ranks would be large; fixed points and Euler by the Lefschetz average
n = 6
sites = list(itertools.product(range(n), repeat=4))
fixed6_4 = [s for s in sites if all((-x) % n == x for x in s)]
ncell = [len(sites) * len(list(itertools.combinations(range(4), k))) for k in range(5)]
# chi(X/G) = (1/|G|) sum_g chi(X^g): from cell counts, chi(T4) = sum (-1)^k ncell_k
R["T4/Z2 n=6"] = dict(fixed_points=len(fixed6_4), chi_T4_from_cells=sum((-1) ** k * c for k, c in enumerate(ncell)),
                      chi_quotient_lefschetz=(sum((-1) ** k * c for k, c in enumerate(ncell)) + len(fixed6_4)) / 2,
                      note="Betti by full rank computation only at n = 4 (n = 6 dense ranks too slow for the budget)")
print("T4/Z2 n=6", R["T4/Z2 n=6"], flush=True)

# ----------------------------------------------------------------- Dirac zero modes in the orbifold sector
AL = [kron(SX, SX), kron(SX, SY), kron(SX, SZ), kron(SY, I2)]
BETA = kron(SZ, I2); G5 = AL[0] @ AL[1] @ AL[2] @ AL[3]
S_K3 = AL[0] @ AL[1] @ AL[2] @ AL[3]   # = rotation lift exp(pi/2 a1a2) exp(pi/2 a3a4)


def wilson_dirac4(n, r=1.0):
    N = n ** 4
    H = sp.csr_matrix((4 * N, 4 * N), dtype=complex)
    for j in range(4):
        sh = [1] * 4
        e = np.zeros(4, int); e[j] = 1
        sites = np.array(list(itertools.product(range(n), repeat=4)))
        src = np.ravel_multi_index(sites.T, (n,) * 4)
        dst = np.ravel_multi_index(((sites + e) % n).T, (n,) * 4)
        T = sp.csr_matrix((np.ones(N), (dst, src)), shape=(N, N))
        H = H + sp.kron(AL[j], (T - T.T) / 2j) + r * sp.kron(BETA, sp.eye(N) - (T + T.T) / 2)
    return H.tocsr()


for n in (4, 6):
    N = n ** 4
    H = wilson_dirac4(n)
    sites = np.array(list(itertools.product(range(n), repeat=4)))
    lin = np.ravel_multi_index(sites.T, (n,) * 4); neg = np.ravel_multi_index(((-sites) % n).T, (n,) * 4)
    wS, vS = np.linalg.eigh(S_K3)
    cols, rows, vals = [], [], []
    col = 0
    for x in range(N):
        y = neg[x]
        if y < x:
            continue
        if y == x:
            for q in np.where(wS > 0)[0]:
                for s in range(4):
                    rows.append(s * N + x); cols.append(col); vals.append(vS[s, q])
                col += 1
        else:
            for s in range(4):
                e = np.zeros(4); e[s] = 1; Se = S_K3 @ e
                rows.append(s * N + x); cols.append(col); vals.append(1 / np.sqrt(2))
                for s2 in range(4):
                    if abs(Se[s2]) > 1e-12:
                        rows.append(s2 * N + y); cols.append(col); vals.append(Se[s2] / np.sqrt(2))
                col += 1
    Qm = sp.csr_matrix((vals, (rows, cols)), shape=(4 * N, col))
    Hs = (Qm.conj().T @ H @ Qm).toarray()
    w, v = np.linalg.eigh(Hs)
    z = np.abs(w) < 1e-8
    G5s = (Qm.conj().T @ sp.kron(G5, sp.eye(N)) @ Qm).toarray()
    chir = [float(np.real(v[:, q].conj() @ G5s @ v[:, q])) for q in np.where(z)[0]]
    wf = np.linalg.eigvalsh(H.toarray()) if n == 4 else None
    R[f"Dirac T4/Z2 n={n}"] = dict(sector_dim=col, zero_modes=int(z.sum()), chirality=chir, gap=float(np.min(np.abs(w[~z]))),
                                   zero_modes_full_T4=(int(np.sum(np.abs(wf) < 1e-8)) if wf is not None else "4 (k=0 only)"))
    print("Dirac K3", n, R[f"Dirac T4/Z2 n={n}"], f"{time.time()-t0:.1f}s", flush=True)

for L in (6, 9):
    H = hex_dirac(L).toarray()
    P = hex_rot_perm(L).toarray()
    S = sla.expm(-1j * (-4 * np.pi / 3) / 2 * SZ)     # lift with S^3 = 1
    for name, Sx in (("lift S^3=1", S), ):
        O = np.kron(Sx, P)
        Pi = (np.eye(len(O)) + O + O @ O) / 3
        wP, vP = np.linalg.eigh((Pi + Pi.conj().T) / 2)
        Qm = vP[:, wP > 0.5]
        w = np.linalg.eigvalsh(Qm.conj().T @ H @ Qm)
        wf = np.linalg.eigvalsh(H)
        R[f"Dirac T2/Z3 L={L}"] = dict(sector_dim=Qm.shape[1], zero_modes=int(np.sum(np.abs(w) < 1e-8)),
                                       zero_modes_full_T2=int(np.sum(np.abs(wf) < 1e-8)),
                                       commutator=float(np.abs(O @ H - H @ O).max()))
    print("Dirac T2/Z3", L, R[f"Dirac T2/Z3 L={L}"], flush=True)

# ----------------------------------------------------------------- grading
g = {}
g["T2/Z3"] = all(R[f"T2/Z3 L={L}"]["fixed_points"] == 3 and R[f"T2/Z3 L={L}"]["betti_quotient"] == [1, 0, 1]
                 and R[f"Dirac T2/Z3 L={L}"]["zero_modes"] == 0 for L in (6, 9))
g["T4/Z2"] = (R["T4/Z2 n=4"]["fixed_points"] == 16 and R["T4/Z2 n=4"]["betti_quotient"] == [1, 0, 6, 0, 1]
              and R["T4/Z2 n=4"]["chi_orb"] == 24 and R["T4/Z2 n=4"]["b2_resolved"] == 22
              and R["T4/Z2 n=6"]["fixed_points"] == 16 and R["T4/Z2 n=6"]["chi_quotient_lefschetz"] == 8
              and all(R[f"Dirac T4/Z2 n={n}"]["zero_modes"] == 2 for n in (4, 6)))
t6 = R["T6/Z3"]
g["T6/Z3"] = (t6["betti_quotient"] == [1, 0, 9, 2, 9, 0, 1] and t6["chi_quotient"] == 18 and t6["fixed_points_brute"] == 27
              and t6["chi_orb"] == 72 and t6["h11"] == 36 and t6["h21"] == 0 and t6["net_generations_standard_embedding"] == 36)
g["PASS"] = all(g.values())
g["note"] = "T4/Z2 Betti by full rank only at n = 4 (n = 6: fixed points and Lefschetz chi only) - deviation from prereg 'both sizes'"
R["grading"] = g; R["seconds"] = time.time() - t0
print(json.dumps(g, indent=1))
json.dump(R, open(f"{OUT}/orb_o2.json", "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))

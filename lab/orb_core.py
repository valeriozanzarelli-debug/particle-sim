"""orb_core.py - small helpers for the orbifold pre-test (orb_o1/o2/o3).
Lattices: hexagonal (triangular) L x L torus in skew coordinates r = i a1 + j a2, a1 = (1,0), a2 = (1/2, sqrt3/2);
120 deg rotation R: (i, j) -> (-i-j, i)  (a1 -> a2 - a1, a2 -> -a1). Bonds e1 = (1,0), e2 = (-1,1), e3 = (0,-1)
(skew), R e_k = e_{k+1}. Hypercubic n^d torus with inversion x -> -x.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

SX = np.array([[0, 1], [1, 0]], complex); SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2, dtype=complex)
HEX_E = [(1, 0), (-1, 1), (0, -1)]
HEX_N = [np.array([np.cos(2 * np.pi * k / 3), np.sin(2 * np.pi * k / 3)]) for k in range(3)]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def kron(*a):
    r = a[0]
    for b in a[1:]:
        r = np.kron(r, b)
    return r


def hex_index(L):
    return lambda i, j: (i % L) * L + (j % L)


def hex_shift(L, e):
    """sparse T_e: (T psi)(x + e) = psi(x)"""
    idx = hex_index(L)
    rows, cols = [], []
    for i in range(L):
        for j in range(L):
            rows.append(idx(i + e[0], j + e[1])); cols.append(idx(i, j))
    return sp.csr_matrix((np.ones(L * L), (rows, cols)), shape=(L * L, L * L))


def hex_rot_perm(L, center=(0, 0)):
    """permutation matrix P: (P psi)(R x) = psi(x), R = 120 deg rotation about the lattice point `center`"""
    idx = hex_index(L)
    ci, cj = center
    rows, cols = [], []
    for i in range(L):
        for j in range(L):
            a, b = i - ci, j - cj
            rows.append(idx(-a - b + ci, a + cj)); cols.append(idx(i, j))
    return sp.csr_matrix((np.ones(L * L), (rows, cols)), shape=(L * L, L * L))


def hex_fixed_points(L):
    idx = hex_index(L)
    return [(i, j) for i in range(L) for j in range(L) if idx(-i - j, i) == idx(i, j)]


def hex_dirac(L, r=1.0):
    """2-component hexagonal Wilson-Dirac generator of the walker: H = sum_k (n_k.sigma)(T_k - T_k^dag)/(2i)*(2/3)
    + r sigma_z sum_k (1 - (T_k + T_k^dag)/2). Massless; zero modes only at k = 0 (doublers lifted by r)."""
    N = L * L
    H = sp.csr_matrix((2 * N, 2 * N), dtype=complex)
    for e, n in zip(HEX_E, HEX_N):
        T = hex_shift(L, e)
        ns = n[0] * SX + n[1] * SY
        H = H + sp.kron(ns, (T - T.T) / (2j) * (2 / 3)) + r * sp.kron(SZ, sp.eye(N) - (T + T.T) / 2)
    return H.tocsr()


def nullspace(A, tol=1e-9):
    u, s, vh = np.linalg.svd(A, full_matrices=False)
    k = int(np.sum(s > tol * max(1.0, s[0] if len(s) else 1)))
    return vh[k:].conj().T, s


def solve_spinor_operator(apply_W, apply_P, dim_s, nvec, rng, nrand=2):
    """all dim_s x dim_s matrices S with W (S x P) = (S x P) W, tested on random states.
    apply_W / apply_P act on arrays shaped (dim_s, nvec). Returns basis list of S."""
    cols = []
    psis = [rng.normal(size=(dim_s, nvec)) + 1j * rng.normal(size=(dim_s, nvec)) for _ in range(nrand)]
    for a in range(dim_s):
        for b in range(dim_s):
            E = np.zeros((dim_s, dim_s), complex); E[a, b] = 1
            col = []
            for psi in psis:
                Pp = apply_P(psi)
                lhs = apply_W(E @ Pp)
                rhs = E @ apply_P(apply_W(psi))
                col.append((lhs - rhs).ravel())
            cols.append(np.concatenate(col))
    A = np.array(cols).T
    ns, s = nullspace(A)
    return [ns[:, q].reshape(dim_s, dim_s) for q in range(ns.shape[1])]


def in_span(S, basis, tol=1e-8):
    if not basis:
        return False
    B = np.array([b.ravel() for b in basis]).T
    c, *_ = np.linalg.lstsq(B, S.ravel(), rcond=None)
    return np.linalg.norm(B @ c - S.ravel()) < tol * np.linalg.norm(S)


def n_invariant(S, tol=1e-8):
    return int(np.sum(np.abs(np.linalg.eigvals(S) - 1) < tol))

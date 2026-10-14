"""ciclo_p3_core.py - Phase 3 of CICLO_LIVELLI: the on-the-fly CHECK and LOCAL REFINEMENT (step 4 of the loop).
STATUS: DEVELOPMENT (PARZIALE). No production run, no exam yet.

Setting (1-D, reuses phase 2 read-only): base level = 2187 cells of 1.6 lambda_C (the phase-2 base, validated against
the fine zig-zag in A8) with the soft-Coulomb atom; coarse level = the automatic b = 3 chain of phase 2 (B = 3, 9).

Adaptive class below = FIRST design: coarse cells carry the amplitude on their sharp block vector, refined cells their
base cells, coarse-coarse block = learned rule, the rest Galerkin. Development result (ciclo_p3_dev.py): the Galerkin
coupling of sharp blocks is ~6x too stiff (B = 9: diagonal 2.4e-3 vs learned 3.9e-4), and gluing gives spurious
states at -6.7e-4 (20x the atom's binding). Planned design (not yet coded here): smoothed-aggregation prolongator
(the automatic block vectors smoothed s = 3 times by the level's own rule, Loewdin), Galerkin rule of it agrees with the
learned rule within 0.3-1.3 % and the bound states within 0.7 % at B = 9; mixed basis = masked smooth functions +
refined base cells (orthonormal), Crank-Nicolson (exactly unitary); check = base rule on patch + buffer vs the
adaptive evolution, refine + rollback above a threshold to be pre-registered.
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ciclo_p2_core as c
import ciclo_p2_1d as d1

N0 = 2187
TH, B0, WB = d1.TH, d1.B0, d1.WB


HIER_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "ciclo_p3_hier.pkl")   # not committed


def hierarchy(nlev=2):
    """cached wrapper of build_hierarchy (deterministic: same phase-2 code and settings)"""
    import pickle
    if os.path.exists(HIER_CACHE):
        hy = pickle.load(open(HIER_CACHE, "rb"))
        if len(hy["levels"]) >= nlev:
            return hy
    hy = build_hierarchy(nlev)
    pickle.dump(hy, open(HIER_CACHE, "wb"))
    return hy


def build_hierarchy(nlev=2):
    """phase-2 b = 3 chain on the free base (kinetic rule learned), recording the accepted W and H of each level.
    Returns dict: g (base stencil), Hb (sparse base kinetic), V0 (base soft Coulomb), levels [(B, W_tot, H_c)],
    where W_tot maps base -> level (rows orthonormal, block supports)."""
    g, berr = c.make_base(TH, B0, WB)
    Hb = c.circulant(g, N0).tocsr()
    V0 = d1.base_potential(N0)[0]
    rec = []
    orig_cv, orig_fit = c.choose_variables, c.fit_generator
    state = {}

    def cv(*a, **k):
        r = orig_cv(*a, **k)
        state["W"] = r[0]
        return r

    def fit(*a, **k):
        r = orig_fit(*a, **k)
        rec.append((state["W"], r[0]))
        return r
    c.choose_variables, c.fit_generator = cv, fit
    try:
        kk = 2 * np.pi * np.fft.fftfreq(N0)
        nobs = int(np.sum(np.abs(kk) <= d1.KOBS + 1e-12))
        out = c.run_chain(3, Hb, nlev=nlev, eps_level=1e-3, wmax=d1.WMAX, nobs=nobs, V0=V0, nb=3, cmax=2, eps_c=2e-2,
                          log=lambda s: None)
    finally:
        c.choose_variables, c.fit_generator = orig_cv, orig_fit
    levels = []
    Wt = sps.identity(N0, format="csr")
    for L in out["levels"]:
        if "flag" in L:
            break
        cand = [(W, H) for W, H in rec if H.shape[0] == L["ncell"]]
        W, H = cand[-1]
        Wt = (W @ Wt).tocsr()
        levels.append(dict(B=L["B"], W=Wt, H=H.tocsr(), w=L["w"], acc_err=L["acc_err"], bound_rel_err=L["bound_rel_err"]))
    return dict(g=g, base_err=berr, Hb=Hb, V0=V0, levels=levels, chain=out)


class Adaptive:
    """mixed coarse / base representation on top of level `lev` of the hierarchy"""

    def __init__(self, hy, lev, refined=None):
        L = hy["levels"][lev]
        self.W = L["W"].tocsr()            # (nc, N0)
        self.Hc = L["H"].tocsr()
        self.Hb = hy["Hb"]
        self.V0 = hy["V0"]
        self.nc = self.W.shape[0]
        self.B = N0 // self.nc
        self.ref = np.zeros(self.nc, bool) if refined is None else refined.copy()
        self.WV0W = (self.W @ sps.diags(self.V0) @ self.W.conj().T).tocsr()
        self.build()

    def build(self):
        """mixed basis S (N0 x m): columns ordered by coarse cell; coarse cell -> 1 column, refined -> B columns"""
        cols, self.slot = [], []
        rows, cidx, vals = [], [], []
        m = 0
        self.cmap = np.full(self.nc, -1)      # column of the coarse amplitude (or -1)
        self.fmap = {}                        # refined cell -> its base columns
        for I in range(self.nc):
            if not self.ref[I]:
                r = self.W.getrow(I)
                rows += list(r.indices); cidx += [m] * r.nnz; vals += list(np.conj(r.data))
                self.cmap[I] = m
                m += 1
            else:
                bs = np.arange(I * self.B, (I + 1) * self.B)
                rows += list(bs); cidx += list(range(m, m + self.B)); vals += [1.0] * self.B
                self.fmap[I] = np.arange(m, m + self.B)
                m += self.B
        self.m = m
        self.S = sps.csr_matrix((vals, (rows, cidx)), shape=(N0, m))
        self.cc = np.where(self.cmap >= 0)[0]                 # coarse cells
        self.ccol = self.cmap[self.cc]
        # static part of the generator: S^dag Hb S everywhere, then the coarse-coarse block replaced by H_c
        A = (self.S.conj().T @ (self.Hb + sps.diags(self.V0)) @ self.S).tolil()
        Hcc = (self.Hc + self.WV0W)[self.cc][:, self.cc].toarray()
        A = A.toarray()
        A[np.ix_(self.ccol, self.ccol)] = Hcc
        self.H0 = sps.csr_matrix(A)
        self.Sd = self.S.conj().T.tocsr()

    def gen(self, Vd):
        """generator at time t for the extra base potential Vd (array N0, or None)"""
        if Vd is None:
            return self.H0
        return (self.H0 + self.Sd @ sps.diags(Vd) @ self.S).tocsr()

    def to_base(self, a):
        return self.S @ a

    def from_base(self, psi):
        return self.Sd @ psi

    def refine(self, cells, a):
        """refine the listed coarse cells (lossless); returns the state in the new basis"""
        psi = self.to_base(a)
        self.ref[np.asarray(cells, int)] = True
        self.build()
        return self.from_base(psi)


def cn_step(H, a, dt):
    """Crank-Nicolson (Cayley) step: exactly unitary for Hermitian H"""
    n = H.shape[0]
    I = sps.identity(n, format="csc", dtype=complex)
    M1 = (I + 0.5j * dt * H).tocsc()
    M2 = (I - 0.5j * dt * H)
    return spla.spsolve(M1, M2 @ a)


def cn_dense(Hd, a, dt):
    n = Hd.shape[0]
    I = np.eye(n)
    return np.linalg.solve(I + 0.5j * dt * Hd, (I - 0.5j * dt * Hd) @ a)

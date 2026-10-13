"""aggmol_core.py - "Livello delle due atomi" (study aggmol_): the machine for the two-atom scale.

Builds on (read-only imports, nothing of theirs is modified):
  lab/agg_zz_core.py   the agglomerator (Level, GenRule, kcheck, sample_state_k, ...)
  lab/aggfull_core.py  the learned cell-field coupling machinery and the Galerkin field (8 "lanterns" per cell)
  lab/aggfull_l1.py    the Taylor-pinned + parity learner (fit mode "T") and the Markov learner (mode "M")
  lab/aggfull_atom.py  kinetic(), LOBPCG solve()
Learned pieces taken as they are: lab/out/aggfull_relearnT.pkl (kinetic T at a0/16 and a0/8), lab/out/aggfull_learn.pkl
(coupling at a0/16 and a0/8). New levels of this study: lab/out/aggmol_levels.pkl (a0/4, a0/2: kinetic + coupling).

New here:
  * MultiField: the coarse field of SEVERAL nuclei = sum of the Galerkin fields of each (the input is only -Z alpha/r for
    each nucleus and its position; positions anywhere, not snapped to the grid).
  * Mol: the coarse one-electron Hamiltonian H = T(k) + gamma0 V_c with exact parity sectors (sx, sy, sz) about the origin
    (the midpoint of the two protons), real arithmetic.
Units hbar = c = m_e = 1 (length lambda_C, energy m c^2)."""
import os, sys, time, pickle, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import scipy.sparse.linalg as spl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_fz_core as ac
from aggmol_fz_core import fftn, ifftn

OUT = ac.OUT
A0 = ac.A0
MC2 = ac.MC2
EH = ac.zz.ALPHA ** 2 * MC2          # Hartree in eV (m_e, no reduced mass)
ALPHA = ac.zz.ALPHA
TH = ac.TH_A
DEV = os.environ.get("AGGMOL_DEV") == "1"
LEVELS_PKL = os.path.join(OUT, ("dev_" if DEV else "") + "aggmol_levels.pkl")


# ------------------------------------------------------------------------------------------------ levels
def level_rec(B):
    """free generator record (K, offs, phi0, C) of level B (cell B * TH): 256/512 from aggfull (T), 1024/2048 from here"""
    if B in (256, 512):
        return [l for l in pickle.load(open(os.path.join(OUT, "aggfull_relearnT.pkl"), "rb")) if l["B"] == B][0]
    st = pickle.load(open(LEVELS_PKL, "rb"))
    return [l for l in st["free"] if l["B"] == B][0]


def coupling_rec(B):
    if B in (256, 512):
        s = [l for l in pickle.load(open(os.path.join(OUT, "aggfull_learn.pkl"), "rb")) if l["B"] == B][0]
    else:
        s = [l for l in pickle.load(open(LEVELS_PKL, "rb"))["coupling"] if l["B"] == B][0]
    C, offs = s["C"], s["offs"]
    g0 = complex(np.trace(ac.eval_coupling(C, offs, np.zeros((1, 3)), np.zeros((1, 3)))[0]) / 2)
    return C, offs, g0


def kinetic_grid(rec, N, h):
    """spin-blind, parity-even learned generator per unit time on the N^3 grid (as aggfull_atom.kinetic)"""
    r = ac.genrule(rec)
    Gk = r.gen_grid(N)
    g = np.real(0.5 * (Gk[..., 0, 0] + Gk[..., 1, 1]))
    del Gk
    g0 = float(np.real(0.5 * np.trace(r.gen_at(np.zeros((1, 3)))[0])))
    g = g - g0
    gm = np.roll(np.flip(g, axis=(0, 1, 2)), 1, axis=(0, 1, 2))
    return 0.5 * (g + gm) / h


def kinetic_exact(N, h):
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij", sparse=True)
    return 0.5 * (KX ** 2 + KY ** 2 + KZ ** 2)


# ------------------------------------------------------------------------------------------------ field
class MultiField:
    """sum of Galerkin fields (aggfull_core.GalerkinField) of several nuclei; centres in CELLS relative to the grid
    point N//2 (the origin)"""

    def __init__(self, N, h, centres, Z=None, dtype=np.float64):
        Z = Z if Z is not None else [1.0] * len(centres)
        self.N, self.h = N, h
        self.subs = [(sx, sy, sz) for sx in (0, 1) for sy in (0, 1) for sz in (0, 1)]
        self.V = None
        for c, z in zip(centres, Z):
            f = ac.GalerkinField(N, h, Z=z, centre=tuple(c))
            if self.V is None:
                self.V = {s: f.V[s].astype(dtype) for s in self.subs}
            else:
                for s in self.subs:
                    self.V[s] += f.V[s].astype(dtype)
            self.x = f.x
            del f
        self.k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
        self.ph1 = np.exp(1j * self.k1 * h / 2)

    apply_k = ac.GalerkinField.apply_k

    def apply(self, phi):
        return ifftn(self.apply_k(fftn(phi)))


# ------------------------------------------------------------------------------------------------ Hamiltonian
class Mol:
    """H = T(k) + gamma0 V_c on an N^3 periodic grid; exact parity sector (sx, sy, sz) in {+1, -1, 0 (none)}^3 about
    the origin; real arithmetic (T even, V real-symmetric)"""

    def __init__(self, T, field, gamma0=1.0, sym=(1, 1, 1), inv=0):
        """sym[ax] = +-1: reflection about the grid point at the origin; +-2: about the half-cell plane (index N/2 - 1/2,
        for nuclei shifted by half a cell); 0: no symmetry along that axis. inv = +-1: inversion through the origin."""
        self.T, self.f, self.g0, self.sym, self.inv = T, field, gamma0, sym, inv
        self.N = T.shape[0]
        self.pre = 1.0 / (T - T.min() + 0.5 * ALPHA ** 2)
        self.napply = 0

    def project(self, v):
        for ax, s in enumerate(self.sym):
            if abs(s) == 1:
                v = 0.5 * (v + s * np.roll(np.flip(v, axis=ax), 1, axis=ax))
            elif abs(s) == 2:
                v = 0.5 * (v + np.sign(s) * np.flip(v, axis=ax))
        if self.inv:
            v = 0.5 * (v + self.inv * np.roll(np.flip(v, axis=(0, 1, 2)), 1, axis=(0, 1, 2)))
        return v

    def H(self, v):
        v = self.project(v)
        F = fftn(v)
        out = ifftn(self.T * F + self.g0 * self.f.apply_k(F)).real
        return self.project(out)

    def _apply(self, X, fn):
        X = np.asarray(X)
        one = X.ndim == 1
        X = X.reshape(self.N ** 3, -1)
        out = np.empty_like(X)
        for c in range(X.shape[1]):
            out[:, c] = fn(X[:, c].reshape((self.N,) * 3)).ravel()
        return out[:, 0] if one else out

    def A(self, X):
        self.napply += 1
        return self._apply(X, self.H)

    def P(self, X):
        def pre(v):
            return self.project(ifftn(self.pre * fftn(self.project(v))).real)
        return self._apply(X, pre)

    def energy(self, v):
        f = v.reshape((self.N,) * 3)
        Hf = self.H(f)
        nrm = np.vdot(f, f).real
        E = np.vdot(f, Hf).real / nrm
        return float(E), float(np.linalg.norm(Hf - E * f) / np.sqrt(nrm))


def solve(mol, X0, Y=None, tol=1e-10, maxiter=40):
    """LOBPCG for the lowest len(X0.T) states of mol (block), orthogonal to Y"""
    N3 = mol.N ** 3
    Aop = spl.LinearOperator((N3, N3), matvec=mol.A, matmat=mol.A, dtype=float)
    Pop = spl.LinearOperator((N3, N3), matvec=mol.P, matmat=mol.P, dtype=float)
    X0 = np.asarray(X0, float)
    if X0.ndim == 1:
        X0 = X0[:, None]
    X0 = mol._apply(X0, mol.project)
    if Y is not None:
        X0 = X0 - Y @ (Y.T @ X0)
    w, v, hist = spl.lobpcg(Aop, X0, M=Pop, Y=Y, tol=tol, maxiter=maxiter, largest=False, retResidualNormsHistory=True)
    return w, v, len(hist)


def solve_until(mol, X0, budget_s, rtol=2e-8, Y=None, chunk=5, log=print, tag=""):
    """repeated short LOBPCG runs (restartable): returns vectors, energies, residuals, converged flag"""
    t0 = time.time()
    X = X0
    its = 0
    conv = False
    while time.time() - t0 < budget_s:
        w, X, it = solve(mol, X, Y=Y, maxiter=chunk)
        its += it
        res = [mol.energy(X[:, i]) for i in range(X.shape[1])]
        log(f"  {tag} it {its}: " + ", ".join(f"E={e * MC2:.8f} eV r={r:.1e}" for e, r in res) + f" ({time.time() - t0:.0f} s)")
        if max(r for _, r in res) < rtol:
            conv = True
            break
    return X, [e for e, _ in res], [r for _, r in res], conv, its


def grid_x(N, h):
    return (np.arange(N) - N // 2) * h


def guess_mol(N, h, centres_cells, kind="1sg", a=A0):
    """starting vectors (sums/differences of hydrogen-like orbitals on the nuclei): no shape is imposed on the result,
    LOBPCG relaxes them; only the parity sector is exact"""
    x = grid_x(N, h)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij", sparse=True)
    out = 0
    for i, c in enumerate(centres_cells):
        c = np.array(c) * h
        r = np.sqrt((X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2)
        s = 1.0 if (i == 0 or kind in ("1sg", "2sg", "1pu", "1s", "2s", "2p")) else -1.0
        g = {"1sg": np.exp(-r / a), "1su": np.exp(-r / a), "1s": np.exp(-r / a),
             "2sg": (1 - r / (2 * a)) * np.exp(-r / (2 * a)), "2s": (1 - r / (2 * a)) * np.exp(-r / (2 * a)),
             "1pu": (X - c[0]) * np.exp(-r / (2 * a)), "2p": (Z - c[2]) * np.exp(-r / (2 * a))}[kind]
        out = out + s * g
    out = np.broadcast_to(out, (N, N, N)).copy()
    return (out / np.linalg.norm(out)).ravel()

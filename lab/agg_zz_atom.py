"""agg_zz_atom.py - A2: hydrogen at the Bohr scale with the LEARNED coarse zig-zag rule (chain A).
W = e^{-i V tau/2} K e^{-i V tau/2}: K = learned kernel of the deepest level with cell <= a0/16 (1s) or <= a0/8 (n = 2),
V = cell average of -alpha/r (nucleus at a cell corner), periodic box. The learned rest phase is divided out (it is a pure
phase at this level). Levels by preconditioned inverse iteration (LOBPCG) on the Hermitian A = (W^dag - W)/2i, whose
eigenvalues are sin(E tau) (same eigenvectors as W); the energy is then read from the rule itself: E = -arg<x|W|x>/tau.
Checkpointed in lab/out/agg_zz_atom.json (one state per call).
Run: python lab/agg_zz_atom.py"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import json, pickle, sys, time
import numpy as np
import scipy.sparse.linalg as spl
from agg_zz_core import *
from agg_zz_core import GenRule

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
MC2 = zz.MC2_EV
RY = 0.5 * zz.ALPHA ** 2 * MC2                 # 13.605693 eV (Schroedinger, infinite nuclear mass)
DEV = os.environ.get("AGG_DEV") == "1"
JS = os.path.join(OUT, "dev_agg_atom.json" if DEV else "agg_zz_atom.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
t_start = time.time()


def cell_avg_inv_r(N, h, centre="corner", nsub=24, near=4):
    """<1/r> over each cubic cell; nucleus at a cell corner (default) or at a cell centre"""
    off = 0.5 if centre == "corner" else 0.0
    x = (np.arange(N) - N // 2 + off) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij", sparse=True)
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    with np.errstate(divide="ignore"):
        inv = 1 / R
    u = (np.arange(nsub) + 0.5) / nsub - 0.5
    U = np.stack(np.meshgrid(u, u, u, indexing="ij"), -1).reshape(-1, 3)
    idx = np.nonzero(np.abs(x) <= near * h)[0]
    for i in idx:
        for j in idx:
            for k in idx:
                p = (np.array([x[i], x[j], x[k]]) + U * h)
                inv[i, j, k] = np.mean(1 / np.linalg.norm(p, axis=1))
    return inv, x


class Atom:
    """the learned coarse rule with a Coulomb phase. The level stores a learned generator G (Hermitian stencil):
    one step = e^{i phi0} e^{-i G}. Eigenvectors from the Hermitian H = G/tau + V (preconditioned inverse iteration,
    LOBPCG); the energy is then read from the RULE itself, W = e^{-iV tau/2} e^{-iG} e^{-iV tau/2} (rest phase divided
    out): E = -arg<x|W|x>/tau, with the residual |Wx - lambda x| reported."""

    def __init__(self, lvrec, TH, N, h_cells=None, centre="corner"):
        self.N = N
        self.tau = lvrec["B"]
        self.h = lvrec["B"] * TH
        r = GenRule(lvrec["K"], lvrec["offs"], 3, lvrec["C"], lvrec["phi0"])
        Gk = r.gen_grid(N)                                         # (N, N, N, 2, 2)
        self.spin_part = float(np.max(np.abs(Gk[..., 0, 1])) + np.max(np.abs(Gk[..., 0, 0] - Gk[..., 1, 1])))
        self.G = np.real(0.5 * (Gk[..., 0, 0] + Gk[..., 1, 1]))    # spin-blind part (the spin part is reported)
        del Gk
        self.rest = lvrec["phi0"]
        self.Mr = np.exp(-1j * self.G)
        inv, self.x = cell_avg_inv_r(N, self.h, centre)
        self.V = -zz.ALPHA * inv
        self.tau_phys = self.tau * TH
        self.D = np.exp(-0.5j * self.V * self.tau_phys)
        self.eps_max = float(self.G.max())
        self.g_min = float(self.G.min())
        self.pre = 1.0 / (self.G / self.tau_phys + 0.5 * zz.ALPHA ** 2)
        self.napply = 0

    def W(self, f):
        return self.D * ifftn(self.Mr * fftn(self.D * f, 3), 3)

    def H(self, f):
        return ifftn(self.G / self.tau_phys * fftn(f, 3), 3) + self.V * f

    def A(self, X):
        X = np.asarray(X)
        one = X.ndim == 1
        X = X.reshape(self.N ** 3, -1)
        out = np.empty_like(X, dtype=complex)
        for c in range(X.shape[1]):
            out[:, c] = self.H(X[:, c].reshape((self.N,) * 3)).ravel()
            self.napply += 1
        return out[:, 0] if one else out

    def P(self, X):
        X = np.asarray(X)
        one = X.ndim == 1
        X = X.reshape(self.N ** 3, -1)
        out = np.empty_like(X, dtype=complex)
        for c in range(X.shape[1]):
            out[:, c] = ifftn(self.pre * fftn(X[:, c].reshape((self.N,) * 3), 3), 3).ravel()
        return out[:, 0] if one else out

    def guess(self, kind):
        X, Y, Z = np.meshgrid(self.x, self.x, self.x, indexing="ij")
        r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
        g = {"1s": np.exp(-r / A0), "2s": (1 - r / (2 * A0)) * np.exp(-r / (2 * A0)), "2p": Z * np.exp(-r / (2 * A0))}[kind]
        return (g / np.linalg.norm(g)).ravel().astype(complex)

    def energy(self, v):
        f = v.reshape((self.N,) * 3)
        Wf = self.W(f)
        lam = np.vdot(f, Wf) / np.vdot(f, f)
        resid = np.linalg.norm(Wf - lam * f) / np.linalg.norm(f)
        Eg = float(np.real(np.vdot(f, self.H(f)) / np.vdot(f, f)))
        return float(-np.angle(lam) / self.tau_phys), float(resid), Eg


def solve(atom, kind, Y=None, tol=1e-12, maxiter=300):
    N3 = atom.N ** 3
    Aop = spl.LinearOperator((N3, N3), matvec=atom.A, matmat=atom.A, dtype=complex)
    Pop = spl.LinearOperator((N3, N3), matvec=atom.P, matmat=atom.P, dtype=complex)
    x0 = atom.guess(kind)[:, None]
    if Y is not None:
        x0 = x0 - Y @ (Y.conj().T @ x0)
    w, v, hist = spl.lobpcg(Aop, x0, M=Pop, Y=Y, tol=tol, maxiter=maxiter, largest=False, retResidualNormsHistory=True)
    E, resid, Eg = atom.energy(v[:, 0])
    atom.E_generator = Eg
    return E, resid, v[:, 0], len(hist)


if __name__ == "__main__":
    st = pickle.load(open(os.path.join(OUT, "dev_agg_3d_A.pkl" if DEV else "agg_zz_3d_A.pkl"), "rb"))
    TH = A0 / 16 / 256
    levels = st["levels"]
    L1 = max([l for l in levels if l["B"] <= 256], key=lambda l: l["B"])
    L2 = max([l for l in levels if l["B"] <= 512], key=lambda l: l["B"])
    NN = int(os.environ.get("AGG_N", "192"))
    jobs = [("1s", L1, NN), ("1s", L2, NN), ("2s", L2, NN), ("2p", L2, NN)]
    for kind, L, N in jobs:
        key = f"{kind}_B{L['B']}_N{N}"
        if key in res:
            continue
        if time.time() - t_start > 30:
            print("next call continues", flush=True)
            break
        t0 = time.time()
        at = Atom(L, TH, N, None)
        Y = None
        if kind != "1s":
            vecs = []
            for k2 in ("1s", "2s"):
                kk = f"{k2}_B{L['B']}_N{N}"
                fn = os.path.join(OUT, ("dev_" if DEV else "") + f"agg_zz_atom_{kk}.npy")
                if k2 != kind and kk in res and os.path.exists(fn):
                    vecs.append(np.load(fn))
            Y = np.stack(vecs, 1) if vecs else None
        E, resid, v, nit = solve(at, kind, Y)
        np.save(os.path.join(OUT, ("dev_" if DEV else "") + f"agg_zz_atom_{key}.npy"), v)
        n = 1 if kind == "1s" else 2
        ex = -RY / n ** 2
        res[key] = dict(kind=kind, B=L["B"], cell_lambdaC=at.h, cell_over_a0=at.h / A0, N=N, box_a0=N * at.h / A0,
                        E_eV=E * MC2, schroedinger_eV=ex, rel_err=E * MC2 / ex - 1, resid=resid, iterations=nit,
                        rule_applications=at.napply, spin_part=at.spin_part, eps_max=at.eps_max, g_min=at.g_min,
                        rest_phase=float(at.rest), E_generator_eV=at.E_generator * MC2,
                        seconds=time.time() - t0)
        print(f"{key}: E = {E * MC2:.6f} eV vs {ex:.6f} ({100 * (E * MC2 / ex - 1):+.4f} %), residual {resid:.1e}, "
              f"{nit} it, {time.time() - t0:.0f} s, spin part {at.spin_part:.1e}", flush=True)
        json.dump(res, open(JS, "w"), indent=1)
        del at

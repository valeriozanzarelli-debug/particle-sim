"""aggmol_frat_md.py - study aggmol_frat_: classical molecular dynamics of H2 molecule-balls (step 2 and 3) and of the
learned blob-balls (the next level of the loop).
Units: length Angstrom, energy K (k_B = 1), mass amu, time t_u = sqrt(amu A^2 / k_B) = 1.0967 ps.
Integrator: BAOAB Langevin (the friction + noise pair obeys fluctuation-dissipation by construction: <xi^2> = 2 gamma T / m).
Forces: tabulated isotropic pair force (energy shifted to 0 at the cutoff), Verlet list from scipy cKDTree."""
import os, sys, json, time
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
TU_PS = np.sqrt(1.66053906660e-27 * 1e-20 / 1.380649e-23) * 1e12


class Table:
    def __init__(self, r, v, dv):
        self.r0, self.dr = r[0], r[1] - r[0]
        self.v, self.dv = v, dv
        self.n = len(r)
        self.rmax = r[-1]

    def __call__(self, r):
        x = (r - self.r0) / self.dr
        i = np.clip(x.astype(np.int64), 0, self.n - 2)
        f = np.clip(x - i, 0, 1)
        v = self.v[i] * (1 - f) + self.v[i + 1] * f
        dv = self.dv[i] * (1 - f) + self.dv[i + 1] * f
        out = r >= self.rmax
        v[out] = 0; dv[out] = 0
        below = r < self.r0                                  # steep continuation (never reached in practice)
        if below.any():
            v[below] = self.v[0] + self.dv[0] * (r[below] - self.r0); dv[below] = self.dv[0]
        return v, dv


class System:
    def __init__(self, x, m, table, skin=1.0, box=None):
        self.x = np.array(x, float); self.N = len(x)
        self.m = m
        self.tab = table; self.rc = table.rmax; self.skin = skin
        self.box = box
        self.v = np.zeros_like(self.x)
        self.build()
        self.f, self.epot = self.force()

    def build(self):
        if self.box is None:
            tr = cKDTree(self.x)
        else:
            tr = cKDTree(np.mod(self.x, self.box), boxsize=self.box)
        self.pairs = tr.query_pairs(self.rc + self.skin, output_type="ndarray")
        self.x_build = self.x.copy()

    def force(self):
        if np.max(np.sum((self.x - self.x_build) ** 2, 1)) > (0.5 * self.skin) ** 2:
            self.build()
        i, j = self.pairs[:, 0], self.pairs[:, 1]
        d = self.x[j] - self.x[i]
        if self.box is not None:
            d -= self.box * np.round(d / self.box)
        r = np.sqrt(np.sum(d * d, 1))
        v, dv = self.tab(r)
        fij = (dv / r)[:, None] * d                           # force on i = +dv/r d (d = xj - xi)
        f = np.zeros_like(self.x)
        for k in range(3):
            f[:, k] = np.bincount(i, fij[:, k], self.N) - np.bincount(j, fij[:, k], self.N)
        self.pair_r, self.pair_v = r, v
        return f, float(np.sum(v))

    def epot_per(self):
        """potential energy of each molecule (half of its pair energies)"""
        i, j = self.pairs[:, 0], self.pairs[:, 1]
        return 0.5 * (np.bincount(i, self.pair_v, self.N) + np.bincount(j, self.pair_v, self.N))

    def baoab(self, dt, T, gamma, rng, nsteps=1):
        m = self.m if np.ndim(self.m) == 0 else self.m[:, None]
        c1 = np.exp(-gamma * dt); c2 = np.sqrt((1 - c1 ** 2) * T / m)
        for _ in range(nsteps):
            self.v += 0.5 * dt * self.f / m
            self.x += 0.5 * dt * self.v
            if gamma > 0:
                self.v = c1 * self.v + c2 * rng.normal(size=self.v.shape)
            self.x += 0.5 * dt * self.v
            self.f, self.epot = self.force()
            self.v += 0.5 * dt * self.f / m

    def ekin(self):
        m = self.m if np.ndim(self.m) == 0 else self.m[:, None]
        return 0.5 * float(np.sum(m * self.v ** 2))


def mackay(nshell, a):
    """Mackay icosahedron with nshell shells (13, 55, 147, ...), nearest-neighbour distance ~a along the radial rows"""
    p = (1 + np.sqrt(5)) / 2
    V = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            V += [(0, s1, s2 * p), (s1, s2 * p, 0), (s2 * p, 0, s1)]
    V = np.array(V, float); V /= np.linalg.norm(V[0])
    d = np.linalg.norm(V[:, None] - V[None], axis=-1)
    e = np.min(d[d > 1e-9])
    edges = [(i, j) for i in range(12) for j in range(i + 1, 12) if abs(d[i, j] - e) < 1e-6]
    faces = [(i, j, k) for i in range(12) for j in range(i + 1, 12) for k in range(j + 1, 12)
             if abs(d[i, j] - e) < 1e-6 and abs(d[i, k] - e) < 1e-6 and abs(d[j, k] - e) < 1e-6]
    pts = [np.zeros(3)]
    for k in range(1, nshell + 1):
        pts += [k * v for v in V]
        for i, j in edges:
            for mm in range(1, k):
                pts.append(k * V[i] + mm * (V[j] - V[i]))
        for i, j, l in faces:
            for aa in range(1, k):
                for bb in range(1, k - aa):
                    pts.append(k * V[i] + aa * (V[j] - V[i]) + bb * (V[l] - V[i]))
    return np.array(pts) * a


def minimise(x, table):
    from scipy.optimize import minimize
    sys_ = System(x, 1.0, table)
    def fg(z):
        sys_.x = z.reshape(-1, 3)
        f, e = sys_.force()
        return e, -f.ravel()
    r = minimize(fg, x.ravel(), jac=True, method="L-BFGS-B", options=dict(maxiter=20000, gtol=1e-8))
    return r.x.reshape(-1, 3), float(r.fun)


def gyration(x, m=None):
    c = x.mean(0) if m is None else np.average(x, 0, m)
    y = x - c
    G = (y.T @ y) / len(x)
    ev = np.sort(np.linalg.eigvalsh(G))
    rg2 = ev.sum()
    asph = ((ev[2] - ev[1]) ** 2 + (ev[2] - ev[0]) ** 2 + (ev[1] - ev[0]) ** 2) / (2 * rg2 ** 2)
    return float(np.sqrt(rg2)), float(asph)


def radial_density(x, edges):
    c = x.mean(0)
    r = np.linalg.norm(x - c, axis=1)
    h, _ = np.histogram(r, edges)
    vol = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    return h / vol


def gofr_core(x, rcore, rmax=10.0, nb=200):
    """g(r) seen from molecules within rcore of the centre (all partners), normalised by the mean core density"""
    c = x.mean(0)
    rc = np.linalg.norm(x - c, axis=1)
    core = np.where(rc < rcore)[0]
    tr = cKDTree(x)
    edges = np.linspace(0, rmax, nb + 1)
    h = np.zeros(nb)
    for i in core:
        idx = tr.query_ball_point(x[i], rmax)
        d = np.linalg.norm(x[idx] - x[i], axis=1)
        h += np.histogram(d[d > 1e-9], edges)[0]
    rho = len(core) / (4 / 3 * np.pi * rcore ** 3)
    shell = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    return 0.5 * (edges[1:] + edges[:-1]), h / (len(core) * rho * shell)

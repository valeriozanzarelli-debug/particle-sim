"""fluid_core.py - engines for the 'fluid in negative' study (F1-F3).

F1: 2D classical MD of soft repulsive balls (WCA, m = sigma = eps = 1) + one free disk ("hole") that repels the balls
    with a WCA wall shifted by R_w. Verlet neighbour list from scipy cKDTree (periodic box), velocity Verlet.
F2: radial Gross-Pitaevskii vortex, psi = f(r) e^{i theta}, units hbar = m = g n = 1 (xi = 1, c_s = 1, n = 1).
F3: Frenkel-Kontorova chain in sine-Gordon units: phi_n'' = (phi_{n+1} - 2 phi_n + phi_{n-1}) / h^2 - sin(phi_n)
    (+ F - alpha phi_n' when driven). Sound speed of the substrate-free chain c = 1; h = ball spacing / kink width.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np
from scipy.spatial import cKDTree

RC = 2.0 ** (1.0 / 6.0)


# ----------------------------------------------------------------------------------------------------- F1: MD
def wca_f_over_r(r2):
    """-dU/dr / r for WCA (eps = sigma = 1), valid for r < RC."""
    ir2 = 1.0 / r2
    ir6 = ir2 ** 3
    return 24.0 * ir2 * ir6 * (2.0 * ir6 - 1.0)


def wca_u(r2):
    ir6 = (1.0 / r2) ** 3
    return 4.0 * ir6 * (ir6 - 1.0) + 1.0


class FluidMD:
    def __init__(self, L, rho, Rw, M0, T0, seed, dt=0.005, skin=0.4):
        self.L, self.Rw, self.M0, self.dt, self.skin = float(L), float(Rw), float(M0), dt, skin
        rng = np.random.default_rng(seed)
        self.rng = rng
        # triangular lattice at density rho, then cut the disk region
        a = np.sqrt(2.0 / (np.sqrt(3.0) * rho))
        nx = int(round(L / a)); ny = int(round(L / (a * np.sqrt(3) / 2)))
        xs, ys = np.meshgrid(np.arange(nx), np.arange(ny), indexing="xy")
        px = (xs + 0.5 * (ys % 2)) * (L / nx); py = ys * (L / ny)
        pos = np.stack([px.ravel(), py.ravel()], 1)
        self.X = np.array([L / 2, L / 2])
        d = self.mic(pos - self.X)
        keep = np.hypot(d[:, 0], d[:, 1]) > Rw + 1.0
        self.pos = pos[keep].copy()
        self.N = len(self.pos)
        self.vel = rng.normal(0, np.sqrt(T0), (self.N, 2))
        self.vel -= self.vel.mean(0)
        self.V = np.zeros(2)
        self.build()
        self.f, self.F = self.forces()

    def mic(self, d):
        return d - self.L * np.round(d / self.L)

    def build(self):
        self.pos %= self.L
        t = cKDTree(self.pos, boxsize=self.L)
        p = t.query_pairs(RC + self.skin, output_type="ndarray")
        self.i, self.j = p[:, 0], p[:, 1]
        self.ref = self.pos.copy()
        self.Xref = self.X.copy()

    def forces(self):
        d = self.mic(self.pos[self.i] - self.pos[self.j])
        r2 = (d * d).sum(1)
        m = r2 < RC * RC
        fr = np.zeros_like(r2)
        fr[m] = wca_f_over_r(r2[m])
        fij = d * fr[:, None]
        f = np.zeros_like(self.pos)
        f[:, 0] = np.bincount(self.i, fij[:, 0], self.N) - np.bincount(self.j, fij[:, 0], self.N)
        f[:, 1] = np.bincount(self.i, fij[:, 1], self.N) - np.bincount(self.j, fij[:, 1], self.N)
        # disk wall: WCA in s = r - Rw
        dd = self.mic(self.pos - self.X)
        r = np.hypot(dd[:, 0], dd[:, 1])
        s = r - self.Rw
        w = s < RC
        F = np.zeros(2)
        if w.any():
            s2 = s[w] ** 2
            fs = wca_f_over_r(s2) * s[w]  # -dU/ds
            fw = dd[w] * (fs / r[w])[:, None]
            f[w] += fw
            F = -fw.sum(0)
        return f, F

    def step(self):
        dt = self.dt
        self.vel += 0.5 * dt * self.f
        self.V += 0.5 * dt * self.F / self.M0
        self.pos += dt * self.vel
        self.X += dt * self.V
        disp = np.abs(self.mic(self.pos - self.ref)).max() + np.abs(self.mic(self.X - self.Xref)).max()
        if disp > 0.5 * self.skin:
            self.build()
        self.f, self.F = self.forces()
        self.vel += 0.5 * dt * self.f
        self.V += 0.5 * dt * self.F / self.M0

    def rescale(self, T):
        k = (self.vel ** 2).sum() + self.M0 * (self.V ** 2).sum()
        Tn = k / (2 * (self.N + 1))
        s = np.sqrt(T / Tn)
        self.vel *= s; self.V *= s

    def state(self):
        return (self.pos.copy(), self.vel.copy(), self.X.copy(), self.V.copy())

    def set_state(self, st):
        self.pos, self.vel, self.X, self.V = [a.copy() for a in st]
        self.build()
        self.f, self.F = self.forces()

    def radial(self, edges):
        dd = self.mic(self.pos - self.X)
        r = np.hypot(dd[:, 0], dd[:, 1])
        h, _ = np.histogram(r, edges)
        return h

    def energy(self):
        d = self.mic(self.pos[self.i] - self.pos[self.j])
        r2 = (d * d).sum(1)
        m = r2 < RC * RC
        u = wca_u(r2[m]).sum()
        dd = self.mic(self.pos - self.X)
        s = np.hypot(dd[:, 0], dd[:, 1]) - self.Rw
        w = s < RC
        u += wca_u(s[w] ** 2).sum()
        return 0.5 * (self.vel ** 2).sum() + 0.5 * self.M0 * (self.V ** 2).sum() + u


# ----------------------------------------------------------------------------------------------------- F2: GP vortex
def gp_vortex(Rmax=200.0, dr=0.05, iters=None):
    """Radial GP vortex f(r) with winding 1: -1/2 (f'' + f'/r - f/r^2) + (f^2 - 1) f = 0, f(0) = 0, f(Rmax) = 1.
    Solved by Newton on a uniform grid (sparse tridiagonal)."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    r = np.arange(1, int(round(Rmax / dr)) + 1) * dr  # r_1..r_M, f(0) = 0 and f(r_M) = 1 enforced
    M = len(r)
    f = np.tanh(r / np.sqrt(2.0))
    f[-1] = 1.0
    for it in range(60):
        fp = np.r_[f[1:], 1.0]; fm = np.r_[0.0, f[:-1]]
        lap = (fp - 2 * f + fm) / dr ** 2 + (fp - fm) / (2 * dr * r) - f / r ** 2
        R = -0.5 * lap + (f ** 2 - 1) * f
        R[-1] = 0.0
        main = 0.5 * (2 / dr ** 2 + 1 / r ** 2) + 3 * f ** 2 - 1
        up = -0.5 * (1 / dr ** 2 + 1 / (2 * dr * r[:-1]))
        lo = -0.5 * (1 / dr ** 2 - 1 / (2 * dr * r[1:]))
        main[-1] = 1.0; lo[-1] = 0.0
        J = diags([lo, main, up], [-1, 0, 1], format="csc")
        df = spsolve(J, -R)
        f += df
        if np.abs(df).max() < 1e-12:
            break
    return r, f


def gp_energy_parts(r, f, R):
    """Energies inside radius R (per unit length, hbar = m = g n = 1): kinetic swirl, gradient, interaction,
    and missing atoms N_miss = int (1 - f^2) dA."""
    dr = r[1] - r[0]
    m = r <= R
    fp = np.gradient(f, dr)
    w = 2 * np.pi * r[m] * dr
    e_swirl = (0.5 * f[m] ** 2 / r[m] ** 2 * w).sum()
    e_grad = (0.5 * fp[m] ** 2 * w).sum()
    e_int = (0.5 * (f[m] ** 2 - 1) ** 2 * w).sum()
    nmiss = ((1 - f[m] ** 2) * w).sum()
    # add the tiny cell 0 < r < dr (f ~ a r): swirl contributes 0.5 a^2 * pi dr^2 , negligible but included
    a = f[0] / r[0]
    e_swirl += 0.5 * a * a * np.pi * dr ** 2
    nmiss += np.pi * dr ** 2
    e_int += 0.5 * np.pi * dr ** 2
    return dict(swirl=e_swirl, grad=e_grad, inter=e_int, total=e_swirl + e_grad + e_int, nmiss=nmiss)


# ----------------------------------------------------------------------------------------------------- F3: FK chain
def fk_acc(phi, h, F=0.0, alpha=0.0, v=None, bc="free"):
    lap = np.empty_like(phi)
    lap[1:-1] = phi[2:] - 2 * phi[1:-1] + phi[:-2]
    if bc == "free":
        lap[0] = phi[1] - phi[0]
        lap[-1] = phi[-2] - phi[-1]
    a = lap / h ** 2 - np.sin(phi) + F
    if alpha and v is not None:
        a = a - alpha * v
    return a


def fk_energy_density(phi, v, h):
    e = 0.5 * v ** 2 + (1 - np.cos(phi))
    g = 0.5 * ((phi[1:] - phi[:-1]) / h) ** 2
    e = e.copy()
    e[:-1] += 0.5 * g; e[1:] += 0.5 * g
    return e * h  # energy per site


def kink_profile(x, x0, gam, sgn=1.0):
    return 4 * np.arctan(np.exp(sgn * gam * (x - x0)))


def fit_kink(x, phi, sgn=1.0):
    """Fit phi ~ 4 arctan exp(sgn (x - X)/w) (mod nothing; phi in [0, 2pi]). Returns X, w."""
    from scipy.optimize import curve_fit
    k = np.argmin(np.abs(phi - np.pi))
    X0 = x[k]
    sel = np.abs(x - X0) < 12
    def mdl(xx, X, w):
        return 4 * np.arctan(np.exp(sgn * (xx - X) / w))
    p, _ = curve_fit(mdl, x[sel], phi[sel], p0=[X0, 1.0])
    return p[0], abs(p[1])

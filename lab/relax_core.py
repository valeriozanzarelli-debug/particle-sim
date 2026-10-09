"""relax_core.py - de Broglie-Bohm surfers in a 2-D box [0,pi]^2 (hbar = m = 1).

psi(x,y,t) = sum_{n,m} c_nm (2/pi) sin(n x) sin(m y) exp(-i (n^2+m^2) t / 2),  v = Im(grad psi / psi).
Tools: field + gradient, sampling, RK4 with per-surfer sub-cycling near nodes, coarse-grained H,
vortex counting. Used by relax_box.py and relax_slit.py.
"""
import numpy as np

PI = np.pi


class Box:
    def __init__(self, M, seed=1):
        s = int(round(np.sqrt(M)))
        assert s * s == M
        self.nmax = s
        rng = np.random.default_rng(seed)
        ph = rng.uniform(0, 2 * PI, (s, s))
        if M == 1:
            ph[:] = 0.0
        self.c = np.exp(1j * ph) / np.sqrt(M) * (2 / PI)
        n = np.arange(1, s + 1)
        self.n = n
        self.E = 0.5 * (n[:, None] ** 2 + n[None, :] ** 2)

    def A(self, t):
        return self.c * np.exp(-1j * self.E * t)

    def field(self, x, y, t, grad=True):
        # E_nm = (n^2 + m^2)/2 is separable: fold the time phases into the 1-D factors (t may be per-surfer)
        n = self.n[:, None]
        t = np.asarray(t, dtype=float)
        ex, ey = np.exp(-0.5j * n ** 2 * t), np.exp(-0.5j * n ** 2 * t)
        sx, sy = np.sin(n * x) * ex, np.sin(n * y) * ey
        B = self.c @ sy                             # (nmax, N)
        psi = (sx * B).sum(0)
        if not grad:
            return psi
        cx, cy = n * np.cos(n * x) * ex, n * np.cos(n * y) * ey
        dx = (cx * B).sum(0)
        dy = (sx * (self.c @ cy)).sum(0)
        return psi, dx, dy

    def vel(self, x, y, t):
        psi, dx, dy = self.field(x, y, t)
        r2 = (psi.real ** 2 + psi.imag ** 2) + 1e-300
        vx = (psi.real * dx.imag - psi.imag * dx.real) / r2
        vy = (psi.real * dy.imag - psi.imag * dy.real) / r2
        return vx, vy

    def dens_grid(self, t, ng):
        g = (np.arange(ng) + 0.5) * PI / ng
        X, Y = np.meshgrid(g, g, indexing="ij")
        psi = self.field(X.ravel(), Y.ravel(), t, grad=False).reshape(ng, ng)
        return psi


def sample(rng, N, wfun, wmax):
    """rejection sampling in [0,pi]^2 for unnormalised density wfun(x,y) <= wmax"""
    xs, ys = [], []
    got = 0
    while got < N:
        x = rng.uniform(0, PI, 4 * N); y = rng.uniform(0, PI, 4 * N)
        ok = rng.uniform(0, wmax, 4 * N) < wfun(x, y)
        xs.append(x[ok]); ys.append(y[ok]); got += ok.sum()
    return np.concatenate(xs)[:N], np.concatenate(ys)[:N]


def rk4_step(box, x, y, t, dt, ds=0.01, nmax_sub=256):
    """advance all surfers by dt; each surfer takes n_i = ceil(|v| dt / ds) RK4 sub-steps"""
    vx, vy = box.vel(x, y, t)
    sp = np.hypot(vx, vy)
    nsub = np.clip(np.ceil(sp * dt / ds), 1, nmax_sub).astype(int)
    h = dt / nsub
    x = x.copy(); y = y.copy()
    for k in range(nsub.max()):
        a = np.nonzero(nsub > k)[0]
        xa, ya, ha = x[a], y[a], h[a]
        ta = t + k * ha
        k1x, k1y = box.vel(xa, ya, ta)
        k2x, k2y = box.vel(xa + 0.5 * ha * k1x, ya + 0.5 * ha * k1y, ta + 0.5 * ha)
        k3x, k3y = box.vel(xa + 0.5 * ha * k2x, ya + 0.5 * ha * k2y, ta + 0.5 * ha)
        k4x, k4y = box.vel(xa + ha * k3x, ya + ha * k3y, ta + ha)
        x[a] = xa + ha / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        y[a] = ya + ha / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)
    eps = 1e-9
    np.clip(x, eps, PI - eps, out=x); np.clip(y, eps, PI - eps, out=y)
    return x, y


def coarse(counts, ref, K):
    """aggregate a (32,32) counts / reference map to K x K cells"""
    f = counts.shape[0] // K
    c = counts.reshape(K, f, K, f).sum((1, 3))
    r = ref.reshape(K, f, K, f).sum((1, 3))
    return c, r


def Hbar(counts, refmass):
    """coarse-grained H = sum rho ln(rho/|psi|^2) with cell masses (both normalised)"""
    p = counts / counts.sum()
    q = refmass / refmass.sum()
    m = p > 0
    return float((p[m] * np.log(p[m] / q[m])).sum())


def vortices(psi):
    """number of phase-winding plaquettes on a grid of complex values"""
    ph = np.angle(psi)
    def w(a, b):
        d = b - a
        return (d + PI) % (2 * PI) - PI
    s = (w(ph[:-1, :-1], ph[1:, :-1]) + w(ph[1:, :-1], ph[1:, 1:])
         + w(ph[1:, 1:], ph[:-1, 1:]) + w(ph[:-1, 1:], ph[:-1, :-1]))
    return int((np.abs(s) > PI).sum())

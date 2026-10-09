"""dirac_core.py - the "Dirac surfer" rule (frozen before any experiment was run).

Units: hbar = c = m_e = 1  (length = reduced Compton wavelength lambda_C = hbar/(m c) = 3.8616e-13 m,
time = hbar/(m c^2) = 1.2881e-21 s, energy = m c^2 = 510998.95 eV).

THE RULE (one tick of length dt):
  wave balls : every wave ball carries a SPINOR (2 complex numbers in 1+1 D and 2+1 D, 4 in 3+1 D);
               they move with the free Dirac propagator  U(k,dt) = cos(E dt) - i sin(E dt) H(k)/E,
               H(k) = alpha.k + beta m,  E = sqrt(k^2 + m^2)  (spectral split-operator: exact for free
               propagation, no fermion doubling because every k on the FFT grid has its own exact E(k)).
               Representation: 1+1 D  alpha = sigma_x, beta = sigma_z;
                               2+1 D  alpha = (sigma_x, sigma_y), beta = sigma_z.
  surfer     : reads the spinor psi at its own position (cubic interpolation of the spinor components
               from the wave balls around it) and moves with
                     v = psi^dagger alpha psi / psi^dagger psi           (Bohm 1953, Holland 1993 ch.12)
               |v| <= c holds for ANY spinor (Cauchy-Schwarz), so it holds also after interpolation.
               Integrator: classical RK4 in time, the wave read at the exact stage times.
  birth      : surfers placed with probability psi^dagger psi (the only random step).
No other parameter. Frozen values below.
"""
import numpy as np
from scipy.ndimage import map_coordinates, spline_filter

# ---------------- frozen parameters (2026-10-09, before any test) ----------------
DT = 0.02            # surfer time step, in hbar/(m c^2)  (ZBW period is pi -> 157 steps per period)
INTERP_ORDER = 3     # cubic interpolation of the spinor components
SEED = 20261009

# ---------------- SI constants (CODATA 2018) ----------------
HBAR = 1.054571817e-34
C = 299792458.0
ME = 9.1093837015e-31
QE = 1.602176634e-19
ALPHA = 7.2973525693e-3
LAMBDA_C = HBAR / (ME * C)         # 3.8616e-13 m
T_UNIT = HBAR / (ME * C**2)        # 1.2881e-21 s
MC2_EV = ME * C**2 / QE            # 510998.95 eV

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


class Dirac:
    """Free Dirac field on a periodic grid in d = 1 or 2 space dimensions (2-component spinor)."""

    def __init__(self, n, L, dim, m=1.0):
        self.n, self.L, self.dim, self.m = n, L, dim, m
        self.dx = L / n
        x = (np.arange(n) - n // 2) * self.dx
        k = 2 * np.pi * np.fft.fftfreq(n, d=self.dx)
        if dim == 1:
            self.x = (x,)
            self.k = (k,)
        else:
            X, Y = np.meshgrid(x, x, indexing="ij")
            KX, KY = np.meshgrid(k, k, indexing="ij")
            self.x = (X, Y)
            self.k = (KX, KY)
        # H(k) = hx sx + hy sy + hz sz
        self.hx = self.k[0]
        self.hy = self.k[1] if dim == 2 else np.zeros_like(self.k[0])
        self.hz = m * np.ones_like(self.k[0])
        self.E = np.sqrt(self.hx**2 + self.hy**2 + self.hz**2)

    # ---- spectral helpers ----
    def fft(self, psi):
        return np.fft.fftn(psi, axes=tuple(range(1, self.dim + 1)))

    def ifft(self, phi):
        return np.fft.ifftn(phi, axes=tuple(range(1, self.dim + 1)))

    def apply_H(self, phi, s=1.0):
        """phi (2,...) in k space -> (s*H(k)) phi"""
        a, b = phi
        hx, hy, hz = self.hx, self.hy, self.hz
        return s * np.array([hz * a + (hx - 1j * hy) * b, (hx + 1j * hy) * a - hz * b])

    def project(self, psi, sign=+1):
        """positive (sign=+1) or negative energy part: P = (1 + sign H/E)/2"""
        phi = self.fft(psi)
        phi = 0.5 * (phi + sign * self.apply_H(phi) / self.E)
        return self.ifft(phi)

    def evolve(self, psi0_k, t):
        """exact free evolution from the k-space spinor psi0_k to time t (returns x-space spinor)"""
        c, s = np.cos(self.E * t), np.sin(self.E * t)
        phi = c * psi0_k - 1j * s * self.apply_H(psi0_k) / self.E
        return self.ifft(phi)

    # ---- observables ----
    def density_current(self, psi):
        a, b = psi
        rho = (np.abs(a)**2 + np.abs(b)**2).real
        jx = 2 * np.real(np.conj(a) * b)
        jy = 2 * np.imag(np.conj(a) * b)
        return rho, jx, jy

    def mean(self, psi, f):
        rho = self.density_current(psi)[0]
        return np.sum(rho * f) / np.sum(rho)

    def norm(self, psi):
        return np.sum(self.density_current(psi)[0]) * self.dx**self.dim

    # ---- the surfer rule ----
    def prepare(self, psi):
        """spline coefficients of the 4 real fields (Re, Im of the 2 components). Pure speed-up: identical
        to letting map_coordinates prefilter at every call (checked in the self-test)."""
        return [spline_filter(f, order=INTERP_ORDER, mode="grid-wrap") for comp in psi for f in (comp.real, comp.imag)]

    def _read_spinor(self, coef, pos):
        """cubic interpolation of the spinor components at surfer positions pos (dim, N)"""
        idx = [(p / self.dx + self.n // 2) for p in pos]
        v = [map_coordinates(c, idx, order=INTERP_ORDER, mode="grid-wrap", prefilter=False) for c in coef]
        return v[0] + 1j * v[1], v[2] + 1j * v[3]

    def surfer_velocity(self, psi, pos, prepared=False):
        a, b = self._read_spinor(psi if prepared else self.prepare(psi), pos)
        rho = np.abs(a)**2 + np.abs(b)**2
        vx = 2 * np.real(np.conj(a) * b) / rho
        if self.dim == 1:
            return np.array([vx])
        vy = 2 * np.imag(np.conj(a) * b) / rho
        return np.array([vx, vy])

    def birth(self, psi, N, rng):
        """N surfers distributed with psi^dagger psi (sampled on the grid + uniform jitter in the cell)"""
        rho = self.density_current(psi)[0].ravel()
        p = rho / rho.sum()
        cells = rng.choice(rho.size, size=N, p=p)
        if self.dim == 1:
            pos = np.array([self.x[0][cells]]) + (rng.random((1, N)) - 0.5) * self.dx
        else:
            i, j = np.unravel_index(cells, (self.n, self.n))
            pos = np.array([self.x[0][i, j], self.x[1][i, j]]) + (rng.random((2, N)) - 0.5) * self.dx
        return pos

    def run_surfers(self, psi0, pos0, T, dt=DT, record_every=1, on_step=None):
        """RK4 surfers. Returns times, trajectory array (nrec, dim, N), max speed seen."""
        psi0_k = self.fft(psi0)
        nsteps = int(round(T / dt))
        pos = pos0.copy()
        traj, times = [pos.copy()], [0.0]
        vmax = 0.0
        psi_t = psi0
        c_t = self.prepare(psi0)
        # exact half-step propagator U(dt/2) as a 2x2 field (group property: U(t+dt/2) = U(dt/2) U(t),
        # verified to 1e-16 in the self-test) - a speed-up only, same numbers as evolve(psi0_k, t)
        ch, sh = np.cos(self.E * dt / 2), np.sin(self.E * dt / 2)
        f = -1j * sh / self.E
        U = (ch + f * self.hz, f * (self.hx - 1j * self.hy), f * (self.hx + 1j * self.hy), ch - f * self.hz)
        half = lambda p: np.array([U[0] * p[0] + U[1] * p[1], U[2] * p[0] + U[3] * p[1]])
        phik = psi0_k
        for s in range(nsteps):
            t = s * dt
            phik = half(phik)
            c_h = self.prepare(self.ifft(phik))
            phik = half(phik)
            psi_1 = self.ifft(phik)
            c_1 = self.prepare(psi_1)
            k1 = self.surfer_velocity(c_t, pos, True)
            k2 = self.surfer_velocity(c_h, pos + 0.5 * dt * k1, True)
            k3 = self.surfer_velocity(c_h, pos + 0.5 * dt * k2, True)
            k4 = self.surfer_velocity(c_1, pos + dt * k3, True)
            for kk in (k1, k2, k3, k4):
                vmax = max(vmax, float(np.max(np.sqrt(np.sum(kk**2, axis=0)))))
            pos = pos + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            psi_t, c_t = psi_1, c_1
            if on_step is not None:
                on_step(s + 1, (s + 1) * dt, psi_t, pos)
            if (s + 1) % record_every == 0:
                traj.append(pos.copy())
                times.append((s + 1) * dt)
        return np.array(times), np.array(traj), vmax


def gaussian_packet(D, sigma, k0, spinor):
    """psi = exp(-r^2/(4 sigma^2) + i k0.r) * spinor, normalized (so |psi|^2 has std sigma per axis)"""
    if D.dim == 1:
        r2 = D.x[0]**2
        ph = k0[0] * D.x[0]
    else:
        r2 = D.x[0]**2 + D.x[1]**2
        ph = k0[0] * D.x[0] + k0[1] * D.x[1]
    g = np.exp(-r2 / (4 * sigma**2) + 1j * ph)
    s = np.asarray(spinor, complex)
    s = s / np.linalg.norm(s)
    psi = np.array([g * s[0], g * s[1]])
    return psi / np.sqrt(D.norm(psi))


if __name__ == "__main__":
    # self-test of the machinery (not of the physics)
    import time
    t0 = time.time()
    D = Dirac(512, 200.0, 1)
    psi = gaussian_packet(D, 5.0, [0.5], [1, 1])
    pk = D.fft(psi)
    # unitarity and group property of the exact propagator
    p1 = D.evolve(pk, 3.7)
    p2 = D.evolve(D.fft(D.evolve(pk, 1.2)), 2.5)
    print("norm drift", abs(D.norm(p1) - 1), " U(1.2)U(2.5)-U(3.7)", np.max(np.abs(p1 - p2)), flush=True)
    # projector idempotent and positive+negative = 1
    pp = D.project(psi, +1); pm = D.project(psi, -1)
    print("P+ + P- - 1:", np.max(np.abs(pp + pm - psi)), " P+P+ - P+:", np.max(np.abs(D.project(pp, 1) - pp)), flush=True)
    # |v|<=1 for random spinors
    rng = np.random.default_rng(SEED)
    a = rng.normal(size=100000) + 1j * rng.normal(size=100000)
    b = rng.normal(size=100000) + 1j * rng.normal(size=100000)
    v = np.hypot(2 * np.real(np.conj(a) * b), 2 * np.imag(np.conj(a) * b)) / (abs(a)**2 + abs(b)**2)
    print("max |v| over 1e5 random spinors:", v.max(), flush=True)
    # incremental propagation in run_surfers == direct evolve
    pos0 = D.birth(psi, 50, rng)
    got = []
    D.run_surfers(psi, pos0, 1.0, on_step=lambda s, t, ps, p: got.append(np.max(np.abs(ps - D.evolve(pk, t)))))
    print("incremental half-step propagator vs direct, max diff over 50 steps:", max(got), flush=True)
    # prefiltered interpolation == map_coordinates with its own prefilter (speed-up only)
    D2 = Dirac(128, 40.0, 2)
    p2 = gaussian_packet(D2, 3.0, [0.4, 0.1], [1, 1j])
    P = rng.uniform(-8, 8, (2, 500))
    idx = [(q / D2.dx + D2.n // 2) for q in P]
    ref = [map_coordinates(f, idx, order=INTERP_ORDER, mode="grid-wrap") for c in p2 for f in (c.real, c.imag)]
    new = D2.prepare(p2)
    new = [map_coordinates(c, idx, order=INTERP_ORDER, mode="grid-wrap", prefilter=False) for c in new]
    print("prefilter speed-up, max difference:", max(np.max(np.abs(a - b)) for a, b in zip(ref, new)), flush=True)
    print("SI: lambda_C = %.4e m, hbar/2mc = %.4e m, 2mc^2/hbar = %.4e rad/s, t_unit = %.4e s"
          % (LAMBDA_C, LAMBDA_C / 2, 2 / T_UNIT, T_UNIT))
    print("time %.2f s" % (time.time() - t0))

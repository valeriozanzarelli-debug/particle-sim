"""emiss_core.py - study "Emissione vera" (docs/EMISSIONE.md): the zig-zag rule for the electron (3-D Dirac walk)
and for the light (3-D spin-1 walk, F = E + iB), coupled with the TRUE charge e = sqrt(4 pi alpha) (no K boost).

Same rule as lab/zz_core.py (checked against it in the self-test: python lab/emiss_core.py) and as
src/atomo/zz.rs + accoppiato.rs (same coupling: the light's E at the atom, Gaussian average over the light cells,
enters the electron's phase as e E.r; the electron's conserved current -e d<r>/dt, deposited with the same
Gaussian, is a source of the light, dF/dt = -i curl F - J).  Only numba kernels are new (one thread).

Units hbar = c = m_e = 1 (length lambda_C, time hbar/mc^2, energy mc^2).  Heaviside-Lorentz: e^2 = 4 pi alpha,
field energy 1/2 sum |F|^2 b^3.

Two ways to couple (both use exactly the same walks):
  * SEMICLASSICAL ("sc"): one electron wave psi, one REAL light field (one walk).  As in accoppiato.rs.
  * ONE-EXCITATION QUANTUM ("qe"): the state is |psi(x)> (x) |no photon>  +  |1s> (x) |one photon, amplitude Phi>.
    The photon amplitude is a COMPLEX field: two light walks (real and imaginary part), the same rule.
    Coupling (rotating-wave, dipole, one excitation):  i d psi/dt = H psi + e (Phi_E(atom) . r) phi_1s(t)
    and the photon amplitude is sourced by the TRANSITION current  J = -e d<phi_1s(t)| r |psi>/dt.
    phi_1s(t) is the lattice 1s (listened) evolved by the same walk.  A pure 2p (psi = 2p, Phi = 0) decays by itself.
"""
import os
import numpy as np
import numba as nb
from zz_core import ALPHA, order_of_tick, Walk, RadialWalk, dirac_formula, T_UNIT, MC2_EV, SEED

E_CHARGE = np.sqrt(4 * np.pi * ALPHA)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
CK = os.path.join(OUT, "emiss_ck")
os.makedirs(CK, exist_ok=True)
S2 = 1 / np.sqrt(2)

# eigenvectors of sigma_j (u+ : sigma_j = +1, u- : sigma_j = -1); the shift does not depend on their phases
UP = np.array([[S2, S2], [S2, 1j * S2], [1, 0]], np.complex128)
DN = np.array([[S2, -S2], [S2, -1j * S2], [0, 1]], np.complex128)


# ------------------------------------------------------------------ Dirac walk (electron) ---------------------------
@nb.njit(cache=True)
def _dirac_axis(psi, out, j, up, dn):
    """translate along axis j: L block (comp 0,1): sigma=+1 part moves -1, -1 part +1; R block (2,3) the opposite."""
    n = psi.shape[1]
    u0, u1, d0, d1 = up[j, 0], up[j, 1], dn[j, 0], dn[j, 1]
    cu0, cu1, cd0, cd1 = np.conj(u0), np.conj(u1), np.conj(d0), np.conj(d1)
    for i in range(n):
        for k in range(n):
            for l in range(n):
                # neighbours: s_m = site - e_j, s_p = site + e_j
                if j == 0:
                    im, km, lm, ip, kp, lp = (i - 1) % n, k, l, (i + 1) % n, k, l
                elif j == 1:
                    im, km, lm, ip, kp, lp = i, (k - 1) % n, l, i, (k + 1) % n, l
                else:
                    im, km, lm, ip, kp, lp = i, k, (l - 1) % n, i, k, (l + 1) % n
                # L: u+ part moves -1 -> comes from s_p ; u- part moves +1 -> comes from s_m
                a = psi[0, ip, kp, lp]; b = psi[1, ip, kp, lp]
                cp = cu0 * a + cu1 * b
                a = psi[0, im, km, lm]; b = psi[1, im, km, lm]
                cm = cd0 * a + cd1 * b
                out[0, i, k, l] = u0 * cp + d0 * cm
                out[1, i, k, l] = u1 * cp + d1 * cm
                # R: u+ part moves +1 -> from s_m ; u- part moves -1 -> from s_p
                a = psi[2, im, km, lm]; b = psi[3, im, km, lm]
                cp = cu0 * a + cu1 * b
                a = psi[2, ip, kp, lp]; b = psi[3, ip, kp, lp]
                cm = cd0 * a + cd1 * b
                out[2, i, k, l] = u0 * cp + d0 * cm
                out[3, i, k, l] = u1 * cp + d1 * cm


@nb.njit(cache=True)
def _dirac_coin(psi, ct, st, pc, px, py, pz):
    """mass flip + phase: phase = pc[site] * px[i] * py[k] * pz[l] (Coulomb and e E.r)"""
    n = psi.shape[1]
    for i in range(n):
        for k in range(n):
            for l in range(n):
                ph = pc[i, k, l] * px[i] * py[k] * pz[l]
                for s in range(2):
                    L = psi[s, i, k, l]
                    R = psi[2 + s, i, k, l]
                    psi[s, i, k, l] = ph * (ct * L - 1j * st * R)
                    psi[2 + s, i, k, l] = ph * (ct * R - 1j * st * L)


class Dirac3D:
    """3-D zig-zag electron in a periodic box, nucleus (point, charge Z) between the 8 central sites.
    tick(psi, eE) advances one tick with the field e*E (vector, uniform over the atom: dipole form)."""

    def __init__(self, n, a, Z=60.0, mass=1.0):
        self.n, self.a, self.Z = n, a, Z
        self.x = (np.arange(n) - n / 2 + 0.5) * a
        X, Y, Zc = np.meshgrid(self.x, self.x, self.x, indexing="ij")
        r = np.sqrt(X**2 + Y**2 + Zc**2)
        self.r = r
        self.pc = np.exp(1j * Z * ALPHA / r * a)            # e^{-i V a}, V = -Z alpha / r
        self.theta = mass * a
        self.buf = np.empty((4, n, n, n), np.complex128)
        self.one = np.ones(n, np.complex128)

    def tick(self, psi, tick_count, eE=None):
        b = self.buf
        cur = psi
        for j in order_of_tick(tick_count):
            _dirac_axis(cur, b, j, UP, DN)
            cur, b = b, cur
        # after 3 axes cur is the buffer (odd number of swaps): copy back
        psi[...] = cur
        if eE is None:
            px = py = pz = self.one
        else:
            px = np.exp(-1j * eE[0] * self.x * self.a)
            py = np.exp(-1j * eE[1] * self.x * self.a)
            pz = np.exp(-1j * eE[2] * self.x * self.a)
        _dirac_coin(psi, np.cos(self.theta), np.sin(self.theta), self.pc, px, py, pz)
        return psi

    def position(self, psi):
        rho = np.sum(np.abs(psi)**2, axis=0)
        x = self.x
        return np.array([np.einsum("ijk,i->", rho, x), np.einsum("ijk,j->", rho, x), np.einsum("ijk,k->", rho, x)])

    def dip_matrix(self, phi, psi):
        """<phi| r |psi> (complex 3-vector)"""
        w = np.sum(np.conj(phi) * psi, axis=0)
        x = self.x
        return np.array([np.einsum("ijk,i->", w, x), np.einsum("ijk,j->", w, x), np.einsum("ijk,k->", w, x)])


def _add_rphi(self, psi, coef, phi):
    """psi += (coef . r) phi   (coef complex 3-vector; r = lattice position)"""
    x = self.x
    f = coef[0] * x[:, None, None] + coef[1] * x[None, :, None] + coef[2] * x[None, None, :]
    psi += f[None] * phi


Dirac3D.add_rphi = _add_rphi


def _tick_phase(self, psi, phi_vec):
    """psi *= exp(-i phi_vec . r)  (a uniform-field kick: e E dt . r)"""
    x = self.x
    f = np.exp(-1j * (phi_vec[0] * x[:, None, None] + phi_vec[1] * x[None, :, None] + phi_vec[2] * x[None, None, :]))
    psi *= f[None]


Dirac3D.tick_phase = _tick_phase


def norm2(psi):
    return float(np.sum(np.abs(psi)**2))


def ovl(a, b):
    return complex(np.vdot(a, b))


# ------------------------------------------------------------------ light walk (spin 1) -----------------------------
@nb.njit(cache=True)
def _light_axis(F, out, j):
    """translate along axis j: F_j stays; c+ = (F_k - i F_l)/sqrt2 moves +1, c- = (F_k + i F_l)/sqrt2 moves -1
    ((j,k,l) cyclic; (S_j)_{kl} = -i eps_jkl: eigenvector +1 is (1, i)/sqrt2 on (F_k, F_l))."""
    n = F.shape[1]
    kk = (j + 1) % 3
    ll = (j + 2) % 3
    s2 = 0.7071067811865476
    for i in range(n):
        for k in range(n):
            for l in range(n):
                if j == 0:
                    im, km, lm, ip, kp, lp = (i - 1) % n, k, l, (i + 1) % n, k, l
                elif j == 1:
                    im, km, lm, ip, kp, lp = i, (k - 1) % n, l, i, (k + 1) % n, l
                else:
                    im, km, lm, ip, kp, lp = i, k, (l - 1) % n, i, k, (l + 1) % n
                cp = (F[kk, im, km, lm] - 1j * F[ll, im, km, lm]) * s2      # c+ comes from site - e_j
                cm = (F[kk, ip, kp, lp] + 1j * F[ll, ip, kp, lp]) * s2      # c- comes from site + e_j
                out[j, i, k, l] = F[j, i, k, l]
                out[kk, i, k, l] = (cp + cm) * s2
                out[ll, i, k, l] = 1j * (cp - cm) * s2


@nb.njit(cache=True)
def _mask_energy(F, m):
    """F *= m (sponge); returns (energy removed, energy left) as sums of |F|^2"""
    n = F.shape[1]
    tolta = 0.0
    resta = 0.0
    for i in range(n):
        for k in range(n):
            for l in range(n):
                w = m[i, k, l]
                e = 0.0
                for c in range(3):
                    v = F[c, i, k, l]
                    e += v.real * v.real + v.imag * v.imag
                if w < 1.0:
                    for c in range(3):
                        F[c, i, k, l] *= w
                tolta += (1.0 - w * w) * e
                resta += w * w * e
    return tolta, resta


class Light3D:
    """spin-1 zig-zag light in a periodic box of n^3 cells of size b, with a sponge (cubic profile, as in
    accoppiato.rs: one crossing attenuates by ~e^{-7.5}), a Gaussian average (sigma cells) around the centre
    for the field at the atom and for the deposit of the atom's current."""

    def __init__(self, n, b, sponge=12, sigma=1.0, strength=2.5):
        self.n, self.b, self.sigma = n, b, sigma
        self.F = np.zeros((3, n, n, n), np.complex128)
        self.buf = np.empty_like(self.F)
        q = np.arange(n)
        d = np.maximum(sponge - np.minimum(q, n - 1 - q), 0) / max(sponge, 1)
        d3 = d**3
        self.mask = np.exp(-strength * (d3[:, None, None] + d3[None, :, None] + d3[None, None, :]))
        h = n // 2
        r = int(np.ceil(3 * sigma))
        idx, w = [], []
        for i in range(-r, r + 1):
            for j in range(-r, r + 1):
                for k in range(-r, r + 1):
                    ww = np.exp(-(i * i + j * j + k * k) / (2 * sigma**2))
                    if ww < 1e-4:
                        continue
                    idx.append((h + i, h + j, h + k))
                    w.append(ww)
        self.idx = tuple(np.array(idx).T)
        self.offs = np.array(idx) - h
        self.w = np.array(w) / np.sum(w)
        self.absorbed = 0.0
        self.tick_count = 0

    def tick(self):
        cur, b = self.F, self.buf
        for j in order_of_tick(self.tick_count):
            _light_axis(cur, b, j)
            cur, b = b, cur
        self.F, self.buf = cur, b
        self.tick_count += 1

    def field_at_atom(self):
        """Gaussian-averaged E = Re F at the atom"""
        return np.array([np.sum(self.w * self.F[c][self.idx].real) for c in range(3)])

    def deposit(self, I):
        """add the source of one light tick: dF = -J dt, J = I w / b^3, dt = b  (I = charge x velocity, real)"""
        for c in range(3):
            if I[c] != 0.0:
                self.F[c][self.idx] -= I[c] * self.w / self.b**2

    def make_transverse(self):
        """Deposit only the TRANSVERSE part of the source (Coulomb gauge: the electron does not feel its own
        static Coulomb field through the light).  On the lattice the two-tick light walk W(k) is a real rotation;
        its axis n(k) is the lattice's longitudinal (static) mode, exactly.  Patterns T_c = P_perp (w e_c),
        P_perp(k) = 1 - n n^T, computed by FFT once and cached.  M = <w| T> is the instantaneous self-field matrix."""
        from zz_core import bloch_batch
        n = self.n
        f = os.path.join(CK, f"transverse_{n}_{self.sigma:g}.npy")
        if os.path.exists(f):
            self.T = np.load(f)
        else:
            kk = 2 * np.pi * np.fft.fftfreq(n)
            K = np.stack([q.ravel() for q in np.meshgrid(kk, kk, kk, indexing="ij")], -1)
            P = np.empty((K.shape[0], 3, 3))
            for b0 in range(0, K.shape[0], 100000):
                sl = slice(b0, b0 + 100000)
                R = bloch_batch(K[sl], spin=1).real
                ax = np.stack([R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], -1)
                nr = np.linalg.norm(ax, axis=1)
                bad = nr < 1e-9
                if np.any(bad):              # rotation by 0 or pi: axis from the eigenvector of eigenvalue 1
                    for i in np.where(bad)[0]:
                        w_, v_ = np.linalg.eig(R[i])
                        j = np.argmin(np.abs(w_ - 1))
                        ax[i] = np.real(v_[:, j]); nr[i] = np.linalg.norm(ax[i])
                ax /= nr[:, None]
                P[sl] = np.eye(3)[None] - ax[:, :, None] * ax[:, None, :]
            P[0] = 2.0 / 3.0 * np.eye(3)          # k = 0
            w = np.zeros((n, n, n))
            w[self.idx] = self.w
            wk = np.fft.fftn(w).ravel()
            T = np.empty((3, 3, n, n, n), np.float32)
            for c in range(3):
                Tc = np.fft.ifftn((P[:, :, c] * wk[:, None]).T.reshape(3, n, n, n), axes=(1, 2, 3))
                print(f"  transverse pattern {c}: max |imag|/max |real| = {np.max(np.abs(Tc.imag)) / np.max(np.abs(Tc.real)):.1e}")
                T[c] = Tc.real
            np.save(f, T)
            self.T = T
        self.M = np.array([[np.sum(self.w * self.T[c][cp][self.idx]) for c in range(3)] for cp in range(3)])

    def deposit_T(self, I):
        """transverse deposit of one light tick: F -= sum_c I_c T_c / b^2 (I real 3-vector)"""
        for c in range(3):
            if I[c] != 0.0:
                self.F.real -= (I[c] / self.b**2) * self.T[c]

    def self_field(self, I):
        """field at the atom created by the deposit of I: -M I / b^2"""
        return -self.M @ np.asarray(I) / self.b**2

    def sponge(self):
        t, r = _mask_energy(self.F, self.mask)
        v = 0.5 * self.b**3
        self.absorbed += v * t
        return v * r

    def energy(self):
        return 0.5 * self.b**3 * float(np.sum(np.abs(self.F)**2))

    def form_factor2(self, k, nang=64):
        """<|W(k n)|^2 sin^2>/<sin^2> over directions n (dipole along z): the power reduction of the Gaussian
        average at wave number k (lambda_C^-1), with the discrete weights (computed, not fitted)."""
        u = (np.arange(nang) + 0.5) / nang * 2 - 1          # cos theta (Gauss-Legendre would be finer)
        ph = (np.arange(2 * nang) + 0.5) / (2 * nang) * 2 * np.pi
        U, P = np.meshgrid(u, ph, indexing="ij")
        S = np.sqrt(1 - U**2)
        nx, ny, nz = S * np.cos(P), S * np.sin(P), U
        kb = k * self.b
        o = self.offs
        Wk = np.zeros_like(U, complex)
        for q in range(len(self.w)):
            Wk += self.w[q] * np.exp(1j * kb * (o[q, 0] * nx + o[q, 1] * ny + o[q, 2] * nz))
        s2 = 1 - U**2
        return float(np.sum(np.abs(Wk)**2 * s2) / np.sum(s2))


# ------------------------------------------------------------------ states from the radial rule ---------------------
def radial_state(kappa, n_ring, a, Za, rmax=60.0):
    """level (omega, g(r), f(r) on r_i = (i+1/2) a) of partial wave kappa from the radial rule (zz_core.RadialWalk)"""
    rw = RadialWalk(kappa, a, rmax, Zalpha=Za)
    l = kappa if kappa > 0 else -kappa - 1
    guess = dirac_formula(n_ring, kappa, Za)          # only chooses which ring to listen to (as zz.rs)
    om, u, res = rw.eigen(guess - 0.002)
    g, f = rw.gf(u)
    return om, rw.r, g.real, f.real


def spin_angle(kappa, mj, X, Y, Z):
    """two-component spin-angle function Omega_{kappa, mj}(r_hat) on a grid (Condon-Shortley, l = kappa or -kappa-1)"""
    from scipy.special import sph_harm_y
    r = np.sqrt(X**2 + Y**2 + Z**2)
    th = np.arccos(np.clip(Z / np.maximum(r, 1e-12), -1, 1))
    phi = np.arctan2(Y, X)
    j = abs(kappa) - 0.5
    l = int(kappa if kappa > 0 else -kappa - 1)
    out = np.zeros((2,) + X.shape, complex)
    for s, ms in ((0, 0.5), (1, -0.5)):
        ml = mj - ms
        if abs(ml) > l:
            continue
        # Clebsch-Gordan <l ml 1/2 ms | j mj>
        if np.isclose(j, l + 0.5):
            c = np.sqrt((l + 0.5 + (mj if ms > 0 else -mj)) / (2 * l + 1))
        else:
            c = (-1 if ms > 0 else 1) * np.sqrt((l + 0.5 - (mj if ms > 0 else -mj)) / (2 * l + 1))
        out[s] = c * sph_harm_y(l, int(round(ml)), th, phi)
    return out


def dirac_3d_from_radial(r_grid, g, f, kappa, mj, x):
    """chiral spinor (L0, L1, R0, R1) on the cubic grid x^3 from radial (g, f):  psi = (g/r Omega_k, i f/r Omega_-k)"""
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    gi = np.interp(r, r_grid, g) / r
    fi = np.interp(r, r_grid, f) / r
    up = spin_angle(kappa, mj, X, Y, Z) * gi
    lo = 1j * spin_angle(-kappa, mj, X, Y, Z) * fi
    L = (up - lo) * S2
    R = (up + lo) * S2
    psi = np.concatenate([L, R]).astype(np.complex128)
    return psi / np.sqrt(norm2(psi))


def zero_cross_freq(t, y, t0=-np.inf, t1=np.inf):
    """frequency from zero crossings (as accoppiato.rs frequenza)"""
    m = (t >= t0) & (t <= t1)
    t, y = t[m], y[m]
    z = []
    for i in range(len(y) - 1):
        if (y[i] <= 0 < y[i + 1]) or (y[i] >= 0 > y[i + 1]):
            z.append(t[i] + (t[i + 1] - t[i]) * y[i] / (y[i] - y[i + 1]))
    if len(z) < 4:
        return np.nan
    return np.pi * (len(z) - 1) / (z[-1] - z[0])


def line_amp(t, y, om, t0, t1):
    """Hann-weighted amplitude of the component at om (as accoppiato.rs ampiezza_a)"""
    m = (t >= t0) & (t <= t1)
    t, y = t[m], y[m]
    h = np.sin(np.pi * (t - t[0]) / (t[-1] - t[0]))**2
    c = np.sum(h * y * np.exp(1j * om * t))
    return 2 * abs(c) / np.sum(h)


if __name__ == "__main__":
    import time
    rng = np.random.default_rng(SEED)
    n = 10
    # 1. Dirac kernel vs zz_core.Walk (Coulomb phase only)
    D = Dirac3D(n, 0.5, Z=60)
    V = -60 * ALPHA / D.r
    Wr = Walk("dirac", n, 0.5, theta=0.5, V=V)
    psi = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n))
    a, b = psi.copy(), psi.copy()
    for t in range(5):
        a = Wr.tick(a)
        D.tick(b, t)
    print("Dirac kernel vs zz_core walk: max diff %.2e" % np.max(np.abs(a - b)))
    # 2. light kernel vs zz_core.Walk('light')
    L = Light3D(n, 1.0, sponge=0, sigma=1.0)
    WL = Walk("light", n, 1.0)
    F = rng.normal(size=(3, n, n, n)) + 1j * rng.normal(size=(3, n, n, n))
    a = F.copy()
    L.F[...] = F
    for t in range(5):
        a = WL.tick(a)
        L.tick()
    print("light kernel vs zz_core walk: max diff %.2e" % np.max(np.abs(a - L.F)))
    # 3. field-phase coin: e E.r phase with a uniform field equals multiplying by exp(-i eE z a)
    eE = np.array([0.0, 0.0, 0.01])
    a, b = psi.copy(), psi.copy()
    D.tick(a, 0, eE)
    Wr.tick_count = 0
    Wr2 = Walk("dirac", n, 0.5, theta=0.5, V=V + eE[2] * D.x[None, None, :])
    b = Wr2.tick(b)
    print("field coin vs walk with V + eE z: max diff %.2e" % np.max(np.abs(a - b)))
    # 4. speed
    for nn in (48, 96):
        D = Dirac3D(nn, 0.5)
        p = np.ones((4, nn, nn, nn), complex)
        D.tick(p, 0)
        t0 = time.time()
        for t in range(3):
            D.tick(p, t, eE)
        print(f"Dirac {nn}^3 tick: {(time.time() - t0) / 3 * 1e3:.1f} ms")
    for nn in (96, 128):
        L = Light3D(nn, 4.0)
        L.tick()
        t0 = time.time()
        for t in range(3):
            L.tick()
            L.sponge()
        print(f"light {nn}^3 tick + sponge: {(time.time() - t0) / 3 * 1e3:.1f} ms")

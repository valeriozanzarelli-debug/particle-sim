"""zz_core.py - the ZIG-ZAG rule (frozen 2026-10-09, before any test was run).

One rule for every particle: every wave ball moves at c and rotates its internal state.
Discretization = split-step quantum walk / quantum cellular automaton on a cubic lattice
(Bialynicki-Birula PRD 49, 6920 (1994); Meyer J. Stat. Phys. 85, 551 (1996); Arrighi, Nesme, Forets
J. Phys. A 47, 465302 (2014); Di Molfetta, Strauch et al. on Dirac quantum walks).

Units: hbar = c = m_e = 1 (length = lambda_C = hbar/(m_e c) = 3.8616e-13 m, time = hbar/(m_e c^2),
energy = m_e c^2 = 510998.95 eV). Lattice spacing a, tick dt = a / c = a.

THE RULE (one tick):
 1. TRANSLATE. Along each axis j in turn (order x,y,z on even ticks, z,y,x on odd ticks) every ball
    splits its internal state on the eigenvectors of its spin matrix S_j:
        eigenvalue +1 -> moves one cell along +j,  eigenvalue -1 -> one cell along -j,  0 -> stays.
    spin-1/2 right-handed (R): S_j = sigma_j ; left-handed (L): S_j = -sigma_j  (Weyl walkers)
    spin-1 (light, Riemann-Silberstein F = E + i c B): (S_k)_ij = -i eps_kij      (Bialynicki-Birula)
    In k-space one axis is exactly exp(-i a k_j S_j): the ball always moves at c.
 2. ROTATE (the "flip"). After the translation, at every site
        L' = e^{-i V dt} (cos th L - i sin th R),   R' = e^{-i V dt} (cos th R - i sin th L)
    th = theta_m = m c^2 dt / hbar = m c a / hbar  (Higgs-like mixing of the two handednesses).
    theta_m = 0: massless. V = charge x potential (static Coulomb of a point nucleus, -Z alpha / r)
    applied as a phase per tick. Light has no flip (theta = 0, no charge).
 3. SURFER: v = c psi^dag alpha psi / psi^dag psi, alpha = diag(-sigma, +sigma) on (L, R):
        v = c (R^dag sigma R - L^dag sigma L) / (|L|^2 + |R|^2)    -> |v| <= c for any spinor.
INPUTS: c, hbar (units), the electron mass (only through theta_m = m c a / hbar), the charge
alpha = e^2/(4 pi eps0 hbar c) and Z of the nucleus. The lattice spacing a is NOT physics:
every number is reported for several a and extrapolated to a -> 0.

Continuum limit (second order, the alternated axis order makes two ticks palindromic):
   i d psi/dt = (c alpha . p + beta m c^2 + V) psi      (Dirac, chiral representation)
   i dF/dt = c curl F                                    (Maxwell in vacuum, F = E + icB)

Radial form (for hydrogen): the same rule written in spherical coordinates for a state of given
kappa (Dirac quantum number). Two radial walkers per site, chi+ (moves outwards at c) and chi- (moves
inwards at c); rotation = exp(-i dt (m sigma_z + kappa/r sigma_x + V)) acting on (g, f) = radial
large/small components. At r = 0 the inward walker turns into the outward one (g(0) = 0).
"""
import numpy as np

HBAR = 1.054571817e-34
C = 299792458.0
ME = 9.1093837015e-31
QE = 1.602176634e-19
ALPHA = 7.2973525693e-3
LAMBDA_C = HBAR / (ME * C)
T_UNIT = HBAR / (ME * C**2)
MC2_EV = ME * C**2 / QE
H_EVS = 4.135667696e-15
HC_EVNM = 1239.84198433
MP_ME = 1836.15267343
SEED = 20261009

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)
PAULI = [SX, SY, SZ]
EPS = np.zeros((3, 3, 3))
for (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    EPS[i, j, k], EPS[i, k, j] = 1, -1
SPIN1 = [-1j * EPS[k] for k in range(3)]          # (S_k)_{ij} = -i eps_{kij}


def eigbasis(S):
    """eigenvalues (rounded: +1, 0, -1) and eigenvector columns of a Hermitian spin matrix"""
    w, v = np.linalg.eigh(S)
    o = np.argsort(-w)
    return np.round(w[o]).astype(int), v[:, o]


BASIS_HALF = [eigbasis(S) for S in PAULI]
BASIS_ONE = [eigbasis(S) for S in SPIN1]


def order_of_tick(n):
    return (0, 1, 2) if n % 2 == 0 else (2, 1, 0)


# =============================== lattice walks ===============================
class Walk:
    """Cubic-lattice zig-zag walk in d = 1, 2 or 3 dimensions, periodic box, n sites per axis.
    kind = 'dirac' : psi shape (4, n,...,n) = (L0, L1, R0, R1)  (d = 3 only; or d = 2 with k_z = 0)
    kind = 'light' : psi shape (3, n,...,n) = F = (Fx, Fy, Fz)
    kind = 'dirac1': psi shape (2, n) = (R, L)  (1-D: R moves +x, L moves -x; = the 3-D walk for fields
                     constant in y, z and spinors along sigma_x)"""

    def __init__(self, kind, n, a, dim=3, theta=0.0, V=None, axes=None, axmap=None):
        self.kind, self.n, self.a, self.dim, self.theta = kind, n, a, dim, theta
        # physical axis j -> array axis of a component field (default: j)
        self.axmap = axmap if axmap is not None else {0: 0, 1: 1, 2: 2}
        self.V = V
        self.axes = axes if axes is not None else tuple(range(dim))   # lattice axes that exist
        x = (np.arange(n) - n / 2 + 0.5) * a
        self.x1 = x
        self.tick_count = 0
        if V is not None:
            self.phase = np.exp(-1j * V * a)
        else:
            self.phase = None

    # ---- one axis translation ----
    def _shift_block(self, blk, j, h, basis):
        lam, P = basis[j]
        Pd = P.conj().T
        c = np.tensordot(Pd, blk, axes=(1, 0))
        for q in range(len(lam)):
            s = h * lam[q]
            if s != 0:
                c[q] = np.roll(c[q], s, axis=self.axmap[j])
        return np.tensordot(P, c, axes=(1, 0))

    def tick(self, psi):
        if self.kind == "dirac1":
            R = np.roll(psi[0], 1)
            L = np.roll(psi[1], -1)
            ct, st = np.cos(self.theta), np.sin(self.theta)
            out = np.array([ct * R - 1j * st * L, ct * L - 1j * st * R])
            if self.phase is not None:
                out *= self.phase
            self.tick_count += 1
            return out
        order = [j for j in order_of_tick(self.tick_count) if j in self.axes]
        if self.kind == "light":
            for j in order:
                psi = self._shift_block(psi, j, +1, BASIS_ONE)
            self.tick_count += 1
            return psi
        L, R = psi[0:2], psi[2:4]
        for j in order:
            L = self._shift_block(L, j, -1, BASIS_HALF)
            R = self._shift_block(R, j, +1, BASIS_HALF)
        ct, st = np.cos(self.theta), np.sin(self.theta)
        out = np.concatenate([ct * L - 1j * st * R, ct * R - 1j * st * L])
        if self.phase is not None:
            out *= self.phase
        self.tick_count += 1
        return out


# =============================== Bloch analysis ===============================
def T_axis(kj, j, h, spin):
    """one-axis translation in k-space: exp(-i h a k_j S_j) (a k_j = kj)"""
    lam, P = (BASIS_HALF if spin == 0.5 else BASIS_ONE)[j]
    return P @ np.diag(np.exp(-1j * h * lam * kj)) @ P.conj().T


def bloch_two_ticks(ka, theta=0.0, spin=0.5):
    """two-tick walk matrix W(k) for dimensionless k a (3-vector). spin 0.5 -> 4x4 (L,R); 1 -> 3x3"""
    mats = []
    for n in (0, 1):
        if spin == 0.5:
            TL, TR = I2.copy(), I2.copy()
            for j in order_of_tick(n):
                TL = T_axis(ka[j], j, -1, 0.5) @ TL
                TR = T_axis(ka[j], j, +1, 0.5) @ TR
            T = np.zeros((4, 4), complex)
            T[:2, :2], T[2:, 2:] = TL, TR
            Cm = np.cos(theta) * np.eye(4) - 1j * np.sin(theta) * np.block([[0 * I2, I2], [I2, 0 * I2]])
            mats.append(Cm @ T)
        else:
            T = np.eye(3, dtype=complex)
            for j in order_of_tick(n):
                T = T_axis(ka[j], j, +1, 1) @ T
            mats.append(T)
    return mats[1] @ mats[0]


def bloch_batch(KA, theta=0.0, spin=0.5):
    """vectorized two-tick walk matrices for many k (KA shape (N,3), dimensionless k a)"""
    KA = np.atleast_2d(KA)
    N = KA.shape[0]

    def Tax(j, h, sp):
        lam, P = (BASIS_HALF if sp == 0.5 else BASIS_ONE)[j]
        ph = np.exp(-1j * h * lam[None, :] * KA[:, j:j + 1])            # (N, d)
        return np.einsum("ab,nb,cb->nac", P, ph, P.conj())
    mats = []
    for nt in (0, 1):
        if spin == 0.5:
            TL = np.broadcast_to(I2, (N, 2, 2)).copy()
            TR = TL.copy()
            for j in order_of_tick(nt):
                TL = Tax(j, -1, 0.5) @ TL
                TR = Tax(j, +1, 0.5) @ TR
            Wt = np.zeros((N, 4, 4), complex)
            ct, st = np.cos(theta), np.sin(theta)
            Wt[:, :2, :] = np.concatenate([ct * TL, -1j * st * TR], axis=2)
            Wt[:, 2:, :] = np.concatenate([-1j * st * TL, ct * TR], axis=2)
            mats.append(Wt)
        else:
            T = np.broadcast_to(np.eye(3, dtype=complex), (N, 3, 3)).copy()
            for j in order_of_tick(nt):
                T = Tax(j, +1, 1) @ T
            mats.append(T)
    return mats[1] @ mats[0]


def quasi_energies(ka, theta=0.0, spin=0.5):
    """eigen-frequencies omega (in units of 1/dt) of the two-tick walk, sorted, and eigenvectors"""
    W = bloch_two_ticks(ka, theta, spin)
    w, v = np.linalg.eig(W)
    om = -np.angle(w) / 2.0
    o = np.argsort(om)
    return om[o], v[:, o]


# =============================== radial walk ===============================
U_CHI = np.array([[1, 1], [1j, -1j]], complex) / np.sqrt(2)   # columns u+ (sigma_y=+1), u- (sigma_y=-1)


class RadialWalk:
    """The rule for one Dirac partial wave kappa in a Coulomb field (point nucleus, charge Z).
    Sites r_i = (i + 1/2) a, i = 0..N-1; state (2, N) = (chi+, chi-)."""

    def __init__(self, kappa, a, rmax, Zalpha=ALPHA, mass=1.0, refl0=-1.0, refl1=-1.0):
        self.kappa, self.a, self.Za, self.m = kappa, a, Zalpha, mass
        self.N = int(np.ceil(rmax / a))
        self.r = (np.arange(self.N) + 0.5) * a
        self.theta = mass * a
        self.refl0, self.refl1 = refl0, refl1
        dt = a
        r = self.r
        V = -Zalpha / r
        hx, hz = kappa / r, mass * np.ones_like(r)
        nrm = np.sqrt(hx**2 + hz**2)
        phi = dt * nrm
        cph, sph = np.cos(phi), np.sin(phi)
        e = np.exp(-1j * V * dt)
        # M = e (cos - i sin (nx sx + nz sz)) in (g,f) basis, then to chi basis: U^dag M U
        M = np.zeros((self.N, 2, 2), complex)
        M[:, 0, 0] = cph - 1j * sph * hz / nrm
        M[:, 1, 1] = cph + 1j * sph * hz / nrm
        M[:, 0, 1] = -1j * sph * hx / nrm
        M[:, 1, 0] = -1j * sph * hx / nrm
        M *= e[:, None, None]
        self.M = np.einsum("ab,nbc,cd->nad", U_CHI.conj().T, M, U_CHI)
        # half coin (for the palindromic state): same with dt/2
        phih = 0.5 * phi
        Mh = np.zeros((self.N, 2, 2), complex)
        Mh[:, 0, 0] = np.cos(phih) - 1j * np.sin(phih) * hz / nrm
        Mh[:, 1, 1] = np.cos(phih) + 1j * np.sin(phih) * hz / nrm
        Mh[:, 0, 1] = -1j * np.sin(phih) * hx / nrm
        Mh[:, 1, 0] = -1j * np.sin(phih) * hx / nrm
        Mh *= np.exp(-0.5j * V * dt)[:, None, None]
        self.Mh = np.einsum("ab,nbc,cd->nad", U_CHI.conj().T, Mh, U_CHI)

    def shift(self, chi):
        p, q = chi
        pn = np.empty_like(p)
        qn = np.empty_like(q)
        pn[1:] = p[:-1]
        pn[0] = self.refl0 * q[0]
        qn[:-1] = q[1:]
        qn[-1] = self.refl1 * p[-1]
        return pn, qn

    def tick(self, chi):
        p, q = self.shift(chi)
        M = self.M
        return np.array([M[:, 0, 0] * p + M[:, 0, 1] * q, M[:, 1, 0] * p + M[:, 1, 1] * q])

    def matrix(self):
        """sparse one-tick operator W = M S in interleaved ordering index = 2 i + s"""
        import scipy.sparse as sp
        N = self.N
        rows, cols, vals = [], [], []
        # shifted components: P_i = chi+_{i-1} (i>=1), refl0 chi-_0 (i=0); Q_i = chi-_{i+1}, refl1 chi+_{N-1}
        src_p = np.where(np.arange(N) >= 1, 2 * (np.arange(N) - 1), 1)        # index of source of P_i
        fac_p = np.where(np.arange(N) >= 1, 1.0, self.refl0)
        src_q = np.where(np.arange(N) <= N - 2, 2 * (np.arange(N) + 1) + 1, 2 * (N - 1))
        fac_q = np.where(np.arange(N) <= N - 2, 1.0, self.refl1)
        i = np.arange(N)
        for s in (0, 1):
            rows += [2 * i + s, 2 * i + s]
            cols += [src_p, src_q]
            vals += [self.M[:, s, 0] * fac_p, self.M[:, s, 1] * fac_q]
        return sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * N, 2 * N))

    def eigen(self, omega_guess, v0=None, iters=8, lu=None):
        """inverse iteration on W at lambda0 = exp(-i omega_guess a) = the infinite-time lock-in.
        returns omega (energy, units m c^2), eigenvector chi (2, N) at the palindromic point."""
        import scipy.sparse as sp
        import scipy.sparse.linalg as spl
        W = self.matrix()
        lam0 = np.exp(-1j * omega_guess * self.a)
        if lu is None:
            lu = spl.splu((W - lam0 * sp.identity(2 * self.N, format="csc")).tocsc())
        x = np.ones(2 * self.N, complex) if v0 is None else v0.copy()
        x /= np.linalg.norm(x)
        lam = lam0
        for _ in range(iters):
            y = lu.solve(x)
            x = y / np.linalg.norm(y)
            lam = np.vdot(x, W @ x)
        om = -np.angle(lam) / self.a
        resid = np.linalg.norm(W @ x - lam * x)
        chi = x.reshape(self.N, 2).T.copy()
        # palindromic state u = Mh^dag v
        Mh = self.Mh
        u = np.array([np.conj(Mh[:, 0, 0]) * chi[0] + np.conj(Mh[:, 1, 0]) * chi[1],
                      np.conj(Mh[:, 0, 1]) * chi[0] + np.conj(Mh[:, 1, 1]) * chi[1]])
        return om, u, resid

    def palindromic(self, chi):
        """state at the symmetric point of the tick: u = Mh^dag chi (chi = eigenvector of W = M S)"""
        Mh = self.Mh
        return np.array([np.conj(Mh[:, 0, 0]) * chi[0] + np.conj(Mh[:, 1, 0]) * chi[1],
                         np.conj(Mh[:, 0, 1]) * chi[0] + np.conj(Mh[:, 1, 1]) * chi[1]])

    def gf(self, chi):
        """large and small radial components (g, f) from (chi+, chi-), made real by a global phase"""
        g = (chi[0] + chi[1]) / np.sqrt(2)
        f = 1j * (chi[0] - chi[1]) / np.sqrt(2)
        k = np.argmax(np.abs(g))
        ph = np.exp(-1j * np.angle(g[k]))
        g, f = g * ph, f * ph
        nrm = np.sqrt(np.sum(np.abs(g)**2 + np.abs(f)**2) * self.a)
        return g / nrm, f / nrm


def dirac_formula(n, kappa, Za=ALPHA):
    """Sommerfeld-Dirac energy (units m c^2, including rest energy). CHECK ONLY, never an input."""
    k = abs(kappa)
    d = k - np.sqrt(k**2 - Za**2)
    return 1 / np.sqrt(1 + (Za / (n - d))**2)


def gaussian3(walk, sigma, spinor, k0=(0, 0, 0), center=(0, 0, 0)):
    x = walk.x1
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    g = np.exp(-((X - center[0])**2 + (Y - center[1])**2 + (Z - center[2])**2) / (4 * sigma**2)
               + 1j * (k0[0] * X + k0[1] * Y + k0[2] * Z))
    s = np.asarray(spinor, complex)
    s = s / np.linalg.norm(s)
    psi = s[:, None, None, None] * g[None]
    return psi / np.sqrt(np.sum(np.abs(psi)**2))


def dirac_velocity(psi):
    """surfer velocity field v = (R^dag sigma R - L^dag sigma L)/|psi|^2 on the lattice (3-D), and |psi|^2
    (bug fixed after the first zz_free run: v was returned without the division by |psi|^2)"""
    L, R = psi[0:2], psi[2:4]
    rho = np.sum(np.abs(psi)**2, axis=0)
    v = []
    for S in PAULI:
        jr = np.real(np.einsum("a...,ab,b...->...", R.conj(), S, R))
        jl = np.real(np.einsum("a...,ab,b...->...", L.conj(), S, L))
        v.append(jr - jl)
    return np.array(v) / np.maximum(rho, 1e-300), rho


if __name__ == "__main__":
    import time
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    # 1. unitarity of the lattice walks (norm after 20 ticks)
    for kind, nc, th in (("dirac", 4, 0.3), ("light", 3, 0.0)):
        W = Walk(kind, 16, 0.5, theta=th, V=(rng.random((16, 16, 16)) if kind == "dirac" else None))
        psi = rng.normal(size=(nc, 16, 16, 16)) + 1j * rng.normal(size=(nc, 16, 16, 16))
        n0 = np.sum(np.abs(psi)**2)
        for _ in range(20):
            psi = W.tick(psi)
        print(f"{kind}: |norm(20 ticks) - norm| / norm = {abs(np.sum(np.abs(psi)**2) - n0) / n0:.2e}")
    # 2. lattice walk == Bloch matrix (one plane wave)
    n, a = 16, 1.0
    kk = 2 * np.pi * np.array([1, 2, 3]) / n
    W = Walk("dirac", n, a, theta=0.2)
    x = np.arange(n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s0 = rng.normal(size=4) + 1j * rng.normal(size=4)
    pw = np.exp(1j * (kk[0] * X + kk[1] * Y + kk[2] * Z))
    psi = s0[:, None, None, None] * pw
    psi = W.tick(W.tick(psi))
    pred = bloch_two_ticks(kk, 0.2) @ s0
    print("plane wave: walk vs Bloch matrix", np.max(np.abs(psi[:, 0, 0, 0] - pred)))
    WL = Walk("light", n, a)
    s1 = rng.normal(size=3) + 0j
    F = WL.tick(WL.tick(s1[:, None, None, None] * pw))
    print("plane wave light: walk vs Bloch", np.max(np.abs(F[:, 0, 0, 0] - bloch_two_ticks(kk, 0, 1) @ s1)))
    # 3. radial walk unitary
    rw = RadialWalk(-1, 0.5, 300.0)
    Wm = rw.matrix()
    v = rng.normal(size=2 * rw.N) + 1j * rng.normal(size=2 * rw.N)
    print("radial walk: | |Wv| - |v| | / |v| =", abs(np.linalg.norm(Wm @ v) - np.linalg.norm(v)) / np.linalg.norm(v))
    chi = v.reshape(rw.N, 2).T
    print("radial walk: tick() == matrix", np.max(np.abs(rw.tick(chi) - (Wm @ v).reshape(rw.N, 2).T)))
    print("time %.1f s" % (time.time() - t0))

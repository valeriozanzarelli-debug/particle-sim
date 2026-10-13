"""locale_core.py - tools of the study "Regola locale pulita" (docs/REGOLA_LOCALE.md).

A generic translation-invariant SPLIT-STEP walk = a product of local factors. Everything is analysed through its
Bloch matrix W(k) (exact for a translation-invariant linear rule = what the running walk does to a plane wave),
plus a real-space runner used only for checks (unitarity, charge, Bloch == running walk).

Factors (all LOCAL, finite stencil):
  ('shift', v, M) : the component with eigenvalue lam of the Hermitian matrix M (integer spectrum) moves by lam*v
                    (v = integer vector in LATTICE coordinates).  k-space: exp(-i (k.v_real) M).
                    This is the zig-zag "corri a c" step along v.
  ('coin', C)     : the same unitary C on every site (the "poi ruota" step; the mass flip is one of these).
  ('link', v, Q, C): rotation of a LINK: the part Q of the state at x+v and the part (1-Q) at x rotate together with C.
                    = S_v^dag C S_v, S_v = (1-Q) + Q exp(-i k.v)   (Q = projector).  Local: radius |v|.
A period = list of factors applied left to right (first factor acts first). T_phys = physical duration of the
period in units a/c (fixed by the Weyl/Maxwell continuum limit: sum over shifts |v| n n^T = T_phys * I).

Units: c = 1, lattice length 1 (the spacing of the integer coordinates), hbar = 1. Mass enters only through the
flip angle theta = m c a / hbar (as in the frozen rule lab/zz_core.py, which is NOT modified here).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

PI = np.pi
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)
PAULI = np.array([SX, SY, SZ])
_EPS = np.zeros((3, 3, 3))
for (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    _EPS[i, j, k], _EPS[i, k, j] = 1, -1
SPIN1 = np.array([-1j * _EPS[k] for k in range(3)])          # (S_k)_ij = -i eps_kij  (zz_core: light F = E + icB)
XCHI = np.kron(SX, I2)                                       # L <-> R (chirality) flip on (L0, L1, R0, R1)
SEED = 20261013


def sdot(n, spin):
    """spin matrix along the unit vector n (3-vector): sigma.n (spin 1/2) or S.n (spin 1)"""
    n = np.asarray(n, float)
    S = PAULI if spin == 0.5 else SPIN1
    return np.einsum("i,ijk->jk", n, S)


# ------------------------------------------------------------------ lattices
LATTICES = {
    # name: (basis B (columns = primitive vectors in real space, 3x3), dim, sites per unit volume)
    "cubic": (np.eye(3), 3),
    "square": (np.eye(3), 2),
    "hex": (np.array([[1, 0.5, 0], [0, np.sqrt(3) / 2, 0], [0, 0, 1]]), 2),
}


class SplitWalk:
    """factors: list of tuples (see module doc). dimint: internal dimension. basis: 3x3 real-space basis of the
    integer coordinates. T_phys: physical time of the whole period (c = 1)."""

    def __init__(self, factors, dimint, basis=np.eye(3), T_phys=None, name="", dim=3):
        self.factors, self.n, self.B, self.name, self.dim = factors, dimint, np.asarray(basis, float), name, dim
        self.pre = []
        self.ivec = []
        for f in factors:
            self.ivec.append(np.asarray(f[1], float) if f[0] in ("shift", "link") else None)
            if f[0] == "shift":
                v = self.B @ np.asarray(f[1], float)
                w, P = np.linalg.eigh(f[2])
                lam = np.round(w).astype(int)
                assert np.allclose(w, lam, atol=1e-9), "shift matrix must have integer spectrum"
                self.pre.append(("shift", v, lam, P))
            elif f[0] == "coin":
                self.pre.append(("coin", np.asarray(f[1], complex)))
            elif f[0] == "link":
                v = self.B @ np.asarray(f[1], float)
                self.pre.append(("link", v, np.asarray(f[2], complex), np.asarray(f[3], complex)))
            else:
                raise ValueError(f[0])
        self.T = T_phys if T_phys is not None else self.frame_time()

    def frame_time(self):
        """T from sum |v| n n^T over the shift factors (must be isotropic in the used dimensions)"""
        A = np.zeros((3, 3))
        # the generator of a shift factor is (k.v) M with M = s.n, n = v/|v|  ->  |v| (k.n)(s.n): weight |v| n n^T
        for f, f0 in zip(self.pre, self.factors):
            if f[0] == "shift":
                v = f[1]
                # signed duration: the spin axis n is stored in the shift matrix; time = (v . n) with n = axis of M
                n = self.axis_of(f0[2])
                tau = float(v @ n) if n is not None else np.linalg.norm(v)
                A += tau * np.outer(n, n) if n is not None else np.outer(v, v) / np.linalg.norm(v)
        d = self.dim
        sub = A[:d, :d]
        T = np.trace(sub) / d
        self.frame_aniso = float(np.max(np.abs(sub - T * np.eye(d))) / T)
        return T

    @staticmethod
    def axis_of(M):
        """unit axis n of a spin matrix M = s.n (spin 1/2: 2x2; spin 1: 3x3; Dirac: R block of a 4x4)"""
        M = np.asarray(M)
        if M.shape == (4, 4):
            M = M[2:, 2:]
        S = PAULI if M.shape == (2, 2) else SPIN1 if M.shape == (3, 3) else None
        if S is None:
            return None
        nrm = np.real(np.einsum("ijk,kj->i", S, M)) / np.real(np.trace(S[0] @ S[0]))
        return nrm / np.linalg.norm(nrm)

    def bloch(self, K):
        """W(k) for K shape (N,3) real-space wave vectors -> (N, n, n). Product in time order."""
        K = np.atleast_2d(np.asarray(K, float))
        N = K.shape[0]
        W = np.broadcast_to(np.eye(self.n, dtype=complex), (N, self.n, self.n)).copy()
        for f in self.pre:
            if f[0] == "shift":
                _, v, lam, P = f
                ph = np.exp(-1j * np.outer(K @ v, lam))                  # (N, m)
                F = np.einsum("ab,nb,cb->nac", P, ph, P.conj())
            elif f[0] == "coin":
                F = f[1][None]
            else:
                _, v, Q, C = f
                e = np.exp(-1j * (K @ v))[:, None, None]
                S = (np.eye(self.n) - Q)[None] + Q[None] * e
                F = np.einsum("nba,bc,ncd->nad", S.conj(), C, S)       # S^dag C S
            W = F @ W
        return W

    def phases(self, K):
        """eigen-phases eps in (-pi, pi] (W v = e^{-i eps} v) and eigenvectors"""
        w, v = np.linalg.eig(self.bloch(K))
        return -np.angle(w), v


# ------------------------------------------------------------------ builders
def make_walk(kind, ticks, basis=np.eye(3), theta=0.0, dim=3, name="", wilson=None):
    """kind: 'weyl' (2 comps, R), 'dirac' (4 comps L,R + mass flip theta after every tick), 'light' (spin 1).
    ticks: list of ticks, each a list of integer direction vectors (lattice coords) in time order.
    wilson: None or dict(r=phi, vecs=[...], scheme='palin') -> local link rotations L(x) <-> R(x+v) after each tick."""
    B = np.asarray(basis, float)
    fac = []
    for it, tk in enumerate(ticks):
        for d in tk:
            m = 1.0
            if isinstance(d, tuple) and len(d) == 2 and np.ndim(d[1]) == 0:
                d, m = d                     # ('direction', integer multiple m): m < 0 = a BACKWARD step
            d0 = np.asarray(d, float)
            n = B @ d0
            n = n / np.linalg.norm(n)
            d = d0 * m
            if kind == "weyl":
                fac.append(("shift", d, sdot(n, 0.5)))
            elif kind == "light":
                fac.append(("shift", d, sdot(n, 1)))
            else:
                M = np.zeros((4, 4), complex)
                M[:2, :2] = -sdot(n, 0.5)
                M[2:, 2:] = sdot(n, 0.5)
                fac.append(("shift", d, M))
        if kind == "dirac":
            sym = wilson is not None and wilson.get("place", "sym") == "sym"
            if wilson is not None:
                wl = dict(wilson); wl["h"] = wilson["h"] / 2 if sym else wilson["h"]
                if wilson.get("alt", True) and it % 2 == 1:      # odd ticks: link sequence mirrored (kills the drift)
                    wl["vecs"] = [-np.asarray(v, float) for v in wilson["vecs"]]
                fac += wilson_factors(wl)
            if theta != 0.0:
                fac.append(("coin", np.cos(theta) * np.eye(4) - 1j * np.sin(theta) * XCHI))
            if sym:
                fac += wilson_factors(wl)
    nint = {"weyl": 2, "dirac": 4, "light": 3}[kind]
    return SplitWalk(fac, nint, B, name=name, dim=dim)


QR = np.diag([0, 0, 1, 1]).astype(complex)        # projector on the R block


YCHI = np.kron(SY, I2)                             # "twisted" chirality rotation (i gamma5-like mass direction)


def wilson_factors(wl):
    """LOCAL Wilson term (a rotation that depends on the neighbours): for every link vector v in wl['vecs'] and -v,
    the L part at x and the R part at x+v rotate into each other by the angle h = wl['h'] about the chirality axis
    wl['axis'] ('X' = same direction as the mass flip, 'Y' = twisted, orthogonal to it); palindromic order so that
    the second-order (BCH) terms cancel; then the on-site part is undone so that NOTHING happens at k = 0:
        prod_v [S_v^dag e^{-i h A/2} S_v] (palindrome)  x  e^{+i h n_links A}
    k-space generator to first order: -h sum_v (1 - cos k.v) A  = a Wilson mass, zero at k = 0, large at the corners.
    Radius = max |v|. Parameters: h (the Wilson strength 'r'), the link set, the axis."""
    h, vecs = wl["h"], [np.asarray(v, float) for v in wl["vecs"]]
    A = XCHI if wl.get("axis", "Y") == "X" else YCHI
    allv = []
    for v in vecs:
        allv += [v, -v]
    out = []
    C = np.cos(h / 2) * np.eye(4) - 1j * np.sin(h / 2) * A
    for v in allv + allv[::-1]:
        out.append(("link", v, QR, C))
    nl = len(allv)
    out.append(("coin", np.cos(h * nl) * np.eye(4) + 1j * np.sin(h * nl) * A))
    return out


CUBE_AXES = [np.array(e, float) for e in np.eye(3)]


def zigzag_ticks(dim=3):
    """the frozen rule: tick 0 order x,y,z ; tick 1 order z,y,x (lab/zz_core.order_of_tick)"""
    ax = CUBE_AXES[:dim]
    return [ax, ax[::-1]]


# ------------------------------------------------------------------ sampling
def fib_sphere(n, half=False):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = PI * (1 + 5 ** 0.5) * i
    p = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], -1)
    if half:
        p = p[p[:, 2] >= 0]
    return p


def dirs_for(dim, n=300):
    if dim == 2:
        a = np.linspace(0, PI / 2, 46) if n >= 46 else np.linspace(0, PI / 2, n)
        return np.stack([np.cos(a), np.sin(a), 0 * a], -1)
    return fib_sphere(2 * n, half=True)


# ------------------------------------------------------------------ light metrics
def light_band(walk, K, ref=None, steps=None, sign=+1):
    """frequency (per unit physical time) of the moving helicity-(+) band, UNWRAPPED by continuation along the ray
    0 -> K (the period can be long, so quasi-energies fold many times).  ref = (omega, vec) from a nearby point:
    then the band is picked by eigenvector overlap and unwrapped against ref (used for finite differences)."""
    K = np.atleast_2d(K)
    kn = np.linalg.norm(K, axis=1)
    if ref is not None:
        eps, v = walk.phases(K)
        ov = np.abs(np.einsum("nai,na->ni", v.conj(), ref[1])) / np.linalg.norm(v, axis=1)
        j = np.argmax(ov, axis=1)
        o = eps[np.arange(len(j)), j]
        w0 = ref[0] * walk.T
        o = w0 + ((o - w0) + PI) % (2 * PI) - PI
        return o / walk.T, v[np.arange(len(j)), :, j]
    M = steps or max(4, int(np.ceil(np.max(kn) * walk.T * 1.5 / 0.4)))
    s0 = 1e-3
    eps, v = walk.phases(K * s0)
    target = sign * kn * s0 * walk.T
    j = np.argmin(np.abs(eps - target[:, None]), axis=1)
    o = eps[np.arange(len(j)), j]
    vec = v[np.arange(len(j)), :, j]
    for s in np.linspace(s0, 1, M + 1)[1:]:
        o, vec = light_band(walk, K * s, ref=(o / walk.T, vec), sign=sign)
        o = o * walk.T
    return o / walk.T, vec


def band_derivs(walk, K, h=1e-4, dim=3):
    """omega, gradient and Hessian (dim x dim) of the tracked light band by central differences"""
    K = np.atleast_2d(K)
    w0, v0 = light_band(walk, K)
    fun = lambda Q: light_band(walk, Q, ref=(w0, v0))[0]
    g = np.zeros((len(K), dim))
    Hs = np.zeros((len(K), dim, dim))
    E = np.eye(3)
    for a in range(dim):
        wp, wm = fun(K + h * E[a]), fun(K - h * E[a])
        g[:, a] = (wp - wm) / (2 * h)
        Hs[:, a, a] = (wp - 2 * w0 + wm) / h**2
        for b in range(a + 1, dim):
            wpp = fun(K + h * (E[a] + E[b])); wmm = fun(K - h * (E[a] + E[b]))
            wpm = fun(K + h * (E[a] - E[b])); wmp = fun(K - h * (E[a] - E[b]))
            Hs[:, a, b] = Hs[:, b, a] = (wpp - wpm - wmp + wmm) / (4 * h * h)
    return w0, g, Hs


def light_metrics(walk, N, dirs):
    """at N cells per lambda (lattice length 1): errors over directions.
    phase: |omega/k - 1|; group: |v_g| - 1 and walk-off angle; diffraction: k * transverse curvature vs 1."""
    k = 2 * PI / N
    K = k * dirs
    d = walk.dim
    w, g, H = band_derivs(walk, K, h=min(1e-3, k * 1e-2), dim=d)
    nd = dirs[:, :d]
    phase = w / k - 1
    vg = np.linalg.norm(g, axis=1)
    walkoff = np.arccos(np.clip(np.sum(g * nd, axis=1) / vg, -1, 1))
    # transverse curvature: Hessian restricted to the plane orthogonal to k
    diff = []
    for i in range(len(K)):
        n = nd[i]
        if d == 2:
            t = np.array([[-n[1], n[0]]])
        else:
            a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
            t1 = np.cross(n, a); t1 /= np.linalg.norm(t1)
            t2 = np.cross(n, t1)
            t = np.stack([t1, t2])
        Ht = t @ H[i] @ t.T
        diff.append(np.linalg.eigvalsh(Ht) * k - 1)
    diff = np.array(diff)
    return dict(phase=phase, group=vg - 1, walkoff=walkoff, diffraction=diff)


def worst(m):
    return dict(phase=float(np.max(np.abs(m["phase"]))), group=float(np.max(np.abs(m["group"]))),
                walkoff=float(np.max(m["walkoff"])), diffraction=float(np.max(np.abs(m["diffraction"]))))


def n_for_tol(walk, dirs, tol=0.01, key="all", lo=3.0, hi=400.0):
    """smallest N (cells/lambda, lattice length units) where the worst error (of `key`, or of all three) <= tol.
    bisection in log N (errors decrease monotonically with N in the range used; checked on a grid by the caller)"""
    def err(N):
        wm = worst(light_metrics(walk, N, dirs))
        if key == "all":
            return max(wm["phase"], wm["group"], wm["diffraction"])
        return wm[key]
    if err(hi) > tol:
        return np.inf
    if err(lo) <= tol:
        return lo
    a, b = np.log(lo), np.log(hi)
    for _ in range(22):
        m = 0.5 * (a + b)
        if err(np.exp(m)) <= tol:
            b = m
        else:
            a = m
    return float(np.exp(b))


# ------------------------------------------------------------------ fermion metrics
def dirac_positive(walk, K):
    """positive per-time energies of the period operator (two largest of eps/T in (-pi/T, pi/T])"""
    eps, _ = walk.phases(K)
    om = np.sort(eps / walk.T, axis=1)
    return om


def su2_coords(U):
    """U (N,2,2) in SU(2) -> a (N,4) with U = a0 + i (a1 sx + a2 sy + a3 sz)"""
    a0 = np.real(np.trace(U, axis1=1, axis2=2)) / 2
    a = [np.real(np.trace(U @ s, axis1=1, axis2=2) / 2j) for s in PAULI]
    return np.stack([a0] + a, -1)


def su2_degree(walk, n=48):
    """degree of k -> U(k) in SU(2) = S^3 over the cube [-pi,pi)^3 (the integer-coordinate torus; for a walk on a
    sublattice this counts all the sheets): (1/2pi^2) Int det[a, da/dkx, da/dky, da/dkz] d^3k (midpoint rule,
    analytic-enough derivatives by central differences)."""
    g = (np.arange(n) + 0.5) / n * 2 * PI - PI
    KX, KY, KZ = np.meshgrid(g, g, g, indexing="ij")
    K = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], -1)
    h = 1e-5
    a = su2_coords(walk.bloch(K))
    da = []
    for q in range(3):
        e = np.zeros(3); e[q] = h
        da.append((su2_coords(walk.bloch(K + e)) - su2_coords(walk.bloch(K - e))) / (2 * h))
    M = np.stack([a] + da, -1)                     # (N, 4, 4)
    det = np.linalg.det(M)
    return float(np.sum(det) * (2 * PI / n) ** 3 / (2 * PI ** 2))


def find_nodes(walk, n=24, tol=1e-7):
    """points where the 2x2 period operator is +-1 (Weyl nodes): grid search + Newton-free polish (Nelder-Mead).
    returns list of (k, sign (+1 for U=+1, -1 for U=-1), chirality)"""
    from scipy.optimize import minimize
    g = np.arange(n) / n * 2 * PI - PI
    KX, KY, KZ = np.meshgrid(g, g, g, indexing="ij")
    K = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], -1)
    a = su2_coords(walk.bloch(K))
    res = []
    for s in (+1, -1):
        dist = np.linalg.norm(a[:, 1:], axis=1) + (1 - s * a[:, 0])
        cand = np.argsort(dist)[: 64]
        for c in cand:
            if dist[c] > 0.6:
                continue
            f = lambda k: float(np.sum(su2_coords(walk.bloch(k[None]))[0, 1:] ** 2) + (1 - s * su2_coords(walk.bloch(k[None]))[0, 0]))
            r = minimize(f, K[c], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-16, maxiter=4000))
            if r.fun < tol:
                k0 = (r.x + PI) % (2 * PI) - PI
                if not any(np.linalg.norm(((k0 - q[0]) + PI) % (2 * PI) - PI) < 1e-4 for q in res):
                    # chirality = sign det d(a1,a2,a3)/dk
                    h = 1e-6
                    J = np.zeros((3, 3))
                    for q in range(3):
                        e = np.zeros(3); e[q] = h
                        J[:, q] = (su2_coords(walk.bloch((k0 + e)[None]))[0, 1:] - su2_coords(walk.bloch((k0 - e)[None]))[0, 1:]) / (2 * h)
                    res.append((k0, s, int(np.sign(np.linalg.det(J))), float(np.linalg.det(J))))
    return res


# ------------------------------------------------------------------ real-space runner (checks only)
class RealWalk:
    """runs the same factor list on a periodic box of n^dim sites in LATTICE coordinates (psi shape (n_int, n, ..)).
    V: optional on-site phase field applied after every period (any on-site unitary; used for the charge check)."""

    def __init__(self, walk, n, dim=3, V=None):
        self.w, self.n, self.dim, self.V = walk, n, dim, V

    def _roll(self, f, v):
        for ax in range(self.dim):
            s = int(round(v[ax]))
            if s:
                f = np.roll(f, s, axis=ax)
        return f

    def period(self, psi):
        for f, iv in zip(self.w.pre, self.w.ivec):
            if f[0] == "shift":
                _, v, lam, P = f
                c = np.tensordot(P.conj().T, psi, axes=(1, 0))
                c = np.array([self._roll(c[q], lam[q] * iv) if lam[q] else c[q] for q in range(len(lam))])
                psi = np.tensordot(P, c, axes=(1, 0))
            elif f[0] == "coin":
                psi = np.tensordot(f[1], psi, axes=(1, 0))
            else:
                _, v, Q, C = f
                qi = np.real(np.diag(Q)) > 0.5
                psi = np.array([self._roll(psi[q], iv) if qi[q] else psi[q] for q in range(len(qi))])
                psi = np.tensordot(C, psi, axes=(1, 0))
                psi = np.array([self._roll(psi[q], -iv) if qi[q] else psi[q] for q in range(len(qi))])
        if self.V is not None:
            psi = psi * np.exp(-1j * self.V)[None]
        return psi


# ------------------------------------------------------------------ integer-step compositions (4th order)
S4_SEQ = (1, 1, 1, 1, -2, 1, 1, 1, 1)        # sum m^3 = 0, palindrome -> 4th order (Yoshida/Suzuki with integers)


def composed_walk(kind, seq, dim=3, theta=0.0, name=""):
    """composition of palindromic 2-tick blocks P(m) = C(m th/2) [shifts by m: x,y,z] C(m th) [shifts by m: z,y,x] C(m th/2)
    (the frozen 2-tick rule is P(1) up to a conjugation by C(th/2)). m < 0 = the block run BACKWARDS in time
    (balls move against their spin, the flip turns back)."""
    ax = CUBE_AXES[:dim]
    fac = []
    nint = {"weyl": 2, "dirac": 4, "light": 3}[kind]

    def coin(a):
        return ("coin", np.cos(a) * np.eye(4) - 1j * np.sin(a) * XCHI)

    def shifts(order, m):
        out = []
        for e in order:
            if kind == "weyl":
                M = sdot(e, 0.5)
            elif kind == "light":
                M = sdot(e, 1)
            else:
                M = np.zeros((4, 4), complex); M[:2, :2] = -sdot(e, 0.5); M[2:, 2:] = sdot(e, 0.5)
            out.append(("shift", e * m, M))
        return out

    for m in seq:
        if kind == "dirac":
            fac.append(coin(m * theta / 2))
        fac += shifts(ax, m)
        if kind == "dirac":
            fac.append(coin(m * theta))
        fac += shifts(ax[::-1], m)
        if kind == "dirac":
            fac.append(coin(m * theta / 2))
    w = SplitWalk(fac, nint, np.eye(3), name=name, dim=dim)
    return w


def group_velocities(walk, K, h=1e-6, deg_tol=1e-7):
    """all band group velocities |grad eps|/T at K: degenerate-subspace Hellmann-Feynman (dW/dk numerical, the
    derivative of the phase is the eigenvalue of the Hermitian matrix i V^dag dW V / lambda inside each cluster of
    degenerate eigenvalues). gap = distance of each eigenvalue to the nearest eigenvalue of ANOTHER cluster
    (velocities near band crossings are unreliable: the caller filters on gap)."""
    K = np.atleast_2d(K)
    W = walk.bloch(K)
    N, n = W.shape[0], W.shape[1]
    dWs = []
    for a in range(walk.dim):
        e = np.zeros(3); e[a] = h
        dWs.append((walk.bloch(K + e) - walk.bloch(K - e)) / (2 * h))
    vel = np.zeros((N, n)); gap = np.zeros((N, n))
    for i in range(N):
        lam, V = np.linalg.eig(W[i])
        Q, _ = np.linalg.qr(V)
        # cluster
        lab = -np.ones(n, int); c = 0
        for p in range(n):
            if lab[p] < 0:
                lab[np.abs(lam - lam[p]) < deg_tol] = c; c += 1
        for cl in range(c):
            idx = np.where(lab == cl)[0]
            # orthonormal basis of the cluster subspace
            B = np.linalg.qr(V[:, idx])[0]
            G = []
            for a in range(walk.dim):
                M = 1j * (B.conj().T @ dWs[a][i] @ B) / lam[idx[0]]
                G.append(0.5 * (M + M.conj().T))
            # common eigenbasis: diagonalize a generic combination, then read each component
            Mc = sum((0.37 + 0.11 * a) * G[a] for a in range(walk.dim))
            _, U = np.linalg.eigh(Mc)
            gv = np.stack([np.real(np.diag(U.conj().T @ G[a] @ U)) for a in range(walk.dim)], -1)
            vel[i, idx] = np.linalg.norm(gv, axis=1) / walk.T
            other = np.abs(lam[idx[0]] - lam[lab != cl])
            gap[i, idx] = other.min() if other.size else 9.0
    return vel, gap

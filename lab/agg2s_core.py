"""agg2s_core.py - "2s all'agglomeratore" (study agg2s_): the coarse hydrogen atom of AGGLOMERATORE COMPLETO
(lab/aggfull_*, imported read-only, never modified) with TWISTED BOUNDARY CONDITIONS, to take the box error apart.

Idea (geometric). The periodic box of side L turns the single atom into a CRYSTAL of its images. The coarse eigenstate
is then the Bloch state of that crystal at a given twist theta (psi(x + L e_a) = e^{i theta_a} psi(x)). Periodic box
= theta = 0 = the bottom of the s band: the 2s "bonds" with its 6 images and is pulled DOWN (the old +1.2 % at +-12 a0).
The isolated atom is the AVERAGE of the band over the twists (the on-site energy of the Wannier state). For a band
E(theta) = E_a + sum_R t_R cos(theta.R):
  * theta = (0,0,0) and (pi,pi,pi) averaged: cancels the 6 faces (R = L) and the 8 corners (R = sqrt3 L) exactly;
    edges (R = sqrt2 L) survive: suppressed by e^{-(sqrt2 - 1) kappa L} relative to the face term.
  * theta = (pi/2, pi/2, pi/2) alone (Baldereschi-like): cancels faces, edges and corners (needs complex arithmetic).
Both are real-arithmetic friendly in the first form (periodic and anti-periodic real states), used here.

Units hbar = c = m_e = 1 (lambda_C, mc^2), as aggfull_core.
"""
import os, sys, json, time, pickle
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import warnings
import numpy as np
warnings.filterwarnings("ignore", message=".*not reaching the requested tolerance.*")
warnings.filterwarnings("ignore", message=".*Exited.*")
import scipy.fft as sfft
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_core as ac
import aggfull_atom as aa

OUT = ac.OUT
TH = ac.TH_A
ALPHA = ac.zz.ALPHA
A0 = ac.A0
MC2 = ac.MC2
RY = ac.RY


def fftn(a):
    return sfft.fftn(a, axes=(-3, -2, -1), workers=1)


def ifftn(a):
    return sfft.ifftn(a, axes=(-3, -2, -1), workers=1)


# ------------------------------------------------------------------ fields (Galerkin, aggfull_core) with a cache
class CompactField:
    """the alias-free Galerkin field of aggfull_core stored by PARITY: the 8 sub-lattice fields are the samples of one field
    on the 2N grid of half cells, which is even about the nucleus (W(-x) = W(x)); only |p| = 0..N half cells per axis are kept
    ((N+1)^3 numbers instead of 8 N^3). V_s is rebuilt on the fly (a copy, no arithmetic)."""

    def __init__(self, N, h, x, Vget):
        self.N, self.h, self.x = N, h, x
        self.subs = [(sx, sy, sz) for sx in (0, 1) for sy in (0, 1) for sz in (0, 1)]
        self.W = None
        self.parity_err = 0.0
        for sub in self.subs:
            V = Vget(sub)
            dt = V.dtype
            if self.W is None:
                self.W = np.zeros((N + 1,) * 3, dt)
            idx = [np.abs(2 * np.arange(N) + s - N) for s in sub]        # |p| of each grid point of this sub-lattice
            # keep the samples with p >= 0 (their mirror images are the same numbers up to round-off)
            sel = [np.where((2 * np.arange(N) + s - N >= 0) | (2 * np.arange(N) + s - N == -N))[0] for s in sub]
            self.W[np.ix_(idx[0][sel[0]], idx[1][sel[1]], idx[2][sel[2]])] = V[np.ix_(sel[0], sel[1], sel[2])]
            del V
        for sub in self.subs[:2] + self.subs[-1:]:                         # parity check on 3 sub-lattices
            V = Vget(sub)
            self.parity_err = max(self.parity_err, float(np.max(np.abs(V - self.Vsub(sub))) / np.max(np.abs(V))))
            del V
        self.idx = {}

    def Vsub(self, sub):
        idx = [np.abs(2 * np.arange(self.N) + s - self.N) for s in sub]
        return self.W[np.ix_(idx[0], idx[1], idx[2])]


def galerkin_field(N, h, cache=True):
    """aggfull_core.GalerkinField, stored compactly (CompactField); reuses the aggfull_ cache file if present (read only),
    else builds and caches under agg2s_ (never commits .npz)"""
    names = [os.path.join(OUT, f"aggfull_field_N{N}_h{h:.6f}.npz"), os.path.join(OUT, f"agg2s_field_N{N}_h{h:.6f}.npz")]
    for fn in names:
        if os.path.exists(fn):
            d = np.load(fn)
            return CompactField(N, h, d["x"], lambda s: d[f"V{s[0]}{s[1]}{s[2]}"])
    t0 = time.time()
    f = ac.GalerkinField(N, h)
    if cache:
        np.savez(names[1], x=f.x, **{f"V{s[0]}{s[1]}{s[2]}": f.V[s] for s in f.subs})
    print(f"  Galerkin field N={N} built in {time.time() - t0:.0f} s", flush=True)
    cf = CompactField(N, h, f.x, lambda s: f.V[s])
    del f
    return cf


# ------------------------------------------------------------------ kinetic energies on a twisted grid
def kgrid_axes(N, h, tw):
    """physical momenta per axis on the twisted grid: k_m + theta_a / L, theta_a = pi tw_a"""
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
    L = N * h
    return [k1 + np.pi * t / L for t in tw]


def kin_schroedinger(N, h, tw, m=1.0):
    ka = kgrid_axes(N, h, tw)
    return (ka[0][:, None, None] ** 2 + ka[1][None, :, None] ** 2 + ka[2][None, None, :] ** 2) / (2 * m)


def kin_learned(rec, N, h, tw):
    """learned spin-blind generator per unit time at the twisted momenta, parity-even part, rest energy removed:
    g(k) = sum_j Re(s_j) e^{i k.j}, s_j = tr(K_j)/2 (Hermiticity: s_{-j} = conj s_j, so Re s_j is the even part).
    Twist: k -> k + delta (cell units delta_a = pi tw_a / N) = coefficients times e^{i delta.j}."""
    K, offs = rec["K"], np.asarray(rec["offs"])
    s = 0.5 * (K[0, 0] + K[1, 1])
    se = s.real.astype(complex)
    g0 = float(np.sum(se.real))
    delta = np.pi * np.asarray(tw, float) / N
    c = se * np.exp(1j * offs @ delta)
    A = np.zeros((N, N, N), complex)
    np.add.at(A, (offs[:, 0] % N, offs[:, 1] % N, offs[:, 2] % N), c)
    g = (ifftn(A) * N ** 3).real
    spin = float(np.max(np.abs(K[0, 1])) + np.max(np.abs(K[0, 0] - K[1, 1])))
    return (g - g0) / h, spin


# ------------------------------------------------------------------ the twisted atom
class TwistAtom:
    """H = T(k + theta/L) + gamma0 V_G on real Bloch states psi(x + L e_a) = (+-1) psi(x) (tw_a = 0: periodic,
    tw_a = 1: anti-periodic), parity sector even / oddz about the nucleus (at the grid point N//2)."""

    def __init__(self, T, field, tw=(0, 0, 0), gamma0=1.0, sym="even"):
        self.T, self.f, self.tw, self.g0, self.sym = T, field, tuple(int(t) for t in tw), gamma0, sym
        self.N = N = T.shape[0]
        self.real = True
        h = field.h
        n = np.arange(N) - N // 2
        self.ax = [np.exp(-1j * np.pi * t * n / N) for t in self.tw]            # u = psi e^{-i theta x / L}
        ka = kgrid_axes(N, h, self.tw)
        self.ph = [np.exp(1j * k * h / 2) for k in ka]                       # interpolation to the half-cell sub-lattice
        self.pre = (1.0 / (T - T.min() + 0.5 * ALPHA ** 2)).astype(np.float32)
        self.napply = 0

    def _to_u(self, v):
        """u = psi e^{-i theta x / L} (complex copy)"""
        u = v.astype(complex)
        if any(self.tw):
            u *= self.ax[0][:, None, None]
            u *= self.ax[1][None, :, None]
            u *= self.ax[2][None, None, :]
        return u

    def _from_u(self, u):
        """psi = Re(u e^{+i theta x / L}) (u is overwritten)"""
        if any(self.tw):
            u *= np.conj(self.ax[0])[:, None, None]
            u *= np.conj(self.ax[1])[None, :, None]
            u *= np.conj(self.ax[2])[None, None, :]
        return u.real.copy()

    def refl(self, v, axis):
        r = np.roll(np.flip(v, axis=axis), 1, axis=axis)
        if self.tw[axis]:
            sl = [slice(None)] * 3
            sl[axis] = 0
            r[tuple(sl)] *= -1.0                       # psi(+L/2) = -psi(-L/2) for an anti-periodic axis
        return r

    def project(self, v, inplace=False):
        if self.sym is None:
            return v
        if not inplace:
            v = v.copy()
        for ax in range(3):
            r = self.refl(v, ax)
            if ax == 2 and self.sym == "oddz":
                v -= r
            else:
                v += r
            del r
            v *= 0.5
        return v

    def apply_field_k(self, F, out, coef=1.0):
        """out += coef * FFT(V_G psi) given F = FFT(u) (twisted interpolation phases); one complex128 temporary
        (a complex64 temporary was tried in development: it shifts the 1s by +3e-6 eV, i.e. 2e-7 relative; not used)"""
        one = np.ones(self.N)
        c = coef / 8
        tmp = np.empty_like(F)
        for sub in self.f.subs:
            p = [self.ph[i] if sub[i] else one for i in range(3)]
            np.multiply(F, p[0][:, None, None], out=tmp)
            tmp *= p[1][None, :, None]
            tmp *= p[2][None, None, :]
            tmp = sfft.ifftn(tmp, overwrite_x=True, workers=1)
            tmp *= self.f.Vsub(sub) if hasattr(self.f, "Vsub") else self.f.V[sub]
            tmp = sfft.fftn(tmp, overwrite_x=True, workers=1)
            tmp *= (c * np.conj(p[0]))[:, None, None]
            tmp *= np.conj(p[1])[None, :, None]
            tmp *= np.conj(p[2])[None, None, :]
            out += tmp
        del tmp
        return out

    def H(self, v, field=True, kin=True, T=None):
        v = self.project(v)
        F = sfft.fftn(self._to_u(v), overwrite_x=True, workers=1)
        TT = self.T if T is None else T
        out = TT * F if kin else np.zeros_like(F)
        if field:
            self.apply_field_k(F, out, self.g0)
        del F
        o = self._from_u(sfft.ifftn(out, overwrite_x=True, workers=1))
        del out
        return self.project(o, inplace=True)

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
            F = sfft.fftn(self._to_u(self.project(v)), overwrite_x=True, workers=1)
            F *= self.pre
            return self.project(self._from_u(sfft.ifftn(F, overwrite_x=True, workers=1)), inplace=True)
        return self._apply(X, pre)

    def energy(self, v):
        f = v.reshape((self.N,) * 3)
        Hf = self.H(f)
        nrm = float(np.vdot(f, f).real)
        E = float(np.vdot(f, Hf).real / nrm)
        return E, float(np.linalg.norm(Hf - E * f) / np.sqrt(nrm))

    def parts(self, v):
        """<T> and <V> separately (for the term-by-term budget)"""
        f = v.reshape((self.N,) * 3)
        nrm = float(np.vdot(f, f).real)
        return (float(np.vdot(f, self.H(f, field=False)).real / nrm), float(np.vdot(f, self.H(f, kin=False)).real / nrm))


def guess(atom, kind):
    x = atom.f.x
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij", sparse=True)
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    g = {"1s": np.exp(-r / A0), "2s": (1 - r / (2 * A0)) * np.exp(-r / (2 * A0)), "2p": Z * np.exp(-r / (2 * A0))}[kind]
    g = atom.project(np.broadcast_to(g, (atom.N,) * 3).copy())
    return (g / np.linalg.norm(g)).ravel()


def solve_state(atom, x0, Y=None, tol_res=2e-8, budget=480, chunk=5, log=None):
    """LOBPCG (aggfull_atom.solve) in chunks until the residual is below tol_res or the time budget is spent"""
    t0 = time.time()
    it_tot = 0
    while True:
        v, _, it = aa.solve(atom, x0, Y=Y, maxiter=chunk)
        it_tot += it
        x0 = v
        E, r = atom.energy(v)
        if log:
            log(f"    it {it_tot}: E = {E * MC2:.9f} eV, residual {r:.2e} ({time.time() - t0:.0f} s)")
        if r < tol_res or time.time() - t0 > budget:
            return v, E, r, it_tot, r < tol_res


# ------------------------------------------------------------------ lean single-vector LOBPCG (memory: ~7 grid vectors)
def lobpcg1(atom, x0, Y=None, tol_res=2e-8, budget=480, maxit=400, log=None, every=4):
    """smallest eigenpair of atom.H (real symmetric) with x orthogonal to the columns of Y; preconditioner atom.P.
    Same algorithm as scipy's LOBPCG for one vector (Rayleigh-Ritz on [x, w, p]); all updates in place
    (memory: x, Ax, w, Aw, p, Ap)."""
    import scipy.linalg as sla
    from scipy.linalg.blas import daxpy
    t0 = time.time()
    sh = (atom.N,) * 3
    H = lambda v: atom.H(v.reshape(sh)).ravel()

    def orthY(v):
        if Y is not None:
            for j in range(Y.shape[1]):
                daxpy(Y[:, j], v, a=-float(Y[:, j] @ v))
        return v
    x = orthY(np.array(x0, float).ravel())
    x /= np.linalg.norm(x)
    Ax = H(x)
    p = Ap = None
    it = 0
    while True:
        lam = float(x @ Ax)
        r = Ax.copy()
        daxpy(x, r, a=-lam)
        rn = float(np.linalg.norm(r))
        it += 1
        if log and (it % every == 0):
            log(f"    it {it}: E = {lam * MC2:.10f} eV, residual {rn:.2e} ({time.time() - t0:.0f} s)")
        if rn < tol_res or time.time() - t0 > budget or it > maxit:
            del r
            break
        w = atom.P(r)
        del r
        w = orthY(w)
        daxpy(x, w, a=-float(x @ w))
        w /= np.linalg.norm(w)
        Aw = H(w)
        S = [x, w] + ([p] if p is not None else [])
        AS = [Ax, Aw] + ([Ap] if p is not None else [])
        n = len(S)
        GA = np.array([[S[i] @ AS[j] for j in range(n)] for i in range(n)])
        GB = np.array([[S[i] @ S[j] for j in range(n)] for i in range(n)])
        GA = 0.5 * (GA + GA.T)
        if n == 3 and np.linalg.cond(GB) > 1e10:
            n = 2
            GA, GB = GA[:2, :2], GB[:2, :2]
        del S, AS
        ev, V = sla.eigh(GA, GB)
        cv = V[:, 0]
        # new direction p <- c1 w + c2 p (in the buffer of w), A p likewise
        w *= cv[1]
        Aw *= cv[1]
        if n == 3:
            daxpy(p, w, a=cv[2])
            daxpy(Ap, Aw, a=cv[2])
        p, Ap = w, Aw
        del w, Aw
        x *= cv[0]
        daxpy(p, x, a=1.0)
        Ax *= cv[0]
        daxpy(Ap, Ax, a=1.0)
        nx = np.linalg.norm(x)
        x /= nx
        Ax /= nx
        npn = np.linalg.norm(p)
        p /= npn
        Ap /= npn
        if it % 10 == 0:                     # refresh A x (round-off drift of the recurrence)
            orthY(x)
            x /= np.linalg.norm(x)
            del Ax
            Ax = H(x)
    del Ax, p, Ap
    E, res = atom.energy(x)
    return x, E, res, it, res < tol_res


# ------------------------------------------------------------------ geometry of the Brillouin cube (Born tail)
def born_constants():
    """I = int_{outside cube [-pi/h, pi/h]^3} d^3k/((2pi)^3 k^4) = I_coef h,  J = int d^3k/((2pi)^3 k^6) = J_coef h^3.
    Spherical coordinates: the cube surface is at k = (pi/h)/max|n_i|, so I = (h/pi)(2pi)^-3 int dOmega max|n_i| and
    J = (h/pi)^3 (2pi)^-3 int dOmega max|n_i|^3 / 3. The angular averages by quadrature on the 1/48 wedge of the sphere."""
    from scipy import integrate
    # wedge: z >= x >= y >= 0 ; n = (sin t cos p, sin t sin p, cos t), 0 <= p <= pi/4, cos t >= sin t cos p
    def wedge(power):
        f = lambda t, p: np.cos(t) ** power * np.sin(t)
        tmax = lambda p: np.arctan(1 / np.cos(p))       # z >= x = sin t cos p
        val, err = integrate.dblquad(f, 0, np.pi / 4, 0, tmax, epsabs=1e-13, epsrel=1e-12)
        return 48 * val / (4 * np.pi)        # average over the sphere
    m1, m3 = wedge(1), wedge(3)
    I_coef = m1 * 4 * np.pi / (np.pi * (2 * np.pi) ** 3)
    J_coef = m3 * 4 * np.pi / (3 * np.pi ** 3 * (2 * np.pi) ** 3)
    return dict(mean_max=m1, mean_max3=m3, I_coef=I_coef, J_coef=J_coef)

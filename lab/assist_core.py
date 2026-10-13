"""assist_core.py - pairs from the vacuum in a time-dependent homogeneous field, 1+1 D, three ways.

Study: dynamically assisted Schwinger effect (docs/SCHWINGER_ASSISTITO.md). Method copied from lab/redo_r1.py (R1):
a homogeneous field E(t) keeps every sea plane wave a plane wave, so the whole Dirac sea evolves mode by mode
(one 2x2 problem per canonical momentum p). Units hbar = c = m = 1, e absorbed in the field: E_crit = m^2/e = 1.

1. THE RULE (zig-zag walk, 1+1 D, same as redo_r1.py / zz_core.py in one dimension). One tick of length a:
       R(x) <- e^{-i a A} R(x-a),  L(x) <- e^{+i a A} L(x+a)          (Peierls phase, A = vector potential)
       R' = cos th R - i sin th L,  L' = cos th L - i sin th R,       th = m a  (the flip = the mass)
   in k-space (k = p a): psi <- C diag(e^{-i(k + aA)}, e^{+i(k + aA)}) psi.  A is taken at the middle of the tick.
   Vacuum = lower Floquet band of the static walk at the start; pairs = weight in the upper band at the end.
   The doubler gap at k = pi is never reached: every mode used keeps |k + aA| < pi - margin (checked).
2. CONTINUUM REFERENCE (exact 1+1 D Dirac, what the quantum kinetic equation solves): i dpsi/dt = H psi with
       H(t) = (p + A(t)) sigma_z + m sigma_x
   integrated with the 4th-order Magnus (2 Gauss points) scheme, exact 2x2 exponentials. Independent of the walk.
3. SEMICLASSICS (WKB / Dykhne, = worldline instanton for time-dependent fields, Dumlu & Dunne 2010-11):
       f(p) ~ exp(-4 Im K),  K = int_{Re t*}^{t*} omega(t) dt,  omega^2 = m^2 + (p + A(t))^2,
   t* = complex turning point (omega = 0) in the upper half plane that minimises Im K.
Fields are given through A(t) (analytic in t), so E(t) = dA/dt switches off smoothly:
   strong slow : A1 = (E1/Om) tanh(Om t)                         [E1 sech^2(Om t)]
   weak Sauter : A2 = (E2/w) tanh(w t)                           [E2 sech^2(w t)]
   weak cos    : A2 = (E2/w) sin(w t) sech^2(Om t)               [~ E2 cos(w t) sech^2(Om t)]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from numba import njit

PI = np.pi
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


# ------------------------------------------------------------------ fields
class Field:
    """E1 strong slow sech^2 pulse (gamma_s = Om/E1) + optional weak pulse ('sauter' | 'cos' | None)."""

    def __init__(self, E1, Om, E2=0.0, w=1.0, shape=None, E1_shape="sauter"):
        self.E1, self.Om, self.E2, self.w, self.shape = E1, Om, E2, w, shape
        self.E1_shape = E1_shape          # 'sauter' (E1 sech^2) or 'none' (weak pulse alone)

    def A(self, t):
        t = np.asarray(t)
        a = np.zeros_like(t, dtype=complex if np.iscomplexobj(t) else float)
        if self.E1_shape == "sauter" and self.E1 != 0:
            a = a + self.E1 / self.Om * np.tanh(self.Om * t)
        if self.shape == "sauter" and self.E2 != 0:
            a = a + self.E2 / self.w * np.tanh(self.w * t)
        elif self.shape == "cos" and self.E2 != 0:
            a = a + self.E2 / self.w * np.sin(self.w * t) / np.cosh(self.Om * t) ** 2
        return a

    def dA(self, t):
        t = np.asarray(t)
        d = np.zeros_like(t, dtype=complex if np.iscomplexobj(t) else float)
        if self.E1_shape == "sauter" and self.E1 != 0:
            d = d + self.E1 / np.cosh(self.Om * t) ** 2
        if self.shape == "sauter" and self.E2 != 0:
            d = d + self.E2 / np.cosh(self.w * t) ** 2
        elif self.shape == "cos" and self.E2 != 0:
            ch = np.cosh(self.Om * t)
            d = d + self.E2 * (np.cos(self.w * t) / ch ** 2
                               - 2 * self.Om / self.w * np.sin(self.w * t) * np.tanh(self.Om * t) / ch ** 2)
        return d

    def T(self):
        """half time window: strong pulse tails down to ~e^-24"""
        if (self.E1_shape == "sauter" and self.E1 != 0) or self.shape == "cos":
            return 12.0 / self.Om
        return 14.0 / self.w

    def amax(self):
        return float(np.max(np.abs(self.A(np.linspace(-self.T(), self.T(), 20001)))))


# ------------------------------------------------------------------ 1. the rule
def walk_bands(k, th):
    """lower / upper Floquet eigenvectors of U(k) = C diag(e^{-ik}, e^{ik}), vectorised over k."""
    c, s = np.cos(th), np.sin(th)
    k = np.asarray(k, float)
    U = np.empty(k.shape + (2, 2), complex)
    U[..., 0, 0] = c * np.exp(-1j * k); U[..., 0, 1] = -1j * s * np.exp(1j * k)
    U[..., 1, 0] = -1j * s * np.exp(-1j * k); U[..., 1, 1] = c * np.exp(1j * k)
    lam, v = np.linalg.eig(U)
    w = -np.angle(lam)
    o = np.argsort(w, axis=-1)
    idx = np.arange(k.size)
    lo = v.reshape(-1, 2, 2)[idx, :, o.reshape(-1, 2)[:, 0]]
    up = v.reshape(-1, 2, 2)[idx, :, o.reshape(-1, 2)[:, 1]]
    return lo, up


@njit(cache=True)
def _walk_run(ka, Aa, th, lo, up):
    c, s = np.cos(th), np.sin(th)
    n = ka.size
    out = np.empty(n)
    for j in range(n):
        R = lo[j, 0]; L = lo[j, 1]
        for i in range(Aa.size):
            q = ka[j] + Aa[i]
            R = R * (np.cos(q) - 1j * np.sin(q))
            L = L * (np.cos(q) + 1j * np.sin(q))
            R, L = c * R - 1j * s * L, c * L - 1j * s * R
        amp = np.conj(up[j, 0]) * R + np.conj(up[j, 1]) * L
        out[j] = amp.real ** 2 + amp.imag ** 2
    return out


def walk_pairs(field, p, th):
    """pair number per mode f(p) from the zig-zag walk with lattice spacing a = th (m = 1)."""
    a = th
    T = field.T()
    nt = int(np.ceil(2 * T / a))
    tmid = -T + (np.arange(nt) + 0.5) * a
    Aa = a * field.A(tmid)
    p = np.atleast_1d(np.asarray(p, float))
    k0 = p * a + a * field.A(-T)
    k1 = p * a + a * field.A(-T + nt * a)
    kmax = np.max(np.abs(p[:, None] * a + Aa[None, ::max(1, nt // 2000)]))
    if kmax > PI - 0.3:
        raise ValueError(f"mode reaches the doubler region: |k|max = {kmax:.3f}")
    lo, _ = walk_bands(k0, th)
    _, up = walk_bands(k1, th)
    return _walk_run(p * a, Aa, th, lo, up)


# ------------------------------------------------------------------ 2. continuum Dirac
@njit(cache=True)
def _magnus_run(p, A1g, A2g, dt, m, v0, v1):
    n = p.size
    out = np.empty(n)
    r3 = np.sqrt(3.0) / 6.0
    for j in range(n):
        u0 = v0[j, 0]; u1 = v0[j, 1]
        for i in range(A1g.size):
            a1 = p[j] + A1g[i]; a2 = p[j] + A2g[i]
            nz = 0.5 * dt * (a1 + a2)
            nx = dt * m
            ny = r3 * dt * dt * (a2 - a1) * m
            nn = np.sqrt(nx * nx + ny * ny + nz * nz)
            cc = np.cos(nn); ss = np.sin(nn) / nn
            # exp(-i n.sigma) = cc - i ss (nx sx + ny sy + nz sz)
            m00 = cc - 1j * ss * nz
            m11 = cc + 1j * ss * nz
            m01 = -1j * ss * (nx - 1j * ny)
            m10 = -1j * ss * (nx + 1j * ny)
            u0, u1 = m00 * u0 + m01 * u1, m10 * u0 + m11 * u1
        amp = np.conj(v1[j, 0]) * u0 + np.conj(v1[j, 1]) * u1
        out[j] = amp.real ** 2 + amp.imag ** 2
    return out


def dirac_bands(q, m=1.0):
    """lower / upper eigenvectors of q sigma_z + m sigma_x"""
    q = np.atleast_1d(np.asarray(q, float))
    e = np.sqrt(q * q + m * m)
    up = np.stack([m * np.ones_like(q), e - q], 1); up /= np.linalg.norm(up, axis=1)[:, None]
    lo = np.stack([-(e - q), m * np.ones_like(q)], 1); lo /= np.linalg.norm(lo, axis=1)[:, None]
    return lo.astype(complex), up.astype(complex)


def dirac_pairs(field, p, dt=0.01, m=1.0):
    T = field.T()
    nt = int(np.ceil(2 * T / dt))
    dt = 2 * T / nt
    t0 = -T + np.arange(nt) * dt
    g = 0.5 - np.sqrt(3) / 6
    A1g = field.A(t0 + g * dt); A2g = field.A(t0 + (1 - g) * dt)
    p = np.atleast_1d(np.asarray(p, float))
    lo, _ = dirac_bands(p + field.A(-T), m)
    _, up = dirac_bands(p + field.A(T), m)
    return _magnus_run(p, A1g, A2g, dt, m, lo, up)


# ------------------------------------------------------------------ exact single Sauter pulse (Narozhny-Nikishov)
def sauter_pulse_exact(p, E, tau, m=1.0):
    """spinor pair number per mode for E(t) = E sech^2(t/tau), A = E tau tanh(t/tau), 1+1 D (m -> m_perp in 3+1)."""
    p = np.asarray(p, float)
    lam = E * tau * tau
    wm = np.sqrt(m * m + (p - E * tau) ** 2)   # t -> -inf
    wp = np.sqrt(m * m + (p + E * tau) ** 2)   # t -> +inf
    # f = [cosh(2 pi lam) - cosh(pi tau (wp - wm))] / [2 sinh(pi tau wp) sinh(pi tau wm)]  (log-safe)
    x1 = 2 * PI * lam; x2 = PI * tau * np.abs(wp - wm)
    num = np.exp(x1 - PI * tau * (wp + wm)) * (1 - np.exp(-(x1 - x2))) * (1 - np.exp(-(x1 + x2)))
    den = (1 - np.exp(-2 * PI * tau * wp)) * (1 - np.exp(-2 * PI * tau * wm))
    return num / den


# ------------------------------------------------------------------ 3. WKB / worldline instanton
def _imK(field, p, ts, m=1.0, n=3000):
    """Im of int_{Re ts}^{ts} omega dt along the vertical path, continuous branch, endpoint sqrt handled."""
    tr, ti = ts.real, ts.imag
    s = (np.arange(n) + 0.5) / n                     # midpoint in s, tau = ti (1 - (1-s)^2): dense at the end
    tau = ti * (1 - (1 - s) ** 2)
    dtau = ti * 2 * (1 - s) / n
    q = p + field.A(tr + 1j * tau)
    w2 = m * m + q * q
    w = np.sqrt(w2.astype(complex))
    # continuity of the branch, starting from the positive root on the real axis
    for i in range(1, n):
        if abs(w[i] + w[i - 1]) < abs(w[i] - w[i - 1]):
            w[i:] = -w[i:]
    if w[0].real < 0:
        w = -w
    K = np.sum(w * 1j * dtau)
    return K.imag


def turning_points(field, p, m=1.0, grid=24):
    """roots of p + A(t) = +- i m in the upper half plane (Newton from a grid), below the first pole line."""
    E1 = field.E1 if field.E1_shape == "sauter" else max(field.E2, 1e-9)
    tmax = 1.6 * m / E1
    imax = tmax
    if field.shape == "sauter":
        imax = min(imax, PI / (2 * field.w) * (1 - 1e-9))
    if field.E1_shape == "sauter":
        imax = min(imax, PI / (2 * field.Om) * 0.999)
    roots = []
    reach = 3 * m / E1 + abs(p) / E1
    ys = np.unique(np.r_[np.linspace(0.05, 0.95, grid // 2), 1 - np.geomspace(1e-7, 0.05, 8)]) * imax
    xs = np.r_[np.linspace(-reach, reach, grid), -p / E1]
    for sgn in (1, -1):
        for x0 in xs:
            for y0 in ys:
                t = complex(x0, y0)
                ok = False
                for it in range(80):
                    g = p + field.A(np.array([t]))[0] - sgn * 1j * m
                    d = field.dA(np.array([t]))[0]
                    if d == 0 or not np.isfinite(d):
                        break
                    dtn = g / d
                    lam = 1.0
                    while lam > 1e-6:
                        tn = t - lam * dtn
                        if 0 < tn.imag < imax and abs(tn.real) < 2 * reach:
                            break
                        lam *= 0.5
                    if lam <= 1e-6:
                        break
                    t = tn
                    if abs(lam * dtn) < 1e-13 * (1 + abs(t)) and lam == 1.0:
                        ok = True
                        break
                if ok and all(abs(t - r) > 1e-7 * (1 + abs(t)) for r in roots):
                    roots.append(t)
    return roots


def wkb_pairs(field, p, m=1.0, return_all=False):
    """f_WKB(p) = exp(-4 Im K) with the dominant turning point."""
    roots = turning_points(field, p, m)
    if not roots:
        return (0.0, []) if return_all else 0.0
    ims = [(_imK(field, p, t, m), t) for t in roots]
    ims = [x for x in ims if x[0] > 0]
    ims.sort(key=lambda x: x[0])
    f = float(np.exp(-4 * ims[0][0]))
    return (f, ims) if return_all else f


# ------------------------------------------------------------------ worldline-instanton closed forms (literature)
def s_wi_sauter(gamma):
    """Schuetzhold, Gies, Dunne PRL 101, 130404 (2008) eq. 15, divided by pi m^2/E: exponent / Schwinger exponent."""
    gamma = np.asarray(gamma, float)
    u = np.minimum(PI / (2 * gamma), 1.0)
    return 2 / PI * (np.arcsin(u) + u * np.sqrt(1 - u * u))


def s_linder_cos(gamma, eps, p=0.0):
    """Linder et al. PRD 92, 085009 (2015) eq. 37 (first order in eps), divided by the Schwinger exponent."""
    from scipy.special import iv
    return 1 - 2 * eps * np.cos(p * gamma) * iv(1, gamma) / gamma


def s_sauter_single(gs):
    """exponent of a single Sauter pulse E sech^2(Om t) at p = 0, / (pi m^2/E): 2/(1 + sqrt(1 + gs^2))."""
    return 2 / (1 + np.sqrt(1 + gs * gs))


def S_wi(gamma, gs, m=1.0):
    """worldline-instanton exponent for eps -> 0 (times E1/m^2), slow field E1 sech^2(Om t), gs = m Om / E1,
    weak Sauter pulse of combined Keldysh gamma = m w / E1: the instanton is cut at tau_p = pi/(2w)
    (Schuetzhold-Gies-Dunne 2008; for gs -> 0 this is s_wi_sauter * pi).  Units E1 = m = 1."""
    from scipy.integrate import quad
    tau_s = np.arctan(gs) / gs if gs > 0 else 1.0          # (1/gs) tan(gs tau_s) = 1
    tau_e = min(PI / (2 * gamma), tau_s)
    f = lambda tau: np.sqrt(max(0.0, 1 - (np.tan(gs * tau) / gs if gs > 0 else tau) ** 2))
    return 4 * quad(f, 0, tau_e, epsabs=1e-13, epsrel=1e-12, limit=200)[0]


# ------------------------------------------------------------------ Floquet (one oscillating field, p = 0)
@njit(cache=True)
def _walk_floquet(Aa, th):
    c, s = np.cos(th), np.sin(th)
    U = np.eye(2).astype(np.complex128)
    for i in range(Aa.size):
        q = Aa[i]
        e1 = np.cos(q) - 1j * np.sin(q); e2 = np.conj(e1)
        # tick = C diag(e1, e2)
        T00 = c * e1; T01 = -1j * s * e2; T10 = -1j * s * e1; T11 = c * e2
        U = np.array([[T00 * U[0, 0] + T01 * U[1, 0], T00 * U[0, 1] + T01 * U[1, 1]],
                      [T10 * U[0, 0] + T11 * U[1, 0], T10 * U[0, 1] + T11 * U[1, 1]]])
    return U


def floquet_split_walk(omega, xi, N, p=0.0, m=1.0):
    """zig-zag walk, A(t) = xi m sin(omega t), N ticks per period (a = 2 pi/(omega N), th = m a).
    Returns the Floquet quasi-energy splitting delta in [0, pi] (eigenphase difference of the period operator)."""
    a = 2 * PI / (omega * N)
    t = (np.arange(N) + 0.5) * a
    Aa = a * (p + xi * m * np.sin(omega * t))
    U = _walk_floquet(Aa, m * a)
    lam = np.linalg.eigvals(U)
    return float(abs(np.angle(lam[0] / lam[1])))


@njit(cache=True)
def _magnus_floquet(a1g, a2g, dt, m):
    r3 = np.sqrt(3.0) / 6.0
    U = np.eye(2).astype(np.complex128)
    for i in range(a1g.size):
        a1 = a1g[i]; a2 = a2g[i]
        nz = 0.5 * dt * (a1 + a2); nx = dt * m; ny = r3 * dt * dt * (a2 - a1) * m
        nn = np.sqrt(nx * nx + ny * ny + nz * nz)
        cc = np.cos(nn); ss = np.sin(nn) / nn
        m00 = cc - 1j * ss * nz; m11 = cc + 1j * ss * nz
        m01 = -1j * ss * (nx - 1j * ny); m10 = -1j * ss * (nx + 1j * ny)
        U = np.array([[m00 * U[0, 0] + m01 * U[1, 0], m00 * U[0, 1] + m01 * U[1, 1]],
                      [m10 * U[0, 0] + m11 * U[1, 0], m10 * U[0, 1] + m11 * U[1, 1]]])
    return U


def floquet_split_dirac(omega, xi, M=4000, p=0.0, m=1.0):
    P = 2 * PI / omega
    dt = P / M
    t0 = np.arange(M) * dt
    g = 0.5 - np.sqrt(3) / 6
    a1g = p + xi * m * np.sin(omega * (t0 + g * dt)); a2g = p + xi * m * np.sin(omega * (t0 + (1 - g) * dt))
    U = _magnus_floquet(a1g, a2g, dt, m)
    lam = np.linalg.eigvals(U)
    return float(abs(np.angle(lam[0] / lam[1])))


def resonance(n, xi, split, m=1.0, scan=400):
    """locate the n-photon resonance: minimum of the Floquet splitting near omega = 2 qbar / n.
    split(omega) -> delta.  Returns (omega_n, Delta_n = delta_min * omega / (2 pi))."""
    from scipy.optimize import minimize_scalar
    q0 = qbar(xi, m)
    lo, hi = 2 * q0 / (n + 0.45), 2 * q0 / (n - 0.45)
    ws = np.linspace(lo, hi, scan)
    d = np.array([split(w) for w in ws])
    i = int(np.argmin(d))
    a, b = ws[max(i - 1, 0)], ws[min(i + 1, scan - 1)]
    r = minimize_scalar(split, bounds=(a, b), method="bounded", options=dict(xatol=1e-14))
    return float(r.x), float(r.fun * r.x / (2 * PI))


def qbar(xi, m=1.0):
    """adiabatic dressed energy at p = 0: time average of sqrt(m^2 + A^2), A = xi m sin(phi)."""
    ph = np.linspace(0, 2 * PI, 20001)[:-1]
    return float(np.mean(np.sqrt(m * m + (xi * m * np.sin(ph)) ** 2)))

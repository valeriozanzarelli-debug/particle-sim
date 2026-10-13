"""supercrit_core.py - ATOMO SUPERCRITICO (docs/SUPERCRITICO.md): the frozen zig-zag radial rule
(lab/zz_core.py RadialWalk) with a general central potential V(r), plus the tools of the study:

  * SWalk: RadialWalk with V(r) given (point nucleus, cutoff/shell nucleus, uniform sphere). The rule is
    UNCHANGED: shift at c (chi+ out, chi- in, reflection at r = 0), then the coin
    exp(-i a (m sigma_z + kappa/r sigma_x + V)) on (g, f). For V = -Z alpha / r it reproduces RadialWalk
    to machine precision (checked in __main__).
  * level(): quasi-energy of the localized eigenstate of the one-tick operator (shift-invert on W).
  * phase_scan(): the stationary equation W chi = lambda chi solved as an outward recursion from r = 0
    (regular solution), decomposed at a far site on the two local Bloch waves; arg(c+/c-) = 2 delta + bg.
  * evolve(): numba time evolution with an absorbing sponge (non-unitary only inside the sponge).
  * Continuum Dirac reference (CHECK ONLY): bound state / critical charge / resonance phase shift from the
    radial Dirac ODE with the same V(r) (scipy solve_ivp).
Units hbar = c = m_e = 1 (length lambda_C = 386.159 fm, energy 510.999 keV).
"""
import os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import ALPHA, U_CHI, RadialWalk, dirac_formula, MC2_EV

LAMBDA_C_FM = 386.15926796
FM = 1.0 / LAMBDA_C_FM            # 1 fm in units of lambda_C


# ----------------------------------------------------------------- potentials
def V_point(r, Za):
    return -Za / r


def V_shell(r, Za, R):
    """'rectangular cutoff' = charge on a shell of radius R (Godunov-Machet-Vysotsky 2017 model)"""
    return np.where(r < R, -Za / R, -Za / r)


def V_sphere(r, Za, R):
    """uniformly charged sphere of radius R"""
    return np.where(r < R, -Za / (2 * R) * (3 - (r / R)**2), -Za / r)


def R_sphere_fm(Z, A=None):
    """DEBT: R = 1.2 A^(1/3) fm; for Z >~ 120 there is no measured nucleus, A = 2.5 Z (Popov's R_min(Z))"""
    if A is None:
        A = 2.5 * Z
    return 1.2 * A**(1.0 / 3.0)


# ----------------------------------------------------------------- the rule with V(r)
def coins(kappa, a, r, V, mass=1.0, frac=1.0):
    """chi-basis coin exp(-i frac a (m sz + kappa/r sx + V)) per site, shape (N, 2, 2) -- same algebra as zz_core"""
    dt = frac * a
    hx, hz = kappa / r, mass * np.ones_like(r)
    nrm = np.sqrt(hx**2 + hz**2)
    phi = dt * nrm
    c, s = np.cos(phi), np.sin(phi)
    M = np.zeros((len(r), 2, 2), complex)
    M[:, 0, 0] = c - 1j * s * hz / nrm
    M[:, 1, 1] = c + 1j * s * hz / nrm
    M[:, 0, 1] = -1j * s * hx / nrm
    M[:, 1, 0] = -1j * s * hx / nrm
    M *= np.exp(-1j * V * dt)[:, None, None]
    return np.einsum("ab,nbc,cd->nad", U_CHI.conj().T, M, U_CHI)


class SWalk(RadialWalk):
    """RadialWalk (inherits shift, tick, matrix, eigen, gf) with a general potential V(r) = Vfun(r)."""

    def __init__(self, kappa, a, rmax, Vfun, mass=1.0, refl0=-1.0, refl1=-1.0):
        self.kappa, self.a, self.m = kappa, a, mass
        self.Za = None
        self.N = int(np.ceil(rmax / a))
        self.r = (np.arange(self.N) + 0.5) * a
        self.theta = mass * a
        self.refl0, self.refl1 = refl0, refl1
        self.V = Vfun(self.r)
        self.M = coins(kappa, a, self.r, self.V, mass)
        self.Mh = coins(kappa, a, self.r, self.V, mass, frac=0.5)


def level(w, omega_guess, iters=12, region=None):
    """eigen quasi-energy closest to omega_guess (shift-invert), localization diagnostics.
    returns omega, density |chi|^2 per site (normalised), residual, staggering fraction"""
    W = w.matrix()
    lam0 = np.exp(-1j * omega_guess * w.a)
    lu = spl.splu((W - lam0 * sp.identity(2 * w.N, format="csc")).tocsc())
    x = np.ones(2 * w.N, complex)
    x /= np.linalg.norm(x)
    I = sp.identity(2 * w.N, format="csc")
    for rnd in range(3):
        for _ in range(iters):
            y = lu.solve(x)
            x = y / np.linalg.norm(y)
        lam = np.vdot(x, W @ x)
        if np.linalg.norm(W @ x - lam * x) < 1e-12 or rnd == 2:
            break
        lu = spl.splu((W - lam * I).tocsc())       # shift update (Rayleigh)
    om = -np.angle(lam) / w.a
    resid = np.linalg.norm(W @ x - lam * x)
    chi = x.reshape(w.N, 2).T
    rho = np.abs(chi[0])**2 + np.abs(chi[1])**2
    g = (chi[0] + chi[1]) / np.sqrt(2)
    # staggering: weight of the k ~ pi/a part of g (doubler copies are staggered)
    stag = np.sum(np.abs(g[1:] + g[:-1])**2) / max(np.sum(np.abs(g[1:] - g[:-1])**2) + np.sum(np.abs(g[1:] + g[:-1])**2), 1e-300)
    return om, rho, resid, 1.0 - stag


# ----------------------------------------------------------------- stationary recursion (phase scan)
def regular_solution(w, omega, i_stop=None):
    """outward recursion of W chi = lambda chi from r = 0. s_i = (chi+_{i-1}, chi-_i).
    (p_i, q_i) = lambda M_i^dag chi_i with p_i = chi+_{i-1} (p_0 = refl0 chi-_0), q_i = chi-_{i+1}.
    returns chip, chim arrays up to i_stop (normalised so that chi-_0 = 1); rescaled to avoid overflow."""
    lam = np.exp(-1j * omega * w.a)
    A = lam * np.conj(np.transpose(w.M, (0, 2, 1)))
    n = w.N if i_stop is None else i_stop
    return _recur(A[:n].copy(), w.refl0)


try:
    from numba import njit
except ImportError:            # pragma: no cover
    def njit(*a, **k):
        def d(f):
            return f
        return d if not (len(a) == 1 and callable(a[0])) else a[0]


@njit(cache=True)
def _recur(A, refl0):
    n = A.shape[0]
    chip = np.zeros(n, np.complex128)
    chim = np.zeros(n + 1, np.complex128)
    chim[0] = 1.0
    prev = refl0 * chim[0]                      # chi+_{-1}
    for i in range(n):
        chip[i] = (prev - A[i, 0, 1] * chim[i]) / A[i, 0, 0]
        chim[i + 1] = A[i, 1, 0] * chip[i] + A[i, 1, 1] * chim[i]
        prev = chip[i]
        s = abs(chip[i]) + abs(chim[i + 1])
        if s > 1e100:
            chip[:i + 1] /= s
            chim[:i + 2] /= s
            prev /= s
    return chip, chim[:n]


def local_transfer(w, omega, i):
    lam = np.exp(-1j * omega * w.a)
    A = lam * np.conj(w.M[i].T)
    return np.array([[1 / A[0, 0], -A[0, 1] / A[0, 0]],
                     [A[1, 0] / A[0, 0], A[1, 1] - A[1, 0] * A[0, 1] / A[0, 0]]])


def bloch_phase(w, omega, i_m):
    """arg(c1/c2) of the regular solution decomposed on the two local Bloch waves at site i_m (+1 site)"""
    chip, chim = regular_solution(w, omega, i_stop=i_m + 2)
    s = np.array([chip[i_m - 1], chim[i_m]])
    T = local_transfer(w, omega, i_m)
    mu, v = np.linalg.eig(T)
    # order by arg(mu) so the labelling is continuous in omega
    o = np.argsort(np.angle(mu))
    v = v[:, o]
    # fix eigenvector phases by a convention (first component real positive)
    v = v / (v[0] / np.abs(v[0]))[None, :]
    c = np.linalg.solve(v, s)
    return np.angle(c[1] / c[0]), np.abs(mu[o])


def phase_curve(w, omegas, i_m):
    ph = np.array([bloch_phase(w, om, i_m)[0] for om in omegas])
    return np.unwrap(ph)


def fit_breit_wigner(E, ph, sign=None):
    """fit ph(E) = b0 + b1 x + b2 x^2 + s*2*atan2(G/2, Er - E)  ->  (Er, G, residual rms)."""
    from scipy.optimize import least_squares
    d = np.gradient(ph, E)
    k = np.argmax(np.abs(d))
    s = np.sign(d[k]) if sign is None else sign
    # background slope estimate from the ends
    G0 = 4.0 / max(abs(d[k]), 1e-300)
    x0 = [ph[k] - s * np.pi, (ph[-1] - ph[0] - s * 2 * np.pi) / (E[-1] - E[0]), 0.0, E[k], G0]
    Ec = E[k]

    def res(p):
        b0, b1, b2, Er, G = p
        return b0 + b1 * (E - Ec) + b2 * (E - Ec)**2 + s * 2 * np.arctan2(abs(G) / 2, Er - E) - ph
    r = least_squares(res, x0, x_scale=[1, abs(x0[1]) + 1, 1, G0, G0], xtol=1e-14, ftol=1e-14)
    b0, b1, b2, Er, G = r.x
    return Er, abs(G), np.sqrt(np.mean(r.fun**2))


def resonance_scan(w, E_guess, i_m, G_guess=None, n=161, zooms=6):
    """locate the 2*pi jump of the Bloch phase near E_guess by successive zooms, then fit Breit-Wigner"""
    width = 0.3 if G_guess is None else max(20 * G_guess, 1e-9)
    Ec = E_guess
    for z in range(zooms):
        E = np.linspace(Ec - width, Ec + width, n)
        ph = phase_curve(w, E, i_m)
        d = np.abs(np.gradient(ph, E))
        # remove smooth background slope (median)
        dd = d - np.median(d)
        k = np.argmax(dd)
        Ec = E[k]
        Gest = 4.0 / max(dd[k], 1e-300)
        new = max(15 * Gest, 1e-10)
        if new > 0.5 * width:
            break
        width = new
    E = np.linspace(Ec - width, Ec + width, n)
    ph = phase_curve(w, E, i_m)
    Er, G, rms = fit_breit_wigner(E, ph)
    return Er, G, rms, (E, ph)


# ----------------------------------------------------------------- time evolution
@njit(cache=True)
def _evolve(chip, chim, M, refl0, refl1, mask, nticks, rin_idx, rec_every, rec):
    N = chip.shape[0]
    p = np.empty(N, np.complex128)
    q = np.empty(N, np.complex128)
    k = 0
    for t in range(nticks):
        p[1:] = chip[:-1]
        p[0] = refl0 * chim[0]
        q[:-1] = chim[1:]
        q[N - 1] = refl1 * chip[N - 1]
        for i in range(N):
            a = M[i, 0, 0] * p[i] + M[i, 0, 1] * q[i]
            b = M[i, 1, 0] * p[i] + M[i, 1, 1] * q[i]
            chip[i] = a * mask[i]
            chim[i] = b * mask[i]
        if (t + 1) % rec_every == 0:
            s_in = 0.0
            s_tot = 0.0
            for i in range(N):
                d = chip[i].real**2 + chip[i].imag**2 + chim[i].real**2 + chim[i].imag**2
                s_tot += d
                if i < rin_idx:
                    s_in += d
            rec[k, 0] = t + 1
            rec[k, 1] = s_in
            rec[k, 2] = s_tot
            k += 1
    return chip, chim


def sponge(r, r0, r1, strength=None, a=None):
    """smooth absorbing mask per tick: 1 for r < r0, decreasing as exp(-s ((r-r0)/(r1-r0))^2 a) beyond"""
    x = np.clip((r - r0) / (r1 - r0), 0, None)
    s = 3.0 if strength is None else strength
    return np.exp(-s * x**2 * (a if a is not None else 1.0))


def evolve(w, chi0, nticks, mask, r_in, rec_every=100):
    chip, chim = chi0[0].astype(complex).copy(), chi0[1].astype(complex).copy()
    rec = np.zeros((nticks // rec_every, 3))
    rin_idx = int(np.searchsorted(w.r, r_in))
    chip, chim = _evolve(chip, chim, w.M, w.refl0, w.refl1, mask, nticks, rin_idx, rec_every, rec)
    return np.array([chip, chim]), rec


# ----------------------------------------------------------------- continuum Dirac reference (CHECK ONLY)
def _rhs(kappa, E, Vf):
    def f(r, y):
        V = Vf(r)
        G, F = y
        return [-kappa / r * G + (E - V + 1) * F, -(E - V - 1) * G + kappa / r * F]
    return f


def dirac_regular(kappa, E, Vf, r_end, r0=1e-9, rtol=1e-11, dense=False, max_step=np.inf, t_eval=None):
    from scipy.integrate import solve_ivp
    # kappa = -1: G ~ r, F ~ (E - V0 - 1) r^2 / 3 near 0 (finite V0)
    V0 = Vf(r0)
    y0 = [r0, (E - V0 - 1) * r0**2 / 3.0]
    return solve_ivp(_rhs(kappa, E, Vf), (r0, r_end), y0, rtol=rtol, atol=1e-30, method="DOP853",
                     dense_output=dense, max_step=max_step, t_eval=t_eval)


def dirac_critical_Z(Vmaker, Zlo, Zhi, r_end=40.0, kappa=-1, tol=1e-6):
    """continuum Z_cr: at E = -1 the regular solution's G(r_end) changes sign when the 1s reaches -mc^2.
    Vmaker(Z) -> V(r) callable. Bisection between Zlo (subcritical) and Zhi."""
    def gend(Z):
        s = dirac_regular(kappa, -1.0, Vmaker(Z), r_end, max_step=0.05)
        return s.y[0, -1]
    glo = gend(Zlo)
    # scan for the FIRST sign change (higher states dive at larger Z)
    Zs = np.linspace(Zlo, Zhi, 41)
    for Zn in Zs[1:]:
        gn = gend(Zn)
        if np.sign(gn) != np.sign(glo):
            Zhi, ghi = Zn, gn
            break
        Zlo, glo = Zn, gn
    else:
        raise ValueError("no critical charge in range")
    while Zhi - Zlo > tol:
        Zm = 0.5 * (Zlo + Zhi)
        gm = gend(Zm)
        if np.sign(gm) == np.sign(glo):
            Zlo, glo = Zm, gm
        else:
            Zhi, ghi = Zm, gm
    return 0.5 * (Zlo + Zhi)


def dirac_bound(Vf, Elo, Ehi, r_end=None, kappa=-1, tol=1e-12):
    """continuum bound-state energy in the gap: sign change of G(r_end) of the regular solution."""
    def gend(E):
        re = r_end if r_end is not None else 30.0 / np.sqrt(max(1 - E * E, 1e-3))
        s = dirac_regular(kappa, E, Vf, min(re, 80.0), max_step=0.05)
        return s.y[0, -1]
    from scipy.optimize import brentq
    return brentq(gend, Elo, Ehi, xtol=tol)


def dirac_phase(Vf, E, r_m, kappa=-1):
    """continuum phase (in the lower continuum E < -1): tan of the angle of (G, F) at r_m mixed with
    the local free waves; returned as the angle of G + i c F with c = sqrt((1 - E')/(..)) normalised so that a free
    wave rotates uniformly. Smooth background, resonance = jump of 2 pi (delta jumps by pi)."""
    s = dirac_regular(kappa, E, Vf, r_m, max_step=0.02)
    G, F = s.y[0, -1], s.y[1, -1]
    V = Vf(r_m)
    Ee = E - V
    p = np.sqrt(max(Ee * Ee - 1, 1e-300))
    # free solution: G ~ sin(pr+..), F ~ (p/(Ee+1)) cos(..) -> angle of (G, (Ee+1)/p F)
    return np.angle(G * 1j + (Ee + 1) / p * F) * 2      # 2*(phase of the wave) -> jumps by 2 pi


def dirac_resonance(Vf, E_guess, r_m, width=0.3, n=161, zooms=8):
    Ec = E_guess
    w = width
    for z in range(zooms):
        E = np.linspace(Ec - w, Ec + w, n)
        ph = np.unwrap([dirac_phase(Vf, e, r_m) for e in E])
        d = np.abs(np.gradient(ph, E))
        dd = d - np.median(d)
        k = np.argmax(dd)
        Ec = E[k]
        Gest = 4.0 / max(dd[k], 1e-300)
        new = max(15 * Gest, 1e-10)
        if new > 0.5 * w:
            break
        w = new
    E = np.linspace(Ec - w, Ec + w, n)
    ph = np.unwrap([dirac_phase(Vf, e, r_m) for e in E])
    return fit_breit_wigner(E, ph)


if __name__ == "__main__":
    # identity with the frozen rule
    Z = 60
    rw = RadialWalk(-1, 0.05, 30.0, Zalpha=Z * ALPHA)
    sw = SWalk(-1, 0.05, 30.0, lambda r: V_point(r, Z * ALPHA))
    print("max |M - M_zz| =", np.max(np.abs(rw.M - sw.M)), " Mh:", np.max(np.abs(rw.Mh - sw.Mh)))
    om, rho, res, st = level(sw, dirac_formula(1, -1, Z * ALPHA))
    print("Z=60 a=0.05 omega_1s =", om, "Dirac", dirac_formula(1, -1, Z * ALPHA), "resid", res, "smooth", st)
    # recursion: at the eigen-omega the regular solution must decay (bounded) far out
    chip, chim = regular_solution(sw, om)
    print("regular solution at eigenvalue: |chi| at r=5, 15, 25:",
          [float(abs(chip[int(x / 0.05)])) for x in (5, 15, 25)])
    # continuum reference: point-like check with a tiny shell, Z = 60
    E = dirac_bound(lambda r: V_shell(r, Z * ALPHA, 1e-6), 0.5, 0.99)
    print("continuum shell R=1e-6 Z=60 E =", E, " Dirac point", dirac_formula(1, -1, Z * ALPHA))

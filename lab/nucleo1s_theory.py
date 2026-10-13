"""nucleo1s_theory.py - the analytic convergence law of the zig-zag radial rule for s states (study NUCLEO_1S).

For a state with Dirac quantum number kappa near a point nucleus, the relative binding error of the rule with cell a is
        e(a) = K a^(2 gamma) + B a^2 + (higher: a^(2 gamma + 1), a^4, ...),      gamma = sqrt(kappa^2 - (Z alpha)^2)
  * B  ("bulk"): the second-order effective Hamiltonian of the palindromic split step (symmetric BCH),
        H_eff = H + a^2/24 [A,[A,Bk]] - a^2/12 [Bk,[Bk,A]],  A = m sigma_z + kappa/r sigma_x + V,  Bk = sigma_y p,
    averaged on the exact Dirac 1s; the radial integrals int r^s e^{-2 lambda r} that diverge at r = 0 are continued
    analytically (Gamma(s+1)/(2 lambda)^(s+1)). No fit.
  * K  ("contact"): the cells next to the nucleus see the scale-free problem (kappa/r and Z alpha/r only, rotation angles
    kappa/(i+1/2) and Z alpha/(i+1/2) per tick, independent of a). Its stationary solution, far from the nucleus, is
    the regular solution r^gamma plus beta x the irregular one r^(-gamma) (in units of the cell). beta is one number
    of the rule (computed here by running the rule outward for 10^5-10^6 cells); the energy shift follows from the
    Wronskian of the Dirac radial equation:  dE = -2 gamma Z alpha c^2 beta a^(2 gamma),  c = coefficient of
    r^gamma (Z alpha, gamma - 1) in the normalised 1s.  No fit (beta is extrapolated in the cell index, not fitted to
    energies).
For Z alpha -> 0 (gamma -> 1) the two terms have opposite poles 1/(1 - gamma) that cancel, leaving a^2 (const - q ln a):
this is the "theta^2 ln theta" contact term seen for hydrogen in ZIGZAG.md section 5.4.
"""
import numpy as np
from scipy.special import gamma as Gf
from numba import njit
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_core import rot_chi, gamma_of, ALPHA  # noqa: F401


def dirac1s_coeffs(Za, m=1.0):
    """normalised Dirac 1s, (G, F) = r (g, f) = C r^gamma e^{-Za m r} (sqrt(1+gamma), -sqrt(1-gamma))"""
    g = gamma_of(-1, Za)
    lam = Za * m
    C = (2 * lam)**(g + 0.5) / np.sqrt(2 * Gf(2 * g + 1))
    return g, lam, C


@njit(cache=True)
def _rec(M, refl0, N):
    p = np.zeros(N, np.complex128)
    q = np.zeros(N + 1, np.complex128)
    q[0] = 1.0
    pin = refl0 * q[0]
    for i in range(N):
        q[i + 1] = (q[i] - M[i, 1, 0] * pin) / M[i, 1, 1]
        p[i] = M[i, 0, 0] * pin + M[i, 0, 1] * q[i + 1]
        pin = p[i]
    return p, q[:N]


def inner_solution(Za, kappa=-1, Nmax=400_000, refl0=-1.0, M_override=None):
    """the scale-free inner problem (cell units, m -> 0, omega -> 0): stationary amplitudes run outward from r = 0.
    M_override: dict {site index: 2x2 chi-basis rotation} to change single cells (decomposition studies)."""
    r = np.arange(Nmax) + 0.5
    M = rot_chi(np.ones(Nmax), kappa, r, Za, mass=0.0)
    if M_override:
        for i, Mi in M_override.items():
            M[i] = Mi
    p, q = _rec(M, refl0, Nmax)
    G = (p + q) / np.sqrt(2)
    F = 1j * (p - q) / np.sqrt(2)
    return r, G, F


def inner_beta(Za, kappa=-1, Nmax=400_000, refl0=-1.0, M_override=None, lo=100, hi=None, nterms=2, return_c1=False):
    """beta_inf: coefficient of the irregular solution, per unit regular one, in cell units (real part; the imaginary part
    is the half-tick offset of the post-rotation amplitudes and carries no energy)."""
    g = gamma_of(kappa, Za)
    r, G, F = inner_solution(Za, kappa, Nmax, refl0, M_override)
    vreg = np.array([Za, g + kappa])          # (G, F) direction of r^gamma   (s gamma: (s + kappa) G = Za F)
    virr = np.array([Za, -g + kappa])
    det = vreg[0] * virr[1] - vreg[1] * virr[0]
    a_ = (G * virr[1] - F * virr[0]) / det
    b_ = (vreg[0] * F - vreg[1] * G) / det
    beta = (b_ * r**g / (a_ / r**g)).real
    hi = hi or Nmax // 2
    idx = np.unique(np.logspace(np.log10(lo), np.log10(hi), 80).astype(int))
    x = r[idx]
    terms = [2 * g - 2, -1.0, 2 * g - 3, -2.0][:nterms]
    A = np.column_stack([np.ones_like(x)] + [x**t for t in terms])
    c, *_ = np.linalg.lstsq(A, beta[idx], rcond=None)
    return (c[0], c[1]) if return_c1 else c[0]


def contact_K(Za, beta=None, **kw):
    """relative 1s binding error coefficient of a^(2 gamma) (positive = too bound)"""
    g, lam, C = dirac1s_coeffs(Za)
    if beta is None:
        beta = inner_beta(Za, **kw)
    c1n = np.sqrt(1 + g) * C / Za
    dE = -2 * g * Za * c1n**2 * beta
    return -dE / (1 - g)


def bulk_B(Za, kappa=-1, m=1.0, parts=("all",)):
    """relative 1s binding error coefficient of a^2 from the symmetric BCH (analytic continuation).
    parts: which pieces of A to keep in the commutators: 'm' (mass), 'k' (kappa/r), 'V' (Coulomb), or 'all'."""
    g, lam, C = dirac1s_coeffs(Za, m)
    psi = [{g: np.sqrt(1 + g) * C}, {g: -np.sqrt(1 - g) * C}]
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], complex)

    def add(d, p, c):
        d[p] = d.get(p, 0) + c

    def mat(Mx, f):
        out = [{}, {}]
        for i in range(2):
            for j in range(2):
                if Mx[i, j] != 0:
                    for p, c in f[j].items():
                        add(out[i], p, Mx[i, j] * c)
        return out

    def rpow(q, c, f):
        return [{p + q: c * v for p, v in comp.items()} for comp in f]

    def plus(*fs):
        out = [{}, {}]
        for f in fs:
            for i in range(2):
                for p, c in f[i].items():
                    add(out[i], p, c)
        return out

    def deriv(f):
        out = [{}, {}]
        for i in range(2):
            for p, c in f[i].items():
                add(out[i], p - 1, p * c)
                add(out[i], p, -lam * c)
        return out

    use = set(parts)
    allp = "all" in use

    def Aop(f):
        terms = []
        if allp or "m" in use:
            terms.append(mat(m * sz, f))
        if allp or "k" in use:
            terms.append(mat(sx, rpow(-1, kappa, f)))
        if allp or "V" in use:
            terms.append(rpow(-1, -Za, f))
        return plus(*terms)

    def Bop(f):
        return mat(sy, [{p: -1j * c for p, c in comp.items()} for comp in deriv(f)])

    def sc(c, f):
        return [{p: c * v for p, v in comp.items()} for comp in f]

    AAB = plus(Aop(Aop(Bop(psi))), sc(-2, Aop(Bop(Aop(psi)))), Bop(Aop(Aop(psi))))
    BBA = plus(Bop(Bop(Aop(psi))), sc(-2, Bop(Aop(Bop(psi)))), Aop(Bop(Bop(psi))))
    dH = plus(sc(1 / 24, AAB), sc(-1 / 12, BBA))
    val = 0
    for i in range(2):
        for p1, c1 in psi[i].items():
            for p2, c2 in dH[i].items():
                P = p1 + p2
                val += np.conj(c1) * c2 * Gf(P + 1) / (2 * lam)**(P + 1)
    return -val.real / (1 - g)


def law(a, Za, K=None, B=None):
    g = gamma_of(-1, Za)
    K = contact_K(Za) if K is None else K
    B = bulk_B(Za) if B is None else B
    return K * a**(2 * g) + B * a**2


if __name__ == "__main__":
    for Z in (1, 20, 40, 60, 80, 92):
        Za = Z * ALPHA
        b, c1 = inner_beta(Za, return_c1=True)
        print(f"Z={Z}: gamma={gamma_of(-1, Za):.6f} beta={b:.6f} c1={c1:.6f} K={contact_K(Za, beta=b):.5f} B={bulk_B(Za):.5f}")

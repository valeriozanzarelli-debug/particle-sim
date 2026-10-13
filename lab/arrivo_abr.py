"""arrivo_abr.py - REFERENCE: Tumulka's absorbing boundary rule (ABR) for the Das-Duerr state (NR, 1-D).

This is an exam/comparison tool, not the rule. Schroedinger (hbar = m = d = 1) on 0 < z < L with
Dirichlet at z = 0 (closed end of the guide) and the absorbing boundary condition at the detector
    d psi / dz = i kappa psi     at z = L       (Tumulka, arXiv:1601.03715; J. Stat. Phys. 2022)
Detection-time density  rho_T(t) = (hbar kappa / m) |psi(L, t)|^2  (= flux at L, always >= 0).
For the waveguide the transverse factor drops out (spinless/scalar ABC: no spin dependence; see
Jozani & Tumulka arXiv:2603.22044).
A plane wave e^{ikz} is reflected with |(k - kappa)/(k + kappa)|^2; kappa = k absorbs it completely.
Tumulka's ideal Dirac ABCD (n.alpha psi = psi) = semi-ideal with theta = 0  ->  kappa = 2 m c / hbar
(arXiv:1601.04571 Eq. 66): in DD units kappa = 2 d / lambda_C.

Two solvers: Crank-Nicolson (any kappa) and the modal formula for kappa >> k (cavity modes sin(n pi z/L)
with decay rates Gamma_n = 2 k_n^2 / (kappa L) to first order in 1/kappa).
"""
import numpy as np
from scipy.linalg import solve_banded


def psi0_grid(z):
    return np.where(z < 1.0, np.sqrt(2.0) * np.sin(np.pi * np.clip(z, 0, 1)), 0.0)


def cn_abr(L, kappa, t_end, h=0.005, dt=5e-4, rec_every=1):
    """Crank-Nicolson with the Robin (ABR) condition at z = L. returns t, rho_T(t), norm(t)."""
    M = int(round(L / h))
    z = np.arange(1, M + 1) * h          # unknowns psi_1..psi_M (psi_0 = 0)
    psi = psi0_grid(z).astype(complex)
    psi /= np.sqrt(np.sum(np.abs(psi)**2) * h)
    # H = -1/2 Lap; last row with ghost psi_{M+1} = psi_{M-1} + 2 i h kappa psi_M
    main = np.full(M, 1.0 / h**2, dtype=complex)
    up = np.full(M - 1, -0.5 / h**2)
    lo = np.full(M - 1, -0.5 / h**2)
    main[-1] = (1.0 - 1j * h * kappa) / h**2
    lo[-1] = -1.0 / h**2
    # A = 1 + i dt/2 H ; B = 1 - i dt/2 H
    ab = np.zeros((3, M), complex)
    ab[0, 1:] = 0.5j * dt * up
    ab[1, :] = 1 + 0.5j * dt * main
    ab[2, :-1] = 0.5j * dt * lo
    bm, bu, bl = 1 - 0.5j * dt * main, -0.5j * dt * up, -0.5j * dt * lo
    nsteps = int(round(t_end / dt))
    ts, rs, ns = [], [], []
    for n in range(nsteps):
        rhs = bm * psi
        rhs[:-1] += bu * psi[1:]
        rhs[1:] += bl * psi[:-1]
        psi = solve_banded((1, 1), ab, rhs)
        if (n + 1) % rec_every == 0:
            ts.append((n + 1) * dt)
            rs.append(kappa * abs(psi[-1])**2)
            ns.append(np.sum(np.abs(psi[:-1])**2) * h + 0.5 * abs(psi[-1])**2 * h)
    return np.array(ts), np.array(rs), np.array(ns)


def modal_cn(L, n_max):
    """expansion coefficients of sqrt2 sin(pi z) (0<z<1) on sin(n pi z / L), 0<z<L (Dirichlet both ends).
    returns k_n, weights w_n = |c_n|^2 L/2 (sum w_n = 1)."""
    n = np.arange(1, n_max + 1)
    k = n * np.pi / L
    # int_0^1 sqrt2 sin(pi z) sin(k z) dz
    with np.errstate(divide="ignore", invalid="ignore"):
        I = np.sqrt(2.0) * 0.5 * (np.sin(np.pi - k) / (np.pi - k) - np.sin(np.pi + k) / (np.pi + k))
    I = np.where(np.isclose(k, np.pi), np.sqrt(2.0) * 0.5, I)
    c = (2.0 / L) * I
    w = np.abs(c)**2 * L / 2
    return k, c, w


def modal_large_kappa(L, kappa, n_max=200000):
    """first order in 1/kappa: Gamma_n = 2 k_n^2/(kappa L); incoherent detection CDF and mean time."""
    k, c, w = modal_cn(L, n_max)
    G = 2 * k**2 / (kappa * L)
    mean = np.sum(w / G)
    return dict(k=k, w=w, Gamma=G, mean=mean, wsum=np.sum(w))


def modal_cdf(mod, t):
    t = np.atleast_1d(t)
    return np.array([np.sum(mod["w"] * (1 - np.exp(-mod["Gamma"] * tt))) for tt in t])


def modal_quantile(mod, q, t_lo=1e-3, t_hi=1e30):
    """t with CDF(t) = q (bisection in log t)"""
    lo, hi = np.log(t_lo), np.log(t_hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if modal_cdf(mod, np.exp(mid))[0] < q:
            lo = mid
        else:
            hi = mid
    return np.exp(0.5 * (lo + hi))

"""aggmol_exact.py - EXAM ONLY: the exact non-relativistic one-electron two-centre problem (H2+ with clamped protons),
solved by separation of variables in prolate spheroidal coordinates (a textbook mathematical fact, no chemistry code,
no tables). Used only as the answer key of the two-centre test of the agglomerated rule (study aggmol_).

Coordinates: xi = (r_a + r_b)/R in [1, inf), eta = (r_a - r_b)/R in [-1, 1], phi.  psi = X(xi) Y(eta) e^{i m phi}.
With p^2 = -E R^2 / 2 (E = electronic energy, Eh, atomic units), the two separated equations (Z_a = Z_b = 1) read
   (xi^2-1) u'' + 2(m+1) xi u' + (m(m+1) + 2 R xi - p^2 xi^2) u = A u,        X = (xi^2-1)^{m/2} u   (radial)
   (1-eta^2) v'' - 2(m+1) eta v' + (p^2 eta^2 - m(m+1)) v = -A v,              Y = (1-eta^2)^{m/2} v (angular)
Both are solved by Chebyshev collocation (radial on xi - 1 in [0, S] with a quadratic map, Dirichlet at S; angular on
[-1, 1], regular solutions picked by collocating at the end points). The energy is the root p of A_rad(p) = A_ang(p)
for the chosen quantum numbers (n_xi radial nodes, n_eta angular nodes, m).
Check: R = 2, 1 sigma_g -> E = -1.1026342144949 Eh (Bates, Ledsham & Stewart 1953 and every later paper)."""
import numpy as np
from scipy.optimize import brentq


def cheb(N):
    """Chebyshev-Lobatto differentiation matrix (Trefethen), points x_j = cos(pi j / N)"""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(1))
    return D, x


def radial_A(p, R, m=0, N=90, S=None):
    """eigenvalues A of the radial operator (largest first) at given p"""
    S = S if S is not None else 46.0 / p
    D, x = cheb(N)
    t = (1 - x) / 2                       # t in [0, 1], t = 0 at x = 1
    s = S * t ** 2                        # xi - 1, clustered near xi = 1
    dsdt = 2 * S * t
    # d/ds = (1/(ds/dt)) d/dt, d/dt = -2 d/dx ; avoid t = 0 singularity of the map: use the map s = S t^2 only away
    # from 0 by switching to s = S (t^2 + t) / 2 (smooth, monotone, ds/dt > 0 everywhere)
    s = S * (t ** 2 + t) / 2
    dsdt = S * (2 * t + 1) / 2
    d2sdt2 = S * np.ones_like(t)
    Dt = -2 * D
    Ds = np.diag(1 / dsdt) @ Dt
    Ds2 = np.diag(1 / dsdt ** 2) @ (Dt @ Dt) - np.diag(d2sdt2 / dsdt ** 3) @ Dt
    xi = 1 + s
    L = np.diag(xi ** 2 - 1) @ Ds2 + np.diag(2 * (m + 1) * xi) @ Ds + np.diag(m * (m + 1) + 2 * R * xi - p ** 2 * xi ** 2)
    # Dirichlet at s = S (index of t = 1 is x = -1: last point)
    keep = np.arange(N)                   # drop the last point (u(S) = 0)
    ev = np.linalg.eigvals(L[np.ix_(keep, keep)])
    ev = ev[np.abs(ev.imag) < 1e-6 * (1 + np.abs(ev.real))].real
    return np.sort(ev)[::-1]


def angular_A(p, m=0, N=40):
    """separation constants A from the angular equation, smallest first (each A = - eigenvalue)"""
    D, x = cheb(N)
    L = np.diag(1 - x ** 2) @ (D @ D) - np.diag(2 * (m + 1) * x) @ D + np.diag(p ** 2 * x ** 2 - m * (m + 1))
    ev = np.linalg.eigvals(L)
    ev = ev[np.abs(ev.imag) < 1e-6 * (1 + np.abs(ev.real))].real
    return np.sort(-ev)                   # A = -eigenvalue, smallest first


def energy(R, n_xi=0, n_eta=0, m=0, N=90):
    """electronic energy (Eh) of the state with n_xi radial nodes, n_eta angular nodes, |m|"""
    def f(p):
        return radial_A(p, R, m, N)[n_xi] - angular_A(p, m)[n_eta]
    # bracket p in (0.05, 5 R)
    ps = np.linspace(0.05, 3.0 + 1.5 * R, 300)
    fs = np.array([f(p) for p in ps])
    roots = [brentq(f, ps[i], ps[i + 1], xtol=1e-15, rtol=1e-15) for i in range(len(ps) - 1)
             if np.sign(fs[i]) != np.sign(fs[i + 1]) and abs(fs[i]) < 50 and abs(fs[i + 1]) < 50]
    if not roots:
        raise RuntimeError("no root")
    p = roots[0]
    return -2 * p ** 2 / R ** 2, p


def radial_node(R, n_xi=1, n_eta=0, m=0, N=90):
    """xi of the radial node(s) of the exact state (n_xi radial nodes): the node surfaces are the confocal ellipsoids
    xi = const with the two protons at the foci"""
    E, p = energy(R, n_xi, n_eta, m, N)
    S = 46.0 / p
    D, x = cheb(N)
    t = (1 - x) / 2
    s = S * (t ** 2 + t) / 2
    dsdt = S * (2 * t + 1) / 2
    Dt = -2 * D
    Ds = np.diag(1 / dsdt) @ Dt
    Ds2 = np.diag(1 / dsdt ** 2) @ (Dt @ Dt) - np.diag(S / dsdt ** 3) @ Dt
    xi = 1 + s
    Lm = np.diag(xi ** 2 - 1) @ Ds2 + np.diag(2 * (m + 1) * xi) @ Ds + np.diag(m * (m + 1) + 2 * R * xi - p ** 2 * xi ** 2)
    keep = np.arange(N)
    ev, V = np.linalg.eig(Lm[np.ix_(keep, keep)])
    A = angular_A(p, m)[n_eta]
    i = np.argmin(np.abs(ev - A))
    u = V[:, i].real
    xs = xi[keep]
    o = np.argsort(xs)
    xs, u = xs[o], u[o]
    roots = []
    for j in range(len(xs) - 1):
        if np.sign(u[j]) != np.sign(u[j + 1]) and xs[j] < 30:
            # Chebyshev interpolation would be exact; local cubic root is enough (grid is dense near xi = 1)
            jj = slice(max(0, j - 2), min(len(xs), j + 3))
            c = np.polyfit(xs[jj], u[jj], min(3, len(xs[jj]) - 1))
            rr = [r.real for r in np.roots(c) if abs(r.imag) < 1e-12 and xs[j] - 1e-9 <= r.real <= xs[j + 1] + 1e-9]
            roots.append(rr[0] if rr else xs[j])
    return roots, E


STATES = {"1sg": (0, 0, 0), "1su": (0, 1, 0), "2sg": (1, 0, 0), "1pu": (0, 0, 1)}


if __name__ == "__main__":
    for N in (70, 90, 120):
        E, p = energy(2.0, N=N)
        print(f"N={N}: R=2 1sg E = {E:.13f} (ref -1.1026342144949)")
    for nm, q in STATES.items():
        print(nm, "R=2:", f"{energy(2.0, *q)[0]:.12f}")
    for R in (1.0, 1.4, 2.0, 3.0, 4.0):
        print(R, f"{energy(R)[0]:.12f}", f"total {energy(R)[0] + 1 / R:.12f}")

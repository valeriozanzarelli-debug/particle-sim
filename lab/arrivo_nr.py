"""arrivo_nr.py - REFERENCE (Schroedinger/Pauli, non-relativistic) for the Das-Duerr set-up.

This file is the *exam*, not the rule: the exact non-relativistic solution of the Das-Duerr problem
(Das & Duerr, Sci. Rep. 9, 2242 (2019), arXiv:1802.07141), used to check what the zig-zag rule gives.

Units of Das-Duerr: hbar = m = d = 1 (d = length of the initial box), time unit m d^2/hbar.
Set-up: semi-infinite waveguide z > 0, hard wall at z = 0, transverse harmonic trap omega,
initial state sqrt(2) sin(pi z) on 0 < z < 1 times the transverse ground state times a spinor
chi = (cos(alpha/2), sin(alpha/2)) (beta = 0). At t = 0 the wall at z = 1 is removed.
Psi(r,t) = phi0(x,y) exp(-i omega t) W(z,t) chi, and W(z,t) is the free evolution of the odd
extension sqrt(2) sin(pi z), |z| < 1 (method of images), written with complex error functions
(Moshinsky-type closed form).

Bohmian (Pauli current with spin term) equations of motion for beta = 0 (derived from
v = Im(grad psi/psi) + Re(grad psi/psi) x s,  s = (sin a, 0, cos a)):
    Xdot = -omega Y cos a
    Ydot = Re(W'/W)(Z,t) sin a + omega X cos a
    Zdot = Im(W'/W)(Z,t) + omega Y sin a
(alpha = 0 : spin along the guide -> Z moves as in 1-D; alpha = pi/2 : spin along x).
"""
import numpy as np
from scipy.special import erf

SQ2 = np.sqrt(2.0)


def W_and_dW(z, t):
    """W(z,t) and dW/dz for the released box state (DD units). z, t broadcastable, t > 0."""
    z = np.asarray(z, float)
    t = np.asarray(t, float)
    s2it = np.sqrt(2j * t)
    W = 0.0
    dW = 0.0
    for sgn in (+1, -1):
        p = sgn * np.pi
        z0 = z - p * t
        w1 = (1.0 - z0) / s2it
        w0 = (-1.0 - z0) / s2it
        ph = np.exp(1j * p * z - 0.5j * p * p * t)
        E = 0.5 * (erf(w1) - erf(w0))
        dE = -(np.exp(-w1 * w1) - np.exp(-w0 * w0)) / (np.sqrt(np.pi) * s2it)
        F = ph * E
        dF = 1j * p * F + ph * dE
        W = W + sgn * F
        dW = dW + sgn * dF
    # sin(pi z') = (e^{i pi z'} - e^{-i pi z'}) / (2 i), times sqrt(2)
    return SQ2 * W / 2j, SQ2 * dW / 2j


def flux(L, t):
    """quantum flux Pi_qf(t) = Im(W* W') at z = L (= spin-up Bohmian arrival density)"""
    W, dW = W_and_dW(L, t)
    return np.imag(np.conj(W) * dW)


def density(z, t):
    W, _ = W_and_dW(z, t)
    return np.abs(W)**2


def sample_initial(N, omega, rng):
    """Z from 2 sin^2(pi z) on (0,1) (rejection), X, Y from |phi0|^2 = Gaussian variance 1/(2 omega)"""
    out = []
    while sum(len(o) for o in out) < N:
        z = rng.random(4 * N)
        u = rng.random(4 * N)
        out.append(z[u < np.sin(np.pi * z)**2])
    Z = np.concatenate(out)[:N]
    s = 1.0 / np.sqrt(2.0 * omega)
    return rng.normal(0, s, N), rng.normal(0, s, N), Z


def vel(X, Y, Z, t, omega, alpha):
    W, dW = W_and_dW(Z, t)
    g = dW / W
    ca, sa = np.cos(alpha), np.sin(alpha)
    if alpha == 0.0:
        return (0.0 * X, 0.0 * Y, np.imag(g))
    return (-omega * Y * ca, np.real(g) * sa + omega * X * ca, np.imag(g) + omega * Y * sa)


def bohm_arrivals(N, L, omega, alpha, t_end, rng, dt_fn=None, t0=1e-6, XYZ0=None):
    """RK4 Bohmian trajectories; returns first-arrival times at Z = L (nan if not arrived by t_end)
    and the number of crossings recorded (first passage only). dt_fn(t) gives the step."""
    if dt_fn is None:
        def dt_fn(t):
            return 2e-4 if t < 0.05 else (2e-3 if t < 2 else 1e-2)
    if XYZ0 is None:
        X, Y, Z = sample_initial(N, omega, rng)
    else:
        X, Y, Z = (np.array(q, float) for q in XYZ0)
        N = X.shape[0]
    T = np.full(N, np.nan)
    if alpha == 0.0:
        # spin along the guide: Z decouples (the transverse rotation omega cos(alpha) does not touch Z)
        Y = np.zeros_like(Y)
        X = np.zeros_like(X)
    t = t0
    alive = np.ones(N, bool)
    while t < t_end and alive.any():
        h = dt_fn(t)
        idx = np.nonzero(alive)[0]
        x, y, zz = X[idx], Y[idx], Z[idx]
        k1 = vel(x, y, zz, t, omega, alpha)
        k2 = vel(x + 0.5 * h * k1[0], y + 0.5 * h * k1[1], zz + 0.5 * h * k1[2], t + 0.5 * h, omega, alpha)
        k3 = vel(x + 0.5 * h * k2[0], y + 0.5 * h * k2[1], zz + 0.5 * h * k2[2], t + 0.5 * h, omega, alpha)
        k4 = vel(x + h * k3[0], y + h * k3[1], zz + h * k3[2], t + h, omega, alpha)
        xn = x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        yn = y + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        zn = zz + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        zn = np.abs(zn)            # hard wall: the velocity field vanishes there; guard against overshoot
        hit = zn >= L
        if hit.any():
            # linear interpolation of the crossing time
            fr = (L - zz[hit]) / np.maximum(zn[hit] - zz[hit], 1e-300)
            T[idx[hit]] = t + h * np.clip(fr, 0, 1)
            alive[idx[hit]] = False
        X[idx], Y[idx], Z[idx] = xn, yn, zn
        t += h
    return T


if __name__ == "__main__":
    import time
    # quick checks: norm conservation and t->0 limit
    z = np.linspace(0, 400, 400001)
    for t in (0.01, 1.0, 10.0):
        print("norm(t=%g) on [0,400] = %.8f" % (t, np.trapezoid(density(z, t), z)))
    W, _ = W_and_dW(np.array([0.25, 0.5, 0.75]), 1e-7)
    print("t->0:", np.abs(W), "expected", SQ2 * np.sin(np.pi * np.array([0.25, 0.5, 0.75])))
    # quick, small Bohm test (exploratory)
    rng = np.random.default_rng(1)
    t0 = time.time()
    T = bohm_arrivals(2000, 100.0, 1000.0, np.pi / 2, 80.0, rng)
    print("spin x: arrived %.4f, max T %.3f, mean %.3f  (%.0f s)" % (np.mean(np.isfinite(T)), np.nanmax(T), np.nanmean(T), time.time() - t0))

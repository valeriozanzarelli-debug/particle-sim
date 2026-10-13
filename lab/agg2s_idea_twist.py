"""agg2s_idea_twist.py - IDEA TEST (post hoc, declared, not pre-registered): ONE complex twist theta = (pi/2, pi/2, pi/2)
instead of the pair (0, pi). For a band E(theta) = E_a + sum_R t_R cos(theta.R) the point (pi/2)^3 cancels the 6 faces,
the 12 edges and the 8 corners of the image crystal at once (cos(pi/2) = 0; edges: pairs cos(pi) + cos(0) cancel; corners:
cos(3pi/2) = cos(pi/2) = 0). In a time-dependent engine psi is complex anyway: the twist is free.
Exact Schroedinger kinetic + Galerkin field, cells a0/4, 2s, boxes L = 20, 24, 28, 32 a0 (N = 80..128).
Writes lab/out/agg2s_idea_twist.json"""
import os, sys, json, time
import numpy as np
import scipy.sparse.linalg as spl
import scipy.fft as sfft
import agg2s_core as c

JS = os.path.join(c.OUT, "agg2s_idea_twist.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}


class CAtom:
    """complex Bloch states psi(x + L e_a) = e^{i theta_a} psi(x), general theta (no parity sector)"""

    def __init__(self, N, h, th, f):
        self.N, self.h, self.f = N, h, f
        L = N * h
        k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
        ka = [k1 + t / L for t in th]
        self.T = (ka[0][:, None, None] ** 2 + ka[1][None, :, None] ** 2 + ka[2][None, None, :] ** 2) / 2
        n = np.arange(N) - N // 2
        self.ax = [np.exp(-1j * t * n / N) for t in th]
        self.ph = [np.exp(1j * k * h / 2) for k in ka]
        self.pre = 1.0 / (self.T - self.T.min() + 0.5 * c.ALPHA ** 2)

    def H(self, v):
        u = v * self.ax[0][:, None, None] * self.ax[1][None, :, None] * self.ax[2][None, None, :]
        F = sfft.fftn(u, workers=1)
        out = self.T * F
        one = np.ones(self.N)
        for sub in self.f.subs:
            p = [self.ph[i] if sub[i] else one for i in range(3)]
            tmp = sfft.ifftn(F * p[0][:, None, None] * p[1][None, :, None] * p[2][None, None, :], workers=1)
            tmp *= self.f.Vsub(sub)
            out += sfft.fftn(tmp, workers=1) * (np.conj(p[0])[:, None, None] * np.conj(p[1])[None, :, None] * np.conj(p[2])[None, None, :]) / 8
        o = sfft.ifftn(out, workers=1)
        return o * np.conj(self.ax[0])[:, None, None] * np.conj(self.ax[1])[None, :, None] * np.conj(self.ax[2])[None, None, :]

    def Pc(self, v):
        u = v * self.ax[0][:, None, None] * self.ax[1][None, :, None] * self.ax[2][None, None, :]
        o = sfft.ifftn(self.pre * sfft.fftn(u, workers=1), workers=1)
        return o * np.conj(self.ax[0])[:, None, None] * np.conj(self.ax[1])[None, :, None] * np.conj(self.ax[2])[None, None, :]


def mat(fn, N):
    def m(X):
        X = np.asarray(X).reshape(N ** 3, -1)
        return np.stack([fn(X[:, j].reshape((N,) * 3)).ravel() for j in range(X.shape[1])], 1)
    return m


h = c.A0 / 4
for L in (20, 24, 28, 32):
    N = 4 * L
    key = f"L{L}"
    if key in res:
        continue
    f = c.galerkin_field(N, h, cache=False)
    at = c.CAtom(N, h, (np.pi / 2,) * 3, f) if hasattr(c, "CAtom") else CAtom(N, h, (np.pi / 2,) * 3, f)
    x = f.x
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij", sparse=True)
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2) / c.A0
    A = spl.LinearOperator((N ** 3,) * 2, matvec=mat(at.H, N), matmat=mat(at.H, N), dtype=complex)
    M = spl.LinearOperator((N ** 3,) * 2, matvec=mat(at.Pc, N), matmat=mat(at.Pc, N), dtype=complex)
    g1 = np.broadcast_to(np.exp(-r), (N,) * 3).astype(complex).ravel()
    w, v1 = spl.lobpcg(A, g1[:, None] / np.linalg.norm(g1), M=M, tol=1e-9, maxiter=60, largest=False)
    g2 = np.broadcast_to((1 - r / 2) * np.exp(-r / 2), (N,) * 3).astype(complex).ravel()
    w2, v2 = spl.lobpcg(A, g2[:, None] / np.linalg.norm(g2), M=M, Y=v1, tol=1e-9, maxiter=80, largest=False)
    E2 = float(w2[0].real) * c.MC2
    res[key] = dict(L_a0=L, E2s_rel=E2 / (-c.RY / 4) - 1, E1s_rel=float(w[0].real) * c.MC2 / (-c.RY) - 1)
    json.dump(res, open(JS, "w"), indent=1)
    print(key, res[key], flush=True)

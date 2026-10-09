"""kk_core.py - the zig-zag rule with ONE hidden compact direction (Kaluza-Klein on the lattice).

Built on lab/zz_core.py (same rule, same conventions): every ball moves at c along each axis in turn
(split on the eigenvectors of sigma_j, L moves opposite to R) and then rotates L<->R by the flip
angle theta = m c a / hbar. Here the three walk axes are
    axis 0 = x (real), axis 1 = y (real), axis 2 = w (HIDDEN circle of N sites, periodic).
So the ball also hops around the hidden circle with the same translate+rotate rule, nothing new.
Lattice units a = dt = c = hbar = 1. Hidden radius R = N a / (2 pi). The "charge" of a ball is its
momentum around the hidden circle, k_w = n / R with n an integer (the circle closes on itself).

Twisted geometry (the only thing added, and it is geometry, not a force): when a ball hops along a
real axis j it also slides around the hidden circle by a distance D_j (an off-diagonal metric
g_{jw}, i.e. a background gauge field A_j = D_j / R written as geometry). D_j is stored on the link
that ARRIVES at the site. The slide by a non-integer distance is the exact band-limited translation
on the circle (multiply hidden Fourier mode n by exp(-i n D / R)).
psi shape: (4, nx, ny, N) = (L0, L1, R0, R1). For 1+1 D use ny = 1.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from zz_core import BASIS_HALF, order_of_tick, bloch_batch

SEED = 20261009


def nvals(N):
    """integer hidden momenta of the N lattice modes (FFT order)"""
    return np.rint(np.fft.fftfreq(N) * N).astype(int)


class KKWalk:
    def __init__(self, nx, ny, N, theta, V=None, D=None, Dsym=None):
        self.shape = (nx, ny, N)
        self.N, self.theta = N, theta
        self.R = N / (2 * np.pi)
        self.kw = 2 * np.pi * nvals(N) / N                     # hidden wavenumbers n / R
        self.phase = None if V is None else np.exp(-1j * V)    # V shape broadcastable to (nx,ny,N)
        self.D = D or {}                                       # {axis: array (nx, ny)} link slide
        self.t = 0

    def _hshift(self, c, d):
        """slide the field c (nx,ny,N) around the hidden circle by d (nx,ny)"""
        F = np.fft.fft(c, axis=2)
        F *= np.exp(-1j * self.kw[None, None, :] * d[..., None])
        return np.fft.ifft(F, axis=2)

    def _axis(self, blk, j, h):
        lam, P = BASIS_HALF[j]
        c = np.tensordot(P.conj().T, blk, axes=(1, 0))
        for q in range(2):
            s = h * lam[q]
            c[q] = np.roll(c[q], s, axis=j)
            if j in self.D:
                Dj = self.D[j]
                d = Dj if s > 0 else -np.roll(Dj, -1, axis=j)
                c[q] = self._hshift(c[q], d)
        return np.tensordot(P, c, axes=(1, 0))

    def tick(self, psi):
        L, R = psi[0:2], psi[2:4]
        for j in order_of_tick(self.t):
            if self.shape[j] == 1 and j not in self.D:
                L = self._axis(L, j, -1); R = self._axis(R, j, +1)   # trivial roll, keeps the phase k=0
                continue
            L = self._axis(L, j, -1)
            R = self._axis(R, j, +1)
        ct, st = np.cos(self.theta), np.sin(self.theta)
        out = np.concatenate([ct * L - 1j * st * R, ct * R - 1j * st * L])
        if self.phase is not None:
            out *= self.phase
        self.t += 1
        return out


def charge(psi):
    """<n> = mean hidden momentum (the 'charge') and the weight on each integer n"""
    N = psi.shape[-1]
    F = np.fft.fft(psi, axis=-1)
    w = np.sum(np.abs(F)**2, axis=tuple(range(psi.ndim - 1)))
    w = w / w.sum()
    n = nvals(N)
    return float(np.sum(n * w)), dict(zip(n.tolist(), w.tolist()))


def packet(nx, ny, N, theta, n, k0=(0.0, 0.0), sigma=5.0, center=None, spinor=(1, 0, 1, 0)):
    """positive-energy wave packet (lattice-exact projector of the two-tick Bloch matrix) with hidden
    momentum n: Gaussian in x,y times exp(i n w / R) on the hidden circle"""
    kx = 2 * np.pi * np.fft.fftfreq(nx)
    ky = 2 * np.pi * np.fft.fftfreq(ny)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    kw = 2 * np.pi * n / N
    cx, cy = center if center is not None else (nx // 2, ny // 2)
    # Gaussian amplitude in k-space centred on k0, phase puts it at (cx, cy)
    dkx = (KX - k0[0] + np.pi) % (2 * np.pi) - np.pi
    dky = (KY - k0[1] + np.pi) % (2 * np.pi) - np.pi
    g = np.exp(-sigma**2 * (dkx**2 + (dky**2 if ny > 1 else 0 * dky)) - 1j * (KX * cx + KY * cy))
    KA = np.stack([KX.ravel(), KY.ravel(), np.full(KX.size, kw)], axis=1)
    W = bloch_batch(KA, theta)
    lam, V = np.linalg.eig(W)
    om = -np.angle(lam) / 2
    s0 = np.asarray(spinor, complex) / np.linalg.norm(spinor)
    coef = np.empty((KA.shape[0], 4), complex)
    for i in range(KA.shape[0]):
        Vp = V[i][:, om[i] > 0]
        Pp = Vp @ np.linalg.solve(Vp.conj().T @ Vp, Vp.conj().T)
        coef[i] = Pp @ s0
    coef = coef.reshape(nx, ny, 4).transpose(2, 0, 1) * g[None]
    psi2 = np.fft.ifft2(coef, axes=(1, 2))
    w = np.arange(N)
    psi = psi2[..., None] * np.exp(1j * kw * w)[None, None, None, :]
    return psi / np.sqrt(np.sum(np.abs(psi)**2))


def density(psi):
    return np.sum(np.abs(psi)**2, axis=(0, 3))


def ring_ab(L, N, theta, n, Dx, p=0.5, sigma=10.0, ticks=None, window=30, init_phase=None):
    """Aharonov-Bohm ring (1+1 D, ring of L sites, ny = 1): two packets with +p and -p start together
    at x = 0, run around the two halves of the ring and meet at x = L/2. Dx (L, 1): hidden slide on the
    x-links. Linear rule -> evolve the two halves separately (psi = a + b) and read the interference term
    2 Re(a^dag b) directly; returns (density, phase of sum a^dag b in the window, fringe phase from the density)"""
    a = packet(L, 1, N, theta, n, k0=(p, 0), sigma=sigma, center=(0, 0))
    b = packet(L, 1, N, theta, n, k0=(-p, 0), sigma=sigma, center=(0, 0))
    if init_phase is not None:
        a = a * init_phase[None, :, None, None]; b = b * init_phase[None, :, None, None]
    W1 = KKWalk(L, 1, N, theta, D={0: Dx}); W2 = KKWalk(L, 1, N, theta, D={0: Dx})
    for _ in range(ticks):
        a = W1.tick(a); b = W2.tick(b)
    psi = a + b
    rho = density(psi)[:, 0]
    xa = np.arange(L); win = np.abs(xa - L // 2) < window
    cross = np.sum(np.conj(a) * b, axis=(0, 3))[:, 0]
    ph_overlap = np.angle(np.sum(cross[win]))
    ph_fringe = np.angle(np.sum(rho[win] * np.exp(-2j * p * xa[win])))
    return rho, ph_overlap, ph_fringe

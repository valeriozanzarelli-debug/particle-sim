"""zz_fast.py - numba versions of the zig-zag rule (same rule as zz_core.py, only faster).
Each kernel is checked against the numpy reference in the self-test (python lab/zz_fast.py)."""
import numpy as np
import numba as nb
from zz_core import BASIS_HALF, order_of_tick


def _bases():
    up = np.zeros((3, 2), np.complex128)
    dn = np.zeros((3, 2), np.complex128)
    for j in range(3):
        lam, P = BASIS_HALF[j]
        up[j], dn[j] = P[:, 0], P[:, 1]        # sigma_j = +1, -1
    return up, dn


UP, DN = _bases()


@nb.njit(parallel=True, cache=True)
def _axis_block(psi, off, j, h, up, dn):
    """translate one Weyl block psi[off:off+2] along axis j (periodic): the sigma_j = +1 part moves h cells,
    the -1 part -h cells (h = +1 right-handed, -1 left-handed)"""
    n = psi.shape[1]
    cp = np.empty((n, n, n), np.complex128)
    cm = np.empty((n, n, n), np.complex128)
    u0, u1, d0, d1 = up[j, 0], up[j, 1], dn[j, 0], dn[j, 1]
    for i in nb.prange(n):
        for k in range(n):
            for l in range(n):
                a = psi[off, i, k, l]
                b = psi[off + 1, i, k, l]
                cp[i, k, l] = np.conj(u0) * a + np.conj(u1) * b
                cm[i, k, l] = np.conj(d0) * a + np.conj(d1) * b
    for i in nb.prange(n):
        for k in range(n):
            for l in range(n):
                if j == 0:
                    p = cp[(i - h) % n, k, l]; m = cm[(i + h) % n, k, l]
                elif j == 1:
                    p = cp[i, (k - h) % n, l]; m = cm[i, (k + h) % n, l]
                else:
                    p = cp[i, k, (l - h) % n]; m = cm[i, k, (l + h) % n]
                psi[off, i, k, l] = u0 * p + d0 * m
                psi[off + 1, i, k, l] = u1 * p + d1 * m


@nb.njit(parallel=True, cache=True)
def _coin(psi, ct, st, phase):
    n = psi.shape[1]
    for i in nb.prange(n):
        for k in range(n):
            for l in range(n):
                ph = phase[i, k, l]
                for s in range(2):
                    L = psi[s, i, k, l]
                    R = psi[2 + s, i, k, l]
                    psi[s, i, k, l] = ph * (ct * L - 1j * st * R)
                    psi[2 + s, i, k, l] = ph * (ct * R - 1j * st * L)


class FastDirac3D:
    """in-place 3-D Dirac walk (periodic box), identical to zz_core.Walk('dirac', ...)"""

    def __init__(self, n, a, theta, V=None):
        self.n, self.a, self.theta = n, a, theta
        self.phase = np.exp(-1j * V * a) if V is not None else np.ones((n, n, n), complex)
        self.tick_count = 0

    def tick(self, psi):
        for j in order_of_tick(self.tick_count):
            _axis_block(psi, 0, j, -1, UP, DN)
            _axis_block(psi, 2, j, +1, UP, DN)
        _coin(psi, np.cos(self.theta), np.sin(self.theta), self.phase)
        self.tick_count += 1
        return psi


@nb.njit(cache=True)
def radial_listen(M, refl0, refl1, chi0, nt, nus, centers, widths, theta, record=True, damp=None):
    """run the radial walk nt ticks from chi0. Records C(n) = <chi0|chi_n> e^{i theta n} (demodulated by
    the rest frequency) and lock-in sums acc_k = sum_n w_k(n) e^{i (theta - nu_k) n} chi_n with Gaussian
    windows w_k (centre, width; cut at 3 widths). nu_k = binding frequency per tick to listen at."""
    N = M.shape[0]
    p = chi0[0].copy()
    q = chi0[1].copy()
    pn = np.empty(N, np.complex128)
    qn = np.empty(N, np.complex128)
    C = np.empty(nt, np.complex128)
    K = nus.shape[0]
    acc = np.zeros((K, 2, N), np.complex128)
    for n in range(nt):
        if record:
            s = 0j
            for i in range(N):
                s += np.conj(chi0[0, i]) * p[i] + np.conj(chi0[1, i]) * q[i]
            C[n] = s * np.exp(1j * theta * n)
        for k in range(K):
            d = n - centers[k]
            if abs(d) < 3.0 * widths[k]:
                w = np.exp(-0.5 * (d / widths[k])**2) * np.exp(1j * (theta - nus[k]) * n)
                for i in range(N):
                    acc[k, 0, i] += w * p[i]
                    acc[k, 1, i] += w * q[i]
        # translate: outward walker one cell out, inward walker one cell in, reflection at the ends
        pn[0] = refl0 * q[0]
        for i in range(1, N):
            pn[i] = p[i - 1]
        for i in range(N - 1):
            qn[i] = q[i + 1]
        qn[N - 1] = refl1 * p[N - 1]
        # rotate (and, if given, let the balls fade out in a sponge near the outer edge)
        for i in range(N):
            a = pn[i]
            b = qn[i]
            p[i] = M[i, 0, 0] * a + M[i, 0, 1] * b
            q[i] = M[i, 1, 0] * a + M[i, 1, 1] * b
            if damp is not None:
                p[i] *= damp[i]
                q[i] *= damp[i]
    return C, acc, np.stack((p, q))


if __name__ == "__main__":
    import time
    from zz_core import Walk, RadialWalk, SEED
    rng = np.random.default_rng(SEED)
    n = 12
    V = rng.random((n, n, n))
    psi = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n))
    Wr = Walk("dirac", n, 0.3, theta=0.3, V=V)
    Wf = FastDirac3D(n, 0.3, 0.3, V)
    a, b = psi.copy(), psi.copy()
    for _ in range(5):
        a = Wr.tick(a)
        b = Wf.tick(b)
    print("numba 3-D Dirac tick vs numpy walk: max diff", np.max(np.abs(a - b)))
    rw = RadialWalk(-1, 0.5, 200.0)
    chi = rng.normal(size=(2, rw.N)) + 1j * rng.normal(size=(2, rw.N))
    ref = chi.copy()
    for _ in range(7):
        ref = rw.tick(ref)
    C, acc, last = radial_listen(rw.M, rw.refl0, rw.refl1, chi, 7, np.zeros(0), np.zeros(0), np.zeros(0), 0.5)
    print("numba radial listen vs numpy walk: max diff", np.max(np.abs(last - ref)))
    t = time.time()
    n = 100
    W = FastDirac3D(n, 0.24, 0.24)
    p = np.ones((4, n, n, n), complex)
    W.tick(p)
    t = time.time()
    for _ in range(3):
        W.tick(p)
    print("100^3 tick: %.3f s" % ((time.time() - t) / 3))

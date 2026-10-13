"""vuoto_dev_check.py - DEVELOPMENT check (declared, not graded): the k-space linear-response code lr_terms against
the real-space sea (exact diagonalisation of the walk on a ring / small 3-D box, finite differences in V).
Also: the 1-D 'theorem' S U S = -U  ->  the whole one-tick sea density is uniform for ANY potential.
Run: python lab/vuoto_dev_check.py   (~1-3 min)"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from vuoto_core import ueig, lr_terms, species_weight, U1d, U3d, PI

t0 = time.time()


def walk1d_matrix(N, th, V):
    """one tick of zz_core 'dirac1' on a ring: basis index c*N + x, c = 0 (R), 1 (L)"""
    S = np.zeros((2 * N, 2 * N))
    x = np.arange(N)
    S[(x + 1) % N, x] = 1                # R(x+1) <- R(x)
    S[N + (x - 1) % N, N + x] = 1        # L(x-1) <- L(x)
    c, s = np.cos(th), np.sin(th)
    C = np.block([[c * np.eye(N), -1j * s * np.eye(N)], [-1j * s * np.eye(N), c * np.eye(N)]])
    Ph = np.diag(np.exp(-1j * np.concatenate([V, V])))
    return Ph @ C @ S


def sea_density(W, ncomp, N, Q=None):
    Om, lam, V = ueig(W[None])
    Vo = V[0][:, Om[0] < 0]
    P = Vo @ Vo.conj().T
    if Q is not None:
        P = Q @ P @ Q
    return np.real(np.diag(P)).reshape(ncomp, -1).sum(0)


def fourier_Q(N, w, ncomp, dim=1):
    k = 2 * PI * np.fft.fftfreq(N)
    if dim == 1:
        F = np.exp(-1j * np.outer(k, np.arange(N))) / np.sqrt(N)
        q1 = F.conj().T @ np.diag(np.sqrt(w(k))) @ F
    else:
        KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
        Ks = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], 1)
        F1 = np.exp(-1j * np.outer(k, np.arange(N))) / np.sqrt(N)
        F = np.kron(np.kron(F1, F1), F1)
        q1 = F.conj().T @ np.diag(np.sqrt(w(Ks))) @ F
    return np.kron(np.eye(ncomp), q1)


# ---------------- 1-D
N, th, j, v = 64, 0.3, 3, 1e-4
q = 2 * PI * j / N
x = np.arange(N)
kk = 2 * PI * np.arange(N) / N - PI
w = species_weight()
for period in (1, 2):
    for lab, wt in (("whole sea", None), ("physical species Q", w)):
        Q = None if wt is None else fourier_Q(N, wt, 2)
        dn = []
        for sgn in (+1, -1):
            U = walk1d_matrix(N, th, sgn * v * np.cos(q * x))
            W = U if period == 1 else U @ U
            dn.append(sea_density(W, 2, N, Q))
        d = (dn[0] - dn[1]) / 2
        dq = np.sum(np.exp(-1j * q * x) * d) / N
        chi_rs = -dq / (v / 2)
        chi_k = np.mean(lr_terms(kk, q, 1, th, period, wt))
        print(f"1-D period {period}, {lab:20s}: real space chi = {chi_rs.real:+.10e}{chi_rs.imag:+.1e}i   k-space = {chi_k.real:+.10e}{chi_k.imag:+.1e}i")
# theorem: whole one-tick sea uniform for a STRONG potential
Vs = -0.8 / np.sqrt((x - N / 2) ** 2 + 1.0)
n_sea = sea_density(walk1d_matrix(N, th, Vs), 2, N)
print("1-D one tick, strong potential (V(0) = -0.8 per tick): whole-sea density spread", n_sea.max() - n_sea.min(), " mean", n_sea.mean())
n_sea2 = sea_density(walk1d_matrix(N, th, Vs) @ walk1d_matrix(N, th, Vs), 2, N)
print("1-D two-tick window sea, same potential: density spread", n_sea2.max() - n_sea2.min())

# ---------------- 3-D small box
n3, th3, v3 = 8, 0.4, 1e-4
q3 = np.array([2 * PI / n3, 0, 0])


def walk3d_matrix(n, th, V, nt):
    """apply zz_core.Walk('dirac') tick to basis vectors (dense matrix), V on the box"""
    from zz_core import Walk
    Wk = Walk("dirac", n, 1.0, theta=th, V=V)
    Wk.tick_count = nt
    D = 4 * n**3
    M = np.zeros((D, D), complex)
    E = np.eye(D, dtype=complex)
    for i in range(D):
        Wk.tick_count = nt
        M[:, i] = Wk.tick(E[:, i].reshape(4, n, n, n)).ravel()
    return M


X = np.arange(n3)
Vfield = np.cos(q3[0] * X)[:, None, None] * np.ones((n3, n3, n3))
k1 = 2 * PI * np.fft.fftfreq(n3)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
Ks = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], 1)
for lab, wt in (("whole sea", None), ("physical species Q", w)):
    Q = None if wt is None else fourier_Q(n3, wt, 4, dim=3)
    dn = []
    for sgn in (+1, -1):
        Vv = sgn * v3 * Vfield
        W = walk3d_matrix(n3, th3, Vv, 1) @ walk3d_matrix(n3, th3, Vv, 0)
        dn.append(sea_density(W, 4, n3, Q).reshape(n3, n3, n3))
    d = (dn[0] - dn[1]) / 2
    dq = np.sum(np.exp(-1j * q3[0] * X)[:, None, None] * d) / n3**3
    chi_rs = -dq / (v3 / 2)
    chi_k = np.mean(lr_terms(Ks, q3, 3, th3, 2, wt))
    print(f"3-D {n3}^3, {lab:20s}: real space chi = {chi_rs.real:+.10e}{chi_rs.imag:+.1e}i   k-space = {chi_k.real:+.10e}{chi_k.imag:+.1e}i")
print("time %.1f s" % (time.time() - t0))

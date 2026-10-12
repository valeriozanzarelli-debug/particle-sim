"""agg_zz_core.py - the AGGLOMERATOR applied to the zig-zag rule (Valerio: "real hydrogen needs ~1e10 cells:
where did the agglomerator go?").

One automatic rule, the same at every level (docs/ZOOM.md: block -> learn -> step -> check -> repeat):
  * BLOCK   b cells per axis -> 1 super-cell (block mean x b^{d/2}, so the norm is kept), b steps of the level
            below -> 1 step (the coarse tick tau = B in fine units: the super-cell still "moves at c").
  * LEARN   a translation-invariant linear one-step kernel  psi'(X) = sum_j K_j psi(X + j),  K_j = C x C complex
            matrices, |j|_inf <= w, by least squares on (blocked state, blocked state one coarse step later) pairs
            of random band-limited packets evolved by the level below (normal equations built with FFTs).
  * STEP    psi -> FFT -> K(k) psi(k) -> iFFT.
  * CHECK   Gaussian packets: candidate coarse rule vs the level below over H coarse steps, error |dpsi|^2/|psi|^2.
  * CHOOSE  (automatic, Nanite-like): among b in {2,3,4} (only multiples of the lower level's period) and
            w in {1,..,4}, keep the candidates whose check error <= eps_level; take the largest b, then the
            smallest w ("the rule that coarse-grains most for the error it is allowed to lose").
  * REPEAT  the next level is learned from data made by THIS learned rule (not from the fine rule): true repetition.

Fine rules (lab/zz_core.py, frozen): 1-D Dirac walk (R moves +1, L moves -1, flip theta) and the 3-D cubic walk
(4 components, palindromic x,y,z / z,y,x order: period 2 ticks). theta = m a (m = 1 in units hbar = c = m_e = 1).
Descriptors: 1-D: (R, L) block means (both energy bands); 3-D: the large Dirac component phi = (L + R)/sqrt 2 of
positive-energy states (2 components = the Pauli spinor; the coarse rule is learned on phi only).
"""
import itertools, time
import numpy as np
import scipy.fft as sfft
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zz_core as zz

WK = 2                     # FFT workers (<= 2 threads: shared machine)


def fftn(a, d):
    return sfft.fftn(a, axes=tuple(range(-d, 0)), workers=WK)


def ifftn(a, d):
    return sfft.ifftn(a, axes=tuple(range(-d, 0)), workers=WK)


def kgrid(n, d):
    k1 = 2 * np.pi * np.fft.fftfreq(n)
    return np.meshgrid(*([k1] * d), indexing="ij")


def offsets(w, d):
    return np.array(list(itertools.product(range(-w, w + 1), repeat=d)), int)


# ============================== fine rules ==============================
def fine1_tick(psi, th):
    """1-D zig-zag (zz_core Walk 'dirac1'); psi (..., 2, n) = (R, L)"""
    R = np.roll(psi[..., 0, :], 1, -1)
    L = np.roll(psi[..., 1, :], -1, -1)
    c, s = np.cos(th), np.sin(th)
    return np.stack([c * R - 1j * s * L, c * L - 1j * s * R], -2)


def fine1_bloch(k, th):
    """one-tick Bloch matrix of the 1-D walk, k = k a (array) -> (..., 2, 2)"""
    c, s = np.cos(th), np.sin(th)
    T = np.zeros(k.shape + (2, 2), complex)
    T[..., 0, 0] = np.exp(-1j * k)
    T[..., 1, 1] = np.exp(1j * k)
    C = np.array([[c, -1j * s], [-1j * s, c]])
    return C @ T


def _shift_weyl(sp, j, h, ax):
    """sp (..., 2, n, n, n) Weyl spinor; move eigen-parts of sigma_j by h*lambda along spatial axis ax"""
    lam, P = zz.BASIS_HALF[j]
    c = np.einsum("ba,...bxyz->...axyz", P.conj(), sp)
    for q in range(2):
        c[..., q, :, :, :] = np.roll(c[..., q, :, :, :], h * lam[q], axis=ax)
    return np.einsum("ab,...bxyz->...axyz", P, c)


def fine3_tick(psi, th, n_tick):
    """3-D zig-zag tick (same as zz_core.Walk('dirac') without potential); psi (S, 4, n, n, n)"""
    L, R = psi[:, 0:2], psi[:, 2:4]
    for j in zz.order_of_tick(n_tick):
        ax = 2 + j - psi.ndim          # negative spatial axis index of the (S,2,n,n,n) block
        L = _shift_weyl(L, j, -1, ax)
        R = _shift_weyl(R, j, +1, ax)
    c, s = np.cos(th), np.sin(th)
    return np.concatenate([c * L - 1j * s * R, c * R - 1j * s * L], axis=1)


# ============================== block / random data ==============================
def block(psi, b, d):
    sh = psi.shape[:-d]
    n = psi.shape[-1]
    m = n // b
    new = sh + sum(((m, b) for _ in range(d)), ())
    x = psi.reshape(new)
    ax = tuple(len(sh) + 2 * i + 1 for i in range(d))
    return x.mean(axis=ax) * b ** (d / 2)


def band_field(rng, S, C, n, d, kmax):
    """random complex field, flat spectrum for |k| < kmax (k in units of this grid's 1/cell), zero outside"""
    K = kgrid(n, d)
    k2 = sum(k * k for k in K)
    mask = k2 < kmax ** 2
    f = (rng.normal(size=(S, C) + (n,) * d) + 1j * rng.normal(size=(S, C) + (n,) * d)) * mask
    x = ifftn(f, d)
    return x / np.sqrt(np.sum(np.abs(x) ** 2, axis=tuple(range(1, x.ndim)), keepdims=True))


def gauss_packets(rng, S, C, n, d, sigma, kappa_max, spinors=None):
    """Gaussian packets centred in the box, width sigma cells, random mean momentum |k_i| < kappa_max, random spinor"""
    x = np.arange(n) - n / 2 + 0.5
    X = np.meshgrid(*([x] * d), indexing="ij")
    out = np.zeros((S, C) + (n,) * d, complex)
    for s in range(S):
        k0 = rng.uniform(-kappa_max, kappa_max, d)
        g = np.exp(-sum(Xi * Xi for Xi in X) / (4 * sigma ** 2) + 1j * sum(k0[i] * X[i] for i in range(d)))
        sp = (rng.normal(size=C) + 1j * rng.normal(size=C)) if spinors is None else spinors[s]
        sp = sp / np.linalg.norm(sp)
        out[s] = sp.reshape((C,) + (1,) * d) * g
    return out / np.sqrt(np.sum(np.abs(out) ** 2, axis=tuple(range(1, out.ndim)), keepdims=True))


# ============================== learned rule ==============================
class Rule:
    """translation-invariant one-step kernel on a periodic lattice"""

    def __init__(self, K, offs, d, C, unitary=True):
        self.K, self.offs, self.d, self.C = K, offs, d, C     # K (C, C, ntap): psi'_a(X) = sum K[a,b,j] psi_b(X+off_j)
        self.unitary = unitary   # "the balls are conserved": K(k) -> its unitary polar factor at every k
        self._cache = {}

    @staticmethod
    def polar(M):
        """M (..., C, C) -> unitary polar factor U = W V^dag (M = W S V^dag)"""
        W, s, Vh = np.linalg.svd(M)
        return W @ Vh

    def mult(self, n):
        """Fourier multiplier K(k) on an n^d grid: sum_j K_j e^{i k.j} (stencil placed on the grid, inverse FFT)"""
        if n not in self._cache:
            A = np.zeros((self.C, self.C) + (n,) * self.d, complex)
            idx = tuple(self.offs[:, i] % n for i in range(self.d))
            for a in range(self.C):
                for b in range(self.C):
                    A[(a, b) + idx] = self.K[a, b]
            M = ifftn(A, self.d) * n ** self.d
            if self.unitary:
                M = np.moveaxis(self.polar(np.moveaxis(M, (0, 1), (-2, -1))), (-2, -1), (0, 1))
            self._cache[n] = M
        return self._cache[n]

    def mult_at(self, kvec):
        """K(k) at arbitrary k points: kvec (..., d) in units of 1/cell"""
        ph = np.exp(1j * np.tensordot(kvec, self.offs.T, axes=(-1, 0)))                               # (..., ntap)
        M = np.einsum("abj,...j->...ab", self.K, ph)
        return self.polar(M) if self.unitary else M

    def step(self, psi, nsteps=1):
        d = self.d
        n = psi.shape[-1]
        M = self.mult(n)
        f = fftn(psi, d)
        for _ in range(nsteps):
            f = np.einsum("ab...,sb...->sa...", M, f)
        return ifftn(f, d)


class GenRule(Rule):
    """rule stored as a learned GENERATOR (coarse Hamiltonian stencil): K(k) = e^{i phi0} exp(-i G(k)),
    G(k) = sum_j G_j e^{i k.j}, G_{-j} = G_j^dag (Hermitian: the balls are conserved exactly)"""

    def __init__(self, G, offs, d, C, phi0):
        super().__init__(G, offs, d, C, unitary=False)
        self.phi0 = phi0

    @staticmethod
    def expm_herm(Gk):
        Gk = 0.5 * (Gk + np.conj(np.swapaxes(Gk, -1, -2)))
        e, U = np.linalg.eigh(Gk)
        return (U * np.exp(-1j * e)[..., None, :]) @ np.conj(np.swapaxes(U, -1, -2))

    def gen_at(self, kvec):
        ph = np.exp(1j * np.tensordot(kvec, self.offs.T, axes=(-1, 0)))
        return np.einsum("abj,...j->...ab", self.K, ph)

    def mult_at(self, kvec):
        return np.exp(1j * self.phi0) * self.expm_herm(self.gen_at(kvec))

    def gen_grid(self, n):
        """the generator G(k) on an n^d grid, shape (n.., C, C)"""
        A = np.zeros((self.C, self.C) + (n,) * self.d, complex)
        idx = tuple(self.offs[:, i] % n for i in range(self.d))
        for a in range(self.C):
            for b in range(self.C):
                A[(a, b) + idx] = self.K[a, b]
        return np.moveaxis(ifftn(A, self.d) * n ** self.d, (0, 1), (-2, -1))

    def mult(self, n):
        if n not in self._cache:
            A = np.zeros((self.C, self.C) + (n,) * self.d, complex)
            idx = tuple(self.offs[:, i] % n for i in range(self.d))
            for a in range(self.C):
                for b in range(self.C):
                    A[(a, b) + idx] = self.K[a, b]
            Gk = np.moveaxis(ifftn(A, self.d) * n ** self.d, (0, 1), (-2, -1))
            M = np.exp(1j * self.phi0) * self.expm_herm(Gk)
            self._cache[n] = np.moveaxis(M, (-2, -1), (0, 1))
        return self._cache[n]


def learn_h(prev, b, desc, w, ks, rng, nk=None, ridge=1e-14, frac_uniform=0.25, w_uniform=0.03, relative=True):
    """least squares for the GENERATOR on plane-wave data: T(k) = (level below)^b seen through the descriptor,
    phi0 = its phase at k = 0 (the rest phase, mod 2 pi), G(k) = i log(T e^{-i phi0}); fit G_j (|j| <= w)."""
    d = prev.d
    offs = offsets(w, d)
    nt = len(offs)
    nk = nk if nk is not None else max(4000, 6 * nt)
    ng = int(nk * (1 - frac_uniform))
    if ks < 0:                          # negative ks = sample like the observable's state of size 1/|ks| cells
        kg = sample_state_k(rng, ng, d, -ks)
        kg = kg[np.all(np.abs(kg) < np.pi, axis=1)]
        ng = len(kg)
    else:
        kg = np.clip(rng.normal(size=(ng, d)) * ks, -np.pi, np.pi)
    ku = np.pi * (2 * rng.random((int(nk * frac_uniform), d)) - 1)
    kap = np.concatenate([kg, ku, np.zeros((1, d))])
    T = conv_k(prev.mult_k(kap / b, b), prev.desc, desc)
    lam0 = np.linalg.eigvals(T[-1])
    phi0 = float(np.angle(np.mean(lam0)))
    lam, V = np.linalg.eig(T * np.exp(-1j * phi0))
    lg = np.real(1j * np.log(lam / np.abs(lam)))       # motion phase per step, principal branch (-pi, pi]
    lg = np.where(lg < -0.5, lg + 2 * np.pi, lg)        # unwrap: the positive band never moves backwards in energy
    Gs = V @ (lg[..., :, None] * np.linalg.inv(V))                                   # (nk, C, C)
    sw = np.ones(len(kap)); sw[ng:-1] = w_uniform
    Phi = np.exp(1j * kap @ offs.T)
    if relative:                        # relative least squares: every sample's error counted against its own size
        g = np.linalg.norm(Gs, axis=(-2, -1))
        sw = sw / (g + 1e-3 * np.median(g[:ng]))
    Phi = Phi * sw[:, None]
    C = Gs.shape[-1]
    rhs = (Gs * sw[:, None, None]).reshape(len(kap), C * C)
    sol = np.linalg.lstsq(Phi, rhs, rcond=1e-13)[0]                            # QR/SVD: no squared condition number
    G = sol.T.reshape(C, C, nt)
    # Hermitian symmetrisation G_{-j} = G_j^dag
    jidx = {tuple(o): i for i, o in enumerate(offs)}
    neg = np.array([jidx[tuple(-o)] for o in offs])
    G = 0.5 * (G + np.conj(np.swapaxes(G, 0, 1))[:, :, neg])
    return GenRule(G, offs, d, C, phi0)


def learn(X, Y, w, d, ridge=1e-11, unitary=True):
    """least-squares kernel: Y_a(x) = sum_{b,j} K[a,b,j] X_b(x + off_j); X, Y (S, C, n..d)"""
    S, C = X.shape[:2]
    offs = offsets(w, d)
    nt = len(offs)
    Fx = fftn(X, d)
    Fy = fftn(Y, d)
    # A[b,b'](t) = sum_s sum_x conj X_b(x) X_b'(x + t);  B[b,a](t) = sum_s sum_x conj X_b(x) Y_a(x + t)
    A = ifftn(np.einsum("sb...,sc...->bc...", Fx.conj(), Fx), d)
    B = ifftn(np.einsum("sb...,sa...->ba...", Fx.conj(), Fy), d)
    n = X.shape[-1]

    def at(arr, t):          # arr (.., n..d) at integer lag vectors t (m, d)
        idx = tuple((t[:, i] % n) for i in range(d))
        return arr[(Ellipsis,) + idx]
    # G[(b,j),(c,j')] = sum conj X_b(x + j) X_c(x + j') = A[b,c](j' - j)
    lag = offs[None, :, :] - offs[:, None, :]                                   # (nt, nt, d): j' - j
    G = np.zeros((C, nt, C, nt), complex)
    Lf = lag.reshape(-1, d)
    for b in range(C):
        for c in range(C):
            G[b, :, c, :] = at(A[b, c], Lf).reshape(nt, nt)
    G = G.reshape(C * nt, C * nt)
    G = 0.5 * (G + G.conj().T)
    lam = ridge * np.real(np.trace(G)) / G.shape[0]
    Gr = G + lam * np.eye(G.shape[0])
    K = np.zeros((C, C, nt), complex)
    for a in range(C):
        r = np.zeros((C, nt), complex)
        for b in range(C):
            r[b] = at(B[b, a], -offs)                                            # sum conj X_b(x+j) Y_a(x) = B[b,a](-j)
        sol = np.linalg.solve(Gr, r.reshape(-1))
        K[a] = sol.reshape(C, nt)
    return Rule(K, offs, d, C, unitary)


def rel_err(a, b):
    ax = tuple(range(1, a.ndim))
    return np.sum(np.abs(a - b) ** 2, axis=ax) / np.sum(np.abs(b) ** 2, axis=ax)


def fine3_project_positive(psi, th):
    """project a 4-component field (S, 4, n, n, n) on the positive-energy band of the fine 3-D walk
    (eigenvectors of its two-tick Bloch matrix, zz_core.bloch_batch)"""
    n = psi.shape[-1]
    K = kgrid(n, 3)
    KA = np.stack([k.ravel() for k in K], 1)
    out = np.empty_like(psi)
    F = fftn(psi, 3).reshape(psi.shape[0], 4, -1)
    Fo = np.empty_like(F)
    ch = 65536
    for i0 in range(0, KA.shape[0], ch):
        W = zz.bloch_batch(KA[i0:i0 + ch], th)
        lam, V = np.linalg.eig(W)
        om = -np.angle(lam) / 2
        D = (om > 0).astype(float)
        Pm = np.einsum("nab,nb,nbc->nac", V, D, np.linalg.inv(V))
        Fo[:, :, i0:i0 + ch] = np.einsum("nab,sbn->san", Pm, F[:, :, i0:i0 + ch])
    return ifftn(Fo.reshape(psi.shape), 3)


def phi_of(psi4):
    """large Dirac component of a chiral (L, R) 4-spinor field"""
    return (psi4[:, 0:2] + psi4[:, 2:4]) / np.sqrt(2)


def fine3_positive_from_phi(phi, th):
    """a positive-energy fine 4-spinor whose large component is phi (used to build fine data and to prolong)"""
    n = phi.shape[-1]
    K = kgrid(n, 3)
    KA = np.stack([k.ravel() for k in K], 1)
    F = fftn(phi, 3).reshape(phi.shape[0], 2, -1)
    Fo = np.zeros((phi.shape[0], 4, F.shape[-1]), complex)
    ch = 65536
    for i0 in range(0, KA.shape[0], ch):
        W = zz.bloch_batch(KA[i0:i0 + ch], th)
        lam, V = np.linalg.eig(W)
        om = -np.angle(lam) / 2
        o = np.argsort(-om, axis=1)[:, :2]                     # the two positive-band eigenvectors
        U = np.take_along_axis(V, o[:, None, :], axis=2)       # (N, 4, 2)
        Phi = (U[:, 0:2, :] + U[:, 2:4, :]) / np.sqrt(2)       # (N, 2, 2) large component of the band basis
        coef = np.linalg.solve(Phi, F[:, :, i0:i0 + ch].transpose(2, 1, 0))   # (N, 2, S)
        Fo[:, :, i0:i0 + ch] = np.einsum("nab,nbs->san", U, coef)
    return ifftn(Fo.reshape((phi.shape[0], 4) + phi.shape[2:]), 3)


def smooth_field(rng, S, C, n, d, b, ks=1.2):
    """random complex field on a grid of n cells whose spectrum lives inside the COARSE Brillouin cube
    |k_i| < pi/b (no aliasing when blocked by b) with Gaussian weight exp(-(|k| b / ks)^2 / 2)"""
    K = kgrid(n, d)
    k2 = sum((k * b) ** 2 for k in K)
    inside = np.all([np.abs(k * b) < np.pi for k in K], axis=0)
    wgt = np.exp(-k2 / (2 * ks ** 2)) * inside
    f = (rng.normal(size=(S, C) + (n,) * d) + 1j * rng.normal(size=(S, C) + (n,) * d)) * wgt
    x = ifftn(f, d)
    return x / np.sqrt(np.sum(np.abs(x) ** 2, axis=tuple(range(1, x.ndim)), keepdims=True))


# ============================== levels and the automatic hierarchy ==============================
class Level:
    """one level of the hierarchy. B = cell in fine cells, tau = B fine ticks per step.
    kind: 'fine1' (1-D walk), 'fine3' (3-D walk), 'learned'. desc: 'RL' (1-D), 'dirac4', 'pauli2'."""

    def __init__(self, kind, d, desc, C, B, theta, rule=None, period=1, b=None, w=None, check=None):
        self.kind, self.d, self.desc, self.C, self.B, self.theta = kind, d, desc, C, B, theta
        self.rule, self.period, self.b, self.w, self.check = rule, period, b, w, check

    def evolve(self, psi, nsteps):
        if self.kind == "fine1":
            for _ in range(nsteps):
                psi = fine1_tick(psi, self.theta)
            return psi
        if self.kind == "fine3":
            for t in range(nsteps):
                psi = fine3_tick(psi, self.theta, t)
            return psi
        return self.rule.step(psi, nsteps)

    def positive(self, psi):
        """projection on the positive-energy band (3-D, 4-component descriptors)"""
        if self.kind == "fine3":
            return fine3_project_positive(psi, self.theta)
        n = psi.shape[-1]
        M = np.moveaxis(self.rule.mult(n), (0, 1), (-2, -1))             # (n.., 4, 4)
        lam, V = np.linalg.eig(M)
        wphi = np.sum(np.abs(V[..., 0:2, :] + V[..., 2:4, :]) ** 2, axis=-2) / 2 / np.sum(np.abs(V) ** 2, axis=-2)
        D = (wphi > 0.5).astype(float)
        Pm = np.einsum("...ab,...b,...bc->...ac", V, D, np.linalg.inv(V))
        F = fftn(psi, 3)
        return ifftn(np.einsum("...ab,sb...->sa...", Pm, F), 3)

    def mult_k(self, k, nsteps):
        """one-step multiplier of this level to the power nsteps at wave vectors k (N, d), units 1/cell of this level"""
        if self.kind == "fine1":
            return np.linalg.matrix_power(fine1_bloch(k[:, 0], self.theta), nsteps)
        if self.kind == "fine3":
            assert nsteps % 2 == 0
            return np.linalg.matrix_power(zz.bloch_batch(k, self.theta), nsteps // 2)
        return np.linalg.matrix_power(self.rule.mult_at(k), nsteps)

    def info(self):
        return dict(kind=self.kind, desc=self.desc, C=self.C, B=self.B, b=self.b, w=self.w, check=self.check,
                    taps=None if self.rule is None else int(len(self.rule.offs)))




def to_desc(prev, psi, desc):
    if desc == prev.desc:
        return psi
    assert desc == "pauli2" and prev.desc == "dirac4"
    return phi_of(psi)


def phi_weight(V):
    """large-component weight of the eigenvector columns V (..., 2 or 4, n): 1-D (R,L) or 3-D (L,R) chiral"""
    h = V.shape[-2] // 2
    return np.sum(np.abs(V[..., :h, :] + V[..., h:, :]) ** 2, axis=-2) / 2 / np.sum(np.abs(V) ** 2, axis=-2)


def conv_k(M, desc_from, desc_to):
    """multiplier in the descriptor desc_to: pauli2 from a 4-component one = positive band seen through phi"""
    if desc_from == desc_to:
        return M
    lam, V = np.linalg.eig(M)
    o = np.argsort(-phi_weight(V), axis=-1)[..., :2]
    U = np.take_along_axis(V, o[..., None, :], axis=-1)
    L = np.take_along_axis(lam, o, axis=-1)
    Phi = (U[..., 0:2, :] + U[..., 2:4, :]) / np.sqrt(2)
    return Phi @ (L[..., :, None] * np.linalg.inv(Phi))


def band_phase(M, desc):
    """phase of the positive band (mean over the spin pair in 3-D)"""
    lam, V = np.linalg.eig(M)
    if desc == "pauli2":
        return np.angle(np.mean(lam, axis=-1)) if True else None
    npos = 1 if desc == "RL" else 2
    o = np.argsort(-phi_weight(V), axis=-1)[..., :npos]
    return np.angle(np.mean(np.take_along_axis(lam, o, axis=-1), axis=-1))


def sample_state_k(rng, n, d, k_s):
    """wave vectors distributed like the momentum density of a hydrogen-like 1s state of size 1/k_s:
    rho(k) ~ 1/(1 + (k/k_s)^2)^4 (d = 3: radial pdf ~ u^2/(1+u^2)^4; d = 1: pdf ~ 1/(1+u^2)^4)"""
    u = np.linspace(0, 40, 200001)
    pdf = (u ** (d - 1)) / (1 + u * u) ** 4
    cdf = np.cumsum(pdf); cdf /= cdf[-1]
    r = np.interp(rng.random(n), cdf, u) * k_s
    v = rng.normal(size=(n, d)); v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v * r[:, None]


def kcheck(prev, cand, b, rng, T_obs, k_obs, sig=4.0, nk=3000, k0_cells=0.5, energy_mode="window"):
    """the CHECK, done exactly in k-space (linear, translation-invariant rules: identical to running packets):
    (i) packet error |dpsi|^2/|psi|^2 after H = ceil(T_obs / B) coarse steps, Gaussian packets sigma = sig cells,
        mean k in {0, +-0.5 per cell along axes and diagonal}, random spinor; worst case;
    (ii) energy error: max |error| of the positive-band motion energy eps(k) - eps(0) for k_obs/10 <= |k| <= k_obs,
        relative to eps(k_obs) (the energy scale of the observable)
        (physical k in 1/fine cell), plus the rest-phase error per step relative to eps(k_obs) tau."""
    d, B = prev.d, cand.B
    H = max(1, int(np.ceil(T_obs / B)))
    dirs = [np.zeros(d)] + [s * np.eye(d)[i] for i in range(d) for s in (1, -1)]
    if d == 3:
        dirs += [np.ones(3) / np.sqrt(3), -np.ones(3) / np.sqrt(3)]
    worst = 0.0
    for dv in dirs:
        k0 = k0_cells * dv
        kap = k0 + rng.normal(size=(nk, d)) / (2 * sig)
        Mc = np.linalg.matrix_power(cand.rule.mult_at(kap), H)
        Mp = conv_k(prev.mult_k(kap / b, b * H), prev.desc, cand.desc)
        D = Mc - Mp
        e = np.mean(np.sum(np.abs(D) ** 2, axis=(-2, -1))) / cand.C
        worst = max(worst, float(e))
    if energy_mode == "state":
        # energy error weighted by the observable's own momentum density (a 1s-like state of size 1/k_obs):
        # <|d eps|> / <eps>, and the rest-phase error in the same units
        kap = sample_state_k(rng, nk, d, k_obs * B)
        kap = kap[np.all(np.abs(kap) < np.pi, axis=1)]
        z = np.zeros((1, d))
        pc_all = band_phase(cand.rule.mult_at(np.concatenate([kap, z])), cand.desc)
        bp = band_phase(conv_k(prev.mult_k(np.concatenate([kap, z]) / b, b), prev.desc, cand.desc), cand.desc)
        pc = np.angle(np.exp(1j * (pc_all[:-1] - pc_all[-1])))
        pp = np.angle(np.exp(1j * (bp[:-1] - bp[-1])))
        mean_e = np.mean(np.abs(pp))
        rel = float(np.mean(np.abs(pc - pp)) / mean_e)
        rest = float(abs(np.angle(np.exp(1j * (pc_all[-1] - bp[-1])))) / mean_e)
        # whole Brillouin cube: the rule must stay tame where the observable does not look (relative error of the
        # motion phase, 90th percentile)
        ku = np.pi * (2 * rng.random((nk, d)) - 1)
        pcu = band_phase(cand.rule.mult_at(ku), cand.desc) - pc_all[-1]
        bpu = band_phase(conv_k(prev.mult_k(ku / b, b), prev.desc, cand.desc), cand.desc) - bp[-1]
        unw = lambda x: np.where(-np.angle(np.exp(1j * x)) < -0.5, -np.angle(np.exp(1j * x)) + 2 * np.pi, -np.angle(np.exp(1j * x)))
        pcu, bpu = unw(pcu), unw(bpu)
        bz = float(np.percentile(np.abs(pcu - bpu) / np.maximum(np.abs(bpu), 1e-12), 90))
        negmin = 0.0
        if isinstance(cand.rule, GenRule):     # no negative motion energy anywhere (no invented bound states)
            ev = np.linalg.eigvalsh(0.5 * (cand.rule.gen_at(ku) + np.conj(np.swapaxes(cand.rule.gen_at(ku), -1, -2))))
            negmin = float(max(0.0, -ev.min()) / mean_e)
        return dict(packet=worst, H=H, energy_rel=rel, rest_rel=rest, kmax_cells=float(k_obs * B), energy_mode="state",
                    bz_rel90=bz, neg_gen=negmin)
    # energy check
    kmax = k_obs * B                       # in cells of the candidate level
    kmax = min(kmax, 0.5 * np.pi)
    u = rng.normal(size=(nk, d)); u /= np.linalg.norm(u, axis=1, keepdims=True)
    r = kmax * (0.1 + 0.9 * rng.random(nk)) ** (1 / d)
    kap = u * r[:, None]
    z = np.zeros((1, d))
    pc = band_phase(cand.rule.mult_at(kap), cand.desc) - band_phase(cand.rule.mult_at(z), cand.desc)
    Mp = conv_k(prev.mult_k(np.concatenate([kap, z]) / b, b), prev.desc, cand.desc)
    bp = band_phase(Mp, cand.desc)
    pp = bp[:-1] - bp[-1]
    pc = np.angle(np.exp(1j * pc)); pp = np.angle(np.exp(1j * pp))
    rel = float(np.max(np.abs(pc - pp)) / np.max(np.abs(pp)))
    rest_c = band_phase(cand.rule.mult_at(z), cand.desc)[0]
    rest = float(abs(np.angle(np.exp(1j * (rest_c - bp[-1])))) / np.max(np.abs(pp)))
    return dict(packet=worst, H=H, energy_rel=rel, rest_rel=rest, kmax_cells=float(kmax))


def make_data(prev, b, desc, rng, Nc, S=None, ks=1.2):
    d = prev.d
    n = Nc * b
    S = S if S is not None else (6 if desc == "pauli2" else prev.C + 2)
    X = smooth_field(rng, S, prev.C, n, d, b, ks=ks)
    if desc != prev.desc:
        X = prev.positive(X)
    Y = prev.evolve(X, b)
    return block(to_desc(prev, X, desc), b, d), block(to_desc(prev, Y, desc), b, d)


def learn_k(prev, b, desc, w, ks, rng, nk=None, ridge=1e-12, frac_uniform=0.25, w_uniform=0.03):
    """least squares on PLANE-WAVE data: the level below applied to e^{ik.x} for many k (its multiplier to the power b,
    seen through the new descriptor) -- the same data as evolved packets (Parseval), but with no box: the k points are
    drawn densely where the observable lives (Gaussian of width ks per new cell) plus a uniform share over the whole
    Brillouin cube (keeps the rule tame where nothing is asked of it)."""
    d = prev.d
    offs = offsets(w, d)
    nt = len(offs)
    nk = nk if nk is not None else max(4000, 6 * nt)
    ng = int(nk * (1 - frac_uniform))
    if ks < 0:                          # negative ks = sample like the observable's state of size 1/|ks| cells
        kg = sample_state_k(rng, ng, d, -ks)
        kg = kg[np.all(np.abs(kg) < np.pi, axis=1)]
        ng = len(kg)
    else:
        kg = np.clip(rng.normal(size=(ng, d)) * ks, -np.pi, np.pi)
    ku = np.pi * (2 * rng.random((int(nk * frac_uniform), d)) - 1)
    kap = np.concatenate([kg, ku, np.zeros((1, d))])
    T = conv_k(prev.mult_k(kap / b, b), prev.desc, desc)                      # (nk, C, C)
    Phi = np.exp(1j * kap @ offs.T)                                           # (nk, nt)
    sw = np.ones(len(kap)); sw[ng:-1] = w_uniform                             # weights (sqrt) of the uniform share
    Phi = Phi * sw[:, None]
    T = T * sw[:, None, None]
    C = T.shape[-1]
    G = Phi.conj().T @ Phi
    lam = ridge * np.real(np.trace(G)) / nt
    rhs = Phi.conj().T @ T.reshape(len(kap), C * C)
    sol = np.linalg.solve(G + lam * np.eye(nt), rhs)                          # (nt, C*C)
    K = sol.T.reshape(C, C, nt)
    return Rule(K, offs, d, C, True)


def auto_level(prev, rng, Nc, eps_packet=1e-3, eps_energy=1e-4, T_obs=1000, k_obs=None, bs=(4, 3, 2),
               ws=(1, 2, 3, 4), descs=None, bmax_total=None, log=print, ks_rule=None, learner="packets",
               packets="cells", energy_mode="window"):
    """THE automatic rule: fewest components, then largest b, then smallest w, with both checks within budget"""
    if k_obs is None:
        k_obs = 0.3 * prev.theta
    if descs is None:
        descs = ["RL"] if prev.d == 1 else (["pauli2", "dirac4"] if prev.desc == "dirac4" else ["pauli2"])
    tried = []
    for desc in descs:
        C = 2 if desc in ("RL", "pauli2") else 4
        for b in bs:
            if (b % prev.period) or (bmax_total is not None and prev.B * b > bmax_total):
                continue
            t0 = time.time()
            ks = 1.2 if ks_rule is None else ks_rule(prev.B * b)
            if learner == "packets":
                Xc, Yc = make_data(prev, b, desc, rng, Nc, ks=ks)
            for w in ws:
                if learner == "packets":
                    if (2 * w + 1) > Nc // 2:
                        continue
                    r = learn(Xc, Yc, w, prev.d)
                    tr = float(rel_err(r.step(Xc), Yc).max())
                elif learner == "k" or (learner == "auto" and desc != "pauli2"):
                    r = learn_k(prev, b, desc, w, ks, rng)
                    tr = float("nan")
                else:
                    r = learn_h(prev, b, desc, w, ks, rng)
                    tr = float("nan")
                if isinstance(r, GenRule):
                    gmax = 1.0
                else:
                    r.unitary = False
                    gmax = float(np.abs(np.linalg.eigvals(r.mult_at(np.pi * (2 * rng.random((4000, prev.d)) - 1)))).max())
                    r.unitary = True
                cand = Level("learned", prev.d, desc, C, prev.B * b, prev.theta, rule=r, period=1, b=b, w=w)
                if packets == "cells":          # sigma = 4 new cells, |mean k| = 0.5 per new cell (A1)
                    ck = kcheck(prev, cand, b, rng, T_obs, k_obs, energy_mode=energy_mode)
                else:                           # physical packets of the observable: mean k = k_obs/2, spread k_obs/4
                    Bn = prev.B * b
                    ck = kcheck(prev, cand, b, rng, T_obs, k_obs, sig=max(1.0, 2.0 / (k_obs * Bn)), k0_cells=0.5 * k_obs * Bn,
                                energy_mode=energy_mode)
                ok = ck["packet"] <= eps_packet and ck["energy_rel"] <= eps_energy and ck["rest_rel"] <= eps_energy                     and ck.get("bz_rel90", 0.0) <= 0.2 and ck.get("neg_gen", 0.0) <= 0.01
                rec = dict(desc=desc, b=b, w=w, ks=ks, train=tr, raw_max_abs_eig=gmax, ok=bool(ok), sec=round(time.time() - t0, 1), **ck)
                tried.append(rec)
                log(f"    cand B={prev.B * b:5d} {desc:6s} b={b} w={w}: train {tr:.1e} packet({ck['H']} st) {ck['packet']:.1e} "
                    f"energy {ck['energy_rel']:.1e} rest {ck['rest_rel']:.1e} bz {ck.get('bz_rel90', 0):.1e} neg {ck.get('neg_gen', 0):.1e} raw|K| {gmax:.4f} {'OK' if ok else ''}")
                if ok:
                    cand.check = ck
                    return cand, tried
    return None, tried

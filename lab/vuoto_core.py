"""vuoto_core.py - the Dirac sea of the zig-zag rule in a static potential (vacuum polarization).

Study 'Vuoto che polarizza' (docs/VUOTO_POLARIZZA.md). Same rule as lab/zz_core.py (frozen), nothing new
in the rule: every ball moves at c, then rotates; the potential is the phase e^{-i V a} after the coin.

WHAT IS COMPUTED
The vacuum = the filled Dirac sea = every Floquet state of the walk with quasi-energy Omega in (-pi, 0)
(per Floquet period: one tick in 1-D, two ticks in 3-D because the 3-D rule alternates xyz / zyx).
A weak static potential V(x) = V_q e^{iqx} changes the sea projector P; the induced density is
dn(q) = Tr[dP e^{-iqx}] = -chi(q) V_q. dP comes from the exact first-order (Riesz) formula for a unitary:
    W(V) = Phi U_{T-1} ... Phi U_0,  Phi = e^{-iV}          (a = dt = 1 lattice units)
    dW_{k+q,k} = -i V_q B(k,q),  B = sum_t U_{T-1}(k+q)..U_{t+1}(k+q) U_t(k)..U_0(k)
    dP = sum_{n occ(k), m emp(k+q)} |m><m|dW|n><n| / (lam_n - lam_m)
       + sum_{m occ(k+q), n emp(k)} |m><m|dW|n><n| / (lam_m - lam_n)
(only occupied-empty pairs: safe for degenerate levels). Density measured at the end of the period
(stroboscopic). Linear in V = first order in Z alpha = exactly the Uehling order.

SPECIES (the 'doublers', declared): the rule has more than one low-energy species.
  1-D (one-tick period):  k = 0 (mass th, normal sea) and k = pi (mass th, gap at quasi-energy pi:
                          its sea is 'upside down' -> it responds with the OPPOSITE sign).  S U S = -U with
                          S = (-1)^x exactly, so the whole sea never polarizes (any V, all orders).
  3-D (two-tick period):  W(k + pi e_j) = W(k); in each pi-cell: 4 massive species at (0,0,0), (pi/2,pi/2,0)
                          and permutations (normal) and 4 MASSLESS cones at quasi-energy pi at (pi/2,0,0) perms
                          and (pi/2,pi/2,pi/2) (upside down -> opposite sign).
'Physical species' (the declared correction, like R1 'count only the k = 0 gap'): the response is computed
with a smooth species projector Q = sum_k sqrt(w(k)) |k><k| around k = 0 only: dn_phys = Tr[Q dP Q ...],
w(|k|) = 1 for |k| <= K1, 0 for |k| >= K2 (C-infinity step), the same w for every theta and every case.

Units inside: lattice units a = 1, so m = theta. Physical: q/m = q_lat/theta.
pit(q) := chi(q)/q^2 (3-D, per species, e = 1). Pi(q) = 4 pi alpha pit(q). Renormalized (on-shell, the same
scheme everywhere): Pi_R(q) = Pi(0) - Pi(q).  Uehling: Pi_R(q) = (2 alpha/pi) int_0^1 x(1-x) ln(1+x(1-x)q^2/m^2).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
from zz_core import BASIS_HALF, order_of_tick, ALPHA, MC2_EV, H_EVS

PI = np.pi
I2 = np.eye(2, dtype=complex)


# ------------------------------------------------------------------ one-tick Bloch matrices
def U1d(k, th):
    """1-D zig-zag (lab/zz_core.py kind 'dirac1'), one tick, basis (R, L): shift (R e^{-ik}, L e^{+ik}) then coin.
    k: array (N,) -> (N, 2, 2)"""
    k = np.asarray(k, float)
    c, s = np.cos(th), np.sin(th)
    U = np.zeros(k.shape + (2, 2), complex)
    eR, eL = np.exp(-1j * k), np.exp(1j * k)
    U[..., 0, 0] = c * eR
    U[..., 0, 1] = -1j * s * eL
    U[..., 1, 0] = -1j * s * eR
    U[..., 1, 1] = c * eL
    return U


def U3d(K, th, nt):
    """3-D zig-zag, ONE tick number nt (0: order x,y,z ; 1: z,y,x), basis (L0, L1, R0, R1) = zz_core.Walk('dirac').
    K: (N, 3) dimensionless k a -> (N, 4, 4). bloch_batch of zz_core = U3d(.,1) @ U3d(.,0)."""
    K = np.atleast_2d(K)
    N = K.shape[0]

    def Tax(j, h):
        lam, P = BASIS_HALF[j]
        ph = np.exp(-1j * h * lam[None, :] * K[:, j:j + 1])
        return np.einsum("ab,nb,cb->nac", P, ph, P.conj())
    TL = np.broadcast_to(I2, (N, 2, 2)).copy()
    TR = TL.copy()
    for j in order_of_tick(nt):
        TL = Tax(j, -1) @ TL
        TR = Tax(j, +1) @ TR
    ct, st = np.cos(th), np.sin(th)
    Wt = np.zeros((N, 4, 4), complex)
    Wt[:, :2, :2] = ct * TL
    Wt[:, :2, 2:] = -1j * st * TR
    Wt[:, 2:, :2] = -1j * st * TL
    Wt[:, 2:, 2:] = ct * TR
    return Wt


def floquet_seq(dim, th, period, start=0, mu=0.0):
    """list of one-tick Bloch functions U_t(k) for the Floquet period; start = 1 measures the density at the
    other time slice (3-D: period begins with the z,y,x tick). mu != 0: EXPLORATION B (not the rule): flip angle
    th - mu on even ticks, th + mu on odd ticks (gives the quasi-energy-pi cones a mass mu)."""
    if dim == 1:
        return [lambda k: U1d(k, th)] * period
    return [(lambda K, nt=nt: U3d(K, th + mu * (1 if (nt + start) % 2 else -1), (nt + start) % 2)) for nt in range(period)]


# ------------------------------------------------------------------ robust eigen-decomposition of a unitary
def ueig(W):
    """eigen-decomposition of (batched) unitary W whose spectrum avoids lambda = 1 (Omega = 0 is inside the mass
    gap). Cayley transform K = i (W + 1)(W - 1)^{-1} is Hermitian with eigenvalues cot(Omega/2), lambda = e^{-i Omega}.
    Returns (Omega in (-pi, pi), lam, V) with V[..., :, n] orthonormal eigenvectors; occupied = Omega < 0."""
    n = W.shape[-1]
    Id = np.broadcast_to(np.eye(n, dtype=complex), W.shape)
    # K = i (W+1)(W-1)^{-1}  ->  solve (W-1)^T X^T = (W+1)^T  (right division)
    A = np.swapaxes(W - Id, -1, -2)
    Bm = np.swapaxes(W + Id, -1, -2)
    K = 1j * np.swapaxes(np.linalg.solve(A, Bm), -1, -2)
    K = 0.5 * (K + np.conj(np.swapaxes(K, -1, -2)))
    t, V = np.linalg.eigh(K)
    Om = 2 * np.arctan2(1.0, -t)           # in (0, 2pi): cot(Om/2) = -t
    Om = np.where(Om > PI, Om - 2 * PI, Om)
    return Om, np.exp(-1j * Om), V


# ------------------------------------------------------------------ the Floquet linear response kernel
def lr_terms(Ks, q, dim, th, period, weight=None, start=0, mu=0.0):
    """integrand F(k) of chi(q) = sum_k F(k) / N (lattice units). Ks: (N,) in 1-D or (N,3) in 3-D; q scalar (1-D)
    or 3-vector. weight: None (whole sea) or function k -> w(k) in [0,1] (species projector, sqrt(w(k) w(k+q))).
    Returns complex array F (N,); chi(q) = mean(F) * (BZ volume factor handled by the caller)."""
    seq = floquet_seq(dim, th, period, start, mu)
    Kq = Ks + q
    Ua = [f(Ks) for f in seq]          # U_t(k)
    Ub = [f(Kq) for f in seq]          # U_t(k+q)
    nb = Ua[0].shape[-1]
    W0 = Ua[0].copy()
    for t in range(1, period):
        W0 = Ua[t] @ W0
    W1 = Ub[0].copy()
    for t in range(1, period):
        W1 = Ub[t] @ W1
    # B = sum_t [U_{T-1}(k+q) .. U_{t+1}(k+q)] [U_t(k) .. U_0(k)]
    B = np.zeros_like(W0)
    for t in range(period):
        left = np.broadcast_to(np.eye(nb, dtype=complex), W0.shape).copy()
        for s in range(t + 1, period):
            left = Ub[s] @ left
        right = np.broadcast_to(np.eye(nb, dtype=complex), W0.shape).copy()
        for s in range(0, t + 1):
            right = Ua[s] @ right
        B = B + left @ right
    Om0, l0, V0 = ueig(W0)
    Om1, l1, V1 = ueig(W1)
    Bmn = np.einsum("nam,nab,nbj->nmj", V1.conj(), B, V0)        # <m,k+q| B |n,k>
    Snm = np.einsum("naj,nam->njm", V0.conj(), V1)                 # <n,k|m,k+q>
    occ0 = Om0 < 0
    occ1 = Om1 < 0
    # pair (n at k, m at k+q): factor 1/(l_n - l_m) if n occ & m emp; 1/(l_m - l_n) if m occ & n emp; else 0
    num = Bmn * np.swapaxes(Snm, -1, -2)                          # [n, m, j] -> want [m, j]: Bmn[m,j] Snm[j,m]
    num = Bmn * np.swapaxes(Snm, -1, -2)
    den = l0[:, None, :] - l1[:, :, None]                           # [m, j] = l_j(k) - l_m(k+q)
    A1 = occ0[:, None, :] & (~occ1[:, :, None])
    A2 = (~occ0[:, None, :]) & occ1[:, :, None]
    fac = np.where(A1, 1.0, 0.0) - np.where(A2, 1.0, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(fac != 0, fac * num / den, 0.0)
    F = 1j * T.sum(axis=(1, 2))                                     # chi = -dn/V,  dn = -i V sum(...)
    if weight is not None:
        F = F * np.sqrt(weight(Ks) * weight(Kq))
    return F


# ------------------------------------------------------------------ species weight (declared, the same everywhere)
K1_DEF, K2_DEF = 0.6, 1.05


def smooth_step(x):
    """C-infinity step: 1 for x <= 0, 0 for x >= 1"""
    x = np.clip(x, 0.0, 1.0)
    with np.errstate(divide="ignore", over="ignore"):
        a = np.where(x < 1, np.exp(-1.0 / np.maximum(1 - x, 1e-300)), 0.0)
        b = np.where(x > 0, np.exp(-1.0 / np.maximum(x, 1e-300)), 0.0)
    return np.where(x <= 0, 1.0, np.where(x >= 1, 0.0, a / (a + b)))


def species_weight(K1=K1_DEF, K2=K2_DEF):
    def w(Ks):
        r = np.abs(Ks) if np.ndim(Ks) == 1 else np.linalg.norm(Ks, axis=-1)
        return smooth_step((r - K1) / (K2 - K1))
    return w


# ------------------------------------------------------------------ quadratures
def sinh_rule(n, kc, beta):
    """Gauss-Legendre on t in [-1,1], k = kc sinh(beta t)/sinh(beta): nodes and weights on [-kc, kc]"""
    t, wt = roots_legendre(n)
    k = kc * np.sinh(beta * t) / np.sinh(beta)
    dk = kc * beta * np.cosh(beta * t) / np.sinh(beta)
    return k, wt * dk


def chi_1d(q, th, weight=None, n=4000, kc=None, beta=None, period=1):
    """chi(q) for the 1-D rule; weight None -> whole Brillouin zone (uniform grid, exact periodic trapezoid)"""
    if weight is None:
        k = -PI + (np.arange(n) + 0.5) * 2 * PI / n
        F = lr_terms(k, q, 1, th, period)
        return np.mean(F)                                         # (1/N) sum_k
    kc = kc if kc is not None else K2_DEF + abs(q)
    k, w = sinh_rule(n, kc, beta if beta is not None else 3.0)
    F = lr_terms(k, q, 1, th, period, weight)
    return np.sum(F * w) / (2 * PI)


def chi_3d(q, th, weight=None, n=64, kc=None, beta=None, chunk=60000, start=0, mu=0.0):
    """chi(q) (q: 3-vector, lattice units) for the 3-D rule (two-tick Floquet period), physical species with weight;
    product sinh-Gauss rule on [-kc, kc]^3. Returns complex chi (per site = per volume a^3)."""
    q = np.asarray(q, float)
    kc = kc if kc is not None else K2_DEF + np.linalg.norm(q) + 0.02
    k1, w1 = sinh_rule(n, kc, beta if beta is not None else 3.5)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    WW = (w1[:, None, None] * w1[None, :, None] * w1[None, None, :]).ravel()
    Ks = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], 1) - 0.5 * q[None, :]   # centre the pair (k, k+q)
    tot = 0.0 + 0.0j
    for s in range(0, len(Ks), chunk):
        F = lr_terms(Ks[s:s + chunk], q, 3, th, 2, weight, start, mu)
        tot += np.sum(F * WW[s:s + chunk])
    return tot / (2 * PI) ** 3


# ------------------------------------------------------------------ continuum references (EXAM ONLY)
def chi1d_cont(q, m=1.0):
    """1+1 D massive Dirac, static density response per species: (1/pi) int x(1-x) q^2/(m^2 + x(1-x) q^2)"""
    return quad(lambda x: x * (1 - x) * q * q / (m * m + x * (1 - x) * q * q), 0, 1, epsabs=1e-15, epsrel=1e-13)[0] / PI


def chi1d_cont_direct(q, m=1.0, kmax=4000.0):
    """same, from the Dirac Hamiltonian h(k) = k sz + m sx (independent check of the formula above)"""
    def E(k):
        return np.sqrt(k * k + m * m)

    def f(k):
        e1, e2 = E(k), E(k + q)
        # |<+,k+q|-,k>|^2 = (1 - (k(k+q) + m^2)/(e1 e2))/2
        ov = 0.5 * (1 - (k * (k + q) + m * m) / (e1 * e2))
        return 2 * ov / (e1 + e2)
    return quad(f, -kmax, kmax, limit=2000, points=[0.0, -q], epsabs=1e-14)[0] / (2 * PI)


def pitR_uehling(q, m=1.0):
    """renormalized polarization per species / (4 pi alpha): (1/(2 pi^2)) int x(1-x) ln(1 + x(1-x) q^2/m^2)"""
    return quad(lambda x: x * (1 - x) * np.log1p(x * (1 - x) * q * q / (m * m)), 0, 1, epsabs=1e-16, epsrel=1e-13)[0] / (2 * PI**2)


LOGCOEF = 1.0 / (6 * PI**2)        # d pit(0) / d(-ln theta)  <->  2 alpha/(3 pi) after 4 pi alpha
C2_UEHLING = 1.0 / (60 * PI**2)    # pitR(q) ~ C2 q^2/m^2 (q -> 0)   <->  alpha/(15 pi)


def uehling_V(r, Za=ALPHA, m=1.0, alpha=ALPHA):
    """Uehling potential energy of the electron (units m c^2, r in hbar/mc): -(Za/r)(2 alpha/3pi) int_1^inf ..."""
    f = quad(lambda t: np.exp(-2 * m * r * t) * (1 + 1 / (2 * t * t)) * np.sqrt(t * t - 1) / (t * t), 1, np.inf,
             epsabs=1e-17, epsrel=1e-12, limit=400)[0]
    return -(Za / r) * (2 * alpha / (3 * PI)) * f


def V_from_pitR(r, qs, pitR, Za, alpha=ALPHA, damp=None):
    """r-space potential from a tabulated pitR(q) (q-grid qs, units m): V(r) = -(2 Za/(pi r)) int sin(qr) Pi_R(q)/q dq,
    Pi_R = 4 pi alpha pitR. Trapezoid on the given grid; optional damping factor array."""
    PiR = 4 * PI * alpha * np.asarray(pitR)
    if damp is not None:
        PiR = PiR * damp
    integ = np.sin(np.outer(r, qs)) * (PiR / qs)[None, :]
    return -(2 * Za / (PI * r)) * np.trapezoid(integ, qs, axis=1)


def MHz(E_mc2):
    return E_mc2 * MC2_EV / H_EVS / 1e6


# ------------------------------------------------------------------ self-tests
def _selftest():
    import time
    t0 = time.time()
    rng = np.random.default_rng(7)
    # (1) bloch_batch identity
    from zz_core import bloch_batch
    K = rng.uniform(-PI, PI, (50, 3))
    print("U3d(1) U3d(0) == zz_core.bloch_batch:", np.max(np.abs(U3d(K, 0.3, 1) @ U3d(K, 0.3, 0) - bloch_batch(K, 0.3))))
    # (2) Riesz first-order projector formula vs finite differences on a random unitary with a gap at Omega = 0
    n = 6
    H = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); H = (H + H.conj().T) / 2
    ev, Vh = np.linalg.eigh(H)
    ev = np.where(np.abs(ev) < 0.3, np.sign(ev) * 0.3, ev); ev = np.clip(ev, -2.8, 2.8)
    U0 = Vh @ np.diag(np.exp(-1j * ev)) @ Vh.conj().T
    U1 = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))[0]
    Vp = np.diag(rng.normal(size=n))
    Wfull = U1 @ U0      # walk Phi U1 Phi U0 ; make W's gap: use W directly
    def P_of(eps):
        Ph = np.diag(np.exp(-1j * eps * np.diag(Vp)))
        W = Ph @ U1 @ Ph @ U0
        Om, lam, V = ueig(W[None])
        V = V[0]; o = Om[0] < 0
        return V[:, o] @ V[:, o].conj().T, Om[0]
    eps = 1e-6
    dP_fd = (P_of(eps)[0] - P_of(-eps)[0]) / (2 * eps)
    Om, lam, V = ueig((U1 @ U0)[None]); Om, lam, V = Om[0], lam[0], V[0]
    B = U1 @ U0 + U1 @ U0 * 0     # placeholder replaced below
    dW = -1j * (Vp @ U1 @ U0 + U1 @ Vp @ U0)
    dP = np.zeros((n, n), complex)
    for a in range(n):
        for b in range(n):
            if Om[a] < 0 and Om[b] >= 0:      # a occ (n), b emp (m)
                dP += np.outer(V[:, b], V[:, b].conj()) @ dW @ np.outer(V[:, a], V[:, a].conj()) / (lam[a] - lam[b])
            if Om[b] < 0 and Om[a] >= 0:      # b occ (m), a emp (n)
                dP += np.outer(V[:, b], V[:, b].conj()) @ dW @ np.outer(V[:, a], V[:, a].conj()) / (lam[b] - lam[a])
    print("Riesz dP vs finite difference (random 6x6, two-step walk):", np.max(np.abs(dP - dP_fd)), " min|Omega|", np.min(np.abs(Om)))
    print("time %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    _selftest()


# ------------------------------------------------------------------ whole sea in 3-D (all species, no projector)
def special_rule(nseg):
    """quadrature on one axis of the pi-cell [-pi/2, pi/2): 4 segments, each concentrated (sinh map) at the special
    coordinate (-pi/2, 0, 0, +pi/2) where the species / cones sit"""
    s, ws = roots_legendre(nseg)
    s = 0.5 * (s + 1); ws = 0.5 * ws               # [0, 1]
    beta = 6.0
    g = np.sinh(beta * s) / np.sinh(beta)
    dg = beta * np.cosh(beta * s) / np.sinh(beta)
    L = PI / 4
    nodes, wts = [], []
    for end, sgn in ((-PI / 2, +1), (0.0, -1), (0.0, +1), (PI / 2, -1)):
        nodes.append(end + sgn * L * g); wts.append(L * dg * ws)
    return np.concatenate(nodes), np.concatenate(wts)


def chi_full_3d(q, th, nseg=16, chunk=60000, mu=0.0):
    """chi per site of the WHOLE two-tick sea (all 32 massive species + 32 massless cones), divided by 32 (= per
    massive species). Integral over one pi-cell (W is pi-periodic): chi_site = 8 int_cell F / (2pi)^3."""
    q = np.asarray(q, float)
    k1, w1 = special_rule(nseg)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    WW = (w1[:, None, None] * w1[None, :, None] * w1[None, None, :]).ravel()
    Ks = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], 1) - 0.5 * q[None, :]
    tot = 0.0 + 0.0j
    for s in range(0, len(Ks), chunk):
        F = lr_terms(Ks[s:s + chunk], q, 3, th, 2, None, 0, mu)
        tot += np.sum(F * WW[s:s + chunk])
    return 8 * tot / (2 * PI) ** 3 / 32

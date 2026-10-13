"""senzadoppi_core.py - tools of the study "Una regola zig-zag senza doppioni?" (docs/SENZA_DOPPIONI.md).

Everything here is a VARIANT of the frozen rule (lab/zz_core.py is imported, never modified):
"corri a c lungo un asse alla volta (la parte con spin +1 lungo l'asse va avanti, -1 indietro), poi ruota".
What can change between variants (and nothing else):
  * the ORDER of the three axis steps, tick by tick, separately for the right-handed part R and the left-handed L;
  * an optional twisted Wilson term (REGOLA_LOCALE F2: link rotations on +-2 e_j around the chirality axis Y,
    with the on-site part removed, so that at k = 0 it is exactly the identity);
  * an optional "spatially staggered" twisted Wilson (link rotation with the sign (-1)^x on the link origin)
    -> it transfers momentum Q = (pi, 0, 0)-type; used only as an exploration;
  * the lattice of the steps: axis steps (cubic) or body-diagonal steps (bcc) as in locale_core.
All analysis is through the exact Bloch matrix of one PERIOD (a translation-invariant linear rule acts on a plane wave
by this matrix; checked against the real-space runner in senzadoppi_r1.py).

Units: c = hbar = 1, lattice spacing a = 1 (so k is k*a, energies are E*a/(hbar c), theta = m c a / hbar).
Basis of the 4-spinor: (L0, L1, R0, R1) as in zz_core.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import itertools
import numpy as np

import zz_core as zz

PI = np.pi
SX, SY, SZ, I2 = zz.SX, zz.SY, zz.SZ, zz.I2
PAULI = np.array(zz.PAULI)
PERMS = list(itertools.permutations((0, 1, 2)))
PNAME = {p: "".join("xyz"[j] for j in p) for p in PERMS}
YCHI = np.kron(SY, I2)          # chirality Y on (L, R)
XCHI = np.kron(SX, I2)          # chirality X = the mass flip axis of zz_core


def axis_step(K, j, s):
    """exp(-i s k_j sigma_j) for a batch K (N,3); s = +1 (R) or -1 (L)"""
    c, sn = np.cos(K[:, j]), np.sin(K[:, j])
    return c[:, None, None] * I2[None] - 1j * s * sn[:, None, None] * PAULI[j][None]


def word(K, order, s):
    """2x2 Bloch matrix of one tick of translations, axes applied in `order` (first acts first)"""
    N = K.shape[0]
    M = np.broadcast_to(I2, (N, 2, 2)).copy()
    for j in order:
        M = axis_step(K, j, s) @ M
    return M


def diag_step(K, d, s):
    """bcc body-diagonal step along d in {(+-1,+-1,+-1)} with spin matrix sigma.d/sqrt3 (eigenvalues +-1):
    the +1 part moves by +d (integer vector), the -1 part by -d.  exp(-i s (k.d) sigma.n)"""
    d = np.asarray(d, float)
    n = d / np.linalg.norm(d)
    S = np.einsum("i,ijk->jk", n, PAULI)
    kd = K @ d
    return np.cos(kd)[:, None, None] * I2[None] - 1j * s * np.sin(kd)[:, None, None] * S[None]


def coin(theta):
    return np.cos(theta) * np.eye(4) - 1j * np.sin(theta) * XCHI


def _blk(L, R):
    N = L.shape[0]
    T = np.zeros((N, 4, 4), complex)
    T[:, :2, :2], T[:, 2:, 2:] = L, R
    return T


def wilson_links(K, h, gen=YCHI, vecs=None, stagger=None):
    """product of link rotations S_v^dag e^{-i a G} S_v e^{+i a G}, v over +-2e_j (6 links), a = h/2 each;
    S_v = P_L + P_R e^{-i k.v} (R part sits at x + v).  First order generator: -h sum_j (1 - cos 2k_j) G.
    stagger: None, or an integer vector Q/pi -> the R partner is taken at momentum k + Q (sign (-1)^{Q.x} on links):
    implemented as S_v = P_L + P_R e^{-i (k+Q).v}... only meaningful in the extended-zone bookkeeping; see doc."""
    N = K.shape[0]
    if vecs is None:
        vecs = [2 * e for e in np.eye(3)] + [-2 * e for e in np.eye(3)]
    a = h / 2.0
    ea = np.cos(a) * np.eye(4) - 1j * np.sin(a) * gen          # e^{-i a G}
    eb = np.cos(a) * np.eye(4) + 1j * np.sin(a) * gen
    M = np.broadcast_to(np.eye(4, dtype=complex), (N, 4, 4)).copy()
    for v in vecs:
        ph = np.exp(-1j * (K @ np.asarray(v, float)))
        S = np.zeros((N, 4, 4), complex)
        S[:, 0, 0] = S[:, 1, 1] = 1
        S[:, 2, 2] = S[:, 3, 3] = ph
        F = np.conj(np.transpose(S, (0, 2, 1))) @ ea[None] @ S @ eb[None]
        M = F @ M
    return M


class Rule:
    """A period of `len(ticks)` ticks. ticks[n] = (orderR, orderL) or ('bcc', dirsR, dirsL).
    After each tick: [half Wilson] coin(theta) [half Wilson]  (Wilson only if h != 0; the order of the links is
    mirrored on odd ticks as in REGOLA_LOCALE F2).  T = number of ticks = physical duration (a/c units)."""

    def __init__(self, ticks, h=0.0, name="", Lsign=-1):
        self.ticks, self.h, self.name, self.Lsign = ticks, h, name, Lsign
        self.T = len(ticks)

    def tick_trans(self, K, n):
        t = self.ticks[n]
        if t[0] == "bcc":
            N = K.shape[0]
            R = np.broadcast_to(I2, (N, 2, 2)).copy()
            L = R.copy()
            for d in t[1]:
                R = diag_step(K, d, +1) @ R
            for d in t[2]:
                L = diag_step(K, d, self.Lsign) @ L
            return _blk(L, R)
        return _blk(word(K, t[1], self.Lsign), word(K, t[0], +1))

    def bloch(self, K, theta):
        K = np.atleast_2d(np.asarray(K, float))
        N = K.shape[0]
        W = np.broadcast_to(np.eye(4, dtype=complex), (N, 4, 4)).copy()
        C = coin(theta)[None]
        for n in range(self.T):
            W = self.tick_trans(K, n) @ W
            if self.h != 0:
                vecs = [2 * e for e in np.eye(3)] + [-2 * e for e in np.eye(3)]
                if n % 2 == 1:
                    vecs = vecs[::-1]
                Wh = wilson_links(K, self.h / 2, vecs=vecs)
                Wh2 = wilson_links(K, self.h / 2, vecs=vecs[::-1])
                W = Wh2 @ C @ Wh @ W
            else:
                W = C @ W
        return W

    def bloch_R(self, K):
        """2x2 massless right-handed part of the period"""
        K = np.atleast_2d(np.asarray(K, float))
        return self.bloch(K, 0.0)[:, 2:, 2:] if self.h == 0 else None


def quasi(W, T):
    """physical quasi-energies (per unit time) of batch W (N,d,d), sorted, in (-pi/T, pi/T]"""
    w = np.linalg.eigvals(W)
    return np.sort(-np.angle(w), axis=1) / T


# --------------------------------------------------------------------- copies (exact sublattice decoupling)
def copy_shifts(rule, theta=0.37, tol=1e-12):
    """G in {0, pi}^3 with W(k + G) = W(k) for all k (tested on random k): exact copies (decoupled sublattices)"""
    rng = np.random.default_rng(1)
    K = rng.uniform(-PI, PI, (6, 3))
    W0 = rule.bloch(K, theta)
    out = []
    for g in itertools.product((0, 1), repeat=3):
        G = PI * np.array(g, float)
        if np.max(np.abs(rule.bloch(K + G, theta) - W0)) < tol:
            out.append(np.array(g))
    return out


# --------------------------------------------------------------------- 2x2 node finder (Weyl, theta = 0)
def _coef(M):
    return np.array([np.real(np.trace(M)) / 2] + [np.real(np.trace(-1j * M @ S)) / 2 for S in PAULI])


def weyl_nodes(rule, n_starts=600, seed=0):
    """nodes of the massless R part (W_R = +-1 per period) in the torus [-pi, pi)^3. Returns list of
    (k, sign(+1 -> quasi-energy 0, -1 -> pi), chirality = local degree)."""
    from scipy.optimize import least_squares
    rng = np.random.default_rng(seed)

    def WR(k):
        return rule.bloch(k[None], 0.0)[0, 2:, 2:]

    def ldeg(k):
        e = 1e-6
        J = np.zeros((4, 3))
        for i in range(3):
            dk = np.zeros(3)
            dk[i] = e
            J[:, i] = (_coef(WR(k + dk)) - _coef(WR(k - dk))) / (2 * e)
        return int(np.sign(np.linalg.det(np.column_stack([_coef(WR(k)), J]))))
    sols = []
    for sgn in (+1, -1):
        tgt = np.array([sgn, 0, 0, 0.0])
        g = PI / 4 * np.arange(-4, 4)
        seeds = [np.array(x) + 1e-3 * rng.normal(size=3) for x in itertools.product(g, g, g)]
        seeds += [rng.uniform(-PI, PI, 3) for _ in range(n_starts)]
        for x0 in seeds:
            r = least_squares(lambda k: _coef(WR(k)) - tgt, x0, xtol=1e-14, ftol=1e-14,
                              gtol=1e-14)
            if np.linalg.norm(r.fun) < 1e-8:
                k = (r.x + PI) % (2 * PI) - PI
                if not any(s[1] == sgn and np.linalg.norm((k - s[0] + PI) % (2 * PI) - PI) < 1e-5 for s in sols):
                    sols.append((k, sgn, ldeg(k)))
    return sols


# --------------------------------------------------------------------- Dirac valley scan
def zone_grid(n):
    g = (np.arange(n) + 0.5) / n * 2 * PI - PI
    return np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T


def valley_scan(rule, theta, n=40, kc=0.45, refine=True):
    """minimum |E| (gap to 0) and minimum distance to the zone edge pi/T, over the torus EXCLUDING balls of radius kc
    around Gamma and its exact copies. Physical energies (units 1/a). Returns dict with E0min, k0, Epimin, kpi."""
    from scipy.optimize import minimize
    shifts = [PI * g for g in copy_shifts(rule)]
    K = zone_grid(n)

    def dist_gamma(K):
        d = np.full(K.shape[0], np.inf)
        for G in shifts:
            D = (K - G + PI) % (2 * PI) - PI
            d = np.minimum(d, np.linalg.norm(D, axis=1))
        return d
    E = quasi(rule.bloch(K, theta), rule.T)
    mask = dist_gamma(K) > kc
    e0 = np.min(np.abs(E), axis=1)
    epi = PI / rule.T - np.max(np.abs(E), axis=1)
    out = {}
    for key, arr in (("E0", e0), ("Epi", epi)):
        a = np.where(mask, arr, np.inf)
        idx = np.argsort(a)[:6]
        best = (a[idx[0]], K[idx[0]])
        if refine:
            for i in idx:
                def f(k):
                    if dist_gamma(k[None])[0] <= kc:
                        return 10.0
                    Ek = quasi(rule.bloch(k[None], theta), rule.T)[0]
                    return np.min(np.abs(Ek)) if key == "E0" else PI / rule.T - np.max(np.abs(Ek))
                r = minimize(f, K[i], method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 3000})
                if r.fun < best[0]:
                    best = (r.fun, (r.x + PI) % (2 * PI) - PI)
        out[key] = float(best[0])
        out["k_" + key] = [float(x) for x in best[1] / PI]
    out["n_copies"] = len(shifts)
    return out


# --------------------------------------------------------------------- isotropy of the massless cone
def fib_dirs(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = PI * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def isotropy(rule, cells_per_lambda=24, ndir=400):
    """phase-velocity error max |omega/(c k) - 1| of the massless bands at |k| = 2 pi / N (all 4 bands, theta = 0)"""
    kk = 2 * PI / cells_per_lambda
    D = fib_dirs(ndir)
    E = quasi(rule.bloch(kk * D, 0.0), rule.T)
    v = np.abs(E) / kk
    return float(np.max(np.abs(v - 1))), float(np.mean(v)), float(np.min(v)), float(np.max(v))


def rest_and_kinetic(rule, theta, dk=1e-3, ndir=26):
    """E(0)/m - 1 and kinetic mass error m_kin/m - 1 (max over directions) of the electron band at Gamma.
    For E(k) = sqrt(m^2 + k^2): E(k) - E(0) = k^2/(2m).  Uses dk small; m = theta."""
    m = theta
    E0 = quasi(rule.bloch(np.zeros((1, 3)), theta), rule.T)[0]
    pos0 = np.min(E0[E0 > 0])
    D = fib_dirs(ndir)
    Ek = quasi(rule.bloch(dk * D, theta), rule.T)
    pos = np.array([np.min(e[e > 0]) for e in Ek])
    mk = dk**2 / (2 * (pos - pos0))
    return float(pos0 / m - 1), float(np.max(np.abs(mk / m - 1))), float(np.mean(mk / m - 1))

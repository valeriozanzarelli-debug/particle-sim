"""fw1vuoto_core.py - vacuum polarization of the FW1 variant WITHOUT the species projector (docs/FW1_VUOTO.md).

What is new here (and what is copied).
  * The period of a rule is written as a list of MOVES. A move is a finite sum of hops: U(k) = sum_u M_u e^{-i k.u}
    (u = integer hop vector; e^{-ik.u} = "moves by +u", the convention of zz_core / senzadoppi_core).
    Frozen rule (zz_core) and FW1 (senzadoppi_r5.FW1: frozen + ONE twisted Wilson, links +-2e_j, h = 0,4, half before
    and half after the coin of the EVEN ticks) are built from the same moves; checked against zz_core.bloch_batch and
    senzadoppi_r5.FW1.bloch (selftest).
  * Static perturbations of wavevector q, always gauge-covariant:
      'V'      scalar potential: the phase e^{-iV} of the rule, one move PHI after every tick (after the coin; for FW1
               tick 0 after the second half of the Wilson: declared choice);
      ('A',j)  vector potential along axis j: Peierls phase e^{i int A.dl} on EVERY hop along j (axis steps AND the
               +-2e_j Wilson links), exact line integral of the plane wave (so a pure gauge is exactly a conjugation).
  * Two responses:
      strobo_F  = the stroboscopic density response of lab/vuoto_core.lr_terms (Debiti, origin/master 48fbe35),
                  copied and generalised to moves (Riesz first-order projector, density at the end of the period);
      energy_F  = NEW: second-order change of the sea quasi-energy  E = sum_{occupied} Omega_n  (Omega in (-pi, 0)),
                  E = Tr f(W), f = i log on the occupied arc, 0 on the empty arc; exact second-order formula for a
                  normal matrix (divided differences of f'). Gives chi (V-V), the magnetic kernel K_jk (A-A) and the
                  mixed V-A kernel with the same bookkeeping, so that e and 1/mu of the vacuum are comparable.
Normalisation (as Debiti): integrals d^3k/(2 pi)^3 per site; 'per species' = one pi-cell of the zone (FW1: 1 light
species per pi-cell; frozen: 4 massive species + 4 massless cones per pi-cell).
Units: a = 1, c = 1, m = theta.  pit(q) = chi(q)/q^2;  Pi = 4 pi alpha pit.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad

PI = np.pi
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
PAULI = [SX, SY, SZ]
I4 = np.eye(4, dtype=complex)
PL = np.diag([1, 1, 0, 0]).astype(complex)
PR = np.diag([0, 0, 1, 1]).astype(complex)
YCHI = np.kron(SY, I2)
XCHI = np.kron(SX, I2)
E3 = np.eye(3, dtype=int)
VEC_A = [2 * E3[j] for j in range(3)] + [-2 * E3[j] for j in range(3)]     # senzadoppi_r1.VEC_A
LOGCOEF = 1.0 / (6 * PI**2)


def blk(L, R):
    T = np.zeros((4, 4), complex)
    T[:2, :2], T[2:, 2:] = L, R
    return T


# ------------------------------------------------------------------ moves
def m_step(j):
    """one axis step along j: R part exp(-i k_j sigma_j), L part exp(+i k_j sigma_j)"""
    Pp, Pm = (I2 + PAULI[j]) / 2, (I2 - PAULI[j]) / 2
    e = tuple(E3[j])
    me = tuple(-E3[j])
    return {"hops": [(e, blk(Pm, Pp)), (me, blk(Pp, Pm))], "phi": False}


def m_const(M):
    return {"hops": [((0, 0, 0), np.asarray(M, complex))], "phi": False}


def m_phi():
    return {"hops": [((0, 0, 0), I4.copy())], "phi": True}


def coin(th):
    return np.cos(th) * I4 - 1j * np.sin(th) * XCHI


def m_link(v, a, G=YCHI):
    """F_v = S_v^dag e^{-iaG} S_v e^{+iaG}, S_v = P_L + P_R T_v (R part hops by v): merged hops 0, +v, -v"""
    Ea = np.cos(a) * I4 - 1j * np.sin(a) * G
    Eb = np.cos(a) * I4 + 1j * np.sin(a) * G
    v = tuple(int(x) for x in v)
    mv = tuple(-x for x in v)
    return {"hops": [((0, 0, 0), PL @ Ea @ PL @ Eb + PR @ Ea @ PR @ Eb),
                     (v, PL @ Ea @ PR @ Eb), (mv, PR @ Ea @ PL @ Eb)], "phi": False}


ORDERS = [(0, 1, 2), (2, 1, 0)]


def frozen_moves(th):
    mv = []
    for n in range(2):
        mv += [m_step(j) for j in ORDERS[n]]
        mv += [m_const(coin(th)), m_phi()]
    return mv


def fw1_moves(th, h=0.4, phi_mid=False):
    """FW1 (senzadoppi_r5): tick 0 = steps xyz, Wilson half (links VEC_A, angle h/4 each), coin, Wilson half (links
    reversed); tick 1 = steps zyx, coin. phi_mid: put the tick-0 potential between the coin and the 2nd Wilson half."""
    a = h / 4.0
    mv = [m_step(j) for j in ORDERS[0]]
    mv += [m_link(v, a) for v in VEC_A]
    mv += [m_const(coin(th))]
    if phi_mid:
        mv += [m_phi()]
    mv += [m_link(v, a) for v in VEC_A[::-1]]
    if not phi_mid:
        mv += [m_phi()]
    mv += [m_step(j) for j in ORDERS[1]]
    mv += [m_const(coin(th)), m_phi()]
    return mv


RULES = {"F0": frozen_moves, "FW1": fw1_moves}
NTICK = 2


# ------------------------------------------------------------------ evaluation of moves and vertices
def g_line(s):
    """(e^{is} - 1)/(is), exact line integral of a plane wave along a hop"""
    s = np.asarray(s, float)
    out = np.ones_like(s, dtype=complex)
    nz = np.abs(s) > 1e-12
    out[nz] = (np.exp(1j * s[nz]) - 1) / (1j * s[nz])
    out[~nz] = 1 + 0.5j * s[~nz]
    return out


def ev_move(mv, K):
    N = K.shape[0]
    U = np.zeros((N, 4, 4), complex)
    for u, M in mv["hops"]:
        if u == (0, 0, 0):
            U += M[None]
        else:
            U += np.exp(-1j * (K @ np.array(u, float)))[:, None, None] * M[None]
    return U


def vertex(mv, K, q, pert, sign):
    """first-order vertex of move mv for perturbation pert with momentum transfer sign*q (sign = +1: k -> k+q),
    evaluated at input momenta K. Returns None if the move does not couple."""
    if pert == "V":
        if not mv["phi"]:
            return None
        return np.broadcast_to(-1j * I4, (K.shape[0], 4, 4))
    j = pert[1]
    out = None
    for u, M in mv["hops"]:
        if u[j] == 0:
            continue
        uu = np.array(u, float)
        gl = g_line(np.full(1, sign * float(q @ uu)))[0]
        ph = np.exp(-1j * ((K + sign * q[None, :]) @ uu))
        term = (1j * u[j] * gl) * ph[:, None, None] * M[None]
        out = term if out is None else out + term
    return out


def contact(mv, K, q, p1, p2):
    """|lambda|^2-type second-order term of one move (field1 at +q, field2 at -q) at input momentum K"""
    if p1 != p2:
        return None
    if p1 == "V":
        if not mv["phi"]:
            return None
        return np.broadcast_to(-I4, (K.shape[0], 4, 4))
    j = p1[1]
    out = None
    for u, M in mv["hops"]:
        if u[j] == 0:
            continue
        uu = np.array(u, float)
        gl = g_line(np.full(1, float(q @ uu)))[0]
        ph = np.exp(-1j * (K @ uu))
        term = (-(u[j] ** 2) * abs(gl) ** 2) * ph[:, None, None] * M[None]
        out = term if out is None else out + term
    return out


def mm(A, B):
    return np.matmul(A, B)


def compose(moves, K, q, p1, p2=None, need_W2=False):
    """W0(k), W0(k+q), W1p(k) [field1, k -> k+q], W1m(k+q) [field2, k+q -> k], and W2(k) [coefficient of
    lambda1 lambda2^* of W(k)] if need_W2."""
    N = K.shape[0]
    p2 = p1 if p2 is None else p2
    Kq, Kmq = K + q[None], K - q[None]
    Ua = [ev_move(m, K) for m in moves]
    Ub = [ev_move(m, Kq) for m in moves]
    nm = len(moves)
    Id = np.broadcast_to(I4, (N, 4, 4))
    # W1p(k) = sum_s Ub_{n-1..s+1} V1+_s(k) Ua_{s-1..0}
    W0a = Id.copy()
    W1p = np.zeros((N, 4, 4), complex)
    for s in range(nm):
        v = vertex(moves[s], K, q, p1, +1)
        W1p = mm(Ub[s], W1p)
        if v is not None:
            W1p = W1p + mm(v, W0a)
        W0a = mm(Ua[s], W0a)
    W0b = Id.copy()
    W1m = np.zeros((N, 4, 4), complex)
    for s in range(nm):
        v = vertex(moves[s], Kq, q, p2, -1)
        W1m = mm(Ua[s], W1m)
        if v is not None:
            W1m = W1m + mm(v, W0b)
        W0b = mm(Ub[s], W0b)
    if not need_W2:
        return W0a, W0b, W1p, W1m, None
    Uc = [ev_move(m, Kmq) for m in moves]
    pre = Id.copy()              # Ua_{s-1..0}
    X1 = np.zeros((N, 4, 4), complex)    # paths k ->(field1,+q) k+q, propagated at k+q
    X2 = np.zeros((N, 4, 4), complex)    # paths k ->(field2,-q) k-q, propagated at k-q
    W2 = np.zeros((N, 4, 4), complex)
    for s in range(nm):
        v1p = vertex(moves[s], K, q, p1, +1)        # k -> k+q
        v2m = vertex(moves[s], K, q, p2, -1)        # k -> k-q
        v2m_b = vertex(moves[s], Kq, q, p2, -1)     # k+q -> k
        v1p_c = vertex(moves[s], Kmq, q, p1, +1)    # k-q -> k
        c = contact(moves[s], K, q, p1, p2)
        new = mm(Ua[s], W2)
        if v2m_b is not None:
            new = new + mm(v2m_b, X1)
        if v1p_c is not None:
            new = new + mm(v1p_c, X2)
        if c is not None:
            new = new + mm(c, pre)
        W2 = new
        X1 = mm(Ub[s], X1) + (mm(v1p, pre) if v1p is not None else 0)
        X2 = mm(Uc[s], X2) + (mm(v2m, pre) if v2m is not None else 0)
        pre = mm(Ua[s], pre)
    return W0a, W0b, W1p, W1m, W2


# ------------------------------------------------------------------ eigen-decomposition (copied from vuoto_core)
def ueig(W):
    """Cayley transform: K = i(W+1)(W-1)^{-1} Hermitian (spectrum must avoid lambda = 1). Omega in (-pi, pi)."""
    n = W.shape[-1]
    Id = np.broadcast_to(np.eye(n, dtype=complex), W.shape)
    A = np.swapaxes(W - Id, -1, -2)
    Bm = np.swapaxes(W + Id, -1, -2)
    K = 1j * np.swapaxes(np.linalg.solve(A, Bm), -1, -2)
    K = 0.5 * (K + np.conj(np.swapaxes(K, -1, -2)))
    t, V = np.linalg.eigh(K)
    Om = 2 * np.arctan2(1.0, -t)
    Om = np.where(Om > PI, Om - 2 * PI, Om)
    return Om, np.exp(-1j * Om), V


# ------------------------------------------------------------------ responses
def strobo_F(moves, K, q, weight=None, pert="V"):
    """Debiti's stroboscopic density response integrand (chi = -dn/V), generalised: dW = W1p (pert 'V': -i B).
    pert = ('A', j) gives the density response to A_j (mixed kernel), same sign convention dn = -F * A."""
    W0a, W0b, W1p, _, _ = compose(moves, K, q, pert, pert, False)
    Om0, l0, V0 = ueig(W0a)
    Om1, l1, V1 = ueig(W0b)
    Bmn = np.einsum("nam,nab,nbj->nmj", V1.conj(), 1j * W1p, V0)        # B = i W1p  (W1p = -i B for 'V')
    Snm = np.einsum("naj,nam->njm", V0.conj(), V1)
    occ0, occ1 = Om0 < 0, Om1 < 0
    num = Bmn * np.swapaxes(Snm, -1, -2)
    den = l0[:, None, :] - l1[:, :, None]
    A1 = occ0[:, None, :] & (~occ1[:, :, None])
    A2 = (~occ0[:, None, :]) & occ1[:, :, None]
    fac = np.where(A1, 1.0, 0.0) - np.where(A2, 1.0, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(fac != 0, fac * num / den, 0.0)
    F = 1j * T.sum(axis=(1, 2))
    if weight is not None:
        F = F * np.sqrt(weight(K) * weight(K + q[None]))
    return F


def energy_F(moves, K, q, p1, p2=None):
    """integrand of S(q) = coefficient of lambda1 lambda2^* in the sea quasi-energy sum per period, for
    field1 = lambda1 e^{iqx} (+c.c.) and field2 = lambda2 e^{iqx} (+c.c.).
    chi_E = -S_VV / T,   K_jk = S_{Aj,Ak} / T,   mixed = S_{V,Aj} / T."""
    W0a, W0b, W1p, W1m, W2 = compose(moves, K, q, p1, p2, True)
    Om0, l0, V0 = ueig(W0a)
    Om1, l1, V1 = ueig(W0b)
    occ0, occ1 = Om0 < 0, Om1 < 0
    fp0 = np.where(occ0, 1j / l0, 0.0)
    fp1 = np.where(occ1, 1j / l1, 0.0)
    # pairs n at k+q, m at k:  <n|W1p|m> <m|W1m|n> D(n, m)
    A = np.einsum("nai,nab,nbm->nim", V1.conj(), W1p, V0)                                     # [i (k+q), m (k)]
    Bm = np.einsum("nam,nab,nbi->nmi", V0.conj(), W1m, V1)              # [m (k), i (k+q)]
    li = l1[:, :, None]
    lm = l0[:, None, :]
    oi = occ1[:, :, None]
    om = occ0[:, None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        Dboth = -1j / (li * lm)
        Dio = (1j / li) / (li - lm)          # i occupied, m empty
        Dmo = -(1j / lm) / (li - lm)         # m occupied, i empty
    D = np.where(oi & om, Dboth, np.where(oi & ~om, Dio, np.where(~oi & om, Dmo, 0.0)))
    pair = np.sum(A * np.swapaxes(Bm, -1, -2) * D, axis=(1, 2))
    W2d = np.einsum("nam,nab,nbm->nm", V0.conj(), W2, V0)
    diag = np.sum(fp0 * W2d, axis=1)
    return pair + diag


# ------------------------------------------------------------------ quadratures
def sinh_rule(n, kc, beta, centre=0.0):
    t, wt = roots_legendre(n)
    k = kc * np.sinh(beta * t) / np.sinh(beta)
    dk = kc * beta * np.cosh(beta * t) / np.sinh(beta)
    return centre + k, wt * dk


def special_rule(nseg, beta=6.0):
    """Debiti's rule on [-pi/2, pi/2): 4 segments concentrated at -pi/2, 0, 0, pi/2 (copied from vuoto_core)"""
    s, ws = roots_legendre(nseg)
    s = 0.5 * (s + 1)
    ws = 0.5 * ws
    g = np.sinh(beta * s) / np.sinh(beta)
    dg = beta * np.cosh(beta * s) / np.sinh(beta)
    L = PI / 4
    nodes, wts = [], []
    for end, sgn in ((-PI / 2, +1), (0.0, -1), (0.0, +1), (PI / 2, -1)):
        nodes.append(end + sgn * L * g)
        wts.append(L * dg * ws)
    return np.concatenate(nodes), np.concatenate(wts)


def grid3(k1, w1, q):
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    WW = (w1[:, None, None] * w1[None, :, None] * w1[None, None, :]).ravel()
    Ks = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], 1) - 0.5 * q[None, :]
    return Ks, WW


def integrate(fun, Ks, WW, chunk=20000):
    tot = 0.0 + 0.0j
    for s in range(0, len(Ks), chunk):
        tot += np.sum(fun(Ks[s:s + chunk]) * WW[s:s + chunk])
    return tot / (2 * PI) ** 3


def quad_cell(kind, n, beta, q):
    """'cellG': one pi-cell [-pi/2, pi/2]^3, sinh-Gauss concentrated at Gamma (FW1 full sea: Gamma is the only light
    point); 'cellS': Debiti's special rule (frozen full sea: species and cones at 0, +-pi/2);
    'ball': Debiti's projector box [-kc, kc]^3 with kc = 1.05 + |q| + 0.02."""
    if kind == "cellG":
        k1, w1 = sinh_rule(n, PI / 2, beta)
    elif kind == "cellS":
        k1, w1 = special_rule(n // 4, beta)
    else:
        k1, w1 = sinh_rule(n, 1.05 + np.linalg.norm(q) + 0.02, beta)
    return grid3(k1, w1, q)


# ------------------------------------------------------------------ species weight (copied from vuoto_core, K1 0.6, K2 1.05)
def smooth_step(x):
    x = np.clip(x, 0.0, 1.0)
    with np.errstate(divide="ignore", over="ignore"):
        a = np.where(x < 1, np.exp(-1.0 / np.maximum(1 - x, 1e-300)), 0.0)
        b = np.where(x > 0, np.exp(-1.0 / np.maximum(x, 1e-300)), 0.0)
    return np.where(x <= 0, 1.0, np.where(x >= 1, 0.0, a / (a + b)))


def species_weight(K1=0.6, K2=1.05):
    def w(Ks):
        return smooth_step((np.linalg.norm(Ks, axis=-1) - K1) / (K2 - K1))
    return w


# ------------------------------------------------------------------ continuum reference (EXAM ONLY, copied)
def pitR_uehling(q, m=1.0):
    return quad(lambda x: x * (1 - x) * np.log1p(x * (1 - x) * q * q / (m * m)), 0, 1,
                epsabs=1e-16, epsrel=1e-13)[0] / (2 * PI**2)


# ------------------------------------------------------------------ selftests
def _supercell_E(moves, Kt, M, lam, pert, qdir=0):
    """exact sea quasi-energy sum on a supercell of M sites along axis qdir, transverse/Bloch momentum Kt (3,),
    with the real field 2 Re(lam e^{iqx}), q = 2 pi/M along qdir. Brute force (for the selftest)."""
    q = 2 * PI / M
    dim = 4 * M
    W = np.eye(dim, dtype=complex)
    for mv in moves:
        U = np.zeros((dim, dim), complex)
        for u, Mu in mv["hops"]:
            for x in range(M):
                xs = (x - u[qdir]) % M                    # source site
                ph = np.exp(-1j * float(np.dot(Kt, u)))
                fac = 1.0 + 0j
                if pert == "V" and mv["phi"]:
                    Vx = 2 * np.real(lam * np.exp(1j * q * x))
                    fac = np.exp(-1j * Vx)
                elif pert != "V" and u[pert[1]] != 0:
                    j = pert[1]
                    uu = np.array(u, float)
                    qv = np.zeros(3)
                    qv[qdir] = q
                    xsrc = x - u[qdir]
                    pos = np.zeros(3)
                    pos[qdir] = xsrc
                    Phi = lam * u[j] * g_line(np.array([qv @ uu]))[0] * np.exp(1j * (qv @ pos))
                    fac = np.exp(1j * 2 * np.real(Phi))
                U[4 * x:4 * x + 4, 4 * xs:4 * xs + 4] += fac * ph * Mu
        W = U @ W
    Om, _, _ = ueig(W[None])
    return np.sum(Om[0][Om[0] < 0])


def _selftest():
    import time
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    t0 = time.time()
    rng = np.random.default_rng(3)
    K = rng.uniform(-PI, PI, (40, 3))
    import zz_core as zz
    th = 0.3
    W = np.broadcast_to(I4, (40, 4, 4)).copy()
    for m in frozen_moves(th):
        W = ev_move(m, K) @ W
    print("frozen moves vs zz_core.bloch_batch:", np.max(np.abs(W - zz.bloch_batch(K, th))))
    import senzadoppi_r5 as r5
    W = np.broadcast_to(I4, (40, 4, 4)).copy()
    for m in fw1_moves(th):
        W = ev_move(m, K) @ W
    print("FW1 moves vs senzadoppi_r5.FW1.bloch:", np.max(np.abs(W - r5.FW1(0.4).bloch(K, th))))
    # energy second order vs brute-force supercell finite differences
    M = 4
    q = np.array([2 * PI / M, 0, 0])
    for name in ("F0", "FW1"):
        mv = RULES[name](th)
        for pert in ("V", ("A", 0), ("A", 1)):
            Kt = np.array([0.37, 0.21, -0.43])
            ks = np.array([Kt + n * q for n in range(M)])
            S = np.sum(energy_F(mv, ks, q, pert))
            lam = 1e-3
            E = [_supercell_E(mv, Kt, M, s * lam, pert) for s in (-1, 0, 1)]
            fd = (E[0] - 2 * E[1] + E[2]) / (2 * lam**2)
            print(f"{name} {pert}: energy formula {S.real:+.10f} (Im {S.imag:+.1e})  supercell FD {fd:+.10f}")
    print("time %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    _selftest()

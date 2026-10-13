"""vuoto_f2_core.py - phase 2 of 'Vuoto che polarizza' (docs/VUOTO_POLARIZZA.md, section 13).

Three engines, all built on the FROZEN rule of lab/zz_core.py (nothing in the rule is changed here):

1. SELF-ENERGY LOOP IN k (mass shift dm).  Electron = two-tick Floquet bands of the 3-D zig-zag (vuoto_core.U3d),
   light = two-tick Floquet bands of the light walk (zz_core.bloch_batch(spin=1)).  Hole-theory / Feynman form,
   electron at rest (p = 0), intermediate electron momentum k, photon momentum -k:
     Coulomb (instantaneous, kernel 4 pi alpha/k^2 = DEBT, an input):
        dE_C(s) = 1/2 int d3k/(2pi)^3 w(k) 4 pi alpha/k^2  sum_n sgn(E_n) |<n,k|s,0>|^2
     transverse (photons of the light walk, polarization projector P_T from its eigenvectors):
        dE_T(s) = int d3k/(2pi)^3 w(k) (2 pi alpha/omega_k) sum_n  J^dag P_T J
                    x [ theta(E_n>0)/(E_s - E_n - omega) + theta(E_n<0)/(E_s - E_n + omega) ]
        J_j = <n,k| dH_eff/dk_j (at the link midpoint k/2) |s,0>,   H_eff = i log W / 2 (per tick)
     dm = [dE(+) - dE(-)]/2  (cancels any state-independent lattice constant; the definition is declared).
   w(k) = the SAME physical-species weight of phase 1 (K1 = 0.6, K2 = 1.05): it is also the UV cutoff (declared).
   The same formulas with the continuum Dirac H = alpha.k + beta m, J = alpha, omega = |k|, P_T = 1 - kk, are the
   CHECK of the formulas (they must give 3 alpha/2 pi ln(1/m) in total).

2. BETHE LOGARITHM from the resolvent of the RADIAL rule (zz_core.RadialWalk, physical units, a = theta lambda_C).
   u = P_+ (r g_n, 0) put into the p channel (kappa' = +1 or -2; the upper component carries l = 1); the velocity
   source is uv = (2/a) T u (= dE <m|r|n> on each state: dipole-velocity form from the dipole-length one);
   P_+ = projector on quasi-energies in (0, pi) = (1 - sign K)/2, K = i(W+1)(W-1)^{-1} (Hermitian, Cayley),
   sign K = (1/pi) int_{-inf}^{inf} ds e^s [ (K + i e^s)^{-1} + (K - i e^s)^{-1} ]  (trapezoid in s, sparse LU).
   T = i(W/lam_n - 1)(W/lam_n + 1)^{-1}: eigenvalues tan(dOmega/2) (= a dE/2 for small dE);
   S(k) = <uv|T (T + a k/2)^{-1}|uv>,  D = (2/a)<uv|T|uv> = sum |<m|v|n>|^2 dE = sum dE^3 |<m|r|n>|^2
   ln k0(K) = ln K - int_0^K S dk / (K S(K))  -> ln k0 for K -> inf (extrapolation in K^{-1/2}, K^{-1}); Rydberg units.
   K S(K) = sum w dE K/(dE + K) is a smooth-cutoff D: states far above K (relativistic / band-top lattice states, whose
   tan-mapped dE^3 weight would dominate the full D) drop out; their effect vanishes like K/dE.

3. BREATHING COIN (exploration, NOT the rule): flip angle theta - mu on even ticks, theta + mu on odd ticks
   (vuoto_core.floquet_seq(mu=...)); here only the band analysis (masses of the species, of the cones, dispersion).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.special import roots_legendre
from zz_core import ALPHA, bloch_batch, RadialWalk, dirac_formula
from vuoto_core import U3d, ueig, species_weight, PI

ALPHA_INV_2PI = ALPHA / (2 * PI)


# ================================================================== 1. self-energy loop
def sphere_rule(nr, nct, nph, kmin, kmax):
    """product rule in spherical coordinates: k log-spaced Gauss in [kmin, kmax], Gauss in cos, trapezoid in phi.
    returns K (N,3), weights for int d3k (N,)"""
    s, ws = roots_legendre(nr)
    lk = np.log(kmin) + 0.5 * (s + 1) * (np.log(kmax) - np.log(kmin))
    k = np.exp(lk)
    wk = 0.5 * ws * (np.log(kmax) - np.log(kmin)) * k * k * k        # k^2 dk = k^3 d(ln k)
    c, wc = roots_legendre(nct)
    ph = (np.arange(nph) + 0.5) * 2 * PI / nph
    wph = np.full(nph, 2 * PI / nph)
    KK, CC, PP = np.meshgrid(k, c, ph, indexing="ij")
    WW = (wk[:, None, None] * wc[None, :, None] * wph[None, None, :]).ravel()
    sn = np.sqrt(1 - CC**2)
    K = np.stack([(KK * sn * np.cos(PP)).ravel(), (KK * sn * np.sin(PP)).ravel(), (KK * CC).ravel()], 1)
    return K, WW


def lattice_electron(K, th):
    """two-tick Floquet eigen of the electron walk at K (N,3): E per tick (N,4), V (N,4,4)"""
    W = U3d(K, th, 1) @ U3d(K, th, 0)
    Om, lam, V = ueig(W)
    return Om / 2, V


def lattice_Heff(K, th):
    E, V = lattice_electron(K, th)
    return np.einsum("nab,nb,ncb->nac", V, E, V.conj())


def lattice_J(K, th, h=1e-4):
    """dH_eff/dk_j at K, central differences -> (N,3,4,4)"""
    out = []
    for j in range(3):
        d = np.zeros(3); d[j] = h
        out.append((lattice_Heff(K + d, th) - lattice_Heff(K - d, th)) / (2 * h))
    return np.stack(out, 1)


def light_modes(Q):
    """two-tick light walk at Q (N,3): transverse frequencies (N,2) per tick and polarization projector P_T (N,3,3)
    (symmetrized with its complex conjugate, declared)"""
    W = bloch_batch(Q, 0.0, spin=1)
    lam, V = np.linalg.eig(W)
    Om = -np.angle(lam) / 2
    o = np.argsort(-np.abs(Om), axis=1)[:, :2]
    om = np.take_along_axis(np.abs(Om), o, 1)
    Vt = np.take_along_axis(V, o[:, None, :], 2)                  # (N,3,2)
    # orthonormalize the two transverse vectors (normal matrix: already orthogonal up to round-off)
    Qm, _ = np.linalg.qr(Vt)
    P = np.einsum("nia,nja->nij", Qm, Qm.conj())
    P = 0.5 * (P + P.conj())
    return om, P


def rest_states(th, lattice=True):
    """positive / negative energy states at p = 0: (E+, V+ (4,2)), (E-, V- (4,2))"""
    if lattice:
        E, V = lattice_electron(np.zeros((1, 3)), th)
        E, V = E[0], V[0]
    else:
        E, V = np.linalg.eigh(dirac_H(np.zeros((1, 3)), th)[0])
    pos, neg = E > 0, E < 0
    return (E[pos], V[:, pos]), (E[neg], V[:, neg])


# ---- continuum Dirac (CHECK of the formulas only)
SX = np.array([[0, 1], [1, 0]], complex); SY = np.array([[0, -1j], [1j, 0]], complex); SZ = np.diag([1.0 + 0j, -1])
ALPH = [np.block([[0 * SX, s], [s, 0 * SX]]) for s in (SX, SY, SZ)]
BETA = np.diag([1.0 + 0j, 1, -1, -1])


def dirac_H(K, m):
    return np.einsum("nj,jab->nab", K, np.array(ALPH)) + m * BETA[None]


def self_energy(th, lattice=True, nr=72, nct=20, nph=40, kmin_fac=1e-4, chunk=20000, K12=(0.6, 1.05)):
    """returns dict: dE_C(+), dE_C(-), dE_T(+), dE_T(-), dm_C/m, dm_T/m (lattice units; m = th)"""
    w = species_weight(*K12)
    K, WW = sphere_rule(nr, nct, nph, kmin_fac * th, K12[1])
    WW = WW * w(K) / (2 * PI) ** 3
    (Ep, Vp), (En, Vn) = rest_states(th, lattice)
    acc = {"C+": 0.0, "C-": 0.0, "T+": 0.0, "T-": 0.0}
    for s0 in range(0, len(K), chunk):
        Kc, Wc = K[s0:s0 + chunk], WW[s0:s0 + chunk]
        kk = np.linalg.norm(Kc, axis=1)
        if lattice:
            E, V = lattice_electron(Kc, th)
            J = lattice_J(0.5 * Kc, th)                              # midpoint between k and 0
            om, PT = light_modes(-Kc)
        else:
            E, V = np.linalg.eigh(dirac_H(Kc, th))
            J = np.broadcast_to(np.array(ALPH)[None], (len(Kc), 3, 4, 4))
            om = np.stack([kk, kk], 1)
            kh = Kc / kk[:, None]
            PT = np.eye(3)[None] - kh[:, :, None] * kh[:, None, :]
        omg = om.mean(1)
        for lab, Es, Vs in (("+", Ep, Vp), ("-", En, Vn)):
            for i in range(Vs.shape[1]):
                s = Vs[:, i]
                Es_ = Es[i]
                ov = np.einsum("nab,a->nb", V.conj(), s)                 # <n,k|s>
                acc["C" + lab] += 0.5 * np.sum(Wc * 4 * PI * ALPHA / kk**2 * np.sum(np.sign(E) * np.abs(ov)**2, 1)) / Vs.shape[1]
                Js = np.einsum("njab,b->nja", J, s)                    # J_j |s>   (N,3,4)
                M = np.einsum("nab,nja->nbj", V.conj(), Js)             # <n,k|J_j|s>  (N,4,3)
                JPJ = np.real(np.einsum("nbi,nij,nbj->nb", M.conj(), PT, M))
                den = np.where(E > 0, Es_ - E - omg[:, None], Es_ - E + omg[:, None])
                acc["T" + lab] += np.sum(Wc * 2 * PI * ALPHA / omg * np.sum(JPJ / den, 1)) / Vs.shape[1]
    res = dict(acc)
    res["dmC_over_m"] = 0.5 * (acc["C+"] - acc["C-"]) / th
    res["dmT_over_m"] = 0.5 * (acc["T+"] - acc["T-"]) / th
    res["dm_over_m"] = res["dmC_over_m"] + res["dmT_over_m"]
    return res


# ================================================================== 2. Bethe logarithm from the radial rule
class RadialResolvent:
    """sparse tools on the one-tick operator W of a RadialWalk (physical units: a = theta)"""

    def __init__(self, rw):
        import scipy.sparse as sp
        self.rw = rw
        self.W = rw.matrix().tocsc()
        self.I = sp.identity(self.W.shape[0], format="csc", dtype=complex)

    def lu(self, A):
        import scipy.sparse.linalg as spl
        return spl.splu(A.tocsc())

    def project_positive(self, x, h=0.5, smin=None, smax=None, tail=19.0):
        # (amendment A-L1: h = 0.35, tail = 28 -> quadrature error ~1e-12 instead of ~6e-9)
        """P_+ x = (x - sign(K) x)/2, sign(K) = (1/pi) int ds e^s [(K + i t)^{-1} + (K - i t)^{-1}], t = e^s.
        (K + i t)^{-1} = (W - 1) [ i(W + 1) + i t (W - 1) ]^{-1} = -i (W - 1)[(1 + t) W + (1 - t)]^{-1};
        (K - i t)^{-1} = [(K + i t)^{-1}]^dagger (K Hermitian) -> adjoint solve with the same LU."""
        th = self.rw.a * self.rw.m
        smin = np.log(th / 2.0) - tail if smin is None else smin
        smax = np.log(2.0 / th) + tail if smax is None else smax
        ss = np.arange(smin, smax + h / 2, h)
        W, I = self.W, self.I
        Wm1 = W - I
        acc = np.zeros_like(x)
        for s in ss:
            t = np.exp(s)
            A = (1 + t) * W + (1 - t) * I
            lu = self.lu(A)
            y1 = -1j * (Wm1 @ lu.solve(x))                         # (K + it)^{-1} x
            # (K - it)^{-1} x = [(K+it)^{-1}]^dag x = (A^{-1})^dag (W-1)^dag (+i) x
            y2 = 1j * lu.solve(Wm1.getH() @ x, trans="H")
            acc += t * (y1 + y2)
        sgn = acc * h / PI
        return 0.5 * (x - sgn), len(ss)

    def band_filter(self, x, omega_c, h=0.5, smin=-4.0, smax=4.0):
        """keep the states with quasi-energy Omega < omega_c (mod 2 pi, on the positive side): (1 + sign K_c) x / 2,
        K_c = Cayley of e^{i omega_c} W (eigenvalues -cot((Omega - omega_c)/2): > 0 for Omega in (omega_c - pi, omega_c)).
        The finite s range makes it a SMOOTH step of width ~ e^{smin} rad around omega_c (declared filter)."""
        ss = np.arange(smin, smax + h / 2, h)
        W = np.exp(1j * omega_c) * self.W
        I = self.I
        Wm1 = W - I
        acc = np.zeros_like(x)
        for s in ss:
            t = np.exp(s)
            lu = self.lu((1 + t) * W + (1 - t) * I)
            acc += t * (-1j * (Wm1 @ lu.solve(x)) + 1j * lu.solve(Wm1.getH() @ x, trans="H"))
        return 0.5 * (x + acc * h / PI), len(ss)


def bethe_log(n, l_kappa, th, Z=1, rmax_bohr=50.0, Ks=None, h=0.5, nk=240, verbose=False, return_raw=False, species=True, omega_cut=0.6, tail=19.0, deflate=0):
    """ln k0 for the s state n (kappa = -1) with intermediate channel kappa' (= +1 or -2). Physical units (m = 1).
    Returns dict with K list, lnk0(K), D, and the reference sum-rule value."""
    import scipy.sparse.linalg as spl
    Za = Z * ALPHA
    a_B = 1.0 / Za
    rmax = rmax_bohr * a_B * n
    rw_s = RadialWalk(-1, th, rmax, Zalpha=Za)
    E0 = dirac_formula(n, -1, Za)
    om, u_s, res = rw_s.eigen(E0, iters=12)
    g, f = rw_s.gf(u_s)
    # dipole source in the p channel: upper component r g (l = 1), lower 0;  put into the chi basis of the walk
    rw_p = RadialWalk(l_kappa, th, rmax, Zalpha=Za)
    gp = rw_p.r * np.real(g)
    if species:
        # the SAME physical-species projector as the vacuum (w(|k a|), K1 = 0.6, K2 = 1.05), DST-II along the radius:
        # without it the source has a tiny content at k a ~ pi (doubler, band top) whose dE^3 weight is huge
        from scipy.fft import dst, idst
        kk = PI * (np.arange(len(gp)) + 1) / len(gp)
        gp = idst(dst(gp, type=2) * np.sqrt(species_weight()(kk)), type=2)
    fp = np.zeros_like(gp)
    # (g, f) -> chi: g = (c+ + c-)/sqrt2, f = i (c+ - c-)/sqrt2  ->  c+ = (g - i f)/sqrt2, c- = (g + i f)/sqrt2
    cp = (gp - 1j * fp) / np.sqrt(2); cm = (gp + 1j * fp) / np.sqrt(2)
    # the walk's eigenvectors live BEFORE the half coin (chi = Mh u):  v = Mh u
    Mh = rw_p.Mh
    v = np.empty(2 * rw_p.N, complex)
    v[0::2] = Mh[:, 0, 0] * cp + Mh[:, 0, 1] * cm
    v[1::2] = Mh[:, 1, 0] * cp + Mh[:, 1, 1] * cm
    v *= np.sqrt(th)                                               # normalization of sum |.|^2 a
    RR = RadialResolvent(rw_p)
    up, nlu = RR.project_positive(v, h=h, tail=tail)
    if omega_cut is not None:
        # keep only the physical band (quasi-energy below Omega_n + omega_cut): the band-top states of the 1-D walk
        # (doubler side, k a ~ pi) would otherwise dominate the dE^3-weighted sums
        up, n2 = RR.band_filter(up, om * th + omega_cut)
        nlu += n2
    lam_n = np.exp(-1j * om * th)
    W, I = RR.W, RR.I
    if deflate:
        # amendment A-L2: the walk's spectrum is mirrored (doubler partner at phase +pi); the mirror of the state
        # degenerate with |n> (2P_1/2 for 2S) sits at W/lam = -1, the pole of the Cayley T: remove those few
        # eigenvectors from the source explicitly (W is unitary: orthogonal projection)
        _, Vd = spl.eigs((W / lam_n).tocsc(), k=deflate, sigma=-1.0, which="LM")
        Vd, _ = np.linalg.qr(Vd)
        up = up - Vd @ (Vd.conj().T @ up)
    # T = i (W/lam - 1)(W/lam + 1)^{-1};  T (T + kap)^{-1} = i(W/lam - 1) [ i (W/lam - 1) + kap (W/lam + 1) ]^{-1}
    Wl = W / lam_n
    luD = spl.splu((Wl + I).tocsc())
    # velocity source: <m|v|n> = i dE <m|r|n>  ->  uv = T P+ u  (T = a dE/2 on each state; the factor (2/a)^2 cancels)
    uv = 1j * ((Wl - I) @ luD.solve(up))
    Tu = 1j * ((Wl - I) @ luD.solve(uv))
    D = (2.0 / th) * np.real(np.vdot(uv, Tu)) * (2.0 / th) ** 2            # = sum dE^3 |<m|r|n>|^2
    up = uv * (2.0 / th)
    norm_neg = np.real(np.vdot(v - up, v - up))
    if Ks is None:
        Ks = np.geomspace(1e-5, 0.3, 61)
    # S(k) on a log grid (trapezoid of k S(k) in ln k) plus the piece [0, k_0] with S(0) = <u|u>
    kg = np.unique(np.concatenate([np.geomspace(1e-8, max(Ks), nk), np.asarray(Ks, float)]))
    S = np.zeros(len(kg))
    for i, k in enumerate(kg):
        kap = th * k / 2.0
        lu = spl.splu((1j * (Wl - I) + kap * (Wl + I)).tocsc())
        y = 1j * ((Wl - I) @ lu.solve(up))
        S[i] = np.real(np.vdot(up, y))
    S0 = np.real(np.vdot(up, up))
    lk = np.log(kg)
    cum = kg[0] * S0 + np.concatenate([[0.0], np.cumsum(0.5 * (kg[1:] * S[1:] + kg[:-1] * S[:-1]) * np.diff(lk))])
    Ry = 0.5 * Za**2                                               # Rydberg (times Z^2), units m c^2
    lnk, lnk_D, DK = [], [], []
    for K in Ks:
        i = int(np.argmin(np.abs(kg - K)))
        IK = cum[i]
        DK.append(K * S[i])                                        # smooth-cutoff D: sum w dE K/(dE + K)
        lnk.append(np.log(K / Ry) - IK / DK[-1])                   # USED: the lattice states far above K drop out
        lnk_D.append(np.log(K / Ry) - IK / D)                      # with the full D (sensitive to the band top)
    sumrule = 2 * PI * Za * (Za**3 / (PI * n**3))                  # NR: sum dE^3 |<m|r|n>|^2 = 2 pi Z alpha |psi(0)|^2
    out = {"n": n, "kappa_p": l_kappa, "theta": th, "Z": Z, "rmax_bohr": rmax_bohr, "N": rw_p.N, "omega_n": om,
           "eig_resid": float(res), "D": D, "D_sumrule_NR": sumrule, "neg_weight": norm_neg, "n_lu_proj": nlu,
           "K": list(map(float, Ks)), "lnk0_K": list(map(float, lnk)), "lnk0_K_fullD": list(map(float, lnk_D)),
           "D_K": list(map(float, DK))}
    if return_raw:
        out["kg"] = kg.tolist(); out["S"] = S.tolist()
    return out


def extrapolate_lnk0(Ks, lnk, kmin=None, kmax=None):
    """fit lnk0(K) = L + b K^{-1/2} + c K^{-1} on the window; returns L, b, c"""
    Ks, lnk = np.asarray(Ks), np.asarray(lnk)
    sel = np.ones(len(Ks), bool)
    if kmin is not None:
        sel &= Ks >= kmin
    if kmax is not None:
        sel &= Ks <= kmax
    A = np.stack([np.ones(sel.sum()), Ks[sel] ** -0.5, Ks[sel] ** -1.0], 1)
    c, *_ = np.linalg.lstsq(A, lnk[sel], rcond=None)
    return c


# ================================================================== 3. breathing coin: bands
def breathing_W(K, th, mu):
    """two-tick W with flip angles th - mu (tick 0, order xyz) and th + mu (tick 1, order zyx)"""
    return U3d(K, th + mu, 1) @ U3d(K, th - mu, 0)


def band_gaps(th, mu):
    """|Omega|/2 at the special points: the 4 massive valleys and the 4 cones (one pi-cell)"""
    P = PI / 2
    pts = {"valle (0,0,0)": (0, 0, 0), "valle (pi/2,pi/2,0)": (P, P, 0), "valle (pi/2,0,pi/2)": (P, 0, P),
           "valle (0,pi/2,pi/2)": (0, P, P), "cono (pi/2,0,0)": (P, 0, 0), "cono (0,pi/2,0)": (0, P, 0),
           "cono (0,0,pi/2)": (0, 0, P), "cono (pi/2,pi/2,pi/2)": (P, P, P)}
    out = {}
    for name, k in pts.items():
        W = breathing_W(np.array([k], float), th, mu)
        lam = np.linalg.eigvals(W[0])
        Om = -np.angle(lam)
        if name.startswith("valle"):
            out[name] = float(np.min(np.abs(Om)) / 2)               # gap around 0 (per tick): mass
        else:
            out[name] = float(np.min(np.abs(np.abs(Om) - PI)) / 2)  # gap around pi: cone mass
    return out


# ================================================================== 1b. g - 2 from the same loop (Breit frame)
# Euclidean (Wick-rotated) form of the same loop: Sigma(E) = int d3k/(2pi)^3 sum_n [4 pi alpha J^dag P_T J g(E-E_n, w)
#   - 4 pi alpha/k^2 |<n|s>|^2 h(b_n)],  g = 1/(2w(D - b w)) with the FEYNMAN prescription by band b = sgn(E_n)
#   (not by the sign of D: a moving external electron has D > 0 for slower intermediate states), h(b) = -b/2.
#   (Checked: gives back the old-fashioned dE_C, dE_T above, and the Ward identity Lambda(p,p) = dSigma/dp.)
# Vertex = first-order change of Sigma under H -> H + V_ext: G -> G + G V G, pair functions
#   f(D1, D2) = [g(D1) - g(D2)]/(D2 - D1),  fC = [h(D1) - h(D2)]/(D2 - D1).
# F2(0) = Lambda_M - Lambda_E: Lambda_E = -dSigma/dE at rest (charge vertex, Ward), Lambda_M = transverse current
# vertex in the Breit frame (p = -Q/2, p' = +Q/2, there the current is pure G_M = F1 + F2) divided by its tree value.
# IR: photon mass lam in omega (numerical regulator, F2 is IR-finite; lam -> 0 checked).
def _g(D, w, b):
    """int dy/2pi D_photon x Feynman electron propagator, band b = +1 (E_n > 0) / -1 (E_n < 0): 1/(2w(D - b w))"""
    return 1.0 / (2 * w * (D - b * w))


def _f(D1, D2, w, b1, b2):
    same = b1 == b2
    with np.errstate(divide="ignore", invalid="ignore"):
        f_same = 1.0 / (2 * w * (D1 - b1 * w) * (D2 - b2 * w))
        f_opp = (_g(D1, w, b1) - _g(D2, w, b2)) / (D2 - D1)
    return np.where(same, f_same, f_opp)


def _fC(D1, D2, b1, b2):
    """Coulomb (instantaneous): h(b) = -b/2"""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(b1 == b2, 0.0, -0.5 * (b1 - b2) / (D2 - D1))


def _electron(K, th, lattice):
    if lattice:
        return lattice_electron(K, th)
    return np.linalg.eigh(dirac_H(K, th))


def _J(K, th, lattice):
    if lattice:
        return lattice_J(K, th)
    return np.broadcast_to(np.array(ALPH)[None], (len(K), 3, 4, 4))


def _photon(Kph, lam, lattice):
    kk = np.linalg.norm(Kph, axis=1)
    if lattice:
        om, PT = light_modes(Kph)
        om = om.mean(1)
    else:
        om = kk.copy()
        kh = Kph / kk[:, None]
        PT = np.eye(3)[None] - kh[:, :, None] * kh[:, None, :]
    return np.sqrt(om**2 + lam**2), PT, kk


def sigma_matrix(E, P, th, lattice, lam, K, WW):
    """4x4 self-energy Sigma(E) at momentum P in the spinor basis (Feynman form, same loop)"""
    om, PT, kk = _photon(K, lam, lattice)
    VC = 4 * PI * ALPHA / kk**2
    En, V = _electron(P[None] - K, th, lattice)
    J = _J(P[None] - 0.5 * K, th, lattice)
    bn = np.sign(En)
    G = _g(E - En, om[:, None], bn)                                  # (N, n)
    # sum_n J_i^dag |n> g <n| J_j P_ij  ->  operator
    JV = np.einsum("nibc,ncm->nibm", J, V)                          # J_i |n>  (N,3,4,n)
    T = np.einsum("nibm,nij,nm,njdm->nbd", JV, PT, G, JV.conj())      # sum_n J_i|n> g <n|J_j  (J Hermitian)
    C = np.einsum("nbm,nm,ndm->nbd", V, -bn / 2, V.conj())
    return np.einsum("n,nbd->bd", WW * 4 * PI * ALPHA, T) - np.einsum("n,nbd->bd", WW * VC, C)


def pgradH(P, th, lattice, h=1e-6):
    """(p . grad_p H)(P): the kinetic operator the velocity counterterm multiplies (continuum: alpha . p)"""
    if not lattice:
        return np.einsum("j,jab->ab", P, np.array(ALPH))
    out = np.zeros((4, 4), complex)
    for j in range(3):
        d = np.zeros(3); d[j] = h
        out += P[j] * (lattice_Heff((P + d)[None], th)[0] - lattice_Heff((P - d)[None], th)[0]) / (2 * h)
    return out


def dH_dtheta(P, th, lattice, h=1e-6):
    if not lattice:
        return BETA.copy()
    return (lattice_Heff(P[None], th + h)[0] - lattice_Heff(P[None], th - h)[0]) / (2 * h)


def g2_all(th, lattice=True, Q=None, lam=None, nr=72, nct=20, nph=40, kmin_fac=1e-3, chunk=8000, K12=(0.6, 1.05),
           pdisp=0.05, pairs=((2, 0), (2, 1), (0, 1), (0, 2), (1, 0), (1, 2))):
    """F2 for each (Q axis, current axis) pair. Shared: Lambda_E = -dSigma/dE at rest, mass counterterm dm and the
    ANISOTROPIC velocity counterterms dv_j (measured on the dispersion along axis j at p = pdisp m).
    Lattice units, m = th. Returns dict with per-pair results and the average."""
    Q = 1e-5 * th if Q is None else Q
    lam = 1e-3 * th if lam is None else lam
    w = species_weight(*K12)
    K, WW = sphere_rule(nr, nct, nph, kmin_fac * lam, K12[1])
    WW = WW * w(K) / (2 * PI) ** 3

    def pos_states(P):
        E, V = _electron(P[None], th, lattice)
        sel = E[0] > 0
        return E[0][sel], V[0][:, sel]
    E0, U0 = pos_states(np.zeros(3))
    # ---- Lambda_E = -dSigma/dE at rest (Ward), spin-averaged; Coulomb part is 0 (instantaneous)
    LamE = 0.0
    for s0 in range(0, len(K), chunk):
        Kc, Wc = K[s0:s0 + chunk], WW[s0:s0 + chunk]
        om, PT, kk = _photon(Kc, lam, lattice)
        E, V = _electron(-Kc, th, lattice)
        J = _J(-0.5 * Kc, th, lattice)
        for i in range(U0.shape[1]):
            s = U0[:, i]
            M = np.einsum("nba,njb->naj", V.conj(), np.einsum("njab,b->nja", J, s))     # <n|J_j|s>
            JPJ = np.real(np.einsum("nai,nij,naj->na", M.conj(), PT, M))
            D = E0[i] - E
            LamE += np.sum(Wc * 4 * PI * ALPHA * np.sum(JPJ * _f(D, D, om[:, None], np.sign(E), np.sign(E)), 1)) / U0.shape[1]
    # ---- counterterms: mass (rest energy) and velocity per axis (dressed E^2 - p^2 = m^2 along each axis)
    S0 = sigma_matrix(E0[0], np.zeros(3), th, lattice, lam, K, WW)
    dm_c = (np.real(np.trace(U0.conj().T @ S0 @ U0)) /
            np.real(np.trace(U0.conj().T @ dH_dtheta(np.zeros(3), th, lattice) @ U0)))
    dv = np.zeros(3)
    for ax in range(3):
        Pt = np.zeros(3); Pt[ax] = pdisp * th
        Et, Ut = pos_states(Pt)
        St = sigma_matrix(Et.mean(), Pt, th, lattice, lam, K, WW)
        dEt = np.real(np.trace(Ut.conj().T @ St @ Ut)) / 2
        mass_part = dm_c * np.real(np.trace(Ut.conj().T @ dH_dtheta(Pt, th, lattice) @ Ut)) / 2
        kin = np.real(np.trace(Ut.conj().T @ pgradH(Pt, th, lattice) @ Ut)) / 2
        dv[ax] = (dEt - mass_part) / kin
    res = {"theta": th, "lattice": bool(lattice), "Q": Q, "lam": lam, "LamE": float(LamE), "dm_c_over_m": float(dm_c / th),
           "dv_over_alpha": (dv / ALPHA).tolist(), "pairs": {}}
    for qax, jax in pairs:
        eq = np.zeros(3); eq[qax] = 1.0
        p, pp = -0.5 * Q * eq, 0.5 * Q * eq
        Ep, Up = pos_states(p)
        Epp, Upp = pos_states(pp)
        Eb = 0.5 * (Ep.mean() + Epp.mean())
        Jt = _J(np.zeros((1, 3)), th, lattice)[0, jax]               # tree current at the midpoint (p+p')/2 = 0
        Tree = Upp.conj().T @ Jt @ Up
        LamM = np.zeros((2, 2), complex)
        for s0 in range(0, len(K), chunk):
            Kc, Wc = K[s0:s0 + chunk], WW[s0:s0 + chunk]
            om, PT, kk = _photon(Kc, lam, lattice)
            VC = 4 * PI * ALPHA / kk**2
            E1, V1 = _electron(p[None] - Kc, th, lattice)             # n' at p - k
            E2, V2 = _electron(pp[None] - Kc, th, lattice)            # n  at p' - k
            Jin = _J(p[None] - 0.5 * Kc, th, lattice)
            Jout = _J(pp[None] - 0.5 * Kc, th, lattice)
            Jext = _J(-Kc, th, lattice)[:, jax]                         # midpoint of p-k and p'-k
            Vext = np.einsum("nab,nbc,ncd->nad", V2.conj().transpose(0, 2, 1), Jext, V1)
            Mi = np.einsum("nab,njbc,cs->najs", V1.conj().transpose(0, 2, 1), Jin, Up)
            Mo = np.einsum("nab,njbc,cs->najs", V2.conj().transpose(0, 2, 1), Jout, Upp)
            D1, D2 = Eb - E2, Eb - E1
            b1, b2 = np.sign(E2)[:, :, None], np.sign(E1)[:, None, :]
            F = _f(D1[:, :, None], D2[:, None, :], om[:, None, None], b1, b2)
            FC = _fC(D1[:, :, None], D2[:, None, :], b1, b2)
            A = np.einsum("naiS,nij,nab,nbjs,nab->nSs", Mo.conj(), PT, Vext, Mi, F, optimize=True)
            LamM += np.einsum("n,nSs->Ss", Wc * 4 * PI * ALPHA, A)
            O2 = np.einsum("nab,bS->naS", V2.conj().transpose(0, 2, 1), Upp)
            O1 = np.einsum("nab,bs->nas", V1.conj().transpose(0, 2, 1), Up)
            C = np.einsum("naS,nab,nbs,nab->nSs", O2.conj(), Vext, O1, FC, optimize=True)
            LamM -= np.einsum("n,nSs->Ss", Wc * VC, C)
        # external legs: on-shell Sigma (minus counterterms) mixes the electron with the negative states at p
        leg = np.zeros((2, 2), complex)
        for P, side in ((p, "in"), (pp, "out")):
            Sren = (sigma_matrix(Eb, P, th, lattice, lam, K, WW) - dm_c * dH_dtheta(P, th, lattice)
                    - dv[qax] * pgradH(P, th, lattice))
            Eall, Vall = _electron(P[None], th, lattice)
            Vneg = Vall[0][:, Eall[0] < 0]; Eneg = Eall[0][Eall[0] < 0]
            for a in range(Vneg.shape[1]):
                va = Vneg[:, a]
                if side == "in":
                    leg += np.outer(Upp.conj().T @ Jt @ va, va.conj() @ Sren @ Up) / (Eb - Eneg[a])
                else:
                    leg += np.outer(Upp.conj().T @ Sren @ va, va.conj() @ Jt @ Up) / (Eb - Eneg[a])
        nT = np.trace(Tree.conj().T @ Tree)
        cV = np.trace(Tree.conj().T @ LamM) / nT
        cL = np.trace(Tree.conj().T @ leg) / nT
        tot = LamM + leg - dv[jax] * Tree
        c = np.trace(Tree.conj().T @ tot) / nT
        resid = np.linalg.norm(tot - c * Tree) / np.linalg.norm(tot)
        F2 = c.real - LamE
        res["pairs"]["%d%d" % (qax, jax)] = {"vertex": float(cV.real), "legs": float(cL.real), "ct": float(-dv[jax]),
                                             "imag": float(c.imag), "resid_struct": float(resid),
                                             "F2_over_alpha2pi": float(F2 / ALPHA_INV_2PI)}
    vals = [v["F2_over_alpha2pi"] for v in res["pairs"].values()]
    res["F2_mean_over_alpha2pi"] = float(np.mean(vals))
    res["F2_spread"] = float(np.max(vals) - np.min(vals))
    return res

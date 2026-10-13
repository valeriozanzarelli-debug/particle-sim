# FROZEN COPY of lab/aggfull_core.py at commit a70a372 (study aggmol_: the aggfull agent works on its file in parallel)
"""aggfull_core.py - AGGLOMERATORE COMPLETO: the coarse zig-zag rule with NO hand-written pieces.

Builds on lab/agg_zz_core.py (read-only import; the learned free hierarchy of chain A lives in
lab/out/agg_zz_3d_A.pkl). New here:

 (1) THE CELL-FIELD COUPLING IS LEARNED. The fine 3-D zig-zag walk (lab/zz_core.py, frozen) is run in a weak potential
     V(x) = eps e^{i q.x} (a plane wave of the field: by linearity and translation invariance this is the same information as
     blocks of the walk in any weak field, e.g. the far field of a nucleus). Its first-order response over one coarse step
     (B ticks), seen through the coarse descriptor (the positive-energy band in the Foldy-Wouthuysen-aligned basis: the
     unitary map whose large component is Hermitian positive), is inverted (Duhamel) into the coupling of the coarse
     GENERATOR:  H_c = G/tau + V_c,   <k'|V_c|k> = Vt(k'-k) Gamma(k', k)   (2x2 in spin).
     Gamma is fitted (least squares) with a local "sandwich" stencil  Gamma(k',k) = sum_{j,l} C_{jl} e^{-i k'.j} e^{i k.l},
     C_{jl} 2x2 complex, C_{lj} = C_{jl}^dag (Hermitian), j,l in a stencil chosen AUTOMATICALLY (smallest that passes the
     check). Nothing about Darwin or spin-orbit is put in: they are read off the learned C afterwards.
 (2) The next level's coupling is learned from the LEARNED level below (repeat), as for the free rule.
 (3) Coarse potential field = the band-limited restriction (coarse Fourier modes) of the fine input potential -Z alpha/r:
     the same plane-wave descriptor convention that the free hierarchy uses (no cell average, no corner/centre choice).
 (4) Recoil: the proton is a ball of the same rule, theta_p = m_p a / hbar on its own fine lattice; its coarse generator is
     learned by the same agglomerator; the pair state in the sector of zero total momentum (exact for translation-invariant
     rules) feels G_e(k) + G_p(-k): the reduced mass must come out of the two learned generators.

Units hbar = c = m_e = 1 (length lambda_C, time hbar/mc^2, energy mc^2).
"""
import os, sys, time, pickle, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
import scipy.fft as sfft
from scipy.special import erf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zz_core as zz
import agg_zz_core as azc
azc.WK = 1                                           # one thread (shared container)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
MC2 = zz.MC2_EV
RY = 0.5 * zz.ALPHA ** 2 * MC2
TH_A = A0 / 16 / 256                                 # chain A fine theta (agg_zz_3d.py)
MP = zz.MP_ME                                        # proton mass / electron mass (input: mass of the nucleus)
I2 = np.eye(2, dtype=complex)
PAULI = [zz.SX, zz.SY, zz.SZ]


def fftn(a):
    return sfft.fftn(a, axes=(-3, -2, -1), workers=1)


def ifftn(a):
    return sfft.ifftn(a, axes=(-3, -2, -1), workers=1)


def load_chain_A():
    st = pickle.load(open(os.path.join(OUT, "agg_zz_3d_A.pkl"), "rb"))
    return st["levels"]


def genrule(rec):
    return azc.GenRule(rec["K"], rec["offs"], 3, rec["C"], rec["phi0"])


# =====================================================================================================
#  fine walk: one-tick Bloch matrices (with the theta flip), positive band in the FW-aligned basis
# =====================================================================================================
def tick_mats(KA, theta):
    """one-tick Bloch matrices M_0(k), M_1(k) (even/odd tick: axis order x,y,z / z,y,x), KA (N,3) = k a"""
    KA = np.atleast_2d(KA)
    N = KA.shape[0]

    def Tax(j, h):
        lam, P = zz.BASIS_HALF[j]
        ph = np.exp(-1j * h * lam[None, :] * KA[:, j:j + 1])
        return np.einsum("ab,nb,cb->nac", P, ph, P.conj())
    out = []
    ct, st = np.cos(theta), np.sin(theta)
    for nt in (0, 1):
        TL = np.broadcast_to(I2, (N, 2, 2)).copy()
        TR = TL.copy()
        for j in zz.order_of_tick(nt):
            TL = Tax(j, -1) @ TL
            TR = Tax(j, +1) @ TR
        Wt = np.zeros((N, 4, 4), complex)
        Wt[:, :2, :] = np.concatenate([ct * TL, -1j * st * TR], axis=2)
        Wt[:, 2:, :] = np.concatenate([-1j * st * TL, ct * TR], axis=2)
        out.append(Wt)
    return out


def fw_band(W2):
    """positive band of the two-tick matrix W2 (N,4,4): orthonormal 4x2 basis U aligned so that its large component
    Phi = (U_L + U_R)/sqrt2 is Hermitian positive (Foldy-Wouthuysen gauge); returns U (N,4,2), eigen-phases (N,2)"""
    lam, V = np.linalg.eig(W2)
    w = azc.phi_weight(V)
    o = np.argsort(-w, axis=-1)[..., :2]
    U = np.take_along_axis(V, o[..., None, :], axis=-1)
    U, _ = np.linalg.qr(U)                                             # orthonormal (unitary W2: bands orthogonal)
    Phi = (U[:, 0:2, :] + U[:, 2:4, :]) / np.sqrt(2)
    Ws, s, Vh = np.linalg.svd(Phi)
    Q = np.conj(np.swapaxes(Vh, -1, -2)) @ np.conj(np.swapaxes(Ws, -1, -2))   # Phi Q = Ws S Ws^dag
    U = U @ Q
    return U


def fine_response(k, kp, theta, B):
    """first-order response of B fine ticks to a unit potential plane wave (V(x) = e^{i(kp-k).x}, phase e^{-iVa} per
    tick, applied after the flip as in zz_core.Walk): returns T(k) (N,4,4), T(kp), dT(kp<-k) (N,4,4).
    k, kp in physical units (1/lambda_C)."""
    a = theta
    Mk = tick_mats(k * a, theta)
    Mp = tick_mats(kp * a, theta)
    N = k.shape[0]
    # prefix products P_t = M_{t}..M_0 (k), suffix S_t = M_{B-1}..M_{t+1} (kp)
    pref = []
    P = np.broadcast_to(np.eye(4, dtype=complex), (N, 4, 4)).copy()
    for t in range(B):
        P = Mk[t % 2] @ P
        pref.append(P)
    suf = [None] * B
    S = np.broadcast_to(np.eye(4, dtype=complex), (N, 4, 4)).copy()
    for t in range(B - 1, -1, -1):
        suf[t] = S
        S = S @ Mp[t % 2]
    Tkp = S
    dT = np.zeros((N, 4, 4), complex)
    for t in range(B):
        dT += suf[t] @ pref[t]
    dT *= -1j * a
    return pref[-1], Tkp, dT


def duhamel_invert(D, Gk, Gkp, phi0):
    """given the first-order change D of one coarse step e^{i phi0} e^{-i G} between k and kp (2x2), return X such that
    D = d/de [e^{i phi0} e^{-i(G + e X)}] (first order, Duhamel). Gk, Gkp: (N,2,2) Hermitian generators."""
    ek, Uk = np.linalg.eigh(Gk)
    ep, Up = np.linalg.eigh(Gkp)
    Dt = np.conj(np.swapaxes(Up, -1, -2)) @ D @ Uk * np.exp(-1j * phi0)
    diff = ep[:, :, None] - ek[:, None, :]
    ssum = ep[:, :, None] + ek[:, None, :]
    fac = -1j * np.exp(-0.5j * ssum) * np.sinc(diff / (2 * np.pi))          # np.sinc(x) = sin(pi x)/(pi x)
    Xt = Dt / fac
    return Up @ Xt @ np.conj(np.swapaxes(Uk, -1, -2))


def coarse_from_fine(k, kp, theta, B):
    """level learned directly from the fine walk: generator G(k), G(kp) (2x2, FW basis), rest phase phi0 and the coupling
    X(kp,k) per unit potential amplitude, in units of the coarse step: H_c tau = G + tau V_c -> returns Gamma = X / tau"""
    Tk, Tkp, dT = fine_response(k, kp, theta, B)
    # two-tick matrices for the band basis (B even: T = W2^(B/2))
    W2k = tick_mats(k * theta, theta)
    W2k = W2k[1] @ W2k[0]
    W2p = tick_mats(kp * theta, theta)
    W2p = W2p[1] @ W2p[0]
    Uk, Up = fw_band(W2k), fw_band(W2p)
    ck = np.conj(np.swapaxes(Uk, -1, -2)) @ Tk @ Uk                       # (N,2,2) coarse free step at k
    cp = np.conj(np.swapaxes(Up, -1, -2)) @ Tkp @ Up
    z = np.zeros((1, 3))
    W20 = tick_mats(z, theta); W20 = W20[1] @ W20[0]
    U0 = fw_band(W20)
    T0 = np.linalg.matrix_power(W20, B // 2)
    phi0 = float(np.angle(np.mean(np.linalg.eigvals(np.conj(U0[0].T) @ T0 @ U0[0]))))
    Gk = gen_from_step(ck, phi0)
    Gp = gen_from_step(cp, phi0)
    D = np.conj(np.swapaxes(Up, -1, -2)) @ dT @ Uk
    X = duhamel_invert(D, Gk, Gp, phi0)
    tau = B * theta
    return Gk, Gp, phi0, X / tau, D


def gen_from_step(M, phi0):
    """Hermitian generator G with M = e^{i phi0} e^{-i G} (positive band: motion phase unwrapped into [-0.5, 2pi-0.5))"""
    lam, V = np.linalg.eig(M * np.exp(-1j * phi0))
    lg = np.real(1j * np.log(lam / np.abs(lam)))
    lg = np.where(lg < -0.5, lg + 2 * np.pi, lg)
    G = V @ (lg[..., :, None] * np.linalg.inv(V))
    return 0.5 * (G + np.conj(np.swapaxes(G, -1, -2)))


# =====================================================================================================
#  the learned coupling: sandwich stencil Gamma(k',k) = sum_{j,l} C_jl e^{-ik'.j} e^{ik.l}
# =====================================================================================================
STENCILS = {
    "point": [(0, 0, 0)],
    "star": [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)],
    "cube": [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)],
}


def stencil(name):
    return np.array(STENCILS[name], int)


def fit_coupling(kc, kpc, Gam, offs, wts=None, rcond=1e-12):
    """least squares for C (2,2,n,n) Hermitian: Gamma(kp,k) ~ sum_jl C_jl e^{-i kp.j} e^{i k.l}; kc, kpc in CELL units.
    Hermiticity C_lj = C_jl^dag imposed by symmetrising the data with the swapped pairs."""
    n = len(offs)
    # data symmetrisation: Gamma(k,kp) = Gamma(kp,k)^dag for a Hermitian coupling (true for the exact one)
    K = np.concatenate([kc, kpc]); KP = np.concatenate([kpc, kc])
    Gm = np.concatenate([Gam, np.conj(np.swapaxes(Gam, -1, -2))])
    w = np.ones(len(K)) if wts is None else np.concatenate([wts, wts])
    Phi = np.einsum("nj,nl->njl", np.exp(-1j * KP @ offs.T), np.exp(1j * K @ offs.T)).reshape(len(K), n * n)
    sol = np.linalg.lstsq(Phi * w[:, None], (Gm.reshape(len(K), 4)) * w[:, None], rcond=rcond)[0]   # (n*n, 4)
    C = sol.T.reshape(2, 2, n, n)
    C = 0.5 * (C + np.conj(np.transpose(C, (1, 0, 3, 2))))
    return C


def eval_coupling(C, offs, kc, kpc):
    ph = np.einsum("nj,nl->njl", np.exp(-1j * kpc @ offs.T), np.exp(1j * kc @ offs.T))
    return np.einsum("abjl,njl->nab", C, ph)


def physics_of_coupling(C, offs, h):
    """read Darwin and spin-orbit off a learned coupling (cell h in lambda_C). Expansion at small k, k' (physical units):
    Gamma = g0 + gk (k^2 + k'^2) + gx k'.k + i gso sigma.(k' x k) + ...   Continuum Dirac (FW): g0 = 1, gk = -1/8,
    gx = +1/4 (together: Darwin -q^2/8), gso = 1/4 (spin-orbit (1/4m^2) sigma.(grad V x p))."""
    eps = 1e-3 / h
    rng = np.random.default_rng(1)
    k = rng.normal(size=(400, 3)) * eps
    kp = rng.normal(size=(400, 3)) * eps
    G = eval_coupling(C, offs, k * h, kp * h)
    s = 0.5 * (G[:, 0, 0] + G[:, 1, 1])
    vec = np.stack([0.5 * np.einsum("nab,ba->n", G, S) for S in PAULI], 1)   # sigma components: tr(G sigma)/2
    A = np.stack([np.ones(len(k)), (k ** 2).sum(1) + (kp ** 2).sum(1), (k * kp).sum(1), (k ** 2).sum(1) - (kp ** 2).sum(1)], 1)
    cs = np.linalg.lstsq(A.astype(complex), s, rcond=None)[0]
    cr = np.cross(kp, k)                                                       # k' x k
    gso = np.linalg.lstsq((1j * cr).reshape(-1, 1), vec.reshape(-1), rcond=None)[0][0]
    resid_so = np.linalg.norm(vec - 1j * gso * cr) / np.linalg.norm(vec)
    return dict(g0=complex(cs[0]), gk=complex(cs[1]), gx=complex(cs[2]), g_asym=complex(cs[3]), gso=complex(gso),
                so_shape_resid=float(resid_so))


# =====================================================================================================
#  coarse potential: band-limited restriction of the fine input -Z alpha / r
# =====================================================================================================
def coulomb_bl(N, h, Z=1.0, centre=(0.0, 0.0, 0.0), sig_cells=3.0):
    """V(X) on an N^3 grid of cell h (lambda_C), nucleus at `centre` (in cells, relative to the grid point N//2):
    band-limited (coarse Fourier modes only) restriction of -Z alpha / r inside the box. Split -1/r = -erf(r/s)/r
    (smooth: its grid samples are band-limited up to e^{-(pi s/2h)^2} ~ 1e-10) + -erfc(r/s)/r (short range: exact
    Fourier transform, -4 pi (1 - e^{-q^2 s^2/4})/q^2)."""
    s = sig_cells * h
    x = (np.arange(N) - N // 2) * h
    c = np.array(centre) * h
    X, Y, Zc = np.meshgrid(x - c[0], x - c[1], x - c[2], indexing="ij", sparse=True)
    r = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        Vl = np.where(r > 0, erf(r / s) / r, 2 / (np.sqrt(np.pi) * s))
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij", sparse=True)
    q2 = KX ** 2 + KY ** 2 + KZ ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        F = np.where(q2 > 0, 4 * np.pi * (1 - np.exp(-q2 * s * s / 4)) / q2, np.pi * s * s)
    # Vs(X) = (1/Vol) sum_q F(q) e^{iq.(X - c)}, X = x0 + n h (periodic images of erfc negligible)
    F = F * np.exp(-1j * (KX * c[0] + KY * c[1] + KZ * c[2])) * np.exp(1j * (KX + KY + KZ) * x[0])
    Vs = sfft.ifftn(F, workers=1) * (N ** 3) / (N * h) ** 3
    V = -Z * zz.ALPHA * (Vl + Vs.real)
    return V, x


def coulomb_cellavg(N, h, Z=1.0):
    """the OLD hand-written restriction (agg_zz_atom.cell_avg_inv_r, nucleus at a cell corner), for comparison only"""
    import agg_zz_atom as aza
    inv, x = aza.cell_avg_inv_r(N, h, "corner")
    return -Z * zz.ALPHA * inv, x


class GalerkinField:
    """the coarse potential operator in the plane-wave descriptor convention, WITHOUT aliasing:
    <k'|V_c|k> = Vt(k' - k) for k, k' in the coarse Brillouin cube, q = k' - k up to twice the cube.
    Exact form: the coarse cell carries the field at its 8 half-cell sub-points (s in {0,1}^3): V_c phi =
    (1/8) sum_s I_s^dag ( V_s . I_s phi ), I_s = band-limited interpolation to the sub-lattice s (a phase in k-space),
    V_s = the field restricted to the Fourier modes of the 2N grid, sampled on sub-lattice s.
    (Splitting -1/r = -erf(r/sig)/r [smooth, sampled] - erfc(r/sig)/r [exact Fourier transform].)"""

    def __init__(self, N, h, Z=1.0, centre=(0.0, 0.0, 0.0), sig_half_cells=3.0):
        self.N, self.h = N, h
        hh = h / 2
        s = sig_half_cells * hh
        x = (np.arange(N) - N // 2) * h
        self.x = x
        c = np.array(centre) * h
        k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
        self.k1 = k1
        # the 2N-grid frequencies as k + G, G in {0, other}: other = k + 2pi/h (k<0) or k - 2pi/h (k>=0)
        k1o = np.where(k1 < 0, k1 + 2 * np.pi / h, k1 - 2 * np.pi / h)
        # Nyquist of the 2N grid (k = -pi/h -> +pi/h representative); both representatives counted once
        self.subs = [(sx, sy, sz) for sx in (0, 1) for sy in (0, 1) for sz in (0, 1)]
        self.V = {}
        vol = (N * h) ** 3
        for sub in self.subs:
            off = np.array(sub) * hh
            X, Y, Zc = np.meshgrid(x + off[0] - c[0], x + off[1] - c[1], x + off[2] - c[2], indexing="ij", sparse=True)
            r = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
            with np.errstate(divide="ignore", invalid="ignore"):
                Vl = np.where(r > 0, erf(r / s) / r, 2 / (np.sqrt(np.pi) * s))
            Vs = np.zeros((N, N, N))
            for m in self.subs:
                ks = [k1o if m[i] else k1 for i in range(3)]
                KX, KY, KZ = np.meshgrid(*ks, indexing="ij", sparse=True)
                q2 = KX ** 2 + KY ** 2 + KZ ** 2
                with np.errstate(divide="ignore", invalid="ignore"):
                    F = np.where(q2 > 0, 4 * np.pi * (1 - np.exp(-q2 * s * s / 4)) / q2, np.pi * s * s)
                # e^{iq.(X + off - c)}, X = x0 + n h: e^{i q x0} e^{i k n h} (e^{i G n h} = 1)
                ph = np.exp(1j * (KX * (x[0] + off[0] - c[0]) + KY * (x[0] + off[1] - c[1]) + KZ * (x[0] + off[2] - c[2])))
                Vs += (sfft.ifftn(F * ph, workers=1) * N ** 3).real / vol
                del F, ph, q2
            self.V[sub] = (-Z * zz.ALPHA * (Vl + Vs)).astype(np.float32 if N > 192 else np.float64)
            del Vl, Vs
        self.ph1 = np.exp(1j * k1 * hh)                    # interpolation phase per axis for s = 1

    def phase(self, sub):
        p = [self.ph1 if sub[i] else np.ones(self.N) for i in range(3)]
        return p[0][:, None, None] * p[1][None, :, None] * p[2][None, None, :]

    def apply_k(self, Fphi):
        """Fphi = FFT of phi (N^3); returns FFT of V_c phi (in-place temporaries: 3 complex N^3 arrays at most)"""
        out = np.zeros_like(Fphi)
        one = np.ones(self.N)
        for sub in self.subs:
            p = [self.ph1 if sub[i] else one for i in range(3)]
            tmp = Fphi * p[0][:, None, None]
            tmp *= p[1][None, :, None]
            tmp *= p[2][None, None, :]
            tmp = sfft.ifftn(tmp, overwrite_x=True, workers=1)
            tmp *= self.V[sub]
            tmp = sfft.fftn(tmp, overwrite_x=True, workers=1)
            tmp *= np.conj(p[0])[:, None, None]
            tmp *= np.conj(p[1])[None, :, None]
            tmp *= np.conj(p[2])[None, None, :]
            out += tmp
            del tmp
        out /= 8
        return out

    def apply(self, phi):
        return ifftn(self.apply_k(fftn(phi)))

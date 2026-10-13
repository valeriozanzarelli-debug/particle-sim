"""emiss_radial.py - the atom in the RADIAL form of the same zig-zag rule (zz_core.RadialWalk: one Dirac partial
wave per walk, same flip, same Coulomb phase), for the one-excitation quantum run of the decay 2p3/2 (m = 3/2) ->
1s (m = 1/2).  Lets us use small cells and any Z (also hydrogen, Z = 1).

Channels: the vacuum sector psi is the kappa = -2 (p3/2), m = 3/2 radial wave; the photon sector carries the
kappa = -1 (s1/2), m = 1/2 radial wave.  By symmetry (the field a sigma+ dipole makes at its own position is along
the dipole) no other channel is reached at first order: exact in the continuum, a declared truncation otherwise.
Angular algebra (geometry, Clebsch-Gordan; checked numerically in the study):
    <1s 1/2 | r | 2p3/2 3/2>  =  R * VEC,   VEC = -(1, i, 0)/sqrt6,   R = Int (g1 g2 + f1 f2) r dr
"""
import numpy as np
import numba as nb
from zz_core import RadialWalk, ALPHA, dirac_formula

VEC = -np.array([1, 1j, 0]) / np.sqrt(6)


@nb.njit(cache=True)
def _radial_tick(chi, M, refl0, refl1, buf):
    N = chi.shape[1]
    p = chi[0]
    q = chi[1]
    pn = buf[0]
    qn = buf[1]
    pn[0] = refl0 * q[0]
    for i in range(1, N):
        pn[i] = p[i - 1]
    for i in range(N - 1):
        qn[i] = q[i + 1]
    qn[N - 1] = refl1 * p[N - 1]
    for i in range(N):
        a = pn[i]
        b = qn[i]
        chi[0, i] = M[i, 0, 0] * a + M[i, 0, 1] * b
        chi[1, i] = M[i, 1, 0] * a + M[i, 1, 1] * b


class RadialAtom:
    def __init__(self, Z, a, rmax, kappa_e=-2, kappa_g=-1):
        self.Z, self.a = Z, a
        Za = Z * ALPHA
        self.we = RadialWalk(kappa_e, a, rmax, Zalpha=Za)
        self.wg = RadialWalk(kappa_g, a, rmax, Zalpha=Za)
        self.r = self.we.r
        self.buf = np.empty((2, self.we.N), np.complex128)

    def tick(self, chi, t, eE=None):
        _radial_tick(chi, self.we.M, self.we.refl0, self.we.refl1, self.buf)
        return chi

    def tick_ground(self, chi, t):
        _radial_tick(chi, self.wg.M, self.wg.refl0, self.wg.refl1, self.buf)
        return chi

    def dip_matrix(self, phi, psi):
        R = np.sum(self.r * (np.conj(phi[0]) * psi[0] + np.conj(phi[1]) * psi[1]))
        return R * VEC

    def add_rphi(self, psi, coef, phi):
        """psi += P_channel (coef . r) phi = (coef . conj(VEC)) r phi"""
        c = np.dot(coef, np.conj(VEC))
        psi += c * self.r[None, :] * phi

    def eigen(self, which, n_ring):
        w = self.we if which == "e" else self.wg
        kap = w.kappa
        om, u, res = w.eigen(dirac_formula(n_ring, kap, self.Z * ALPHA) - 0.1 * (self.Z * ALPHA)**2 / n_ring**2)
        # u is the palindromic state; the walk's eigenvector at the tick boundary is chi (Mh u): use chi directly
        lam0 = np.exp(-1j * om * self.a)
        import scipy.sparse as sp, scipy.sparse.linalg as spl
        Wm = w.matrix()
        lu = spl.splu((Wm - lam0 * sp.identity(2 * w.N, format="csc")).tocsc())
        x = np.ones(2 * w.N, complex)
        for _ in range(6):
            x = lu.solve(x); x /= np.linalg.norm(x)
        lam = np.vdot(x, Wm @ x)
        res = np.linalg.norm(Wm @ x - lam * x)
        chi = x.reshape(w.N, 2).T.copy()
        return -np.angle(lam) / self.a, chi, res


@nb.njit(cache=True)
def _qe_block(psi, phi, Me, r0e, r1e, Mg, r0g, r1g, cr, sub, buf):
    """sub ticks of the one-excitation step with a constant photon field: psi, phi advance by their walks,
    then psi += cr * phi (cr = -i a e (E.conj(VEC)) r, per site)."""
    N = psi.shape[1]
    for s in range(sub):
        _radial_tick(psi, Me, r0e, r1e, buf)
        _radial_tick(phi, Mg, r0g, r1g, buf)
        for i in range(N):
            psi[0, i] += cr[i] * phi[0, i]
            psi[1, i] += cr[i] * phi[1, i]


def _block(self, psi, phis, fields, sub, a, e):
    c = -1j * a * e * np.dot(fields[0], np.conj(VEC))
    _qe_block(psi, phis[0], self.we.M, self.we.refl0, self.we.refl1, self.wg.M, self.wg.refl0, self.wg.refl1,
              (c * self.r).astype(np.complex128), sub, self.buf)


RadialAtom.qe_block = _block

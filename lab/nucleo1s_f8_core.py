"""nucleo1s_f8_core.py - phase 8 (docs/NUCLEO_1S.md section 18): the project of section 16.2.

Inside (r < R_f = (J0 + 1/2) h): the zig-zag radial rule itself on the zoom chain (nucleo1s_core.mesh_zoom), down to
2^-12 of the fine cell of chain A; it resolves the Coulomb cusp, where the Darwin term lives.
Outside: the TAYLOR-PINNED learned level of docs/AGGLOMERATO_COMPLETO.md (lab/out/aggfull_relearnT.pkl, written by the
aggfull agent, READ ONLY here): B = 256 (cell a0/16) or B = 512 (cell a0/8), reduced to the s wave exactly as
nucleo1s_hybrid.CoarseRadial does for the old level (spin-blind part 0.5 (K00 + K11), averaged over 120 directions,
T(k) = (G_iso(k) - G_iso(0)) / tau, 1-D kernel t_n acting on u = r psi with the odd reflection at r = 0).
Glue: the cusp glue of section 12.4 (nucleo1s_glue.solve_glue mode "cusp"): ghosts j < J0 from the inner solution, the
coarse rule's error on the Coulomb cusp removed with u_c = r e^{-lam r}, amplitude from Kato's law (the equation itself).

The alpha^4 shift is measured as a DIFFERENCE between two Z at the same geometry in units of a0/Z:
  (Z = 1, B = 512, cell a0/8)  and  (Z = 2, B = 256, cell a0/16 = (a0/2)/8):
binding B(Z) = Z^2 eps + Z^4 delta (+ Z^6 ...); delta_meas = (B(2)/4 - B(1)) / 3. The non-relativistic part of the error
of the coarse rule (potential sampled at the cell centre, lattice kinetic operator, glue) is a function of the geometry
in a0/Z only, hence the same eps for both Z, and it cancels; what is left is the relativistic content (Z^4) plus the
part of the stencil error that is NOT the same between the two learned levels (they are two different fits).
Reference: the same formula applied to the Dirac levels (exam only).

The solver is a dense linear solve (no LOBPCG, no projection on a sector): the wrong-parity issue of aggfull section 13
(Atom(..., clean=True)) does not arise here; the s-wave reduction keeps only the isotropic part of G(k).
Units hbar = c = m_e = 1.  Inputs: m_e, alpha, Z.  Nothing fitted.
"""
import os, sys, pickle
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_hybrid import CoarseRadial, zz_regular, fib_dirs, TH_A, A0  # noqa: E402
from nucleo1s_core import mesh_zoom, dirac_binding, ALPHA  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


class CoarseT(CoarseRadial):
    """band = 'learned': the aggfull Taylor-pinned level; 'exact': the fine walk's Dirac band sqrt(1+k^2)-1 (control,
    spectral kernel on the same cells); 'old': the old level of agg_zz_3d_A.pkl (comparison)."""

    def __init__(self, B=256, band="learned", ndir=120, nk=1024, nfine=32768, vmode="point"):
        self.vmode, self.band = vmode, band
        self.B, self.h, self.tau = B, B * TH_A, B * TH_A
        ks = np.linspace(0, np.pi, nk)
        self.ks = ks
        if band == "exact":
            self.spin_part = 0.0
            self.aniso_1a0 = 0.0
            kp = ks / self.h
            self.Tk = np.sqrt(1 + kp**2) - 1
            self.G0 = 0.0
        else:
            if band == "learned":
                lv = [l for l in pickle.load(open(os.path.join(OUT, "aggfull_relearnT.pkl"), "rb")) if l["B"] == B][0]
            else:
                lv = [l for l in pickle.load(open(os.path.join(OUT, "agg_zz_3d_A.pkl"), "rb"))["levels"] if l["B"] == B][0]
            offs = np.asarray(lv["offs"], float)
            K = lv["K"]
            Gs = np.real(0.5 * (K[0, 0] + K[1, 1]))
            self.spin_part = float(np.max(np.abs(K[0, 1])) / np.max(np.abs(Gs)))
            dirs = fib_dirs(ndir)
            Giso = np.zeros(nk)
            for c0 in range(0, nk, 128):
                kv = ks[c0:c0 + 128, None, None] * dirs[None, :, :]
                ph = np.exp(1j * np.tensordot(kv, offs.T, axes=(-1, 0)))
                Giso[c0:c0 + 128] = np.real(ph @ Gs).mean(axis=1)
            self.G0 = Giso[0]
            self.Tk = (Giso - self.G0) / self.tau
            kk = (self.h / A0) * dirs
            g1 = np.real(np.exp(1j * kk @ offs.T) @ Gs)
            self.aniso_1a0 = float((g1.max() - g1.min()) / max(abs(g1.mean() - self.G0), 1e-300))
        self._spl = CubicSpline(ks, self.Tk)
        self._kf = np.linspace(0, np.pi, nfine)
        self._Tf = self._spl(self._kf)
        kphys = ks[1:6] / self.h
        self.m_kin = float(np.mean(kphys**2 / (2 * self.Tk[1:6])))
        self._tcache = {}

    def band_dev(self):
        """relative deviation of the s-wave band from the exact Dirac band, at k h = pi/4, pi/2, 3pi/4, pi"""
        kp = self.ks / self.h
        ex = np.sqrt(1 + kp**2) - 1
        out = {}
        for f in (0.25, 0.5, 0.75, 1.0):
            i = int(round(f * (len(self.ks) - 1)))
            out[f"{f:g}pi"] = float(self.Tk[i] / ex[i] - 1)
        # the relativistic k^4 part, from the small-k curvature: T = k^2/2m - c4 k^4 + ...; Dirac c4 = 1/8
        sel = slice(1, 40)
        x = kp[sel]**2
        p = np.polyfit(x, self.Tk[sel] / x, 2)
        out["k4_over_dirac"] = float(-p[1] / 0.125)
        out["m_kin"] = self.m_kin
        return out


class Glue8:
    """cusp glue at energy E; rmax in lambda_C; J0 coarse ghost cells; inner zoom chain to R_f = (J0 + 1/2) h."""

    def __init__(self, c, Z, J0, rmax, zoom_levels=12, npl=16):
        self.c, self.Z, self.J0 = c, Z, J0
        self.Za = Za = Z * ALPHA
        h = c.h
        self.N = N = int(np.ceil(rmax / h))
        self.rj = rj = (np.arange(N) + 0.5) * h
        self.Hm = c.H(N, Za)
        self.K = self.Hm - np.diag(c.V(N, Za))
        self.mesh = mesh_zoom(TH_A, (J0 + 0.5) * h, zoom_levels, npl)
        lam = self.lam = 2 * c.m_kin * Za
        uc = rj * np.exp(-lam * rj)
        self.tauc = self.K @ uc + (-2 * lam + lam**2 * rj) * np.exp(-lam * rj) / (2 * c.m_kin)

    def parts(self, E):
        r, g, f = zz_regular(self.mesh, -1, self.Za, 1 + E)
        G = g.real
        self.imag_ratio = float(np.max(np.abs(g.imag)) / np.max(np.abs(g.real)))
        ug = CubicSpline(r, G)(self.rj[:self.J0])
        sel = (r > 0.5) & (r < 4.0)
        a1, a2 = np.polyfit(r[sel], G[sel] / r[sel], 1)[::-1]
        s = (-2 * self.c.m_kin * self.Za * a1) / (-2 * self.lam)
        return ug, s * self.tauc

    def res(self, E, return_u=False):
        J0, N = self.J0, self.N
        ug, t = self.parts(E)
        A = self.Hm - E * np.eye(N)
        uo = np.linalg.solve(A[J0:, J0:], -A[J0:, :J0] @ ug + t[J0:])
        r = (A[J0 - 1, :J0] @ ug + A[J0 - 1, J0:] @ uo - t[J0 - 1]) / np.abs(ug).max()
        if return_u:
            return r, np.concatenate([ug, uo])
        return r

    def solve(self, nst, span=0.03, npts=25):
        Eref = -dirac_binding(nst, -1, self.Za)
        Es = np.linspace(Eref * (1 + span), Eref * (1 - span), npts)
        rs = [self.res(E) for E in Es]
        br = [i for i in range(npts - 1) if np.sign(rs[i]) != np.sign(rs[i + 1])]
        if not br:
            return None
        i = min(br, key=lambda i: abs(0.5 * (Es[i] + Es[i + 1]) - Eref))
        E = brentq(self.res, Es[i], Es[i + 1], xtol=1e-20, rtol=1e-15, maxiter=200)
        _, u = self.res(E, return_u=True)
        big = np.abs(u) > 1e-6 * np.abs(u).max()
        v = u[big]
        nodes = int(np.sum(np.sign(v[1:]) * np.sign(v[:-1]) < 0))
        wall = float(np.abs(u[-20:]).max() / np.abs(u).max())
        return dict(B=-E, rel=-E / dirac_binding(nst, -1, self.Za) - 1, nodes=nodes, wall=wall)


def coarse_only_levels(c, Z, rmax, nmax=3):
    N = int(np.ceil(rmax / c.h))
    v = np.linalg.eigvalsh(c.H(N, Z * ALPHA))
    return [-float(x) for x in v[:nmax]]


def delta_pair(B1, B2, Z1=1, Z2=2):
    """B(Z) = Z^2 eps + Z^4 delta: delta from two Z"""
    return (B2 / Z2**2 - B1 / Z1**2) / (Z2**2 - Z1**2)


def delta_dirac(n, Z1=1, Z2=2):
    return delta_pair(dirac_binding(n, -1, Z1 * ALPHA), dirac_binding(n, -1, Z2 * ALPHA), Z1, Z2)

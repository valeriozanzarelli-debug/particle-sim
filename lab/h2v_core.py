"""h2v_core.py - VIDEO H2: everything computed from the two-electron ECG wavefunction of H2 (bases rebuilt by
lab/h2v_wf.py from Debiti's lab/aggmol_h2.py) and from Debiti's Born-Oppenheimer curve.

  WF(R)                  the wavefunction: Psi(r1, r2) = (1/2) sum_k c_k nrm_k sum_ops phi_k(op r)   (normalised to 1)
  WF.psi(X)              Psi at points X (n, 2, 3)
  WF.pairs()             pair terms <phi_k | . | sum_ops phi_k'> as 6-D Gaussians (weight, matrix C, centre R0)
  WF.bins_xz(edges)      exact number of electrons per (x, z) cell (column density integrated over y and over the cell)
  WF.expect()            <T>, <V_ne>, <V_ee>, <1/r12>, norm (internal check against the solver)
  WF.prob_z(...)         probability of a two-electron event on (z1, z2) (Gauss-Legendre quadrature, exact Gaussians)
  WF.abs_z()             <|z|> of one electron
  vib()                  vibrational levels/wavefunctions on Debiti's curve (same recipe as aggmol_h2.py job_vib)
"""
import os, sys, json
import numpy as np
from scipy.special import erf, ndtr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_ecg as E

OUT = E.OUT
OPS = [(np.eye(2), 1.0), (np.eye(2), -1.0), (np.array([[0., 1.], [1., 0.]]), 1.0), (np.array([[0., 1.], [1., 0.]]), -1.0)]
RGRID = [1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7]
MP = 1836.15267343
AU_FS = 2.4188843265857e-2      # atomic unit of time in fs


class WF:
    def __init__(self, R, path=None):
        z = np.load(path or os.path.join(OUT, f"h2v_basis_R{R:.2f}.npz"))
        self.R = R
        self.A, self.s, self.w = z["A"], z["s"], z["c"] * z["nrm"]
        self.E_el = float(z["E_el"])
        if np.sum(self.w) < 0:
            self.w = -self.w
        Aa, sa, wa = [], [], []
        for P, inv in OPS:
            Aa.append(np.einsum("ji,kjl,lm->kim", P, self.A, P))
            sa.append(inv * np.einsum("ij,kjx->kix", P, self.s))
            wa.append(self.w)
        self.Aa, self.za, self.wa = np.concatenate(Aa), np.concatenate(sa)[..., 2], 0.5 * np.concatenate(wa)
        self._pairs = None

    def psi(self, X):
        """X (n, 2, 3) -> Psi (n). The exponent is linear in 6 features of the walker: one matrix product"""
        if not hasattr(self, "_G"):
            a11, a12, a22 = self.Aa[:, 0, 0], self.Aa[:, 0, 1], self.Aa[:, 1, 1]
            s1, s2 = self.za[:, 0], self.za[:, 1]
            self._G = np.stack([a11, 2 * a12, a22, -2 * (a11 * s1 + a12 * s2), -2 * (a12 * s1 + a22 * s2),
                                a11 * s1 ** 2 + 2 * a12 * s1 * s2 + a22 * s2 ** 2])
        x1, x2 = X[:, 0], X[:, 1]
        F = np.stack([np.sum(x1 * x1, 1), np.sum(x1 * x2, 1), np.sum(x2 * x2, 1), x1[:, 2], x2[:, 2], np.ones(len(X))], 1)
        out = np.empty(len(X))
        for i in range(0, len(X), 4096):
            out[i:i + 4096] = np.exp(-(F[i:i + 4096] @ self._G)) @ self.wa
        return out

    # ---------------------------------------------------------------- analytic pieces
    def pairs(self):
        if self._pairs is not None:
            return self._pairs
        K = len(self.A)
        out = []
        for P, inv in OPS:
            Ap = np.einsum("ji,kjl,lm->kim", P, self.A, P)
            sp = inv * np.einsum("ij,kjx->kix", P, self.s)
            A1, s1 = np.repeat(self.A, K, 0), np.repeat(self.s, K, 0)
            A2, s2 = np.tile(Ap, (K, 1, 1)), np.tile(sp, (K, 1, 1))
            C = A1 + A2
            Ci = np.linalg.inv(C)
            det = np.linalg.det(C)
            v = A1 @ s1 + A2 @ s2
            R0 = Ci @ v
            c0 = np.einsum("nix,nij,njx->n", s1, A1, s1) + np.einsum("nix,nij,njx->n", s2, A2, s2)
            ex = np.einsum("nix,nij,njx->n", v, Ci, v) - c0
            ww = np.repeat(self.w, K) * np.tile(self.w, K)
            S = (np.pi ** 2 / det) ** 1.5 * np.exp(ex)
            out.append(dict(A1=A1, s1=s1, A2=A2, s2=s2, C=C, Ci=Ci, det=det, R0=R0, ex=ex, ww=ww, S=S))
        self._pairs = {k: np.concatenate([o[k] for o in out]) for k in out[0]}
        return self._pairs

    def expect(self):
        p = self.pairs()
        S, ww = p["S"], p["ww"]
        norm = np.sum(ww * S)
        Sig = 0.5 * p["Ci"]
        M = p["A1"] @ p["A2"]
        d1, d2 = p["R0"] - p["s1"], p["R0"] - p["s2"]
        T = 4 * (3 * np.einsum("nij,nji->n", M, Sig) + np.einsum("nix,nij,njx->n", d1, M, d2)) / (2 * E.M_KIN)
        Vne = 0.0
        for i in range(2):
            for X in ((0, 0, self.R / 2), (0, 0, -self.R / 2)):
                Vne = Vne - E.coul(p["R0"][:, i, :] - np.array(X), Sig[:, i, i])
        iee = E.coul(p["R0"][:, 0, :] - p["R0"][:, 1, :], Sig[:, 0, 0] + Sig[:, 1, 1] - 2 * Sig[:, 0, 1])
        f = lambda a: float(np.sum(ww * S * a) / norm)
        o = dict(norm_vs_solver=float(norm), T=f(T), Vne=f(Vne), inv_r12=f(iee))
        o["E_el"] = o["T"] + o["Vne"] + o["inv_r12"]
        o["E_el_solver"] = self.E_el
        return o

    def _gauss1(self):
        """one-electron marginal Gaussians: list of (weight, alpha, z-centre) for both electrons and both inversion images"""
        p = self.pairs()
        C, R0, wt = p["C"], p["R0"], p["ww"] * np.exp(p["ex"])
        det = p["det"]
        a1 = det / C[:, 1, 1]
        a2 = det / C[:, 0, 0]
        w1 = wt * (np.pi / C[:, 1, 1]) ** 1.5
        w2 = wt * (np.pi / C[:, 0, 0]) ** 1.5
        norm = np.sum(p["ww"] * p["S"])
        W = np.concatenate([w1, w2]) / norm
        A = np.concatenate([a1, a2])
        Z = np.concatenate([R0[:, 0, 2], R0[:, 1, 2]])
        # gerade: average with the inverted image
        return np.concatenate([W, W]) / 2, np.concatenate([A, A]), np.concatenate([Z, -Z])

    def bins_xz(self, edges):
        """exact expected number of electrons (of 2) per cell of the (x, z) grid (integrated over y and over the cell)"""
        W, A, Z = self._gauss1()
        out = np.zeros((len(edges) - 1, len(edges) - 1))
        for i in range(0, len(W), 50000):
            w, a, z = W[i:i + 50000], A[i:i + 50000], Z[i:i + 50000]
            sa = np.sqrt(a)
            # int dy exp(-a y^2) = sqrt(pi/a); int over a cell of exp(-a (x-c)^2) = sqrt(pi/a)/2 [erf]
            ex = 0.5 * np.sqrt(np.pi / a)[:, None] * np.diff(erf(sa[:, None] * edges[None]), axis=1)
            ez = 0.5 * np.sqrt(np.pi / a)[:, None] * np.diff(erf(sa[:, None] * (edges[None] - z[:, None])), axis=1)
            out += ((w * np.sqrt(np.pi / a))[:, None] * ex).T @ ez
        return out          # [ix, iz]

    def abs_z(self):
        W, A, Z = self._gauss1()
        sig = np.sqrt(0.5 / A)
        tot = np.sum(W * (np.pi / A) ** 1.5)
        m = sig * np.sqrt(2 / np.pi) * np.exp(-Z ** 2 / (2 * sig ** 2)) + Z * (1 - 2 * ndtr(-Z / sig))
        return float(np.sum(W * (np.pi / A) ** 1.5 * m) / tot)

    def prob_z(self, events, ngl=96):
        """events: list of (lo1, hi1, lo2, hi2) on (z1, z2); returns the average probability over the list
        (pass the list of images of an event under swap and inversion so the observable is symmetric)"""
        p = self.pairs()
        Sig = 0.5 * p["Ci"]
        m1, m2 = p["R0"][:, 0, 2], p["R0"][:, 1, 2]
        s1, s2 = np.sqrt(Sig[:, 0, 0]), np.sqrt(Sig[:, 1, 1])
        rho = Sig[:, 0, 1] / (s1 * s2)
        wt = p["ww"] * p["S"]
        norm = np.sum(wt)
        t, gw = np.polynomial.legendre.leggauss(ngl)
        tot = 0.0
        for lo1, hi1, lo2, hi2 in events:
            u0 = np.clip((lo1 - m1) / s1, -9, 9)
            u1 = np.clip((hi1 - m1) / s1, -9, 9)
            u = u0[:, None] + (u1 - u0)[:, None] * (t[None] + 1) / 2
            jac = (u1 - u0)[:, None] / 2
            phi = np.exp(-u ** 2 / 2) / np.sqrt(2 * np.pi)
            cm = m2[:, None] + rho[:, None] * s2[:, None] * u
            cs = (s2 * np.sqrt(np.maximum(1 - rho ** 2, 1e-300)))[:, None]
            pc = ndtr((hi2 - cm) / cs) - ndtr((lo2 - cm) / cs)
            P = np.sum(gw[None] * jac * phi * pc, axis=1)
            tot += np.sum(wt * P) / norm
        return float(tot / len(events))


def images(lo1, hi1, lo2, hi2):
    """the event, its swap, its inversion, both"""
    ev = [(lo1, hi1, lo2, hi2), (lo2, hi2, lo1, hi1)]
    return ev + [(-b, -a, -d, -c) for a, b, c, d in ev]


def vib():
    """same recipe as lab/aggmol_h2.py job_vib (spline on Debiti's 18 points, x in [0.8, 6], 6000 points)"""
    from scipy.interpolate import CubicSpline
    from scipy.linalg import eigh_tridiagonal
    cur = json.load(open(os.path.join(OUT, "aggmol_h2.json")))["curve"]
    Rs = np.array(sorted(v["R"] for v in cur.values()))
    Et = np.array([cur[f"{R:.2f}"]["E_tot"] for R in Rs])
    cs = CubicSpline(Rs, Et)
    mu = MP / 2
    x = np.linspace(0.8, 6.0, 6000); dx = x[1] - x[0]
    d = 1 / (mu * dx ** 2) + cs(x); o = -0.5 / (mu * dx ** 2) * np.ones(len(x) - 1)
    lev, vec = eigh_tridiagonal(d, o, select="i", select_range=(0, 1))
    vec = vec / np.sqrt(dx)
    for k in range(2):
        if vec[np.argmax(np.abs(vec[:, k])), k] < 0:
            vec[:, k] *= -1
    return dict(x=x, dx=dx, E=lev, chi=vec, spline=cs)


def packet(V, t_au):
    """|chi(R, t)|^2 for (chi0 + chi1)/sqrt(2); t in atomic units"""
    w = V["E"][1] - V["E"][0]
    c0, c1 = V["chi"][:, 0], V["chi"][:, 1]
    return 0.5 * (c0 ** 2 + c1 ** 2 + 2 * c0 * c1 * np.cos(w * t_au))


def grid_weights(V, dens):
    """weights of the R grid points (cells between mid-points; tails lumped into the end points, declared)"""
    x = V["x"]
    mids = np.array([(a + b) / 2 for a, b in zip(RGRID[:-1], RGRID[1:])])
    idx = np.searchsorted(mids, x)
    w = np.bincount(idx, dens * V["dx"], len(RGRID))
    return w / w.sum()


def rho1_table(wf, rmax=6.0, nr=241, nz=481):
    """one-electron marginal rho1(rho_cyl, z) (integrates to 1) on a grid; exact sum of the Gaussians"""
    W, A, Z = wf._gauss1()
    r = np.linspace(0, rmax, nr); z = np.linspace(-rmax, rmax, nz)
    T = np.zeros((nr, nz))
    for i in range(0, len(W), 20000):
        w, a, c = W[i:i + 20000], A[i:i + 20000], Z[i:i + 20000]
        T += (w[:, None] * np.exp(-a[:, None] * r[None] ** 2)).T @ np.exp(-a[:, None] * (z[None] - c[:, None]) ** 2)
    return r, z, T / 2


class Gmap:
    """pair correlation of the other electron given one at r1: g(r2 | r1) = |Psi(r1, r2)|^2 / (rho1(r1) rho1(r2));
    the average of g over balls sampled from rho1 is 1 (check)"""
    def __init__(self, wf, cache=None):
        from scipy.interpolate import RegularGridInterpolator
        if cache and os.path.exists(cache):
            z = np.load(cache); r, zz, T = z["r"], z["z"], z["T"]
        else:
            r, zz, T = rho1_table(wf)
            if cache:
                np.savez(cache, r=r, z=zz, T=T)
        self.wf = wf
        self.ip = RegularGridInterpolator((r, zz), T, bounds_error=False, fill_value=None)

    def rho1(self, P):
        rc = np.minimum(np.hypot(P[:, 0], P[:, 1]), 6.0); z = np.clip(P[:, 2], -6, 6)
        return np.maximum(self.ip(np.stack([rc, z], 1)), 1e-300)

    def g(self, r1, P):
        X = np.concatenate([np.broadcast_to(r1, (len(P), 1, 3)), P[:, None, :]], 1)
        return self.wf.psi(X) ** 2 / (self.rho1(np.asarray(r1)[None])[0] * self.rho1(P))

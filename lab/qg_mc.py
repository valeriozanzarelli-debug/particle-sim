"""qg_mc.py - PART 1 engine: SU(3) pure gauge theory in 3+1 D imaginary time (Wilson action), numpy.

U[mu] has shape (3, 3, T, L, L, L); mu = 0 is time, 1..3 = x, y, z.  S = beta sum_P (1 - 1/3 Re tr U_P).
Update: Cabibbo-Marinari (3 SU(2) subgroups) with Kennedy-Pendleton heat bath, plus over-relaxation, on a
checkerboard (links of the same direction on sites of the same parity do not share plaquettes).

Measurements (all gauge invariant):
  plaquette; APE smearing (spatial for the quark lines, 4-D for the measurement field: knobs, declared);
  planar Wilson loops W(R, T) -> static q-qbar potential -> string tension;
  baryonic (3-quark) Wilson loop with the junction at the origin and the 3 quarks at d*s_k*e_k, contracted with
  eps_abc eps_a'b'c' (colour singlet) -> 3Q energy, and its correlation with the plaquette field (action density)
  on the middle time slice -> the flux distribution around three static quarks.
"""
import numpy as np

from qg_su3 import (SUB, dag, eps_eps, ident, kp, left_sub, mm, mm3, qdag, qmul, rand_su2, reunit, rtr, sh,
                    sub_quat)


class Lat:
    def __init__(self, T, L, beta, seed=0):
        self.T, self.L, self.beta = T, L, beta
        self.shape = (T, L, L, L)
        self.rng = np.random.default_rng(seed)
        self.U = np.stack([ident(self.shape) for _ in range(4)])
        g = np.indices(self.shape).sum(0) % 2
        self.par = [g == 0, g == 1]

    # ------------------------------------------------------------ staples
    def staple(self, mu, U=None):
        """V with local action of U_mu(x) = -(beta/3) Re tr(U_mu(x) V)."""
        U = self.U if U is None else U
        V = np.zeros_like(U[0])
        for nu in range(4):
            if nu == mu:
                continue
            V += mm3(sh(U[nu], mu), dag(sh(U[mu], nu)), dag(U[nu]))
            Unm = sh(U[nu], nu, -1)
            V += mm3(dag(sh(Unm, mu)), dag(sh(U[mu], nu, -1)), Unm)
        return V

    # ------------------------------------------------------------ updates
    def _update(self, heat):
        b = self.beta
        for mu in range(4):
            for p in (0, 1):
                m = self.par[p]
                V = self.staple(mu)[:, :, m]
                Um = self.U[mu][:, :, m]
                W = mm(Um, V)
                for (i, j) in SUB:
                    r = sub_quat(W, i, j)
                    k = np.sqrt((r * r).sum(0))
                    k = np.maximum(k, 1e-12)
                    v = r / k
                    if heat:
                        x = rand_su2(kp(2.0 * b * k / 3.0, self.rng), self.rng)
                        a = qmul(x, qdag(v))
                    else:
                        vd = qdag(v)
                        a = qmul(vd, vd)
                    left_sub(Um, a, i, j)
                    left_sub(W, a, i, j)
                self.U[mu][:, :, m] = Um

    def heatbath(self):
        self._update(True)

    def overrelax(self):
        self._update(False)

    def reunit(self):
        for mu in range(4):
            self.U[mu] = reunit(self.U[mu])

    # ------------------------------------------------------------ observables
    def plaq_field(self, U, mu, nu):
        return rtr(mm(mm(U[mu], sh(U[nu], mu)), mm(dag(sh(U[mu], nu)), dag(U[nu])))) / 3.0

    def plaquette(self, U=None):
        U = self.U if U is None else U
        s, n = 0.0, 0
        for mu in range(4):
            for nu in range(mu + 1, 4):
                s += self.plaq_field(U, mu, nu).mean()
                n += 1
        return s / n

    def ape(self, U, alpha, n, dirs):
        """APE smearing of the links in `dirs` using staples in `dirs` only:
        U' = Proj_SU3[(1 - alpha) U + alpha/(2(nd-1)) sum staples^+]  (projection = Gram-Schmidt, declared)."""
        U = U.copy()
        nd = len(dirs)
        for _ in range(n):
            new = U.copy()
            for mu in dirs:
                V = np.zeros_like(U[0])
                for nu in dirs:
                    if nu == mu:
                        continue
                    V += mm3(U[nu], sh(U[mu], nu), dag(sh(U[nu], mu)))
                    Unm = sh(U[nu], nu, -1)
                    V += mm3(dag(Unm), sh(U[mu], nu, -1), sh(Unm, mu))
                new[mu] = reunit((1 - alpha) * U[mu] + alpha / (2 * (nd - 1)) * V)
            U = new
        return U

    def site_plaq(self, U):
        """per-site averages of Re tr/3 of the 12 temporal and the 12 spatial plaquettes touching the site."""
        pt = np.zeros(self.shape)
        ps = np.zeros(self.shape)
        for mu in range(4):
            for nu in range(mu + 1, 4):
                P = self.plaq_field(U, mu, nu)
                # the plaquette at x (plane mu,nu) touches x, x+mu, x+nu, x+mu+nu
                Pm = np.roll(P, 1, axis=mu)
                acc = (P + Pm + np.roll(P, 1, axis=nu) + np.roll(Pm, 1, axis=nu)) / 4.0
                if mu == 0:
                    pt += acc / 3.0
                else:
                    ps += acc / 3.0
        return pt, ps

    # ------------------------------------------------------------ lines
    def line(self, U, mu, n, sign=1):
        """product of n links from x along sign*e_mu, for all x: (3,3,*lat)."""
        out = ident(self.shape)
        if sign > 0:
            for k in range(n):
                out = mm(out, np.roll(U[mu], -k, axis=2 + mu))
        else:
            for k in range(1, n + 1):
                out = mm(out, dag(np.roll(U[mu], k, axis=2 + mu)))
        return out

    def wilson_planar(self, Us, Ut, Rmax, Tmax):
        """<W(R,T)> averaged over positions and the 3 spatial directions; Us = spatially smeared links (lines),
        Ut = links used for the temporal lines."""
        out = np.zeros((Rmax + 1, Tmax + 1))
        tl = [ident(self.shape)]
        for t in range(Tmax):
            tl.append(mm(tl[-1], np.roll(Ut[0], -t, axis=2)))
        for mu in (1, 2, 3):
            S = ident(self.shape)
            for R in range(1, Rmax + 1):
                S = mm(S, np.roll(Us[mu], -(R - 1), axis=2 + mu))
                for T in range(1, Tmax + 1):
                    # W = S(x,t) Tt(x+R,t) S(x,t+T)^+ Tt(x,t)^+
                    A = mm(S, np.roll(tl[T], -R, axis=2 + mu))
                    B = mm(tl[T], np.roll(S, -T, axis=2))
                    out[R, T] += (rtr(mm(A, dag(B))) / 3.0).mean() / 3.0
        return out

    def baryon(self, Us, Ut, d, Tlist, sgn):
        """baryonic Wilson loop with junction at x, quarks at x + d*sgn[k]*e_{k+1}; returns {T: W field (complex)}."""
        sp = [self.line(Us, k + 1, d, sgn[k]) for k in range(3)]
        tl = {0: ident(self.shape)}
        cur = ident(self.shape)
        for t in range(max(Tlist)):
            cur = mm(cur, np.roll(Ut[0], -t, axis=2))
            tl[t + 1] = cur
        out = {}
        for T in Tlist:
            M = []
            for k in range(3):
                q = np.roll(tl[T], -d * sgn[k], axis=3 + k)          # temporal line at x + q_k
                M.append(mm(mm(sp[k], q), dag(np.roll(sp[k], -T, axis=2))))
            out[T] = eps_eps(*M)
        return out

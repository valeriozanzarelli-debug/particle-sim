"""qg_ym.py - PART 2 engine: classical SU(3) Yang-Mills in real time on a 3-D periodic lattice (Hamiltonian,
temporal gauge A_0 = 0), leapfrog.

Fields: links U[i] (3,3,X,Y,Z) in SU(3), electric fields E[i] (3,3,X,Y,Z) hermitian traceless (E = E^a T^a).
H = sum_links tr(E^2) + 2 sum_plaq (3 - Re tr U_P)   (-> 1/2 E^a E^a + 1/2 B^a B^a for smooth fields; g = 1, a = 1).
Equations:  dU/dt = i E U ;   dE/dt = - [traceless hermitian part of (U S - (U S)^+)/(2i)],  S = staple sum.
Gauss law:  G(x) = sum_i [E_i(x) - U_i(x-i)^+ E_i(x-i) U_i(x-i)] = 0  (conserved exactly by this leapfrog, because the
link update exp(i dt E) commutes with E).

Gauge-invariant diagnostics per site: energy density (electric + magnetic part), Poynting vector
S_k = eps_kij 2 tr(E_i B_j) (E and B both moved to the site, B from the clover), used to carry the drawn balls.
"""
import numpy as np

from qg_su3 import TA, dag, expi, herm_tl, ident, mm, mm3, reunit, rtr, tr

DIRS = (0, 1, 2)


def shift(a, d, n=1):
    return np.roll(a, -n, axis=2 + d)


class YM:
    def __init__(self, shape):
        self.shape = tuple(shape)
        self.U = np.stack([ident(self.shape) for _ in DIRS])
        self.E = np.zeros_like(self.U)

    # -------------------------------------------------------------- forces
    def staple(self, i):
        U = self.U
        V = np.zeros_like(U[0])
        for j in DIRS:
            if j == i:
                continue
            V += mm3(shift(U[j], i), dag(shift(U[i], j)), dag(U[j]))
            Ujm = shift(U[j], j, -1)
            V += mm3(dag(shift(Ujm, i)), dag(shift(U[i], j, -1)), Ujm)
        return V

    def force(self):
        return np.stack([-herm_tl(mm(self.U[i], self.staple(i))) for i in DIRS])

    def step(self, dt, F=None):
        """one leapfrog step: E is at half steps (E(t - dt/2) -> E(t + dt/2)), U at integer steps."""
        F = self.force() if F is None else F
        self.E += dt * F
        for i in DIRS:
            self.U[i] = mm(expi(dt * self.E[i], order=10), self.U[i])

    def reunit(self):
        for i in DIRS:
            self.U[i] = reunit(self.U[i])

    # -------------------------------------------------------------- diagnostics
    def plaq(self, i, j):
        U = self.U
        return mm(mm(U[i], shift(U[j], i)), mm(dag(shift(U[i], j)), dag(U[j])))

    def gauss(self):
        G = np.zeros_like(self.E[0])
        for i in DIRS:
            Um = shift(self.U[i], i, -1)
            G += self.E[i] - mm3(dag(Um), shift(self.E[i], i, -1), Um)
        return G

    def energy_fields(self, E=None):
        """site energy densities (electric, magnetic); sums are exactly H's two parts."""
        E = self.E if E is None else E
        eE = np.zeros(self.shape)
        for i in DIRS:
            e = rtr(mm(E[i], E[i]))
            eE += 0.5 * (e + np.roll(e, 1, axis=i))
        eB = np.zeros(self.shape)
        for i, j in ((0, 1), (0, 2), (1, 2)):
            p = 2 * (3 - rtr(self.plaq(i, j)))
            pm = np.roll(p, 1, axis=i)
            eB += 0.25 * (p + pm + np.roll(p, 1, axis=j) + np.roll(pm, 1, axis=j))
        return eE, eB

    def site_E(self, E=None):
        E = self.E if E is None else E
        out = []
        for i in DIRS:
            Um = shift(self.U[i], i, -1)
            out.append(0.5 * (E[i] + mm3(dag(Um), shift(E[i], i, -1), Um)))
        return out

    def site_B(self):
        """clover: F_ij(x) = traceless hermitian part of (Q - Q^+)/(2i) / 4; B_k = F_ij, (i,j,k) cyclic."""
        U = self.U
        F = {}
        for i, j in ((1, 2), (2, 0), (0, 1)):
            Ui, Uj = U[i], U[j]
            Uim = shift(Ui, i, -1)              # U_i(x-i)
            Ujm = shift(Uj, j, -1)              # U_j(x-j)
            # four leaves, all from x back to x with the same orientation (i then j)
            Q = mm(mm(Ui, shift(Uj, i)), mm(dag(shift(Ui, j)), dag(Uj)))
            Q += mm(mm(Uj, dag(shift(Uim, j))), mm(dag(shift(Uj, i, -1)), Uim))
            Q += mm(mm(dag(Uim), dag(shift(Ujm, i, -1))), mm(shift(Uim, j, -1), Ujm))
            Q += mm(mm(dag(Ujm), shift(Ui, j, -1)), mm(shift(Ujm, i), dag(Ui)))
            F[(i, j)] = herm_tl(Q) / 4.0
        # B_x = F_yz, B_y = F_zx, B_z = F_xy
        return [F[(1, 2)], F[(2, 0)], F[(0, 1)]]

    def poynting(self):
        Es = self.site_E()
        B = self.site_B()
        S = np.zeros((3,) + self.shape)
        for k, (i, j) in enumerate(((1, 2), (2, 0), (0, 1))):
            # minus sign: our E is +dA/dt (physical E is -dA/dt); checked: a packet moving +x gives S_x > 0
            S[k] = -2 * (rtr(mm(Es[i], B[j])) - rtr(mm(Es[j], B[i])))
        return S


def comp_field(c, gen):
    """c: real field (*lat) -> c * T^gen as (3,3,*lat)."""
    return np.einsum("...,ij->ij...", c, TA[gen])

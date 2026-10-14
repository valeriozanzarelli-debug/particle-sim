"""tavg_radial.py - livello fine "a zoom radiale" per gli atomi da B in su (studio tavg_).

Hartree-Fock ristretto con energia del termine (ROHF di Froese Fischer) per il determinante di peso massimo del termine
di Hund (max S, poi max L): per B..Kr quel determinante e' un termine puro, quindi la sua energia e' l'energia del
termine. Operatore: quello non relativistico della regola agglomerata, p^2/(2 m_kin) - Z/r + 1/r_ij (m_kin -> 1 per
theta -> 0, scala esatta E(m) = m E(1)). Niente orbitali da tabelle: si parte da orbitali idrogenoidi nudi (carica Z)
e il campo medio si costruisce da solo.

Griglia: x uniforme, r = e^x, P(r) = sqrt(r) u(x). Cinetico con differenze finite di ordine 8; integrali con la regola
dei trapezi in x (spettrale per funzioni che decadono ai due estremi); potenziali di Hartree/scambio Y^k risolvendo
w'' - (k+1/2)^2 w = -(2k+1) sqrt(r) rho (stesse differenze finite: soluzione liscia, niente spigolo del nucleo r</r>).

Unita' atomiche (Eh, a0). Interfaccia principale:
    atom = Atom(Z, config)  # config: lista (n, l, occ)
    res = atom.scf()        # dict con E, orbitali, eps, ...
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import math
from functools import lru_cache
import numpy as np
import scipy.linalg as sla

L_LETTER = "spdf"


# ------------------------------------------------------------------------------------------------ angular algebra
@lru_cache(maxsize=None)
def wigner3j(j1, j2, j3, m1, m2, m3):
    if m1 + m2 + m3 != 0 or j3 < abs(j1 - j2) or j3 > j1 + j2 or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    f = math.factorial
    tri = f(j1 + j2 - j3) * f(j1 - j2 + j3) * f(-j1 + j2 + j3) / f(j1 + j2 + j3 + 1)
    pre = math.sqrt(tri * f(j1 + m1) * f(j1 - m1) * f(j2 + m2) * f(j2 - m2) * f(j3 + m3) * f(j3 - m3))
    s = 0.0
    for t in range(0, j1 + j2 + j3 + 1):
        d = [t, j3 - j2 + t + m1, j3 - j1 + t - m2, j1 + j2 - j3 - t, j1 - t - m1, j2 - t + m2]
        if min(d) < 0:
            continue
        s += (-1) ** t / (f(d[0]) * f(d[1]) * f(d[2]) * f(d[3]) * f(d[4]) * f(d[5]))
    return (-1) ** (j1 - j2 - m3) * pre * s


@lru_cache(maxsize=None)
def ck(k, l, m, lp, mp):
    """Condon-Shortley c^k(lm, l'm') = <lm|C^k_{m-m'}|l'm'>"""
    return (-1) ** m * math.sqrt((2 * l + 1) * (2 * lp + 1)) * wigner3j(l, k, lp, 0, 0, 0) * wigner3j(l, k, lp, -m, m - mp, mp)


def hund_determinant(config):
    """spin-orbitals (shell index, m, spin) of the highest-weight determinant of the Hund term: closed shells full;
    open shells: spin-up first (m = l, l-1, ...), then spin-down (m = l, l-1, ...)"""
    so = []
    for a, (n, l, q) in enumerate(config):
        ms = list(range(l, -l - 1, -1))
        up = min(q, 2 * l + 1); dn = q - up
        so += [(a, m, +1) for m in ms[:up]] + [(a, m, -1) for m in ms[:dn]]
    return so


def energy_coefficients(config, det=None):
    """E = sum_a q_a I_a + sum_{a<=b} sum_k [d^k_ab F^k(a,b) - e^k_ab G^k(a,b)], with pairs i,j counted as
    1/2 sum_{i,j} (self pairs included: they cancel). Returns d[a][b][k], e[a][b][k] for a<=b."""
    det = det if det is not None else hund_determinant(config)
    ns = len(config)
    d = {}; e = {}
    for i in det:
        for j in det:
            a, b = i[0], j[0]
            la, lb = config[a][1], config[b][1]
            key = (min(a, b), max(a, b))
            w = 0.5 if a == b else (1.0 if a < b else 0.0)   # a<b counted once (from the (i in a, j in b) ordering)
            if w == 0:
                continue
            for k in range(0, 2 * min(la, lb) + 1, 2):        # direct: c^k(la ma la ma) c^k(lb mb lb mb)
                v = ck(k, la, i[1], la, i[1]) * ck(k, lb, j[1], lb, j[1])
                if abs(v) > 1e-14:
                    d.setdefault(key, {}); d[key][k] = d[key].get(k, 0.0) + w * v
            if i[2] == j[2]:                                   # exchange
                for k in range(abs(la - lb), la + lb + 1, 2):
                    v = ck(k, la, i[1], lb, j[1]) ** 2
                    if abs(v) > 1e-14:
                        e.setdefault(key, {}); e[key][k] = e[key].get(k, 0.0) + w * v
    return d, e


# ------------------------------------------------------------------------------------------------ radial grid
FD8 = np.array([-1 / 560, 8 / 315, -1 / 5, 8 / 5, -205 / 72, 8 / 5, -1 / 5, 8 / 315, -1 / 560])


def geig_low(H, Sd, nlow, sigma):
    """lowest nlow eigenpairs of H u = lam diag(Sd) u, well conditioned: A = H - sigma S (SPD for sigma below the
    spectrum) = L L^T; C = L^-1 S L^-T; lam = sigma + 1/mu for the largest mu of C"""
    A = 0.5 * (H + H.T) - sigma * np.diag(Sd)
    L = np.linalg.cholesky(A)
    Li = sla.solve_triangular(L, np.eye(len(Sd)), lower=True)
    C = (Li * Sd[None, :]) @ Li.T
    M = len(Sd)
    mu, V = sla.eigh(0.5 * (C + C.T), subset_by_index=[M - nlow, M - 1])
    mu = mu[::-1]; V = V[:, ::-1]
    U = Li.T @ V
    return sigma + 1.0 / mu, U


class Grid:
    def __init__(self, Z, rmin_Z=1e-12, rmax=45.0, h=0.03):
        x0 = math.log(rmin_Z / Z); x1 = math.log(rmax)
        self.M = int(math.ceil((x1 - x0) / h)) + 1
        self.x = x0 + h * np.arange(self.M); self.h = h
        self.r = np.exp(self.x)
        M = self.M
        D2 = np.zeros((M, M))
        for o, c in zip(range(-4, 5), FD8):
            D2 += np.diag(np.full(M - abs(o), c), o)
        self.D2 = D2 / h ** 2
        self._L = {}

    def integ(self, f):
        """int f(r) dr = sum f r h (trapezoid in x)"""
        return float(np.sum(f * self.r) * self.h)

    def Lk(self, k):
        """dense map rho(r) (values on the grid) -> V(r) = int rho(r') r<^k / r>^(k+1) dr'.
        w = sqrt(r) V... solved from w'' - c^2 w = -(2k+1) sqrt(r) rho with c = k + 1/2, where w := r^(1/2) V(r)?
        (derivation: with U = r V, U'' - k(k+1) U / r^2 = -(2k+1) rho / r; U = sqrt(r) w gives the line above with
        w = U / sqrt(r) = sqrt(r) V.) Right boundary: ghost values from the multipole tail V = Q_k / r^(k+1)."""
        if k in self._L:
            return self._L[k]
        M, h, r = self.M, self.h, self.r
        c = k + 0.5
        A = self.D2 - c * c * np.eye(M)
        # RHS map: rho -> -(2k+1) sqrt(r) rho ; tail: Q = sum rho r^k r h ; ghost w_g = sqrt(r_g) Q / r_g^(k+1)
        R = np.diag(-(2 * k + 1) * np.sqrt(r))
        q = r ** k * r * h                                   # Q = q . rho
        G = np.zeros(M)                                      # contribution of ghost values to rows (moved to RHS)
        for g in range(1, 5):                                # ghost index M-1+g
            rg = r[-1] * math.exp(g * h)
            wg = math.sqrt(rg) / rg ** (k + 1)               # w_ghost = wg * Q
            for row in range(M - 4, M):
                off = (M - 1 + g) - row
                if off <= 4:
                    G[row] += FD8[4 + off] / h ** 2 * wg
        rhs = R - np.outer(G, q)
        W = np.linalg.solve(A, rhs)                          # w = W rho
        Lmat = W / np.sqrt(r)[:, None]                       # V = w / sqrt(r)
        self._L[k] = Lmat
        return Lmat


# ------------------------------------------------------------------------------------------------ atom
class Atom:
    def __init__(self, Z, config, mkin=1.0, grid=None, h=0.03, rmax=45.0, extV=None, det=None, nuc=None, extNL=None):
        """config: list of (n, l, occ). extV: optional dict l -> extra local radial potential (pseudo-atoms);
        nuc: optional array, the nuclear potential on the grid (default -Z/r)"""
        self.Z, self.config, self.mkin = Z, list(config), mkin
        self.g = grid if grid is not None else Grid(max(Z, 1), h=h, rmax=rmax)
        g = self.g
        self.d, self.e = energy_coefficients(self.config, det)
        self.q = np.array([c[2] for c in self.config], float)
        self.ls = sorted(set(c[1] for c in self.config))
        self.Vnuc = nuc if nuc is not None else -Z / g.r
        self.extV = extV or {}
        self.extNL = extNL or {}                         # l -> list of (chi(r) on the grid, E): sum |chi><chi| / E
        self.kmax = 2 * max(self.ls)
        for k in range(self.kmax + 1):
            g.Lk(k)
        self.P = None
        self.frozen = set()                              # shells kept fixed during scf (frozen-core tests)

    # one-electron operator in the u-basis: (1/mkin)(-1/2 D2 + 1/8) + r^2 (V + l(l+1)/(2 r^2)), overlap r^2
    def h_matrix(self, l):
        if not hasattr(self, "_hcache"):
            self._hcache = {}
        if l in self._hcache:
            return self._hcache[l]
        g = self.g; r = g.r
        T = (-0.5 * g.D2 + 0.125 * np.eye(g.M)) / self.mkin
        V = r * r * (self.Vnuc + self.extV.get(l, 0.0)) + 0.5 * l * (l + 1) / self.mkin
        H = T + np.diag(V)
        for chi, E in self.extNL.get(l, []):
            v = r ** 1.5 * chi
            H = H + np.outer(v, v) * g.h / E
        self._hcache[l] = H
        return H

    def J(self, b, k):
        return self.g.Lk(k) @ (self.P[b] ** 2)

    def K_matrix(self, b, k):
        """matrix (u-basis, already with the r^(3/2) factors) of phi -> P_b(r) V[P_b phi](r)"""
        g = self.g; r = g.r
        Lk = g.Lk(k)
        # phi(r) = sqrt(r) u ; output (K phi)(r) ; weak form row i: r_i^2 * (K phi)(r_i)/sqrt(r_i)
        Pb = self.P[b]
        Mtx = (r ** 1.5 * Pb)[:, None] * Lk * (Pb * np.sqrt(r))[None, :]
        return 0.5 * (Mtx + Mtx.T)

    def fock(self, a):
        """F_a in the u-basis (matrix), from the energy expression (see energy_coefficients)"""
        g = self.g; r = g.r
        la = self.config[a][1]
        F = self.h_matrix(la)
        loc = np.zeros(g.M)
        for (x, y), dk in self.d.items():
            if a not in (x, y):
                continue
            b = y if x == a else x
            f = 2.0 if x == y else 1.0
            for k, v in dk.items():
                loc += f * v / self.q[a] * self.J(b, k)
        F = F + np.diag(r * r * loc)
        for (x, y), ek in self.e.items():
            if a not in (x, y):
                continue
            b = y if x == a else x
            f = 2.0 if x == y else 1.0
            for k, v in ek.items():
                F = F - f * v / self.q[a] * self.K_matrix(b, k)
        return F

    def energy(self):
        g = self.g; r = g.r
        E1 = 0.0
        parts = {}
        for a, (n, l, q) in enumerate(self.config):
            Pa = self.P[a]
            u = Pa / np.sqrt(r)
            hu = self.h_matrix(l) @ u
            I = float(u @ hu) * g.h
            E1 += q * I
        E2 = 0.0
        for (x, y), dk in self.d.items():
            for k, v in dk.items():
                E2 += v * g.integ(self.P[x] ** 2 * self.J(y, k))
        for (x, y), ek in self.e.items():
            for k, v in ek.items():
                E2 -= v * g.integ(self.P[x] * self.P[y] * (g.Lk(k) @ (self.P[x] * self.P[y])))
        return E1 + E2, E1, E2

    def sigma(self):
        return -0.6 * self.Z ** 2 - 5.0

    def start(self):
        """bare-nucleus hydrogen-like orbitals (charge Z: no screening put by hand), from diagonalising h"""
        g = self.g; r = g.r
        S = r * r
        self.P = {}
        for l in self.ls:
            shells = [a for a, c in enumerate(self.config) if c[1] == l]
            H = self.h_matrix(l)
            w, U = geig_low(H, S, len(shells) + 2, self.sigma())
            for rank, a in enumerate(sorted(shells, key=lambda a: self.config[a][0])):
                u = U[:, rank]
                self.P[a] = self._norm(u * np.sqrt(r))
        return self

    def _norm(self, P):
        P = P / math.sqrt(self.g.integ(P * P))
        i = int(np.argmax(np.abs(P)))
        # sign: positive near the origin (first lobe)
        j = int(np.argmax(np.abs(P) > 1e-3 * abs(P[i])))
        return P if P[j] > 0 else -P

    def scf(self, tol=1e-9, maxit=300, mix=0.5, verbose=False):
        if self.P is None:
            self.start()
        g = self.g; r = g.r; S = np.diag(r * r)
        Eold = 0.0
        self.eps = {}
        hist = []
        for it in range(maxit):
            newP = {}
            for l in self.ls:
                shells = [a for a, c in enumerate(self.config) if c[1] == l]
                shells.sort(key=lambda a: self.config[a][0])
                closed = [a for a in shells if self.config[a][2] == 2 * (2 * l + 1)]
                opens = [a for a in shells if a not in closed]
                groups = ([closed] if closed else []) + [[a] for a in opens]
                for grp in groups:
                    F = self.fock(grp[0])
                    others = [b for b in shells if b not in grp]
                    if others:
                        # project out the other shells of the same l (orthogonality); u-basis with metric S
                        C = np.array([self.P[b] / np.sqrt(r) for b in others]).T      # M x m
                        SC = S @ C
                        Q = np.eye(g.M) - C @ np.linalg.solve(C.T @ SC * g.h, SC.T * g.h)
                        Fp = Q.T @ F @ Q
                        shift = 1e4
                        Fp = Fp + shift * (SC @ np.linalg.solve(C.T @ SC, SC.T))
                    else:
                        Fp = F
                    nlow = len(grp)
                    w, U = geig_low(Fp, S.diagonal(), nlow, self.sigma())
                    # the group takes the lowest eigenvectors of its (projected) Fock operator
                    for idx, a in enumerate(sorted(grp, key=lambda a: self.config[a][0])):
                        P = self._norm(U[:, idx] * np.sqrt(r))
                        newP[a] = P
                        self.eps[a] = float(w[idx])
            # mixing + re-orthonormalisation per l (Gram-Schmidt in n order: closed inner shells first)
            for a in list(newP):
                if a in self.frozen:
                    continue
                Pm = (1 - mix) * self.P[a] + mix * newP[a] if it > 0 else newP[a]
                self.P[a] = Pm
            self._orthonormalise()
            self._brillouin()
            E, E1, E2 = self.energy()
            hist.append(E)
            if verbose:
                print(it, f"{E:.10f}")
            if abs(E - Eold) < tol and it > 5:
                break
            Eold = E
        self.E = E
        self.hist = hist
        return dict(E=E, it=it, eps={self._name(a): self.eps[a] for a in self.eps})

    def _orthonormalise(self):
        for l in self.ls:
            shells = sorted([a for a, c in enumerate(self.config) if c[1] == l], key=lambda a: self.config[a][0])
            done = []
            for a in shells:
                P = self.P[a]
                for b in done:
                    P = P - self.g.integ(P * self.P[b]) * self.P[b]
                self.P[a] = self._norm(P)
                done.append(a)

    def _brillouin(self):
        """rotate pairs (closed, open) of the same l to the energy minimum (off-diagonal Lagrange condition)"""
        for l in self.ls:
            shells = [a for a, c in enumerate(self.config) if c[1] == l]
            closed = [a for a in shells if self.config[a][2] == 2 * (2 * l + 1)]
            opens = [a for a in shells if a not in closed]
            pairs = [(a, b) for a in closed for b in opens] + [(a, b) for i, a in enumerate(opens) for b in opens[i + 1:]]
            for a, b in pairs:
                if a in self.frozen or b in self.frozen:
                    continue
                Pa, Pb = self.P[a].copy(), self.P[b].copy()
                def Et(t):
                    c, s = math.cos(t), math.sin(t)
                    self.P[a] = c * Pa + s * Pb; self.P[b] = -s * Pa + c * Pb
                    return self.energy()[0]
                ts = np.array([-2e-3, 0.0, 2e-3])
                Es = np.array([Et(t) for t in ts])
                cur = (Es[2] - 2 * Es[1] + Es[0]) / (ts[2] ** 2)
                grad = (Es[2] - Es[0]) / (2 * ts[2])
                t = -grad / cur if cur > 0 else 0.0
                t = float(np.clip(t, -0.05, 0.05))
                if Et(t) > Es[1]:
                    Et(0.0)

    def _name(self, a):
        n, l, q = self.config[a]
        return f"{n}{L_LETTER[l]}"

    # ------------------------------------------------------------------ observables
    def rho(self):
        """total electron density (spherically averaged) rho(r), electrons / a0^3"""
        r = self.g.r
        return sum(q * self.P[a] ** 2 for a, (n, l, q) in enumerate(self.config)) / (4 * np.pi * r * r)

    def rho0(self):
        """density at the nucleus: only s shells, P/r -> R(0); extrapolated from the first grid points"""
        r = self.g.r
        tot = 0.0
        for a, (n, l, q) in enumerate(self.config):
            if l == 0:
                R = self.P[a] / r
                i = int(np.searchsorted(r, 1e-5 / self.Z))
                tot += q * R[i] ** 2 / (4 * np.pi)
        return tot

    def mean_r(self, a, p=1):
        return self.g.integ(self.P[a] ** 2 * self.g.r ** p)

    def kato(self):
        """-(1/(2 rho)) d rho/dr at r -> 0 (Kato cusp: = Z)"""
        r = self.g.r; rho = self.rho()
        i = int(np.searchsorted(r, 1e-4 / self.Z))
        dr = np.gradient(rho, r)
        return float(-dr[i] / (2 * rho[i]))


def config_from_string(s):
    """'1s2 2s2 2p3' -> [(1,0,2),(2,0,2),(2,1,3)]"""
    out = []
    for t in s.split():
        n = int(t[0]); l = L_LETTER.index(t[1]); q = int(t[2:])
        out.append((n, l, q))
    return out


if __name__ == "__main__":
    import time, sys
    for Z, c, ref in [(1, "1s1", -0.5), (2, "1s2", -2.861679996), (4, "1s2 2s2", -14.573023168),
                      (3, "1s2 2s1", -7.432726931), (6, "1s2 2s2 2p2", -37.688618963),
                      (7, "1s2 2s2 2p3", -54.400934209), (10, "1s2 2s2 2p6", -128.547098109)]:
        t0 = time.time()
        a = Atom(Z, config_from_string(c))
        res = a.scf()
        print(Z, c, f"{res['E']:.9f}", f"ref {ref}", f"diff {res['E'] - ref:.2e}", f"it {res['it']}",
              f"{time.time() - t0:.1f}s", {k: round(v, 6) for k, v in res["eps"].items()}, flush=True)

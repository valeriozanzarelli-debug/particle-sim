"""aggmol_frat_pot.py - study aggmol_frat_ ("test frattale con H2"), STEP 1: molecule -> molecule-ball.

Fine level: rigid H2 molecules with the full anisotropic H2-H2 force of aggmol phase 3 (rank-2 invariant expansion,
5 orientations H, X, T, L, S):
  V(R, u1, u2) = a0 + a1 [(u1.n)^2 + (u2.n)^2] + a2 (u1.u2)^2 + a3 (u1.n)(u2.n)(u1.u2) + a4 (u1.n)^2 (u2.n)^2
  H = a0 + a2, X = a0, T = a0 + a1, L = a0 + 2a1 + a2 + a3 + a4, S = a0 + a1 + a2 + a3/2 + a4/4.
Coarse level ("ball"): the spherical average V0 = a0 + 2a1/3 + a2/3 + (a3 + a4)/9.
What is lost, measured on dimers and trimers:
  - quantum rotors (para-H2, even J; B = 1/(m_p r^2) from the proton mass and r = 1.4 a0): ground energy of
    B (J1^2 + J2^2 [+ J3^2]) + sum of pair forces, against the ball (V0) - the physically right comparison at 20 K
    (the J = 2 level is 6B = 500 K above J = 0);
  - classical rotors at T (what a classical fine level would see): Boltzmann free energy over orientations;
  - three-body (Axilrod-Teller, C9 from the phase-3 pseudo-states), which the ball can carry.
Units: Eh, a0 inside; MD tables in K and Angstrom. Output lab/out/aggmol_frat_step1.json."""
import os, sys, json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import sph_harm_y
from scipy.linalg import eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
EH_K = 315775.02480407
A0_A = 0.529177210903
MP = 1836.15267343
ME_AMU = 5.48579909065e-4
M_H2_AMU = (2 * MP + 2) * ME_AMU            # molecule mass (protons + electrons), from constants only
RB = 1.4
BROT = 1.0 / (MP * RB ** 2)                 # rigid rotor constant 1/(2I), I = m_p r^2 / 2
RCUT = 24.0                                 # a0, pair cutoff (energy shifted)


def coefficients():
    """a_k(R) on the phase-3 R grid from the 5 orientations"""
    g = json.load(open(os.path.join(OUT, "aggmol_h2h2.json")))["geom"]
    g.update(json.load(open(os.path.join(OUT, "aggmol_frat_aniso.json")))["geom"])
    Rs = np.array(sorted({v["R"] for v in g.values()}))
    O = {o: np.array([g[f"{o}_{R:.1f}"]["E_int"] for R in Rs]) for o in "HXTLS"}
    a0 = O["X"]; a2 = O["H"] - O["X"]; a1 = O["T"] - O["X"]
    p = O["L"] - a0 - 2 * a1 - a2
    q = O["S"] - a0 - a1 - a2
    a3 = 4 * q - p; a4 = 2 * p - 4 * q
    return Rs, np.stack([a0, a1, a2, a3, a4])


class Pair:
    """a_k(R) for any R: spline on [5, 12]; below 5 exponential continuation (common decay of the isotropic wall);
    above 12: isotropic -C6/R^6 - C8/R^8 (tail fit of phase 3), anisotropic parts ~ R^-5 (quadrupole-quadrupole)."""

    def __init__(self):
        self.Rs, self.a = coefficients()
        self.sp = [CubicSpline(self.Rs, ak) for ak in self.a]
        t = json.load(open(os.path.join(OUT, "aggmol_h2h2.json")))["tail"]
        self.C6, self.C8 = t["C6_fit"], t["C8_fit"]
        v0 = self.iso_grid()
        self.b = np.log(v0[0] / v0[1]) / (self.Rs[1] - self.Rs[0])
        self.C9 = json.load(open(os.path.join(OUT, "aggmol_frat_aniso.json")))["C9"]

    def iso_grid(self):
        a0, a1, a2, a3, a4 = self.a
        return a0 + 2 * a1 / 3 + a2 / 3 + (a3 + a4) / 9

    def coef(self, R):
        R = np.atleast_1d(np.asarray(R, float))
        out = np.zeros((5, len(R)))
        lo, hi = R < self.Rs[0], R > self.Rs[-1]
        mid = ~(lo | hi)
        for k in range(5):
            out[k, mid] = self.sp[k](R[mid])
            out[k, lo] = self.sp[k](self.Rs[0]) * np.exp(-self.b * (R[lo] - self.Rs[0]))
            out[k, hi] = self.sp[k](self.Rs[-1]) * (self.Rs[-1] / R[hi]) ** 5
        # isotropic tail: a0 + 2a1/3 + a2/3 + (a3+a4)/9 = -C6/R^6 - C8/R^8 beyond 12 (keep a1..a4 as R^-5, fix a0)
        if hi.any():
            rest = 2 * out[1, hi] / 3 + out[2, hi] / 3 + (out[3, hi] + out[4, hi]) / 9
            v12 = -self.C6 / self.Rs[-1] ** 6 - self.C8 / self.Rs[-1] ** 8
            iso12 = float(self.iso_grid()[-1])
            out[0, hi] = (-self.C6 / R[hi] ** 6 - self.C8 / R[hi] ** 8) * (iso12 / v12) - rest
        return out

    def V0(self, R):
        a = self.coef(R)
        return a[0] + 2 * a[1] / 3 + a[2] / 3 + (a[3] + a[4]) / 9

    def V(self, R, u1, u2, n):
        """full anisotropic pair energy; R (...), unit vectors (..., 3)"""
        a = self.coef(np.ravel(R)).reshape((5,) + np.shape(R))
        c1 = np.sum(u1 * n, -1); c2 = np.sum(u2 * n, -1); c12 = np.sum(u1 * u2, -1)
        return a[0] + a[1] * (c1 ** 2 + c2 ** 2) + a[2] * c12 ** 2 + a[3] * c1 * c2 * c12 + a[4] * c1 ** 2 * c2 ** 2


# ------------------------------------------------------------------------------------------- quantum rotors
def rotor_basis(Jmax, nth=28, nph=56):
    """real spherical harmonics, even J <= Jmax, on a Gauss-Legendre x uniform grid; returns Y (nb, ng), w, u (ng, 3), J"""
    x, w = np.polynomial.legendre.leggauss(nth)
    th = np.arccos(x); ph = 2 * np.pi * np.arange(nph) / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    W = (w[:, None] * np.ones(nph)[None] * 2 * np.pi / nph).ravel()
    u = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1).reshape(-1, 3)
    Ys, Js = [], []
    for J in range(0, Jmax + 1, 2):
        for m in range(-J, J + 1):
            y = sph_harm_y(J, abs(m), TH, PH).ravel()
            if m < 0:
                y = np.sqrt(2) * (-1) ** m * y.imag
            elif m > 0:
                y = np.sqrt(2) * (-1) ** m * y.real
            else:
                y = y.real
            Ys.append(y); Js.append(J)
    Y = np.array(Ys)
    return Y, W, u, np.array(Js)


def one_body_ops(Jmax):
    Y, W, u, J = rotor_basis(Jmax)
    S = (Y * W) @ Y.T
    assert np.allclose(S, np.eye(len(J)), atol=1e-10), "rotor basis not orthonormal"
    O = np.array([[(Y * W * u[:, i] * u[:, j]) @ Y.T for j in range(3)] for i in range(3)])   # <Y|u_i u_j|Y'>
    return O, J


def pair_operator(a, n, O):
    """V for one pair as a list of (coef, op_on_1, op_on_2); n unit vector 1 -> 2"""
    nb = O.shape[-1]
    I = np.eye(nb)
    N = np.einsum("i,ijab,j->ab", n, O, n)
    P = np.einsum("k,kiab->iab", n, O)             # (u.n) u_i
    terms = [(a[0], I, I), (a[1], N, I), (a[1], I, N), (a[4], N, N)]
    for i in range(3):
        terms.append((a[3], P[i], P[i]))
        for j in range(3):
            terms.append((a[2], O[i, j], O[i, j]))
    return terms


def dimer_quantum(pair, R, Jmax=6):
    O, J = one_body_ops(Jmax)
    nb = len(J)
    a = pair.coef(R)[:, 0]
    Hm = BROT * (np.kron(np.diag(J * (J + 1.0)), np.eye(nb)) + np.kron(np.eye(nb), np.diag(J * (J + 1.0))))
    for c, A, B in pair_operator(a, np.array([0, 0, 1.0]), O):
        Hm = Hm + c * np.kron(A, B)
    return float(eigh(Hm, eigvals_only=True, subset_by_index=[0, 0])[0])


def trimer_quantum(pair, R, Jmax=4):
    """equilateral trimer of side R, quantum rotors"""
    O, J = one_body_ops(Jmax)
    nb = len(J); I = np.eye(nb)
    pos = np.array([[0, 0, 0], [R, 0, 0], [R / 2, R * np.sqrt(3) / 2, 0]])
    Jd = np.diag(J * (J + 1.0))
    Hm = BROT * (np.kron(np.kron(Jd, I), I) + np.kron(np.kron(I, Jd), I) + np.kron(np.kron(I, I), Jd))
    a = pair.coef(R)[:, 0]
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        n = pos[j] - pos[i]; n /= np.linalg.norm(n)
        for c, A, B in pair_operator(a, n, O):
            ops = [I, I, I]; ops[i] = A; ops[j] = B
            Hm = Hm + c * np.kron(np.kron(ops[0], ops[1]), ops[2])
    return float(eigh(Hm, eigvals_only=True, subset_by_index=[0, 0])[0])


# ------------------------------------------------------------------------------------------- classical rotors
def sphere_grid(nth=24, nph=48):
    x, w = np.polynomial.legendre.leggauss(nth)
    th = np.arccos(x); ph = 2 * np.pi * np.arange(nph) / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    W = (w[:, None] * np.ones(nph)[None] * 2 * np.pi / nph).ravel() / (4 * np.pi)
    u = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1).reshape(-1, 3)
    return u, W


def dimer_classical(pair, R, T_K):
    u, W = sphere_grid()
    n = np.array([0, 0, 1.0])
    V = pair.V(np.full((len(u), len(u)), R), u[:, None], u[None], n)
    beta = EH_K / T_K
    Vm = V.min()
    Z = np.sum(W[:, None] * W[None] * np.exp(-beta * (V - Vm)))
    return dict(F=float(Vm - np.log(Z) / beta), Vmin=float(Vm),
                Emean=float(np.sum(W[:, None] * W[None] * V * np.exp(-beta * (V - Vm))) / Z))


def trimer_classical(pair, R, T_K, nsamp=400000, seed=5):
    rng = np.random.default_rng(seed)
    pos = np.array([[0, 0, 0], [R, 0, 0], [R / 2, R * np.sqrt(3) / 2, 0]])
    u = rng.normal(size=(nsamp, 3, 3)); u /= np.linalg.norm(u, axis=-1, keepdims=True)
    E = 0.0
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        n = pos[j] - pos[i]; n /= np.linalg.norm(n)
        E = E + pair.V(np.full(nsamp, R), u[:, i], u[:, j], n)
    beta = EH_K / T_K
    Em = E.min()
    Z = np.mean(np.exp(-beta * (E - Em)))
    # minimum: local polish of the best samples
    from scipy.optimize import minimize
    def f(x):
        th, ph = x[:3], x[3:]
        uu = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], -1)
        e = 0.0
        for (i, j) in [(0, 1), (0, 2), (1, 2)]:
            n = pos[j] - pos[i]; n /= np.linalg.norm(n)
            e += float(pair.V(np.array([R]), uu[i][None], uu[j][None], n)[0])
        return e
    best = np.inf
    for k in np.argsort(E)[:8]:
        th = np.arccos(np.clip(u[k, :, 2], -1, 1)); ph = np.arctan2(u[k, :, 1], u[k, :, 0])
        r = minimize(f, np.concatenate([th, ph]), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-12, maxiter=4000))
        best = min(best, r.fun)
    return dict(F=float(Em - np.log(Z) / beta), Emin=float(best))


def atm(C9, r12, r13, r23):
    """Axilrod-Teller-Muto triple-dipole energy (isotropic C9)"""
    c1 = (r12 ** 2 + r13 ** 2 - r23 ** 2) / (2 * r12 * r13)
    c2 = (r12 ** 2 + r23 ** 2 - r13 ** 2) / (2 * r12 * r23)
    c3 = (r13 ** 2 + r23 ** 2 - r12 ** 2) / (2 * r13 * r23)
    return C9 * (1 + 3 * c1 * c2 * c3) / (r12 * r13 * r23) ** 3


def job_step1():
    pair = Pair()
    res = dict(units="Eh, a0", B_rot=BROT, B_rot_K=BROT * EH_K, m_H2_amu=M_H2_AMU, C9=pair.C9, C6_tail=pair.C6,
               C8_tail=pair.C8, wall_decay=pair.b)
    res["coef_grid"] = dict(R=pair.Rs.tolist(), a=pair.a.tolist())
    x = np.linspace(5.0, 12.0, 7001)
    v0 = pair.V0(x); i = int(np.argmin(v0))
    res["ball_well"] = dict(R=float(x[i]), V=float(v0[i]), V_K=float(v0[i] * EH_K))
    # dimers
    Rd = [5.5, 6.0, 6.25, 6.5, 6.75, 7.0, 7.5, 8.0, 9.0, 10.0]
    dim = []
    for R in Rd:
        q4 = dimer_quantum(pair, R, 4); q6 = dimer_quantum(pair, R, 6)
        row = dict(R=R, V0=float(pair.V0(R)[0]), quantum_J4=q4, quantum_J6=q6)
        for T in (10.0, 20.0, 50.0):
            row[f"classical_{int(T)}K"] = dimer_classical(pair, R, T)
        dim.append(row)
        print(f"dimer R={R}: V0 {row['V0']*EH_K:.3f} K, quantum {q6*EH_K:.3f} K (J4 {q4*EH_K:.3f}), "
              f"classical F(20K) {row['classical_20K']['F']*EH_K:.3f} K, Vmin {row['classical_20K']['Vmin']*EH_K:.3f} K",
              flush=True)
    res["dimer"] = dim
    # trimers (equilateral)
    tri = []
    for R in [6.0, 6.5, 7.0, 7.5]:
        q = trimer_quantum(pair, R, 4)
        c = trimer_classical(pair, R, 20.0)
        row = dict(R=R, ball=float(3 * pair.V0(R)[0]), quantum_J4=q, classical_20K_F=c["F"], classical_min=c["Emin"],
                   ATM=float(atm(pair.C9, R, R, R)))
        tri.append(row)
        print(f"trimer R={R}: ball {row['ball']*EH_K:.3f} K, quantum {q*EH_K:.3f} K, classical F20 {c['F']*EH_K:.3f} K,"
              f" min {c['Emin']*EH_K:.3f} K, ATM {row['ATM']*EH_K:.3f} K", flush=True)
    res["trimer"] = tri
    # second virial coefficient (classical translation): ball vs classical anisotropic rotors
    R = np.linspace(3.0, 40.0, 1500)
    bv = {}
    u, W = sphere_grid(12, 24)
    for T in (20.0, 50.0, 100.0, 300.0):
        beta = EH_K / T
        f0 = np.exp(-beta * pair.V0(R)) - 1
        fa = []
        for r in R:
            V = pair.V(np.full((len(u), len(u)), r), u[:, None], u[None], np.array([0, 0, 1.0]))
            fa.append(np.sum(W[:, None] * W[None] * np.exp(-beta * V)) - 1)
        fa = np.array(fa)
        conv = (A0_A * 1e-8) ** 3 * 6.02214076e23                     # a0^3 per molecule -> cm^3/mol
        bv[f"{int(T)}K"] = dict(ball=float(-2 * np.pi * np.trapezoid(f0 * R ** 2, R) * conv),
                                aniso=float(-2 * np.pi * np.trapezoid(fa * R ** 2, R) * conv))
        print("B2", T, bv[f"{int(T)}K"], flush=True)
    res["B2_cm3mol"] = bv
    json.dump(res, open(os.path.join(OUT, "aggmol_frat_step1.json"), "w"), indent=1)


def table_K_A(pair=None, rmin_A=2.2, rmax_A=None, n=20001):
    """isotropic ball potential and its derivative as tables in K and Angstrom (shifted to 0 at the cutoff)"""
    pair = pair or Pair()
    rmax_A = rmax_A or RCUT * A0_A
    r = np.linspace(rmin_A, rmax_A, n)
    v = pair.V0(r / A0_A) * EH_K
    v = v - v[-1]
    dv = np.gradient(v, r)
    return r, v, dv


if __name__ == "__main__":
    job_step1()

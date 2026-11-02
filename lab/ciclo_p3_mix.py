"""ciclo_p3_mix.py - Phase 3 of CICLO_LIVELLI: mixed coarse / base representation with SMOOTHED automatic variables,
the on-the-fly check and local refinement (step 4 of the loop).

Ingredients (all from phase 2, read-only, plus the smoothing found in development, ciclo_p3_dev.py):
  base     : 2187 cells of 1.6 lambda_C, kinetic stencil of the fine zig-zag (phase 2 base), soft-Coulomb atom
  variables: the AUTOMATIC block vectors of the phase-2 b = 3 chain (2187 -> 729 -> 243), each level SMOOTHED s = 3
             times by its own rule (smoothed aggregation, step 4/(3 rho)), Loewdin -> prolongator P (2187 x 243);
             the coarse rule is P^dag H P (local; agrees with the phase-2 LEARNED rule within 0.3-1.3 %)
  mixed    : coarse cell K -> one number on its smooth function P_K masked off the refined base sites (Loewdin among
             the coarse ones); refined cell -> its 9 base cells.  S (2187 x m) orthonormal; H_mix = S^dag H S
             (variational: no ghost terms at the seam); Crank-Nicolson steps = exactly unitary
  check    : every NCHK steps, each patch of 3 coarse cells holding weight: run the BASE rule on patch +- BUF base cells
             from the current state for the same interval; e = || psi_base - psi_mix ||_patch.  e > EPS_REF ->
             refine the patch (+ MARGIN coarse cells per side), roll back to the start of the interval, redo
  truth    : the base rule everywhere (2187 cells), same integrator, same start
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, pickle, time
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import scipy.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ciclo_p3_core as p3
import ciclo_p3_dev as dev

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
N0 = p3.N0
S_SMOOTH = 3
SV_TOL = 1e-8          # drop masked-tail directions whose weight is below this (numerical floor)
PROL_CACHE = os.path.join(OUT, "ciclo_p3_prol.pkl")      # not committed


def prolongators():
    """P1 (2187 x 729) and P (2187 x 243): composite smoothed aggregation, s = 3 per level"""
    if os.path.exists(PROL_CACHE):
        return pickle.load(open(PROL_CACHE, "rb"))
    hy = p3.hierarchy(2)
    Hb = hy["Hb"].toarray()
    W1 = hy["levels"][0]["W"].toarray(); W2 = hy["levels"][1]["W"].toarray()
    P1 = dev.SA(W1, Hb, S_SMOOTH)
    G1 = P1.conj().T @ Hb @ P1
    P = P1 @ dev.SA(W2 @ W1.conj().T, G1, S_SMOOTH)
    out = dict(P1=P1, P=P, Hb=hy["Hb"], V0=hy["V0"])
    pickle.dump(out, open(PROL_CACHE, "wb"))
    return out


class Mixed:
    def __init__(self, P, Hs, refined):
        """P: (N0, nc) smooth coarse functions; Hs: sparse static base generator (kinetic + V0)"""
        self.P, self.Hs = P, Hs
        self.nc = P.shape[1]
        self.B = N0 // self.nc
        self.ref = refined.copy()
        self.build()

    def build(self):
        fine = np.zeros(N0, bool)
        for I in np.where(self.ref)[0]:
            fine[I * self.B:(I + 1) * self.B] = True
        self.fine_sites = np.where(fine)[0]
        # the smooth functions of ALL cells, masked off the refined sites (a refined cell keeps its tail outside), then
        # canonical orthogonalisation (drop directions with overlap eigenvalue < SV_TOL): old span inside new span,
        # so refining is lossless up to the dropped directions
        Phi = self.P.copy()
        Phi[fine] = 0
        Sv = Phi.conj().T @ Phi
        e, U = np.linalg.eigh(Sv)
        k = e > SV_TOL
        Phi = Phi @ (U[:, k] / np.sqrt(e[k]))
        self.ncoarse = int(k.sum())
        E = np.zeros((N0, len(self.fine_sites)))
        E[self.fine_sites, np.arange(len(self.fine_sites))] = 1
        self.S = np.concatenate([Phi, E], 1)
        self.m = self.S.shape[1]
        self.A0 = self.S.conj().T @ (self.Hs @ self.S)

    def set_drive(self, Vd):
        self.Ad = None if Vd is None else self.S.conj().T @ (Vd[:, None] * self.S)

    def H(self, f):
        return self.A0 if (self.Ad is None or f == 0) else self.A0 + f * self.Ad


def cn(Hd, a, dt):
    n = Hd.shape[0]
    M = np.eye(n) + 0.5j * dt * Hd
    return sla.solve(M, a - 0.5j * dt * (Hd @ a))


class Scenario:
    def __init__(self, amp, xd, sigma, omega, Tp, T, dt):
        self.amp, self.xd, self.sigma, self.omega, self.Tp, self.T, self.dt = amp, xd, sigma, omega, Tp, T, dt
        x = np.arange(N0) + 0.5
        self.Vd = None if amp == 0 else amp * np.exp(-((x - xd) ** 2) / (2 * sigma ** 2))
        self.nsteps = int(round(T / dt))

    def f(self, t):
        """drive factor at time t (midpoint rule)"""
        if self.Vd is None or t > self.Tp:
            return 0.0
        return np.sin(np.pi * t / self.Tp) ** 2 * np.sin(self.omega * t)


def truth(sc, Hs, psi0, save_every):
    """base everywhere, Crank-Nicolson with sparse LU (factor reused while the drive is off)"""
    I = sps.identity(N0, format="csc", dtype=complex)
    psi = psi0.astype(complex)
    snaps = [psi.copy()]
    lu_static = None
    t0 = time.time()
    for n in range(sc.nsteps):
        t = (n + 0.5) * sc.dt
        fv = sc.f(t)
        H = Hs if fv == 0 else (Hs + sps.diags(fv * sc.Vd)).tocsc()
        rhs = psi - 0.5j * sc.dt * (H @ psi)
        if fv == 0:
            if lu_static is None:
                lu_static = spla.splu((I + 0.5j * sc.dt * Hs).tocsc())
            psi = lu_static.solve(rhs)
        else:
            psi = spla.splu((I + 0.5j * sc.dt * H).tocsc()).solve(rhs)
        if (n + 1) % save_every == 0:
            snaps.append(psi.copy())
    return np.array(snaps), time.time() - t0


def run_mixed(sc, P, Hs, a0_base, nchk, eps_ref, buf, margin, patch=3, check=True, maxit=6, log=print):
    """the adaptive run; returns dict with snapshots (base, every nchk steps), refinement history, costs, energies"""
    nc = P.shape[1]; B = N0 // nc
    # a patch whose window holds weight w cannot fail: unitary evolutions keep ||y||, ||psi_mix|| <= sqrt(w) on the
    # window, so e <= 2 sqrt(w) up to what flows in through the buffer; skip if 2 sqrt(w) < eps_ref (no knob)
    wmin = (eps_ref / 2) ** 2
    mx = Mixed(P, Hs, np.zeros(nc, bool))
    mx.set_drive(sc.Vd)
    a = mx.S.conj().T @ a0_base
    Hd_static = Hs.toarray() if check else None
    snaps = [mx.S @ a]
    hist, events, cost_steps, cost_check = [], [], 0, 0
    npatch = nc // patch
    norm_loss = 0.0
    t_start = time.time()
    n = 0
    while n < sc.nsteps:
        nn = min(nchk, sc.nsteps - n)
        emax = 0.0
        for it in range(maxit):
            b = a.copy()
            for q in range(nn):
                t = (n + q + 0.5) * sc.dt
                b = cn(mx.H(sc.f(t)), b, sc.dt)
            cost_steps += nn * mx.m
            if not check:
                break
            psi_s = mx.S @ a
            psi_e = mx.S @ b
            flags, emax = [], 0.0
            for pp in range(npatch):
                cells = np.arange(pp * patch, (pp + 1) * patch)
                if mx.ref[cells].all():
                    continue
                lo, hi = pp * patch * B, (pp + 1) * patch * B
                win = np.arange(lo - buf, hi + buf) % N0
                w0 = psi_s[win]
                if np.sum(np.abs(w0) ** 2) < wmin:
                    continue
                Hw = Hd_static[np.ix_(win, win)]
                Vw = None if sc.Vd is None else sc.Vd[win]
                y = w0
                for q in range(nn):
                    t = (n + q + 0.5) * sc.dt
                    fv = sc.f(t)
                    y = cn(Hw if fv == 0 else Hw + np.diag(fv * Vw), y, sc.dt)
                cost_check += nn * len(win)
                e = float(np.linalg.norm(y[buf:buf + hi - lo] - psi_e[lo:hi]))
                emax = max(emax, e)
                if e > eps_ref:
                    flags.append(pp)
            if not flags:
                break
            new = np.zeros(nc, bool)
            for pp in flags:
                new[max(0, pp * patch - margin):min(nc, (pp + 1) * patch + margin)] = True
            new &= ~mx.ref
            if not new.any():
                break
            psi = mx.S @ a
            mx.ref |= new
            mx.build(); mx.set_drive(sc.Vd)
            a = mx.S.conj().T @ psi
            lost = float(np.linalg.norm(psi) ** 2 - np.linalg.norm(a) ** 2)
            norm_loss += lost
            events.append(dict(step=n, t=n * sc.dt, cells=np.where(new)[0].tolist(), lost=lost, emax=emax, m=mx.m))
            log(f"    t={n * sc.dt:8.0f} refine {int(new.sum())} cells (m = {mx.m}), e_max {emax:.2e}, lost {lost:.1e}")
        a = b
        n += nn
        snaps.append(mx.S @ a)
        E_own = float(np.real(np.vdot(a, mx.H(0.0) @ a)))
        hist.append(dict(t=n * sc.dt, m=mx.m, emax=emax, nref=int(mx.ref.sum()), norm=float(np.linalg.norm(a)), E_own=E_own))
    return dict(snaps=np.array(snaps), hist=hist, events=events, cost_steps=cost_steps, cost_check=cost_check,
                norm_loss=norm_loss, refined=mx.ref.tolist(), wall=time.time() - t_start)

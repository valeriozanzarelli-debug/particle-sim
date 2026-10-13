"""surf_core.py - "Surfisti esatti": exactly equivariant surfers for the zig-zag walk (study SURFISTI_ESATTI).

The zig-zag rule (zz_core.py, frozen) moves probability ONLY in the translation sub-steps, and there it moves it
across links in a perfectly known way: along axis j, the part of the spinor heading + (eigenvalue +1 of the
translation) at site x crosses the link (x, x+e_j), the part heading - at x+e_j crosses it backwards.
So the walk has EXACT link currents:
        gross:  A_j(x) = mass at x heading +,  B_j(x+e_j) = mass at x+e_j heading -
        net:    F_j(x) = A_j(x) - B_j(x+e_j)                    (flux from cell x to cell x+e_j in one sub-step)
and   rho_after = rho_before - div_j F_j   exactly (the coin acts inside a site and does not move |psi|^2).

Surfer rules implemented here (all in LATTICE units: cell = 1, tick = 1):
  S0 'vecchio'   : the old reader (Rust zz.rs::muovi_surfisti): v = c psi^dag alpha psi / psi^dag psi from a
                   trilinear interpolation of the spinor at the tick boundary, midpoint rule.   NOT equivariant.
  S1 'balistico' : Bell jump process with the GROSS flows of every sub-step (lambda = 1 of pre_core): a surfer at
                   x jumps + with prob A/rho, - with prob B/rho. Exactly equivariant, maximally random.
  S2 'bell'      : Bell minimal jump process per sub-step (lambda = 0 of pre_core, Vink 1993): only the NET flow
                   is carried, prob max(F,0)/rho. Exactly equivariant, minimal randomness per sub-step.
  S2b 'bellflux' : Bell continuous-time jump process over the whole tick in FLOW TIME (net flux of the tick,
                   F = sum of the sub-step fluxes, rho(s) linear in s inside the tick, rates F+/rho(s)).
                   Exactly equivariant (master equation solved exactly), randomness only where the tick net
                   current is non-zero.
  S3 'cdfasse'   : deterministic monotone CDF map along each line for every sub-step (IDEE_RIPRESE R3).
                   Exactly equivariant for the piecewise-uniform measure.
  S4 'flusso'    : NEW. Deterministic flow in "flow time" (Pollock 1988 semi-analytic particle tracking) of the
                   net flux of the tick: inside each cell the flux is the lowest-order Raviart-Thomas field (each
                   component linear in its own coordinate, matching the face fluxes), rho(s) = rho0 + s (rho1-rho0)
                   uniform in the cell. Continuity holds exactly in every cell and normal flux is continuous across
                   faces, so the flow pushes the piecewise-uniform measure rho_t EXACTLY onto rho_{t+1}. With the
                   flow time tau = int ds / rho(s) every coordinate moves as q_j(tau) = q_j0 exp(d_j tau): solved
                   analytically. Deterministic, reversible (run it with -F from rho1 to rho0), zero random numbers
                   after the initial Born sampling. In 1-D it coincides with S3 (the monotone map is unique).
Positions: integer cell index + local coordinate xi in [0,1]^d (cell centre = 0.5, as in Rust x_of()).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import numba as nb
from zz_core import BASIS_HALF, order_of_tick, RadialWalk, ALPHA

SEED = 20261013
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


# ======================================================================================================
#  WALKS THAT REPORT THEIR LINK CURRENTS
# ======================================================================================================
class Walk1D:
    """1-D zig-zag (zz_core 'dirac1'): psi (2, n) = (R, L); R moves +1, L moves -1, then the flip, then the phase.
    periodic=True: ring. tick() returns (psi, [(0, A, B)]) with A = |R|^2 (heading +), B = |L|^2 (heading -)."""

    def __init__(self, n, theta, V=None, periodic=True):
        self.n, self.theta, self.periodic = n, theta, periodic
        self.phase = None if V is None else np.exp(-1j * V)      # V in units 1/tick
        self.tick_count = 0

    def tick(self, psi):
        R, L = psi
        A, B = np.abs(R)**2, np.abs(L)**2
        R2, L2 = np.roll(R, 1), np.roll(L, -1)
        ct, st = np.cos(self.theta), np.sin(self.theta)
        out = np.array([ct * R2 - 1j * st * L2, ct * L2 - 1j * st * R2])
        if self.phase is not None:
            out *= self.phase
        self.tick_count += 1
        return out, [(0, A, B)]


class Radial:
    """The radial form of the rule (zz_core.RadialWalk) with its link currents. chi+ (index 0) moves outwards,
    chi- inwards; at the first site chi- turns into chi+ (no link crossed), at the last site the reverse.
    Link i = (i, i+1), i = 0..N-2; A = |chi+|^2, B = |chi-|^2."""

    def __init__(self, rw):
        self.rw = rw
        self.tick_count = 0

    def tick(self, chi):
        A, B = np.abs(chi[0])**2, np.abs(chi[1])**2
        out = self.rw.tick(chi)
        self.tick_count += 1
        return out, [(0, A, B)]


class Walk3D:
    """3-D zig-zag (zz_core.Walk 'dirac', dim 3) on a periodic n^3 box; psi (4, n, n, n) = (L0, L1, R0, R1).
    V (n,n,n) in units of 1/tick (phase e^{-iV} per tick, as zz_core with V*a). tick() returns
    (psi, [(j, A_j, B_j) for the three sub-steps in the order they are applied])."""

    def __init__(self, n, theta, V=None):
        self.n, self.theta = n, theta
        self.phase = None if V is None else np.exp(-1j * V)
        self.tick_count = 0
        self.P = [BASIS_HALF[j][1] for j in range(3)]

    def tick(self, psi):
        subs = []
        L, R = psi[0:2], psi[2:4]
        for j in order_of_tick(self.tick_count):
            P = self.P[j]
            Pd = P.conj().T
            cL = np.tensordot(Pd, L, axes=(1, 0))
            cR = np.tensordot(Pd, R, axes=(1, 0))
            # BASIS_HALF order: q = 0 eigenvalue +1, q = 1 eigenvalue -1. L moves -lam, R moves +lam.
            A = np.abs(cL[1])**2 + np.abs(cR[0])**2
            B = np.abs(cL[0])**2 + np.abs(cR[1])**2
            subs.append((j, A, B))
            cL = np.stack([np.roll(cL[0], -1, axis=j), np.roll(cL[1], 1, axis=j)])
            cR = np.stack([np.roll(cR[0], 1, axis=j), np.roll(cR[1], -1, axis=j)])
            L = np.tensordot(P, cL, axes=(1, 0))
            R = np.tensordot(P, cR, axes=(1, 0))
        ct, st = np.cos(self.theta), np.sin(self.theta)
        out = np.concatenate([ct * L - 1j * st * R, ct * R - 1j * st * L])
        if self.phase is not None:
            out *= self.phase
        self.tick_count += 1
        return out, subs


def net_flux(j, A, B, periodic=True):
    """net flux across link (x, x+e_j): A(x) - B(x+e_j). Non-periodic (radial / walls): last link = 0."""
    F = A - np.roll(B, -1, axis=j)
    if not periodic:
        idx = [slice(None)] * A.ndim
        idx[j] = -1
        F[tuple(idx)] = 0.0
    return F


def density(psi):
    return np.sum(np.abs(psi)**2, axis=0)


# ======================================================================================================
#  S4: DETERMINISTIC FLOW IN FLOW TIME (Pollock / RT0), numba, any d = 1, 2, 3
# ======================================================================================================
@nb.njit(cache=True)
def _phi(x):
    """log1p(x)/x, -> 1 at 0; +inf at x <= -1"""
    if x <= -1.0:
        return np.inf
    if abs(x) < 1e-8:
        return 1.0 - 0.5 * x + x * x / 3.0
    return np.log1p(x) / x


@nb.njit(cache=True)
def _E(x):
    """expm1(x)/x, -> 1 at 0"""
    if abs(x) < 1e-8:
        return 1.0 + 0.5 * x + x * x / 6.0
    return np.expm1(x) / x


@nb.njit(cache=True)
def _flat(idx, shape):
    k = 0
    for a in range(idx.shape[0]):
        k = k * shape[a] + idx[a]
    return k


@nb.njit(cache=True)
def pollock_move(cell, xi, rho0, rho1, F, shape, periodic, maxcross):
    """move surfers (cell (N,d) int64, xi (N,d) float64) through one tick.
    rho0, rho1: flattened densities (C order of `shape`); F: (d, ncell) net flux of the tick, F[j, c] = flux from
    cell c to c+e_j. periodic: wrap; else a link out of the box has (must have) zero flux.
    Returns (number of face crossings per surfer, number of surfers that hit maxcross)."""
    N, d = cell.shape
    ncross = np.zeros(N, np.int64)
    stuck = 0
    nbr = np.zeros(d, np.int64)
    for p in range(N):
        s = 0.0
        for it in range(maxcross + 1):
            c = _flat(cell[p], shape)
            r0, r1 = rho0[c], rho1[c]
            D = r1 - r0
            rs = r0 + s * D
            if rs <= 0.0:
                rs = 1e-300
            # flow time left in the tick
            tau_end = (1.0 - s) / rs * _phi(D * (1.0 - s) / rs)
            tau_best = tau_end
            jbest, side = -1, 0
            for j in range(d):
                b = F[j, c]
                for a_ in range(d):
                    nbr[a_] = cell[p, a_]
                nbr[j] -= 1
                if nbr[j] < 0:
                    if periodic:
                        nbr[j] += shape[j]
                        a = F[j, _flat(nbr, shape)]
                    else:
                        a = 0.0
                else:
                    a = F[j, _flat(nbr, shape)]
                dj = b - a
                x0 = xi[p, j]
                q0 = a + dj * x0
                if q0 > 0.0 and b > 0.0:
                    t = (1.0 - x0) / q0 * _phi(dj * (1.0 - x0) / q0)
                    if t < tau_best:
                        tau_best, jbest, side = t, j, 1
                elif q0 < 0.0 and a < 0.0:
                    t = (-x0 / q0) * _phi(-dj * x0 / q0)
                    if t < tau_best:
                        tau_best, jbest, side = t, j, -1
            # advance every coordinate by tau_best
            if tau_best == np.inf:
                tau_best = 0.0     # (cannot happen for rho1 > 0; rho1 = 0 cells are emptied by an exit)
            for j in range(d):
                b = F[j, c]
                for a_ in range(d):
                    nbr[a_] = cell[p, a_]
                nbr[j] -= 1
                if nbr[j] < 0:
                    if periodic:
                        nbr[j] += shape[j]
                        a = F[j, _flat(nbr, shape)]
                    else:
                        a = 0.0
                else:
                    a = F[j, _flat(nbr, shape)]
                dj = b - a
                x0 = xi[p, j]
                q0 = a + dj * x0
                xn = x0 + q0 * tau_best * _E(dj * tau_best)
                xi[p, j] = min(max(xn, 0.0), 1.0)
            if jbest < 0:
                break                                    # reached the end of the tick inside the cell
            s = s + rs * tau_best * _E(D * tau_best)
            if s > 1.0:
                s = 1.0
            # cross the face
            ncross[p] += 1
            if side > 0:
                cell[p, jbest] += 1
                xi[p, jbest] = 0.0
                if cell[p, jbest] >= shape[jbest]:
                    cell[p, jbest] -= shape[jbest]
            else:
                cell[p, jbest] -= 1
                xi[p, jbest] = 1.0
                if cell[p, jbest] < 0:
                    cell[p, jbest] += shape[jbest]
            if it == maxcross:
                stuck += 1
    return ncross, stuck


def flow_tick(cell, xi, rho0, rho1, Flist, periodic=True, maxcross=100000):
    """S4 wrapper. Flist: list (one per axis) of net fluxes of the tick (summed over sub-steps)."""
    shape = np.array(rho0.shape, np.int64)
    F = np.stack([f.ravel() for f in Flist]).astype(np.float64)
    return pollock_move(cell, xi, rho0.ravel().astype(np.float64), rho1.ravel().astype(np.float64), F, shape,
                        periodic, maxcross)


# ======================================================================================================
#  S2b: BELL CONTINUOUS-TIME JUMPS IN FLOW TIME (net flux of the tick)
# ======================================================================================================
@nb.njit(cache=True)
def _bellflux(cell, rho0, rho1, F, shape, periodic, u, maxj):
    N, d = cell.shape
    njump = np.zeros(N, np.int64)
    nbr = np.zeros(d, np.int64)
    out = np.zeros(2 * d)
    k = 0
    for p in range(N):
        s = 0.0
        for it in range(maxj + 1):
            c = _flat(cell[p], shape)
            r0, r1 = rho0[c], rho1[c]
            D = r1 - r0
            rs = r0 + s * D
            if rs <= 0.0:
                rs = 1e-300
            tau_end = (1.0 - s) / rs * _phi(D * (1.0 - s) / rs)
            O = 0.0
            for j in range(d):
                b = F[j, c]
                for a_ in range(d):
                    nbr[a_] = cell[p, a_]
                nbr[j] -= 1
                a = 0.0
                if nbr[j] < 0:
                    if periodic:
                        nbr[j] += shape[j]
                        a = F[j, _flat(nbr, shape)]
                else:
                    a = F[j, _flat(nbr, shape)]
                out[2 * j] = b if b > 0 else 0.0
                out[2 * j + 1] = -a if a < 0 else 0.0
                O += out[2 * j] + out[2 * j + 1]
            if O <= 0.0:
                break
            E = -np.log(u[k % u.shape[0]]); k += 1
            tj = E / O
            if tj >= tau_end:
                break
            s = s + rs * tj * _E(D * tj)
            if s > 1.0:
                s = 1.0
            w = u[k % u.shape[0]] * O; k += 1
            acc = 0.0
            pick = 0
            for q in range(2 * d):
                acc += out[q]
                if w < acc:
                    pick = q
                    break
                pick = q
            j = pick // 2
            if pick % 2 == 0:
                cell[p, j] += 1
                if cell[p, j] >= shape[j]:
                    cell[p, j] -= shape[j]
            else:
                cell[p, j] -= 1
                if cell[p, j] < 0:
                    cell[p, j] += shape[j]
            njump[p] += 1
    return njump


def bellflux_tick(cell, rho0, rho1, Flist, rng, periodic=True, maxj=200):
    shape = np.array(rho0.shape, np.int64)
    F = np.stack([f.ravel() for f in Flist]).astype(np.float64)
    u = rng.random(4 * cell.shape[0] + 1024)
    u = np.maximum(u, 1e-300)
    return _bellflux(cell, rho0.ravel().astype(np.float64), rho1.ravel().astype(np.float64), F, shape, periodic,
                     u, maxj)


# ======================================================================================================
#  S1 / S2: DISCRETE BELL JUMPS PER SUB-STEP (pre_core family, w = 1)
# ======================================================================================================
def bell_substep(cell, j, A, B, rho, lam, rng, periodic=True):
    """one sub-step along axis j. lam = 1: ballistic (gross flows), lam = 0: minimal net flows.
    Returns (cell, number of jumps, Shannon entropy (bits) of the choices summed over surfers)."""
    Bn = np.roll(B, -1, axis=j)
    if lam == 1:
        Ff, Fb = A.copy(), Bn.copy()
    else:
        Ff, Fb = np.maximum(A - Bn, 0.0), np.maximum(Bn - A, 0.0)
    if not periodic:
        idx = [slice(None)] * A.ndim
        idx[j] = -1
        Ff[tuple(idx)] = 0.0
        Fb[tuple(idx)] = 0.0
    r = np.maximum(rho, 1e-300)
    pf = np.minimum(Ff / r, 1.0)
    pb = np.minimum(np.roll(Fb, 1, axis=j) / r, 1.0)
    t = tuple(cell[:, a] for a in range(cell.shape[1]))
    p1, p2 = pf[t], pb[t]
    p0 = np.clip(1.0 - p1 - p2, 0.0, 1.0)
    u = rng.random(len(cell))
    fwd = u < p1
    bwd = (~fwd) & (u < p1 + p2)
    step = fwd.astype(np.int64) - bwd.astype(np.int64)
    cell = cell.copy()
    n = rho.shape[j]
    if periodic:
        cell[:, j] = (cell[:, j] + step) % n
    else:
        cell[:, j] = np.clip(cell[:, j] + step, 0, n - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -sum(np.where(q > 0, q * np.log2(q), 0.0) for q in (p0, p1, p2))
    return cell, int(fwd.sum() + bwd.sum()), float(H.sum())


# ======================================================================================================
#  S3: MONOTONE CDF MAP ALONG LINES, PER SUB-STEP (R3), any d, numba
# ======================================================================================================
@nb.njit(cache=True)
def _cdf_line_move(cell, xi, j, R0, R1, shape, periodic):
    """R0, R1: (nlines, n) densities along axis j for every line; line index = flat index of the other coords.
    Non-periodic only (lines with walls) or periodic treated with a cut at the cell 0 face (exact as well,
    for the measure: the monotone map on the cut circle is a valid transport, but NOT the minimal one)."""
    N, d = cell.shape
    n = shape[j]
    nl = R0.shape[0]
    C0 = np.zeros((nl, n + 1)); C1 = np.zeros((nl, n + 1))
    for l in range(nl):
        for i in range(n):
            C0[l, i + 1] = C0[l, i] + R0[l, i]
            C1[l, i + 1] = C1[l, i] + R1[l, i]
    for p in range(N):
        # line index
        l = 0
        for a in range(d):
            if a != j:
                l = l * shape[a] + cell[p, a]
        i = cell[p, j]
        c = C0[l, i] + xi[p, j] * (C0[l, i + 1] - C0[l, i])
        # invert C1 (scale: totals equal up to rounding)
        tot0, tot1 = C0[l, n], C1[l, n]
        if tot0 <= 0 or tot1 <= 0:
            continue
        c = c * tot1 / tot0
        lo, hi = 0, n - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if C1[l, mid] <= c:
                lo = mid
            else:
                hi = mid - 1
        k = lo
        while k < n - 1 and C1[l, k + 1] <= c:
            k += 1
        den = C1[l, k + 1] - C1[l, k]
        fr = (c - C1[l, k]) / den if den > 0 else 0.5
        cell[p, j] = k
        xi[p, j] = min(max(fr, 0.0), 1.0)


def cdf_substep(cell, xi, j, rho0, rho1, periodic=False):
    d = rho0.ndim
    R0 = np.moveaxis(rho0, j, -1).reshape(-1, rho0.shape[j])
    R1 = np.moveaxis(rho1, j, -1).reshape(-1, rho1.shape[j])
    _cdf_line_move(cell, xi, j, np.ascontiguousarray(R0), np.ascontiguousarray(R1),
                   np.array(rho0.shape, np.int64), periodic)


# ======================================================================================================
#  S0: THE OLD READER (Rust zz.rs::muovi_surfisti, trilinear spinor + midpoint, 3-D lattice units)
# ======================================================================================================
SXm = np.array([[0, 1], [1, 0]], complex)
SYm = np.array([[0, -1j], [1j, 0]], complex)
SZm = np.array([[1, 0], [0, -1]], complex)


def trilinear(psi, X):
    """psi (4, n, n, n), X (N,3) continuous lattice coords (cell centre at i + 0.5, periodic). Returns (N, 4)."""
    n = psi.shape[1]
    f = X - 0.5
    i0 = np.floor(f).astype(np.int64)
    t = f - i0
    out = np.zeros((X.shape[0], psi.shape[0]), complex)
    for dx in (0, 1):
        wx = t[:, 0] if dx else 1 - t[:, 0]
        ix = (i0[:, 0] + dx) % n
        for dy in (0, 1):
            wy = t[:, 1] if dy else 1 - t[:, 1]
            iy = (i0[:, 1] + dy) % n
            for dz in (0, 1):
                wz = t[:, 2] if dz else 1 - t[:, 2]
                iz = (i0[:, 2] + dz) % n
                out += (wx * wy * wz)[:, None] * psi[:, ix, iy, iz].T
    return out


def bohm_dirac_v(sp):
    """sp (N, 4) chiral (L0, L1, R0, R1) -> v (N, 3) = (R^dag sigma R - L^dag sigma L)/rho (units of c)"""
    L, R = sp[:, 0:2], sp[:, 2:4]
    rho = np.sum(np.abs(sp)**2, 1)
    v = []
    for S in (SXm, SYm, SZm):
        jr = np.real(np.einsum("na,ab,nb->n", R.conj(), S, R))
        jl = np.real(np.einsum("na,ab,nb->n", L.conj(), S, L))
        v.append(jr - jl)
    return np.array(v).T / np.maximum(rho, 1e-300)[:, None]


def old_reader_move(psi, X, verso=1.0):
    """Rust muovi_surfisti in lattice units (one tick = one cell at c)."""
    v1 = bohm_dirac_v(trilinear(psi, X))
    q = X + 0.5 * verso * v1
    v2 = bohm_dirac_v(trilinear(psi, q))
    n = psi.shape[1]
    return (X + verso * v2) % n


# ======================================================================================================
#  SAMPLING AND STATISTICS
# ======================================================================================================
def born_sample(rng, rho, N):
    """N surfers from the piecewise-uniform measure rho (cells): cell index + uniform xi."""
    p = rho.ravel() / rho.sum()
    k = rng.choice(p.size, size=N, p=p)
    cell = np.stack(np.unravel_index(k, rho.shape), axis=1).astype(np.int64)
    xi = rng.random(cell.shape)
    return cell, xi


def chi2_z(cell, rho, block=1, min_expected=5.0):
    """chi^2 of the surfer histogram against rho on blocks of block^d cells; bins with expected < min_expected
    merged into one. Returns (z = (chi2 - dof)/sqrt(2 dof), chi2, dof, TV distance)."""
    d = rho.ndim
    N = cell.shape[0]
    shp = tuple(s // block for s in rho.shape)
    cb = cell // block
    ok = np.all(cb < np.array(shp), axis=1)
    cnt = np.zeros(shp)
    np.add.at(cnt, tuple(cb[ok, a] for a in range(d)), 1.0)
    rr = rho[tuple(slice(0, s * block) for s in shp)]
    for a in range(d):
        rr = rr.reshape(rr.shape[:a] + (shp[a], block) + rr.shape[a + 1:]).sum(a + 1)
    q = rr.ravel() / rho.sum()
    o = cnt.ravel()
    e = q * N
    big = e >= min_expected
    O = np.concatenate([o[big], [o[~big].sum()]])
    Ee = np.concatenate([e[big], [e[~big].sum() + N * (1 - q.sum())]])
    m = Ee > 0
    chi2 = float(np.sum((O[m] - Ee[m])**2 / Ee[m]))
    dof = int(m.sum() - 1)
    tv = 0.5 * float(np.abs(o / N - q).sum())
    return (chi2 - dof) / np.sqrt(2 * dof), chi2, dof, tv


def tv_null(rng, rho, N, block=1, reps=20):
    """mean TV distance of an IDEAL sampler (N iid surfers from rho), for scale"""
    out = []
    for _ in range(reps):
        c, _ = born_sample(rng, rho, N)
        out.append(chi2_z(c, rho, block)[3])
    return float(np.mean(out)), float(np.std(out))


# ======================================================================================================
#  STATES
# ======================================================================================================
def radial_state(kappa, n_level, a, Za, rmax, guess=None):
    """eigenstate of the radial rule (inverse iteration = infinite listening). Returns (omega, chi at tick
    boundary, RadialWalk)."""
    rw = RadialWalk(kappa, a, rmax, Zalpha=Za)
    if guess is None:
        k = abs(kappa)
        nr = n_level - k
        guess = 1.0 / np.sqrt(1 + (Za / (nr + np.sqrt(k * k - Za * Za)))**2)
    import scipy.sparse as sp
    import scipy.sparse.linalg as spl
    W = rw.matrix()
    lam0 = np.exp(-1j * guess * a)
    lu = spl.splu((W - lam0 * sp.identity(2 * rw.N, format="csc")).tocsc())
    x = np.ones(2 * rw.N, complex)
    x /= np.linalg.norm(x)
    for _ in range(12):
        y = lu.solve(x)
        x = y / np.linalg.norm(y)
    lam = np.vdot(x, W @ x)
    om = -np.angle(lam) / a
    chi = x.reshape(rw.N, 2).T.copy()
    return om, chi, rw


def spinor3d_from_radial(rw, chi, kappa, X, Y, Z):
    """put a radial state (kappa = -1: s1/2, +1: p1/2; m = +1/2) on 3-D points (same construction as Rust
    zz.rs::spinore + dirac_a_chirale). Returns chiral (L0, L1, R0, R1)."""
    g, f = rw.gf(rw.palindromic(chi))
    r = np.sqrt(X**2 + Y**2 + Z**2)
    r = np.maximum(r, 1e-9)
    G = np.interp(r, rw.r, g, right=0.0) / r
    Fr = np.interp(r, rw.r, f, right=0.0) / r
    y00 = 0.5 / np.sqrt(np.pi)
    om_s = (y00 * np.ones_like(r) + 0j, 0j * r)
    om_p = (-y00 * Z / r + 0j, -y00 * (X + 1j * Y) / r)
    up, lo = (om_s, om_p) if kappa < 0 else (om_p, om_s)
    d = [G * up[0], G * up[1], 1j * Fr * lo[0], 1j * Fr * lo[1]]
    s = 1 / np.sqrt(2)
    Lc = [(d[0] - d[2]) * s, (d[1] - d[3]) * s]
    Rc = [(d[0] + d[2]) * s, (d[1] + d[3]) * s]
    return np.array([Lc[0], Lc[1], Rc[0], Rc[1]])


def listen3d(walk, psi0, omega, ticks, a):
    """Hann-windowed lock-in at frequency omega (units mc^2/hbar) = filtered eigenstate of the 3-D rule."""
    acc = np.zeros_like(psi0)
    psi = psi0.copy()
    ov = []
    for t in range(ticks):
        w = np.sin(np.pi * (t + 0.5) / ticks)**2
        acc += w * np.exp(1j * omega * t * a) * psi
        ov.append(np.vdot(psi0, psi))
        psi, _ = walk.tick(psi)
    return acc / np.sqrt(np.sum(np.abs(acc)**2)), np.array(ov)


if __name__ == "__main__":
    import time
    from zz_core import Walk
    rng = np.random.default_rng(SEED)
    t0 = time.time()
    # 1) Walk3D == zz_core.Walk('dirac') (with a random potential)
    n = 12
    V = rng.random((n, n, n))
    a = 0.5
    W0 = Walk("dirac", n, a, theta=0.3, V=V)
    W1 = Walk3D(n, 0.3, V=V * a)
    psi = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n))
    p0, p1 = psi.copy(), psi.copy()
    err_flux = 0.0
    for _ in range(5):
        r0 = density(p1)
        p0 = W0.tick(p0)
        p1, subs = W1.tick(p1)
        Fs = [net_flux(j, A, B) for (j, A, B) in subs]
        div = sum(F - np.roll(F, 1, axis=j) for (j, _, _), F in zip(subs, Fs))
        err_flux = max(err_flux, np.max(np.abs(r0 - div - density(p1))) / np.max(r0))
    print("Walk3D vs zz_core.Walk: %.2e ; rho1 = rho0 - div F: %.2e" % (np.max(np.abs(p0 - p1)), err_flux))
    # 2) Pollock 1-D == monotone map (quantile kept) on a random tick
    n = 64
    W = Walk1D(n, 0.2)
    psi = rng.normal(size=(2, n)) + 1j * rng.normal(size=(2, n))
    psi /= np.linalg.norm(psi)
    r0 = density(psi)
    psi2, subs = W.tick(psi)
    r1 = density(psi2)
    F = net_flux(0, subs[0][1], subs[0][2])
    cell, xi = born_sample(rng, r0, 5000)
    C0 = np.concatenate([[0], np.cumsum(r0)])
    q0 = C0[cell[:, 0]] + xi[:, 0] * r0[cell[:, 0]]
    c4, x4 = cell.copy(), xi.copy()
    ncr, st = flow_tick(c4, x4, r0, r1, [F], periodic=True)
    C1 = np.concatenate([[0], np.cumsum(r1)])
    q1 = C1[c4[:, 0]] + x4[:, 0] * r1[c4[:, 0]]
    dq = (q1 - q0 + 0.5) % 1.0 - 0.5       # periodic: quantile defined up to a common shift (flux through the cut)
    print("1-D flow: quantile shift spread %.2e (a common shift = flux through cell 0's face)" % (dq.max() - dq.min()))
    print("time %.1f s" % (time.time() - t0))

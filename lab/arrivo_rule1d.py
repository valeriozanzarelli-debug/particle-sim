"""arrivo_rule1d.py - the Das-Duerr set-up on the zig-zag rule, 1+1 D (longitudinal motion).

THE RULE (frozen, lab/zz_core.py, kind 'dirac1'): psi = (R, L) on sites z_j = (j + 1/2) a, j = 0..N-1;
one tick:  R moves one cell to +z, L one cell to -z (at c), then the flip
           R' = cos th R - i sin th L,  L' = cos th L - i sin th R,   th = m c a / hbar = a  (hbar = c = m = 1).
Additions of this study (declared, nothing else):
  * hard wall at z = 0 (the closed end of the guide): the L amplitude that would leave site 0 to the left
    comes back as R on site 0 with a factor -1 (the same reflection as RadialWalk at r = 0: the large
    component (R + L)/sqrt2 vanishes at the wall);
  * the initial box 0 < z < d: a second wall of the same kind at z = d (only for the preparation);
    the initial state is the lowest positive-energy eigenvector of the one-tick operator in the box
    (dense diagonalisation, this is the 'ground state of the rule in the box');
  * numerical edge of the grid (NOT physics): either the exact light cone (grid longer than c t_end: nothing
    can reach the edge), or a smooth complex absorbing ramp far behind the detector plane (checked against
    the light cone);
  * the rule's own detector: a cell that removes everything that enters it (= nothing comes back from
    beyond z = L: L(z_L) = 0 at the boundary = Tumulka's ideal ABCD  alpha_z psi = psi).
Read-outs: flux across the bond at z = L (|R_{jL-1}|^2 - |L_{jL}|^2 per tick, exactly conserved), absorbed
probability per tick, density, and 'hybrid' Bohmian trajectories: the longitudinal surfer velocity of the
rule v = (|R|^2 - |L|^2)/rho plus the NR spin term of the transverse factor (debt, see docs/TEMPI_ARRIVO.md):
    Zdot = v(Z) + omega_r Y sin(alpha),   Ydot = rho'(Z)/(2 rho) sin(alpha) + omega_r X cos(alpha),
    Xdot = -omega_r Y cos(alpha).
"""
import numpy as np
import numba as nb


# ------------------------------------------------------------------ box ground state of the rule
def box_state(D, a):
    """lowest positive quasi-energy eigenvector of the one-tick operator of the rule in the box (0, D)
    with the two walls. returns (R, L) on the box sites (n = round(D/a)), and the quasi-energy omega."""
    n = int(round(D / a))
    th = a
    c, s = np.cos(th), np.sin(th)
    # shift S: R'_j = R_{j-1} (j >= 1), R'_0 = -L_0 ; L'_j = L_{j+1} (j <= n-2), L'_{n-1} = -R_{n-1}
    # index: R_j -> j, L_j -> n + j
    Wm = np.zeros((2 * n, 2 * n), complex)
    S = np.zeros((2 * n, 2 * n))
    for j in range(n):
        if j >= 1:
            S[j, j - 1] = 1.0
        else:
            S[0, n + 0] = -1.0
        if j <= n - 2:
            S[n + j, n + j + 1] = 1.0
        else:
            S[n + n - 1, n - 1] = -1.0
    C = np.zeros((2 * n, 2 * n), complex)
    for j in range(n):
        C[j, j] = c
        C[j, n + j] = -1j * s
        C[n + j, n + j] = c
        C[n + j, j] = -1j * s
    Wm = C @ S
    lam, V = np.linalg.eig(Wm)
    om = -np.angle(lam) / a
    # positive branch near m + pi^2/(2 m D^2)
    k = np.pi / D
    target = np.arccos(np.cos(th) * np.cos(k * a)) / a
    i = int(np.argmin(np.abs(om - target)))
    v = V[:, i]
    v = v / np.linalg.norm(v)
    R, L = v[:n].copy(), v[n:].copy()
    # global phase: make the large component real positive at the centre
    g = R[n // 2] + L[n // 2]
    ph = np.exp(-1j * np.angle(g))
    return R * ph, L * ph, om[i], target


# ------------------------------------------------------------------ kernels
@nb.njit(cache=True, fastmath=False)
def _tick(R, L, Rn, Ln, c, s, nact, absorber_end, damp, jflux):
    """one tick on sites 0..nact-1 (beyond nact everything is zero). returns (flux across bond jflux-1|jflux,
    probability absorbed at the end cell if absorber_end). damp: per-site factor (complex absorbing ramp) or
    empty array."""
    N = R.shape[0]
    fl = 0.0
    if jflux > 0 and jflux < N:
        fl = abs(R[jflux - 1])**2 - abs(L[jflux])**2
    ab = 0.0
    # shift
    Rn[0] = -L[0]
    for j in range(1, nact):
        Rn[j] = R[j - 1]
    for j in range(0, nact - 1):
        Ln[j] = L[j + 1]
    if absorber_end:
        ab = abs(R[nact - 1])**2
        Ln[nact - 1] = 0.0
    else:
        Ln[nact - 1] = L[nact] if nact < N else 0.0
        if nact < N:
            Rn[nact] = R[nact - 1]
    m = nact + 1 if (not absorber_end and nact < N) else nact
    # coin
    for j in range(m):
        r = Rn[j]
        l = Ln[j]
        Rn[j] = c * r - 1j * s * l
        Ln[j] = c * l - 1j * s * r
    if damp.shape[0] > 0:
        for j in range(m):
            Rn[j] *= damp[j]
            Ln[j] *= damp[j]
    return fl, ab, m


@nb.njit(cache=True)
def _fields(R, L, a, m, rho, J, drho):
    # J = conserved lattice current, centred on the site: mean of the two bond fluxes
    # F_{j+1/2} = |R_j|^2 - |L_{j+1}|^2 (what crosses the bond in the next tick); F_{-1/2} = 0 (wall).
    # (the single-site formula |R_j|^2 - |L_j|^2 of zz_core.dirac_velocity is off by O(a) grad rho: see doc)
    n = R.shape[0]
    for j in range(m):
        rho[j] = abs(R[j])**2 + abs(L[j])**2
    for j in range(m):
        fr = abs(R[j])**2 - (abs(L[j + 1])**2 if j + 1 < n else 0.0)
        fl = (abs(R[j - 1])**2 - abs(L[j])**2) if j >= 1 else 0.0
        J[j] = 0.5 * (fr + fl)
    for j in range(m, rho.shape[0]):
        rho[j] = 0.0
        J[j] = 0.0
    for j in range(m):
        if j == 0:
            # wall at z = 0: rho is even about the wall image (rho(-z) = rho(z)) -> one-sided
            drho[j] = (rho[1] - rho[0]) / (2 * a)     # image site -1 has rho_{-1} = rho_0
        elif j == rho.shape[0] - 1:
            drho[j] = (rho[j] - rho[j - 1]) / a
        else:
            drho[j] = (rho[j + 1] - rho[j - 1]) / (2 * a)


@nb.njit(cache=True)
def _interp(f, z, a):
    x = z / a - 0.5
    if x <= 0.0:
        return f[0]
    i = int(x)
    if i >= f.shape[0] - 1:
        return f[f.shape[0] - 1]
    w = x - i
    return (1 - w) * f[i] + w * f[i + 1]


@nb.njit(cache=True)
def _vel(z, y, x, a, om, ca, sa, rho, J, drho):
    r = _interp(rho, z, a)
    if r < 1e-300:
        r = 1e-300
    v = _interp(J, z, a) / r
    g = 0.5 * _interp(drho, z, a) / r
    return -om * y * ca, g * sa + om * x * ca, v + om * y * sa


@nb.njit(cache=True)
def _traj_step(X, Y, Z, T, alive, t0, h, nsub, a, om, ca, sa, zL,
               rho0, J0, d0, rho1, J1, d1, tmp_r, tmp_J, tmp_d):
    """advance trajectories from time t0 to t0 + h (fields linearly interpolated in time between the two
    stored snapshots), nsub Heun substeps."""
    hs = h / nsub
    for k in range(nsub):
        w0 = k / nsub
        w1 = (k + 1) / nsub
        for i in range(X.shape[0]):
            if not alive[i]:
                continue
            # fields at w0 and w1 interpolated on the fly (pointwise)
            x, y, z = X[i], Y[i], Z[i]
            # k1 at time w0
            r = (1 - w0) * _interp(rho0, z, a) + w0 * _interp(rho1, z, a)
            jj = (1 - w0) * _interp(J0, z, a) + w0 * _interp(J1, z, a)
            dd = (1 - w0) * _interp(d0, z, a) + w0 * _interp(d1, z, a)
            if r < 1e-300:
                r = 1e-300
            ax = -om * y * ca
            ay = 0.5 * dd / r * sa + om * x * ca
            az = jj / r + om * y * sa
            xp, yp, zp = x + hs * ax, y + hs * ay, abs(z + hs * az)
            r = (1 - w1) * _interp(rho0, zp, a) + w1 * _interp(rho1, zp, a)
            jj = (1 - w1) * _interp(J0, zp, a) + w1 * _interp(J1, zp, a)
            dd = (1 - w1) * _interp(d0, zp, a) + w1 * _interp(d1, zp, a)
            if r < 1e-300:
                r = 1e-300
            bx = -om * yp * ca
            by = 0.5 * dd / r * sa + om * xp * ca
            bz = jj / r + om * yp * sa
            xn = x + 0.5 * hs * (ax + bx)
            yn = y + 0.5 * hs * (ay + by)
            zn = abs(z + 0.5 * hs * (az + bz))
            if zn >= zL:
                fr = (zL - z) / max(zn - z, 1e-300)
                if fr < 0:
                    fr = 0.0
                if fr > 1:
                    fr = 1.0
                T[i] = t0 + k * hs + fr * hs
                alive[i] = False
            X[i], Y[i], Z[i] = xn, yn, zn


def damp_ramp(N, a, z_c, ell, W0):
    """complex absorbing ramp exp(-W(z) a), W = W0 ((z - z_c)/ell)^2 for z > z_c (numerical edge only)"""
    z = (np.arange(N) + 0.5) * a
    W = np.where(z > z_c, W0 * ((z - z_c) / ell)**2, 0.0)
    return np.exp(-W * a)


class Run1D:
    """Released box on the rule. geometry in units of lambda_C; times in hbar/(m c^2).
    mode 'free'     : grid long enough (light cone) or with absorbing ramp; flux recorded at zL.
    mode 'detector' : the rule's cell detector at zL (grid ends there)."""

    def __init__(self, D, a, zL, mode, t_end, edge="cone", ramp=None):
        self.D, self.a, self.zL, self.mode, self.t_end = D, a, zL, mode, t_end
        self.th = a
        self.c, self.s = np.cos(a), np.sin(a)
        self.jL = int(round(zL / a))
        R0, L0, om, target = box_state(D, a)
        self.omega_box, self.omega_target = om, target
        nb_ = R0.shape[0]
        if mode == "detector":
            N = self.jL
            self.absorber = True
            self.damp = np.zeros(0)
        else:
            self.absorber = False
            if edge == "cone":
                N = nb_ + int(np.ceil(t_end / a)) + 4
                self.damp = np.zeros(0)
            else:
                z_c, ell, W0 = ramp
                N = int(np.ceil((z_c + ell) / a))
                self.damp = damp_ramp(N, a, z_c, ell, W0)
        self.N = N
        self.R = np.zeros(N, complex)
        self.L = np.zeros(N, complex)
        self.R[:nb_], self.L[:nb_] = R0, L0
        self.Rn = np.zeros(N, complex)
        self.Ln = np.zeros(N, complex)
        self.nact = nb_ if (mode != "detector" and edge == "cone") else N
        self.tick = 0
        self.flux = []
        self.absorbed = []

    def advance(self, nticks):
        fl = np.empty(nticks)
        ab = np.empty(nticks)
        for k in range(nticks):
            f, b, m = _tick(self.R, self.L, self.Rn, self.Ln, self.c, self.s, self.nact,
                            self.absorber, self.damp, self.jL)
            self.R, self.Rn = self.Rn, self.R
            self.L, self.Ln = self.Ln, self.L
            self.nact = m
            fl[k], ab[k] = f, b
            self.tick += 1
        self.flux.append(fl)
        self.absorbed.append(ab)

    def norm(self):
        return float(np.sum(np.abs(self.R)**2 + np.abs(self.L)**2))


@nb.njit(cache=True)
def _advance_with_traj(R, L, Rn, Ln, c, s, nact, absorber, damp, jL, nticks, a, tick0,
                       stride_ticks, nsub, X, Y, Z, T, alive, om, ca, sa, zL,
                       rho0, J0, d0, rho1, J1, d1, fl_out, ab_out, last_t, do_traj):
    """advance the field nticks, updating trajectories every stride_ticks ticks (Heun, nsub substeps).
    The state ends in (R, L) (an odd number of internal swaps is copied back)."""
    m = nact
    tsnap = last_t
    Ra, La, Rb, Lb = R, L, Rn, Ln
    for k in range(nticks):
        f, b, m = _tick(Ra, La, Rb, Lb, c, s, m, absorber, damp, jL)
        Ra, Rb = Rb, Ra
        La, Lb = Lb, La
        fl_out[k] = f
        ab_out[k] = b
        if do_traj and (tick0 + k + 1) % stride_ticks == 0:
            tnow = (tick0 + k + 1) * a
            _fields(Ra, La, a, Ra.shape[0], rho1, J1, d1)
            _traj_step(X, Y, Z, T, alive, tsnap, tnow - tsnap, nsub, a, om, ca, sa, zL,
                       rho0, J0, d0, rho1, J1, d1, rho1, J1, d1)
            for j in range(rho0.shape[0]):
                rho0[j] = rho1[j]
                J0[j] = J1[j]
                d0[j] = d1[j]
            tsnap = tnow
    if nticks % 2 == 1:
        for j in range(R.shape[0]):
            R[j] = Rn[j]       # after an odd number of swaps the newest state is in Rn, Ln
            L[j] = Ln[j]
    return m, tsnap


# ------------------------------------------------------------------ driver: hybrid trajectories + flux
SCHEDULE = ((0.05, 2e-4), (2.0, 2e-3), (1e9, 1e-2))      # (until tau, step in DD units), as the NR reference


def run_hybrid(D, a, Lz, omega, alpha, XYZ0, tau_end, edge="cone", ramp=None, schedule=SCHEDULE,
               do_traj=True, verbose=False):
    """Das-Duerr on the rule. Lz, omega, tau in DD units (d = 1). XYZ0: initial (X, Y, Z) in DD units.
    returns dict: tau (bond times), flux density (DD units), arrival times T (DD units, nan = not arrived)."""
    import time
    t0 = time.time()
    zL = Lz * D
    t_end = tau_end * D**2
    run = Run1D(D, a, zL, "free", t_end, edge=edge, ramp=ramp)
    N = run.N
    X = XYZ0[0] * D
    Y = XYZ0[1] * D
    Z = XYZ0[2] * D
    T = np.full(X.shape[0], np.nan)
    alive = np.ones(X.shape[0], np.bool_)
    om = omega / D**2
    ca, sa = np.cos(alpha), np.sin(alpha)
    rho0, J0, d0 = np.zeros(N), np.zeros(N), np.zeros(N)
    rho1, J1, d1 = np.zeros(N), np.zeros(N), np.zeros(N)
    _fields(run.R, run.L, a, N, rho0, J0, d0)
    fl_all, tick = [], 0
    tsnap = 0.0
    m = run.nact
    n_end = int(np.ceil(t_end / a))
    for (tau_until, hDD) in schedule:
        n_until = min(n_end, int(np.ceil(tau_until * D**2 / a)))
        if n_until <= tick:
            continue
        stride = max(1, int(round(hDD * D**2 / a)))
        nsub = max(1, int(np.ceil(a / (hDD * D**2))))
        # make the segment a multiple of the stride
        nt = n_until - tick
        fl = np.empty(nt)
        ab = np.empty(nt)
        m, tsnap = _advance_with_traj(run.R, run.L, run.Rn, run.Ln, run.c, run.s, m, False, run.damp, run.jL,
                                      nt, a, tick, stride, nsub, X, Y, Z, T, alive, om, ca, sa, zL,
                                      rho0, J0, d0, rho1, J1, d1, fl, ab, tsnap, do_traj)
        tick += nt
        fl_all.append(fl)
        if verbose:
            print("  seg to tau %.3g: ticks %d stride %d nsub %d alive %d  (%.0f s)" %
                  (tick * a / D**2, tick, stride, nsub, alive.sum(), time.time() - t0), flush=True)
        if tick >= n_end:
            break
    fl = np.concatenate(fl_all) / a * D**2
    tau = (np.arange(fl.shape[0]) + 0.5) * a / D**2
    norm = float(np.sum(np.abs(run.R)**2 + np.abs(run.L)**2))
    return dict(tau=tau, flux=fl, T=T / D**2, norm=norm, omega_box=run.omega_box, N=N, secs=time.time() - t0)

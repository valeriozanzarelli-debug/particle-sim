"""arrivo_rule2d.py - the Das-Duerr set-up on the zig-zag rule in 2+1 D (y, z): a parallel-plate waveguide.

THE RULE (frozen, lab/zz_core.py, kind 'dirac', axes (y, z), k_x = 0): psi = (L0, L1, R0, R1) on sites
(y_i, z_j) = ((i - Ny/2 + 1/2) a, (j + 1/2) a). One tick: translate along y then z (even ticks) or z then y
(odd ticks): R splits on the eigenvectors of sigma_y (sigma_z) and the +1 part moves +1 cell, the -1 part -1
cell; L the opposite (S = -sigma). Then the flip L' = e^{-iVa}(cos th L - i sin th R), R' = e^{-iVa}(cos th R
- i sin th L), th = m c a/hbar = a.
Additions of this study (declared):
  * the guide: V(y) = m omega_r^2 y^2 / 2 (an external potential applied as a phase, like the Coulomb phase of
    the atom; the plates of the guide are a hand-put potential = debt), y periodic far outside;
  * hard wall at z = 0 (and at z = d for the preparation): a -z mover that would leave the box comes back
    as the +z mover with the same spin and the other handedness, factor -1 (R1 <-> L1, L0 <-> R0: the large
    component (L + R)/sqrt2 vanishes at the wall; spin is not touched);
  * initial state: spectral filter of the guess phi0(y) sin(pi z/d) chi with the rule itself in the box
    (Gaussian window in time centred on the Rayleigh quasi-energy): the ground state of the rule;
  * numerical edge of the grid (complex absorbing ramp far behind the detector plane), or the rule's cell
    detector at z = L (what moves beyond z_L is removed; nothing comes back: alpha_z psi = psi there).
Read-outs: bond flux through z = L, absorbed probability per tick, and surfers (Bohm-Dirac trajectories) with
the conserved lattice current, bond-centred: v_j = (F_{j-1/2} + F_{j+1/2}) / (2 rho) along each axis.
NO spin term is added by hand here: whatever spin dependence the surfers show comes from the rule.
"""
import numpy as np
import numba as nb

S2 = 1.0 / np.sqrt(2.0)


@nb.njit(cache=True)
def _shift_y(P, Q):
    """y translation of P into Q (periodic in y). components 0,1 = L, 2,3 = R."""
    Ny, Nz = P.shape[1], P.shape[2]
    for i in range(Ny):
        ip = (i + 1) % Ny       # source for a +y mover arriving at i is i-1; we push instead
        im = (i - 1) % Ny
        for j in range(Nz):
            # R: c+ = (R0 - i R1)/sqrt2 moves +y ; c- = (R0 + i R1)/sqrt2 moves -y
            # value at i: c+ from i-1, c- from i+1
            r0a, r1a = P[2, im, j], P[3, im, j]
            r0b, r1b = P[2, ip, j], P[3, ip, j]
            cp = (r0a - 1j * r1a) * S2
            cm = (r0b + 1j * r1b) * S2
            Q[2, i, j] = (cp + cm) * S2
            Q[3, i, j] = (1j * cp - 1j * cm) * S2
            # L: c+ (sigma_y = +1) moves -y (from i+1), c- moves +y (from i-1)
            l0a, l1a = P[0, ip, j], P[1, ip, j]
            l0b, l1b = P[0, im, j], P[1, im, j]
            cp = (l0a - 1j * l1a) * S2
            cm = (l0b + 1j * l1b) * S2
            Q[0, i, j] = (cp + cm) * S2
            Q[1, i, j] = (1j * cp - 1j * cm) * S2


@nb.njit(cache=True)
def _shift_z(P, Q, top_mode, jL):
    """z translation. +z movers: R0, L1 ; -z movers: R1, L0. wall at z=0 always.
    top_mode 0: open (lost), 1: wall (box), 2: cell detector (absorbed, returned).
    returns (absorbed, flux across bond jL-1|jL)."""
    Ny, Nz = P.shape[1], P.shape[2]
    ab = 0.0
    fl = 0.0
    for i in range(Ny):
        if jL > 0 and jL < Nz:
            fl += (abs(P[2, i, jL - 1])**2 + abs(P[1, i, jL - 1])**2
                   - abs(P[3, i, jL])**2 - abs(P[0, i, jL])**2)
        for j in range(1, Nz):
            Q[2, i, j] = P[2, i, j - 1]
            Q[1, i, j] = P[1, i, j - 1]
        for j in range(0, Nz - 1):
            Q[3, i, j] = P[3, i, j + 1]
            Q[0, i, j] = P[0, i, j + 1]
        # bottom wall: -z movers at j=0 come back as +z movers, same spin, other handedness, factor -1
        Q[2, i, 0] = -P[0, i, 0]       # L0 (-z) -> R0 (+z)
        Q[1, i, 0] = -P[3, i, 0]       # R1 (-z) -> L1 (+z)
        t = Nz - 1
        if top_mode == 1:
            Q[3, i, t] = -P[1, i, t]   # L1 (+z) -> R1 (-z)
            Q[0, i, t] = -P[2, i, t]   # R0 (+z) -> L0 (-z)
        else:
            if top_mode == 2:
                ab += abs(P[2, i, t])**2 + abs(P[1, i, t])**2
            Q[3, i, t] = 0.0
            Q[0, i, t] = 0.0
    return ab, fl


@nb.njit(cache=True)
def _coin(P, c, s, ph, damp):
    Ny, Nz = P.shape[1], P.shape[2]
    for i in range(Ny):
        e = ph[i]
        for j in range(Nz):
            d = e * damp[j]
            l0, l1, r0, r1 = P[0, i, j], P[1, i, j], P[2, i, j], P[3, i, j]
            P[0, i, j] = d * (c * l0 - 1j * s * r0)
            P[1, i, j] = d * (c * l1 - 1j * s * r1)
            P[2, i, j] = d * (c * r0 - 1j * s * l0)
            P[3, i, j] = d * (c * r1 - 1j * s * l1)


@nb.njit(cache=True)
def _ticks(P, Q, n0, nticks, c, s, ph, damp, top_mode, jL, fl_out, ab_out):
    """nticks of the rule; the state is returned in P (copied back if needed)."""
    A, B = P, Q
    for k in range(nticks):
        n = n0 + k
        ab = 0.0
        fl = 0.0
        if n % 2 == 0:
            _shift_y(A, B)
            ab, fl = _shift_z(B, A, top_mode, jL)
        else:
            ab, fl = _shift_z(A, B, top_mode, jL)
            _shift_y(B, A)
        _coin(A, c, s, ph, damp)
        fl_out[k] = fl
        ab_out[k] = ab
    return 0


@nb.njit(cache=True)
def _fields2(P, rho, Jy, Jz):
    """density and bond-centred conserved currents (per site) along y and z."""
    Ny, Nz = P.shape[1], P.shape[2]
    # movers along z: +: |R0|^2 + |L1|^2 ; -: |R1|^2 + |L0|^2
    # movers along y: R: +y |(R0 - i R1)|^2/2, -y |(R0 + i R1)|^2/2 ; L: -y |(L0 - i L1)|^2/2, +y |(L0 + i L1)|^2/2
    for i in range(Ny):
        for j in range(Nz):
            rho[i, j] = abs(P[0, i, j])**2 + abs(P[1, i, j])**2 + abs(P[2, i, j])**2 + abs(P[3, i, j])**2
    for i in range(Ny):
        ip = (i + 1) % Ny
        im = (i - 1) % Ny
        for j in range(Nz):
            # z
            pz = abs(P[2, i, j])**2 + abs(P[1, i, j])**2
            mz_next = (abs(P[3, i, j + 1])**2 + abs(P[0, i, j + 1])**2) if j + 1 < Nz else 0.0
            fr = pz - mz_next
            if j >= 1:
                pz_prev = abs(P[2, i, j - 1])**2 + abs(P[1, i, j - 1])**2
                mz = abs(P[3, i, j])**2 + abs(P[0, i, j])**2
                fl = pz_prev - mz
            else:
                fl = 0.0
            Jz[i, j] = 0.5 * (fr + fl)
            # y
            def_p = 0.5 * (abs(P[2, i, j] - 1j * P[3, i, j])**2 + abs(P[0, i, j] + 1j * P[1, i, j])**2)
            m_next = 0.5 * (abs(P[2, ip, j] + 1j * P[3, ip, j])**2 + abs(P[0, ip, j] - 1j * P[1, ip, j])**2)
            p_prev = 0.5 * (abs(P[2, im, j] - 1j * P[3, im, j])**2 + abs(P[0, im, j] + 1j * P[1, im, j])**2)
            m_here = 0.5 * (abs(P[2, i, j] + 1j * P[3, i, j])**2 + abs(P[0, i, j] - 1j * P[1, i, j])**2)
            Jy[i, j] = 0.5 * ((def_p - m_next) + (p_prev - m_here))


@nb.njit(cache=True)
def _bil(f, y, z, a, y0):
    """bilinear interpolation; y0 = y coordinate of row 0 centre"""
    Ny, Nz = f.shape
    xi = (y - y0) / a
    xj = z / a - 0.5
    if xj < 0.0:
        xj = 0.0
    i = int(np.floor(xi))
    j = int(xj)
    if i < 0:
        i = 0
        xi = 0.0
    if i >= Ny - 1:
        i = Ny - 2
        xi = Ny - 1.0
    if j >= Nz - 1:
        j = Nz - 2
        xj = Nz - 1.0
    u = xi - i
    w = xj - j
    return ((1 - u) * (1 - w) * f[i, j] + u * (1 - w) * f[i + 1, j]
            + (1 - u) * w * f[i, j + 1] + u * w * f[i + 1, j + 1])


@nb.njit(cache=True)
def _traj2(Y, Z, T, alive, t0, h, a, y0, zL, r0, y0f, z0f, r1, y1f, z1f):
    """Heun step from t0 to t0+h with fields at the two ends (rho, Jy, Jz)."""
    for k in range(Y.shape[0]):
        if not alive[k]:
            continue
        y, z = Y[k], Z[k]
        r = _bil(r0, y, z, a, y0)
        if r < 1e-300:
            r = 1e-300
        ay = _bil(y0f, y, z, a, y0) / r
        az = _bil(z0f, y, z, a, y0) / r
        yp, zp = y + h * ay, abs(z + h * az)
        r = _bil(r1, yp, zp, a, y0)
        if r < 1e-300:
            r = 1e-300
        by = _bil(y1f, yp, zp, a, y0) / r
        bz = _bil(z1f, yp, zp, a, y0) / r
        yn = y + 0.5 * h * (ay + by)
        zn = abs(z + 0.5 * h * (az + bz))
        if zn >= zL:
            fr = (zL - z) / max(zn - z, 1e-300)
            fr = min(max(fr, 0.0), 1.0)
            T[k] = t0 + fr * h
            alive[k] = False
        Y[k], Z[k] = yn, zn


class Guide2D:
    def __init__(self, D, a, omega_dd, Lz, ymax_dd, mode, zc_dd=None, ell_dd=None, W0=0.05):
        """D: box length in lambda_C; omega_dd, Lz, ymax in DD units (d = 1). mode 'free' or 'detector'."""
        self.D, self.a = D, a
        self.c, self.s = np.cos(a), np.sin(a)
        self.om_r = omega_dd / D**2
        self.Ny = 2 * int(np.ceil(ymax_dd * D / a))
        self.y = (np.arange(self.Ny) - self.Ny / 2 + 0.5) * a
        self.y0 = self.y[0]
        self.jL = int(round(Lz * D / a))
        self.zL = self.jL * a
        self.mode = mode
        if mode == "detector":
            self.Nz = self.jL
            self.top = 2
            self.damp = np.ones(self.Nz)
        else:
            self.top = 0
            zc, ell = zc_dd * D, ell_dd * D
            self.Nz = int(np.ceil((zc + ell) / a))
            z = (np.arange(self.Nz) + 0.5) * a
            Wc = np.where(z > zc, W0 * ((z - zc) / ell)**2, 0.0)
            self.damp = np.exp(-Wc * a)
        V = 0.5 * self.om_r**2 * self.y**2
        self.ph = np.exp(-1j * V * a)
        self.z = (np.arange(self.Nz) + 0.5) * a
        self.tick = 0

    def prepare(self, spinor, Tf_factor=12.0, verbose=False):
        """ground state of the rule in the box (0, D) x guide: filter of the guess with the rule itself"""
        D, a = self.D, self.a
        nb_ = int(round(D / a))
        chi = np.asarray(spinor, complex)
        chi = chi / np.linalg.norm(chi)
        zb = (np.arange(nb_) + 0.5) * a
        g = np.exp(-0.5 * self.om_r * self.y**2)[:, None] * np.sin(np.pi * zb / D)[None, :]
        P = np.zeros((4, self.Ny, nb_), complex)
        for q in range(2):
            P[q] = S2 * chi[q] * g       # L
            P[2 + q] = S2 * chi[q] * g   # R
        P /= np.sqrt(np.sum(np.abs(P)**2))
        Q = np.zeros_like(P)
        damp = np.ones(nb_)
        fo, ao = np.zeros(1), np.zeros(1)

        def W2(Pin, n0):
            X = Pin.copy()
            Y = np.zeros_like(X)
            _ticks(X, Y, n0, 2, self.c, self.s, self.ph, damp, 1, -1, np.zeros(2), np.zeros(2))
            return X
        def filt(P0, Om, sig_t):
            """Gaussian filter in time centred on quasi-energy Om, width sig_t, truncated at +-4 sig_t"""
            npair = int(np.ceil(4 * sig_t / (2 * a)))
            acc = np.zeros_like(P0)
            X = P0.copy()
            for n in range(2 * npair + 1):
                t = 2 * a * (n - npair)
                acc += np.exp(-0.5 * (t / sig_t)**2) * np.exp(1j * Om * 2 * a * n) * X
                X = W2(X, 0)
            return acc / np.sqrt(np.sum(np.abs(acc)**2)), 4 * npair
        # stage 1: broad filter around the NR estimate (removes the negative-energy part, which biases the
        # Rayleigh quotient), stage 2: Rayleigh quasi-energy, stage 3: narrow filter (twice, re-centred)
        Om0 = 1.0 + np.pi**2 / (2 * D**2) + 0.5 * self.om_r
        P1, nt1 = filt(P, Om0, 5.0)
        lam = np.vdot(P1, W2(P1, 0))
        Om = -np.angle(lam) / (2 * a)
        gap = min(3 * np.pi**2 / (2 * D**2), 2 * self.om_r)
        sig = Tf_factor / (4 * gap)          # energy resolution gap / Tf_factor * 4
        acc, nt2 = filt(P1, Om, sig)
        Om = -np.angle(np.vdot(acc, W2(acc, 0))) / (2 * a)
        acc, nt3 = filt(acc, Om, sig)
        npair = (nt1 + nt2 + nt3) // 4
        acc /= np.sqrt(np.sum(np.abs(acc)**2))
        lam2 = np.vdot(acc, W2(acc, 0))
        resid = np.linalg.norm(W2(acc, 0) - lam2 * acc)
        if verbose:
            print("  filter: Om %.8f -> %.8f, |lam| %.12f, residual %.2e, ticks %d" %
                  (Om, -np.angle(lam2) / (2 * a), abs(lam2), resid, 4 * npair))
        self.P = np.zeros((4, self.Ny, self.Nz), complex)
        self.P[:, :, :nb_] = acc
        self.Q = np.zeros_like(self.P)
        self.tick = 0
        self.prep_resid = float(resid)
        self.prep_Om = float(-np.angle(lam2) / (2 * a))
        return acc

    def spin_expect(self):
        """<sigma> of the large component (L + R)/sqrt2 (approximately the Pauli spin)"""
        G = (self.P[0:2] + self.P[2:4]) * S2
        sx = 2 * np.real(np.sum(np.conj(G[0]) * G[1]))
        sy = 2 * np.imag(np.sum(np.conj(G[0]) * G[1]))
        sz = np.sum(np.abs(G[0])**2 - np.abs(G[1])**2)
        n = np.sum(np.abs(G)**2)
        return np.array([sx, sy, sz]) / n

    def advance(self, nticks):
        fl = np.zeros(nticks)
        ab = np.zeros(nticks)
        _ticks(self.P, self.Q, self.tick, nticks, self.c, self.s, self.ph, self.damp, self.top, self.jL, fl, ab)
        self.tick += nticks
        return fl, ab

    def fields(self):
        rho = np.zeros((self.Ny, self.Nz))
        Jy = np.zeros_like(rho)
        Jz = np.zeros_like(rho)
        _fields2(self.P, rho, Jy, Jz)
        return rho, Jy, Jz

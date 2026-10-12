"""glue_core.py - 2+1 D SU(2) pure gauge theory on a periodic lattice (Wilson action).

Links are SU(2) matrices stored as unit quaternions U = a0 + i a.sigma  ->  array (..., 4).
Lattice array U[mu, t, x, y, 4], mu = 0 (time), 1 (x), 2 (y).  Action S = beta sum_P (1 - 1/2 tr U_P).
The ONLY input is beta = 4/(a g^2): there is no mass anywhere in the rule.

Update: checkerboard Kennedy-Pendleton heat-bath + over-relaxation (exact, vectorised numpy).
Measurements (per configuration, per time slice, zero spatial momentum):
  - APE-smeared spatial links (several levels), spatial Wilson loops of several shapes ->
    0++ (A1 of the square group) and 2++ (B1 = loop(0 deg) - loop(90 deg)) operators;
  - smeared spatial Polyakov loops (torelon = flux tube wound around the box) -> string tension;
  - unsmeared planar Wilson loops W(R,T) (all planes) -> Creutz ratios.
"""
import numpy as np


# ---------------- quaternion algebra ----------------
def qmul(a, b):
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    out = np.empty(np.broadcast_shapes(a.shape, b.shape), dtype=a.dtype)
    out[..., 0] = a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3
    # vector part: a0 b + b0 a - a x b   (U = a0 + i a.sigma)
    out[..., 1] = a0 * b1 + b0 * a1 - (a2 * b3 - a3 * b2)
    out[..., 2] = a0 * b2 + b0 * a2 - (a3 * b1 - a1 * b3)
    out[..., 3] = a0 * b3 + b0 * a3 - (a1 * b2 - a2 * b1)
    return out


def qdag(a):
    out = -a
    out[..., 0] = a[..., 0]
    return out


def qnorm(a):
    return np.sqrt(np.sum(a * a, axis=-1))


def sh(A, d, n):
    """field at x + n*e_d (A has lattice axes 0,1,2 = t,x,y followed by the quaternion axis)."""
    return np.roll(A, -n, axis=d)


def _qmat(q):
    s = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
    return q[0] * np.eye(2) + 1j * sum(q[k + 1] * s[k] for k in range(3))


def selftest():
    rng = np.random.default_rng(1)
    a = rng.normal(size=4); a /= np.linalg.norm(a)
    b = rng.normal(size=4); b /= np.linalg.norm(b)
    assert np.allclose(_qmat(qmul(a, b)), _qmat(a) @ _qmat(b))
    assert np.allclose(_qmat(qdag(a)), _qmat(a).conj().T)
    assert np.isclose(np.trace(_qmat(a)).real / 2, a[0])
    return True


# ---------------- lattice ----------------
class Lattice:
    def __init__(self, T, L, beta, seed=0):
        self.T, self.L, self.beta = T, L, beta
        self.rng = np.random.default_rng(seed)
        U = np.zeros((3, T, L, L, 4)); U[..., 0] = 1.0      # cold start
        self.U = U
        t, x, y = np.meshgrid(np.arange(T), np.arange(L), np.arange(L), indexing="ij")
        self.par = ((t + x + y) % 2).astype(bool)

    def staple(self, mu):
        """V such that the local action of U_mu(x) is -(beta/2) tr(U_mu(x) V)."""
        U = self.U
        V = np.zeros_like(U[0])
        for nu in range(3):
            if nu == mu:
                continue
            # forward: U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
            V += qmul(qmul(sh(U[nu], mu, 1), qdag(sh(U[mu], nu, 1))), qdag(U[nu]))
            # backward: U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu)
            Unm = sh(U[nu], nu, -1)
            V += qmul(qmul(qdag(sh(Unm, mu, 1)), qdag(sh(U[mu], nu, -1))), Unm)
        return V

    def _kp(self, a):
        """Kennedy-Pendleton: x0 with density sqrt(1-x0^2) exp(a x0), vectorised with retries."""
        n = a.size
        x0 = np.empty(n)
        todo = np.arange(n)
        rng = self.rng
        while todo.size:
            aa = a[todo]
            r1 = 1.0 - rng.random(todo.size); r2 = rng.random(todo.size); r3 = 1.0 - rng.random(todo.size)
            r4 = rng.random(todo.size)
            lam2 = -(np.log(r1) + np.cos(2 * np.pi * r2) ** 2 * np.log(r3)) / (2 * aa)
            ok = r4 * r4 <= 1.0 - lam2
            x0[todo[ok]] = 1.0 - 2.0 * lam2[ok]
            todo = todo[~ok]
        return x0

    def _random_vec(self, x0):
        r = np.sqrt(np.clip(1.0 - x0 * x0, 0, None))
        ct = 2 * self.rng.random(x0.size) - 1
        ph = 2 * np.pi * self.rng.random(x0.size)
        st = np.sqrt(1 - ct * ct)
        return np.stack([x0, r * st * np.cos(ph), r * st * np.sin(ph), r * ct], axis=-1)

    def heatbath(self):
        for mu in range(3):
            for p in (False, True):
                m = self.par == p
                V = self.staple(mu)[m]
                k = qnorm(V)
                Vt = V / k[:, None]
                X = self._random_vec(self._kp(self.beta * k))
                self.U[mu][m] = qmul(X, qdag(Vt))

    def overrelax(self):
        for mu in range(3):
            for p in (False, True):
                m = self.par == p
                V = self.staple(mu)[m]
                Vd = qdag(V / qnorm(V)[:, None])
                self.U[mu][m] = qmul(qmul(Vd, qdag(self.U[mu][m])), Vd)

    def reunit(self):
        self.U /= qnorm(self.U)[..., None]

    def plaquette(self):
        U = self.U
        s = 0.0
        for mu in range(3):
            for nu in range(mu + 1, 3):
                P = qmul(qmul(U[mu], sh(U[nu], mu, 1)), qmul(qdag(sh(U[mu], nu, 1)), qdag(U[nu])))
                s += P[..., 0].mean()
        return s / 3

    def action_density_slices(self):
        """sum over each time slice of (1 - 1/2 tr U_P) for the 3 plaquette orientations:
        returns (T, 3) [spatial (xy), tx, ty] -- used for the energy check (G4)."""
        U = self.U
        out = np.zeros((self.T, 3))
        k = 0
        for mu, nu in ((1, 2), (0, 1), (0, 2)):
            P = qmul(qmul(U[mu], sh(U[nu], mu, 1)), qmul(qdag(sh(U[mu], nu, 1)), qdag(U[nu])))
            out[:, k] = (1.0 - P[..., 0]).sum(axis=(1, 2))
            k += 1
        return out


# ---------------- measurements ----------------
def ape_smear(Ux, Uy, alpha, n):
    """APE smearing of spatial links inside each time slice (spatial staples only).
    Ux, Uy: arrays (T, L, L, 4). U' = Proj_SU2[ U + alpha * (sum of the two spatial staples) ].
    For SU(2) the projection of a sum of SU(2) matrices is exact: normalise the quaternion."""
    for _ in range(n):
        nx = Ux + alpha * (qmul(qmul(Uy, sh(Ux, 2, 1)), qdag(sh(Uy, 1, 1)))
                           + qmul(qmul(qdag(sh(Uy, 2, -1)), sh(Ux, 2, -1)), sh(sh(Uy, 2, -1), 1, 1)))
        ny = Uy + alpha * (qmul(qmul(Ux, sh(Uy, 1, 1)), qdag(sh(Ux, 2, 1)))
                           + qmul(qmul(qdag(sh(Ux, 1, -1)), sh(Uy, 1, -1)), sh(sh(Ux, 1, -1), 2, 1)))
        Ux = nx / qnorm(nx)[..., None]
        Uy = ny / qnorm(ny)[..., None]
    return Ux, Uy


def lines(Ud, d, lmax):
    """products of l consecutive links in direction d (axis d of the (T,L,L) arrays), l = 1..lmax."""
    out = {1: Ud}
    cur = Ud
    for l in range(2, lmax + 1):
        cur = qmul(cur, sh(Ud, d, l - 1))
        out[l] = cur
    return out


def loopfield(Lx, Ly, lx, ly):
    """1/2 tr of the lx * ly rectangle in the (x,y) plane with corner at each site -> (T,L,L)."""
    A = qmul(Lx[lx], sh(Ly[ly], 1, lx))      # x -> x+lx e_x -> x+lx e_x+ly e_y
    B = qmul(Ly[ly], sh(Lx[lx], 2, ly))      # x -> x+ly e_y -> same corner
    return qmul(A, qdag(B))[..., 0]


def loop(Lx, Ly, lx, ly):
    return loopfield(Lx, Ly, lx, ly).sum(axis=(1, 2))


def polyakov_spatial(Ux, Uy, L):
    """smeared spatial Polyakov loops (winding once around the box), zero momentum, per time slice:
    sum_y P_x(t,y) + sum_x P_y(t,x)  (rotation-symmetric torelon operator)."""
    Px = Ux.copy()
    for k in range(1, L):
        Px = qmul(Px, sh(Ux, 1, k))
    Py = Uy.copy()
    for k in range(1, L):
        Py = qmul(Py, sh(Uy, 2, k))
    return Px[:, 0, :, 0].sum(axis=1) + Py[:, :, 0, 0].sum(axis=1)


def glue_ops(lat, levels, alpha, sq_sizes, rect_sizes):
    """Operators per time slice.
    Returns (ops, pops): ops[name] = (T,) zero-momentum sums; pops[name] = (T,2) complex, momentum
    p = 2 pi/L along x and along y (A1 operators only, for the dispersion check E(p)^2 = m^2 + p^2).
    A1 (0++): squares l*l and rectangles 2l*l + l*2l; E (2++, B1 irrep): 2l*l - l*2l;
    P: smeared spatial Polyakov loop (torelon)."""
    L = lat.L
    Ux, Uy = lat.U[1], lat.U[2]
    ops, pops = {}, {}
    ph = np.exp(2j * np.pi * np.arange(L) / L)
    def addp(name, f):
        pops[name] = np.stack([(f * ph[None, :, None]).sum(axis=(1, 2)),
                               (f * ph[None, None, :]).sum(axis=(1, 2))], axis=-1)
    done = 0
    lmax = max(max(sq_sizes), 2 * max(rect_sizes))
    for n in levels:
        Ux, Uy = ape_smear(Ux, Uy, alpha, n - done)
        done = n
        Lx = lines(Ux, 1, lmax); Ly = lines(Uy, 2, lmax)
        for l in sq_sizes:
            f = loopfield(Lx, Ly, l, l)
            ops[f"A1_sq{l}_s{n}"] = f.sum(axis=(1, 2)); addp(f"A1_sq{l}_s{n}", f)
        for l in rect_sizes:
            h = loopfield(Lx, Ly, 2 * l, l)
            v = loopfield(Lx, Ly, l, 2 * l)
            ops[f"A1_re{l}_s{n}"] = (h + v).sum(axis=(1, 2)); addp(f"A1_re{l}_s{n}", h + v)
            ops[f"E_re{l}_s{n}"] = (h - v).sum(axis=(1, 2))
        ops[f"P_s{n}"] = polyakov_spatial(Ux, Uy, L)
    return ops, pops


def wilson_loops(lat, Rmax):
    """unsmeared planar Wilson loops W(R,T), averaged over all sites and the 6 ordered planes, R,T = 1..Rmax."""
    U = lat.U
    W = np.zeros((Rmax + 1, Rmax + 1))
    Ls = [lines_numpy(U[mu], mu, Rmax) for mu in range(3)]   # numpy version (any direction)
    for mu in range(3):
        for nu in range(3):
            if nu == mu:
                continue
            for r in range(1, Rmax + 1):
                for t in range(1, Rmax + 1):
                    A = qmul(Ls[mu][r], sh(Ls[nu][t], mu, r))
                    B = qmul(Ls[nu][t], sh(Ls[mu][r], nu, t))
                    W[r, t] += qmul(A, qdag(B))[..., 0].mean()
    return W / 6


# ---------------- numba kernels (same update, ~10x faster) ----------------
import numba as nb


@nb.njit(cache=True)
def _nb_seed(s):
    np.random.seed(s)


@nb.njit(inline="always")
def _m(a0, a1, a2, a3, b0, b1, b2, b3):
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + b0 * a1 - (a2 * b3 - a3 * b2),
            a0 * b2 + b0 * a2 - (a3 * b1 - a1 * b3),
            a0 * b3 + b0 * a3 - (a1 * b2 - a2 * b1))


@nb.njit(cache=True)
def _sweep(U, beta, fw, bw, par, mode):
    """U: (3, V, 4) flat view; fw/bw: (3, V) neighbour tables; par: (V,) site parity.
    mode 0 = Kennedy-Pendleton heat-bath, 1 = over-relaxation. Checkerboard order."""
    V = U.shape[1]
    for mu in range(3):
        for p in range(2):
            for s in range(V):
                if par[s] != p:
                    continue
                v0 = 0.0; v1 = 0.0; v2 = 0.0; v3 = 0.0
                smu = fw[mu, s]
                for nu in range(3):
                    if nu == mu:
                        continue
                    snu = fw[nu, s]
                    # forward: U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
                    a = U[nu, smu]; b = U[mu, snu]; d = U[nu, s]
                    r = _m(a[0], a[1], a[2], a[3], b[0], -b[1], -b[2], -b[3])
                    r = _m(r[0], r[1], r[2], r[3], d[0], -d[1], -d[2], -d[3])
                    v0 += r[0]; v1 += r[1]; v2 += r[2]; v3 += r[3]
                    # backward: U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu)
                    w = bw[nu, s]
                    z = fw[mu, w]
                    a = U[nu, z]; b = U[mu, w]; d = U[nu, w]
                    r = _m(a[0], -a[1], -a[2], -a[3], b[0], -b[1], -b[2], -b[3])
                    r = _m(r[0], r[1], r[2], r[3], d[0], d[1], d[2], d[3])
                    v0 += r[0]; v1 += r[1]; v2 += r[2]; v3 += r[3]
                k = np.sqrt(v0 * v0 + v1 * v1 + v2 * v2 + v3 * v3)
                v0 /= k; v1 /= k; v2 /= k; v3 /= k
                if mode == 0:
                    aa = beta * k
                    while True:
                        r1 = 1.0 - np.random.random(); r2 = np.random.random()
                        r3 = 1.0 - np.random.random(); r4 = np.random.random()
                        c2 = np.cos(2 * np.pi * r2)
                        lam2 = -(np.log(r1) + c2 * c2 * np.log(r3)) / (2 * aa)
                        if r4 * r4 <= 1.0 - lam2:
                            break
                    x0 = 1.0 - 2.0 * lam2
                    rr = np.sqrt(max(0.0, 1.0 - x0 * x0))
                    ct = 2 * np.random.random() - 1
                    ph = 2 * np.pi * np.random.random()
                    st = np.sqrt(1 - ct * ct)
                    n = _m(x0, rr * st * np.cos(ph), rr * st * np.sin(ph), rr * ct, v0, -v1, -v2, -v3)
                else:
                    u = U[mu, s]
                    n = _m(v0, -v1, -v2, -v3, u[0], -u[1], -u[2], -u[3])
                    n = _m(n[0], n[1], n[2], n[3], v0, -v1, -v2, -v3)
                nn = np.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2 + n[3] ** 2)
                U[mu, s, 0] = n[0] / nn
                U[mu, s, 1] = n[1] / nn
                U[mu, s, 2] = n[2] / nn
                U[mu, s, 3] = n[3] / nn


class FastLattice(Lattice):
    """same physics as Lattice, numba sweeps (single thread)."""

    def __init__(self, T, L, beta, seed=0):
        super().__init__(T, L, beta, seed)
        _nb_seed(seed)
        idx = np.arange(T * L * L).reshape(T, L, L)
        self.fw = np.stack([np.roll(idx, -1, axis=d).ravel() for d in range(3)]).astype(np.int64)
        self.bw = np.stack([np.roll(idx, 1, axis=d).ravel() for d in range(3)]).astype(np.int64)
        self.pflat = self.par.ravel().astype(np.int64)

    def _flat(self):
        return self.U.reshape(3, -1, 4)

    def heatbath(self):
        _sweep(self._flat(), self.beta, self.fw, self.bw, self.pflat, 0)

    def overrelax(self):
        _sweep(self._flat(), self.beta, self.fw, self.bw, self.pflat, 1)


@nb.njit(cache=True)
def _ape_step(Ux, Uy, alpha, nx_, ny_):
    """one APE step (same formula as ape_smear) on (T,L,L,4) arrays, writes into nx_, ny_."""
    T, L = Ux.shape[0], Ux.shape[1]
    for t in range(T):
        for x in range(L):
            xp = (x + 1) % L; xm = (x - 1) % L
            for y in range(L):
                yp = (y + 1) % L; ym = (y - 1) % L
                # x-link: Uy(x) Ux(x+y) Uy(x+x)^+  +  Uy(x-y)^+ Ux(x-y) Uy(x-y+x)
                a = Uy[t, x, y]; b = Ux[t, x, yp]; c = Uy[t, xp, y]
                r = _m(a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3])
                s1 = _m(r[0], r[1], r[2], r[3], c[0], -c[1], -c[2], -c[3])
                a = Uy[t, x, ym]; b = Ux[t, x, ym]; c = Uy[t, xp, ym]
                r = _m(a[0], -a[1], -a[2], -a[3], b[0], b[1], b[2], b[3])
                s2 = _m(r[0], r[1], r[2], r[3], c[0], c[1], c[2], c[3])
                u = Ux[t, x, y]
                n0 = u[0] + alpha * (s1[0] + s2[0]); n1 = u[1] + alpha * (s1[1] + s2[1])
                n2 = u[2] + alpha * (s1[2] + s2[2]); n3 = u[3] + alpha * (s1[3] + s2[3])
                nn = np.sqrt(n0 * n0 + n1 * n1 + n2 * n2 + n3 * n3)
                nx_[t, x, y, 0] = n0 / nn; nx_[t, x, y, 1] = n1 / nn
                nx_[t, x, y, 2] = n2 / nn; nx_[t, x, y, 3] = n3 / nn
                # y-link: Ux(x) Uy(x+x) Ux(x+y)^+  +  Ux(x-x)^+ Uy(x-x) Ux(x-x+y)
                a = Ux[t, x, y]; b = Uy[t, xp, y]; c = Ux[t, x, yp]
                r = _m(a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3])
                s1 = _m(r[0], r[1], r[2], r[3], c[0], -c[1], -c[2], -c[3])
                a = Ux[t, xm, y]; b = Uy[t, xm, y]; c = Ux[t, xm, yp]
                r = _m(a[0], -a[1], -a[2], -a[3], b[0], b[1], b[2], b[3])
                s2 = _m(r[0], r[1], r[2], r[3], c[0], c[1], c[2], c[3])
                u = Uy[t, x, y]
                n0 = u[0] + alpha * (s1[0] + s2[0]); n1 = u[1] + alpha * (s1[1] + s2[1])
                n2 = u[2] + alpha * (s1[2] + s2[2]); n3 = u[3] + alpha * (s1[3] + s2[3])
                nn = np.sqrt(n0 * n0 + n1 * n1 + n2 * n2 + n3 * n3)
                ny_[t, x, y, 0] = n0 / nn; ny_[t, x, y, 1] = n1 / nn
                ny_[t, x, y, 2] = n2 / nn; ny_[t, x, y, 3] = n3 / nn


def ape_smear_fast(Ux, Uy, alpha, n):
    Ux = np.ascontiguousarray(Ux); Uy = np.ascontiguousarray(Uy)
    for _ in range(n):
        nx_ = np.empty_like(Ux); ny_ = np.empty_like(Uy)
        _ape_step(Ux, Uy, alpha, nx_, ny_)
        Ux, Uy = nx_, ny_
    return Ux, Uy


@nb.njit(cache=True)
def _poly(Ux, Uy):
    T, L = Ux.shape[0], Ux.shape[1]
    out = np.zeros(T)
    for t in range(T):
        for y in range(L):
            p = (1.0, 0.0, 0.0, 0.0)
            for x in range(L):
                u = Ux[t, x, y]
                p = _m(p[0], p[1], p[2], p[3], u[0], u[1], u[2], u[3])
            out[t] += p[0]
        for x in range(L):
            p = (1.0, 0.0, 0.0, 0.0)
            for y in range(L):
                u = Uy[t, x, y]
                p = _m(p[0], p[1], p[2], p[3], u[0], u[1], u[2], u[3])
            out[t] += p[0]
    return out


# use the fast versions in glue_ops
ape_smear_numpy, polyakov_spatial_numpy = ape_smear, polyakov_spatial
ape_smear = ape_smear_fast
def polyakov_spatial(Ux, Uy, L):
    return _poly(np.ascontiguousarray(Ux), np.ascontiguousarray(Uy))


# numba versions of lines/loopfield (identical results; replaces the numpy ones used by glue_ops)
@nb.njit(cache=True)
def _lines_nb(Ud, d, lmax):
    T, L = Ud.shape[0], Ud.shape[1]
    out = np.empty((lmax + 1, T, L, L, 4))
    for t in range(T):
        for x in range(L):
            for y in range(L):
                p = (Ud[t, x, y, 0], Ud[t, x, y, 1], Ud[t, x, y, 2], Ud[t, x, y, 3])
                out[1, t, x, y, 0] = p[0]; out[1, t, x, y, 1] = p[1]; out[1, t, x, y, 2] = p[2]; out[1, t, x, y, 3] = p[3]
                for l in range(2, lmax + 1):
                    if d == 1:
                        u = Ud[t, (x + l - 1) % L, y]
                    else:
                        u = Ud[t, x, (y + l - 1) % L]
                    p = _m(p[0], p[1], p[2], p[3], u[0], u[1], u[2], u[3])
                    out[l, t, x, y, 0] = p[0]; out[l, t, x, y, 1] = p[1]; out[l, t, x, y, 2] = p[2]; out[l, t, x, y, 3] = p[3]
    return out


@nb.njit(cache=True)
def _loop_nb(Lx, Ly, lx, ly):
    T, L = Lx.shape[1], Lx.shape[2]
    f = np.empty((T, L, L))
    for t in range(T):
        for x in range(L):
            xr = (x + lx) % L
            for y in range(L):
                yu = (y + ly) % L
                a = Lx[lx, t, x, y]; b = Ly[ly, t, xr, y]
                A = _m(a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3])
                a = Ly[ly, t, x, y]; b = Lx[lx, t, x, yu]
                B = _m(a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3])
                f[t, x, y] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2] + A[3] * B[3]   # 1/2 tr(A B^+)
    return f


class _LinesNB:
    def __init__(self, arr):
        self.arr = arr
    def __getitem__(self, l):
        return self.arr[l]


lines_numpy, loopfield_numpy = lines, loopfield
def lines(Ud, d, lmax):
    return _LinesNB(_lines_nb(np.ascontiguousarray(Ud), d, lmax))
def loopfield(Lx, Ly, lx, ly):
    if isinstance(Lx, _LinesNB):
        return _loop_nb(Lx.arr, Ly.arr, lx, ly)
    return loopfield_numpy(Lx, Ly, lx, ly)

"""fw1atomo_core.py - tools of the study fw1atomo_ (docs/FW1_ATOMO.md): FW1 AS A VARIANT ON TRIAL next to the
frozen rule.

* The rules (copied from lab/senzadoppi_real.py, same factors, checked bit for bit in fw1atomo_check.py):
  F0  = the frozen rule of lab/zz_core.py (tick 0: steps x, y, z; tick 1: steps z, y, x; coin theta each tick);
  FW1 = F0 + ONE twisted Wilson per period: around the coin of the EVEN ticks only, for each of the 6 links
        v = +-2 e_j, the L part at x and the R part at x + v rotate into each other by h/4 about the chirality axis
        Y, and the same rotation is removed on the spot (k = 0 untouched); half before the coin, half after
        (mirror order).  h = 0,4 is a DEBT (a regulator chosen by hand, senzadoppi_ amendment A3).
  Coulomb of a point nucleus as the phase exp(-i V a) applied with every coin (as zz_core / the Laboratorio).
* Lattice: periodic box of n^3 sites, site i at x_i = (i - n/2) a + s, nucleus at x = 0.  The offset s per axis
  chooses where the nucleus sits: s = a/2 on the three axes = corner of 8 cells (Laboratorio, BATTIMENTO_1S);
  (a/2, a/2, 0) = middle of an edge (4 cells); (a/2, 0, 0) = middle of a face (2 cells); (0, 0, 0) = on a site
  (cell centre): there V(0) = -Z alpha <1/r>_cell = -Z alpha * 2,380077 / a (average of 1/r over the cell: a
  DEBT, a regularisation chosen by hand, only for that configuration).
* Ritz: filter diagonalisation as in src/atomo/zz.rs::ritz (BATTIMENTO_1S): evolve the start state(s), keep on
  every period (even ticks) the Hann-windowed sums sum_p w(p) e^{i om p T} psi(p) on a grid of om, orthonormalise,
  restrict the period map U2 to that space, diagonalise; every level with its residual |U2 v - lam v|.
* Diagnostics of a level: valley content (Voronoi in k-space: Gamma valley, the 3 N valleys (pi/2,pi/2,0)-type,
  the 4 pi valleys (pi/2,0,0)/(pi/2,pi/2,pi/2)-type, all modulo the copies {0,pi}^3), the 4 colours
  ((i+j)%2, (j+k)%2) and the 8 period-sublattices (i%2, j%2, k%2 at even ticks), spin vector <Sigma> (L and R
  summed), <r>, weight inside r < 1 lambda_C.
Units: hbar = c = m = 1, lengths in lambda_C, theta = a."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import itertools
import numpy as np
import zz_core as zz

BH = zz.BASIS_HALF
TICKS = [((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))]
VEC = [2 * e for e in np.eye(3)] + [-2 * e for e in np.eye(3)]
YCHI = np.kron(zz.SY, zz.I2)
XCHI = np.kron(zz.SX, zz.I2)
CUBE_INV_R = 2.380077363979557          # <1/r> over a unit cube centred at 0 (numerical, 7e-8)
OFFSETS = {"angolo": (0.5, 0.5, 0.5), "spigolo": (0.5, 0.5, 0.0), "faccia": (0.5, 0.0, 0.0), "sito": (0.0, 0.0, 0.0)}
ZA60 = 60 * zz.ALPHA
SCR = "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad/fw1atomo"
os.makedirs(SCR, exist_ok=True)


def axis_shift(blk, j, h):
    lam, P = BH[j]
    c = np.tensordot(P.conj().T, blk, axes=(1, 0))
    for q in range(2):
        c[q] = np.roll(c[q], h * lam[q], axis=j)
    return np.tensordot(P, c, axes=(1, 0))


def rot4(psi, M):
    return np.tensordot(M, psi, axes=(1, 0))


def link(psi, v, a, G):
    ea = np.cos(a) * np.eye(4) - 1j * np.sin(a) * G
    eb = np.cos(a) * np.eye(4) + 1j * np.sin(a) * G
    psi = rot4(psi, eb)
    v = tuple(int(x) for x in v)
    psi[2:] = np.roll(psi[2:], v, axis=(1, 2, 3))
    psi = rot4(psi, ea)
    psi[2:] = np.roll(psi[2:], tuple(-x for x in v), axis=(1, 2, 3))
    return psi


class Rule:
    """name 'F0' (h ignored) or 'FW1' (h); theta = a; V (n,n,n) potential or None."""

    def __init__(self, name, theta, h=0.4, V=None):
        self.name, self.theta, self.h = name, theta, (h if name == "FW1" else 0.0)
        self.C = np.cos(theta) * np.eye(4) - 1j * np.sin(theta) * XCHI
        self.ph = None if V is None else np.exp(-1j * V)        # V given in lattice units (V * a)

    def _coin(self, psi, n):
        w = self.h and n % 2 == 0
        if w:
            for v in VEC:
                psi = link(psi, v, self.h / 4, YCHI)
        psi = rot4(psi, self.C)
        if self.ph is not None:
            psi = psi * self.ph[None]
        if w:
            for v in VEC[::-1]:
                psi = link(psi, v, self.h / 4, YCHI)
        return psi

    def period(self, psi):
        psi = psi.copy()
        for n, (oR, oL) in enumerate(TICKS):
            L, R = psi[:2], psi[2:]
            for j in oR:
                R = axis_shift(R, j, +1)
            for j in oL:
                L = axis_shift(L, j, -1)
            psi = self._coin(np.concatenate([L, R]), n)
        return psi


class Lattice:
    def __init__(self, n, a, pos="angolo", Z=60):
        self.n, self.a, self.pos, self.ZA = n, a, pos, Z * zz.ALPHA
        s = np.array(OFFSETS[pos]) * a
        i = np.arange(n)
        xs = [(i - n // 2) * a + s[d] for d in range(3)]
        self.X, self.Y, self.Z = np.meshgrid(*xs, indexing="ij")
        self.r = np.sqrt(self.X**2 + self.Y**2 + self.Z**2)
        rr = np.where(self.r > 1e-12, self.r, 1.0)
        inv = np.where(self.r > 1e-12, 1.0 / rr, CUBE_INV_R / a)
        self.V = -self.ZA * inv * a                       # lattice units: phase per tick
        I, J, K = np.meshgrid(i, i, i, indexing="ij")
        self.colour = ((I + J) % 2) * 2 + ((J + K) % 2)
        self.cls = (I % 2) * 4 + (J % 2) * 2 + (K % 2)

    def rule(self, name, h=0.4, coulomb=True):
        return Rule(name, self.a, h, self.V if coulomb else None)

    # ---------------- start states (TEST states, not inputs of the rule) ----------------
    def dirac1s(self, spin=(1, 0)):
        """continuum Dirac 1s (Sommerfeld-Dirac formula: CHECK/TEST state only), spin spinor chi"""
        ZA, r = self.ZA, np.where(self.r > 1e-12, self.r, 0.25 * self.a)
        gam = np.sqrt(1 - ZA**2)
        g = r**(gam - 1) * np.exp(-ZA * r)
        f = -np.sqrt((1 - gam) / (1 + gam)) * g
        nx, ny, nz = self.X / r, self.Y / r, self.Z / r
        chi = np.asarray(spin, complex) / np.linalg.norm(spin)
        u = np.array([g * chi[0], g * chi[1]])
        sn = np.array([nz * chi[0] + (nx - 1j * ny) * chi[1], (nx + 1j * ny) * chi[0] - nz * chi[1]])
        dd = 1j * f * sn
        psi = np.concatenate([(u - dd) / np.sqrt(2), (u + dd) / np.sqrt(2)]).astype(complex)
        return psi / np.linalg.norm(psi)

    def radial1s(self, spin_up=True):
        """1s (kappa = -1, m = +1/2) of the RADIAL form of the same rule (zz_core.RadialWalk, same a, Z alpha), put
        on the cube as src/atomo/zz.rs::con_stato does (linear interpolation of r g, r f; Dirac -> chiral)."""
        W = zz.RadialWalk(-1, self.a, 40.0, self.ZA)
        om, u, res = W.eigen(1 - 0.5 * self.ZA**2 * 1.2)
        g, f = W.gf(u)
        r = np.maximum(self.r, 1e-9)
        x = r / self.a - 0.5
        i = np.clip(np.floor(np.maximum(x, 0)).astype(int), 0, len(g) - 2)
        t = np.clip(x - i, 0, 1)
        G = np.real((g[i] * (1 - t) + g[i + 1] * t)) / r
        F = np.real((f[i] * (1 - t) + f[i + 1] * t)) / r
        y00 = 0.5 / np.sqrt(np.pi)
        cz, cx, cy = self.Z / r, self.X / r, self.Y / r
        up = [y00 * G, 0 * G]
        lo = [-y00 * cz * F, (-y00 * cx - 1j * y00 * cy) * F]
        d = [up[0], up[1], 1j * lo[0], 1j * lo[1]]
        if not spin_up:
            raise NotImplementedError
        s = 1 / np.sqrt(2)
        psi = np.array([(d[0] - d[2]) * s, (d[1] - d[3]) * s, (d[0] + d[2]) * s, (d[1] + d[3]) * s], complex)
        return psi / np.linalg.norm(psi), float(om)

    def blind(self, seed=20261013, width=2):
        """BLIND start: random 4-spinors on the (2 width)^3 sites nearest to the nucleus (all momenta, all valleys,
        both energy signs, all colours): no projector, no choice of k = 0."""
        rng = np.random.default_rng(seed)
        psi = np.zeros((4,) + self.r.shape, complex)
        c = self.n // 2
        sl = tuple(slice(c - width, c + width) for _ in range(3))
        sh = psi[(slice(None),) + sl].shape
        psi[(slice(None),) + sl] = rng.normal(size=sh) + 1j * rng.normal(size=sh)
        return psi / np.linalg.norm(psi)

    # ---------------- diagnostics ----------------
    def diag(self, v):
        v = v / np.linalg.norm(v)
        rho = np.sum(np.abs(v)**2, axis=0)
        out = {"r_mean": float(np.sum(rho * self.r)), "w_r_lt_1": float(rho[self.r < 1.0].sum())}
        out["colours"] = [float(rho[self.colour == c].sum()) for c in range(4)]
        out["classes"] = [float(rho[self.cls == c].sum()) for c in range(8)]
        S = []
        for P in zz.PAULI:
            s = 0.0
            for blk in (v[:2], v[2:]):
                s += np.real(np.sum(np.conj(blk) * np.tensordot(P, blk, axes=(1, 0))))
            S.append(float(s))
        out["spin"] = S
        out.update(valleys(v))
        return out


_VAL_CACHE = {}


def _valley_labels(n):
    if n in _VAL_CACHE:
        return _VAL_CACHE[n]
    k = np.fft.fftfreq(n) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    h = np.pi / 2
    groups = [("Gamma", [(0, 0, 0)]), ("N", [(h, h, 0), (h, 0, h), (0, h, h)]),
              ("pi", [(h, 0, 0), (0, h, 0), (0, 0, h), (h, h, h)])]
    best = np.full(KX.shape, np.inf)
    lab = np.zeros(KX.shape, int)
    for gi, (_, cs) in enumerate(groups):
        for c in cs:
            for sgn in itertools.product((1, -1), repeat=3):
                for G in itertools.product((0, np.pi), repeat=3):
                    d2 = 0
                    for Kc, cc, gg, sg in ((KX, c[0], G[0], sgn[0]), (KY, c[1], G[1], sgn[1]), (KZ, c[2], G[2], sgn[2])):
                        d2 = d2 + ((Kc - sg * cc - gg + np.pi) % (2 * np.pi) - np.pi)**2
                    m = d2 < best
                    best[m] = d2[m]
                    lab[m] = gi
    _VAL_CACHE[n] = (lab, best)
    return lab, best


def valleys(v):
    n = v.shape[1]
    F = np.fft.fftn(v, axes=(1, 2, 3))
    P = np.sum(np.abs(F)**2, axis=0)
    P /= P.sum()
    lab, d2 = _valley_labels(n)
    return {"val_Gamma": float(P[lab == 0].sum()), "val_N": float(P[lab == 1].sum()),
            "val_pi": float(P[lab == 2].sum()), "val_far": float(P[d2 > 0.6**2].sum())}


def ritz(rule, starts, a, lo, hi, dw, nper, ck=None, tmax=None, t0=None, sub=1):
    """filter diagonalisation.  starts: list of start states.  Every `sub` periods the Hann-windowed sums
    sum_p w(p) e^{i om p T} psi(p) are kept on the grid om = lo, lo + dw, ..., hi; orthonormalised (Gram matrix);
    U2 restricted; eigenvectors with their exact residual |U2 v - lam v| (one more period each).
    Returns list of dicts (omega, resid, v, w_start).  With ck (path prefix) and tmax (s): checkpoints and
    raises SystemExit(3) when the time budget is over (rerun to continue)."""
    import time, json
    t0 = t0 or time.time()
    T = 2 * a
    oms = np.arange(lo, hi + 1e-12, dw)
    K1 = len(oms)
    w = np.hanning(nper + 1)
    shp = starts[0].shape
    st = {"s": 0, "p": 0}
    if ck and os.path.exists(ck + "_st.json"):
        st = json.load(open(ck + "_st.json"))
        acc = np.load(ck + "_acc.npy")
        psi = np.load(ck + "_psi.npy")
    else:
        acc = np.zeros((len(starts) * K1,) + shp, complex)
        psi = None
    while st["s"] < len(starts):
        if st["p"] == 0:
            psi = starts[st["s"]].copy()
        base = st["s"] * K1
        while st["p"] <= nper:
            if tmax and time.time() - t0 > tmax:
                np.save(ck + "_acc.npy", acc)
                np.save(ck + "_psi.npy", psi)
                json.dump(st, open(ck + "_st.json", "w"))
                print("checkpoint", ck, st, flush=True)
                raise SystemExit(3)
            p = st["p"]
            if p % sub == 0:
                ph = np.exp(1j * oms * p * T) * w[p]
                for q in range(K1):
                    acc[base + q] += ph[q] * psi
            psi = rule.period(psi)
            st["p"] += 1
        st["s"] += 1
        st["p"] = 0
    A = acc.reshape(acc.shape[0], -1)
    Gm = A.conj() @ A.T
    ev, U = np.linalg.eigh(Gm)
    keep = ev > ev.max() * 1e-13
    Cq = (U[:, keep] / np.sqrt(ev[keep])).T                # Q = Cq @ A, rows orthonormal
    Q = Cq @ A
    del acc, A
    m = Q.shape[0]
    M = np.zeros((m, m), complex)
    for j in range(m):
        M[:, j] = Q.conj() @ rule.period(Q[j].reshape(shp)).ravel()
    lam, C = np.linalg.eig(M)
    out = []
    for k in range(len(lam)):
        v = C[:, k] @ Q
        v /= np.linalg.norm(v)
        Uv = rule.period(v.reshape(shp)).ravel()
        l = np.vdot(v, Uv)
        res = float(np.linalg.norm(Uv - l * v))
        om = float(-np.angle(l) / T)
        wts = [float(abs(np.vdot(s_.ravel(), v))**2) for s_ in starts]
        out.append({"omega": om, "resid": res, "v": v.reshape(shp), "w_start": wts})
    out.sort(key=lambda d: d["omega"])
    if ck:
        for f in ("_acc.npy", "_psi.npy", "_st.json"):
            if os.path.exists(ck + f):
                os.remove(ck + f)
    return out


def survival(rule, psi0, nper, every=1):
    psi = psi0.copy()
    c = np.zeros(nper + 1, complex)
    c[0] = 1
    for p in range(1, nper + 1):
        psi = rule.period(psi)
        c[p] = np.vdot(psi0, psi)
    return c


def lines(c, T, Emin=0.6, Emax=1.05, dE=2e-5, thr=0.02):
    t = np.arange(len(c)) * T
    w = np.hanning(len(c))
    Es = np.arange(Emin, Emax, dE)
    S = np.abs(np.exp(1j * np.outer(Es, t)) @ (c * w)) / np.sum(w)
    pk = [i for i in range(1, len(Es) - 1) if S[i] >= S[i - 1] and S[i] >= S[i + 1] and S[i] >= thr]
    return [(float(Es[i]), float(S[i])) for i in pk]

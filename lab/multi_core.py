"""MULTI-PHOTON PILOT RULE — the single-photon rule of REGOLA_FOTONE.md extended to n photons.
(2-D lab: x transverse, z along the beams, λ = 1.)  The rule is defined ONCE, here; every experiment
script (multi_*.py) only builds a geometry (sources, objects, detectors) and an initial state.

WAVE BALLS (unchanged from the single-photon rule)
  * fly straight, are re-emitted as regular fans (origins every DS, all forward directions) when they
    touch an object; amplitude at an origin = incoming field there (read from the balls that arrive)
    × the object's local matrix.  An object is a Jones-type matrix acting on the ball's internal
    vector = (path label) ⊗ (polarization), times a local phase ramp that sets the outgoing tilt
    (beam splitter BS, polarizing beam splitter PBS, mirror/prism, wave plate, phase shifter).
  * a set of wave balls that left one source beam and took one route through the objects is a
    CHANNEL: (label, tilt, fan of origins).  Channels with the same label interfere, channels with
    different labels do not (exactly like orthogonal polarizations in the single-photon rule).

THE n-PHOTON WAVE (what the wave balls represent)
  Ψ(x1..xn) = Σ_k c_k  ψ_k1(x1) ⊗ ψ_k2(x2) ⊗ … ⊗ ψ_kn(xn)          (a SMALL sum of products)
  ψ_ki = the single-photon wave of photon i in term k = Σ over the channels descending from its source
  beam, each with internal vector J_channel · pol_ki, and a "tag" vector τ_ki (spectral / timing
  mode, untouched by the objects).  Unentangled photons: 1 term.  Entangled pair: Schmidt terms
  (2 for |Φ+>).  Identical photons in the same arm: the terms are symmetrized (bosons).
  Cost = (number of terms) × (number of photons) single-photon waves, never a grid in 2n dimensions.

PHOTON BALLS (the only new thing; reduces exactly to the single-photon rule for one photon)
  At each tick every photon ball i reads its own wave balls at its own position x_i — for every
  term k — and the GLOBAL configuration enters only through the other photons' factors:
      Q^j_kl = flux overlap of photon j between terms k and l, read where photon j IS
      X^i_kl = the same for the transverse current of photon i
      dx_i/dz = Re Σ_kl c_k* c_l X^i_kl Π_{j≠i} Q^j_kl  /  Re Σ_kl c_k* c_l Π_j Q^j_kl
  For a product state this is exactly j_x/j_z of the photon's own wave (single-photon rule).
  For an entangled state the velocity of i depends on where the OTHER photons are: this is the one
  nonlocal ingredient, required by Bell's theorem.  It is exactly equivariant (the transition
  currents are divergence-free), so photons born with the Born density stay Born-distributed.
  * birth: the source places all n photon balls at once with density ρ (the source event);
  * non-absorbing objects (BS, PBS, mirrors): no coin; the first step behind the plane (length DELTA)
    is taken along the outgoing current (implicit Euler, as in the single-photon rule);
  * detector plane: the ball stops there; its click = the output port / position where it lands; a
    detected ball keeps its (frozen) position in the guidance of the others.
  * absorbing objects (not used in these experiments): local coin with p = Φ_out/Φ_in of the
    photon's own conditional flux (same formula with X replaced by the flux before/after).

ACCELERATION (does not change the rule): the neighbour sum of a channel is the same function for
every photon that reads it, so it is tabulated once on a grid that follows the beam (Hermite in x,
linear in z).  `selftest()` checks it against the direct ball sum.
"""
import os, sys, time, hashlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from optics import K                        # noqa: E402  (also reconfigures stdout to utf-8)
from photon_rule import kernel, Fan, field_lattice  # noqa: E402

DS = 0.1            # origin spacing on smooth objects (REGOLA §3.3)
NDIR = 31           # directions per origin (dθ = 0.1, REGOLA variant H) — only for counting balls
DELTA = 0.5         # first step behind an object (= 5 DS)
DZMIN, DZMAX, GROW = 0.05, 2.0, 0.2
DU, DZG = 0.25, 0.5  # cache grid (transverse, along z)
AMP_CUT = 1e-6
ALPHA = 0.25         # standard deflection of BS / PBS outputs (rad)
CACHE = os.path.join(HERE, "out", "multi_cache")
EXACT = False        # True: photon balls read the wave balls directly (no table)

T50 = 1 / np.sqrt(2)
R50 = 1j / np.sqrt(2)
I2 = np.eye(2, dtype=complex)


def lin(a):
    return np.array([np.cos(a), np.sin(a)], complex)


H, V = lin(0.0), lin(np.pi / 2)
RC = np.array([1, 1j]) / np.sqrt(2)      # right circular
LC = np.array([1, -1j]) / np.sqrt(2)


# ─────────────────────────────── direct ball sum ───────────────────────────────────────
def fan_sum(x, z, s, a, ze, chunk=768):
    """Field and gradient of a lattice fan (origins s, amplitudes a incl. ds) at points (x, z>ze).
    Equal to the hat-window neighbour reading of REGOLA §3.2 (weights sum to 1)."""
    x = np.atleast_1d(np.asarray(x, float))
    z = np.broadcast_to(np.asarray(z, float), x.shape)
    out = np.zeros((3, x.size), complex)
    for i in range(0, x.size, chunk):
        xx, zz = x[i:i + chunk], z[i:i + chunk]
        G, Gx, Gz = kernel(xx[:, None] - s[None, :], (zz - ze)[:, None])
        out[0, i:i + chunk] = G @ a
        out[1, i:i + chunk] = Gx @ a
        out[2, i:i + chunk] = Gz @ a
    return out


# ─────────────────────────────── channels (sets of wave balls) ─────────────────────────
class Channel:
    """Wave balls of one route: fan emitted at plane ze from origins s with scalar amplitudes a
    (ds included), label, nominal tilt theta, accumulated Jones matrix, source beam."""

    def __init__(self, label, theta, ze, s, a, jones, src):
        self.label, self.theta, self.ze = label, float(theta), float(ze)
        self.s, self.a = np.asarray(s, float), np.asarray(a, complex)
        self.J, self.src = np.asarray(jones, complex), src
        self.z_end = None
        p = np.abs(self.a) ** 2
        self.x0 = (p * self.s).sum() / p.sum()
        self.sig0 = np.sqrt((p * (self.s - self.x0) ** 2).sum() / p.sum())
        # angular spectrum of the origins -> mean direction and angular spread
        n = 1 << int(np.ceil(np.log2(8 * self.s.size)))
        A = np.fft.fft(self.a, n)
        kx = 2 * np.pi * np.fft.fftfreq(n, DS)
        ok = np.abs(kx) < K * 0.999
        w = np.abs(A[ok]) ** 2
        tn = kx[ok] / np.sqrt(K ** 2 - kx[ok] ** 2)
        self.tan = (w * tn).sum() / w.sum()
        self.sig_t = np.sqrt((w * (tn - self.tan) ** 2).sum() / w.sum())

    @property
    def n_origins(self):
        return self.s.size

    @property
    def n_balls(self):
        return self.s.size * NDIR

    def center(self, z):
        return self.x0 + self.tan * (z - self.ze)

    def half_width(self, z):
        return 7 * (self.sig0 + abs(z - self.ze) * self.sig_t) + 2

    def exact(self, x, z):
        return fan_sum(x, z, self.s, self.a, self.ze)

    def key(self):
        h = hashlib.sha1()
        for v in (self.s, self.a, np.array([self.theta, self.ze, self.z_end, DU, DZG, DELTA])):
            h.update(np.ascontiguousarray(v).tobytes())
        return h.hexdigest()[:20]

    def build(self, memo):
        k = self.key()
        if k in memo:
            self.__dict__.update({n: memo[k][n] for n in ("zg", "Uh", "E", "Eu", "Ez")})
            return 0.0
        fn = os.path.join(CACHE, k + ".npz")
        t0 = time.time()
        if os.path.exists(fn):
            d = np.load(fn)
            self.zg, self.Uh, self.E, self.Eu, self.Ez = d["zg"], float(d["Uh"]), d["E"], d["Eu"], d["Ez"]
        else:
            nz = max(2, int(np.ceil((self.z_end - self.ze - DELTA) / DZG)) + 1)
            self.zg = np.linspace(self.ze + DELTA, self.z_end, nz)
            self.Uh = DU * np.ceil(self.half_width(self.z_end) / DU)
            u = np.arange(-self.Uh, self.Uh + DU / 2, DU)
            X = self.center(self.zg)[:, None] + u[None, :]
            Z = np.broadcast_to(self.zg[:, None], X.shape)
            f, fx, fz = self.exact(X.ravel(), Z.ravel())
            f, fx, fz = (q.reshape(X.shape) for q in (f, fx, fz))
            st, ct = np.sin(self.theta), np.cos(self.theta)
            e = np.exp(-1j * K * (st * X + ct * (Z - self.ze)))
            self.E = f * e
            self.Eu = fx * e - 1j * K * st * self.E
            self.Ez = fz * e - 1j * K * ct * self.E
            os.makedirs(CACHE, exist_ok=True)
            np.savez(fn, zg=self.zg, Uh=self.Uh, E=self.E, Eu=self.Eu, Ez=self.Ez)
        memo[k] = {n: getattr(self, n) for n in ("zg", "Uh", "E", "Eu", "Ez")}
        return time.time() - t0

    def eval(self, x, z):
        """φ, ∂xφ, ∂zφ at points x (array) and z (scalar or array) — from the table."""
        if EXACT:
            return self.exact(x, z)
        x = np.asarray(x, float)
        z = np.broadcast_to(np.asarray(z, float), x.shape)
        nz, nu = self.E.shape
        fu = (x - self.center(z) + self.Uh) / DU
        fz = (z - self.zg[0]) / (self.zg[1] - self.zg[0])
        inside = (fu >= 0) & (fu <= nu - 1)
        j = np.clip(np.floor(fu).astype(int), 0, nu - 2)
        t = np.clip(fu - j, 0, 1)
        iz = np.clip(np.floor(fz).astype(int), 0, nz - 2)
        sz = np.clip(fz - iz, 0, 1)
        t2, t3 = t * t, t * t * t
        h00, h10, h01, h11 = 2 * t3 - 3 * t2 + 1, t3 - 2 * t2 + t, -2 * t3 + 3 * t2, t3 - t2
        d00, d10, d01, d11 = 6 * t2 - 6 * t, 3 * t2 - 4 * t + 1, -6 * t2 + 6 * t, 3 * t2 - 2 * t
        val = der = ezz = 0
        for row, wz in ((iz, 1 - sz), (iz + 1, sz)):
            E0, E1 = self.E[row, j], self.E[row, j + 1]
            D0, D1 = self.Eu[row, j] * DU, self.Eu[row, j + 1] * DU
            val = val + wz * (h00 * E0 + h10 * D0 + h01 * E1 + h11 * D1)
            der = der + wz * (d00 * E0 + d10 * D0 + d01 * E1 + d11 * D1) / DU
            ezz = ezz + wz * ((1 - t) * self.Ez[row, j] + t * self.Ez[row, j + 1])
        st, ct = np.sin(self.theta), np.cos(self.theta)
        e = np.exp(1j * K * (st * x + ct * (z - self.ze))) * inside
        return np.array([e * val, e * (der + 1j * K * st * val), e * (ezz + 1j * K * ct * val)])


def tilted_gaussian(s, x0, theta, w0, zw):
    """Paraxial 2-D Gaussian beam on the plane z = 0: axis through (x0, 0) at angle theta, waist w0
    at z = zw (measured along z)."""
    zR = K * w0 ** 2 / 2
    xi = (s - x0) * np.cos(theta)
    eta = (s - x0) * np.sin(theta)
    q = eta - zw / np.cos(theta) - 1j * zR          # e^{+ikz} convention: 1/q = 1/R + i/zR·...
    return np.sqrt(-1j * zR / q) * np.exp(1j * K * eta + 1j * K * xi ** 2 / (2 * q))


# ─────────────────────────────── arms (one 2-D plane each) ─────────────────────────────
class Arm:
    """One 2-D plane with its source beams, objects and detector plane.  Photons that share an arm
    share its wave balls (channels)."""

    def __init__(self, name, z_det):
        self.name, self.z_det = name, float(z_det)
        self.sources, self.objects, self.channels = {}, [], []

    def source(self, beam, label, x0=0.0, theta=0.0, w0=5.0, zw=0.0):
        zR = K * w0 ** 2 / 2
        wpl = w0 * np.sqrt(1 + (zw / zR) ** 2)
        half = 5 * wpl / np.cos(theta) + 1
        n = int(round(2 * half / DS))
        s = x0 - half + (np.arange(n) + 0.5) * DS
        g = tilted_gaussian(s, x0, theta, w0, zw)
        keep = np.abs(g) > AMP_CUT * np.abs(g).max()
        s, g = s[keep], g[keep]
        g = g / np.sqrt((np.abs(g) ** 2).sum() * DS * np.cos(theta))    # unit flux
        self.sources[beam] = Channel(label, theta, 0.0, s, g * DS, I2, beam)
        return self

    def obj(self, z, rules, name=""):
        """rules: {in_label: [(out_label, out_theta, M), ...]}, M a scalar or 2x2 Jones matrix.
        Labels not in rules fly through untouched (the object does not cover them)."""
        self.objects.append((float(z), rules, name))
        return self

    def build(self, memo):
        live = list(self.sources.values())
        self.channels = []
        for z, rules, _ in sorted(self.objects, key=lambda o: o[0]):
            nxt = []
            for ch in live:
                if ch.label not in rules:
                    nxt.append(ch); continue
                ch.z_end = z
                self.channels.append(ch)
                c, hw = ch.center(z), ch.half_width(z)
                s = np.arange(c - hw, c + hw, DS) + DS / 2
                f = ch.exact(s, z)[0]
                keep = np.abs(f) > AMP_CUT * np.abs(f).max()
                s, f = s[keep], f[keep]
                for lab, th, M in rules[ch.label]:
                    M = np.asarray(M, complex)
                    Mj = M * I2 if M.ndim == 0 else M
                    ramp = np.exp(1j * K * (np.sin(th) - np.sin(ch.theta)) * s)
                    flux = np.sqrt(np.cos(ch.theta) / np.cos(th))
                    a = f * ramp * flux * DS
                    child = Channel(lab, th, z, s, a, Mj @ ch.J, ch.src)
                    nxt.append(child)
            live = nxt
        for ch in live:
            ch.z_end = self.z_det
            self.channels.append(ch)
        t = 0.0
        for ch in self.channels:
            t += ch.build(memo)
        self.labels = sorted({ch.label for ch in self.channels})
        self.lab_idx = {l: i for i, l in enumerate(self.labels)}
        self.out = [ch for ch in self.channels if ch.z_end == self.z_det]
        return t

    def planes(self):
        return sorted({o[0] for o in self.objects})

    def port(self, x):
        """Output port where a ball lands: label of the nearest outgoing beam centre."""
        labs = [ch.label for ch in self.out]
        cen = np.array([ch.center(self.z_det) for ch in self.out])
        return np.array(labs, dtype=object)[np.argmin(np.abs(np.asarray(x)[:, None] - cen[None, :]), 1)]

    @property
    def n_balls(self):
        return sum(ch.n_balls for ch in self.channels)


# ─────────────────────────────── the n-photon state and its guidance ───────────────────
class System:
    """photons: list of arm names (photon i lives in arm photons[i]).
    terms: list of (c_k, [(beam, pol, tag) for each photon]); pol: Jones 2-vector; tag: any vector
    (spectral/timing mode) or None."""

    _memo = {}

    def __init__(self, arms, photons, terms):
        self.arms = {a.name: a for a in arms}
        self.ph = list(photons)
        self.n = len(self.ph)
        self.c = np.array([t[0] for t in terms], complex)
        self.K = len(terms)
        t0 = time.time()
        self.build_time = sum(a.build(System._memo) for a in arms)
        self.spec = [[terms[k][1][i] for k in range(self.K)] for i in range(self.n)]
        self.P = np.conj(self.c)[:, None] * self.c[None, :]
        self.T = []
        for i in range(self.n):
            tags = [np.atleast_1d(np.asarray(1.0 if s[2] is None else s[2], complex)) for s in self.spec[i]]
            self.T.append(np.array([[np.vdot(a, b) for b in tags] for a in tags]))
        self.zdet = np.array([self.arms[a].z_det for a in self.ph])
        planes = set()
        for a in arms:
            planes |= set(a.planes())
        self.planes = sorted(planes)
        self.wall = time.time() - t0

    # -- single-photon waves of photon i for all terms at positions x (N,), plane z (scalar)
    def waves(self, i, x, z):
        arm = self.arms[self.ph[i]]
        N, L = x.size, len(arm.labels)
        psi = np.zeros((3, N, self.K, 2 * L), complex)
        for ch in arm.channels:
            if not (ch.ze + DELTA - 1e-7 <= z <= ch.z_end + 1e-7):
                continue
            ks = [k for k in range(self.K) if self.spec[i][k][0] == ch.src]
            if not ks:
                continue
            f = ch.eval(x, z)                                   # (3, N)
            li = arm.lab_idx[ch.label]
            for k in ks:
                v = ch.J @ np.asarray(self.spec[i][k][1], complex)
                psi[:, :, k, 2 * li:2 * li + 2] += f[:, :, None] * v[None, None, :]
        return psi

    def _pairs(self, i, z):
        """Active channels of photon i at plane z and the constant tensor
        M[p, k, l] = <J_c pol_ik | J_c' pol_il> for every pair p = (c, c') of channels with the same label
        (channels with different labels never interfere)."""
        arm = self.arms[self.ph[i]]
        act = [ch for ch in arm.channels if ch.ze + DELTA - 1e-7 <= z <= ch.z_end + 1e-7
               and any(self.spec[i][k][0] == ch.src for k in range(self.K))]
        key = (i, tuple(id(ch) for ch in act))
        if key not in self._pc:
            vec = np.zeros((len(act), self.K, 2), complex)
            for a, ch in enumerate(act):
                for k in range(self.K):
                    if self.spec[i][k][0] == ch.src:
                        vec[a, k] = ch.J @ np.asarray(self.spec[i][k][1], complex)
            pairs, M = [], []
            for a, ca in enumerate(act):
                for b, cb in enumerate(act):
                    if ca.label == cb.label:
                        m = np.conj(vec[a]) @ vec[b].T             # (K, K)
                        if np.abs(m).max() > 0:
                            pairs.append((a, b)); M.append(m.reshape(-1))
            M = np.array(M).reshape(len(pairs), self.K, self.K)
            Mt = np.conj(np.swapaxes(M, 1, 2))
            # Q = -i/(2K) (G - G^H) folded into one matrix: Q = [A, conj A] @ Mq  (tags included)
            Mq = (-0.5j / K) * np.concatenate([M, -Mt]) * self.T[i][None]
            self._pc[key] = (act, np.array(pairs, int).reshape(-1, 2), Mq.reshape(2 * len(pairs), -1))
        return self._pc[key]

    def overlaps(self, i, x, z, want_x=True):
        """Q^i (flux overlap) and X^i (transverse-current overlap), both (N, K, K), tag-weighted.
        G_kl = Σ_pairs conj(φ_c) ∂φ_c' M[p,k,l]: one matrix product per photon."""
        act, pr, M = self._pc_get(i, z)
        N = x.size
        if len(act) == 0:
            Z = np.zeros((N, self.K, self.K), complex)
            return Z, (Z if want_x else None)
        f = np.array([ch.eval(x, z) for ch in act])              # (nch, 3, N)
        fc = np.conj(f[pr[:, 0], 0])                              # (npairs, N)
        A = fc * f[pr[:, 1], 2]
        Q = (np.concatenate([A, np.conj(A)]).T @ M).reshape(N, self.K, self.K)
        if not want_x:
            return Q, None
        A = fc * f[pr[:, 1], 1]
        X = (np.concatenate([A, np.conj(A)]).T @ M).reshape(N, self.K, self.K)
        return Q, X

    def _pc_get(self, i, z):
        if not hasattr(self, "_pc"):
            self._pc = {}
        return self._pairs(i, z)

    def density(self, Xp, zs):
        Qs = [self.overlaps(i, Xp[:, i], zs[i])[0] for i in range(self.n)]
        prod = np.prod(np.stack(Qs), 0)
        return np.real((self.P[None] * prod).sum((1, 2)))

    def slopes(self, Xp, zs, frozen, prev):
        QX = [self.overlaps(i, Xp[:, i], zs[i], not frozen[i]) for i in range(self.n)]
        Qs = [q for q, _ in QX]
        # products of the other photons' factors (prefix/suffix, no division)
        pre = [None] * (self.n + 1); suf = [None] * (self.n + 1)
        pre[0] = self.P[None]
        for i in range(self.n):
            pre[i + 1] = pre[i] * Qs[i]
        suf[self.n] = 1.0
        for i in range(self.n - 1, -1, -1):
            suf[i] = Qs[i] * suf[i + 1]
        rho = np.real(pre[self.n].sum((1, 2)))
        out = np.zeros_like(Xp)
        good = rho > 1e-12 * np.abs(rho).max(initial=1e-300)
        for i in range(self.n):
            if frozen[i]:
                continue
            num = np.real((pre[i] * QX[i][1] * suf[i + 1]).sum((1, 2)))
            out[:, i] = np.where(good, num / np.where(good, rho, 1), prev[:, i])
        return out

    # -- the source event: place all photon balls at once with density ρ (rejection sampling)
    def sample(self, N, rng, z0=DELTA):
        prof = {}
        for i in range(self.n):
            arm = self.arms[self.ph[i]]
            for k in range(self.K):
                b = self.spec[i][k][0]
                if (arm.name, b) in prof:
                    continue
                ch = arm.sources[b]
                x = ch.center(z0) + np.arange(-ch.Uh, ch.Uh + DU / 2, DU / 4)
                f, fx, fz = ch.eval(x, z0)
                j = np.clip(np.imag(np.conj(f) * fz), 0, None)
                cdf = np.cumsum(j); cdf /= cdf[-1]
                dens = j / (j.sum() * DU / 4)
                prof[(arm.name, b)] = (x, cdf, dens)

        def draw(m):
            ks = rng.integers(0, self.K, m)
            Xp = np.zeros((m, self.n))
            for i in range(self.n):
                for k in range(self.K):
                    sel = ks == k
                    x, cdf, _ = prof[(self.ph[i], self.spec[i][k][0])]
                    idx = np.searchsorted(cdf, rng.random(sel.sum()))
                    Xp[sel, i] = x[np.clip(idx, 0, x.size - 1)] + (rng.random(sel.sum()) - 0.5) * DU / 4
            q = np.zeros(m)
            for k in range(self.K):
                qk = np.ones(m)
                for i in range(self.n):
                    x, _, dens = prof[(self.ph[i], self.spec[i][k][0])]
                    qk *= np.interp(Xp[:, i], x, dens, left=0, right=0)
                q += qk / self.K
            return Xp, q

        zs = np.full(self.n, z0)
        Xp, q = draw(20000)
        r = self.density(Xp, zs) / np.maximum(q, 1e-300)
        M = 1.3 * r.max()
        out = []
        got = 0
        while got < N:
            Xp, q = draw(max(2 * (N - got), 2000))
            rho = self.density(Xp, zs)
            acc = rng.random(Xp.shape[0]) * M * q < rho
            out.append(Xp[acc]); got += acc.sum()
        return np.concatenate(out)[:N]

    # -- the dynamics: all photon balls move together (common clock z), RK4 between planes
    def run(self, X0, z0=DELTA, log=None):
        Xp = X0.copy()
        N = Xp.shape[0]
        frozen = np.zeros(self.n, bool)
        prev = np.zeros_like(Xp)
        t, last = z0, 0.0
        stops = sorted(set(self.planes) | set(self.zdet.tolist()))
        nsteps = 0
        for P in stops:
            if P <= t + 1e-9:
                continue
            while t < P - 1e-9:
                dz = min(DZMAX, max(DZMIN, GROW * (t - last)), P - t)
                zs = lambda tt: np.where(frozen, self.zdet, tt)
                k1 = self.slopes(Xp, zs(t), frozen, prev)
                k2 = self.slopes(Xp + 0.5 * dz * k1, zs(t + 0.5 * dz), frozen, k1)
                k3 = self.slopes(Xp + 0.5 * dz * k2, zs(t + 0.5 * dz), frozen, k2)
                k4 = self.slopes(Xp + dz * k3, zs(t + dz), frozen, k3)
                prev = (k1 + 2 * k2 + 2 * k3 + k4) / 6
                Xp = Xp + dz * prev
                t += dz; nsteps += 1
            t = P
            frozen = frozen | (np.abs(self.zdet - P) < 1e-9)
            if P in self.planes and not frozen.all():
                # first step behind the object along the outgoing current (implicit, 2 iterations)
                zs1 = np.where(frozen, self.zdet, P + DELTA)
                X1 = Xp.copy()
                for _ in range(2):
                    X1 = Xp + DELTA * self.slopes(X1, zs1, frozen, prev)
                prev = (X1 - Xp) / DELTA
                Xp = X1
                t, last = P + DELTA, P
        self.nsteps = nsteps
        return Xp

    def ports(self, Xf):
        return [self.arms[self.ph[i]].port(Xf[:, i]) for i in range(self.n)]

    @property
    def n_balls(self):
        return sum(a.n_balls for a in self.arms.values())


# ─────────────────────────────── standard objects ──────────────────────────────────────
def bs_rules(inA="a", inB="b", outC="c", outD="d", thA=ALPHA, thB=-ALPHA, t=T50, r=R50):
    """50:50 beam splitter: input A (arriving tilted thA) and B; output C leaves at +ALPHA, D at −ALPHA.
    Transmission keeps the side, reflection swaps it."""
    return {inA: [(outC, ALPHA, t), (outD, -ALPHA, r)], inB: [(outC, ALPHA, r), (outD, -ALPHA, t)]}


def pbs_rules(e1, inlab="in", plus="+", minus="-"):
    """Polarizing beam splitter projecting on the orthonormal basis (e1, e2): e1 -> +ALPHA, e2 -> −ALPHA."""
    e1 = np.asarray(e1, complex); e1 = e1 / np.linalg.norm(e1)
    e2 = np.array([-np.conj(e1[1]), np.conj(e1[0])])
    return {inlab: [(plus, ALPHA, np.outer(e1, np.conj(e1))), (minus, -ALPHA, np.outer(e2, np.conj(e2)))]}


def pbs_arm(name, e1, zP=60.0, zD=160.0, w0=5.0, zw=60.0):
    return Arm(name, zD).source("s", "in", 0.0, 0.0, w0, zw).obj(zP, pbs_rules(e1), "PBS")


def bs_arm(name, zB=60.0, zD=160.0, w0=5.0):
    x0 = zB * np.tan(ALPHA)
    return (Arm(name, zD).source("a", "a", -x0, ALPHA, w0, zB).source("b", "b", x0, -ALPHA, w0, zB)
            .obj(zB, bs_rules(), "BS"))


def basis(name):
    """Measurement basis vector e1 (outcome +1 = port '+')."""
    if name == "R":
        return RC
    return lin(np.deg2rad(float(name)))


# ─────────────────────────────── self test ─────────────────────────────────────────────
def selftest():
    global EXACT
    t0 = time.time()
    print("multi_core self-test", flush=True)
    arm = pbs_arm("A", lin(np.deg2rad(30)))
    sysA = System([arm], ["A"], [(1.0, [("s", H, None)])])
    print(f"  arm: {len(arm.channels)} channels, {arm.n_balls} wave balls, tables built in {sysA.build_time:.1f}s", flush=True)
    rng = np.random.default_rng(1)
    # 1. table vs direct ball sum (and vs the REGOLA hat-window reading of the same fan)
    worst = 0.0
    for ch in arm.channels:
        z = rng.uniform(ch.ze + DELTA, ch.z_end, 400)
        x = ch.center(z) + rng.normal(0, 1, 400) * (ch.sig0 + (z - ch.ze) * ch.sig_t) * 1.5
        d = ch.exact(x, z); c = ch.eval(x, z)
        scale = np.abs(d[0]).max()
        err = max(np.abs(c[0] - d[0]).max() / scale, np.abs(c[1] - d[1]).max() / (K * scale),
                  np.abs(c[2] - d[2]).max() / (K * scale))
        worst = max(worst, err)
    fan = Fan(arm.channels[0].s, arm.channels[0].a, 0.0, dth=0.1)
    xs = rng.normal(0, 5, 50)
    U, Ux, Uz = field_lattice(xs, 30.0, fan)
    d = arm.channels[0].exact(xs, 30.0)
    hat = np.abs(U[:, 0] - d[0]).max() / np.abs(d[0]).max()
    print(f"  table vs direct ball sum: max rel. error {worst:.2e} (field and gradients); "
          f"hat-window reading vs direct sum {hat:.1e}", flush=True)
    # 2. reduction: a product state of 2 photons moves exactly like 2 single photons
    arm2 = pbs_arm("B", lin(np.deg2rad(-20)))
    sys2 = System([arm, arm2], ["A", "B"], [(1.0, [("s", lin(0.3), None), ("s", lin(1.1), None)])])
    s1 = System([arm], ["A"], [(1.0, [("s", lin(0.3), None)])])
    X = np.column_stack([rng.normal(0, 4, 300), rng.normal(0, 4, 300)])
    red = 0.0
    for z in (20.0, 70.0, 120.0):
        v2 = sys2.slopes(X, np.array([z, z]), np.zeros(2, bool), np.zeros_like(X))
        v1 = s1.slopes(X[:, :1], np.array([z]), np.zeros(1, bool), np.zeros((300, 1)))
        red = max(red, np.abs(v2[:, 0] - v1[:, 0]).max())
    print(f"  product state: |v(2-photon rule) - v(1-photon rule)| max {red:.1e}", flush=True)
    # 3. table vs direct: trajectories
    X0 = sysA.sample(40, np.random.default_rng(3))
    xf = sysA.run(X0)
    EXACT = True
    xe = sysA.run(X0)
    EXACT = False
    print(f"  40 trajectories, table vs direct ball sum: max |Δx| at the detector {np.abs(xf - xe).max():.1e} λ "
          f"(beam width ~8 λ)  [{time.time() - t0:.0f}s]", flush=True)
    ok = worst < 2e-3 and red < 1e-10 and np.abs(xf - xe).max() < 0.05
    print("  SELFTEST", "PASS" if ok else "FAIL", flush=True)
    return ok


if __name__ == "__main__":
    selftest()

"""THE BASE RULE — one fixed rule, two roles (wave balls + ONE photon ball), in 2-D and in 3-D.
Spec in plain Italian: docs/REGOLA_BASE.md.  Changelog of the rule versions: RULE_LOG below.

WAVE BALLS (every object does the same thing: source, aperture, polarizer, beam splitter, mirror, plate)
  * an object re-emits the wave that reaches it as a regular FAN of wave balls:
      origins on a regular grid over the illuminated part of the object (cell centroids, amplitude
      a = B(origin) · cell area, with B = T·ψ_in(+side) + R·ψ_in(−side) the boundary value of the
      outgoing wave, T, R the object's transmission / reflection Jones matrices);
      directions on a regular grid covering the whole forward half-space (Θ = 90°).
  * HOW MANY balls is not chosen: it is the field's number of degrees of freedom (Shannon number)
      level 0 ("propagating" balls):  ds0 = λ/(2η)   (Nyquist of the propagating band |f| ≤ 1/λ,
                                       oversampling η)  →  N0 = 2ηW/λ origins per transverse axis,
                                       d(sinθ) = λ/(ηW)  →  2ηW/λ directions per axis,
      level ℓ ≥ 1 ("evanescent" balls): the evanescent band that is still above ε at distance R is
                                       f(R) = ln(1/ε)/(2πR); Nyquist with the same η gives
                                       ds_ℓ = ds0/2^ℓ, heard only for R < R0·2^(2−ℓ), R0 = ln(1/ε)λ/(2π).
      W = extent of the object where |B| ≥ ε·max|B|.   ℓ_max is fixed by the closest distance at which a
      photon ever reads a fan: δ (the crossing distance) → ℓ_max = ceil(log2(R0/δ)) + 1.
  * a wave ball flies straight at c; it dies when it meets another object (it is re-emitted there).

PHOTON BALL (one per photon)
  * born on the source with probability |B_src|² (where the light is); local coin at the source (below);
  * at every step it hears the wave balls that pass near it — the balls of every fan whose ray from
    its origin reaches the photon unblocked (2-D: exact ray test; 3-D: planes ⊥ z, so: the last
    object crossed) — and adds their exact wavelets (2-D: (ik/2)(R·n)H1(kR)/R; 3-D: vector
    Rayleigh–Sommerfeld, E_t = −2 Σ B_t ∂_n g, E_n = 2 Σ B_t·∇_t g, g = e^{ikR}/(4πR)), level
    weights w_ℓ(R) (partition of unity in log R);
  * j = Σ_pol Im(ψ*∇ψ) (3-D: Σ over E_x, E_y, E_z); it moves at c along j (Dormand–Prince 5(4),
    fixed tolerance); no wave heard → straight on;
  * at an object: opaque part → dies (absorbed);  transparent part → flies straight to δ beyond;
    ABSORBING object (aperture walls, polarizer, source): survives with p = min(1, Φ_out/Φ_in),
    Φ = n·j read at its own position (in: just before, out: at δ);  NON-absorbing (beam splitter,
    phase plate): no coin, it simply continues with the current on the other side; screen → click.

Everything that depends on the experiment (W, number of balls, ℓ_max cut by ε) is an OUTPUT.
The sums over origins run in Rust (lab/base_kernel, same arithmetic; numpy reference below,
tested by base_test.py).
"""
import ctypes
import os
import subprocess
import sys
import time

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

LAM = 1.0
K = 2 * np.pi / LAM
NTHREADS = int(os.environ.get("BASE_THREADS", "4"))

# ───────────────────────────── THE RULE (frozen; change = new version) ─────────────────────────────
RULE = dict(
    version="v2",
    eta=2.0,       # oversampling of the Shannon/Nyquist count (one number for every object, 2-D and 3-D)
    eps=1e-3,      # amplitude below which a field component is "absent" (defines W and the evanescent range)
    delta=0.02,    # distance (λ) behind an object at which the photon resumes / reads Φ_out
    atol=1e-4,     # trajectory integrator: position error per step (λ)
    hmax=2.0,      # largest step (λ)
    hmin=1e-6,
    edge_zone=0.5,  # v2: width (λ) of the zone next to an edge where the crossing keeps the photon's rank in the flux
)
RULE_LOG = """
v1  (2026-10-03) eta=2, eps=1e-3, delta=0.02, DP45 atol=1e-4, hmax=2.  First frozen version.
v2  (2026-10-03) + edge zone (λ/2): a photon crossing an absorbing object closer than λ/2 to an edge is not
    flown straight + coin; it keeps its rank in the flux across the edge (monotone map of the incident flux
    onto the outgoing flux j_z at δ, both read on the local cut through the edge) and survives with the
    zone's flux ratio.  Reason: the Arago spot (v1 FAIL) is fed by the ~2 % of the flux that the sharp-edge
    (Rayleigh–Sommerfeld/Kirchhoff) wave sends into the geometric shadow within δ, at ANY δ; a coin can only
    remove photons and never put them there.
"""


def derived(rule=RULE):
    eta, eps, delta = rule["eta"], rule["eps"], rule["delta"]
    r0 = np.log(1 / eps) * LAM / (2 * np.pi)
    ds0 = LAM / (2 * eta)
    lmax = int(np.ceil(np.log2(r0 / delta))) + 1
    return dict(r0=r0, ds0=ds0, lmax=lmax, sds=ds0 / 4, dtouch=delta / 4)


# ───────────────────────────── Rust bridge ──────────────────────────────────────────────────────
_KDIR = os.path.join(HERE, "base_kernel")
_DLL = os.path.join(_KDIR, "target", "release", "base_kernel.dll" if os.name == "nt" else "libbase_kernel.so")
_lib = None


def lib():
    global _lib
    if _lib is None:
        src = os.path.join(_KDIR, "src", "lib.rs")
        if not os.path.exists(_DLL) or os.path.getmtime(_DLL) < os.path.getmtime(src):
            subprocess.run(["cargo", "build", "--release"], cwd=_KDIR, check=True)
        L = ctypes.CDLL(_DLL)
        P = ctypes.POINTER(ctypes.c_double)
        L.scene_new.restype = ctypes.c_void_p
        L.scene_new.argtypes = [P, ctypes.c_size_t]
        L.scene_free.argtypes = [ctypes.c_void_p]
        L.fan_origins.restype = ctypes.c_size_t
        L.fan_origins.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t] + [ctypes.c_double] * 4 + [P, ctypes.c_size_t]
        L.field2.argtypes = [ctypes.c_void_p, ctypes.c_size_t, P, P, P, ctypes.c_int64, P, ctypes.c_size_t]
        L.field3.argtypes = [ctypes.c_void_p, ctypes.c_size_t, P, P, P, ctypes.c_size_t]
        _lib = L
    return _lib


def _p(a):
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


# ───────────────────────────── Gaussian beams (the sources) ─────────────────────────────────────
def gauss_beam(r2, w0, zrel, dim):
    """Paraxial Gaussian beam on a plane at distance zrel from its waist (zrel < 0: before the waist).
    r2 = transverse distance², dim = 2 (cylindrical beam) or 3.  w0 = inf → plane wave."""
    if not np.isfinite(w0):
        return np.ones_like(r2, dtype=complex)
    zr = np.pi * w0 ** 2 / LAM
    q0 = -1j * zr
    q = zrel - 1j * zr
    amp = (q0 / q) ** (0.5 if dim == 2 else 1.0)
    return amp * np.exp(1j * K * r2 / (2 * q))


# ───────────────────────────── 2-D objects and scene ────────────────────────────────────────────
class Obj2:
    """A segment object in the x–z plane.  Points are (x, z).
    c: centre, t: unit tangent, half: half-length, normal n = (−t_z, t_x) is its '+' side.
    kind: 'source' | 'aperture' (absorbing walls, transmission T on `ints`) | 'polarizer' (absorbing)
          | 'bs' | 'mirror' | 'phase' (non-absorbing) | 'stop' (opaque) | 'screen' (click) | 'plane' (passive record)
    ints: transmitting intervals in the tangent coordinate (default: the whole segment)."""

    def __init__(self, name, kind, c, t, half, ints=None, T=None, R=None, absorbing=None, src=None):
        self.name, self.kind = name, kind
        self.c = np.asarray(c, float)
        t = np.asarray(t, float)
        self.t = t / np.hypot(*t)
        self.n = np.array([-self.t[1], self.t[0]])
        self.half = float(half)
        self.ints = [(-half, half)] if ints is None else [tuple(map(float, iv)) for iv in ints]
        I2 = np.eye(2, dtype=complex)
        self.T = None if T is None else np.asarray(T, complex) * I2 if np.ndim(T) == 0 else np.asarray(T, complex)
        self.R = None if R is None else np.asarray(R, complex) * I2 if np.ndim(R) == 0 else np.asarray(R, complex)
        if kind in ("aperture", "polarizer", "phase", "source") and self.T is None:
            self.T = I2
        self.absorbing = kind in ("aperture", "polarizer", "source") if absorbing is None else absorbing
        self.src = src                # dict(w0, zw, jones): Gaussian beam with waist at distance zw along +n
        self.blocks = kind not in ("source", "plane")

    def p1(self):
        return self.c - self.half * self.t

    def p2(self):
        return self.c + self.half * self.t

    def acts(self, s):
        """the part of the object where it transmits or reflects (its intervals)"""
        s = np.atleast_1d(s)
        ok = np.zeros(s.shape, bool)
        for a, b in self.ints:
            ok |= (s >= a) & (s <= b)
        return ok

    def transmits(self, s):
        s = np.atleast_1d(s)
        ok = np.zeros(s.shape, bool)
        if self.T is None:
            return ok
        for a, b in self.ints:
            ok |= (s >= a) & (s <= b)
        return ok

    def src_B(self, s):
        """Boundary value of a source: Gaussian beam × polarization (rows: s, cols: Jones)."""
        u = gauss_beam(s ** 2, self.src["w0"], -self.src.get("zw", 0.0), 2)
        return u[:, None] * np.asarray(self.src.get("jones", (1, 0)), complex)[None, :]


class Scene2D:
    def __init__(self, rule=RULE):
        self.rule = dict(rule)
        self.d = derived(rule)
        self.objs = []
        self.fans = []                # python records (dict) in the same order as in Rust
        self.tog_obj, self.tog_init = -1, True
        self.ptr = None
        self.build_time = 0.0

    # --- packing -----------------------------------------------------------------------------
    def _pack(self):
        d = self.d
        buf = [2, K, d["r0"], d["ds0"], d["lmax"], len(self.fans),
               sum(o.blocks for o in self.objs), self.tog_obj, 1.0 if self.tog_init else 0.0]
        for i, o in enumerate(self.objs):
            if o.blocks:
                buf += [i, *o.p1(), *o.p2()]
        for f in self.fans:
            buf += [f["obj"], *f["c"], *f["t"], *f["n"], f["glo"], f["ghi"], len(f["ints"])]
            for iv in f["ints"]:
                buf += list(iv)
            buf += [f["sp0"], f["sds"], len(f["samp"])]
            buf += list(np.column_stack([f["samp"][:, 0].real, f["samp"][:, 0].imag,
                                         f["samp"][:, 1].real, f["samp"][:, 1].imag]).ravel())
            buf += [f["bthr"]]
        arr = np.ascontiguousarray(buf, float)
        if self.ptr:
            lib().scene_free(self.ptr)
        self.ptr = lib().scene_new(_p(arr), arr.size)
        self._buf = arr

    def __del__(self):
        try:
            if self.ptr:
                lib().scene_free(self.ptr)
        except Exception:
            pass

    def set_toggle(self, obj_index, init_present):
        self.tog_obj, self.tog_init = obj_index, bool(init_present)
        self._pack()

    # --- field --------------------------------------------------------------------------------
    def field(self, P, t=None, tog=None, excl=-999):
        P = np.ascontiguousarray(np.atleast_2d(P), float)
        n = P.shape[0]
        t = np.zeros(n) if t is None else np.ascontiguousarray(np.broadcast_to(t, (n,)), float)
        tog = np.full(n, np.inf) if tog is None else np.ascontiguousarray(np.broadcast_to(tog, (n,)), float)
        out = np.zeros((n, 12))
        if n:
            lib().field2(self.ptr, n, _p(P), _p(t), _p(tog), int(excl), _p(out), NTHREADS)
        z = out[:, 0::2] + 1j * out[:, 1::2]
        return z[:, 0:2], z[:, 2:4], z[:, 4:6]          # ψ, ∂xψ, ∂zψ  (n × 2 Jones)

    def current(self, P, t=None, tog=None, excl=-999):
        U, Ux, Uz = self.field(P, t, tog, excl)
        jx = np.imag(np.conj(U) * Ux).sum(-1)
        jz = np.imag(np.conj(U) * Uz).sum(-1)
        return np.column_stack([jx, jz]), (np.abs(U) ** 2).sum(-1), U

    # --- objects ------------------------------------------------------------------------------
    def add(self, o):
        """Add an object; build its fan(s) from the wave that reaches it (in order of the light)."""
        t0 = time.time()
        i = len(self.objs)
        self.objs.append(o)
        if self.ptr is None or o.blocks:
            self._pack()                       # the new object already blocks rays (for its incident field)
        if o.kind in ("stop", "screen", "plane"):
            return i
        d = self.d
        ds0, sds = d["ds0"], d["sds"]
        hull = (max(-o.half, min(a for a, b in o.ints)), min(o.half, max(b for a, b in o.ints)))
        if o.kind == "source":
            Bfun = {+1: o.src_B}
        else:
            def psi(side):
                def f(s):
                    P = o.c[None, :] + s[:, None] * o.t[None, :] + side * 1e-9 * o.n[None, :]
                    return self.field(P, excl=i)[0]
                return f
            pm, pp = psi(-1), psi(+1)
            Bfun = {}
            for side in (+1, -1):            # fan emitted into the `side` half-space
                inc_from_other = pm if side == +1 else pp
                inc_same = pp if side == +1 else pm
                def B(s, a=inc_from_other, b=inc_same):
                    out = 0
                    if o.T is not None:
                        out = out + a(s) @ o.T.T
                    if o.R is not None:
                        out = out + b(s) @ o.R.T
                    return out if np.ndim(out) else np.zeros((s.size, 2), complex)
                Bfun[side] = B
        for side, B in Bfun.items():
            sc = np.arange(hull[0] + ds0 / 2, hull[1], ds0 / 4)
            if sc.size == 0:
                continue
            Bc = B(sc)
            mag = np.sqrt((np.abs(Bc) ** 2).sum(-1)) * o.acts(sc) if o.kind != "source" else \
                np.sqrt((np.abs(Bc) ** 2).sum(-1)) * o.transmits(sc)
            if mag.max(initial=0) <= 1e-300:
                continue
            on = np.nonzero(mag >= RULE_EPS(self) * mag.max())[0]
            glo = max(hull[0], sc[on[0]] - ds0)
            ghi = min(hull[1], sc[on[-1]] + ds0)
            sp0 = glo - 3 * sds
            ns = int(np.ceil((ghi - glo) / sds)) + 7
            ss = sp0 + np.arange(ns) * sds
            samp = B(ss)
            inside = o.acts(ss)
            bthr = RULE_EPS(self) * np.sqrt((np.abs(samp[inside]) ** 2).sum(-1)).max(initial=0)
            n = o.n * side
            self.fans.append(dict(obj=i, side=side, c=o.c, t=o.t, n=n, glo=glo, ghi=ghi,
                                  ints=[(max(a, glo), min(b, ghi)) for a, b in o.ints if min(b, ghi) > max(a, glo)],
                                  sp0=sp0, sds=sds, samp=np.ascontiguousarray(samp), bthr=bthr))
        self._pack()
        self.build_time += time.time() - t0
        return i

    # --- bookkeeping: how many balls the rule made -----------------------------------------------
    def origins(self, fi, level=0):
        cap = 2_000_000
        buf = np.zeros(cap * 6)
        n = lib().fan_origins(self.ptr, fi, level, -1e300, 1e300, 0, 0, _p(buf), cap)
        b = buf[:6 * min(n, cap)].reshape(-1, 6)
        return b[:, :2], np.column_stack([b[:, 2] + 1j * b[:, 3], b[:, 4] + 1j * b[:, 5]])

    def ball_count(self):
        """Per fan: (object, N origins level 0, W, directions/origin, balls level 0, balls all levels)."""
        eta, d = self.rule["eta"], self.d
        rows = []
        for fi, f in enumerate(self.fans):
            P, _ = self.origins(fi, 0)
            if len(P) == 0:
                continue
            s = (P - f["c"]) @ f["t"]
            W = s.max() - s.min() + d["ds0"]
            ndir = int(np.ceil(2 * eta * W / LAM))
            n0 = len(P)
            nall = sum(n0 * 2 ** l for l in range(d["lmax"] + 1))
            rows.append(dict(obj=self.objs[f["obj"]].name, side=f["side"], n0=n0, W=W, ndir=ndir,
                             balls0=n0 * ndir, balls=nall * ndir))
        return rows


def RULE_EPS(scene):
    return scene.rule["eps"]


# ───────────────────────────── 3-D objects and scene (planes ⊥ z) ──────────────────────────────
class Mask3:
    """'full' | ('circ_in', cx, cy, r) | ('circ_out', cx, cy, r) | ('rects', [(x0,x1,y0,y1), ...])"""

    def __init__(self, spec="full"):
        self.spec = spec

    def code(self):
        s = self.spec
        if s == "full":
            return [0, 0]
        if s[0] == "circ_in":
            return [1, 3, *s[1:4]]
        if s[0] == "circ_out":
            return [2, 3, *s[1:4]]
        rects = s[1]
        return [3, 1 + 4 * len(rects), len(rects)] + [v for r in rects for v in r]

    def inside(self, x, y):
        s = self.spec
        if s == "full":
            return np.ones(np.shape(x), bool)
        if s[0] == "circ_in":
            return (x - s[1]) ** 2 + (y - s[2]) ** 2 < s[3] ** 2
        if s[0] == "circ_out":
            return (x - s[1]) ** 2 + (y - s[2]) ** 2 >= s[3] ** 2
        ok = np.zeros(np.shape(x), bool)
        for x0, x1, y0, y1 in s[1]:
            ok |= (x >= x0) & (x < x1) & (y >= y0) & (y < y1)
        return ok

    def hull(self):
        s = self.spec
        if s == "full" or s[0] == "circ_out":
            return None
        if s[0] == "circ_in":
            return (s[1] - s[3], s[1] + s[3], s[2] - s[3], s[2] + s[3])
        r = np.array(s[1])
        return (r[:, 0].min(), r[:, 1].max(), r[:, 2].min(), r[:, 3].max())


class Obj3:
    def __init__(self, name, kind, z, mask="full", T=None, absorbing=None, src=None):
        self.name, self.kind, self.z = name, kind, float(z)
        self.mask = mask if isinstance(mask, Mask3) else Mask3(mask)
        self.T = (np.eye(2, dtype=complex) if T is None else np.asarray(T, complex)) \
            if kind in ("aperture", "polarizer", "phase", "source") else None
        self.absorbing = kind in ("aperture", "polarizer", "source") if absorbing is None else absorbing
        self.src = src

    def src_B(self, x, y):
        """Gaussian beam (w0, waist at distance zw downstream) × Jones; w0y: different waist along y
        (separable paraxial beam = product of two cylindrical beams)."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        zw = -self.src.get("zw", 0.0)
        if "w0y" in self.src:
            u = gauss_beam(x ** 2, self.src["w0"], zw, 2) * gauss_beam(y ** 2, self.src["w0y"], zw, 2)
        else:
            u = gauss_beam(x ** 2 + y ** 2, self.src["w0"], zw, 3)
        return u[..., None] * np.asarray(self.src.get("jones", (1, 0)), complex)


class Scene3D:
    def __init__(self, rule=RULE):
        self.rule = dict(rule)
        self.d = derived(rule)
        self.objs = []
        self.fans = []
        self.fan_of = {}            # object index -> fan index
        self.ptr = None
        self.build_time = 0.0

    def _pack(self):
        d = self.d
        buf = [3, K, d["r0"], d["ds0"], d["lmax"], len(self.fans)]
        for f in self.fans:
            buf += [f["obj"], f["z0"], *f["g"], *f["mask"].code(), f["sx0"], f["sy0"], f["sds"], f["nx"], f["ny"]]
            s = f["samp"].reshape(-1, 2)
            buf += list(np.column_stack([s[:, 0].real, s[:, 0].imag, s[:, 1].real, s[:, 1].imag]).ravel())
            buf += [f["bthr"]]
        arr = np.ascontiguousarray(buf, float)
        if self.ptr:
            lib().scene_free(self.ptr)
        self.ptr = lib().scene_new(_p(arr), arr.size)

    def __del__(self):
        try:
            if self.ptr:
                lib().scene_free(self.ptr)
        except Exception:
            pass

    def field(self, P, hear):
        P = np.ascontiguousarray(np.atleast_2d(P), float)
        n = P.shape[0]
        h = np.ascontiguousarray(np.broadcast_to(hear, (n,)), float)
        out = np.zeros((n, 24))
        if n:
            lib().field3(self.ptr, n, _p(P), _p(h), _p(out), NTHREADS)
        z = out[:, 0::2] + 1j * out[:, 1::2]
        return z[:, 0:3], z[:, 3:12].reshape(n, 3, 3)       # E_a, ∂_b E_a  [n, a, b]

    def current(self, P, hear):
        E, dE = self.field(P, hear)
        j = np.imag(np.conj(E)[:, :, None] * dE).sum(1)
        return j, (np.abs(E) ** 2).sum(-1), E

    def last_fan(self):
        return len(self.fans) - 1

    def add(self, o):
        t0 = time.time()
        i = len(self.objs)
        self.objs.append(o)
        if o.kind in ("stop", "screen", "plane"):
            return i
        d, eps = self.d, self.rule["eps"]
        ds0, sds = d["ds0"], d["sds"]
        if o.kind == "source":
            def B(x, y):
                return o.src_B(x, y)
            def ext(w0):
                if not np.isfinite(w0):
                    return 50.0
                zr = np.pi * w0 ** 2
                return 4 * w0 * np.sqrt(1 + (o.src.get("zw", 0.0) / zr) ** 2)
            bx, by = ext(o.src["w0"]), ext(o.src.get("w0y", o.src["w0"]))
            box = [-bx, bx, -by, by]
        else:
            prev = self.last_fan()
            pf = self.fans[prev]

            def B(x, y):
                P = np.column_stack([np.ravel(x), np.ravel(y), np.full(np.size(x), o.z)])
                E, _ = self.field(P, prev)
                return (E[:, :2] @ o.T.T).reshape(np.shape(x) + (2,))
            m = 0.5 * (o.z - pf["z0"]) + 2.0
            box = [pf["g"][0] - m, pf["g"][1] + m, pf["g"][2] - m, pf["g"][3] + m]
        hull = o.mask.hull()
        if hull is not None:
            box = [max(box[0], hull[0]), min(box[1], hull[1]), max(box[2], hull[2]), min(box[3], hull[3])]
        # find the region where |B| >= eps·max (grow the box if the light reaches its border)
        for _ in range(8):
            xs = np.arange(box[0] + ds0 / 2, box[1], ds0)
            ys = np.arange(box[2] + ds0 / 2, box[3], ds0)
            X, Y = np.meshgrid(xs, ys)
            mag = np.sqrt((np.abs(B(X, Y)) ** 2).sum(-1)) * o.mask.inside(X, Y)
            mx = mag.max(initial=0)
            if mx <= 1e-300:
                return i
            on = mag >= eps * mx
            edge = on[0, :].any() or on[-1, :].any() or on[:, 0].any() or on[:, -1].any()
            if not edge or hull is not None and all(abs(a - b) < 1e-9 for a, b in zip(box, hull)):
                break
            cx, cy = 0.5 * (box[0] + box[1]), 0.5 * (box[2] + box[3])
            hx, hy = 0.75 * (box[1] - box[0]), 0.75 * (box[3] - box[2])
            box = [cx - hx, cx + hx, cy - hy, cy + hy]
            if hull is not None:
                box = [max(box[0], hull[0]), min(box[1], hull[1]), max(box[2], hull[2]), min(box[3], hull[3])]
        iy, ix = np.nonzero(on)
        g = [max(box[0], xs[ix.min()] - ds0), min(box[1], xs[ix.max()] + ds0),
             max(box[2], ys[iy.min()] - ds0), min(box[3], ys[iy.max()] + ds0)]
        sx0, sy0 = g[0] - 3 * sds, g[2] - 3 * sds
        nx = int(np.ceil((g[1] - g[0]) / sds)) + 7
        ny = int(np.ceil((g[3] - g[2]) / sds)) + 7
        SX, SY = np.meshgrid(sx0 + np.arange(nx) * sds, sy0 + np.arange(ny) * sds)
        samp = np.zeros((ny, nx, 2), complex)
        for r0 in range(0, ny, 64):                     # in slabs (memory)
            samp[r0:r0 + 64] = B(SX[r0:r0 + 64], SY[r0:r0 + 64])
        ins = o.mask.inside(SX, SY)
        bthr = eps * np.sqrt((np.abs(samp[ins]) ** 2).sum(-1)).max(initial=0)
        self.fans.append(dict(obj=i, z0=o.z, g=g, mask=o.mask, sx0=sx0, sy0=sy0, sds=sds, nx=nx, ny=ny,
                              samp=np.ascontiguousarray(samp), bthr=bthr))
        self.fan_of[i] = len(self.fans) - 1
        self._pack()
        self.build_time += time.time() - t0
        return i

    def origins(self, fi, level=0):
        cap = 4_000_000
        buf = np.zeros(cap * 6)
        n = lib().fan_origins(self.ptr, fi, level, -1e300, 1e300, -1e300, 1e300, _p(buf), cap)
        b = buf[:6 * min(n, cap)].reshape(-1, 6)
        return b[:, :2], np.column_stack([b[:, 2] + 1j * b[:, 3], b[:, 4] + 1j * b[:, 5]])

    def ball_count(self):
        eta, d = self.rule["eta"], self.d
        rows = []
        for fi, f in enumerate(self.fans):
            P, _ = self.origins(fi, 0)
            if len(P) == 0:
                continue
            Wx = P[:, 0].max() - P[:, 0].min() + d["ds0"]
            Wy = P[:, 1].max() - P[:, 1].min() + d["ds0"]
            ndir = int(np.ceil(np.pi / 4 * (2 * eta * Wx / LAM) * (2 * eta * Wy / LAM)))
            n0 = len(P)
            nall = sum(n0 * 4 ** l for l in range(d["lmax"] + 1))
            rows.append(dict(obj=self.objs[f["obj"]].name, n0=n0, W=(Wx, Wy), ndir=ndir,
                             balls0=n0 * ndir, balls=nall * ndir))
        return rows


def balls_summary(scene):
    rows = scene.ball_count()
    b0 = sum(r["balls0"] for r in rows)
    ba = sum(r["balls"] for r in rows)
    txt = "; ".join(f"{r['obj']}{'' if r.get('side', 1) == 1 else '(-)'}: {r['n0']} origins x {r['ndir']} dir"
                    for r in rows)
    return dict(per_fan=txt, balls_level0=int(b0), balls_all_levels=int(ba), lmax=scene.d["lmax"])


# ───────────────────────────── the photon ball ────────────────────────────────────────────────
# Dormand–Prince 5(4)
_C = np.array([0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1])
_A = [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9],
      [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
      [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
      [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]]
_B5 = np.array([35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0])
_B4 = np.array([5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40])

ALIVE, CLICK, ABSORBED, COIN, LOST, STUCK = 0, 1, 2, 3, 4, 5
STATUS = {CLICK: "click", ABSORBED: "absorbed (opaque part)", COIN: "absorbed (local coin)",
          LOST: "left the lab", STUCK: "step limit"}


class Photons:
    """Runs photon balls through a scene. dim = 2 (Scene2D) or 3 (Scene3D)."""

    def __init__(self, scene, box, rng, tog=None, max_steps=20000):
        self.sc = scene
        self.dim = 2 if isinstance(scene, Scene2D) else 3
        self.box = np.asarray(box, float)          # [lo0, hi0, lo1, hi1(, lo2, hi2)]
        self.rng = rng
        self.r = scene.rule
        self.dd = scene.d
        self.max_steps = max_steps
        self.tog = tog
        self._planes = [i for i, o in enumerate(scene.objs) if o.kind == "plane"]

    # velocity field: unit vector along j (own direction where no wave is heard)
    def vel(self, X, t, dirs, hear, tog):
        if self.dim == 2:
            j, rho, _ = self.sc.current(X, t, tog)
        else:
            j, rho, _ = self.sc.current(X, hear)
        nj = np.sqrt((j ** 2).sum(-1))
        ok = (rho > 1e-300) & (nj > 1e-12 * K * np.maximum(rho, 1e-300))
        v = dirs.copy()
        v[ok] = j[ok] / nj[ok, None]
        return v

    # ---- births -----------------------------------------------------------------------------
    def birth(self, src_index, n):
        """n photon balls born on the source object with probability |B|² (where the light is);
        returns positions X (after the source coin, at δ in front of the source), directions, mask of survivors,
        plus the fraction of |B|² that falls on opaque parts of the source mask (absorbed before birth)."""
        o = self.sc.objs[src_index]
        rng = self.rng
        sds = self.dd["sds"]
        if self.dim == 2:
            f = [f for f in self.sc.fans if f["obj"] == src_index][0]
            lo, hi = -o.half, o.half
            if np.isfinite(o.src["w0"]):
                zr = np.pi * o.src["w0"] ** 2
                w = o.src["w0"] * np.sqrt(1 + (o.src.get("zw", 0) / zr) ** 2)
                lo, hi = max(lo, -5 * w), min(hi, 5 * w)
            s = np.arange(lo + sds / 2, hi, sds)
            p = (np.abs(o.src_B(s)) ** 2).sum(-1)
            pin = p * o.transmits(s)
            wall = 1 - pin.sum() / p.sum()
            # sample |B|² on every cell that touches the opening, then keep the points inside it (exact near edges)
            touch = o.transmits(s) | o.transmits(s - sds / 2) | o.transmits(s + sds / 2)
            pt = p * touch
            cdf = np.cumsum(pt) / pt.sum()
            k = np.searchsorted(cdf, rng.random(n))
            ss = s[k] + (rng.random(n) - 0.5) * sds
            keep = o.transmits(ss)
            ss = ss[keep]
            X0 = o.c[None, :] + ss[:, None] * o.t[None, :]
            nrm = np.repeat(o.n[None, :], ss.size, 0)
            Bm = (np.abs(o.src_B(ss)) ** 2).sum(-1)
            phi_in = K * Bm
        else:
            fi = self.sc.fan_of[src_index]
            f = self.sc.fans[fi]
            g = f["g"]
            xs = np.arange(g[0] + sds / 2, g[1], sds)
            ys = np.arange(g[2] + sds / 2, g[3], sds)
            X, Y = np.meshgrid(xs, ys)
            p = (np.abs(o.src_B(X, Y)) ** 2).sum(-1)
            touch = np.zeros(X.shape, bool)
            for ddx in (-0.5, 0, 0.5):
                for ddy in (-0.5, 0, 0.5):
                    touch |= o.mask.inside(X + ddx * sds, Y + ddy * sds)
            pt = p * touch
            cdf = np.cumsum(pt.ravel()) / pt.sum()
            k = np.searchsorted(cdf, rng.random(n))
            x = X.ravel()[k] + (rng.random(n) - 0.5) * sds
            y = Y.ravel()[k] + (rng.random(n) - 0.5) * sds
            keep = o.mask.inside(x, y)
            x, y = x[keep], y[keep]
            X0 = np.column_stack([x, y, np.full(x.size, o.z)])
            nrm = np.repeat(np.array([[0, 0, 1.0]]), x.size, 0)
            Bm = (np.abs(o.src_B(x, y)) ** 2).sum(-1)
            phi_in = K * Bm
            # fraction of the beam on the opaque part (over the whole beam, not only the fan's box)
            if np.isfinite(o.src["w0"]):
                def ext(w0):
                    zr = np.pi * w0 ** 2
                    return 4 * w0 * np.sqrt(1 + (o.src.get("zw", 0.0) / zr) ** 2)
                bx, by = ext(o.src["w0"]), ext(o.src.get("w0y", o.src["w0"]))
                hx, hy = min(sds, bx / 200), min(sds, by / 200)
                XX, YY = np.meshgrid(np.arange(-bx + hx / 2, bx, hx), np.arange(-by + hy / 2, by, hy))
                pp = (np.abs(o.src_B(XX, YY)) ** 2).sum(-1)
                wall = 1 - (pp * o.mask.inside(XX, YY)).sum() / pp.sum()
            else:
                wall = float("nan")             # plane wave: only the aperture is lit (not meaningful)
        # the source coin: fly δ along the normal, read Φ_out there
        delta = self.r["delta"]
        X1 = X0 + delta * nrm
        tg = None if self.tog is None else np.asarray(self.tog, float)[: len(keep)][keep]
        if self.dim == 2:
            j, _, _ = self.sc.current(X1, np.full(len(X1), delta), tg)
        else:
            j, _, _ = self.sc.current(X1, self.sc.fan_of[src_index])
        phi_out = (j * nrm).sum(-1)
        pr = np.clip(phi_out / np.maximum(phi_in, 1e-300), 0, 1)
        u, E, dv, W, z, cvv = self._edge_z(o, X0)
        if z.any():
            if self.dim == 2:
                tgz = None if tg is None else tg[z]

                def fin(P, M):
                    return K * (np.abs(o.src_B((P - o.c) @ o.t)) ** 2).sum(-1)

                def fout(P, M):
                    return self.sc.current(P, np.full(len(P), delta), None if tgz is None else np.repeat(tgz, len(P) // M))[0]
            else:
                def fin(P, M):
                    return K * (np.abs(o.src_B(P[:, 0], P[:, 1])) ** 2).sum(-1)

                def fout(P, M):
                    return self.sc.current(P, self.sc.fan_of[src_index])[0]
            X1[z], pr[z] = self._transport(u[z], E[z], dv[z], W[z], nrm[z], fin, fout, cvv[z])
        surv = rng.random(len(X1)) < pr
        self._keep = keep
        return X1, nrm, surv, wall, (n - keep.sum())


    # ---- v2: edges --------------------------------------------------------------------------
    def _edge(self, o, xc):
        """Nearest edge of the transmitting part of object o for crossing points xc.
        Returns u (distance from the edge into the opening; inf = not in an edge zone), E (edge point),
        dv (unit vector from the edge into the opening), W (zone half-width)."""
        n = len(xc)
        Ue = self.r.get("edge_zone", 0.0)
        u = np.full(n, np.inf)
        E = np.zeros_like(xc)
        dv = np.zeros_like(xc)
        W = np.zeros(n)
        self._curv = np.zeros(n)                 # 1/R of a circular edge (+: opening outside), for the area weight
        if Ue <= 0 or n == 0:
            return u, E, dv, W
        if self.dim == 2:
            s = (xc - o.c) @ o.t
            for a, b in o.ints:
                inside = (s >= a) & (s <= b)
                for end, sg in ((a, 1.0), (b, -1.0)):
                    if abs(abs(end) - o.half) < 1e-9:
                        continue                         # end of the object, not an edge
                    d = sg * (s - end)
                    m = inside & (d < u)
                    u[m] = d[m]
                    E[m] = o.c + end * o.t
                    dv[m] = sg * o.t
                    W[m] = min(Ue, 0.5 * (b - a))
        else:
            spec = o.mask.spec
            if spec == "full":
                return u, E, dv, W
            x, y = xc[:, 0], xc[:, 1]
            if spec[0] in ("circ_in", "circ_out"):
                cx, cy, R = spec[1:4]
                dx, dy = x - cx, y - cy
                r = np.maximum(np.hypot(dx, dy), 1e-12)
                ex, ey = dx / r, dy / r
                if spec[0] == "circ_in":
                    u = R - r
                    dv = np.column_stack([-ex, -ey, np.zeros(n)])
                    W[:] = min(Ue, R)
                    self._curv[:] = -1.0 / R
                else:
                    u = r - R
                    dv = np.column_stack([ex, ey, np.zeros(n)])
                    W[:] = Ue
                    self._curv[:] = 1.0 / R
                E = np.column_stack([cx + R * ex, cy + R * ey, xc[:, 2]])
            else:
                for x0, x1, y0, y1 in spec[1]:
                    inside = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
                    for d, ev, dvec, wid in ((x - x0, (x0, None), (1, 0), x1 - x0), (x1 - x, (x1, None), (-1, 0), x1 - x0),
                                             (y - y0, (None, y0), (0, 1), y1 - y0), (y1 - y, (None, y1), (0, -1), y1 - y0)):
                        m = inside & (d < u)
                        u[m] = d[m]
                        Ef = np.column_stack([x if ev[0] is None else np.full(n, ev[0]),
                                              y if ev[1] is None else np.full(n, ev[1]), xc[:, 2]])
                        E[m] = Ef[m]
                        dv[m] = np.array([dvec[0], dvec[1], 0.0])
                        W[m] = min(Ue, 0.5 * wid)
        u = np.where(u < W, u, np.inf)
        return u, E, dv, W

    def _edge_z(self, o, xc):
        u, E, dv, W = self._edge(o, xc)
        z = np.isfinite(u)
        return u, E, dv, W, z, self._curv

    def _transport(self, u, E, dv, W, nside, phi_in_fn, field_out, curv=None):
        """Monotone (rank-keeping) map across an edge zone.  For each photon: the cut E + v·dv, v in [-W, W];
        F_in = cumulative incident flux on the opening side (v in [0, W]); F_out = cumulative outgoing j_z at δ
        (both sides of the edge).  Survival p = F_out(W)/F_in(W) (<= 1); new position: F_out(v*) = p·F_in(u)."""
        delta = self.r["delta"]
        M = len(u)
        tt = np.linspace(-1, 1, 161)
        g = np.sinh(6 * tt) / np.sinh(6)                      # dense near the edge
        v = W[:, None] * g[None, :]                            # (M, NC)
        Pout = E[:, None, :] + v[..., None] * dv[:, None, :] + delta * nside[:, None, :]
        jout = field_out(Pout.reshape(-1, self.dim), M).reshape(M, -1, self.dim)
        jz = np.maximum((jout * nside[:, None, :]).sum(-1), 0)
        cv = np.zeros(M) if curv is None else curv
        jz = jz * (1 + cv[:, None] * v)                       # area element of a circular edge
        Fo = np.concatenate([np.zeros((M, 1)), np.cumsum(0.5 * (jz[:, 1:] + jz[:, :-1]) * np.diff(v, axis=1), 1)], 1)
        pos = g >= 0
        vp = v[:, pos]
        Pin = E[:, None, :] + vp[..., None] * dv[:, None, :]
        fin = np.maximum(phi_in_fn(Pin.reshape(-1, self.dim), M).reshape(M, -1), 0) * (1 + cv[:, None] * vp)
        Fi = np.concatenate([np.zeros((M, 1)), np.cumsum(0.5 * (fin[:, 1:] + fin[:, :-1]) * np.diff(vp, axis=1), 1)], 1)
        p = np.clip(Fo[:, -1] / np.maximum(Fi[:, -1], 1e-300), 0, 1)
        newpos = np.zeros_like(E)
        for k in range(M):
            target = p[k] * np.interp(u[k], vp[k], Fi[k])
            vs = np.interp(target, Fo[k], v[k])
            newpos[k] = E[k] + vs * dv[k] + delta * nside[k]
        return newpos, p

    def _tog(self, n):
        if self.tog is None:
            return np.full(n, np.inf)
        return np.broadcast_to(self.tog, (n,)).astype(float)

    # ---- objects hit along a chord ----------------------------------------------------------
    def _hits(self, Xa, Xb, t, only=None, tg=None):
        """First object crossed by the chord Xa→Xb (sources and passive planes excluded, unless only=i).
        Returns (idx or −1, fraction of the chord)."""
        n = len(Xa)
        best = np.full(n, np.inf)
        idx = np.full(n, -1)
        for i, o in enumerate(self.sc.objs):
            if only is None and o.kind in ("source", "plane"):
                continue
            if only is not None and i != only:
                continue
            if self.dim == 2:
                p1, p2 = o.p1(), o.p2()
                r = Xb - Xa
                s = p2 - p1
                den = r[:, 0] * s[1] - r[:, 1] * s[0]
                q = p1[None, :] - Xa
                with np.errstate(divide="ignore", invalid="ignore"):
                    tt = (q[:, 0] * s[1] - q[:, 1] * s[0]) / den
                    uu = (q[:, 0] * r[:, 1] - q[:, 1] * r[:, 0]) / den
                ok = (np.abs(den) > 1e-300) & (tt >= 0) & (tt <= 1) & (uu >= 0) & (uu <= 1)
                if self.tog is not None and i == self.sc.tog_obj:
                    tc = t + tt * np.sqrt(((Xb - Xa) ** 2).sum(-1))
                    pres = np.where(tc < (np.full(n, np.inf) if tg is None else tg), self.sc.tog_init, not self.sc.tog_init)
                    ok &= pres
            else:
                dz = Xb[:, 2] - Xa[:, 2]
                with np.errstate(divide="ignore", invalid="ignore"):
                    tt = (o.z - Xa[:, 2]) / dz
                ok = (dz > 0) & (tt >= 0) & (tt <= 1)
            better = ok & (tt < best)
            best[better] = tt[better]
            idx[better] = i
        return idx, best

    def _ray_dist(self, X, D, t, tg=None):
        """Distance along the current direction to the first object ahead (for the step cap)."""
        far = X + 1e6 * D
        idx, fr = self._hits(X, far, t, tg=tg)
        return np.where(idx >= 0, fr * 1e6, np.inf)

    # ---- the run ------------------------------------------------------------------------------
    def run(self, X, D, t0=None, hear=None, planes=(), on_event=None):
        """Move photon balls until each clicks / dies.  X, D: (n, dim).  Returns dict with status,
        click positions, object of the click, time, recorded crossings of passive planes."""
        r, dd = self.r, self.dd
        n = len(X)
        X = X.copy()
        D = D.copy()
        t = np.zeros(n) if t0 is None else np.asarray(t0, float).copy()
        hear = np.full(n, -1) if hear is None else np.asarray(hear).copy()
        tog = self._tog(n)
        self._togrun = tog
        st = np.full(n, ALIVE)
        where = np.full(n, -1)
        h = np.full(n, dd["dtouch"])
        steps = np.zeros(n, int)
        rec = []                                   # (photon, plane obj, position)
        events = []                                # (photon, obj, kind, info)
        K1 = np.zeros_like(X)
        k1ok = np.zeros(n, bool)
        while True:
            act = np.nonzero(st == ALIVE)[0]
            if act.size == 0:
                break
            Xa, Da, ta, ha = X[act], D[act], t[act], h[act]
            hz = hear[act]
            tg = tog[act]
            # cap: do not fly more than (distance to the next object + δ) in one step
            cap = self._ray_dist(Xa, Da, ta, tg) + r["delta"]
            ha = np.minimum(np.minimum(ha, cap), r["hmax"])
            need = ~k1ok[act]
            if need.any():
                K1[act[need]] = self.vel(Xa[need], ta[need], Da[need], hz[need], tg[need])
            ks = [K1[act]]
            for s in range(1, 7):
                Y = Xa + ha[:, None] * sum(a * kk for a, kk in zip(_A[s], ks))
                ks.append(self.vel(Y, ta + _C[s] * ha, Da, hz, tg))
            y5 = Xa + ha[:, None] * sum(b * kk for b, kk in zip(_B5, ks))
            y4 = Xa + ha[:, None] * sum(b * kk for b, kk in zip(_B4, ks))
            err = np.sqrt(((y5 - y4) ** 2).sum(-1))
            acc = (err <= r["atol"]) | (ha <= r["hmin"] * 1.0001)
            fac = np.clip(0.9 * (r["atol"] / np.maximum(err, 1e-300)) ** 0.2, 0.2, 5.0)
            hnew = np.clip(ha * fac, r["hmin"], r["hmax"])
            # accepted steps: check crossings of the chord
            ai = act[acc]
            Xo = Xa[acc]
            Xn = y5[acc]
            tn = ta[acc] + ha[acc]
            idx, fr = self._hits(Xo, Xn, ta[acc], tg=tg[acc])
            # passive planes: record the crossing, do not touch the photon
            for pl in self._planes:
                hp, fp = self._hits(Xo, Xn, ta[acc], only=pl, tg=tg[acc])
                for q in np.nonzero((hp >= 0) & (fp <= fr))[0]:
                    rec.append((ai[q], pl, Xo[q] + fp[q] * (Xn[q] - Xo[q])))
            free = idx < 0
            ii = ai[free]
            X[ii] = Xn[free]
            D[ii] = ks[-1][acc][free]
            t[ii] = tn[free]
            steps[ii] += 1
            K1[ii] = ks[-1][acc][free]
            k1ok[ii] = True
            h[act] = hnew
            # crossings → events
            hit = np.nonzero(~free)[0]
            if hit.size:
                pi = ai[hit]
                Xc = Xo[hit] + fr[hit, None] * (Xn[hit] - Xo[hit])
                chord = Xn[hit] - Xo[hit]
                din = chord / np.maximum(np.sqrt((chord ** 2).sum(-1)), 1e-300)[:, None]
                tc = ta[acc][hit] + fr[hit] * ha[acc][hit]
                k1ok[pi] = False
                self._events(pi, idx[hit], Xc, din, tc, X, D, t, st, where, hear, h, events)
            # bookkeeping
            out = np.zeros(n, bool)
            for a in range(self.dim):
                out |= (X[:, a] < self.box[2 * a]) | (X[:, a] > self.box[2 * a + 1])
            st[(st == ALIVE) & out] = LOST
            st[(st == ALIVE) & (steps > self.max_steps)] = STUCK
        return dict(status=st, X=X, where=where, t=t, rec=rec, events=events, steps=steps)

    def _events(self, pi, oi, Xc, din, tc, X, D, t, st, where, hear, h, events):
        delta = self.r["delta"]
        for o_idx in np.unique(oi):
            sel = oi == o_idx
            ps, xc, dn, tt = pi[sel], Xc[sel], din[sel], tc[sel]
            o = self.sc.objs[o_idx]
            if o.kind == "screen":
                X[ps] = xc
                t[ps] = tt
                st[ps] = CLICK
                where[ps] = o_idx
                continue
            if o.kind == "stop":
                X[ps] = xc
                st[ps] = ABSORBED
                where[ps] = o_idx
                continue
            if self.dim == 2:
                nrm = o.n
                s = (xc - o.c) @ o.t
                tr = o.transmits(s)
                cosn = dn @ nrm
            else:
                nrm = np.array([0, 0, 1.0])
                tr = o.mask.inside(xc[:, 0], xc[:, 1]) & (o.T is not None)
                cosn = dn[:, 2]
            # opaque part (or mirror: nothing on the other side) → the photon dies
            dead = ~tr
            X[ps[dead]] = xc[dead]
            st[ps[dead]] = ABSORBED
            where[ps[dead]] = o_idx
            ps, xc, dn, tt, cosn = ps[tr], xc[tr], dn[tr], tt[tr], cosn[tr]
            if ps.size == 0:
                continue
            sgn = np.sign(cosn)
            Xd = xc + (delta / np.maximum(np.abs(cosn), 1e-3))[:, None] * dn
            if self.dim == 3:
                newhear = self.sc.fan_of.get(o_idx, -1)
            if o.absorbing:
                if self.dim == 2:
                    jin, _, _ = self.sc.current(xc - (sgn * 1e-7)[:, None] * nrm[None, :], tt, self._togrun[ps], excl=o_idx)
                    jout, _, _ = self.sc.current(Xd, tt + delta, self._togrun[ps])
                    nn = sgn[:, None] * nrm[None, :]
                else:
                    jin, _, _ = self.sc.current(xc - np.array([0, 0, 1e-7]), hear[ps])
                    jout, _, _ = self.sc.current(Xd, newhear)
                    nn = np.repeat(nrm[None, :], len(ps), 0)
                fin = (jin * nn).sum(-1)
                fout = (jout * nn).sum(-1)
                p = np.where(fin > 0, np.clip(fout / np.maximum(fin, 1e-300), 0, 1), (fout > 0) * 1.0)
                u, E, dvv, W, zz, cvv = self._edge_z(o, xc)
                if zz.any():
                    zi = np.nonzero(zz)[0]
                    if self.dim == 2:
                        def fi_(P, M, zi=zi):
                            r_ = len(P) // M
                            nn_ = np.repeat(nn[zi], r_, 0)
                            tg_ = np.repeat(self._togrun[ps[zi]], r_)
                            j_, _, _ = self.sc.current(P - 1e-7 * nn_, np.repeat(tt[zi], r_), tg_, excl=o_idx)
                            return (j_ * nn_).sum(-1)

                        def fo_(P, M, zi=zi):
                            r_ = len(P) // M
                            tg_ = np.repeat(self._togrun[ps[zi]], r_)
                            return self.sc.current(P, np.repeat(tt[zi] + delta, r_), tg_)[0]
                    else:
                        def fi_(P, M, zi=zi):
                            P = P.copy()
                            P[:, 2] = o.z
                            j_, _, _ = self.sc.current(P, np.repeat(hear[ps[zi]], len(P) // M))
                            return j_[:, 2]

                        def fo_(P, M):
                            return self.sc.current(P, newhear)[0]
                    Xd[zz], p[zz] = self._transport(u[zz], E[zz], dvv[zz], W[zz], nn[zz], fi_, fo_, cvv[zz])
                ok = self.rng.random(len(ps)) < p
                st[ps[~ok]] = COIN
                where[ps[~ok]] = o_idx
                X[ps[~ok]] = xc[~ok]
                events.append((ps, o_idx, p))
                ps, Xd, dn, tt = ps[ok], Xd[ok], dn[ok], tt[ok]
            X[ps] = Xd
            D[ps] = dn
            t[ps] = tt + delta
            h[ps] = self.dd["dtouch"]
            if self.dim == 3:
                hear[ps] = newhear


# ───────────────────────────── numpy reference of the field sums (for base_test.py) ───────────
def _blend_np(l, R, r0, lmax):
    u = np.log2(r0 / R)
    c = np.ceil(u)
    f = u - (c - 1)
    s = f * f * (3 - 2 * f)
    a = np.minimum(c, lmax)
    b = np.minimum(c + 1, lmax)
    w = np.where(a == l, 1 - s, 0) + np.where(b == l, s, 0)
    return np.where(u <= -1, 1.0 * (l == 0), w)


def _cr(f):
    return [(-f ** 3 + 2 * f ** 2 - f) / 2, (3 * f ** 3 - 5 * f ** 2 + 2) / 2, (-3 * f ** 3 + 4 * f ** 2 + f) / 2, (f ** 3 - f ** 2) / 2]


def origins2_np(f, ds, sa=-1e300, sb=1e300):
    """Same as Fan2::origins (numpy): cell centroids inside the transmitting intervals, amplitude B·length."""
    lo, hi = max(sa, f["glo"]), min(sb, f["ghi"])
    if hi <= lo:
        return np.zeros((0, 2)), np.zeros((0, 2), complex)
    k = np.arange(max(0, int(np.floor((lo - f["glo"]) / ds))), int(np.ceil((hi - f["glo"]) / ds)))
    a = f["glo"] + k * ds
    b = np.minimum(a + ds, f["ghi"])
    ln = np.zeros(a.size)
    mom = np.zeros(a.size)
    for x0, x1 in f["ints"]:
        L_ = np.maximum(a, x0)
        R_ = np.minimum(b, x1)
        ov = np.clip(R_ - L_, 0, None)
        ln += ov
        mom += ov * 0.5 * (L_ + R_)
    ok = ln > 0
    s = mom[ok] / ln[ok]
    u = (s - f["sp0"]) / f["sds"]
    j = np.floor(u).astype(int)
    w = _cr(u - j)
    B = sum(w[m][:, None] * f["samp"][np.clip(j - 1 + m, 0, len(f["samp"]) - 1)] for m in range(4))
    good = np.sqrt((np.abs(B) ** 2).sum(-1)) >= f["bthr"]
    s, B, L_ = s[good], B[good], ln[ok][good]
    return f["c"][None, :] + s[:, None] * f["t"][None, :], B * L_[:, None]


def field2_numpy(scene, P):
    """Numpy reference of Scene2.field with everything present (no gating); used by base_test.py."""
    from scipy.special import hankel1
    d = scene.d
    out = np.zeros((len(P), 6), complex)
    for f in scene.fans:
        for l in range(d["lmax"] + 1):
            O, a = origins2_np(f, d["ds0"] / 2 ** l)
            if len(O) == 0:
                continue
            r = P[:, None, :] - O[None, :, :]
            zn = r @ f["n"]
            R = np.sqrt((r ** 2).sum(-1))
            w = _blend_np(l, R, d["r0"], d["lmax"]) * (zn > 1e-12)
            # blocking
            for i, o in enumerate(scene.objs):
                if not o.blocks or i == f["obj"]:
                    continue
                p1, p2 = o.p1(), o.p2()
                rr = r                                  # segment O→P
                s = p2 - p1
                den = rr[..., 0] * s[1] - rr[..., 1] * s[0]
                q = p1[None, None, :] - O[None, :, :]
                with np.errstate(divide="ignore", invalid="ignore"):
                    tt = (q[..., 0] * s[1] - q[..., 1] * s[0]) / den
                    uu = (q[..., 0] * rr[..., 1] - q[..., 1] * rr[..., 0]) / den
                blk = (np.abs(den) > 1e-300) & (tt > 1e-12) & (tt < 1 - 1e-12) & (uu >= 0) & (uu <= 1)
                w = np.where(blk, 0, w)
            Rs = np.where(w > 0, R, 1.0)
            H0 = hankel1(0, K * Rs)
            H1 = hankel1(1, K * Rs)
            F = H1 / Rs
            Fp = (K * H0 - 2 * H1 / Rs) / Rs
            c = 0.5j * K * w
            G = c * zn * F
            Gx = c * (f["n"][0] * F + zn * Fp * r[..., 0] / Rs)
            Gz = c * (f["n"][1] * F + zn * Fp * r[..., 1] / Rs)
            out[:, 0:2] += G @ a
            out[:, 2:4] += Gx @ a
            out[:, 4:6] += Gz @ a
    return out[:, 0:2], out[:, 2:4], out[:, 4:6]


def field3_numpy(scene, P, fi):
    """Numpy reference of Scene3.field (levels, blend, vector kernel) using Rust's origin lists."""
    d = scene.d
    f = scene.fans[fi]
    E = np.zeros((len(P), 3), complex)
    dE = np.zeros((len(P), 3, 3), complex)
    for l in range(d["lmax"] + 1):
        cap = 4_000_000
        buf = np.zeros(cap * 6)
        nn = lib().fan_origins(scene.ptr, fi, l, -1e300, 1e300, -1e300, 1e300, _p(buf), cap)
        if nn == 0 or nn > cap:
            continue
        b = buf[:6 * nn].reshape(-1, 6)
        O = b[:, :2]
        a = np.column_stack([b[:, 2] + 1j * b[:, 3], b[:, 4] + 1j * b[:, 5]])
        r = np.concatenate([P[:, None, :2] - O[None, :, :], np.broadcast_to((P[:, 2] - f["z0"])[:, None, None], (len(P), len(O), 1))], -1)
        R = np.sqrt((r ** 2).sum(-1))
        w = _blend_np(l, R, d["r0"], d["lmax"])
        e = np.exp(1j * K * R) / (4 * np.pi * R) * w
        h = (1j * K - 1 / R) * e
        h2 = (-K ** 2 - 2j * K / R + 2 / R ** 2) * e
        hr = h / R
        q = (h2 - hr) / R ** 2
        g = hr[..., None] * r
        H = q[..., None, None] * r[..., :, None] * r[..., None, :] + hr[..., None, None] * np.eye(3)
        bx, by = a[:, 0], a[:, 1]
        E[:, 0] += -2 * g[..., 2] @ bx
        E[:, 1] += -2 * g[..., 2] @ by
        E[:, 2] += 2 * (g[..., 0] @ bx + g[..., 1] @ by)
        for bb in range(3):
            dE[:, 0, bb] += -2 * H[..., bb, 2] @ bx
            dE[:, 1, bb] += -2 * H[..., bb, 2] @ by
            dE[:, 2, bb] += 2 * (H[..., bb, 0] @ bx + H[..., bb, 1] @ by)
    return E, dE


# ───────────────────────────── small helpers for the experiment scripts ─────────────────────────
def log(*a):
    print(*a, flush=True)


def result_line(name, measured, theory, tol, ok, extra=""):
    return dict(name=name, measured=measured, theory=theory, tol=tol, ok=bool(ok), extra=extra)


def save_results(tag, rows, meta):
    import json
    path = os.path.join(OUT, f"base_{tag}_{RULE['version']}.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(dict(rule=RULE, derived=derived(), meta=meta, rows=rows), fh, indent=1, default=str)
    log(f"saved {path}")


def simulate(scene, src, n, box, seed, batch=2000, tog_fn=None, max_steps=20000, verbose=True):
    """Emit n photons from the source object `src`, run them, return what happened.
    tog_fn(rng, m) -> per-photon toggle times (Wheeler), optional."""
    rng = np.random.default_rng(seed)
    t0 = time.time()
    res = dict(n=0, src_coin=0, clicks=[], where=[], status=[], rec=[], t=[], tog=[], wall_frac=None, steps=[])
    done = 0
    while done < n:
        m = min(batch, n - done)
        ph = Photons(scene, box, rng, max_steps=max_steps)
        tog = None if tog_fn is None else tog_fn(rng, m)
        ph.tog = tog
        X, D, surv, wall, nout = ph.birth(src, m)
        res["wall_frac"] = wall
        res["n"] += len(X)
        res["src_coin"] += int((~surv).sum())
        X, D = X[surv], D[surv]
        tg = None if tog is None else np.asarray(tog)[ph._keep][surv]
        ph.tog = tg
        hear = None
        if ph.dim == 3:
            hear = np.full(len(X), scene.fan_of[src])
        out = ph.run(X, D, t0=np.full(len(X), scene.rule["delta"]), hear=hear)
        res["clicks"].append(out["X"])
        res["where"].append(out["where"])
        res["status"].append(out["status"])
        res["t"].append(out["t"])
        res["steps"].append(out["steps"])
        res["tog"].append(tg if tg is not None else np.full(len(X), np.inf))
        base = sum(len(a) for a in res["status"][:-1])
        res["rec"] += [(base + a, b, c) for a, b, c in out["rec"]]
        done += m
        if verbose:
            st = np.concatenate(res["status"])
            log(f"    ... {done}/{n} photons, {int((st == CLICK).sum())} clicks, {time.time() - t0:.0f}s")
    for k in ("clicks", "where", "status", "t", "tog", "steps"):
        res[k] = np.concatenate(res[k]) if res[k] else np.zeros(0)
    res["time"] = time.time() - t0
    st = res["status"]
    res["counts"] = {STATUS[s]: int((st == s).sum()) for s in STATUS}
    res["counts"]["source coin"] = res["src_coin"]
    return res

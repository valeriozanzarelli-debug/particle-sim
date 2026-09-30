"""THE PHOTON RULE — one local rule per ball (2-D lab version, x transverse, z along the beam).

A photon = N wave balls + 1 photon ball.

WAVE BALL (every ball, always the same rule)
  * flies straight at c; its phase grows by 2π per λ of path;
  * carries: origin (point of its last emission), direction u = tanθ, amplitude a (Jones vector,
    complex) fixed at emission;
  * when it touches an object (source, slit wall, polarizer/PBS) it is (re)emitted: the object
    turns it into a fan of balls with a regular set of directions, amplitude = local field × transmission.
    Balls that hit an opaque part of an object die.

PHOTON BALL (one of the balls, flagged at emission)
  * at each tick it reads ONLY the wave balls around it: for each origin, the two balls of that origin
    whose rays bracket it (a "hat" window whose width = local spacing of that fan);
  * each such ball contributes its wavelet at the photon's position
      E2 ("wavelet"):  a · G(r − origin)            (exact 2-D Huygens kernel, exact gradient)
      E1 ("plane"):    a · G(ℓ) · e^{ik n·d}        (ball's own phase + plane-wave extrapolation; ∇ = ik n)
  * ψ, ∇ψ = sums of these (Jones components kept separate);
    current j = Σ_pol Im(ψ* ∇ψ);  the photon ball moves at c along j  (dx/dz = j_x / j_z);
  * if there is no neighbour (isolated ball) it keeps its own direction;
  * it dies where it hits an opaque object; the click is where it hits the screen.

Everything below is vectorized over many independent photons (each photon still only uses its own
neighbourhood). Emission of the wave balls is a regular lattice (see REGOLA_FOTONE.md for why).
"""
import numpy as np, os, time, sys
from multiprocessing import Pool
from optics import *

HERE = os.path.dirname(os.path.abspath(__file__))
TH_FULL = np.pi / 2 - 1e-3            # a fan covers the whole forward half-plane


# ─────────────────────────────── wave balls ────────────────────────────────────────────
class Fan:
    """All wave balls emitted from the plane z = ze.

    lattice mode: origins s_i (n_src) with amplitude a_i (n_src, ncomp); every origin emits the same
                  regular set of directions θ_j = −TH + (j+½)dth, j < nu  (TH ≈ π/2: the whole
                  forward half-plane).  N_balls = n_src·nu.
    mc mode     : N random balls (origin, direction) — the old pilot / current Rust emission — each with
                  its own Monte-Carlo amplitude.
    """

    def __init__(self, s, a, ze=0.0, TH=TH_FULL, dth=0.01, mode="lattice", u=None):
        self.s = np.asarray(s, float); self.a = np.asarray(a, complex)
        if self.a.ndim == 1:
            self.a = self.a[:, None]
        self.ze, self.TH, self.dth, self.mode = ze, TH, dth, mode
        self.nu = int(round(2 * TH / dth)) if mode == "lattice" else 1
        self.u = u                                         # mc: direction of each ball

    @property
    def n_balls(self):
        return self.s.size * (self.nu if self.mode == "lattice" else 1)


def lattice_fan(s_lo, s_hi, ds, amp_fn, ze=0.0, TH=None, dth=0.01, jones=None):
    """Regular origins on [s_lo, s_hi) (several intervals allowed), amplitude amp_fn(s)*ds*jones."""
    ss, aa = [], []
    for k, (lo, hi) in enumerate(zip(np.atleast_1d(s_lo), np.atleast_1d(s_hi))):
        n = max(1, int(round((hi - lo) / ds)))
        s = lo + (np.arange(n) + 0.5) * (hi - lo) / n
        a = amp_fn(s) * (hi - lo) / n
        J = np.array([1.0 + 0j]) if jones is None else np.asarray(jones[k], complex)
        ss.append(s); aa.append(a[:, None] * J[None, :])
    return Fan(np.concatenate(ss), np.concatenate(aa), ze, TH_FULL if TH is None else TH, dth)


def mc_fan(n, rng, theta_max=0.45, jones=None):
    """Old pilot / Rust emission at the double slit: balls sampled from |U_inc|^2, those in a slit
    re-emitted with θ uniform in ±theta_max, amplitude obliquity/inc (+ normalization)."""
    x = rng.normal(0, W0 / 2, n)
    s = x[in_slit(x)]
    th = rng.uniform(-theta_max, theta_max, s.size)
    u = np.tan(th)
    # unbiased weights: a = U0(s) / (N p(s) p_u(u)),  p(s) ∝ inc^2 on the whole line, p_u = cos^2θ/(2θm)
    ps = np.exp(-2 * s ** 2 / W0 ** 2) / (W0 / 2 * np.sqrt(2 * np.pi))
    pu = np.cos(th) ** 2 / (2 * theta_max)
    a = inc_amp(s) / (n * ps * pu)
    if jones is None:
        A = a[:, None].astype(complex)
    else:
        J = np.where((s > 0)[:, None], np.asarray(jones[1])[None, :], np.asarray(jones[0])[None, :])
        A = a[:, None] * J
    return Fan(s, A, 0.0, theta_max, None, mode="mc", u=u)


# ─────────────────────── what the photon ball reads from its neighbours ─────────────────
def kernel(dx, dz):
    """Exact 2-D Rayleigh-Sommerfeld wavelet G and its gradient, for offsets (dx, dz) from the origin."""
    R = np.sqrt(dx ** 2 + dz ** 2)
    kR = K * R
    H0 = hank(0, kR); H1 = hank(1, kR)
    F = H1 / R
    Fp = (K * H0 - 2 * H1 / R) / R
    c = 0.5j * K
    return c * dz * F, c * dz * Fp * dx / R, c * (F + dz * Fp * dz / R)


def ball_value(th, dz):
    """What a ball carries on its own ray: G at its own position (distance ℓ from origin, angle θ)."""
    cos = np.cos(th)
    ell = dz / cos
    return 0.5j * K * cos * hank(1, K * ell), cos


def field_lattice(X, z, fan, est="E2"):
    """Hat window: for every origin, the two balls whose rays bracket the photon."""
    dz = z - fan.ze
    t = np.arctan2(X[:, None] - fan.s[None, :], dz)       # direction from origin to photon
    q = (t + fan.TH) / fan.dth - 0.5
    j = np.floor(q); f = q - j
    w0 = np.where((j >= 0) & (j < fan.nu), 1 - f, 0.0)
    w1 = np.where((j + 1 >= 0) & (j + 1 < fan.nu), f, 0.0)
    if est == "E2":
        G, Gx, Gz = kernel(X[:, None] - fan.s[None, :], dz)
        W = w0 + w1
        return (W * G) @ fan.a, (W * Gx) @ fan.a, (W * Gz) @ fan.a
    U = Ux = Uz = 0
    for jj, ww in ((j, w0), (j + 1, w1)):
        th = -fan.TH + (jj + 0.5) * fan.dth
        g, cos = ball_value(th, dz)
        sin = np.sin(th)
        d = X[:, None] - (fan.s[None, :] + dz * np.tan(th))  # offset photon − ball (same z)
        e = ww * g * np.exp(1j * K * sin * d)
        U = U + e @ fan.a
        Ux = Ux + (1j * K * sin * e) @ fan.a
        Uz = Uz + (1j * K * cos * e) @ fan.a
    return U, Ux, Uz


def field_mc(X, z, fan, est="E2", h=1.0):
    """Gaussian window of width h over random balls (old pilot style)."""
    dz = z - fan.ze
    h = h(dz) if callable(h) else h
    xr = fan.s + dz * fan.u
    order = np.argsort(xr); xs = xr[order]
    lo = np.searchsorted(xs, X - 3 * h); hi = np.searchsorted(xs, X + 3 * h)
    cnt = hi - lo; tot = cnt.sum()
    nc = fan.a.shape[1]
    out = [np.zeros((X.size, nc), complex) for _ in range(3)]
    if tot == 0:
        return out
    pid = np.repeat(np.arange(X.size), cnt)
    b = order[np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt) + np.repeat(lo, cnt)]
    d = X[pid] - xr[b]
    w = np.exp(-0.5 * (d / h) ** 2) * dz / (np.sqrt(2 * np.pi) * h)
    if est == "E2":
        G, Gx, Gz = kernel(X[pid] - fan.s[b], dz)
    else:
        th = np.arctan(fan.u[b])
        g, cos = ball_value(th, dz)
        sin = np.sin(th)
        G = g * np.exp(1j * K * sin * d)
        Gx = 1j * K * sin * G; Gz = 1j * K * cos * G
    for o, gg in zip(out, (G, Gx, Gz)):
        v = (w * gg)[:, None] * fan.a[b]
        for c in range(nc):
            o[:, c] = np.bincount(pid, v[:, c].real, X.size) + 1j * np.bincount(pid, v[:, c].imag, X.size)
    return out


def slope(X, z, fan, est, own_u, h=None):
    """dx/dz of the photon ball = j_x / j_z; own direction where no neighbour is heard."""
    if fan.mode == "lattice":
        U, Ux, Uz = field_lattice(X, z, fan, est)
    else:
        U, Ux, Uz = field_mc(X, z, fan, est, h)
    jx = np.imag(np.conj(U) * Ux).sum(-1)
    jz = np.imag(np.conj(U) * Uz).sum(-1)
    ok = jz > 1e-12 * (np.abs(U) ** 2).sum(-1).max(initial=1e-300) * K
    return np.where(ok, jx / np.where(ok, jz, 1), own_u)


def guide(X, fan, z0, z1, est="E2", dzmax=2.0, grow=0.1, dzmin=0.01, h=None, rk=4):
    """Move photon balls from z0 to z1. Returns final x and final own direction."""
    X = X.copy(); z = z0
    u = np.zeros_like(X)
    while z < z1 - 1e-9:
        dz = min(dzmax, max(dzmin, grow * (z - fan.ze)), z1 - z)
        if rk == 1:
            u = slope(X, z, fan, est, u, h); X = X + dz * u
        else:
            k1 = slope(X, z, fan, est, u, h)
            k2 = slope(X + 0.5 * dz * k1, z + 0.5 * dz, fan, est, k1, h)
            k3 = slope(X + 0.5 * dz * k2, z + 0.5 * dz, fan, est, k2, h)
            k4 = slope(X + dz * k3, z + dz, fan, est, k3, h)
            u = (k1 + 2 * k2 + 2 * k3 + k4) / 6
            X = X + dz * u
        z += dz
    return X, u


# ─────────────────────────────── experiments ──────────────────────────────────────────
def slit_fan(cfg):
    jones = cfg.get("jones")
    if cfg.get("mode", "lattice") == "mc":
        return mc_fan(cfg["n_mc"], np.random.default_rng(cfg["seed"]), cfg.get("theta_max", 0.45), jones)
    return lattice_fan([c - A / 2 for c in SLITS], [c + A / 2 for c in SLITS], cfg["ds"], inc_amp,
                       TH=cfg["TH"], dth=cfg["dth"], jones=jones)


def local_field(X, z, fan, est="E2", h=None):
    return field_lattice(X, z, fan, est) if fan.mode == "lattice" else field_mc(X, z, fan, est, h)


def edge_rule(X, fan, delta, inc, est, rng, h=None):
    """Photon ball just crossed the aperture at X (flown straight to z = delta).
    It survives with probability min(1, j_z / (k |psi_inc|^2)), both read locally:
    j_z from the neighbour balls of the fan, |psi_inc|^2 from the incident balls at the wall."""
    U, Ux, Uz = local_field(X, fan.ze + delta, fan, est, h)
    jz = np.imag(np.conj(U) * Uz).sum(-1)
    r = jz / (K * inc ** 2)
    return rng.random(X.size) < np.minimum(1.0, r), r


def first_step(x0, fan, delta, est="E2", h=None, iters=2):
    """First step of length delta behind an object: along the OUTGOING current, evaluated at the end
    of the step (implicit Euler). 'straight' (old) = keep the incoming direction."""
    x1 = x0.copy()
    for _ in range(iters):
        x1 = x0 + delta * slope(x1, fan.ze + delta, fan, est, np.zeros_like(x0), h)
    return x1


def _slit_worker(args):
    x0, cfg = args
    fan = slit_fan(cfg)
    rng = np.random.default_rng(cfg["seed"] + 17)
    x1 = first_step(x0, fan, cfg["z0"], cfg["est"], cfg.get("h")) if cfg.get("first") == "implicit" else x0
    if cfg.get("edge", True):
        keep, r = edge_rule(x1, fan, cfg["z0"], inc_amp(x0), cfg["est"], rng, cfg.get("h"))
    else:
        keep = np.ones(x0.size, bool)
    xl, _ = guide(x1[keep], fan, cfg["z0"], L, cfg["est"], cfg["dzmax"], cfg["grow"], cfg["dzmin"],
                  cfg.get("h"), cfg.get("rk", 4))
    return xl, (~keep).sum()


DEFAULT = dict(est="E2", mode="lattice", ds=0.02, TH=TH_FULL, dth=0.01, z0=0.1, dzmax=2.0, grow=0.1, dzmin=0.01,
               nph=6000, seed=7, rk=4, edge=True)


def double_slit(cfg, pool, chunks=32):
    """Returns click positions, fraction absorbed (wall + edge rule), wall-time."""
    c = dict(DEFAULT); c.update(cfg)
    rng = np.random.default_rng(c["seed"])
    # the source emits the photon ball with |U_inc|^2 (quantum equilibrium); the wall absorbs it unless in a slit
    n_enter = int(c["nph"] * (1.1 if c["edge"] else 1.0))
    x = rng.normal(0, W0 / 2, int(n_enter / 0.2) + 2000)
    through = in_slit(x)
    x0 = x[through][:n_enter]
    n_emit = np.searchsorted(np.cumsum(through), n_enter) + 1
    t = time.time()
    subs = []
    for i, part in enumerate(np.array_split(x0, chunks)):
        ci = dict(c); ci["seed"] = c["seed"] * 1000 + i          # mc: each chunk has its own random wave
        subs.append((part, ci))
    res = pool.map(_slit_worker, subs)
    xl = np.concatenate([r[0] for r in res])
    edge_abs = sum(r[1] for r in res)
    absorbed = 1 - xl.size / n_emit
    return xl, absorbed, time.time() - t


if __name__ == "__main__":
    # the final rule on the double slit (6000 photons entering, ~2 min on 8 cores)
    refs = screen_reference("flux")
    b = ideal_baseline(refs, 6000, 200)
    print(f"ideal sampler 6000 clicks: corr {b[0]:.4f}+-{b[1]:.4f}  fringe rms {100*b[2]:.1f}%  max {100*b[3]:.1f}%", flush=True)
    with Pool(8) as pool:
        xl, ab, dt = double_slit({}, pool)
        report("double slit, final rule (lattice fan, hat window, E2, edge rule)", xl, refs,
               f"absorbed {100*ab:.1f}% (theory: geometric 76.2%, RS model 78.1%)  [{dt:.0f}s]")

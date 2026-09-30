"""Shared optics toolkit for the photon-rule lab (2-D problem: x transverse, z along the beam).

Everything here is REFERENCE physics or MEASUREMENT, never part of the ball rule:
  * exact 2-D Rayleigh-Sommerfeld field behind the double slit (Hankel kernel)
  * reference screen patterns (flux j_z = what a flat detector counts, |U|^2, paraxial Fresnel)
  * scoring: coarse correlation, fringe-height errors, ideal-sampler baseline
"""
import sys
import numpy as np
from scipy.special import hankel1

try:                                    # Windows console: allow λ, ψ, … in prints
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

LAM = 1.0
K = 2 * np.pi / LAM
W0, D, A, L = 7.0, 8.0, 2.0, 160.0          # beam waist, slit separation, slit width, slit->screen
BINS, HALF, COARSE = 240, 60.0, 8           # screen histogram; coarse score = 8 fine bins merged
EDGES = np.linspace(-HALF, HALF, BINS + 1)
XS = 0.5 * (EDGES[1:] + EDGES[:-1])
SLITS = (-D / 2, D / 2)


def inc_amp(x):
    """Incident Gaussian amplitude on the slit plane (intensity sigma = W0/2)."""
    return np.exp(-(x ** 2) / W0 ** 2)


def in_slit(x):
    return (np.abs(x - D / 2) <= A / 2) | (np.abs(x + D / 2) <= A / 2)


def sources(m=200, jones=None):
    """Midpoint samples of the aperture field U0(s) = inc(s)·rect, weight ds.
    jones: optional pair (J_left, J_right) of Jones vectors applied per slit (which-path)."""
    ds = A / m
    s = np.concatenate([c - A / 2 + (np.arange(m) + 0.5) * ds for c in SLITS])
    w = inc_amp(s) * ds
    if jones is None:
        return s, w[:, None].astype(complex)          # scalar field: one component
    J = np.concatenate([np.repeat(np.asarray(jones[i], complex)[None, :], m, 0) for i in range(2)])
    return s, w[:, None] * J


# ── Hankel functions: exact near, 4-term asymptotic far (rel. error < 2e-5 for u > 12) ──
_A0 = [1.0, -1 / 8, 9 / 128, -225 / 3072]
_A1 = [1.0, 3 / 8, -15 / 128, 315 / 3072]


def _hank_asym(n, u):
    a = _A0 if n == 0 else _A1
    iu = 1j / u
    ser = a[0] + iu * (a[1] + iu * (a[2] + iu * a[3]))
    return np.sqrt(2 / (np.pi * u)) * np.exp(1j * (u - n * np.pi / 2 - np.pi / 4)) * ser


def hank(n, u):
    out = np.empty(u.shape, complex)
    far = u > 12
    out[far] = _hank_asym(n, u[far])
    if (~far).any():
        out[~far] = hankel1(n, u[~far])
    return out


def rs_field(x, z, s, w):
    """Exact 2-D Rayleigh-Sommerfeld (1st kind) field and gradient at points (x, z>0).
    U = (ik/2) z ∫ U0(s) H1(kR)/R ds.  Returns U, dU/dx, dU/dz with shape (n, ncomp)."""
    x = np.atleast_1d(x); z = np.broadcast_to(z, x.shape)
    dx = x[:, None] - s[None, :]
    zz = z[:, None]
    R = np.sqrt(dx ** 2 + zz ** 2)
    u = K * R
    H0 = hank(0, u); H1 = hank(1, u)
    F = H1 / R
    Fp = (K * H0 - 2 * H1 / R) / R                   # dF/dR
    c = 0.5j * K
    U = c * zz * F
    Ux = c * zz * Fp * dx / R
    Uz = c * (F + zz * Fp * zz / R)
    return U @ w, Ux @ w, Uz @ w


def flux(U, Ux, Uz):
    """Probability current (per unit k), summed over polarization components."""
    jx = np.imag(np.conj(U) * Ux).sum(-1)
    jz = np.imag(np.conj(U) * Uz).sum(-1)
    return jx, jz


# ── reference screen patterns ──────────────────────────────────────────────────
def screen_reference(kind="flux", jones=None, sub=8, m=400):
    """Probability per fine bin on the screen z = L (normalized over the screen)."""
    s, w = sources(m, jones)
    xf = np.linspace(-HALF, HALF, BINS * sub + 1)
    xf = 0.5 * (xf[1:] + xf[:-1])
    U, Ux, Uz = rs_field(xf, L, s, w)
    if kind == "flux":
        p = flux(U, Ux, Uz)[1]
    elif kind == "abs2":
        p = (np.abs(U) ** 2).sum(-1)
    elif kind == "paraxial":
        from rules import fresnel_double_slit
        return fresnel_double_slit(XS, W0, D, A, L)
    p = p.reshape(BINS, sub).sum(1)
    return p / p.sum()


# ── scoring ─────────────────────────────────────────────────────────────────────
def hist(clicks):
    h, _ = np.histogram(clicks, EDGES)
    return h


def coarse_corr(clicks, ref):
    h = hist(clicks)
    return np.corrcoef(h.reshape(-1, COARSE).sum(1), ref.reshape(-1, COARSE).sum(1))[0, 1]


def fringe_windows(ref):
    """Bright fringes with |x| < 45: window = +-quarter period around each maximum."""
    period = LAM * L / D
    fine = np.convolve(ref, np.ones(3) / 3, "same")
    peaks = [i for i in range(1, BINS - 1) if fine[i] >= fine[i - 1] and fine[i] > fine[i + 1]
             and abs(XS[i]) < 45 and fine[i] > 0.05 * fine.max()]
    return [(XS[i] - period / 4, XS[i] + period / 4) for i in peaks]


def fringe_errors(clicks, ref):
    """Relative error of the click fraction in each bright fringe vs theory.
    Returns (array of relative errors, rms, max abs)."""
    clicks = np.asarray(clicks)
    n = (np.abs(clicks) < HALF).sum()          # reference is normalized over the screen
    errs = []
    for a, b in fringe_windows(ref):
        pth = ref[(XS >= a) & (XS < b)].sum()
        pm = ((clicks >= a) & (clicks < b)).sum() / max(n, 1)
        errs.append(pm / pth - 1)
    e = np.array(errs)
    return e, np.sqrt(np.mean(e ** 2)), np.max(np.abs(e))


def visibility(clicks_or_prob, is_prob=False, half=2.5):
    """Central fringe visibility V = (I0 - Id)/(I0 + Id): I0 = counts in |x| < half (central bright
    fringe), Id = mean counts in |x -+ P/2| < half (the two first dark fringes, P = λL/d)."""
    P = LAM * L / D
    if is_prob:
        p = np.asarray(clicks_or_prob)
        cnt = lambda c: p[np.abs(XS - c) < half].sum()
    else:
        x = np.asarray(clicks_or_prob)
        cnt = lambda c: float((np.abs(x - c) < half).sum())
    i0 = cnt(0.0); idk = 0.5 * (cnt(-P / 2) + cnt(P / 2))
    return (i0 - idk) / max(i0 + idk, 1e-30)


def ideal_baseline(ref, n, trials=200, seed=0):
    """What a perfect sampler of `ref` scores with n clicks (mean, std of corr; mean rms / max fringe err)."""
    rng = np.random.default_rng(seed)
    cdf = np.cumsum(ref); cdf /= cdf[-1]
    cs, rm, mx = [], [], []
    for _ in range(trials):
        idx = np.searchsorted(cdf, rng.random(n))
        x = XS[idx] + (rng.random(n) - 0.5) * (2 * HALF / BINS)
        cs.append(coarse_corr(x, ref))
        _, r, m = fringe_errors(x, ref)
        rm.append(r); mx.append(m)
    return np.mean(cs), np.std(cs), np.mean(rm), np.mean(mx)


def chi2(clicks, ref):
    """Pearson chi^2 per degree of freedom on the coarse bins (expected >= 5); ~1 for a perfect sampler."""
    h = hist(clicks).reshape(-1, COARSE).sum(1).astype(float)
    e = ref.reshape(-1, COARSE).sum(1) * h.sum()
    ok = e >= 5
    return ((h[ok] - e[ok]) ** 2 / e[ok]).sum() / max(ok.sum() - 1, 1)


def report(name, clicks, ref, extra=""):
    c = coarse_corr(clicks, ref)
    e, r, m = fringe_errors(clicks, ref)
    print(f"{name}: {len(clicks)} clicks  corr {c:.4f}  chi2/dof {chi2(clicks, ref):.2f}  fringe err rms {100*r:.1f}% max {100*m:.1f}%  "
          f"per-fringe {np.round(100*e).astype(int).tolist()} {extra}", flush=True)
    return c, r, m

"""Search for a single local rule per ball that produces a real photon.

Transverse (x) dynamics of a paraxial photon: z plays the role of time,
the photon behaves like a particle of 'mass' k = 2π/λ with ħ = 1:
    i ∂ψ/∂z = −(1/2k) ∂²ψ/∂x²
Balls move along z at c; the rule decides how their x changes.
"""
import numpy as np
from scipy.special import ndtri
import sys, time

LAM = 1.0
K = 2 * np.pi / LAM


# ── exact references ─────────────────────────────────────────────────────────
def gaussian_width(w0, z):
    zr = np.pi * w0**2 / LAM
    return w0 * np.sqrt(1 + (z / zr) ** 2)


def fresnel_double_slit(xs, w0, d, a, L, m=4000):
    """Paraxial Fresnel integral over the two slits, |U|² normalized."""
    s = np.concatenate([np.linspace(c - a / 2, c + a / 2, m) for c in (-d / 2, d / 2)])
    ds = a / m
    inc = np.exp(-(s**2) / w0**2)
    U = (inc[None, :] * np.exp(1j * K * (xs[:, None] - s[None, :]) ** 2 / (2 * L))).sum(1) * ds
    I = np.abs(U) ** 2
    return I / I.sum()


# ── rules ────────────────────────────────────────────────────────────────────
def accel_miw(x, stencil=1):
    """Many-Interacting-Worlds force from neighbour gaps (Hall, Deckert,
    Wiseman 2014), U = ħ²/8m Σ (1/d_n − 1/d_{n−1})².  With stencil s > 1 the
    local density is read over s gaps (smoother, less noisy)."""
    N = len(x)
    Kc = 1.0 / (8 * K)
    if stencil == 1:
        a = np.zeros(N + 1)
        a[1:N] = 1.0 / np.diff(x)
        g = np.zeros(N + 1)
        g[1:N] = 2 * Kc * (2 * a[1:N] - a[0:N - 1] - a[2:N + 1])
        dU = g[1:] * a[1:] ** 2 - g[:-1] * a[:-1] ** 2
        return -dU / K
    # smoothed: density over s gaps, quantum potential Q = −(1/2m) √ρ''/√ρ
    s = stencil
    xp = np.concatenate([x[:1] - (x[s:0:-1] - x[0]), x, x[-1] + (x[-1] - x[-2:-s - 2:-1])])
    # local spacing from s-neighbour span
    span = (xp[2 * s:] - xp[:-2 * s]) / (2 * s)
    rho = 1.0 / (N * span)
    r = np.sqrt(rho)
    # derivatives on the irregular grid
    d1 = np.gradient(r, x)
    d2 = np.gradient(d1, x)
    Q = -(1.0 / (2 * K)) * d2 / r
    Q[:s] = Q[s]; Q[-s:] = Q[-s - 1]
    F = -np.gradient(Q, x)
    return F / K


def evolve(x, v, L, rule, dz, stencil=1, max_steps=400000):
    """Velocity-Verlet with an adaptive step: never let two neighbours close
    their gap in one step, never let the fastest change exceed a fraction."""
    x = x.copy(); v = v.copy()
    if rule == "rays":
        return x + v * L
    acc = accel_miw(x, stencil)
    z = 0.0; crossings = 0; steps = 0
    while z < L and steps < max_steps:
        gap = np.diff(x)
        rel_v = np.abs(np.diff(v)) + 1e-12
        rel_a = np.abs(np.diff(acc)) + 1e-12
        h = min(dz, L - z,
                0.2 * np.min(gap / rel_v),
                0.2 * np.sqrt(np.min(gap / rel_a)))
        v += 0.5 * h * acc
        x += h * v
        if np.any(np.diff(x) <= 0):
            crossings += 1
            order = np.argsort(x); x = x[order]; v = v[order]
        acc = accel_miw(x, stencil)
        v += 0.5 * h * acc
        z += h; steps += 1
    return x, (crossings, steps)


# ── experiments ──────────────────────────────────────────────────────────────
def gaussian_balls(N, w0, rng=None):
    sig = w0 / 2
    q = (np.arange(N) + 0.5) / N
    x = sig * ndtri(q)
    if rng is None:
        return x, np.zeros(N)
    # rays: random angles from the Wigner function (initial condition)
    rng.shuffle(x)
    th = rng.normal(0, LAM / (2 * np.pi * w0), N)
    return np.sort(x), th


def free_beam(rule, N=2000, w0=3.0, stencil=1, dz=0.02):
    zr = np.pi * w0**2 / LAM
    out = []
    for z in (zr, 2 * zr, 5 * zr):
        if rule == "rays":
            rng = np.random.default_rng(1)
            sig = w0 / 2
            x0 = rng.normal(0, sig, N); th = rng.normal(0, LAM / (2 * np.pi * w0), N)
            xz = x0 + th * z
        else:
            x0, v0 = gaussian_balls(N, w0)
            xz, _ = evolve(x0, v0, z, rule, dz, stencil)
        w = 2 * xz.std()
        out.append((z, w, gaussian_width(w0, z)))
    return out


def double_slit(rule, N=20000, w0=7.0, d=8.0, a=2.0, L=160.0, stencil=1, dz=0.01, bins=240, half=60):
    edges = np.linspace(-half, half, bins + 1)
    xs = 0.5 * (edges[1:] + edges[:-1])
    ref = fresnel_double_slit(xs, w0, d, a, L)
    if rule == "rays":
        rng = np.random.default_rng(2)
        x0 = rng.normal(0, w0 / 2, N); th = rng.normal(0, LAM / (2 * np.pi * w0), N)
        keep = (np.abs(x0 - d / 2) <= a / 2) | (np.abs(x0 + d / 2) <= a / 2)
        xl = x0[keep] + th[keep] * L
        cr = 0
    else:
        x0, v0 = gaussian_balls(N, w0)
        keep = (np.abs(x0 - d / 2) <= a / 2) | (np.abs(x0 + d / 2) <= a / 2)
        xl, cr = evolve(x0[keep], v0[keep], L, rule, dz, stencil)
    h, _ = np.histogram(xl, edges)
    h = h / max(h.sum(), 1)
    corr = np.corrcoef(h, ref)[0, 1]
    return corr, keep.sum(), cr, h, ref, xs


if __name__ == "__main__":
    rules = [("R1 rays", "rays", 1), ("R2 MIW nearest", "miw", 1), ("R3 MIW smooth s=8", "miw", 8)]
    for name, rule, s in rules:
        t = time.time()
        fb = free_beam(rule, stencil=s)
        txt = "  ".join(f"z={z:5.1f}: w={w:5.2f} (th {th:5.2f})" for z, w, th in fb)
        print(f"{name:22s} free beam  {txt}")
        corr, n, cr, *_ = double_slit(rule, stencil=s)
        print(f"{name:22s} double slit: corr with Fresnel {corr:.3f}  ({n} balls through, reorders/steps {cr})  [{time.time()-t:.1f}s]")
        sys.stdout.flush()

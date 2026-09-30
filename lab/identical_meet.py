"""Where does a photon of identical balls break? Two SMOOTH half-photons that meet and interfere.

  crossing : two Gaussian beams (w = 2) at x = ±6 tilted towards each other (angle ±0.1) -> fringes of period 5
  soft     : 'soft double slit' = two Gaussian slits (w = 1) at x = ±4, no tilt; they diffract and overlap
  phase    : like soft, but the right half starts with an extra phase π (central DARK fringe)

Every ball carries ψ (rule C, log carrier, consistent Crank-Nicolson). Start: quantiles (or random) of |ψ|².
Reference: exact paraxial evolution (FFT). Scores: correlation of the ball histogram with |ψ(L)|² (bins of 1),
fraction of balls in the dark zones, ball order swaps, and the maximum |ψ| (blow-up detector).

  python lab/identical_meet.py CASE N [seed|q] [m=5 P=4 dz=0.02 vmax=1]
"""
import sys, time, numpy as np
import identical_core as ic
from identical_core import K, run_psi, ndtri

CASES = {
    "crossing": dict(c=6.0, w=2.0, u=0.1, ph=0.0, L=120.0),
    "soft": dict(c=4.0, w=1.0, u=0.0, ph=0.0, L=60.0),
    "phase": dict(c=4.0, w=1.0, u=0.0, ph=np.pi, L=60.0),
}


def psi0(x, c, w, u, ph):
    return (np.exp(-(x + c) ** 2 / w ** 2 + 1j * K * u * x)
            + np.exp(-(x - c) ** 2 / w ** 2 - 1j * K * u * x + 1j * ph))


def soft_exact(x, z, ph=0.0, c=4.0, w=1.0):
    """analytic paraxial field of the 'soft' double slit (two Gaussians, waist w, at ±c)."""
    G = lambda y: np.sqrt(w * w / (w * w + 2j * z / K)) * np.exp(-y ** 2 / (w * w + 2j * z / K))
    return G(x + c) + np.exp(1j * ph) * G(x - c)


def meeting_error(N, zc=(4.0, 6.0, 8.0), m=3, P=2, dz=0.01, case="soft"):
    """Follow the balls up to the moment the two halves meet and compare each ball with the exact field:
    returns, at each zc, the worst relative error of the carried ψ, the worst velocity error, and where the
    innermost ball sits (the balls closest to the centre are the ones that meet first)."""
    ic.DEMOD["vmax"] = 0
    p = CASES[case]
    x, ps, _ = start(case, N)
    x_in = np.abs(x).min()
    t = time.time(); out = []
    z = 0.0; v = np.zeros(N); ps = ps[:, None]
    for zt in zc:
        sn, pp, info = run_psi(x, ps, (zt,), m, P, lambda x, z: dz, z0=z, scheme="cn2", v0=v, tmax=500)
        if sn[0] is None:
            out.append((zt, np.nan, np.nan)); break
        x, ps, z = sn[0], pp[0], zt
        v, *_ = ic.psi_fit(x, ps, v, m, P)
        U = soft_exact(x, z, p["ph"]); h = 1e-5
        Ux = (soft_exact(x + h, z, p["ph"]) - soft_exact(x - h, z, p["ph"])) / (2 * h)
        vex = np.imag(np.conj(U) * Ux) / (K * np.abs(U) ** 2)
        out.append((zt, np.max(np.abs(ps[:, 0] - U) / np.abs(U)), np.max(np.abs(v - vex))))
    return x_in, out, time.time() - t


def reference(case, n=2 ** 15, X=600.0):
    p = CASES[case]
    xg = (np.arange(n) - n // 2) * (X / n)
    kx = 2 * np.pi * np.fft.fftfreq(n, X / n)
    U = np.fft.ifft(np.fft.fft(psi0(xg, p["c"], p["w"], p["u"], p["ph"])) * np.exp(-1j * kx ** 2 * p["L"] / (2 * K)))
    return xg, np.abs(U) ** 2


def start(case, N, seed=None, floor=0.0, span=5.0):
    """floor = 0: balls distributed as |ψ|² (all weights equal: the strict requirement).
    floor = f > 0: a fraction f of the balls is spread uniformly over |x| < span (they carry the wave in the
    dark); every ball then carries a weight μ ∝ |ψ|²/(sampling density), constant along its path."""
    p = CASES[case]
    xg = np.linspace(-20, 20, 40001)
    I = np.abs(psi0(xg, p["c"], p["w"], p["u"], p["ph"])) ** 2
    I /= np.trapezoid(I, xg)
    q = (1 - floor) * I + floor * (np.abs(xg) < span) / (2 * span)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (q[1:] + q[:-1]))]); cdf /= cdf[-1]
    u = (np.arange(N) + 0.5) / N if seed is None else np.sort(np.random.default_rng(seed).random(N))
    x = np.interp(u, cdf, xg)
    mu = np.interp(x, xg, I) / np.interp(x, xg, q); mu /= mu.sum()
    return x, psi0(x, p["c"], p["w"], p["u"], p["ph"]), mu


def run(case, N, seed=None, m=5, P=4, dz=0.02, vmax=1.0, tmax=520, floor=0.0):
    ic.DEMOD["vmax"] = vmax
    p = CASES[case]
    x, ps, mu = start(case, N, seed, floor)
    t = time.time()
    try:
        sn, pp, info = run_psi(x, ps, (p["L"],), m, P, lambda x, z: dz, scheme="cn2", tmax=tmax)
    except Exception as e:                                         # linear algebra on NaNs
        return None, dict(fail=True, z=np.nan, swaps=-1, steps=-1, time=time.time() - t, err=str(e)), None, None
    if sn[0] is None:
        return None, info, None, None
    # the weight travels with its ball: find each ball's final rank (balls are re-sorted when they swap)
    mu_end = mu[info["lab"]]
    return sn[0], info, (pp[0] if pp else None), mu_end


def score(case, xl, mu=None):
    xg, I = reference(case)
    p = CASES[case]
    half = 30.0 if case != "crossing" else 20.0
    edges = np.arange(-half, half + 1e-9, 1.0)
    ref = np.array([I[(xg >= a) & (xg < b)].sum() for a, b in zip(edges[:-1], edges[1:])]); ref /= ref.sum()
    h, _ = np.histogram(xl, edges, weights=mu)
    corr = np.corrcoef(h, ref)[0, 1]
    dark = ref < 0.1 * ref.max()
    return corr, h[dark].sum() / max(h.sum(), 1), ref[dark].sum()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "meeting":
        for N in [int(a) for a in sys.argv[2].split(",")]:
            x_in, out, dt = meeting_error(N)
            txt = "  ".join(f"z={z:.0f}: ψ err {100*e:.2g}%  v err {dv:.2g}" for z, e, dv in out)
            print(f"meeting (soft double slit, exact start, rule C m=3 P=2) N={N:6d}: innermost ball at |x|={x_in:.3f}  {txt}  [{dt:.0f}s]", flush=True)
        sys.exit()
    case = sys.argv[1] if len(sys.argv) > 1 else "soft"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    seed = None if len(sys.argv) <= 3 or sys.argv[3] == "q" else int(sys.argv[3])
    kw = {k: float(v) for k, v in (a.split("=") for a in sys.argv[4:])}
    xl, info, pl, mu = run(case, N, seed, int(kw.get("m", 5)), int(kw.get("P", 4)), kw.get("dz", 0.02), kw.get("vmax", 1.0),
                           floor=kw.get("floor", 0.0))
    name = f"{case} N={N} start={'quantile' if seed is None else seed} {kw}"
    if xl is None or info.get("fail") or not np.isfinite(xl).all():
        print(f"{name}: FAILED at z={info['z']:.2f} swaps {info['swaps']} [{info['time']:.0f}s] {info.get('err','')}", flush=True)
    else:
        c, df, dth = score(case, xl, mu if kw.get("floor", 0) > 0 else None)
        print(f"{name}: corr {c:.4f}  balls in dark zones {100*df:.2f}% (theory {100*dth:.2f}%)  swaps {info['swaps']}  "
              f"max|ψ| {np.abs(pl).max():.2f}  steps {info['steps']} [{info['time']:.0f}s]", flush=True)

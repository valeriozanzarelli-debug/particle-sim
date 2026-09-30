"""TASK 2 — double slit (L = 160, d = 8, a = 2, w0 = 7, λ = 1) with ONE photon made of N identical balls.

The balls start just behind the slits (flux-based start of the previous study: positions distributed as
the exact Rayleigh-Sommerfeld flux j_z(x, z0); each ball receives the exact local field there, i.e. the
slit is the object that emits them). From z0 on, every ball follows the same local rule. The ball
positions on the screen are the photon's pattern; they are scored against
   * the exact RS flux on the screen (the target, as in REGOLA_FOTONE.md),
   * the paraxial evolution of the same start field (FFT) = what a PERFECT paraxial rule would give.

  python lab/identical_slit.py RULE N [seed|q] [key=value ...]
     RULE: A (HDW, Verlet)  Ai (HDW, implicit midpoint)  BH (Hamiltonian SG, Verlet)  P (polar R,S, RK4)  C (ψ carried, CN)  Ce (ψ, RK4)
     keys: z0 (1.0), m, P (degree), cfl, dzmax, grow, tmax, save=1
"""
import sys, os, time, numpy as np
from identical_core import *

HERE = os.path.dirname(os.path.abspath(__file__))


def run_slit(rule, N, seed=None, z0=1.0, m=None, deg=None, cfl=None, dzmax=None, grow=None, tmax=540, log=None):
    x, U, Ux = slit_start(N, z0, seed)
    v = np.imag(np.conj(U) * Ux) / (K * np.abs(U) ** 2)          # paraxial Bohm velocity of the start field
    gmin = lambda x: np.diff(x).min()
    t = time.time()
    if rule == "Ai":
        dm = dzmax or 0.05
        snaps, info = implicit_midpoint(x, v, (L,), acc_hdw, 2, lambda x: dm, z0=z0, tmax=tmax)
    elif rule in ("A", "BH"):
        if rule == "A":
            acc, span = acc_hdw, gmin
        else:
            mm = int(m or 2); acc = lambda x: acc_sgH(x, mm, int(deg or 2)); span = lambda x: mm * gmin(x)
        c = cfl or 0.2
        snaps, info = verlet(x, v, (L,), acc, lambda x: c * K * span(x) ** 2, z0=z0, tmax=tmax, log=log)
    elif rule == "P":
        mm = int(m or 3); c = cfl or 0.2
        S = np.unwrap(np.angle(U))
        snaps, info = run_polar(x, np.abs(U), S, (L,), mm, int(deg or 2), lambda x: c * K * (mm * gmin(x)) ** 2, z0=z0, tmax=tmax)
    else:
        mm = int(m or 3); deg = deg or 2
        if rule == "Ce":
            c = cfl or 0.2
            dzf = lambda x, z: c * K * (mm * gmin(x)) ** 2
            snaps, _, info = run_psi(x, U, (L,), mm, int(deg), dzf, z0=z0, scheme="rk4", tmax=tmax, v0=v, log=log)
        else:
            dm = dzmax or 0.5; gr = grow or 0.02
            dzf = lambda x, z: min(dm, gr * z)
            snaps, _, info = run_psi(x, U, (L,), mm, int(deg), dzf, z0=z0, scheme="cn2", tmax=tmax, v0=v, log=log)
    info["time"] = time.time() - t
    return snaps[0], info


def dark_fraction(xl, ref, thr=0.1):
    """Fraction of balls landing in dark regions (reference < thr·max, |x| < 45) vs the same for theory."""
    dark = (ref < thr * ref[np.abs(XS) < 45].max()) & (np.abs(XS) < 45)
    h, _ = np.histogram(xl, EDGES)
    return h[dark].sum() / max(len(xl), 1), ref[dark].sum()


def evaluate(xl, z0, N, seed, info, name):
    ref, par = slit_refs(z0)
    c, r, mx = score(xl, ref)
    cp, rp, mp = score(xl, par)
    df, dth = dark_fraction(xl, ref)
    q = quantile_sampler(par, N)
    cq, rq, mq = score(q, ref)
    line = (f"{name}: corr {c:.4f}  fringe err rms {100*r:.1f}% max {100*mx:.1f}%   | vs paraxial-FFT corr {cp:.4f} "
            f"rms {100*rp:.1f}%  | balls in dark zones {100*df:.2f}% (theory {100*dth:.2f}%)  swaps {info['swaps']} "
            f"steps {info['steps']} [{info['time']:.0f}s]")
    base = f"   perfect paraxial photon of {N} balls (quantiles of the FFT pattern): corr {cq:.4f} rms {100*rq:.1f}% max {100*mq:.1f}%"
    if seed is not None:
        b = ideal_baseline(ref, N, 100)
        base += f";  ideal random sampler {N}: corr {b[0]:.4f}±{b[1]:.4f} rms {100*b[2]:.1f}%"
    return (c, r, mx, cp, df), line, base


if __name__ == "__main__":
    rule = sys.argv[1] if len(sys.argv) > 1 else "C"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    seed = None if len(sys.argv) <= 3 or sys.argv[3] == "q" else int(sys.argv[3])
    kw = {k: float(v) for k, v in (a.split("=") for a in sys.argv[4:])}
    z0 = kw.get("z0", 1.0)
    ref, par = slit_refs(z0)

    # paraxial ceiling: FFT pattern scored as an (infinitely many balls) histogram against the exact flux
    cc = np.corrcoef(par.reshape(-1, 8).sum(1), ref.reshape(-1, 8).sum(1))[0, 1]
    print(f"paraxial ceiling (FFT from z0={z0} vs exact RS flux): corr {cc:.4f}", flush=True)
    xl, info = run_slit(rule, N, seed, z0, kw.get("m"), kw.get("P"), kw.get("cfl"), kw.get("dzmax"), kw.get("grow"),
                        kw.get("tmax", 540), int(kw["log"]) if "log" in kw else None)
    name = f"rule {rule} N={N} start={'quantile' if seed is None else 'seed %d' % seed} {kw}"
    if xl is None or info.get("fail"):
        print(f"{name}: FAILED at z={info['z']:.2f} (steps {info['steps']}, swaps {info['swaps']}) [{info['time']:.0f}s]", flush=True)
    else:
        _, line, base = evaluate(xl, z0, N, seed, info, name)
        print(line, flush=True); print(base, flush=True)
        if kw.get("save"):
            os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
            np.save(os.path.join(HERE, "out", f"identical_slit_{rule}_{N}_{'q' if seed is None else seed}.npy"), xl)

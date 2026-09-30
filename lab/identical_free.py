"""TASK 1 — free Gaussian beam (w0 = 3, λ = 1) with identical balls.

Width = 2·std(x) of the N balls of ONE photon at z = zR, 2zR, 5zR vs w0·sqrt(1+(z/zR)²).

  python lab/identical_free.py RULE N [seed|q] [key=value ...]
     RULE: A  (HDW, explicit Verlet)     Ai (HDW, implicit midpoint)
           B  (rank-space SG density)    Bi (same, implicit midpoint)
           BH (Hamiltonian SG = smoothed HDW, P=2)  BHi (implicit midpoint)
           P  (polar R,S carried, RK4)
           C  (ψ carried, Crank-Nicolson) Ce (ψ carried, explicit RK4)
           H  (ψ carried, unitary nearest-neighbour links + move-and-reread; tried, worse)
     seed: 'q' (default) = regular quantile start, integer = random start from |ψ|² with that seed
     keys: m, P(deg), cfl, dzmax, grow, tmax
"""
import sys, time, numpy as np
from identical_core import *


def parse(argv):
    rule = argv[1] if len(argv) > 1 else "C"
    N = int(argv[2]) if len(argv) > 2 else 500
    seed = None if len(argv) <= 3 or argv[3] == "q" else int(argv[3])
    kw = dict(a.split("=") for a in argv[4:])
    return rule, N, seed, {k: float(v) for k, v in kw.items()}


def run(rule, N, seed=None, m=None, deg=None, cfl=None, dzmax=None, grow=None, tmax=540, log=None):
    x, psi = free_start(N, seed)
    import identical_core as _ic
    _ic.DEMOD["theta"] = THETA; _ic.DEMOD["vsmooth"] = VSMOOTH
    v = np.zeros(N)
    gmin = lambda x: np.diff(x).min()
    if rule in ("A", "Ai"):
        acc = acc_hdw; band = 2; span = lambda x: gmin(x)
    elif rule in ("B", "Bi"):
        mm = int(m or 4); acc = lambda x: acc_sg(x, mm, 4); band = 2 * mm
        span = lambda x: mm * gmin(x)
    elif rule in ("BH", "BHi"):
        mm = int(m or 2); acc = lambda x: acc_sgH(x, mm, int(deg or 2)); band = 2 * mm
        span = lambda x: mm * gmin(x)
    if rule in ("A", "B", "BH"):
        c = cfl or 0.2
        snaps, info = verlet(x, v, ZF, acc, lambda x: c * K * span(x) ** 2, tmax=tmax, log=log)
    elif rule in ("Ai", "Bi", "BHi"):
        dm = dzmax or 0.2
        snaps, info = implicit_midpoint(x, v, ZF, acc, band, lambda x: dm, tmax=tmax)
    elif rule == "P":
        mm = int(m or 3); c = cfl or 0.2
        snaps, info = run_polar(x, np.abs(psi), np.zeros(N), ZF, mm, int(deg or 2), lambda x: c * K * (mm * gmin(x)) ** 2, tmax=tmax)
    elif rule == "H":
        mm = int(m or 3); dm = dzmax or 0.1
        snaps, _, info = run_psi(x, psi, ZF, mm, int(deg or 2), lambda x, z: dm, scheme="herm", tmax=tmax)
    elif rule in ("C", "Ce"):
        mm = int(m or 3); deg = deg or 2
        if rule == "Ce":
            c = cfl or 0.2
            dzf = lambda x, z: c * K * (mm * gmin(x)) ** 2
            snaps, _, info = run_psi(x, psi, ZF, mm, int(deg), dzf, scheme="rk4", tmax=tmax, log=log)
        else:
            dm = dzmax or 0.5; gr = grow or 0.02
            dzf = lambda x, z: min(dm, gr * (z + ZR))
            snaps, _, info = run_psi(x, psi, ZF, mm, int(deg), dzf, scheme="cn2", tmax=tmax, log=log)
    res = []
    s0 = 2 * x.std() / W0F            # the start sample's own width (random starts: sampling noise ~ 1/sqrt(2N))
    info["w_start_ratio"] = s0
    for z, X in zip(ZF, snaps):
        if X is None or not np.all(np.isfinite(X)):
            res.append((z, np.nan, w_theory(z), np.nan, np.nan))
        else:
            w = 2 * X.std()
            res.append((z, w, w_theory(z), w / w_theory(z) - 1, w / (s0 * w_theory(z)) - 1))
    return res, info


THETA = 0.5
VSMOOTH = False

if __name__ == "__main__":
    rule, N, seed, kw = parse(sys.argv)
    THETA = kw.get("theta", 0.5); VSMOOTH = bool(kw.get("vs", 0))
    t = time.time()
    res, info = run(rule, N, seed, m=kw.get("m"), deg=kw.get("P"), cfl=kw.get("cfl"), dzmax=kw.get("dzmax"),
                    grow=kw.get("grow"), tmax=kw.get("tmax", 540), log=int(kw["log"]) if "log" in kw else None)
    txt = "  ".join(f"z={z/ZR:.0f}zR: w={w:.3f} (th {wt:.3f}, err {100*e:+.2f}%, vs own start {100*e2:+.2f}%)" for z, w, wt, e, e2 in res)
    print(f"rule {rule} N={N} start={'quantile' if seed is None else 'random seed %d' % seed} {kw}: {txt}  "
          f"start width/w0 {info['w_start_ratio']:.4f} steps {info['steps']} swaps {info['swaps']} {'FAIL(z=%.2f)' % info['z'] if info['fail'] else ''} [{time.time()-t:.1f}s]",
          flush=True)

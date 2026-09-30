"""TASK 2 — how good is the NEIGHBOUR estimate of the field that guides the photon ball?

Variants (all: photon enters with |U_inc|^2, local edge rule, dx/dz = j_x/j_z, RK4 adaptive).
Random-ball variants A-C run WITHOUT the edge rule: with a fixed window h=1 their near field
is so underestimated (the fan is narrower than the window) that the rule would kill ~every photon.
  A  random balls (old pilot / current Rust emission), Gaussian window h=1, plane-wave extrapolation (E1)
  B  random balls, Gaussian window h=1, exact wavelet from the ball's origin (E2)
  D  lattice balls, hat window (2 balls per origin), E1
  E  lattice balls, hat window, E2, ds=0.01, delta=0.05
  F  as E without the edge rule                      (shows why the edge rule is needed)
  G  as E, coarse lattice ds=0.04, delta=0.2         (fewer balls)
  H  as E, dth=0.1 (only 31 directions per origin)  (few balls; E2 does not care about dth)
  I  as E, ds=0.02, delta=0.1                        (half the origins, 2x faster) <- DEFAULT of photon_rule.py
  J  as I, first step behind the slit along the outgoing current (slightly worse: see doc)

Run:  python lab/neighbour_study.py            (variant I only = final rule, ~2 min)
      python lab/neighbour_study.py A D E      (any subset; one variant per command keeps each < 5 min:
                                                A,B ~4-5 min, D,E,F,H ~3-4 min, G ~1 min, I ~2 min)
"""
import sys, os, numpy as np
from multiprocessing import Pool
from optics import *
from photon_rule import double_slit, Fan, DEFAULT

VARIANTS = {
    "A": dict(mode="mc", n_mc=20000, h=1.0, est="E1", nph=2000, edge=False, rk=2, dzmin=0.05),
    "B": dict(mode="mc", n_mc=20000, h=1.0, est="E2", nph=1400, edge=False, rk=2, dzmin=0.05),
    "D": dict(est="E1", ds=0.01, z0=0.05),
    "E": dict(est="E2", ds=0.01, z0=0.05),
    "F": dict(est="E2", ds=0.01, z0=0.05, edge=False),
    "G": dict(est="E2", ds=0.04, z0=0.2),
    "H": dict(est="E2", ds=0.01, z0=0.05, dth=0.1),
    "I": dict(est="E2", ds=0.02, z0=0.1, first="straight"),
    "J": dict(est="E2", ds=0.02, z0=0.1, first="implicit"),
}


def n_balls(v):
    c = dict(DEFAULT); c.update(v)
    if c["mode"] == "mc":
        return int(c["n_mc"] * 0.245)
    return int(round(2 * A / c["ds"])) * int(round(2 * c["TH"] / c["dth"]))


if __name__ == "__main__":
    which = sys.argv[1:] or ["I"]
    ref = screen_reference("flux")
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out"), exist_ok=True)
    for n in (3000, 6000):
        b = ideal_baseline(ref, n, 200)
        print(f"ideal sampler {n} clicks: corr {b[0]:.4f}+-{b[1]:.4f}  fringe rms {100*b[2]:.1f}%  max {100*b[3]:.1f}%", flush=True)
    with Pool(8) as pool:
        for k in which:
            v = VARIANTS[k]
            xl, ab, dt = double_slit(v, pool)
            np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", f"neigh_{k}.npy"), xl)
            report(f"{k} {v} balls/photon~{n_balls(v)}", xl, ref, f"absorbed {100*ab:.1f}% [{dt:.0f}s]")

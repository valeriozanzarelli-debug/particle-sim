"""pre_p1_posthoc.py - POST-HOC additions to P1 (declared as such in docs/PRE_ESPERIMENTI.md).
 (1) P1d redo: the pre-registered 'lam = 1' run used a scalar lam = 1 also on the box WALL links (w = 0), where
     it lets surfers jump across the wrap-around link (x = n-1 -> 0): a teleport bug. Redo with lam = 'w'
     (1 inside, 0 at walls) in the box, and on a TORUS (no walls), where the sheets cannot be glued.
 (2) P1e refinement: is the M = 1 'fake mixing' a lattice artefact that dies in the continuum limit (Vink 1993:
     Bell's lattice jump process -> deterministic Bohm as a -> 0)? theta = m a in {0.5, 0.25, 0.125} with
     n = 24, 48, 96 (same box in physical units), omega0 and T scaled with theta; tau reported in units of 1/m.
Writes lab/out/pre_p1_posthoc.json (~2 min).
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from pre_p1p5 import run

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

if __name__ == "__main__":
    t0 = time.time()
    res = {}
    for name, kw in (("box_lamw", dict(periodic=False, lam="w")), ("torus_lamw", dict(periodic=True, lam="w")),
                     ("torus_lam0", dict(periodic=True, lam=0.0))):
        r = run(24, 16, "sheet0", kw["lam"], 400, sheets=True, periodic=kw["periodic"])
        res["P1d_" + name] = dict(Hs0=float(r["Hs"][0]), HsT=float(r["Hs"][-1]), ratio=float(r["Hs"][-1] / r["Hs"][0]),
                                  Fs=r["Fs"], H_cells_T=float(r["H"][-1]), Hs_curve=r["Hs"][::8].tolist())
        print(name, {k: v for k, v in res["P1d_" + name].items() if k != "Hs_curve"}, flush=True)
    ref = {}
    for th, n in ((0.5, 24), (0.25, 48), (0.125, 96)):
        f = 0.5 / th
        for M in (1, 16):
            for lam in ("w", 0.0):
                T = int(400 * f)
                r = run(n, M, "uniform", lam, T, cell=int(4 * f), omega0=0.9 * th / 0.5, every=int(5 * f))
                key = f"theta{th}_M{M}_lam{lam}"
                ref[key] = dict(n=n, tau_ticks=r["tau"], tau_phys=r["tau"] * th, ratio=r["ratio"],
                                H0=float(r["H"][0]))
                print(key, ref[key], f"{time.time()-t0:.0f}s", flush=True)
    res["P1e_refinement"] = ref
    res["time_s"] = time.time() - t0
    with open(os.path.join(OUT, "pre_p1_posthoc.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print("time", res["time_s"])

"""nucleo1s_posthoc.py - AFTER the pre-registered runs (not pre-registered, declared in docs/NUCLEO_1S.md):
(1) automatic zoom driven by ALL ten levels of the orbitali scene (Z = 60), budget 0.5 %: what the Laboratorio Atomo
    would need for every level within ~1 %; (2) the Z = 20 L1 failure: R(a) along the whole a sequence."""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_run import *

part = sys.argv[1]
if part == "allzoom":
    Za = 60 * ALPHA
    rmax = rmax_of(4, Za)
    z = AutoZoom(0.5, rmax, STATES, Za, 0.005, log=log).run()
    lv = {NAMES[s]: float(z.E[i] / dirac_binding(s[0], s[1], Za) - 1) for i, s in enumerate(STATES)}
    save("posthoc_allzoom", dict(budget=0.005, cells=len(z.edges) - 1, uniform_cells_a05=int(np.ceil(rmax / 0.5)),
                                 levels=lv, est=z.history[-1]["est_err"], edges=z.edges.tolist(),
                                 iters=len(z.history)))
    log("all-level zoom:", len(z.edges) - 1, lv)
elif part == "geo":
    # geometric idea: cells proportional to r (h = eps r, between h_min and h_max): a scale-free lattice for a
    # scale-free singularity. Prediction: no contact term (error independent of h_min), only eps^2 and h_max^2.
    Za = 60 * ALPHA
    rmax = rmax_of(1, Za)
    ref = dirac_binding(1, -1, Za)

    def geo(eps, hmin, hmax):
        e = [0.0]
        while e[-1] < rmax:
            e.append(e[-1] + min(hmax, max(hmin, eps * (e[-1] + 0.5 * hmin))))
        return Mesh(np.array(e))
    out = []
    for hmax in (0.5, 4.0):
        for eps in (0.4, 0.2, 0.1, 0.05, 0.025, 0.0125):
            for hmin in (1e-3, 1e-5, 1e-7):
                m = geo(eps, hmin, hmax)
                b, _ = level(m, 1, -1, Za, guess_rel=1.05)
                out.append(dict(hmax=hmax, eps=eps, hmin=hmin, cells=m.N, err=b / ref - 1))
                log(out[-1])
    save("posthoc_geo", out)
elif part == "R":
    r = load()["law"]
    out = {}
    for Z in [z for z in r if z.isdigit()]:
        pz = PRE["predictions_law"]["values"][Z]
        d = r[Z]["1,-1"]
        a = np.array(d["a"]); e = np.array(d["err"])
        R = (e - pz["K"] * a**pz["two_gamma"]) / (pz["B"] * a**2)
        out[Z] = dict(a=a.tolist(), R=R.tolist())
        log(Z, np.round(R, 5))
    save("posthoc_R", out)

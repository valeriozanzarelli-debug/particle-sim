"""supercrit_run_g1.py - exam G1/G3: resonance energy and width of the dived 1s (shell R = 0.031) from the
phase scan of the stationary rule (Bloch phase at r_m = 25 jumps by 2 pi). Seed energy from the quasi-bound
box state (closed at the middle of the barrier). Checkpointed: lab/out/supercrit_g1.json.
Usage: python supercrit_run_g1.py [Z ...]"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supercrit_core import SWalk, V_shell, ALPHA, resonance_scan2, quasi_bound
from supercrit_prereg import P, GODUNOV

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "supercrit_g1.json")
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
cfg = P["configs"]["G1"]
Zlist = [int(x) for x in sys.argv[1:]] or cfg["Z"]
alist = cfg["a"]
t0 = time.time()
rm = cfg["phase_site_r"]
for Z in Zlist:
    Za = Z * ALPHA
    d = res.setdefault(str(Z), {})
    E0 = -1.0 - 0.052 * (Z - 177.46)          # rough seed (slope of the shell track near Z_cr, dev)
    for a in alist:
        key = f"{a:.6g}"
        if key in d or time.time() - t0 > 450:
            continue
        Vf = lambda r: V_shell(r, Za, 0.031)
        om, chi, rw = quasi_bound(lambda rmax: SWalk(-1, a, rmax, Vf), E0)
        w = SWalk(-1, a, rm + 2.0, Vf)
        Er, G, rms, (E, ph), GF = resonance_scan2(w, om, int(rm / a), int(rw / a))
        d[key] = dict(a=a, E_box=om, r_w=rw, E_res=Er, Gamma=G, Gamma_peak=GF, fit_rms=rms)
        E0 = Er
        print(Z, a, "box", round(om, 6), "res", Er, "Gamma", G, "GF", GF, "rms", rms, f"{time.time() - t0:.0f}s", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)
    runs = sorted([v for k, v in d.items() if k[0].isdigit()], key=lambda x: x["a"])
    if len(runs) == len(alist):
        A = np.array([x["a"] for x in runs])
        ex = {}
        for q in ("E_res", "Gamma"):
            y = np.array([x[q] for x in runs])
            ex[q + "_quad3"] = float(np.polyval(np.polyfit(A[:3], y[:3], 2), 0))
            ex[q + "_lin2"] = float(np.polyval(np.polyfit(A[:2], y[:2], 1), 0))
        if Z in GODUNOV:
            Eg, Gg = GODUNOV[Z]
            ex["dE_rel"] = (ex["E_res_quad3"] - Eg) / abs(Eg)
            ex["dG_rel"] = (ex["Gamma_quad3"] - Gg) / Gg
        d["extrap"] = ex
        print(Z, ex, flush=True)
        json.dump(res, open(OUT, "w"), indent=1)

"""box_formation.py - A4: the 'formation length' of a pair. 1+1 D zig-zag rule (lab/redo_r1.py), Sauter step
V(x) = (V0/2)(1 + tanh(x/ell)) with FIXED peak field E = V0/(2 ell) = 2 th^2, ell = 2..160 sites.
Yield Gamma(ell) = int dw/2pi N(w) over the Klein zone (exact lattice transfer matrix, lattice_N of redo_r1), vs the
local-constant-field estimate Gamma_LCFA = int dx E(x)/2pi exp(-pi m^2/E(x)). Formation length ell_f = 2m/E.
Also (not graded) the same with m -> m_eff of the strip d = 16 (lowest transverse mode treated as a 1+1 D particle).
Writes lab/out/box_formation.json. (~3 min)"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys, time
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from redo_r1 import lattice_N, sauter_formula

PI = np.pi
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
TH = 0.1
E = 2 * TH**2
ELLS = [2, 3, 5, 7, 10, 14, 20, 28, 40, 56, 80, 113, 160]


def lcfa(m, V0, ell):
    f = lambda x: (V0 / (2 * ell)) / np.cosh(x / ell) ** 2
    g = lambda x: f(x) / (2 * PI) * np.exp(-PI * m * m / f(x)) if f(x) > 1e-300 else 0.0
    return 2 * quad(g, 0, 30 * ell, limit=400)[0]


def yield_(th, ell, nw=120):
    V0 = 2 * E * ell
    lo, hi = max(th, V0 - PI + th), min(PI - th, V0 - th)
    if hi <= lo:
        return 0.0, 0.0, V0, (lo, hi)
    ws = np.linspace(lo, hi, nw + 2)[1:-1]
    N = np.array([lattice_N(w, th, V0, ell) for w in ws])
    N = np.nan_to_num(N)
    G = float(np.trapezoid(N, ws) / (2 * PI))
    # continuum Sauter (closed form, checked against the ODE in R1c), only meaningful while the window is not clipped
    Gc = float(np.trapezoid([sauter_formula(w, th, V0, ell) for w in ws], ws) / (2 * PI)) if V0 < PI - th else float("nan")
    return G, Gc, V0, (lo, hi)


if __name__ == "__main__":
    t0 = time.time()
    res = dict(th=TH, E=E, ell_f=2 * TH / E, rows={})
    for th, tag in ((TH, "free"), (0.16432048228126003, "strip_d16")):
        rows = []
        for ell in ELLS:
            G, Gc, V0, win = yield_(th, ell)
            Gl = lcfa(th, V0, ell)
            rows.append(dict(ell=ell, V0=V0, window=win, Gamma_lat=G, Gamma_sauter_exact=Gc, Gamma_LCFA=Gl,
                             ratio=G / Gl if Gl > 0 else float("nan"), band_wrap=bool(V0 > PI - th)))
            print(tag, ell, round(V0, 3), G, Gc, Gl, round(G / Gl, 4) if Gl > 0 else None, f"{time.time()-t0:.0f}s", flush=True)
        res["rows"][tag] = rows
    lf = res["ell_f"]
    fr = res["rows"]["free"]
    a = [r for r in fr if r["ell"] >= 8 * lf]
    b = [r for r in fr if r["ell"] <= lf]
    res["A4a_PASS"] = bool(all(abs(r["ratio"] - 1) <= 0.15 for r in a))
    res["A4b_PASS"] = bool(all(r["ratio"] <= 0.5 for r in b))
    res["A4_PASS"] = res["A4a_PASS"] and res["A4b_PASS"]
    res["A4a_rows"] = [(r["ell"], r["ratio"], r["band_wrap"]) for r in a]
    res["A4b_rows"] = [(r["ell"], r["ratio"]) for r in b]
    print("A4a", res["A4a_PASS"], res["A4a_rows"], "A4b", res["A4b_PASS"], res["A4b_rows"])
    res["time_s"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "box_formation.json"), "w"), indent=1, default=float)

"""locale_f3_stagger.py - EXPLORATION after the phase-3 run (NOT pre-registered; a proposal for a next phase).
Finding of the run: the 4 pi-valleys of one copy are massless for the frozen rule AND for the twisted Wilson F2, because a
chirality-mixing term (the mass flip, the Wilson rotation) changes sign between the two ticks at those momenta and
cancels over a period. Idea: a Wilson term whose sign ALTERNATES in time (tick 0: h_u + h_s, tick 1: h_u - h_s).
At k = 0 every Wilson term vanishes, so the low-energy electron should not notice it.
Measures, at theta = 0.1: half-gaps at the 7 doubler valleys, the zone minimum outside |k| < 0.4 (from 0 and from pi),
E(0), dispersion and kinetic mass. Writes lab/out/locale_f3_stagger.json."""
import os, sys, json
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import SplitWalk, sdot, wilson_factors, zigzag_ticks, XCHI, CUBE_AXES

PI = np.pi
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f3_stagger.json")
VALLEYS = {"N_xy": ((.5, .5, 0), 0.0), "N_xz": ((.5, 0, .5), 0.0), "N_yz": ((0, .5, .5), 0.0),
           "X_x": ((.5, 0, 0), PI), "X_y": ((0, .5, 0), PI), "X_z": ((0, 0, .5), PI), "R_xyz": ((.5, .5, .5), PI)}


def walk(theta, hu, hs):
    fac = []
    for it, tk in enumerate(zigzag_ticks(3)):
        for d in tk:
            M = np.zeros((4, 4), complex); M[:2, :2] = -sdot(d, 0.5); M[2:, 2:] = sdot(d, 0.5)
            fac.append(("shift", d, M))
        h = hu + (hs if it == 0 else -hs)
        vecs = [2 * e for e in CUBE_AXES] if it == 0 else [-2 * e for e in CUBE_AXES]
        wl = dict(h=h / 2, vecs=vecs, axis="Y")
        if h != 0:
            fac += wilson_factors(wl)
        fac.append(("coin", np.cos(theta) * np.eye(4) - 1j * np.sin(theta) * XCHI))
        if h != 0:
            fac += wilson_factors(wl)
    return SplitWalk(fac, 4, np.eye(3), dim=3)


def rel(eps, c):
    return np.angle(np.exp(-1j * (eps - c)))


def measure(w, th):
    out = {}
    for n, (kv, c) in VALLEYS.items():
        e, _ = w.phases(PI * np.array(kv, float)[None])
        out[n] = float(np.min(np.abs(rel(e[0], c))) / w.T)
    n = 24
    g = (np.arange(n) + 0.5) / n * PI - PI / 2
    K = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    K = K[np.linalg.norm(K, axis=1) >= 0.4]
    e, _ = w.phases(K)
    out["zone_g0"] = float(np.min(np.abs(e)) / w.T)
    out["zone_gpi"] = float(np.min(PI - np.abs(e)) / w.T)
    El = lambda Q: np.min(np.abs(w.phases(np.atleast_2d(Q))[0]), axis=1) / w.T
    out["E0_rel"] = float(El(np.zeros((1, 3)))[0] / th - 1)
    errs, mk = [], []
    for d in (np.array([1, 0, 0.]), np.array([1, 1, 0]) / 2 ** .5, np.array([1, 1, 1]) / 3 ** .5):
        P = np.linspace(0.1, 1, 10)[:, None] * th * d[None]
        errs.append(float(np.max(np.abs(El(P) / np.hypot(np.linalg.norm(P, axis=1), th) - 1))))
        ps = np.linspace(0.02, 0.1, 9) * th
        E = El(ps[:, None] * d[None]); E0 = El(np.zeros((1, 3)))[0]
        cc = np.linalg.lstsq(np.stack([ps ** 2, ps ** 4], 1), E ** 2 - E0 ** 2, rcond=None)[0]
        mk.append(float(1 / cc[0] - 1))
    out["disp"] = max(errs); out["mkin_rel"] = mk
    out["min_doubler_gap"] = min(out["zone_g0"], out["zone_gpi"], *[out[k] for k in VALLEYS])
    return out


R = {}
th = 0.1
for hu in (0.0, 0.1, 0.6):
    for hs in (0.0, 0.1, 0.3, 0.6):
        r = measure(walk(th, hu, hs), th)
        R[f"{hu}|{hs}"] = r
        print(f"hu={hu} hs={hs}: N {r['N_xy']:.3f} X {r['X_x']:.3f} R {r['R_xyz']:.3f} zone0 {r['zone_g0']:.3f} "
              f"zonepi {r['zone_gpi']:.3f} | E0 {r['E0_rel']:.1e} disp {r['disp']:.2e} mkin {r['mkin_rel'][0]:+.2e}", flush=True)
with open(OUT, "w") as f:
    json.dump(R, f, indent=1)

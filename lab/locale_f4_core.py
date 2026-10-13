"""locale_f4_core.py - phase 4 of REGOLA LOCALE: the zig-zag 3-D Dirac walk with extra LOCAL terms that act differently
on the two ticks of the period (to give mass to the pi-valleys):
  coin of tick t:  exp(-i [(theta + s_t mu_x) X + s_t mu_y Y]),  s_0 = -1, s_1 = +1   ('breathing' coin; mu_x = VUOTO's)
  twisted link Wilson (axis Y) of tick t with strength h_u + s'_t h_s,  s'_0 = +1, s'_1 = -1 (uniform + alternating).
All factors are on-site or nearest-link unitaries: unitary, local (radius 2), charge conserving."""
import numpy as np
from locale_core import SplitWalk, sdot, wilson_factors, zigzag_ticks, XCHI, YCHI, CUBE_AXES

PI = np.pi
VALLEYS = {"N_xy": ((.5, .5, 0), 0.0), "N_xz": ((.5, 0, .5), 0.0), "N_yz": ((0, .5, .5), 0.0),
           "X_x": ((.5, 0, 0), PI), "X_y": ((0, .5, 0), PI), "X_z": ((0, 0, .5), PI), "R_xyz": ((.5, .5, .5), PI)}


def rot(a, b):
    """exp(-i (a X + b Y)) on the chirality index (X, Y anticommute, square to 1)"""
    n = np.hypot(a, b)
    if n == 0:
        return np.eye(4, dtype=complex)
    return np.cos(n) * np.eye(4) - 1j * np.sin(n) * (a * XCHI + b * YCHI) / n


def walk(theta, mu_x=0.0, mu_y=0.0, hu=0.0, hs=0.0):
    fac = []
    for it, tk in enumerate(zigzag_ticks(3)):
        s = -1 if it == 0 else +1
        for d in tk:
            M = np.zeros((4, 4), complex); M[:2, :2] = -sdot(d, 0.5); M[2:, 2:] = sdot(d, 0.5)
            fac.append(("shift", d, M))
        h = hu + (hs if it == 0 else -hs)
        vecs = [2 * e for e in CUBE_AXES] if it == 0 else [-2 * e for e in CUBE_AXES]
        wl = dict(h=h / 2, vecs=vecs, axis="Y")
        if h != 0:
            fac += wilson_factors(wl)
        fac.append(("coin", rot(theta + s * mu_x, s * mu_y)))
        if h != 0:
            fac += wilson_factors(wl)
    return SplitWalk(fac, 4, np.eye(3), dim=3)


def rel(eps, c):
    return np.angle(np.exp(-1j * (eps - c)))


def E_low(w, K):
    return np.min(np.abs(w.phases(np.atleast_2d(K))[0]), axis=1) / w.T


DIRS = [np.array([1, 0, 0.]), np.array([1, 1, 0]) / 2 ** .5, np.array([1, 1, 1]) / 3 ** .5, np.array([1, 2, 3]) / 14 ** .5]


def measure(w, th, ref=None, n=24):
    out = {}
    for name, (kv, c) in VALLEYS.items():
        e, _ = w.phases(PI * np.array(kv, float)[None])
        out[name] = float(np.min(np.abs(rel(e[0], c))) / w.T)
    g = (np.arange(n) + 0.5) / n * PI - PI / 2
    K = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    K = K[np.linalg.norm(K, axis=1) >= 0.4]
    e, _ = w.phases(K)
    out["zone_g0"] = float(np.min(np.abs(e)) / w.T)
    out["zone_gpi"] = float(np.min(PI - np.abs(e)) / w.T)
    out["min_doubler_gap"] = min(out["zone_g0"], out["zone_gpi"], *[out[k] for k in VALLEYS])
    out["E0_rel"] = float(abs(E_low(w, np.zeros((1, 3)))[0] / th - 1))
    mk, disp = [], []
    for d in DIRS:
        ps = np.linspace(0.02, 0.1, 9) * th
        E = E_low(w, ps[:, None] * d[None]); E0 = E_low(w, np.zeros((1, 3)))[0]
        cc = np.linalg.lstsq(np.stack([ps ** 2, ps ** 4], 1), E ** 2 - E0 ** 2, rcond=None)[0]
        mk.append(float(1 / cc[0] - 1))
        P = np.linspace(0.1, 1, 10)[:, None] * th * d[None]
        disp.append(float(np.max(np.abs(E_low(w, P) / np.hypot(np.linalg.norm(P, axis=1), th) - 1))))
    out["mkin_rel"] = mk; out["disp"] = max(disp)
    # drift
    p = 1e-5
    gr = []
    for j in range(3):
        ev = np.zeros(3); ev[j] = p
        a = np.sort(w.phases(ev[None])[0][0]); b = np.sort(w.phases(-ev[None])[0][0])
        gr.append(float(np.max(np.abs(a - b)) / w.T / (2 * p)))
    out["drift"] = max(gr)
    U = w.bloch(K[:50]); out["unitarity"] = float(np.max(np.abs(np.einsum("nba,nbc->nac", U.conj(), U) - np.eye(4))))
    if ref is not None:
        # deviation from the reference walk (the frozen rule) near Gamma: dE(p) and its power law
        ps = np.geomspace(0.1, 2.0, 12) * th
        dev = []
        for d in DIRS:
            dev.append(np.abs(E_low(w, ps[:, None] * d[None]) - E_low(ref, ps[:, None] * d[None])))
        dev = np.max(np.array(dev), axis=0)
        ok = dev > 1e-14
        slope = float(np.polyfit(np.log(ps[ok]), np.log(dev[ok]), 1)[0]) if ok.sum() >= 3 else None
        out["dE_vs_frozen"] = dict(p_over_m=(ps / th).tolist(), dE=dev.tolist(), slope=slope,
                                   rel_at_p_eq_m=float(dev[np.argmin(np.abs(ps - th))] / th))
        out["mkin_vs_frozen"] = [float((1 + a) / (1 + b) - 1) for a, b in zip(mk, ref_mk(ref, th))]
    return out


_cache = {}


def ref_mk(ref, th):
    key = (id(ref), th)
    if key not in _cache:
        mk = []
        for d in DIRS:
            ps = np.linspace(0.02, 0.1, 9) * th
            E = E_low(ref, ps[:, None] * d[None]); E0 = E_low(ref, np.zeros((1, 3)))[0]
            cc = np.linalg.lstsq(np.stack([ps ** 2, ps ** 4], 1), E ** 2 - E0 ** 2, rcond=None)[0]
            mk.append(float(1 / cc[0] - 1))
        _cache[key] = mk
    return _cache[key]

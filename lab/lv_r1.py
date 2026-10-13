"""lv_r1.py - esami E1-E8 della pre-registrazione (controllo numerico della dispersione dalla regola congelata).
Scrive lab/out/lv_r1.json. ~1 minuto, 1 core."""
import json
import os
import time
import itertools
import numpy as np
import lv_core as lv

t0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261013)
res = {}


def verdict(x, thr):
    return "PASS" if x < thr else "FAIL"


# E1 luce, forma chiusa
K = rng.uniform(-np.pi / 2, np.pi / 2, (2000, 3))
omL = lv.omega_numeric(K, 0, 1)                      # (N,3): -w, 0, +w
e1 = np.max(np.abs(omL[:, 2] - lv.omega_light_closed(K)))
res["E1 forma chiusa luce"] = {"max_err": e1, "esito": verdict(e1, 1e-12)}
# E4 birifrangenza: |omega_+| = |omega_-|, modo fermo a 0
e4 = np.max(np.abs(omL[:, 2] + omL[:, 0]))
res["E4 birifrangenza"] = {"max | |w+|-|w-| |": e4, "modo longitudinale max|w|": np.max(np.abs(omL[:, 1])),
                          "esito": verdict(e4, 1e-12)}

# E2 Dirac, forma chiusa + degenerazione di spin
th = rng.uniform(0, 1, 2000)
e2, deg = 0.0, 0.0
for i in range(2000):
    om = lv.omega_numeric(K[i:i + 1] / 2, th[i], 0.5)[0]    # k in [-pi/4, pi/4] (rami di 2k senza ripiegamento)
    e2 = max(e2, abs(om[3] - lv.omega_dirac_closed(K[i:i + 1] / 2, th[i])[0]))
    deg = max(deg, abs(om[3] - om[2]), abs(om[0] - om[1]))
res["E2 forma chiusa Dirac"] = {"max_err": e2, "spin splitting max": deg,
                                "esito": verdict(max(e2, deg), 1e-12)}

# E7 identita' omega_e(k) = omega_gamma(2k)/2
Ke = K / 2
oe = lv.omega_numeric(Ke, 0.0, 0.5)[:, 3]
og = lv.omega_numeric(2 * Ke, 0.0, 1)[:, 2]
e7 = np.max(np.abs(oe - og / 2))
res["E7 identita' elettrone = luce a 2k"] = {"max_err": e7, "esito": verdict(e7, 1e-12)}

# E3 / E5 coefficienti dello sviluppo dalla regola numerica
dirs = {"asse (1,0,0)": [1, 0, 0], "faccia (1,1,0)": [1, 1, 0], "corpo (1,1,1)": [1, 1, 1],
        "generica (1,2,3)": [1, 2, 3], "generica (3,-1,2)": [3, -1, 2]}
eps = np.array([1e-2, 5e-3, 2.5e-3])


def coeffs(fun, n):
    """fit omega/eps - 1 = b1 eps + b2 eps^2 + b4 eps^4 (3 punti piu' Richardson)"""
    n = np.asarray(n, float) / np.linalg.norm(n)
    g = np.array([fun(e * n) / e - 1 for e in eps])
    M = np.vstack([eps, eps ** 2, eps ** 4]).T
    b = np.linalg.solve(M, g)
    return b


E3, E5 = {}, {}
worst3, worst5, worst5s = 0.0, 0.0, 0.0
for name, n in dirs.items():
    A = float(lv.A_of(n))
    bl = coeffs(lambda k: lv.omega_numeric(k, 0, 1)[0, 2], n)
    be = coeffs(lambda k: lv.omega_numeric(k, 0, 0.5)[0, 3], n)
    nn = np.asarray(n, float) / np.linalg.norm(n)
    bsp = coeffs(lambda k: lv.single_tick_omega(k)[0, 2], n)
    bsm = coeffs(lambda k: lv.single_tick_omega(-k)[0, 2], n)
    pred_s = abs(nn[0] * nn[1] * nn[2]) / 2
    E3[name] = {"A": A, "luce c3 num": bl[1], "luce c3 analitico": -A / 6, "elettrone c3 num": be[1],
                "elettrone c3 analitico": -2 * A / 3}
    E5[name] = {"luce b1 (lineare in E)": bl[0], "elettrone b1": be[0],
                "un tick: b1(+k), b1(-k)": [bsp[0], bsm[0]], "un tick previsto |b1|": pred_s}
    worst3 = max(worst3, abs(bl[1] + A / 6), abs(be[1] + 2 * A / 3))
    worst5 = max(worst5, abs(bl[0]), abs(be[0]))
    worst5s = max(worst5s, abs(abs(bsp[0]) - pred_s), abs(abs(bsm[0]) - pred_s), abs(bsp[0] + bsm[0]))
symk = np.max(np.abs(lv.omega_numeric(K, 0, 1)[:, 2] - lv.omega_numeric(-K, 0, 1)[:, 2]))
res["E3 coefficiente k^3"] = {"direzioni": E3, "max |c3_num - c3_an|": worst3, "esito": verdict(worst3, 1e-6)}
res["E5 niente termine n=1"] = {"direzioni": E5, "max |b1| palindromo": worst5, "omega(k)-omega(-k) max": symk,
                                "un tick: max scarto da -+nxnynz/2": worst5s,
                                "esito": "PASS" if (worst5 < 1e-8 and symk < 1e-12 and worst5s < 1e-6) else "FAIL"}

# E6 simmetria cubica (48 operazioni)
Ks = K[:300] / 2
base_l = lv.omega_numeric(Ks, 0, 1)[:, 2]
base_d = lv.omega_numeric(Ks, 0.3, 0.5)[:, 3]
e6 = 0.0
for perm in itertools.permutations(range(3)):
    for sg in itertools.product((1, -1), repeat=3):
        Kp = Ks[:, perm] * np.array(sg)
        e6 = max(e6, np.max(np.abs(lv.omega_numeric(Kp, 0, 1)[:, 2] - base_l)),
                 np.max(np.abs(lv.omega_numeric(Kp, 0.3, 0.5)[:, 3] - base_d)))
res["E6 simmetria cubica"] = {"max_err": e6, "esito": verdict(e6, 1e-12)}

# E8 soglie (unita' scalate: m = 1, E_a = 1e4..1e6 -> la legge di scala porta ai valori fisici)
E8 = {}
xs = np.linspace(1e-4, 0.5, 20001)
worst8 = 0.0
for name in ("faccia (1,1,0)", "corpo (1,1,1)", "generica (1,2,3)"):
    n = dirs[name]
    A = float(lv.A_of(n))
    for Ea in (1e4, 1e5, 1e6):
        kan = lv.kth_photon_decay_analytic(A, Ea, m=1.0)
        f = lambda k: np.max(lv.photon_decay_excess(k, n, Ea, xs, m=1.0))
        lo, hi = 0.5 * kan, 2.0 * kan
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        knum = 0.5 * (lo + hi)
        xbest = xs[np.argmax(lv.photon_decay_excess(1.2 * knum, n, Ea, xs, m=1.0))]
        E8["%s E_a=%g" % (name, Ea)] = {"k_th numerica": knum, "k_th analitica": kan,
                                         "scarto rel": knum / kan - 1, "x migliore a 1.2 k_th": xbest}
        worst8 = max(worst8, abs(knum / kan - 1))
# non collineare: a 1.2 k_th, corpo, E_a = 1e5: aprire un angolo riduce l'eccesso?
n = np.array([1, 1, 1]) / np.sqrt(3)
Ea = 1e5
k = 1.2 * lv.kth_photon_decay_analytic(1 / 3, Ea, 1.0)
perp = np.array([1, -1, 0]) / np.sqrt(2)
best = -np.inf
best_ang = None
for x in np.linspace(0.05, 0.5, 46):
    for ang in np.linspace(0, 2e-3, 41):
        p1 = x * k * (np.cos(ang) * n + np.sin(ang) * perp)
        p2 = k * n - p1
        e1_ = lv.omega_dirac_closed(p1[None] / Ea, 1 / Ea)[0] * Ea
        e2_ = lv.omega_dirac_closed(p2[None] / Ea, 1 / Ea)[0] * Ea
        d = lv.omega_light_closed((k * n)[None] / Ea)[0] * Ea - e1_ - e2_
        if d > best:
            best, best_ang = d, (x, ang)
E8["non collineare (corpo, E_a=1e5, k=1.2 k_th)"] = {"max eccesso": best, "a (x, angolo)": best_ang}
# e -> e gamma e gamma -> gamma gamma: devono essere vietati
cher, split = -np.inf, -np.inf
for name, nd in dirs.items():
    nd = np.asarray(nd, float) / np.linalg.norm(nd)
    for pa in np.linspace(0.01, 0.5, 50):            # p a fino a 0.5
        for th_ in (1e-3, 1e-2, 0.1):
            p = pa
            Ein = lv.omega_dirac_closed((p * nd)[None], th_)[0]
            x = np.linspace(1e-3, 1, 400)
            Eout = lv.omega_dirac_closed(((1 - x) * p)[:, None] * nd[None], th_) + \
                lv.omega_light_closed((x * p)[:, None] * nd[None])
            cher = max(cher, np.max(Ein - Eout))
        Eg = lv.omega_light_closed((pa * nd)[None])[0]
        x = np.linspace(1e-3, 0.5, 400)
        Eo = lv.omega_light_closed((x * pa)[:, None] * nd[None]) + lv.omega_light_closed(((1 - x) * pa)[:, None] * nd[None])
        split = max(split, np.max(Eg - Eo))
E8["e -> e gamma max Delta E (deve essere <0)"] = cher
E8["gamma -> gamma gamma max Delta E (deve essere <=0)"] = split
res["E8 soglie"] = {"casi": E8, "max scarto rel soglia": worst8,
                    "esito": "PASS" if (worst8 < 0.01 and cher < 0 and split <= 1e-15) else "FAIL"}
res["tempo_s"] = time.time() - t0


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    return o


with open(os.path.join(HERE, "out", "lv_r1.json"), "w") as f:
    json.dump(clean(res), f, indent=1, ensure_ascii=False)
for k_, v in res.items():
    if isinstance(v, dict):
        print(k_, v.get("esito"), {kk: vv for kk, vv in v.items() if kk not in ("direzioni", "casi", "esito")})
print(json.dumps(clean(res["E3 coefficiente k^3"]["direzioni"]), indent=0))
print(json.dumps(clean(res["E5 niente termine n=1"]["direzioni"]), indent=0))
print(json.dumps(clean(res["E8 soglie"]["casi"]), indent=0))
print("tempo %.1f s" % res["tempo_s"])

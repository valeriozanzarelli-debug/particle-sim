"""geo_inverse.py - Part B1: from the masses back to the hidden shape, with a train/test split.

Hidden shape (pre-registered): a flat oblique torus, sides a = 1, b = 1.37, angle 63 deg, base mass
m0 = 0.35 * 2 pi. Its 'masses' are the Kaluza-Klein tower E = sqrt(m0^2 + |k|^2) (the rest energies of
the zig-zag ball on the hidden space, continuum limit checked in CARICA_NASCOSTA test 2), k in the
dual lattice. 40 levels with multiplicity, relative noise 1e-3.
We learn the shape from the first k levels inside a family with p parameters and PREDICT levels
k+1..40. Families: F1 square torus (p=2), F2 rectangle (p=3), F3 oblique torus (p=4, the true one),
F4 shape-free polynomial in the level index (p = P). Multistart least squares (60 starts per fit).
Non-identifiability: among all starts that fit the training data within the noise, we record the
spread of the predictions (if several shapes fit, the data did not decide).
Writes lab/out/geo_inverse.json and .log.
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from scipy.optimize import least_squares

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SEED = 20261009
NLEV = 40
NOISE = 1e-3
TRUE = dict(a=1.0, b=1.37, phi=np.deg2rad(63.0), m0=0.35 * 2 * np.pi)
NN = np.array([(i, j) for i in range(-12, 13) for j in range(-12, 13)], float)


def gauss_reduce(E):
    """Lagrange-Gauss reduction of a 2-D lattice basis (same lattice, shortest basis). Bug fix after
    run 1 (lab/out/geo_inverse_run1_bug.log): a very skewed basis of the TRUE lattice (b' = b + 5a) was
    found by the fit at k = 20, and the fixed enumeration window |n| <= 12 then missed modes."""
    u, v = E[0].copy(), E[1].copy()
    for _ in range(100):
        if u @ u > v @ v:
            u, v = v, u
        mu = np.round((u @ v) / (u @ u))
        if mu == 0:
            break
        v = v - mu * u
    return np.array([u, v])


def tower(a, b, phi, m0, K=NLEV):
    """first K masses of a flat torus with lattice vectors a(1,0), b(cos phi, sin phi)"""
    E = gauss_reduce(np.array([[a, 0.0], [b * np.cos(phi), b * np.sin(phi)]]))
    F = 2 * np.pi * np.linalg.inv(E).T            # dual vectors (rows), f_i . e_j = 2 pi delta_ij
    k2 = np.sum((NN @ F)**2, axis=1)
    return np.sqrt(m0**2 + np.sort(k2)[:K])


FAMILIES = {
    "F1 toro quadrato": (2, lambda p, K: tower(abs(p[0]), abs(p[0]), np.pi / 2, abs(p[1]), K),
                         lambda r: [r.uniform(0.3, 3), r.uniform(0.1, 6)]),
    "F2 rettangolo": (3, lambda p, K: tower(abs(p[0]), abs(p[1]), np.pi / 2, abs(p[2]), K),
                      lambda r: [r.uniform(0.3, 3), r.uniform(0.3, 3), r.uniform(0.1, 6)]),
    "F3 toro obliquo": (4, lambda p, K: tower(abs(p[0]), abs(p[1]), p[2], abs(p[3]), K),
                        lambda r: [r.uniform(0.3, 3), r.uniform(0.3, 3), r.uniform(0.3, np.pi - 0.3),
                                   r.uniform(0.1, 6)]),
}


def fit_family(name, data_tr, k, rng, starts=60):
    p, model, init = FAMILIES[name]
    sols = []
    for _ in range(starts):
        x0 = init(rng)
        try:
            r = least_squares(lambda x: np.log(model(x, k)) - np.log(data_tr), x0, method="trf",
                              max_nfev=400 * (p + 1))
        except Exception:
            continue
        rms = np.sqrt(np.mean(r.fun**2))
        sols.append((rms, r.x))
    sols.sort(key=lambda s: s[0])
    return sols, model


def poly_fit(data_tr, k, P):
    """shape-free control: E^2 = sum_{q<P} c_q n^(q/2), n = level index (Weyl-like)"""
    n = np.arange(1, k + 1, dtype=float)
    A = np.stack([n**(q / 2) for q in range(P)], axis=1)
    c, *_ = np.linalg.lstsq(A, data_tr**2, rcond=None)
    nn = np.arange(1, NLEV + 1, dtype=float)
    pred2 = np.stack([nn**(q / 2) for q in range(P)], axis=1) @ c
    return np.sqrt(np.maximum(pred2, 1e-12))


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    exact = tower(**TRUE)
    data = exact * (1 + NOISE * rng.normal(size=NLEV))
    log = []
    say = lambda s: (print(s), log.append(s))
    say(f"true masses (first 12): {np.round(exact[:12], 3).tolist()}")
    ks = [3, 4, 5, 6, 8, 10, 14, 20]
    res = {"true": {k: float(v) if np.isscalar(v) else v for k, v in TRUE.items()},
           "data": data.tolist(), "families": {}, "poly": {}}
    for name in FAMILIES:
        p = FAMILIES[name][0]
        res["families"][name] = []
        for k in ks:
            sols, model = fit_family(name, data[:k], k, rng)
            best_rms, xb = sols[0]
            pred = model(xb, NLEV)
            test = np.abs(pred[k:] / data[k:] - 1)
            # every start that fits the training data within the noise: how much do they disagree?
            good = [x for (rms, x) in sols if rms < 2 * NOISE]
            preds = np.array([model(x, NLEV) for x in good]) if good else pred[None]
            spread = float(np.median(np.std(preds[:, k:], axis=0) / data[k:])) if len(good) > 1 else 0.0
            row = dict(k=k, p=p, train_rms=float(best_rms), test_median=float(np.median(test)),
                       test_max=float(test.max()), n_good=len(good), pred_spread=spread,
                       params=np.abs(xb).tolist())
            res["families"][name].append(row)
            say(f"{name:18s} p={p} k={k:2d}  train rms {best_rms:.1e}  test median {np.median(test):.2e} "
                f"max {test.max():.2e}  fits-within-noise {len(good):2d} spread {spread:.1e}")
    for P in (2, 3, 4, 6, 8):
        res["poly"][P] = []
        for k in ks:
            if k < P:
                continue
            pred = poly_fit(data[:k], k, P)
            test = np.abs(pred[k:] / data[k:] - 1)
            train = np.abs(pred[:k] / data[:k] - 1)
            res["poly"][P].append(dict(k=k, p=P, train_max=float(train.max()),
                                       test_median=float(np.median(test)), test_max=float(test.max())))
            say(f"F4 polinomio      p={P} k={k:2d}  train max {train.max():.1e}  test median "
                f"{np.median(test):.2e}  max {test.max():.2e}")
    # pre-registered expectation
    f3 = {r["k"]: r["test_median"] for r in res["families"]["F3 toro obliquo"]}
    res["prereg_F3_pass"] = all(f3[k] < 1e-2 for k in ks if k >= 6)
    say(f"PRE-REG: F3 test median < 1e-2 for k >= 6: {res['prereg_F3_pass']}  ({time.time() - t0:.0f} s)")
    with open(os.path.join(OUT, "geo_inverse.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "geo_inverse.log"), "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()

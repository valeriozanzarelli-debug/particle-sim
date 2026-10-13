"""qg_analisi1.py - PART 1 analysis (jackknife over blocks of measurements).

    python3 lab/qg_analisi1.py <tag> [--skip N]    -> <SCR>/qg_<tag>_analisi.npz and printed summary
    (as a module: load(tag), analizza(data, cfg) used by verdict / figure / render)
"""
import glob
import json
import os
import sys

import numpy as np
from scipy.ndimage import map_coordinates

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_run import SCR, cfg  # noqa: E402


def load(tag, skip=0):
    fs = sorted(glob.glob(os.path.join(SCR, f"qg_{tag}_mis_*.npz")))
    parts = [dict(np.load(f)) for f in fs]
    keys = parts[0].keys()
    d = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    return {k: v[skip:] for k, v in d.items()}


def jk_blocks(n, nb):
    nb = min(nb, n)
    edges = np.linspace(0, n, nb + 1).astype(int)
    return [np.r_[0:edges[i], edges[i + 1]:n] for i in range(nb)]


def jackknife(fn, data, nb):
    """fn(data_subset) -> array; returns mean (full sample), error."""
    n = len(next(iter(data.values())))
    full = fn({k: v for k, v in data.items()})
    reps = np.array([fn({k: v[idx] for k, v in data.items()}) for idx in jk_blocks(n, nb)])
    m = len(reps)
    err = np.sqrt((m - 1) / m * ((reps - reps.mean(0)) ** 2).sum(0))
    return full, err


def flux_fields(D, d, T):
    """Delta of the temporal / spatial plaquette fields near the 3 quarks (canonical frame), and the two
    combinations: f = (dpt + dps)/2 = minus the action density change (positive where vacuum fluctuations are
    expelled), eps = dpt - dps ~ energy density (Minkowski E^2 + B^2, up to beta/a^4)."""
    W = D[f"Wsum_{d}_{T}"].sum()
    dpt = D[f"Ct_{d}_{T}"].sum(0) / W - D["pt_mean"].mean()
    dps = D[f"Cs_{d}_{T}"].sum(0) / W - D["ps_mean"].mean()
    return np.stack([dpt, dps, 0.5 * (dpt + dps), dpt - dps])


def punti(d):
    q = np.array([[d, 0, 0], [0, d, 0], [0, 0, d]], float)
    c = q.mean(0)
    mids = np.array([(q[0] + q[1]) / 2, (q[1] + q[2]) / 2, (q[0] + q[2]) / 2])
    return q, c, mids


def at(field, pts):
    """periodic trilinear interpolation of a (L,L,L) field at points (n,3)."""
    return map_coordinates(field, np.asarray(pts, float).T, order=1, mode="grid-wrap")


def analizza(D, c, nb=20):
    res = {"n_meas": int(len(D["plaq"]))}
    p, pe = jackknife(lambda x: np.array([x["plaq"].mean()]), D, nb)
    res["plaq"] = [float(p[0]), float(pe[0])]
    # static potential and Cornell fit
    Tf = c.get("T_pot", 3)
    Rs = np.arange(1, c["Rmax"] + 1)
    Rfit = np.array(c.get("R_fit", [2, 3, 4, 5, 6]))

    def pot(x):
        w = x["wil"].mean(0)
        return np.log(w[1:, Tf] / w[1:, Tf + 1])

    V, Ve = jackknife(pot, D, nb)
    res["V"] = V.tolist()
    res["V_err"] = Ve.tolist()

    def cornell(x):
        v = pot(x)[Rfit - 1]
        A = np.stack([np.ones(len(Rfit)), -1.0 / Rfit, Rfit.astype(float)], 1)
        wgt = 1 / np.maximum(Ve[Rfit - 1], 1e-9)
        sol = np.linalg.lstsq(A * wgt[:, None], v * wgt, rcond=None)[0]
        return np.array([sol[0], sol[1], sol[2], np.sqrt(max(sol[2], 0))])

    fp, fe = jackknife(cornell, D, nb)
    res["cornell"] = {"V0": [fp[0], fe[0]], "alpha": [fp[1], fe[1]], "sigma": [fp[2], fe[2]], "a_sqrt_sigma": [fp[3], fe[3]]}
    # 3Q energies
    res["E3Q"] = {}
    for d in c["d_list"]:
        def e3(x, d=d):
            w = x[f"wb_{d}"].real.mean(0)
            return np.log(w[:-1] / w[1:])
        e, ee = jackknife(e3, D, nb)
        res["E3Q"][d] = {"T": c["Tb_list"][:-1], "E": e.tolist(), "err": ee.tolist()}
    # flux
    res["flux"] = {}
    for d in c["corr_d"]:
        for T in c["corr_T"]:
            q, cen, mids = punti(d)
            orig = np.zeros((1, 3))

            def fl(x, d=d, T=T, cen=cen, mids=mids):
                F = flux_fields(x, d, T)
                out = []
                for k in range(4):
                    fc = at(F[k], cen[None])[0]
                    fm = at(F[k], mids).mean()
                    fo = at(F[k], orig)[0]
                    out += [fc, fm, fo]
                f = F[2]
                fc, fm = out[6], out[7]
                return np.array(out + [fc / fm if fm != 0 else np.nan, fc - fm])

            v, ve = jackknife(fl, D, nb)
            Ffull = flux_fields(D, d, T)

            def flf(x, d=d, T=T):
                return flux_fields(x, d, T)

            _, Ferr = jackknife(flf, D, nb)
            names = ["pt", "ps", "f", "eps"]
            r = {}
            for k, nm in enumerate(names):
                r[nm] = {"centroide": [v[3 * k], ve[3 * k]], "medi_lati": [v[3 * k + 1], ve[3 * k + 1]],
                         "origine": [v[3 * k + 2], ve[3 * k + 2]]}
            r["R_Y"] = [v[12], ve[12]]
            r["diff_f"] = [v[13], ve[13]]
            res["flux"][f"{d}_{T}"] = r
            res[f"F_{d}_{T}"] = Ffull
            res[f"Ferr_{d}_{T}"] = Ferr
    return res


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items() if not isinstance(v, np.ndarray) or v.ndim < 2}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


if __name__ == "__main__":
    tag = sys.argv[1]
    skip = int(sys.argv[sys.argv.index("--skip") + 1]) if "--skip" in sys.argv else 0
    c = cfg(tag)
    D = load(tag, skip)
    res = analizza(D, c, nb=min(20, len(D["plaq"])))
    np.savez(os.path.join(SCR, f"qg_{tag}_campi.npz"), **{k: v for k, v in res.items() if isinstance(v, np.ndarray)})
    js = _clean(res)
    json.dump(js, open(os.path.join(SCR, f"qg_{tag}_analisi.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in js.items() if k not in ("flux",)}, indent=None)[:3000])
    for k, v in js["flux"].items():
        print(k, "f:", v["f"], "R_Y", v["R_Y"], "eps", v["eps"]["centroide"], v["eps"]["medi_lati"])

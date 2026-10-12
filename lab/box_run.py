"""box_run.py - A1 (Schwinger in the strip: rate vs m_eff, formula with pi/d, suppression at fixed E) and A3
(crossover width) of lab/out/box_prereg.json + box_amend.json. Writes lab/out/box_run.json. (~3 min, 1 thread)"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys, time
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from box_core import static_part, static_part_reflect, masses, min_gap, schwinger_strip, floquet, PI

OUT = os.path.join(HERE, "out")
TH = 0.1
DS = [4, 6, 8, 12, 16, 24, 32, 64, 128]


def p_infmass(d, m):
    """continuum Dirac strip with infinite-mass walls: p d = pi - arctan(p/m)"""
    return brentq(lambda p: p * d - (PI - np.arctan(p / m)), 1e-9, PI / d * 1.0 + 1e-9)


if __name__ == "__main__":
    t0 = time.time()
    res = dict(th=TH, A1a={}, A1c={}, meff={})
    # ---- A1a / A1b
    for d in DS:
        M = static_part(TH, d)
        m = masses(TH, d, M)
        g, kmin = min_gap(TH, d, M)
        p = p_infmass(d, TH)
        res["meff"][d] = dict(m_eff=float(m[0]), modes=m[:6].tolist(), min_gap_over_kx=g, kx_of_min=kmin,
                              formula_pi_over_d=float(np.sqrt(TH**2 + (PI / d) ** 2)),
                              formula_pi_over_2d=float(np.sqrt(TH**2 + (PI / (2 * d)) ** 2)),
                              formula_infinite_mass=float(np.sqrt(TH**2 + p**2)))
        for x in (1, 2, 4):
            E = x * m[0] ** 2
            P0, Pf, n = schwinger_strip(TH, d, E, M=M)
            pred = float(np.exp(-PI * m**2 / E).sum())
            res["A1a"][f"d{d}_x{x}"] = dict(d=d, x=x, E=E, ticks=n, P_lat=P0, P_pred=pred, ratio=P0 / pred,
                                            P_pred_lowest_only=float(np.exp(-PI * m[0] ** 2 / E)),
                                            P_1D_same_E=float(np.exp(-PI * TH**2 / E)))
            print(d, x, round(E, 4), n, round(P0, 5), round(pred, 5), "ratio", round(P0 / pred, 4), f"{time.time()-t0:.0f}s",
                  flush=True)
    res["A1a_PASS"] = bool(all(abs(v["ratio"] - 1) <= 0.10 for v in res["A1a"].values()))
    res["A1a_fails"] = [k for k, v in res["A1a"].items() if abs(v["ratio"] - 1) > 0.10]
    dev = {d: v["m_eff"] / v["formula_pi_over_d"] - 1 for d, v in res["meff"].items()}
    res["A1b_dev_vs_pi_over_d"] = dev
    res["A1b_PASS"] = bool(all(abs(x) <= 0.10 for x in dev.values()))
    res["A1b_dev_vs_infinite_mass"] = {d: v["m_eff"] / v["formula_infinite_mass"] - 1 for d, v in res["meff"].items()}
    res["A1b_dev_vs_pi_over_2d"] = {d: v["m_eff"] / v["formula_pi_over_2d"] - 1 for d, v in res["meff"].items()}
    print("A1a PASS", res["A1a_PASS"], res["A1a_fails"], "A1b PASS", res["A1b_PASS"], {d: round(x, 3) for d, x in dev.items()},
          flush=True)
    # ---- A1c fixed E = 2 th^2
    E = 2 * TH**2
    P1D = float(np.exp(-PI * TH**2 / E))
    for d in DS:
        M = static_part(TH, d)
        m = masses(TH, d, M)
        P0, Pf, n = schwinger_strip(TH, d, E, M=M)
        res["A1c"][d] = dict(P_lat=P0, P_pred=float(np.exp(-PI * m**2 / E).sum()), P_1D=P1D,
                             P_lat_per_P1D=P0 / P1D)
        print("A1c", d, P0, res["A1c"][d]["P_pred"], P0 / P1D, flush=True)
    Ps = [res["A1c"][d]["P_lat"] for d in DS]
    mono = bool(all(Ps[i] < Ps[i + 1] for i in range(len(Ps) - 1)))
    small = bool(all(res["A1c"][d]["P_lat_per_P1D"] <= 0.1 for d in DS if d <= 8))
    res["A1c_monotonic"], res["A1c_factor10_at_d_le_8"] = mono, small
    res["A1c_PASS"] = mono and small
    # per transverse width: pairs per unit area also (P/d), the rate per unit x-length is P * E/(2 pi)
    # ---- A3 crossover: m_eff(d) fine
    dd = np.arange(2, 81)
    me = np.array([masses(TH, int(d))[0] for d in dd])
    res["A3_curve"] = dict(d=dd.tolist(), m_eff=me.tolist())
    target = np.sqrt(2) * TH
    i = int(np.where(me < target)[0][0])
    dstar = float(dd[i - 1] + (target - me[i - 1]) * (dd[i] - dd[i - 1]) / (me[i] - me[i - 1]))
    res["A3_dstar"] = dstar
    res["A3_dstar_th"] = dstar * TH
    res["A3_PASS"] = bool(0.5 <= dstar * TH <= 5)
    print("A3 d* =", dstar, "d* th =", dstar * TH, res["A3_PASS"], flush=True)
    # ---- the pre-registered reflecting wall (with edge states): A1a numbers for the record
    old = {}
    for d in (16, 64):
        M = static_part_reflect(TH, d)
        w, _ = floquet(M, 0.0)
        m = np.sort(w[(w > 0) & (w < PI / 2)])
        for x in (1, 2, 4):
            E = x * m[0] ** 2
            P0, Pf, n = schwinger_strip(TH, d, E, K=1.2, M=M)
            old[f"d{d}_x{x}"] = dict(ratio=P0 / float(np.exp(-PI * m**2 / E).sum()))
    res["A1a_old_reflecting_wall"] = old
    print("old wall", old)
    res["time_s"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "box_run.json"), "w"), indent=1, default=float)

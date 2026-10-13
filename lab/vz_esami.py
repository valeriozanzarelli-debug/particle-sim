"""vz_esami.py - VIDEO SUPERCRITICO: the pre-registered exams V1..V5 (lab/out/vz_prereg.json) -> lab/out/vz_esami.json

    python3 lab/vz_esami.py      (~1-2 min, 1 core)
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vz_comune import (N, NF, T, DT, RB, RC, RHO, J, QE, SC, RW, ER, ZA, N_NONMONO, palline, raggi)  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "vz_esami.json")


def main():
    t0 = time.time()
    m, _ = palline()
    R = np.empty((NF, N), np.float32)
    for k in range(NF):
        R[k] = raggi(m, T[k])
    res = {"N": N, "celle_Q_non_monotone_ripulite": N_NONMONO}

    # V1: radial histogram vs rho
    edges = np.concatenate([[0.0], np.geomspace(0.01, 40.0, 41)])
    ov = []
    for k in range(NF):
        r = R[k][np.isfinite(R[k])]
        if r.size < 3000:
            ov.append(None)
            continue
        hb = np.histogram(r, edges)[0] / r.size
        q = np.interp(edges, RB, QE[k])
        hd = np.diff(q)
        hd = hd / hd.sum()
        ov.append(float(1 - 0.5 * np.abs(hb - hd).sum()))
    v = np.array([x for x in ov if x is not None])
    res["V1_istogramma"] = {"sovrapposizione_min": float(v.min()), "mediana": float(np.median(v)),
                            "fotogrammi_valutati": int(v.size), "min_al_fotogramma": int(np.nanargmin([x if x is not None else np.nan for x in ov])),
                            "PASS": bool(v.min() >= 0.95 and np.median(v) >= 0.98)}
    res["serie_sovrapposizione"] = [None if x is None else round(x, 5) for x in ov]

    # V2: RK4 on v = j/rho for 2000 balls vs the mass coordinate
    rng = np.random.default_rng(7)
    sel = rng.choice(N, 2000, replace=False)
    r = R[0][sel].astype(np.float64)
    re = RB[1:]

    def vel(rr, t):
        k = int(min(max(np.ceil(t / DT - 1e-12), 1), NF - 1))       # j[k] holds on (t_{k-1}, t_k]
        u = (t - T[k - 1]) / DT
        rho = (1 - u) * RHO[k - 1] + u * RHO[k]
        rh = np.interp(rr, RC, rho)
        jj = np.interp(rr, np.concatenate([[0.0], re]), np.concatenate([[0.0], J[k]]))
        return jj / np.maximum(rh, 1e-30)

    h = 0.05
    t = 0.0
    chk = {}
    for tt in (40.0, 80.0, 160.0):
        while t < tt - 1e-9:
            alive = r < 40.0
            ra = r[alive]
            k1 = vel(ra, t)
            k2 = vel(ra + 0.5 * h * k1, t + 0.5 * h)
            k3 = vel(ra + 0.5 * h * k2, t + 0.5 * h)
            k4 = vel(ra + h * k3, t + h)
            r[alive] = np.maximum(ra + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 1e-6)
            t += h
        kk = int(round(tt / DT))
        rq = R[kk][sel]
        ok = np.isfinite(rq) & (r < 40.0)
        dev = np.abs(r[ok] - rq[ok]) / rq[ok]
        chk[str(int(tt))] = {"palline": int(ok.sum()), "mediana": float(np.median(dev)), "frazione_oltre_0.10": float((dev > 0.10).mean()),
                             "solo_una_delle_due_dentro_40": int((np.isfinite(rq) != (r < 40.0)).sum())}
    res["V2_ode"] = {**chk, "PASS": bool(all(c["mediana"] <= 0.03 and c["frazione_oltre_0.10"] <= 0.10 for c in chk.values()))}

    # V3: charge
    fin = np.isfinite(R)
    p_in = ((R < RW) & fin).mean(1)
    p30 = 1 - ((R < 30.0) & fin).mean(1)
    q40 = fin.mean(1)
    d1 = np.abs(p_in - np.array(SC["P_in_rw"]))
    d2 = np.abs(p30 - np.array(SC["passato_r30"]))
    d3 = np.abs(q40 - np.array(SC["Q_entro_40"]))
    res["V3_carica"] = {"max_P_in_rw": float(d1.max()), "max_passato_r30": float(d2.max()), "max_Q_entro_40": float(d3.max()),
                        "PASS": bool(max(d1.max(), d2.max(), d3.max()) <= 0.01)}
    res["serie_palline"] = {"P_in_rw": [round(float(x), 5) for x in p_in], "passato_r30": [round(float(x), 5) for x in p30],
                            "Q_entro_40": [round(float(x), 5) for x in q40]}

    # V4/V5: positron speed
    vs, ns = 0.0, 0
    sj, sr = 0.0, 0.0
    ce = (RB[1:] >= 20.0) & (RB[1:] <= 30.0)
    cc = (RC >= 20.0) & (RC <= 30.0)
    for k in range(NF - 1):
        if not (40.0 <= T[k] < 160.0):
            continue
        a, b = R[k], R[k + 1]
        ok = np.isfinite(a) & np.isfinite(b) & (a >= 20.0) & (a <= 30.0)
        vs += float(((b[ok] - a[ok]) / DT).sum())
        ns += int(ok.sum())
        sj += float(J[k + 1][ce].sum())
        sr += float((0.5 * (RHO[k] + RHO[k + 1]))[cc].sum())
    vb = vs / ns
    vd = sj / sr
    eps = abs(ER) - ZA / 25.0
    beta = float(np.sqrt(1 - 1 / eps ** 2))
    res["V4_velocita_dati"] = {"palline": vb, "dati_sum_j_su_sum_rho": vd, "scarto": vb / vd - 1, "campioni": ns,
                               "PASS": bool(abs(vb / vd - 1) <= 0.03)}
    res["V5_velocita_cinematica"] = {"palline": vb, "beta_cinematica_r25": beta, "eps": eps, "scarto": vb / beta - 1,
                                     "documento_dice": "~0,8 c", "PASS": bool(abs(vb / beta - 1) <= 0.05)}
    # speeds seen on screen (for captions): max |v| of balls between frames
    vmax = 0.0
    for k in range(NF - 1):
        ok = np.isfinite(R[k]) & np.isfinite(R[k + 1])
        vmax = max(vmax, float(np.abs(R[k + 1][ok] - R[k][ok]).max() / DT))
    res["v_max_palline_fra_fotogrammi"] = vmax
    res["secondi"] = time.time() - t0
    json.dump(res, open(OUT, "w"), indent=1)
    for kk in ("V1_istogramma", "V2_ode", "V3_carica", "V4_velocita_dati", "V5_velocita_cinematica"):
        print(kk, json.dumps(res[kk]))
    print("vmax", vmax, "nonmono", N_NONMONO)


if __name__ == "__main__":
    main()

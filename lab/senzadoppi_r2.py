"""senzadoppi_r2.py - exams E4 (candidates PAW, B2W and references F0, PA: valleys, rest mass, isotropy, kinetic
mass, maximal group velocity, unitarity, Bloch == real space) and E5 (time reversal: Kramers pairs).
Output lab/out/senzadoppi_r2.json.  Run: OMP_NUM_THREADS=2 timeout 590 python lab/senzadoppi_r2.py"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import senzadoppi_core as s
import senzadoppi_real as sr
from senzadoppi_r1 import RULES, VEC_A, VEC_B

OUT = os.path.join(os.path.dirname(__file__), "out", "senzadoppi_r2.json")
T0 = time.time()

REAL = {  # name -> (ticks, mode, h, vecs)
    "F0": ([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], "each", 0.0, []),
    "PA": ([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], "once", 0.0, []),
    "PAW": ([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], "once", 0.3, VEC_A),
    "B2W": ([((0, 1, 2), (2, 1, 0))], "each", 0.3, VEC_B),
}


def vgroup_max(R, theta, n=28, dk=1e-4):
    K = s.zone_grid(n)
    vmax = 0.0
    E0 = s.quasi(R.bloch(K, theta), R.T)
    for i in range(3):
        d = np.zeros(3)
        d[i] = dk
        Ep = s.quasi(R.bloch(K + d, theta), R.T)
        Em = s.quasi(R.bloch(K - d, theta), R.T)
        g = (Ep - Em) / (2 * dk)
        gaps = E0[:, 2] - E0[:, 1]
        ok = (gaps > 1e-2) & (np.abs(np.abs(E0).max(1) - np.pi / R.T) > 1e-2)
        if i == 0:
            G2 = np.zeros_like(g)
            OK = ok
        G2 += g**2
        OK &= ok
    v = np.sqrt(G2)[OK]
    # refine along body diagonal and axes of the worst point
    return float(np.max(v))


def bloch_vs_real(name, theta=0.2, n=8):
    ticks, mode, h, vecs = REAL[name]
    RR = sr.RealRule(ticks, theta, mode, h, vecs)
    rng = np.random.default_rng(5)
    kk = 2 * np.pi * np.array([1, 2, 3]) / n
    x = np.arange(n)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s0 = rng.normal(size=4) + 1j * rng.normal(size=4)
    psi = s0[:, None, None, None] * np.exp(1j * (kk[0] * X + kk[1] * Y + kk[2] * Z))[None]
    out = RR.period(psi)
    pred = RULES[name].bloch(kk[None], theta)[0] @ s0
    err_bloch = float(np.max(np.abs(out[:, 0, 0, 0] - pred)))
    # unitarity with a random potential, 50 periods, 16^3
    m = 16
    V = rng.random((m, m, m))
    RR = sr.RealRule(ticks, theta, mode, h, vecs, V=V)
    psi = rng.normal(size=(4, m, m, m)) + 1j * rng.normal(size=(4, m, m, m))
    n0 = np.sum(np.abs(psi)**2)
    for _ in range(50):
        psi = RR.period(psi)
    err_u = float(abs(np.sum(np.abs(psi)**2) - n0) / n0)
    return err_bloch, err_u


def kramers(name, n=6, theta=0.37):
    ticks, mode, h, vecs = REAL[name]
    rng = np.random.default_rng(11)
    V = rng.uniform(-0.5, 0.5, (n, n, n))
    RR = sr.RealRule(ticks, theta, mode, h, vecs, V=V)
    D = 4 * n**3
    M = np.zeros((D, D), complex)
    for i in range(D):
        e = np.zeros(D, complex)
        e[i] = 1
        M[:, i] = RR.period(e.reshape(4, n, n, n)).ravel()
    w = np.linalg.eigvals(M)
    ph = np.sort(np.angle(w))
    # pair consecutive phases (Kramers doublets), split = max over pairs of the distance inside the pair
    splits = ph[1::2] - ph[0::2]
    return float(np.max(np.abs(splits))), float(np.max(np.abs(np.abs(w) - 1)))


def main():
    res = {}
    for nm in ("F0", "PA", "PAW", "B2W"):
        R = RULES[nm]
        r = {}
        for th in (0.1, 0.02):
            vs = s.valley_scan(R, th, n=32)
            rk = s.rest_and_kinetic(R, th, dk=min(1e-3, th / 50))
            r[f"theta{th}"] = {"E0": vs["E0"], "Epi": vs["Epi"], "k_E0": vs["k_E0"], "k_Epi": vs["k_Epi"],
                               "rest_err": rk[0], "mkin_err_max": rk[1], "mkin_err_mean": rk[2]}
            print(nm, th, r[f"theta{th}"], round(time.time() - T0), flush=True)
        r["copies"] = vs["n_copies"]
        iso = s.isotropy(R, 24, 600)
        r["iso24"] = {"maxerr": iso[0], "mean": iso[1], "vmin": iso[2], "vmax": iso[3]}
        r["vgroup_max_theta0.1"] = vgroup_max(R, 0.1)
        r["bloch_vs_real"], r["unitarity_50"] = bloch_vs_real(nm)
        r["kramers_split"], r["kramers_unit"] = kramers(nm)
        print(nm, {k: v for k, v in r.items() if not k.startswith("theta")}, round(time.time() - T0), flush=True)
        res[nm] = r
    iso_F0 = res["F0"]["iso24"]["maxerr"]
    verdict = {}
    for nm in ("PAW", "B2W"):
        r = res[nm]
        a, b = r["theta0.1"], r["theta0.02"]
        C = {
            "C1": a["E0"] >= 0.3 and b["E0"] >= 0.3,
            "C2": a["Epi"] >= 0.3 and b["Epi"] >= 0.3,
            "C3": abs(a["rest_err"]) <= 1e-12 and abs(b["rest_err"]) <= 1e-12,
            "C4": r["iso24"]["maxerr"] <= 1.1 * iso_F0,
            "C5": a["mkin_err_max"] <= 2e-2 and a["mkin_err_max"] / max(b["mkin_err_max"], 1e-300) >= 15,
            "C6": r["vgroup_max_theta0.1"] <= 1.0,
            "C7": r["unitarity_50"] <= 1e-12,
            "C8": r["bloch_vs_real"] <= 1e-12,
        }
        verdict[nm] = {k: bool(v) for k, v in C.items()}
        verdict[nm]["PROPOSTA_forte"] = all(C[k] for k in ("C1", "C2", "C3", "C4", "C5"))
    E5 = {nm: res[nm]["kramers_split"] for nm in res}
    E5pass = E5["F0"] <= 1e-12 and E5["PA"] <= 1e-12 and E5["PAW"] >= 1e-4 and E5["B2W"] >= 1e-4
    out = {"results": res, "E4_verdict": verdict, "E5": E5, "E5_PASS": bool(E5pass), "time_s": time.time() - T0}
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(verdict, indent=1), E5, E5pass)


if __name__ == "__main__":
    main()

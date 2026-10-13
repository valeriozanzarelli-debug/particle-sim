"""senzadoppi_r4.py - exam E7: does the 1s of a Z = 60 nucleus "beat" (valley multiplet) on the candidate rule?
40^3 cells, a = 0,5 lambda_C (theta = 0,5 as the Laboratorio), nucleus on the corner of 8 cells (origin), point
Coulomb V = -Z alpha / r as a phase with the coin.  Initial state: continuum Dirac 1s spinor, spin up (a TEST state
only: it is not an input of the rule).  1600 periods; survival c(t) = <psi0|psi(t)>, p_exc = 1 - |c|^2, line weights
from a windowed Fourier transform of c(t).
Run: OMP_NUM_THREADS=2 timeout 590 python lab/senzadoppi_r4.py <F0|PA|PAW>; then python lab/senzadoppi_r4.py sum"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import zz_core as zz
import senzadoppi_real as sr
from senzadoppi_r2 import REAL

OUTD = os.path.join(os.path.dirname(__file__), "out")
SCR = "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad"
Z, A, N, NPER = 60, 0.5, 40, 1600
ZA = Z * zz.ALPHA


def initial(n=N, a=A):
    x = (np.arange(n) - n / 2 + 0.5) * a
    X, Y, Zc = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Zc**2)
    gam = np.sqrt(1 - ZA**2)
    g = r**(gam - 1) * np.exp(-ZA * r)
    f = -np.sqrt((1 - gam) / (1 + gam)) * g
    nx, ny, nz = X / r, Y / r, Zc / r
    # chi = spin up (1, 0); sigma.n chi = (nz, nx + i ny)
    u = np.array([g, 0 * g])
    dd = 1j * f * np.array([nz, nx + 1j * ny])
    L = (u - dd) / np.sqrt(2)
    R = (u + dd) / np.sqrt(2)
    psi = np.concatenate([L, R]).astype(complex)
    return psi / np.sqrt(np.sum(np.abs(psi)**2)), r


def run(name):
    if name.startswith("FW1"):
        import senzadoppi_r5  # noqa: F401  (registers FW1 in REAL, amendment A3)
    ticks, mode, h, vecs = REAL[name]
    psi0, r = initial()
    V = -ZA / r * A                          # lattice units (energy x a): the phase per tick is V a
    RR = sr.RealRule(ticks, A, mode, h, vecs, V=V)
    T = len(ticks) * A                       # physical duration of a period (lambda_C / c)
    psi = psi0.copy()
    c = np.zeros(NPER + 1, complex)
    c[0] = 1
    t0 = time.time()
    for p in range(1, NPER + 1):
        psi = RR.period(psi)
        c[p] = np.vdot(psi0, psi)
    norm = float(np.sum(np.abs(psi)**2))
    np.save(os.path.join(SCR, f"senzadoppi_r4_{name}.npy"), c)
    json.dump({"T": T, "norm_end": norm, "time_s": time.time() - t0}, open(os.path.join(SCR, f"senzadoppi_r4_{name}.json"), "w"))
    print(name, "done", time.time() - t0, "norm", norm, "max p_exc", float(np.max(1 - np.abs(c)**2)))


def lines(c, T, Emin=0.6, Emax=1.05, dE=2e-5):
    t = np.arange(len(c)) * T
    w = np.hanning(len(c))
    Es = np.arange(Emin, Emax, dE)
    S = np.abs(np.exp(1j * np.outer(Es, t)) @ (c * w)) / np.sum(w)
    pk = [i for i in range(1, len(Es) - 1) if S[i] >= S[i - 1] and S[i] >= S[i + 1] and S[i] >= 0.02]
    return Es, S, [(float(Es[i]), float(S[i])) for i in pk]


def summary():
    res = {}
    for name in ("F0", "PA", "PAW", "FW1"):
        fn = os.path.join(SCR, f"senzadoppi_r4_{name}.npy")
        if not os.path.exists(fn):
            continue
        c = np.load(fn)
        meta = json.load(open(os.path.join(SCR, f"senzadoppi_r4_{name}.json")))
        pex = 1 - np.abs(c)**2
        Es, S, pk = lines(c, meta["T"])
        pk_sorted = sorted(pk, key=lambda x: -x[1])
        main = pk_sorted[0] if pk_sorted else (None, 0)
        # main line weight, adding a partner within 5e-3 (spin doublet)
        wmain = main[1] + sum(w for e, w in pk_sorted[1:] if abs(e - main[0]) <= 5e-3)
        res[name] = {"max_p_exc": float(pex.max()), "max_p_exc_after_t100": float(pex[int(100 / meta["T"]):].max()),
                     "mean_p_exc": float(pex.mean()), "lines_ge_2pc": pk_sorted, "main_line": main,
                     "main_weight_with_doublet": float(wmain), "norm_end": meta["norm_end"], "T": meta["T"],
                     "dirac_1s": float(zz.dirac_formula(1, -1, ZA))}
        print(name, json.dumps(res[name]))
    if "PAW" in res:
        r = res["PAW"]
        res["E7_PASS"] = bool(r["max_p_exc"] <= 0.15 and r["main_weight_with_doublet"] >= 0.8)
    if "FW1" in res:
        r = res["FW1"]
        res["E7_PASS_FW1"] = bool(r["max_p_exc"] <= 0.15 and r["main_weight_with_doublet"] >= 0.8)
    json.dump(res, open(os.path.join(OUTD, "senzadoppi_r4.json"), "w"), indent=1)
    print("E7", res.get("E7_PASS"))


if __name__ == "__main__":
    if sys.argv[1] == "sum":
        summary()
    else:
        run(sys.argv[1])

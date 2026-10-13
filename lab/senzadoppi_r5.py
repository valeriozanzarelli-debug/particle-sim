"""senzadoppi_r5.py - amendment A3: candidate FW1 = the frozen rule (orders xyz / zyx, coin theta each tick) plus ONE
twisted Wilson per period (chirality axis Y, links +-2 e_j, half before and half after the coin of the EVEN ticks only).
Exams E4 (C1-C8), E5 (Kramers) for h = 0,4 (production, fixed in A3) and h = 0,3 (reported), E6 (pairs over the whole
zone, h = 0,4).  Output lab/out/senzadoppi_r5.json.
Run: OMP_NUM_THREADS=2 timeout 590 python lab/senzadoppi_r5.py   (rerun until exit 0)"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import senzadoppi_core as s
import senzadoppi_r1 as r1
import senzadoppi_r2 as r2
import senzadoppi_r3 as r3

OUT = os.path.join(os.path.dirname(__file__), "out", "senzadoppi_r5.json")
CK = os.path.join(r3.os.path.dirname(r3.CK), "senzadoppi_r5_ck.json")
T0 = time.time()
TICKS = [((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))]


class FW1(s.Rule):
    def __init__(self, h, vecs=r1.VEC_A, name="FW1"):
        super().__init__(TICKS, name=name)
        self.hh, self.vecs = h, vecs

    def bloch(self, K, theta):
        K = np.atleast_2d(np.asarray(K, float))
        C = s.coin(theta)[None]
        W = self.tick_trans(K, 0)
        W = s.wilson_links(K, self.hh / 2, vecs=self.vecs[::-1]) @ C @ s.wilson_links(K, self.hh / 2, vecs=self.vecs) @ W
        W = C @ self.tick_trans(K, 1) @ W
        return W


for h, nm in ((0.4, "FW1"), (0.3, "FW1h03")):
    r1.RULES[nm] = FW1(h, name=nm)
    r2.REAL[nm] = (TICKS, "each_w0", h, r1.VEC_A)
    r3.SPEC[nm] = (TICKS, "each_w0", h, r1.VEC_A, 8)


# r3's evolve/period_matrix know 'each' and 'once'; FW1 needs 'each_w0' -> local versions
def coin_blk(K, n, h):
    if n % 2 == 0:
        return r3.coin_block(K, r3.THETA, h, r1.VEC_A, False)
    return np.broadcast_to(s.coin(r3.THETA)[None], (K.shape[0], 4, 4))


def evolve(name, K0, E, nticks):
    h = r3.SPEC[name][2]
    R = s.Rule(TICKS)
    A = r3.kx_profile(E, nticks)
    U = np.broadcast_to(np.eye(4, dtype=complex), (K0.shape[0], 4, 4)).copy()
    for n in range(nticks):
        K = K0.copy()
        K[:, 0] += A[2 * n + 1]
        U = coin_blk(K, n, h) @ R.tick_trans(K, n % 2) @ U
    return U


def period_matrix(name, K):
    h = r3.SPEC[name][2]
    R = s.Rule(TICKS)
    U = coin_blk(K, 0, h) @ R.tick_trans(K, 0)
    return coin_blk(K, 1, h) @ R.tick_trans(K, 1) @ U


r3.evolve = evolve
r3.period_matrix = period_matrix


def main():
    d = json.load(open(CK)) if os.path.exists(CK) else {}
    for nm in ("FW1", "FW1h03"):
        if nm in d:
            continue
        R = r1.RULES[nm]
        r = {}
        for th in (0.1, 0.02):
            vs = s.valley_scan(R, th, n=32)
            rk = s.rest_and_kinetic(R, th, dk=min(1e-3, th / 50))
            rkF = s.rest_and_kinetic(r1.RULES["F0"], th, dk=min(1e-3, th / 50))
            r[f"theta{th}"] = {"E0": vs["E0"], "Epi": vs["Epi"], "k_E0": vs["k_E0"], "k_Epi": vs["k_Epi"],
                               "rest_err": rk[0], "mkin_err_max": rk[1], "mkin_err_mean": rk[2],
                               "mkin_minus_F0": rk[1] - rkF[1]}
            print(nm, th, r[f"theta{th}"], round(time.time() - T0), flush=True)
        vs5 = s.valley_scan(R, 0.5, n=28)
        r["theta0.5"] = {"E0": vs5["E0"], "Epi": vs5["Epi"], "mkin": s.rest_and_kinetic(R, 0.5)[1]}
        iso = s.isotropy(R, 24, 600)
        r["iso24"] = {"maxerr": iso[0], "mean": iso[1], "vmin": iso[2], "vmax": iso[3]}
        r["vgroup_max_theta0.1"] = r2.vgroup_max(R, 0.1)
        r["bloch_vs_real"], r["unitarity_50"] = r2.bloch_vs_real(nm)
        r["kramers_split"], r["kramers_unit"] = r2.kramers(nm)
        print(nm, {k: v for k, v in r.items() if not k.startswith("theta")}, round(time.time() - T0), flush=True)
        d[nm] = r
        json.dump(d, open(CK, "w"))
    # E6 pairs (h = 0,4)
    nlines = r3.NPERP**2
    for fac in (1.0, 2.0):
        key = f"pairs_E{fac}"
        st = d.setdefault(key, {"sum": 0.0, "done": 0})
        while st["done"] < nlines:
            if time.time() - T0 > 540:
                json.dump(d, open(CK, "w"))
                print("checkpoint", key, st["done"], flush=True)
                sys.exit(3)
            idx = np.arange(st["done"], min(st["done"] + r3.CHUNK, nlines))
            P, Eeff = r3.run_chunk("FW1", fac * r3.THETA**2, idx)
            st["sum"] += float(np.sum(P))
            st["Eeff"] = Eeff
            st["done"] = int(idx[-1] + 1)
            print(key, st["done"], round(time.time() - T0), flush=True)
        json.dump(d, open(CK, "w"))
    dA = (2 * np.pi / r3.NPERP)**2
    for fac in (1.0, 2.0):
        st = d[f"pairs_E{fac}"]
        E = st["Eeff"]
        ref = 2 * E * np.exp(-np.pi * r3.THETA**2 / E)
        st["R"] = st["sum"] * dA / (8 * ref)
    iso_F0 = 0.01543382769313073
    verdict = {}
    for nm in ("FW1", "FW1h03"):
        r = d[nm]
        a, b = r["theta0.1"], r["theta0.02"]
        C = {"C1": a["E0"] >= 0.3 and b["E0"] >= 0.3, "C2": a["Epi"] >= 0.3 and b["Epi"] >= 0.3,
             "C3": abs(a["rest_err"]) <= 1e-12 and abs(b["rest_err"]) <= 1e-12,
             "C4": r["iso24"]["maxerr"] <= 1.1 * iso_F0,
             "C5": a["mkin_err_max"] <= 2e-2 and a["mkin_err_max"] / b["mkin_err_max"] >= 15,
             "C6": r["vgroup_max_theta0.1"] <= 1.0, "C7": r["unitarity_50"] <= 1e-12, "C8": r["bloch_vs_real"] <= 1e-12}
        verdict[nm] = {k: bool(v) for k, v in C.items()}
        verdict[nm]["PROPOSTA_forte"] = all(C[k] for k in ("C1", "C2", "C3", "C4", "C5"))
    out = {"results": {k: d[k] for k in ("FW1", "FW1h03")}, "E4_verdict": verdict,
           "E5_kramers": {nm: d[nm]["kramers_split"] for nm in ("FW1", "FW1h03")},
           "E6_pairs": {k: {"R": d[k]["R"], "E_over_Ec": d[k]["Eeff"] / r3.THETA**2} for k in ("pairs_E1.0", "pairs_E2.0")},
           "E6_PASS": bool(all(abs(d[k]["R"] - 1) <= 0.1 for k in ("pairs_E1.0", "pairs_E2.0")))}
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out["E4_verdict"]), out["E5_kramers"], out["E6_pairs"], out["E6_PASS"])


if __name__ == "__main__":
    main()

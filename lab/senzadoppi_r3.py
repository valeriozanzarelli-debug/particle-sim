"""senzadoppi_r3.py - exam E6: pairs from the vacuum in a uniform electric field, counted over the WHOLE zone with no
valley projector.  Gauge of the vector potential: at tick n every factor of the rule sees k_x -> k_x0 + E (n + 1/2).
Each mode (k_y, k_z) on a grid of the whole transverse zone starts in the two E < 0 bands at k_x0 = -pi; after a sweep
of 2 pi (2 pi / E ticks) we project on the E > 0 bands at the same k.  Sum over the grid -> number of pairs.
Reference: one Dirac species, Landau-Zener-Schwinger per crossing: 2 spins x E exp(-pi m^2 / E) (k_perp integral);
the zone contains n_copies exact copies of Gamma, each crossed once.  R = total / (n_copies x reference).
Output lab/out/senzadoppi_r3.json.  Run: OMP_NUM_THREADS=2 timeout 590 python lab/senzadoppi_r3.py  (rerun until 0)"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import senzadoppi_core as s
from senzadoppi_r1 import VEC_A, VEC_B

OUT = os.path.join(os.path.dirname(__file__), "out", "senzadoppi_r3.json")
CK = os.path.join("/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad",
                  "senzadoppi_r3b_ck.json")
T0 = time.time()
BUDGET = 540
THETA = 0.1
NPERP = 160
CHUNK = 6400
KX0 = -0.75 * np.pi - 0.3   # start/end far from every valley (amendment A2)

SPEC = {  # name -> (ticks, mode, h, vecs, n_copies)
    "F0": ([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], "each", 0.0, [], 8),
    "PA": ([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], "once", 0.0, [], 8),
    "PAW": ([((0, 1, 2), (0, 1, 2)), ((2, 1, 0), (2, 1, 0))], "once", 0.3, VEC_A, 8),
}


def coin_block(K, theta_c, h, vecs, mirror):
    C = s.coin(theta_c)[None]
    if not h:
        return np.broadcast_to(C, (K.shape[0], 4, 4))
    v = vecs[::-1] if mirror else vecs
    return s.wilson_links(K, h / 2, vecs=v[::-1]) @ C @ s.wilson_links(K, h / 2, vecs=v)


RAMP = 120          # ticks of the smooth (sin^2) switch-on and switch-off of the field (amendment A2)


def kx_profile(E, nticks):
    """k_x displacement A(t) at half-integer and integer times: field E f(t), f = sin^2 ramps of RAMP ticks,
    total displacement exactly 2 pi (E is rescaled for that)"""
    t = np.arange(2 * nticks + 1) / 2.0
    f = np.ones_like(t)
    up = t < RAMP
    dn = t > nticks - RAMP
    f[up] = np.sin(np.pi * t[up] / (2 * RAMP))**2
    f[dn] = np.sin(np.pi * (nticks - t[dn]) / (2 * RAMP))**2
    A = np.concatenate([[0], np.cumsum((f[1:] + f[:-1]) / 2 * 0.5)])
    return A * (2 * np.pi / A[-1])


def evolve(name, K0, E, nticks):
    ticks, mode, h, vecs, _ = SPEC[name]
    T = len(ticks)
    R = s.Rule(ticks)
    N = K0.shape[0]
    A = kx_profile(E, nticks)
    U = np.broadcast_to(np.eye(4, dtype=complex), (N, 4, 4)).copy()
    for n in range(nticks):
        K = K0.copy()
        K[:, 0] += A[2 * n + 1]
        U = R.tick_trans(K, n % T) @ U
        if mode == "each":
            U = coin_block(K, THETA, h, vecs, (n % T) == 1) @ U
        elif (n % T) == T - 1:
            Kc = K0.copy()
            Kc[:, 0] += A[2 * n + 2]
            U = coin_block(Kc, T * THETA, h, vecs, False) @ U
    return U


def period_matrix(name, K):
    ticks, mode, h, vecs, _ = SPEC[name]
    R = s.Rule(ticks)
    T = len(ticks)
    N = K.shape[0]
    U = np.broadcast_to(np.eye(4, dtype=complex), (N, 4, 4)).copy()
    for n in range(T):
        U = R.tick_trans(K, n) @ U
        if mode == "each":
            U = coin_block(K, THETA, h, vecs, n == 1) @ U
    if mode == "once":
        U = coin_block(K, T * THETA, h, vecs, False) @ U
    return U


def run_chunk(name, E, idx):
    g = (np.arange(NPERP) + 0.5) / NPERP * 2 * np.pi - np.pi
    KY, KZ = np.meshgrid(g, g, indexing="ij")
    ky, kz = KY.ravel()[idx], KZ.ravel()[idx]
    K0 = np.stack([np.full_like(ky, KX0), ky, kz], 1)
    T = len(SPEC[name][0])
    nper = int(round((2 * np.pi / E + RAMP) / T))     # plateau field = E (ramps of RAMP ticks lose RAMP/2 each)
    A = kx_profile(E, nper * T)
    Eeff = float(np.max(np.diff(A)) * 2)               # plateau field actually used
    W0 = period_matrix(name, K0)
    w, v = np.linalg.eig(W0)
    qe = -np.angle(w)
    lo = np.argsort(qe, axis=1)
    # lower two (E < 0) and upper two eigenvectors; orthonormalize within each pair
    Vlo = np.take_along_axis(v, lo[:, None, :2], axis=2)
    Vhi = np.take_along_axis(v, lo[:, None, 2:], axis=2)
    Vlo, _ = np.linalg.qr(Vlo)
    Vhi, _ = np.linalg.qr(Vhi)
    U = evolve(name, K0, Eeff, nper * T)
    out = U @ Vlo
    P = np.sum(np.abs(np.conj(np.transpose(Vhi, (0, 2, 1))) @ out)**2, axis=(1, 2))
    return P, Eeff


def main():
    d = json.load(open(CK)) if os.path.exists(CK) else {}
    nlines = NPERP * NPERP
    for name in ("PAW", "PA", "F0"):
        for fac in (1.0, 2.0):
            key = f"{name}_E{fac}"
            st = d.setdefault(key, {"sum": 0.0, "done": 0, "Pgrid": []})
            while st["done"] < nlines:
                if time.time() - T0 > BUDGET:
                    json.dump(d, open(CK, "w"))
                    print("checkpoint", key, st["done"], flush=True)
                    sys.exit(3)
                idx = np.arange(st["done"], min(st["done"] + CHUNK, nlines))
                P, Eeff = run_chunk(name, fac * THETA**2, idx)
                st["sum"] += float(np.sum(P))
                st["Eeff"] = Eeff
                st["done"] = int(idx[-1] + 1)
                st["Pgrid"] += [float(x) for x in P]
                print(key, st["done"], round(time.time() - T0), flush=True)
            json.dump(d, open(CK, "w"))
    res = {}
    dA = (2 * np.pi / NPERP)**2
    for key, st in d.items():
        name = key.split("_")[0]
        E = st["Eeff"]
        ref = 2 * E * np.exp(-np.pi * THETA**2 / E)
        total = st["sum"] * dA
        Pg = np.array(st["Pgrid"]).reshape(NPERP, NPERP)
        res[key] = {"E_over_Ec": E / THETA**2, "total": total, "ref_one_species_one_crossing": ref,
                    "n_copies": SPEC[name][4], "R": total / (SPEC[name][4] * ref),
                    "max_P_line": float(Pg.max()), "frac_lines_P_gt_0.5": float(np.mean(Pg > 0.5))}
        print(key, res[key], flush=True)
    passE6 = all(abs(res[f"PAW_E{f}"]["R"] - 1) <= 0.1 for f in (1.0, 2.0))
    out = {"theta": THETA, "nperp": NPERP, "results": res, "E6_PASS": bool(passE6), "time_s_last": time.time() - T0}
    json.dump(out, open(OUT, "w"), indent=1)
    print("E6 PASS", passE6)


if __name__ == "__main__":
    main()

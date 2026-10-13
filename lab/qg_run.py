"""qg_run.py - PART 1 Monte Carlo in pieces (<= budget seconds per call), resumable.

    python3 lab/qg_run.py <tag> [--budget 500]

Config from lab/out/qg_prereg.json ("parte1") or, for the pilot, from the built-in PILOT dict (tag 'pilota').
State (links, rng, counters) in <SCR>/qg_<tag>_state.npz; measurements appended to <SCR>/qg_<tag>_mis_XXXX.npz
(one file per call, never rewritten).
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_mc import Lat  # noqa: E402

LAB = os.path.dirname(os.path.abspath(__file__))
SCR = os.environ.get("QG_SCR", "/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/qg")
os.makedirs(SCR, exist_ok=True)

PILOT = dict(T=12, L=12, beta=5.7, seed=11, therm=60, sweep_or=3, between=2, n_meas=10 ** 6,
             ape4_alpha=0.5, ape4_n=2, apes_alpha=0.5, apes_n=12, Rmax=6, Tmax=5,
             d_list=[2, 3, 4], Tb_list=[1, 2, 3, 4, 5], corr_d=[3, 4], corr_T=[2, 4])

OCT = [(sx, sy, sz) for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]


def cfg(tag):
    if tag == "pilota":
        return PILOT
    pr = json.load(open(os.path.join(LAB, "out", "qg_prereg.json")))
    return pr["parte1"]["config"]


def canon(C, sgn):
    """C_s(y) -> C_can(y) = C_s(S y): reflect the axes with sign -1 (y -> -y mod L)."""
    for k, s in enumerate(sgn):
        if s < 0:
            C = np.roll(np.flip(C, axis=k), 1, axis=k)
    return C


def measure(lat, c):
    U = lat.U
    Um = lat.ape(U, c["ape4_alpha"], c["ape4_n"], (0, 1, 2, 3)) if c["ape4_n"] > 0 else U
    Us = lat.ape(Um, c["apes_alpha"], c["apes_n"], (1, 2, 3))
    out = {"plaq": lat.plaquette(U), "plaq_m": lat.plaquette(Um)}
    out["wil"] = lat.wilson_planar(Us, Um, c["Rmax"], c["Tmax"])
    pt, ps = lat.site_plaq(Um)
    out["pt_mean"], out["ps_mean"] = pt.mean(), ps.mean()
    Fpt = {}
    Fps = {}
    L = lat.L
    for d in c["d_list"]:
        wb = np.zeros(len(c["Tb_list"]), complex)
        for T in c["corr_T"]:
            out[f"Ct_{d}_{T}"] = np.zeros((L, L, L))
            out[f"Cs_{d}_{T}"] = np.zeros((L, L, L))
        for sgn in OCT:
            W = lat.baryon(Us, Um, d, c["Tb_list"], sgn)
            for i, T in enumerate(c["Tb_list"]):
                wb[i] += W[T].mean() / len(OCT)
            if d in c["corr_d"]:
                for T in c["corr_T"]:
                    h = T // 2
                    if T not in Fpt:
                        Fpt[T] = np.fft.fftn(np.roll(pt, -h, axis=0), axes=(1, 2, 3))
                        Fps[T] = np.fft.fftn(np.roll(ps, -h, axis=0), axes=(1, 2, 3))
                    FW = np.conj(np.fft.fftn(W[T].real, axes=(1, 2, 3)))
                    ct = np.fft.ifftn(FW * Fpt[T], axes=(1, 2, 3)).real.sum(0)
                    cs = np.fft.ifftn(FW * Fps[T], axes=(1, 2, 3)).real.sum(0)
                    out[f"Ct_{d}_{T}"] += canon(ct, sgn)
                    out[f"Cs_{d}_{T}"] += canon(cs, sgn)
            if d in c["corr_d"]:
                for T in c["corr_T"]:
                    out.setdefault(f"Wsum_{d}_{T}", 0.0)
                    out[f"Wsum_{d}_{T}"] += W[T].real.sum()
        out[f"wb_{d}"] = wb
    return out


def main():
    tag = sys.argv[1]
    budget = float(sys.argv[sys.argv.index("--budget") + 1]) if "--budget" in sys.argv else 500.0
    c = cfg(tag)
    t0 = time.time()
    st = os.path.join(SCR, f"qg_{tag}_state.npz")
    lat = Lat(c["T"], c["L"], c["beta"], seed=c["seed"])
    nsw, nmeas = 0, 0
    if os.path.exists(st):
        z = np.load(st, allow_pickle=True)
        lat.U = z["U"]
        lat.rng.bit_generator.state = z["rng"].item()
        nsw, nmeas = int(z["nsw"]), int(z["nmeas"])
    recs = []
    plaqs = []

    def unit():
        lat.heatbath()
        for _ in range(c["sweep_or"]):
            lat.overrelax()

    while time.time() - t0 < budget and nmeas < c["n_meas"]:
        if nsw < c["therm"]:
            unit()
            nsw += 1
            if nsw % 10 == 0:
                lat.reunit()
            plaqs.append((nsw, lat.plaquette()))
            continue
        for _ in range(c["between"]):
            unit()
            nsw += 1
        lat.reunit()
        m = measure(lat, c)
        m["nsw"] = nsw
        recs.append(m)
        nmeas += 1
        print(f"meas {nmeas} sw {nsw} plaq {m['plaq']:.5f} t {time.time() - t0:.0f}s", flush=True)
    np.savez(st + ".tmp.npz", U=lat.U, rng=np.array(lat.rng.bit_generator.state, dtype=object), nsw=nsw, nmeas=nmeas)
    if recs:
        k = len([f for f in os.listdir(SCR) if f.startswith(f"qg_{tag}_mis_")])
        keys = recs[0].keys()
        np.savez(os.path.join(SCR, f"qg_{tag}_mis_{k:04d}.npz"), **{kk: np.array([r[kk] for r in recs]) for kk in keys})
    if plaqs:
        with open(os.path.join(SCR, f"qg_{tag}_therm.txt"), "a") as f:
            for a, b in plaqs:
                f.write(f"{a} {b:.6f}\n")
    os.replace(st + ".tmp.npz", st)
    print("done", tag, "nsw", nsw, "nmeas", nmeas, f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()

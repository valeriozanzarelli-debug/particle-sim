"""qg_gluoni.py - PART 2 runs: two gluon wave packets collide (classical SU(3) Yang-Mills, real time).

    python3 lab/qg_gluoni.py <run> [--budget 480]

Runs are defined in lab/out/qg_prereg.json ("parte2" -> "corse"). Resumable in pieces: state in
<SCR>/qg_ym_<run>_state.npz. Diagnostics every `diag_every` steps -> <SCR>/qg_ym_<run>_diag.json.
Final (and, for video runs, per-frame) gauge-invariant fields -> <SCR>/qg_ym_<run>_*.npz.
For video runs the drawn balls are carried INSIDE the run by v = S / eps (Poynting / energy density, trilinear),
so they move with the field's own energy flow; positions saved per frame.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import map_coordinates

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_pacchetti import build  # noqa: E402
from qg_ym import YM  # noqa: E402

LAB = os.path.dirname(os.path.abspath(__file__))
SCR = os.environ.get("QG_SCR", "/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/qg")


def run_cfg(name):
    pr = json.load(open(os.path.join(LAB, "out", "qg_prereg.json")))
    p2 = pr["parte2"]
    c = dict(p2["comune"])
    r = dict(p2["corse"][name])
    if "come" in r:                       # same packets as another run
        r["pacchetti"] = [dict(p) for p in p2["corse"][r["come"]]["pacchetti"]]
    if "gen2" in r:                       # same-colour control: second packet along the same generator
        r["pacchetti"][1]["gen"] = r["gen2"]
    c.update(r)
    return c


def make_packets(c):
    X, Y, Z = c["shape"]
    out = []
    for k, pk in enumerate(c["pacchetti"]):
        d = dict(pk)
        d["amp"] = c["amp"]
        d["k0"] = 2 * np.pi / c["lam"]
        d["sig_par"], d["sig_perp"] = c["sig_par"], c["sig_perp"]
        out.append(d)
    return out


def interp(field, pos):
    """periodic trilinear interpolation of a site field at ball positions pos (3, N)."""
    return map_coordinates(field, pos, order=1, mode="grid-wrap")


def main():
    name = sys.argv[1]
    budget = float(sys.argv[sys.argv.index("--budget") + 1]) if "--budget" in sys.argv else 480.0
    c = run_cfg(name)
    t0 = time.time()
    shape = tuple(c["shape"])
    dt = c["dt"]
    st = os.path.join(SCR, f"qg_ym_{name}_state.npz")
    dg = os.path.join(SCR, f"qg_ym_{name}_diag.json")
    ym = YM(shape)
    rng = np.random.default_rng(c.get("seed", 5))
    video = c.get("video", False)
    if os.path.exists(st):
        z = np.load(st)
        ym.U, ym.E, n = z["U"], z["E"], int(z["n"])
        balls = z["balls"] if video else None
        diag = json.load(open(dg))
    else:
        build(ym, make_packets(c))
        n = 0
        diag = {"t": [], "H": [], "HE": [], "HB": [], "gauss_rel": [], "gauss_abs": [], "unit_err": []}
        balls = None
        if video:
            eE, eB = ym.energy_fields()
            e = (eE + eB).ravel()
            idx = rng.choice(e.size, size=c["n_balls"], p=e / e.sum())
            pos = np.array(np.unravel_index(idx, shape), float) + rng.random((3, c["n_balls"])) - 0.5
            balls = pos % np.array(shape)[:, None]
        ym.E -= 0.5 * dt * ym.force()            # E lives at half steps
    nsteps = c["steps"]
    fe = c.get("frame_every", 0)
    while n < nsteps and time.time() - t0 < budget:
        F = ym.force()
        need_diag = n % c["diag_every"] == 0
        need_frame = video and fe and n % fe == 0
        if need_diag or need_frame or video:
            Eh = ym.E + 0.5 * dt * F               # E at the integer time of U
            eE, eB = ym.energy_fields(Eh)
            e = eE + eB
        if need_diag:
            G = ym.gauss()
            ue = 0.0
            for i in range(3):
                Ui = ym.U[i]
                UU = np.einsum("ji...,jk...->ik...", np.conj(Ui), Ui)
                ue = max(ue, float(np.abs(UU - np.eye(3)[:, :, None, None, None]).max()))
            diag["t"].append(n * dt)
            diag["H"].append(float(e.sum()))
            diag["HE"].append(float(eE.sum()))
            diag["HB"].append(float(eB.sum()))
            diag["gauss_rel"].append(float(np.sqrt((np.abs(G) ** 2).sum() / (np.abs(Eh) ** 2).sum())))
            diag["gauss_abs"].append(float(np.sqrt((np.abs(G) ** 2).sum())))
            diag["unit_err"].append(ue)
            print(f"{name} n {n} t {n * dt:.2f} H {e.sum():.6g} gauss {diag['gauss_rel'][-1]:.2e} {time.time() - t0:.0f}s", flush=True)
        if video:
            Eold = ym.E
            ym.E = Eh
            S = ym.poynting()
            ym.E = Eold
            if need_frame:
                k = n // fe
                ratio = eE / np.maximum(e, 1e-30)
                q = np.round(balls / np.array(shape)[:, None] * 65535).astype(np.uint16)
                # fidelity: balls per coarse cell vs field energy per coarse cell (Pearson r)
                cb = c["fid_cell"]
                hb = np.histogramdd(balls.T, bins=[np.arange(0, s_ + 1, cb) for s_ in shape])[0]
                he = e.reshape(shape[0] // cb, cb, shape[1] // cb, cb, shape[2] // cb, cb).sum((1, 3, 5))
                fid = float(np.corrcoef(hb.ravel(), he.ravel())[0, 1])
                diag.setdefault("frame_t", []).append(n * dt)
                diag.setdefault("fid_r", []).append(fid)
                np.savez(os.path.join(SCR, f"qg_ym_{name}_f{k:04d}.npz"), q=q,
                         rb=interp(ratio, balls).astype(np.float16), eb=interp(e, balls).astype(np.float16),
                         t=n * dt, fid=fid)
            floor = c["eps_floor"] * e.mean()
            v = [interp(S[k] / (e + floor), balls) for k in range(3)]
            v = np.array(v)
            sp = np.sqrt((v ** 2).sum(0))
            v *= np.minimum(1.0, 1.0 / np.maximum(sp, 1e-30))
            # midpoint (RK2) with the same frozen field over one step: dt small against field changes
            mid = (balls + 0.5 * dt * v) % np.array(shape)[:, None]
            v2 = np.array([interp(S[k] / (e + floor), mid) for k in range(3)])
            sp2 = np.sqrt((v2 ** 2).sum(0))
            v2 *= np.minimum(1.0, 1.0 / np.maximum(sp2, 1e-30))
            balls = (balls + dt * v2) % np.array(shape)[:, None]
        ym.step(dt, F)
        n += 1
        if n % 25 == 0:
            ym.reunit()
    done = n >= nsteps
    if done:
        Eh = ym.E + 0.5 * dt * ym.force()
        eE, eB = ym.energy_fields(Eh)
        Eold = ym.E
        ym.E = Eh
        S = ym.poynting()
        ym.E = Eold
        np.savez(os.path.join(SCR, f"qg_ym_{name}_fine.npz"), eE=eE.astype(np.float32), eB=eB.astype(np.float32),
                 S=S.astype(np.float32), t=n * dt)
    json.dump(diag, open(dg + ".tmp", "w"))
    np.savez(st + ".tmp.npz", U=ym.U, E=ym.E, n=n, balls=balls if video else np.zeros(1))
    os.replace(dg + ".tmp", dg)
    os.replace(st + ".tmp.npz", st)
    print("done" if done else "pausa", name, n, f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

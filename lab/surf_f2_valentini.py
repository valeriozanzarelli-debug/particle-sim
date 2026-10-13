"""surf_f2_valentini.py n M - PHASE 2, part 1: relaxation to quantum equilibrium (Valentini) on the 2-D zig-zag lattice
with the EXACT deterministic flow (no random jumps). Exams V1-V7 in lab/out/surf_f2_prereg.json.
usage: python surf_f2_valentini.py 48|96 1|4|16|25   -> lab/out/surf_f2_val_<n>_<M>.json
       python surf_f2_valentini.py summary          -> lab/out/surf_f2_valentini.json
"""
import json
import sys
import time
import numpy as np
from surf_core import *
from pre_core import LatticeWalk, modes_state, hbar, fit_tau

CFG = {48: dict(theta=0.25, cell=4, ticks=4000, every=20), 96: dict(theta=0.125, cell=8, ticks=8000, every=40)}
NS = 50000


def tick_flux(wk, psi):
    F = []
    for j in (0, 1):
        psi, (A, B, Ap, Bp) = wk.axis_step(psi, j, need_flux=True)
        F.append(Ap - B)
    return wk.coin(psi), F


def fits(t, h, floor):
    t = np.asarray(t, float); y = np.asarray(h) - floor
    y0 = y[0]
    m = (y > 0.05 * y0) & (y > 0)
    out = dict(tau=fit_tau(t, h, floor), ratio=float(y[-1] / y0) if y0 > 0 else None)
    if m.sum() >= 4:
        A = np.vstack([t[m], np.ones(m.sum())]).T
        c, res, *_ = np.linalg.lstsq(A, np.log(y[m]), rcond=None)
        ss = np.sum((np.log(y[m]) - np.log(y[m]).mean())**2)
        out["exp_rate"] = float(-c[0]); out["exp_R2"] = float(1 - res[0] / ss) if len(res) and ss > 0 else None
        mm = m & (t > 0)
        if mm.sum() >= 4:
            A = np.vstack([np.log(t[mm]), np.ones(mm.sum())]).T
            c2, res2, *_ = np.linalg.lstsq(A, np.log(y[mm]), rcond=None)
            ss2 = np.sum((np.log(y[mm]) - np.log(y[mm]).mean())**2)
            out["pow_p"] = float(-c2[0]); out["pow_R2"] = float(1 - res2[0] / ss2) if len(res2) and ss2 > 0 else None
    return out


def run(n, M):
    c = CFG[n]
    t0 = time.time()
    wk = LatticeWalk(n, c["theta"], periodic=False)
    psi, om = modes_state(wk, M, omega0=2 * c["theta"])
    rho = density(psi)
    rng = np.random.default_rng(SEED + 100 * M + n)
    ens = {}
    ens["born"] = born_sample(rng, rho, NS)
    ens["rho2"] = born_sample(rng, rho**2, NS)
    ens["uniform"] = born_sample(rng, np.ones_like(rho), NS)
    half = rho.copy(); half[n // 2:, :] = 0
    ens["half"] = born_sample(rng, half, NS)
    K = (n // c["cell"])**2
    floor = (K - 1) / (2 * NS)
    rec = {k: [] for k in ens}
    ts = []
    for t in range(c["ticks"] + 1):
        if t % c["every"] == 0:
            r = density(psi)
            ts.append(t)
            for k, (cl, xi) in ens.items():
                rec[k].append(hbar(cl, r, c["cell"])[0])
        if t == c["ticks"]:
            break
        r0 = density(psi)
        psi, F = tick_flux(wk, psi)
        r1 = density(psi)
        for k, (cl, xi) in ens.items():
            flow_tick(cl, xi, r0, r1, F, periodic=True)
    sig = float(np.std(om) / c["theta"]) if M > 1 else 0.0
    res = dict(n=n, M=M, theta=c["theta"], cell=c["cell"], ticks=c["ticks"], floor=floor, omega=om.tolist(),
               sigma_omega_per_time=sig, t=ts, H=rec, time_s=time.time() - t0,
               fits={k: fits(ts, v, floor) for k, v in rec.items()})
    for k in res["fits"]:
        if res["fits"][k]["tau"] is not None:
            res["fits"][k]["tau_phys"] = res["fits"][k]["tau"] * c["theta"]
            res["fits"][k]["tau_sigma"] = res["fits"][k]["tau_phys"] * sig
    with open(os.path.join(OUT, f"surf_f2_val_{n}_{M}.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    print(n, M, "sigma", round(sig, 4), {k: {q: (round(v, 4) if isinstance(v, float) else v) for q, v in fv.items()}
                                          for k, fv in res["fits"].items()}, "maxH born", max(rec["born"]), "3F", 3 * floor,
          "time", round(time.time() - t0), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "summary":
        pass
    else:
        run(int(sys.argv[1]), int(sys.argv[2]))

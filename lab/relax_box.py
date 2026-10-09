"""relax_box.py - Test 1 and 2 of relax_prereg: relaxation to |psi|^2 in a 2-D box.

usage: python lab/relax_box.py RUN [RUN ...]   (RUN like M16_psi4, M25_unif, M16_eq, M16_cos2, M16_cos12)
       python lab/relax_box.py --quick M25_psi4   (timing probe)
writes lab/out/relax_box_<RUN>.npz
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from relax_core import Box, sample, rk4_step, coarse, Hbar, vortices, PI

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
N = 20000
T_END = 4 * PI       # prereg said 8 pi; cut to 4 pi for the time budget (timing probe: 11 s per unit time)
DT = 0.02
EVERY = 10           # output every 0.2
NG = 32              # finest coarse-graining grid
SUB = 8              # sub-points per cell for |psi|^2 cell mass


def run(name, t_end=T_END, n=N):
    M = int(name.split("_")[0][1:]); start = name.split("_")[1]
    box = Box(M, seed=1)
    rng = np.random.default_rng(100 + M)
    d0 = lambda x, y: np.abs(box.field(x, y, 0.0, grad=False)) ** 2
    gm = np.abs(box.dens_grid(0.0, 256)) ** 2
    if start == "psi4":
        x, y = sample(rng, n, lambda x, y: d0(x, y) ** 2, 1.2 * gm.max() ** 2)
    elif start == "unif":
        x, y = rng.uniform(0, PI, n), rng.uniform(0, PI, n)
    elif start == "eq":
        x, y = sample(rng, n, d0, 1.2 * gm.max())
    elif start.startswith("cos"):
        k = int(start[3:])
        x, y = sample(rng, n, lambda x, y: d0(x, y) * (1 + 0.8 * np.cos(k * x)), 1.2 * 1.8 * gm.max())
    t = 0.0; ts = []; H = {4: [], 8: [], 16: [], 32: []}; vort = []
    t0 = time.time(); step = 0
    while True:
        if step % EVERY == 0:
            cnt, _, _ = np.histogram2d(x, y, bins=NG, range=[[0, PI], [0, PI]])
            psi = box.dens_grid(t, NG * SUB)
            ref = (np.abs(psi) ** 2).reshape(NG, SUB, NG, SUB).sum((1, 3))
            for K in H:
                c, r = coarse(cnt, ref, K)
                H[K].append(Hbar(c, r))
            ts.append(t)
            if step % (EVERY * 6) == 0:
                vort.append(vortices(box.dens_grid(t, 128)))
            if step == 0:
                x0, y0 = x.copy(), y.copy()
        if t >= t_end - 1e-9:
            break
        x, y = rk4_step(box, x, y, t, DT, ds=0.02)
        t += DT; step += 1
    dt = time.time() - t0
    floor = {K: (K * K - 1) / (2 * n) for K in H}
    print(f"{name}: H16 {H[16][0]:.4f} -> {H[16][-1]:.4f} (floor {floor[16]:.4f}), vortices {np.mean(vort):.1f}, {dt:.0f}s",
          flush=True)
    np.savez(os.path.join(OUT, f"relax_box_{name}.npz"), t=np.array(ts),
             **{f"H{K}": np.array(H[K]) for K in H}, vort=np.array(vort), x=x, y=y, x0=x0, y0=y0,
             M=M, n=n, secs=dt, t_end=t_end)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    args = sys.argv[1:]
    if args[0] == "--quick":
        run(args[1], t_end=1.0, n=N)
    else:
        for a in args:
            run(a)

"""zz_slit.py - T1 L5: Young's double slit with the LIGHT walk (spin-1 walkers, F = E + i c B).

Same geometry as docs/REGOLA_BASE.md (lab/optics.py): lambda = 1, two slits 2 wide with centres 8 apart,
Gaussian beam waist W0 = 7 on the slit plane, screen at L = 160. 2-D (x, z), fields uniform in y,
E along y (TE). The wall is one cell thick and absorbs every ball that lands on it (F set to 0);
the box edges are sponges (balls fade out). A source plane at z = -3 adds a forward-moving
(S_z = +1) wave every tick (cos omega t, with a smooth switch-on).
Measured: time-averaged Poynting flux c (E x B)_z on the screen (E = Re F, cB = Im F) over two periods
in the steady state, binned like the photon tests, vs the exact 2-D Rayleigh-Sommerfeld flux.
Tests (tolerances fixed before the first run of THIS script; same thresholds as the photon rule):
  L5a coarse correlation >= 0.98 at lambda = 8a;  L5b bright-fringe rms error <= 8 % at lambda = 8a.
  Convergence: other lambda/a given on the command line (python lab/zz_slit.py 16) -> zz_slit_16.json.
(The test was first written inside zz_light.py: too slow there in numpy (753 s estimated), moved here
with a numba kernel; same rule, same tolerances.)
Run: timeout 590 python lab/zz_slit.py
"""
import json, os, sys, time
import numpy as np
import numba as nb
from zz_core import Walk, BASIS_ONE, SEED
from optics import screen_reference, XS, BINS, HALF, COARSE, LAM, W0, D, A, L, fringe_windows

TOL = dict(corr=0.98, fringe_rms=0.08)
OUT = os.path.join(os.path.dirname(__file__), "out")
res, verdict = {}, []
t_start = time.time()


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


@nb.njit(parallel=True, cache=True)
def axis_step(F, up, um, axis):
    """one axis translation: the S=+1 part moves +1 cell, the S=-1 part -1 cell (zero beyond the edge)"""
    nx, nz = F.shape[1], F.shape[2]
    cp = np.empty((nx, nz), np.complex64)
    cm = np.empty((nx, nz), np.complex64)
    for i in nb.prange(nx):
        for k in range(nz):
            a0, a1, a2 = F[0, i, k], F[1, i, k], F[2, i, k]
            cp[i, k] = np.conj(up[0]) * a0 + np.conj(up[1]) * a1 + np.conj(up[2]) * a2
            cm[i, k] = np.conj(um[0]) * a0 + np.conj(um[1]) * a1 + np.conj(um[2]) * a2
    for i in nb.prange(nx):
        for k in range(nz):
            if axis == 0:
                sp = cp[i - 1, k] if i > 0 else np.complex64(0)
                sm = cm[i + 1, k] if i < nx - 1 else np.complex64(0)
            else:
                sp = cp[i, k - 1] if k > 0 else np.complex64(0)
                sm = cm[i, k + 1] if k < nz - 1 else np.complex64(0)
            dp = sp - cp[i, k]
            dm = sm - cm[i, k]
            for q in range(3):
                F[q, i, k] += up[q] * dp + um[q] * dm


def tick(F, n):
    for j in ((0, 2) if n % 2 == 0 else (2, 0)):
        lam, P = BASIS_ONE[j]
        axis_step(F, P[:, 0].astype(np.complex64), P[:, 2].astype(np.complex64), 0 if j == 0 else 1)
    return F


# machinery check against the generic walk (field kept away from the edges, so periodic == open)
rng = np.random.default_rng(SEED)
G = np.zeros((3, 30, 26), complex)
G[:, 8:22, 8:18] = rng.normal(size=(3, 14, 10)) + 1j * rng.normal(size=(3, 14, 10))
Wg = Walk("light", 30, 1.0, axes=(0, 2), axmap={0: 0, 2: 1})
a_, b_ = G.copy(), G.astype(np.complex64)
for t in range(5):
    a_ = Wg.tick(a_)
    b_ = tick(b_, t)
print(f"numba 2-D tick vs generic walk: max diff {np.max(np.abs(a_ - b_)):.1e} (float32)", flush=True)


def double_slit(per_lam, Ls=L, xhalf=67.0, half=HALF):
    """Ls: slit -> screen distance; xhalf: half width of the box; half: half width of the screen"""
    a = 1.0 / per_lam
    zlo, zhi, sponge = -6.0, Ls + 7.0, 5.0
    nx = int(round(2 * xhalf / a))
    nz = int(round((zhi - zlo) / a))
    x = (np.arange(nx) - nx / 2 + 0.5) * a
    z = zlo + (np.arange(nz) + 0.5) * a
    F = np.zeros((3, nx, nz), np.complex64)
    dx_edge = np.clip((np.abs(x) - (xhalf - sponge)) / sponge, 0, 1)
    dz_edge = np.clip(np.maximum((z - (zhi - sponge)) / sponge, (-3.0 - z) / 3.0), 0, 1)
    mx = np.exp(-8 * dx_edge**2 * a).astype(np.float32)
    mz = np.exp(-8 * dz_edge**2 * a).astype(np.float32)
    ix_s = np.nonzero(mx < 1)[0]
    iz_s = np.nonzero(mz < 1)[0]
    iz_wall = int(np.argmin(np.abs(z - 0.0)))
    blocked = ~((np.abs(x - D / 2) <= A / 2) | (np.abs(x + D / 2) <= A / 2))
    iz_src = int(np.argmin(np.abs(z - (-3.0))))
    iz_scr = int(np.argmin(np.abs(z - Ls)))
    u_plus = BASIS_ONE[2][1][:, 0]
    prof = np.exp(-x**2 / W0**2)
    omega = 2 * np.pi * a / LAM
    path = np.sqrt(Ls**2 + half**2) + 3.0
    n_ramp = int(6 * per_lam)
    n_run = int((path + 6.0) / a) + n_ramp
    n_avg = 2 * per_lam
    Sz = np.zeros(nx)
    t0 = time.time()
    for t in range(n_run + n_avg):
        tick(F, t)
        amp = np.sin(0.5 * np.pi * min(1.0, t / n_ramp))**2 * np.cos(omega * t)
        for q in range(3):
            F[q, :, iz_src] += (amp * u_plus[q] * prof).astype(np.complex64)
        F[:, blocked, iz_wall] = 0
        F[:, ix_s, :] *= mx[ix_s][None, :, None]
        F[:, :, iz_s] *= mz[iz_s][None, None, :]
        if t >= n_run:
            E = F[:, :, iz_scr].real.astype(float)
            cB = F[:, :, iz_scr].imag.astype(float)
            Sz += E[0] * cB[1] - E[1] * cB[0]
        if t == 30:
            print(f"  lambda/a = {per_lam}: grid {nx} x {nz}, {n_run + n_avg} ticks, estimated {(time.time() - t0) / 31 * (n_run + n_avg):.0f} s", flush=True)
    nb_ = int(round(2 * half / 0.5))
    edges = np.linspace(-half, half, nb_ + 1)
    xs = 0.5 * (edges[1:] + edges[:-1])
    if Ls == L and half == HALF:
        ref = screen_reference("flux")
    else:
        from optics import rs_field, sources, flux
        s0, w0 = sources(400)
        xf = np.linspace(-half, half, nb_ * 8 + 1)
        xf = 0.5 * (xf[1:] + xf[:-1])
        ref = flux(*rs_field(xf, Ls, s0, w0))[1].reshape(nb_, 8).sum(1)
        ref = ref / ref.sum()
    p = np.array([Sz[(x >= edges[i]) & (x < edges[i + 1])].sum() for i in range(nb_)])
    neg = float(-p[p < 0].sum() / np.abs(p).sum())
    p = np.clip(p, 0, None)
    p = p / p.sum()
    corr_coarse = np.corrcoef(p.reshape(-1, COARSE).sum(1), ref.reshape(-1, COARSE).sum(1))[0, 1] if nb_ % COARSE == 0 else np.nan
    corr_fine = np.corrcoef(p, ref)[0, 1]
    if Ls == L and half == HALF:
        errs = np.array([p[(XS >= lo) & (XS < hi)].sum() / ref[(XS >= lo) & (XS < hi)].sum() - 1 for lo, hi in fringe_windows(ref)])
    else:
        errs = np.array([np.nan])
    return dict(per_lambda=per_lam, screen=Ls, grid=[nx, nz], ticks=n_run + n_avg, corr_coarse=float(corr_coarse), corr_fine=float(corr_fine),
                fringe_rms=float(np.sqrt(np.mean(errs**2))), fringe_max=float(np.max(np.abs(errs))),
                fringe_errs=errs.tolist(), negative_flux_fraction=neg, seconds=time.time() - t0), p, ref


if __name__ == "__main__":
    res["L5"] = {}
    pats = {}
    PER = [int(v) for v in sys.argv[1:]] or [8]
    TAG = "" if PER == [8] else "_" + "_".join(map(str, PER))
    for per_lam in PER:
        r, p, ref = double_slit(per_lam)
        res["L5"][str(per_lam)] = r
        pats[per_lam] = p
        print(f"  lambda/a = {per_lam}: corr coarse {r['corr_coarse']:.4f} fine {r['corr_fine']:.4f}, bright-fringe err rms "
              f"{100 * r['fringe_rms']:.1f}% max {100 * r['fringe_max']:.1f}%  ({r['seconds']:.0f} s)", flush=True)
    r8 = res["L5"].get("8")
    if r8 is not None:
      check("L5a double slit, correlation with exact RS flux (lambda = 8a)", r8["corr_coarse"] >= TOL["corr"],
          f"{r8['corr_coarse']:.4f} (photon base rule: 0.9962)")
      check("L5b bright-fringe rms error (lambda = 8a)", r8["fringe_rms"] <= TOL["fringe_rms"], f"{100 * r8['fringe_rms']:.1f} % (photon base rule 4.3 %)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 3.4))
        ax.plot(XS, ref, "k-", lw=1.2, label="exact Rayleigh-Sommerfeld flux")
        for pl, c in ((6, "tab:blue"), (8, "r"), (12, "tab:orange"), (16, "tab:green")):
            if pl in pats:
                ax.plot(XS, pats[pl], "-", color=c, lw=0.8, label=f"light walk, lambda = {pl}a")
        ax.set_xlabel("x / lambda on the screen (z = 160 lambda)"); ax.set_ylabel("probability per bin")
        ax.legend(fontsize=7); ax.set_title("Young double slit made of spin-1 walkers")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, f"zz_slit{TAG}.png"), dpi=90)
    except Exception as ex:
        print("figure skipped:", ex)
    res["verdict"] = verdict
    res["seconds"] = time.time() - t_start
    json.dump(res, open(os.path.join(OUT, f"zz_slit{TAG}.json"), "w"), indent=1, default=float)
    npass = sum(v for _, v in verdict)
    print(f"=== {npass} PASS / {len(verdict) - npass} FAIL  ({time.time() - t_start:.0f} s)", flush=True)

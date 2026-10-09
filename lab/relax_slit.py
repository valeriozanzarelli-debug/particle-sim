"""relax_slit.py - Test 3 of relax_prereg: does a wrong start relax by itself in the double slit?

3a: 1 transverse dim (x), z = time, incident Gaussian -> double slit -> screen. Start right behind the
    slits either at |psi|^2 (equilibrium) or uniform over the open slits (wrong, like the old 0.89 case).
3b: 2 transverse dims (x,y). Upstream a square multimode guide [0,pi]^2 (relax_core.Box, M modes) of
    length Z0 where the surfers start uniform (wrong); at z = Z0 a double-slit mask in x; then free FFT
    propagation to the screen. Score = correlation of the x-histogram of clicks with the exact x-marginal.
usage: python lab/relax_slit.py 1d | 2d
"""
import os, sys, time, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from relax_core import Box, rk4_step, PI

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def corr_score(clicks, xs, refx, lo, hi, nfine=400, group=8):
    """histogram the clicks in nfine bins on [lo,hi], group by `group`, correlate with the reference mass"""
    edges = np.linspace(lo, hi, nfine + 1)
    h, _ = np.histogram(clicks, edges)
    idx = np.clip(np.searchsorted(edges, xs) - 1, 0, nfine - 1)
    inside = (xs >= lo) & (xs < hi)
    r = np.bincount(idx[inside], weights=refx[inside], minlength=nfine)
    hc = h.reshape(-1, group).sum(1); rc = r.reshape(-1, group).sum(1)
    return float(np.corrcoef(hc, rc)[0, 1]), r


def ideal_band(r, n, lo, hi, nfine=400, group=8, reps=60, seed=5):
    rng = np.random.default_rng(seed)
    p = r / r.sum()
    edges = np.linspace(lo, hi, nfine + 1)
    out = []
    for _ in range(reps):
        b = rng.choice(nfine, size=n, p=p)
        x = edges[b] + rng.uniform(0, edges[1] - edges[0], n)
        h, _ = np.histogram(x, edges)
        out.append(np.corrcoef(h.reshape(-1, group).sum(1), r.reshape(-1, group).sum(1))[0, 1])
    return float(np.mean(out)), float(np.std(out))


def smooth_slits(x, c1, c2, a, e):
    from scipy.special import erf
    s = lambda c: 0.5 * (erf((x - c + a / 2) / e) - erf((x - c - a / 2) / e))
    return np.clip(s(c1) + s(c2), 0, 1)


# ---------------------------------------------------------------- 3a: one transverse dimension
def test1d(N=6000):
    K = 10.0                                  # wavenumber, i dpsi/dz = -(1/2K) d2psi/dx2
    L, nx = 160.0, 8192
    x = (np.arange(nx) - nx // 2) * (L / nx)
    kx = 2 * PI * np.fft.fftfreq(nx, L / nx)
    mask = smooth_slits(x, -1.0, 1.0, 0.5, 0.03)
    psi0 = np.exp(-x ** 2 / (2 * 3.0 ** 2)) * mask
    psi0 /= np.sqrt((np.abs(psi0) ** 2).sum())
    F0 = np.fft.fft(psi0)

    def fld(z):
        Fz = F0 * np.exp(-1j * kx ** 2 * z / (2 * K))
        return np.fft.ifft(Fz), np.fft.ifft(1j * kx * Fz)

    def vel(xp, z):
        p, d = fld(z)
        v = np.imag(np.conj(p) * d) / (K * (np.abs(p) ** 2 + 1e-300))
        return np.interp(xp, x, v)

    rng = np.random.default_rng(3)
    rho = np.abs(psi0) ** 2
    starts = {
        "eq": rng.choice(x, size=N, p=rho / rho.sum()) + rng.uniform(-0.5, 0.5, N) * (L / nx),
        "wrong": rng.choice(x, size=N, p=(mask > 0.5) / (mask > 0.5).sum()) + rng.uniform(-0.5, 0.5, N) * (L / nx),
        # post-hoc (declared): a much more wrong start, rho0 = |psi|^4, and the KS distance as a finer score
        "wrong4": rng.choice(x, size=N, p=rho ** 2 / (rho ** 2).sum()) + rng.uniform(-0.5, 0.5, N) * (L / nx),
        # post-hoc: wrong slit balance, 2/3 of the surfers in the left slit (|psi|^2 shape inside each slit)
        "wrong_lr": rng.choice(x, size=N, p=rho * np.where(x < 0, 2, 1) / (rho * np.where(x < 0, 2, 1)).sum())
                    + rng.uniform(-0.5, 0.5, N) * (L / nx),
    }
    screens = [0.5, 2.0, 5.0, 10.0, 20.0]
    zs = np.concatenate([np.linspace(0, 0.5, 251), np.linspace(0.5, 20, 976)[1:]])
    res = {"screens": screens, "eq": [], "wrong": [], "wrong4": [], "wrong_lr": [], "band": [], "ks": {"eq": [], "wrong": [], "wrong4": [], "wrong_lr": []}}
    pos = {k: v.copy() for k, v in starts.items()}
    si = 0
    for i in range(len(zs)):
        z = zs[i]
        if si < len(screens) and abs(z - screens[si]) < 1e-9:
            p, _ = fld(z)
            refx = np.abs(p) ** 2
            # window holding 99.8 % of the flux
            c = np.cumsum(refx) / refx.sum()
            lo, hi = x[np.searchsorted(c, 0.001)], x[np.searchsorted(c, 0.999)]
            for k in pos:
                cs, r = corr_score(pos[k], x, refx, lo, hi)
                res[k].append(cs)
                cdf = np.cumsum(refx) / refx.sum()
                emp = np.searchsorted(np.sort(pos[k]), x) / N
                res["ks"][k].append(float(np.abs(cdf - emp).max()))
            res["band"].append(ideal_band(r, N, lo, hi))
            print(f"1D z={z:5.1f}: eq {res['eq'][-1]:.4f}  wrong {res['wrong'][-1]:.4f}  wrong4 {res['wrong4'][-1]:.4f} "
                  f"lr {res['wrong_lr'][-1]:.4f} KS eq/wrong/wrong4/lr {res['ks']['eq'][-1]:.4f}/{res['ks']['wrong'][-1]:.4f}/"
                  f"{res['ks']['wrong4'][-1]:.4f}/{res['ks']['wrong_lr'][-1]:.4f}  ideal "
                  f"{res['band'][-1][0]:.4f}+-{res['band'][-1][1]:.4f}", flush=True)
            if z == 20.0:
                res["screen_ref"] = refx[(x >= lo) & (x < hi)][::4].tolist()
                res["screen_x"] = x[(x >= lo) & (x < hi)][::4].tolist()
                res["screen_eq"] = pos["eq"].tolist(); res["screen_wrong"] = pos["wrong"].tolist(); res["screen_wrong4"] = pos["wrong4"].tolist(); res["screen_lr"] = pos["wrong_lr"].tolist()
            si += 1
        if i == len(zs) - 1:
            break
        h = zs[i + 1] - z
        for k in pos:
            xp = pos[k]
            k1 = vel(xp, z); k2 = vel(xp + 0.5 * h * k1, z + 0.5 * h)
            k3 = vel(xp + 0.5 * h * k2, z + 0.5 * h); k4 = vel(xp + h * k3, z + h)
            pos[k] = xp + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return res


# ---------------------------------------------------------------- 3b: two transverse dimensions
def downstream2d(box, Z0, xs0, ys0, rng, ZS=1.5, nsteps=150):
    nx, ny = 1024, 256
    Lx, Ly = 64.0, 32.0
    x = (np.arange(nx) - nx // 2) * (Lx / nx) + PI / 2
    y = (np.arange(ny) - ny // 2) * (Ly / ny) + PI / 2
    X, Y = np.meshgrid(x, y, indexing="ij")
    inbox = (X > 0) & (X < PI) & (Y > 0) & (Y < PI)
    psi = np.zeros((nx, ny), complex)
    psi[inbox] = box.field(X[inbox], Y[inbox], Z0, grad=False)
    m = smooth_slits(x, PI / 2 - 0.7, PI / 2 + 0.7, 0.8, 0.08)
    psi *= m[:, None]
    kx = 2 * PI * np.fft.fftfreq(nx, Lx / nx); ky = 2 * PI * np.fft.fftfreq(ny, Ly / ny)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    F0 = np.fft.fft2(psi)
    dx, dy = Lx / nx, Ly / ny

    cache = {}

    def vel(xp, yp, z):
        if cache.get("z") != z:
            Fz = F0 * np.exp(-0.5j * (KX ** 2 + KY ** 2) * z)
            p = np.fft.ifft2(Fz); px = np.fft.ifft2(1j * KX * Fz); py = np.fft.ifft2(1j * KY * Fz)
            r2 = np.abs(p) ** 2 + 1e-300
            cache.update(z=z, vx=np.imag(np.conj(p) * px) / r2, vy=np.imag(np.conj(p) * py) / r2)
        vx, vy = cache["vx"], cache["vy"]
        fi = np.clip((xp - x[0]) / dx, 0, nx - 1.001); fj = np.clip((yp - y[0]) / dy, 0, ny - 1.001)
        i, j = fi.astype(int), fj.astype(int); a, b = fi - i, fj - j
        def bil(f):
            return ((1 - a) * (1 - b) * f[i, j] + a * (1 - b) * f[i + 1, j]
                    + (1 - a) * b * f[i, j + 1] + a * b * f[i + 1, j + 1])
        return bil(vx), bil(vy)

    # survivors of the mask: probability m(x)^2 (equivariant absorption)
    keep = rng.uniform(size=xs0.size) < np.interp(xs0, x, m) ** 2
    xp, yp = xs0[keep], ys0[keep]
    # equilibrium start at the same plane, for the method's own ceiling
    rho = np.abs(psi) ** 2
    flat = rng.choice(nx * ny, size=xp.size, p=(rho / rho.sum()).ravel())
    xe = x[flat // ny] + rng.uniform(-0.5, 0.5, xp.size) * dx
    ye = y[flat % ny] + rng.uniform(-0.5, 0.5, xp.size) * dy
    zg = ZS * (np.linspace(0, 1, nsteps + 1) ** 2)       # small steps near the mask
    P = [[xp, yp], [xe, ye]]
    for s in range(nsteps):
        z, h = zg[s], zg[s + 1] - zg[s]
        for q in P:
            a1 = vel(q[0], q[1], z)
            a2 = vel(q[0] + 0.5 * h * a1[0], q[1] + 0.5 * h * a1[1], z + 0.5 * h)
            q[0] = q[0] + h * a2[0]; q[1] = q[1] + h * a2[1]
    Fz = F0 * np.exp(-0.5j * (KX ** 2 + KY ** 2) * ZS)
    refx = (np.abs(np.fft.ifft2(Fz)) ** 2).sum(1)
    c = np.cumsum(refx) / refx.sum()
    lo, hi = x[np.searchsorted(c, 0.001)], x[np.searchsorted(c, 0.999)]
    cw, r = corr_score(P[0][0], x, refx, lo, hi, nfine=200, group=4)
    ce, _ = corr_score(P[1][0], x, refx, lo, hi, nfine=200, group=4)
    band = ideal_band(r, xp.size, lo, hi, nfine=200, group=4)
    return dict(wrong=cw, eq=ce, band=band, n=int(xp.size),
                x=x[(x >= lo) & (x < hi)].tolist(), ref=refx[(x >= lo) & (x < hi)].tolist(),
                clicks_wrong=P[0][0].tolist())


def test2d(N=10000):
    res = {}
    for M in (16, 1):
        box = Box(M, seed=1)
        rng = np.random.default_rng(7)
        xs, ys = rng.uniform(0, PI, N), rng.uniform(0, PI, N)     # wrong start in the guide
        t = 0.0; t0 = time.time()
        for Z0 in (0.0, PI, 4 * PI):
            while t < Z0 - 1e-9:
                xs, ys = rk4_step(box, xs, ys, t, 0.02, ds=0.02); t += 0.02
            r = downstream2d(box, Z0, xs, ys, np.random.default_rng(11))
            res[f"M{M}_Z{Z0:.2f}"] = r
            print(f"2D M={M} Z0={Z0:5.2f}: wrong {r['wrong']:.4f}  eq {r['eq']:.4f}  ideal {r['band'][0]:.4f}"
                  f"+-{r['band'][1]:.4f}  clicks {r['n']}  [{time.time()-t0:.0f}s]", flush=True)
    return res


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    which = sys.argv[1]
    t0 = time.time()
    res = test1d() if which == "1d" else test2d()
    res["secs"] = time.time() - t0
    with open(os.path.join(OUT, f"relax_slit_{which}.json"), "w") as f:
        json.dump(res, f)
    print(f"done {res['secs']:.0f}s")

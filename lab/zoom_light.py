"""ZOOM test c: bundles of photon balls behind a slit -> does the agglomerator find the "ray" level of light
(docs/MULTIFOTONE.md section 6.3) by itself?

Fine level (brute force): a slit of half-width a lit by a plane wave, lambda = 1, screen at distance L,
Fresnel number N_F = a^2 / (lambda L).  The wave is the exact Fresnel solution (the continuum limit of
the wave balls, which reproduce it to 0.02% in MULTIFOTONE section 6.2); N0 = 4096 photon balls start
evenly in the slit and are guided by it.  In 1-D the guided balls never cross, so photon i lands on the
quantile (i + 1/2)/N0 of the screen intensity and owns the light between quantiles i/N0 and (i+1)/N0.

Agglomeration (zoom_core, same code as for springs and spins):
  block(k=2) of neighbouring photon balls -> super-ball with count n, centre X, internal 2nd moment s2.
  learn('transfer'): change of (X, s2) from slit to screen = generic local model (pair terms with the
  neighbouring super-balls + bias), learned on OTHER slits (half-widths 1.5a, 2a, 3a, 4a, same L, same
  ball spacing), never on the test slit.  Level n+1 is learned from the blocked PREDICTIONS of level n
  (true repetition: the fine level is used only for level 1).
  coarse picture at level n: each super-ball deposits its n photons uniformly on X +- sqrt(3 s2).  At the
  top level (one super-ball holding all photons) this is ONE RAY BUNDLE with a photon count.
  check: total-variation distance from the exact wave, over the whole screen and in random small
  windows (spot-checks, refinement flag at 10%).
Nothing in the code knows about Fresnel numbers: the error vs N_F is measured afterwards and compared
with the documented criterion (rays good to 1% for N_F > 4400, 10% for N_F > 13).
Run: timeout 590 python lab/zoom_light.py
"""
import os, sys, time, json
import numpy as np
from scipy.special import fresnel

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zoom_core as zc

T0 = time.time()
OUT = os.path.join(HERE, "out")
A = 1.0
N0 = 4096
H = 2 * A / N0
NF_LIST = 2.0 ** np.arange(-2, 15)            # 0.25 ... 16384
TRAIN = [1.5, 2.0, 3.0, 4.0]                  # training slits (half-width / a)
NG = 800001
rng = np.random.default_rng(31)


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


def fresnel_slit(x, a, Lz):
    """Exact paraxial (Fresnel) intensity behind a 2-D slit |x|<a lit by a unit plane wave (lambda=1)."""
    s = np.sqrt(2 / Lz)
    S1, C1 = fresnel(s * (-a - x)); S2, C2 = fresnel(s * (a - x))
    return 0.5 * ((C2 - C1) ** 2 + (S2 - S1) ** 2)


def screen(a, Lz):
    half = a + 12 * np.sqrt(Lz) + 2 * Lz / a
    x = np.linspace(-half, half, NG)
    I = fresnel_slit(x, a, Lz)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (I[1:] + I[:-1]) * np.diff(x))])
    return x, cdf / cdf[-1]


def fine(a, Lz):
    """Fine level: photon balls in the slit and where they land (non-crossing guided balls = quantiles)."""
    n = int(round(2 * a / H))
    x0 = -a + (np.arange(n) + 0.5) * H
    x, cdf = screen(a, Lz)
    edges = np.interp(np.arange(n + 1) / n, cdf, x)
    xo = 0.5 * (edges[1:] + edges[:-1])
    cell = np.diff(edges)
    st_in = zc.State(x0, np.ones(n), H, ref=x0, s2=np.full(n, H * H / 12))
    st_out = zc.State(xo, np.ones(n), H, ref=x0, s2=cell ** 2 / 12)
    return st_in, st_out, (x, cdf)


def coarse_cdf(X, s2, nI, x):
    """CDF on the grid x of uniform boxes X +- sqrt(3 s2) carrying nI photons each."""
    r = np.sqrt(3 * np.maximum(s2, 0)) + 1e-12
    lo, hi = X - r, X + r
    bp = np.concatenate([lo, hi])
    sl = np.concatenate([nI / (2 * r), -nI / (2 * r)])
    o = np.argsort(bp)
    bp, sl = bp[o], sl[o]
    dens = np.cumsum(sl)[:-1]
    c = np.concatenate([[0], np.cumsum(dens * np.diff(bp))])
    return np.interp(x, bp, c, left=0, right=c[-1]) / nI.sum()


def tv(cdf_a, cdf_b):
    return 0.5 * np.abs(np.diff(cdf_a) - np.diff(cdf_b)).sum()


def main():
    log(f"slit a={A}, lambda=1, N0={N0} photon balls (spacing {H:.5f}), N_F = {NF_LIST[0]} ... {NF_LIST[-1]:.0f}; training slits {TRAIN}")
    nlev = int(np.log2(N0))
    table = []
    for NF in NF_LIST:
        Lz = A * A / NF
        t0 = time.time()
        test_in, test_out, (x, cdf) = fine(A, Lz)
        trains = [fine(A * f, Lz)[:2] for f in TRAIN]
        e_fine = tv(cdf, coarse_cdf(test_out.pos[:, 0], test_out.s2, test_out.n, x))
        rows = []
        tin, tout = test_in, test_out
        tr_in = [t[0] for t in trains]; tr_tgt = [t[1] for t in trains]      # level-n training targets
        dir_in = [t[0] for t in trains]; dir_out = [t[1] for t in trains]    # direct (always from fine)
        for lev in range(1, nlev + 1):
            # block test and training states with the same generic code
            tin, B = zc.block(tin, 2); tout, _ = zc.block(tout, 2, blocking=B)
            new_in, new_tgt, nd_in, nd_out = [], [], [], []
            for si, so in zip(tr_in, tr_tgt):
                ci, Bt = zc.block(si, 2); co, _ = zc.block(so, 2, blocking=Bt)
                new_in.append(ci); new_tgt.append(co)
            for si, so in zip(dir_in, dir_out):
                ci, Bt = zc.block(si, 2); co, _ = zc.block(so, 2, blocking=Bt)
                nd_in.append(ci); nd_out.append(co)
            tr_in, tr_tgt, dir_in, dir_out = new_in, new_tgt, nd_in, nd_out
            m = zc.LocalModel()
            m.learn_transfer(tr_in, [o.pos - i.pos for i, o in zip(tr_in, tr_tgt)], [o.s2 - i.s2 for i, o in zip(tr_in, tr_tgt)])
            # predictions of this level on the training slits become the targets of the next level
            tr_tgt = []
            for si in tr_in:
                dX, ds2 = m.transfer(si)
                tr_tgt.append(si.copy(pos=si.pos + dX, s2=np.maximum(si.s2 + ds2, 0)))
            dX, ds2 = m.transfer(tin)
            Xp = tin.pos[:, 0] + dX[:, 0]; s2p = np.maximum(tin.s2 + ds2, 0)
            e_l = tv(cdf, coarse_cdf(Xp, s2p, tin.n, x))
            e_ray = tv(cdf, coarse_cdf(tin.pos[:, 0], tin.s2, tin.n, x))        # super-ball goes straight
            e_best = tv(cdf, coarse_cdf(tout.pos[:, 0], tout.s2, tout.n, x))   # true descriptors (no rule error)
            md = zc.LocalModel()
            md.learn_transfer(dir_in, [o.pos - i.pos for i, o in zip(dir_in, dir_out)], [o.s2 - i.s2 for i, o in zip(dir_in, dir_out)])
            dXd, ds2d = md.transfer(tin)
            e_dir = tv(cdf, coarse_cdf(tin.pos[:, 0] + dXd[:, 0], np.maximum(tin.s2 + ds2d, 0), tin.n, x))
            # spot checks: random windows 4 super-balls wide, error relative to the light inside
            ccdf = coarse_cdf(Xp, s2p, tin.n, x)
            spots = []
            for _ in range(3):
                c0 = rng.uniform(-A, A)
                w = 2 * tin.h
                sel = (x >= c0 - w) & (x <= c0 + w)
                ii = np.where(sel)[0]
                def fine_fn(R):
                    return np.diff(cdf[R[0]:R[-1] + 1])
                def coarse_fn(R):
                    return np.diff(ccdf[R[0]:R[-1] + 1])
                ck = zc.check(ii, fine_fn, coarse_fn, lambda f, c: 0.5 * np.abs(f - c).sum() / max(f.sum(), 1e-300), 0.10)
                spots.append(ck)
            nbf = (tin.h / 2) ** 2 / Lz
            rows.append(dict(level=lev, balls=tin.N, width=tin.h, NF_block=nbf, err=e_l, err_ray=e_ray,
                             err_best=e_best, err_direct=e_dir, spot=[s["err"] for s in spots],
                             spot_flag=[s["flagged"] for s in spots],
                             dX_edge=float(np.abs(dX).max() / tin.h), ds2_bias=float(m.theta_s[0])))
        top = rows[-1]
        r = dict(NF=float(NF), L=Lz, err_fine=e_fine, levels=rows, t=time.time() - t0)
        for tol in (0.01, 0.10):
            ok = [x_["level"] for x_ in rows if x_["err"] <= tol]
            r[f"level_{int(tol * 100)}"] = max(ok) if ok else (0 if e_fine <= tol else -1)
        table.append(r)
        log(f"N_F={NF:9.2f}: fine err {e_fine * 100:.3f}% | level 2 ({rows[1]['balls']} balls) {rows[1]['err'] * 100:.2f}% (ray {rows[1]['err_ray'] * 100:.2f}%) | "
            f"level 6 ({rows[5]['balls']}) {rows[5]['err'] * 100:.2f}% (ray {rows[5]['err_ray'] * 100:.2f}%) | TOP = 1 bundle: learned {top['err'] * 100:.2f}%, "
            f"pure ray {top['err_ray'] * 100:.2f}%, direct {top['err_direct'] * 100:.2f}% | coarsest level ok @1%: {r['level_1']}, @10%: {r['level_10']} ({time.time() - t0:.1f}s)")
    # documented curve for comparison (MULTIFOTONE 6.3)
    doc = lambda nf: (0.231 + 0.051 * np.log(nf)) / np.sqrt(nf)
    log("one-bundle (top level) error vs N_F, against MULTIFOTONE 6.3 fit (0.231 + 0.051 ln N_F)/sqrt(N_F):")
    for r in table:
        t = r["levels"][-1]
        log(f"   N_F={r['NF']:9.2f}  learned {t['err'] * 100:7.3f}%   pure ray {t['err_ray'] * 100:7.3f}%   doc fit {doc(r['NF']) * 100:7.3f}%   learned width factor {np.exp(t['ds2_bias'] / 2):.3f}   spot-checks {', '.join(f'{s * 100:.1f}%' for s in t['spot'])}")
    # thresholds by log-interpolation of the top-level error
    nf = np.array([r["NF"] for r in table])
    thr = {}
    for key in ("err", "err_ray"):
        e = np.array([r["levels"][-1][key] for r in table])
        for tol in (0.10, 0.01):
            k = np.where((e[:-1] > tol) & (e[1:] <= tol))[0]
            if len(k):
                k = k[-1]
                f = (np.log(e[k]) - np.log(tol)) / (np.log(e[k]) - np.log(e[k + 1]))
                thr[f"{key}_{tol}"] = float(np.exp(np.log(nf[k]) + f * (np.log(nf[k + 1]) - np.log(nf[k]))))
            else:
                thr[f"{key}_{tol}"] = None
    log(f"N_F above which ONE bundle (ray) is good: learned rule  10%: {thr['err_0.1']:.1f}, 1%: {thr['err_0.01']}"
        f" | pure ray 10%: {thr['err_ray_0.1']:.1f}, 1%: {thr['err_ray_0.01']}   (documented: 13 and 4400)")
    # error vs Fresnel number of the super-ball itself (all levels, all N_F)
    log("error of level n vs the super-ball's own Fresnel number N_F,block = (w/2)^2/(lambda L), w = width of one super-ball:")
    bins = 2.0 ** np.arange(-12, 16, 2)
    allp = [(l["NF_block"], l["err"], l["level"]) for r in table for l in r["levels"]]
    for b0, b1 in zip(bins[:-1], bins[1:]):
        e = [p[1] for p in allp if b0 <= p[0] < b1]
        if e:
            log(f"   N_F,block in [{b0:g},{b1:g}): {len(e):3d} cases, error median {np.median(e) * 100:6.2f}%  max {np.max(e) * 100:6.2f}%")
    # cost table: balls needed for 1% and 10%
    log(f"balls needed (coarsest level whose global error <= tol); fine level = {N0} photon balls + the wave:")
    for r in table:
        l1, l10 = r["level_1"], r["level_10"]
        f1 = f"level {l1} ({N0 >> l1} balls)" if l1 >= 0 else "none (even the fine level fails)"
        log(f"   N_F={r['NF']:9.2f}: 1% -> {f1}, 10% -> level {l10} ({N0 >> max(l10, 0)} balls)")
    json.dump(dict(table=table, thresholds=thr), open(os.path.join(OUT, "zoom_light.json"), "w"), indent=1)
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), dpi=90)
        for lev, col in ((2, "C0"), (6, "C1"), (9, "C2"), (nlev, "k")):
            ax[0].loglog(nf, [100 * r["levels"][lev - 1]["err"] for r in table], "o-", color=col, ms=3,
                         label=f"livello {lev} ({N0 >> lev} super-palline)" if lev < nlev else f"livello {nlev} = 1 raggio")
        ax[0].loglog(nf, 100 * doc(nf), "k:", label="MULTIFOTONE 6.3")
        ax[0].axhline(10, color="gray", lw=0.5); ax[0].axhline(1, color="gray", lw=0.5)
        ax[0].set_xlabel("numero di Fresnel N_F"); ax[0].set_ylabel("errore vs onda esatta (%)"); ax[0].legend(fontsize=6)
        ax[1].loglog([p[0] for p in allp], [100 * p[1] for p in allp], ".", ms=3)
        ax[1].set_xlabel("N_F della super-pallina"); ax[1].set_ylabel("errore (%)")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "zoom_light.png")); plt.close(fig)
    except Exception as e:
        log("figure failed:", e)
    log("done")


if __name__ == "__main__":
    main()

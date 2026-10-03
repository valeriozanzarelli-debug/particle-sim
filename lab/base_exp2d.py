"""2-D experiments with the frozen base rule (lab/base_rule.py).  Each one prints measured / theory /
tolerance / PASS-FAIL / photons / clicks / wave balls / runtime and saves lab/out/base_<name>_<ver>.json.

  python lab/base_exp2d.py free | double | which | single | grating | talbot | edge
  ('double' also does the Taylor build-up, experiment 9)

The tolerances are written here BEFORE running (TOL) and are never changed to make a test pass.
"""
import sys
import time
import numpy as np
from base_rule import *
import optics as OP

TOL = dict(
    width_rel=0.02,          # free beam: |w_meas/w_th − 1| ≤ max(3σ, 2%)
    ds_corr=0.98,            # double slit: coarse correlation with the exact flux ≥ 0.98
    ds_rms=0.08,             # double slit: rms error of the bright-fringe populations ≤ 8%
    ds_absorbed=0.02,        # absorbed fraction within ±2 points of the exact RS model value
    vis_par=0.85, vis_orth=0.15,
    pos_rel=0.02,            # positions of minima / maxima / zeros: max(3σ, 2%)
    grating_pos_rel=0.01,    # principal maxima: max(3σ, 1%)
    talbot_rel=0.02,         # Talbot distance: max(plane spacing, 2%)
    edge_abs=0.02,           # Fresnel edge: I(0)/I0 = 1/4 within max(3σ, 0.02)
)


def bootstrap(fn, x, nb=200, seed=3):
    rng = np.random.default_rng(seed)
    v = [fn(x[rng.integers(0, len(x), len(x))]) for _ in range(nb)]
    return float(np.nanstd(v))


def rs_screen(s, B, x, z):
    """Exact 2-D RS field of the boundary values B(s) (fine samples), flux j_z on the points x at z."""
    w = B * np.gradient(s)[:, None] if B.ndim == 2 else (B * np.gradient(s))[:, None]
    U, Ux, Uz = OP.rs_field(x, z, s, w.astype(complex))
    return OP.flux(U, Ux, Uz)


def report(rows, name, meas, th, tol, ok, extra=""):
    rows.append(result_line(name, meas, th, tol, ok, extra))
    log(f"  RESULT {name}: measured {meas}  theory {th}  tol {tol}  -> {'PASS' if ok else 'FAIL'}  {extra}")


# ─────────────────────────────── 1. free Gaussian beam (2-D) ─────────────────────────────────
def exp_free(n=20000):
    w0 = 3.0
    zr = np.pi * w0 ** 2
    sc = Scene2D()
    src = sc.add(Obj2("source", "source", (0, 0), (1, 0), 40.0, src=dict(w0=w0, zw=0.0)))
    zs = [zr, 2 * zr, 5 * zr]
    for z in zs:
        sc.add(Obj2(f"plane z={z:.1f}", "plane", (0, z), (1, 0), 200.0))
    sc.add(Obj2("screen", "screen", (0, 5 * zr + 1), (1, 0), 200.0))
    log("free 2-D:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-100, 100, -1, 200], seed=11)
    rows = []
    for k, z in enumerate(zs):
        x = np.array([p[0] for i, pl, p in res["rec"] if pl == k + 1])
        wm = 2 * x.std()
        wth = w0 * np.sqrt(1 + (z / zr) ** 2)
        sig = 1 / np.sqrt(2 * len(x))
        err = wm / wth - 1
        report(rows, f"2-D free beam width at z={z / zr:.0f} zR", round(wm, 4), round(wth, 4),
               f"max(3σ={3 * sig:.3f}, {TOL['width_rel']})", abs(err) <= max(3 * sig, TOL["width_rel"]),
               f"err {100 * err:+.2f}%  n={len(x)}")
    meta = dict(photons=res["n"], clicks=res["counts"]["click"], counts=res["counts"], balls=balls_summary(sc),
                runtime=res["time"] + sc.build_time)
    save_results("free2d", rows, meta)
    return rows, meta


# ─────────────────────────────── 2. double slit + 9. build-up + 10. which-path ───────────────
ZS = 10.0


def ds_scene(jones=(1, 0), filters=None, rule=RULE):
    """Gaussian beam w0=7 with its waist on the wall (z=0), slits 2 wide at ±4, screen at 160 (as REGOLA_FOTONE)."""
    sc = Scene2D(rule)
    src = sc.add(Obj2("source", "source", (0, -ZS), (1, 0), 40.0, src=dict(w0=OP.W0, zw=ZS, jones=jones)))
    if filters is None:
        sc.add(Obj2("wall", "aperture", (0, 0), (1, 0), 400.0, ints=[(c - OP.A / 2, c + OP.A / 2) for c in OP.SLITS]))
    else:          # a polarizer on each slit: the wall is two objects (left half / right half)
        for k, (c, f) in enumerate(zip(OP.SLITS, filters)):
            f = np.asarray(f, complex)
            T = np.outer(f, f.conj())
            lo, hi = (-400.0, 0.0) if k == 0 else (0.0, 400.0)
            cc = 0.5 * (lo + hi)
            sc.add(Obj2(f"wall{k}+pol", "aperture", (cc, 0), (1, 0), 0.5 * (hi - lo),
                        ints=[(c - OP.A / 2 - cc, c + OP.A / 2 - cc)], T=T))
    sc.add(Obj2("screen", "screen", (0, OP.L), (1, 0), 400.0))
    return sc, src


def run_double(n, seed, jones=(1, 0), filters=None, rule=RULE, verbose=True):
    sc, src = ds_scene(jones, filters, rule)
    res = simulate(sc, src, n, box=[-200, 200, -ZS - 1, OP.L + 1], seed=seed, verbose=verbose)
    cl = res["clicks"][res["status"] == CLICK][:, 0]
    return sc, res, cl


def exp_double(n=27000):
    rows = []
    sc, res, xl = run_double(n, seed=7)
    np.save(os.path.join(OUT, "base_double_clicks.npy"), xl)
    ref = OP.screen_reference("flux")
    c = OP.coarse_corr(xl, ref)
    e, rms, mx = OP.fringe_errors(xl, ref)
    chi = OP.chi2(xl, ref)
    base = OP.ideal_baseline(ref, len(xl), 200)
    absorbed = 1 - len(xl) / res["n"]
    # RS-model absorbed fraction: 1 − (flux through the screen plane of the exact RS field)/(incident flux)
    s0, w0 = OP.sources(800)
    xx = np.linspace(-300, 300, 60001)
    jz = np.concatenate([OP.flux(*OP.rs_field(q, 2.0, s0, w0))[1] for q in np.array_split(xx, 30)])
    passed = jz.sum() * (xx[1] - xx[0]) / (K * OP.W0 * np.sqrt(np.pi / 2))
    log(f"double slit: {len(xl)} clicks of {res['n']} photons; corr {c:.4f} (ideal {base[0]:.4f}±{base[1]:.4f}), chi2/dof {chi:.2f}, "
        f"fringe rms {100 * rms:.1f}% max {100 * mx:.1f}% (ideal rms {100 * base[2]:.1f}%), per fringe {np.round(100 * e).astype(int).tolist()}")
    report(rows, "Young double slit: correlation with exact flux", round(c, 4), f"1 (ideal sampler {base[0]:.4f}±{base[1]:.4f})",
           f">= {TOL['ds_corr']}", c >= TOL["ds_corr"], f"REGOLA_FOTONE: 0.9949")
    report(rows, "Young double slit: bright-fringe population error rms / max", f"{100 * rms:.1f}% / {100 * mx:.1f}%",
           f"0 (ideal sampler rms {100 * base[2]:.1f}%)", f"rms <= {100 * TOL['ds_rms']:.0f}%", rms <= TOL["ds_rms"],
           "REGOLA_FOTONE: 3.9% / 4.7%")
    report(rows, "Young double slit: absorbed fraction", round(absorbed, 4), round(1 - passed, 4),
           f"±{TOL['ds_absorbed']}", abs(absorbed - (1 - passed)) <= TOL["ds_absorbed"], "REGOLA_FOTONE: 0.789")
    np.save(os.path.join(OUT, "base_double_clicks.npy"), xl)
    meta = dict(photons=res["n"], clicks=len(xl), counts=res["counts"], balls=balls_summary(sc),
                runtime=res["time"] + sc.build_time, chi2=chi, per_fringe=e.tolist())
    # 9. build-up: photons arrive one at a time (they are independent): the clicks in arrival order
    log("build-up (Taylor 1909):")
    rng = np.random.default_rng(5)
    order = rng.permutation(len(xl))
    ok_all = True
    txt = []
    for m in (10, 50, 200, 1000, 3000, len(xl)):
        cm = OP.coarse_corr(xl[order[:m]], ref)
        b = OP.ideal_baseline(ref, m, 200)
        good = cm >= b[0] - 3 * b[1]
        ok_all &= good
        txt.append(f"N={m}: {cm:.3f} (ideal {b[0]:.3f}±{b[1]:.3f})")
        log(f"  after {m:5d} clicks: corr {cm:.3f}   ideal sampler {b[0]:.3f} ± {b[1]:.3f}  {'ok' if good else 'LOW'}")
    report(rows, "Taylor build-up: correlation vs number of clicks", "; ".join(txt), "ideal sampler",
           "within 3σ of the ideal sampler at every N", ok_all)
    save_results("double", rows, meta)
    return rows, meta


def exp_which(n=30000):
    rows = []
    H, V = np.array([1, 0], complex), np.array([0, 1], complex)
    D45 = np.array([1, 1], complex) / np.sqrt(2)
    metas = {}
    for name, f in (("parallel", (H, H)), ("orthogonal", (H, V))):
        sc, res, xl = run_double(n, seed=13 if name == "parallel" else 17, jones=D45, filters=f)
        jones = [fi * np.vdot(fi, D45) for fi in f]
        ref = OP.screen_reference("flux", jones=jones)
        V_ = OP.visibility(xl)
        Vth = OP.visibility(ref, True)
        sv = bootstrap(OP.visibility, xl)
        c = OP.coarse_corr(xl, ref)
        if name == "parallel":
            ok = V_ > TOL["vis_par"]
            tol = f"> {TOL['vis_par']}"
        else:
            ok = V_ < TOL["vis_orth"]
            tol = f"< {TOL['vis_orth']}"
        report(rows, f"which-path, {name} polarizers: fringe visibility", f"{V_:.3f} ± {sv:.3f}", round(Vth, 3), tol, ok,
               f"corr with exact {c:.4f}; {len(xl)} clicks")
        metas[name] = dict(photons=res["n"], clicks=len(xl), counts=res["counts"], balls=balls_summary(sc),
                           runtime=res["time"] + sc.build_time)
    save_results("which", rows, metas)
    return rows, metas


# ─────────────────────────────── 3. single slit ─────────────────────────────────────────────
def min_position(x, guess, half, maximum=False, norm=None):
    """Position of a local minimum (or maximum) of the click density: weighted cubic fit of the histogram
    (30 bins in guess±half), stationary point nearest to the guess.  The estimator is checked on ideal
    samples of the exact pattern (ideal_check) because a symmetric parabola is biased on a skew minimum."""
    sel = x[np.abs(x - guess) < half]
    if len(sel) < 30:
        return np.nan
    h, e = np.histogram(sel, np.linspace(guess - half, guess + half, 31))
    c = 0.5 * (e[1:] + e[:-1]) - guess
    nn = np.ones_like(c) if norm is None else norm(c + guess)
    p = np.polyfit(c, h / nn, 3, w=nn / np.sqrt(np.maximum(h, 1)))
    r = np.roots(np.polyder(p))
    r = r[np.isreal(r)].real
    sgn = -1 if maximum else 1
    r = [q for q in r if abs(q) < half and sgn * np.polyval(np.polyder(p, 2), q) > 0]
    return guess + min(r, key=abs) if r else np.nan


def ideal_samples(xx, p, n, seed):
    """n samples of the exact click density p(x) on the grid xx (the 'ideal sampler')."""
    p = np.maximum(p, 0)
    cdf = np.cumsum(p) / p.sum()
    return np.interp(np.random.default_rng(seed).random(n), cdf, xx)


def ideal_check(est, xx, p, n, reps=40):
    v = [est(ideal_samples(xx, p, n, 100 + k)) for k in range(reps)]
    return float(np.nanmean(v)), float(np.nanstd(v))


def exp_single(n=60000):
    a, Ls = 4.0, 400.0
    sc = Scene2D()
    src = sc.add(Obj2("slit", "source", (0, 0), (1, 0), 20.0, ints=[(-a / 2, a / 2)], src=dict(w0=np.inf)))
    sc.add(Obj2("screen", "screen", (0, Ls), (1, 0), 1000.0))
    log("single slit:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-600, 600, -1, Ls + 1], seed=21)
    x = res["clicks"][res["status"] == CLICK][:, 0]
    np.save(os.path.join(OUT, "base_single_clicks.npy"), x)
    rows = []
    # exact RS pattern (same boundary values, fine quadrature): for the ideal-sampler check of the estimator
    s = np.linspace(-a / 2, a / 2, 4001)
    s = 0.5 * (s[1:] + s[:-1])
    xx = np.linspace(-550, 550, 44001)
    jz = np.concatenate([rs_screen(s, np.ones_like(s), q, Ls)[1] for q in np.array_split(xx, 40)])
    for m in (1, 2):
        th = np.arcsin(m * LAM / a)
        xf = Ls * np.tan(th)
        win = 0.3 * LAM * Ls / a
        sel = np.abs(xx - xf) < 15
        xrs = xx[sel][np.argmin(jz[sel])]
        est = lambda q: 0.5 * (min_position(q, xf, win) + min_position(-q, xf, win))
        xm = est(x)
        sig = bootstrap(est, x)
        im, isd = ideal_check(est, xx, jz, len(x))
        err = xm / xf - 1
        report(rows, f"single slit a={a}: minimum m={m} (mean of ±)", f"{xm:.2f} ± {sig:.2f}",
               f"{xf:.2f} = L·tan(asin(mλ/a)) [paraxial mλL/a = {m * LAM * Ls / a:.1f}; exact RS {xrs:.2f}]",
               f"max(3σ, {100 * TOL['pos_rel']:.0f}%)", abs(xm - xf) <= max(3 * sig, TOL["pos_rel"] * xf),
               f"err {100 * err:+.2f}%; same estimator on ideal samples of the exact pattern: {im:.2f} ± {isd:.2f}")
    meta = dict(photons=res["n"], clicks=len(x), counts=res["counts"], balls=balls_summary(sc),
                runtime=res["time"] + sc.build_time)
    save_results("single", rows, meta)
    return rows, meta


# ─────────────────────────────── 4. N-slit grating ──────────────────────────────────────────
def exp_grating(n=40000):
    N, d, a, Lg = 5, 8.0, 2.0, 200.0
    wg = 60.0
    w0f = LAM * Lg / (np.pi * wg)            # converging Gaussian: focus (Fraunhofer plane) on the screen
    cs = (np.arange(N) - (N - 1) / 2) * d
    sc = Scene2D()
    src = sc.add(Obj2("grating", "source", (0, 0), (1, 0), 60.0, ints=[(c - a / 2, c + a / 2) for c in cs],
                      src=dict(w0=w0f, zw=Lg)))
    sc.add(Obj2("screen", "screen", (0, Lg), (1, 0), 1000.0))
    log("grating:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-600, 600, -1, Lg + 1], seed=31)
    x = res["clicks"][res["status"] == CLICK][:, 0]
    rows = []
    for m in (1, 2):
        th = np.arcsin(m * LAM / d)
        xf = Lg * np.tan(th)
        half = 0.5 * LAM * Lg / (N * d)
        cen = lambda q, s_: q[np.abs(q - s_ * xf) < half].mean() * s_
        xm = 0.5 * (cen(x, 1) + cen(x, -1))
        sig = bootstrap(lambda q: 0.5 * (cen(q, 1) + cen(q, -1)), x)
        report(rows, f"grating N={N}, d={d}: principal maximum m={m} (centroid, mean of ±)", f"{xm:.3f} ± {sig:.3f}",
               f"{xf:.3f} = L·tan(asin(mλ/d)) [paraxial {m * LAM * Lg / d:.2f}]",
               f"max(3σ, {100 * TOL['grating_pos_rel']:.0f}%)", abs(xm - xf) <= max(3 * sig, TOL["grating_pos_rel"] * xf),
               f"err {100 * (xm / xf - 1):+.2f}%")
    # secondary maxima between the orders 0 and ±1: significant local maxima of the histogram
    P1 = Lg * np.tan(np.arcsin(LAM / d))
    step = LAM * Lg / (N * d)

    def count_peaks(q, side):
        e = np.arange(0.5 * step, P1 - 0.5 * step + 1e-9, step / 5)
        h, _ = np.histogram(side * q, e)
        cnt = 0
        for i in range(1, len(h) - 1):
            if not (h[i] >= h[i - 1] and h[i] > h[i + 1]):
                continue
            # prominence: drop to the lowest point before meeting a higher bin, on both sides
            l = i - 1
            lmin = h[i]
            while l >= 0 and h[l] <= h[i]:
                lmin = min(lmin, h[l]); l -= 1
            r_ = i + 1
            rmin = h[i]
            while r_ < len(h) and h[r_] <= h[i]:
                rmin = min(rmin, h[r_]); r_ += 1
            if h[i] - max(lmin, rmin) > 3 * np.sqrt(h[i]):
                cnt += 1
        return cnt
    found = [count_peaks(x, +1), count_peaks(x, -1)]
    # the same counter on ideal samples of the exact pattern (exact RS of the same boundary values)
    s_ = np.concatenate([np.linspace(c - a / 2, c + a / 2, 401)[:-1] + a / 800 for c in cs])
    B_ = sc.objs[src].src_B(s_)[:, 0]
    xx = np.linspace(-3 * P1, 3 * P1, 6001)
    jz = np.concatenate([rs_screen(s_, B_ * (a / 400) / np.gradient(s_), q, Lg)[1] for q in np.array_split(xx, 12)])
    ideal = [count_peaks(ideal_samples(xx, jz, len(x), 200 + k), +1) for k in range(20)]
    report(rows, f"grating N={N}: secondary maxima between order 0 and ±1", f"{found[0]} (+), {found[1]} (−)", N - 2,
           "exact", all(f == N - 2 for f in found),
           f"zeros every λL/(Nd) = {step:.1f}; same counter on 20 ideal samples of the exact pattern: {np.mean(ideal):.2f} (min {min(ideal)})")
    meta = dict(photons=res["n"], clicks=len(x), counts=res["counts"], balls=balls_summary(sc),
                runtime=res["time"] + sc.build_time)
    np.save(os.path.join(OUT, "base_grating_clicks.npy"), x)
    save_results("grating", rows, meta)
    return rows, meta


# ─────────────────────────────── 7. Talbot ───────────────────────────────────────────────────
def exp_talbot(n=20000):
    d, a, w = 8.0, 2.0, 40.0
    zt = 2 * d ** 2 / LAM
    cs = np.arange(-15, 16) * d
    sc = Scene2D()
    src = sc.add(Obj2("grating", "source", (0, 0), (1, 0), 130.0, ints=[(c - a / 2, c + a / 2) for c in cs],
                      src=dict(w0=w, zw=0.0)))
    zp = zt * np.arange(0.40, 1.2001, 0.01)
    for z in zp:
        sc.add(Obj2(f"plane {z:.1f}", "plane", (0, z), (1, 0), 400.0))
    sc.add(Obj2("screen", "screen", (0, zp[-1] + 1), (1, 0), 400.0))
    log("talbot:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-400, 400, -1, zp[-1] + 2], seed=41)
    rec = res["rec"]
    pl = np.array([r[1] for r in rec])
    xs = np.array([r[2][0] for r in rec])
    half = 2.0 * d                                   # central periods only (uniform part of the beam)

    def contrast(x):
        x = x[np.abs(x) < half]
        ph = np.mod(x + d / 2, d) - d / 2              # position relative to the nearest slit centre
        return np.mean(np.abs(ph) < a / 2) if len(x) else np.nan
    C = np.array([contrast(xs[pl == i + 1]) for i in range(len(zp))])
    # exact RS reference of the same quantity on the same planes
    s = np.concatenate([np.linspace(c - a / 2, c + a / 2, 81)[:-1] + a / 160 for c in cs])
    B = np.exp(-s ** 2 / w ** 2)
    Cr = []
    xx = np.linspace(-half, half, 801)
    for z in zp:
        jz = rs_screen(s, B * (a / 80) / np.gradient(s), xx, z)[1]
        ph = np.mod(xx + d / 2, d) - d / 2
        Cr.append(jz[np.abs(ph) < a / 2].sum() / jz.sum())
    Cr = np.array(Cr)
    near = zp > 0.8 * zt
    zm = zp[near][np.nanargmax(C[near])]
    zr_ = zp[near][np.argmax(Cr[near])]
    zh = zp[~near][np.nanargmin(C[~near])] if (~near).any() else np.nan
    log("  z/zT, contrast clicks, contrast exact: " + " ".join(f"{z / zt:.2f}:{c:.3f}/{r:.3f}" for z, c, r in zip(zp[::5], C[::5], Cr[::5])))
    tol = max(0.01 * zt, TOL["talbot_rel"] * zt)
    report(rows := [], "Talbot self-image: z of best image (clicks)", round(zm, 1),
           f"{zt:.1f} = 2d²/λ  [exact RS best plane {zr_:.1f}]", f"±{tol:.1f}", abs(zm - zt) <= tol,
           f"image fraction in the slit windows at best plane {np.nanmax(C[near]):.3f} (open fraction a/d = {a / d:.2f}; exact {Cr[near].max():.3f})")
    report(rows, "Talbot half-distance: the image is shifted by d/2 (contrast minimum near zT/2)", round(zh, 1),
           f"{zt / 2:.1f}", f"±{tol:.1f}", abs(zh - zt / 2) <= tol, f"contrast there {np.nanmin(C[~near]):.3f}")
    meta = dict(photons=res["n"], clicks=res["counts"]["click"], counts=res["counts"], balls=balls_summary(sc),
                runtime=res["time"] + sc.build_time, z=zp.tolist(), C=C.tolist(), C_exact=Cr.tolist())
    save_results("talbot", rows, meta)
    return rows, meta


# ─────────────────────────────── 8. Fresnel edge ─────────────────────────────────────────────
def fresnel_edge(u):
    from scipy.special import fresnel
    S, Cf = fresnel(u)
    return 0.5 * ((Cf + 0.5) ** 2 + (S + 0.5) ** 2)


def exp_edge(n=60000):
    w, Le = 20.0, 100.0
    sc = Scene2D()
    src = sc.add(Obj2("edge", "source", (0, 0), (1, 0), 200.0, ints=[(0.0, 200.0)], src=dict(w0=w, zw=0.0)))
    sc.add(Obj2("screen", "screen", (0, Le), (1, 0), 1000.0))
    log("edge:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-400, 400, -1, Le + 1], seed=51)
    x = res["clicks"][res["status"] == CLICK][:, 0]
    np.save(os.path.join(OUT, "base_edge_clicks.npy"), x)
    rows = []
    nb = res["n"]                          # photons born on the open half (= the incident flux there)
    zr = np.pi * w ** 2
    wL = w * np.sqrt(1 + (Le / zr) ** 2)
    # density of clicks the FREE beam (no edge) would give, with the same incident flux per unit length
    rho_free = lambda xx: 2 * nb * np.exp(-2 * xx ** 2 / wL ** 2) * np.sqrt(2 / np.pi) / wL
    F = np.sqrt(LAM * Le / 2)
    # exact RS pattern of the same boundary values (for the ideal-sampler check of the estimators)
    s = np.arange(0.00625, 4 * w, 0.0125)
    xx = np.linspace(-90, 90, 6001)
    jz = np.concatenate([rs_screen(s, np.exp(-s ** 2 / w ** 2), q, Le)[1] for q in np.array_split(xx, 30)])
    b = 0.5
    est0 = lambda q: (np.abs(q) < b / 2).sum() / b / rho_free(0.0)
    i0, si = est0(x), bootstrap(est0, x)
    ic = ideal_check(lambda q: (np.abs(q) < b / 2).sum() / b / (2 * len(q) * np.sqrt(2 / np.pi) / wL), xx, jz, len(x))
    report(rows, "Fresnel edge: I(0)/I0 at the geometric shadow edge", f"{i0:.4f} ± {si:.4f}", 0.25,
           f"max(3σ, {TOL['edge_abs']})", abs(i0 - 0.25) <= max(3 * si, TOL["edge_abs"]),
           f"(window |x|<{b / 2}; ideal sampler of the exact pattern, normalised to its own clicks: {ic[0]:.4f} ± {ic[1]:.4f})")
    for name, u0, I0, mx in (("first maximum", 1.2172, 1.3704, True), ("first minimum", 1.8725, 0.7784, False),
                             ("second maximum", 2.3445, 1.1993, True)):
        xf = u0 * F
        shape = lambda xx_: np.exp(-2 * xx_ ** 2 / wL ** 2)
        est = lambda q: min_position(q, xf, 0.25 * F, maximum=mx, norm=shape)
        xm, sx = est(x), bootstrap(est, x)
        im = ideal_check(est, xx, jz, len(x))
        # theory for THIS beam (Gaussian w, not a plane wave): exact RS ratio to the free beam; the plane-wave
        # Fresnel numbers (u0, I0) are quoted for reference (they differ by up to 18% here: w = 20 is only 2.8 Fresnel scales)
        xs_ = np.linspace(xf - 0.3 * F, xf + 0.3 * F, 121)
        sfull = np.arange(-4 * w + 0.00625, 4 * w, 0.0125)
        j_e = rs_screen(s, np.exp(-s ** 2 / w ** 2), xs_, Le)[1]
        j_f = rs_screen(sfull, np.exp(-sfull ** 2 / w ** 2), xs_, Le)[1]
        ratio = j_e / j_f
        x_ex = xs_[np.argmax(ratio) if mx else np.argmin(ratio)]
        report(rows, f"Fresnel edge: {name} position", f"{xm:.3f} ± {sx:.3f}",
               f"{x_ex:.3f} (exact RS, this beam) [plane-wave Fresnel u={u0}: {xf:.3f}]",
               f"max(3σ, {100 * TOL['pos_rel']:.0f}%)", abs(xm - x_ex) <= max(3 * sx, TOL["pos_rel"] * x_ex),
               f"ideal sampler of the exact pattern {im[0]:.3f} ± {im[1]:.3f}")
        sel_ = np.abs(xs_ - xf) <= 0.1 * F
        hth = ratio[sel_].mean()
        hfr = np.mean(fresnel_edge(np.linspace(xf - 0.1 * F, xf + 0.1 * F, 201) / F))
        est_h = lambda q: (np.abs(q - xf) < 0.1 * F).sum() / (0.2 * F) / rho_free(xf)
        hm, sh = est_h(x), bootstrap(est_h, x)
        report(rows, f"Fresnel edge: {name} height I/I0 (window ±0.1·sqrt(λL/2))", f"{hm:.4f} ± {sh:.4f}",
               f"{hth:.4f} (exact RS, this beam) [plane-wave Fresnel {hfr:.4f}]",
               f"max(3σ, {TOL['edge_abs']})", abs(hm - hth) <= max(3 * sh, TOL["edge_abs"]))
    meta = dict(photons=res["n"], clicks=len(x), counts=res["counts"], balls=balls_summary(sc),
                runtime=res["time"] + sc.build_time)
    save_results("edge", rows, meta)
    return rows, meta


if __name__ == "__main__":
    what = sys.argv[1:] or ["free"]
    for w_ in what:
        t = time.time()
        rows, meta = {"free": exp_free, "double": exp_double, "which": exp_which, "single": exp_single,
                      "grating": exp_grating, "talbot": exp_talbot, "edge": exp_edge}[w_]()
        log(f"== {w_}: {'PASS' if all(r['ok'] for r in rows) else 'FAIL'}  [{time.time() - t:.0f}s]")

"""3-D experiments with the frozen base rule (lab/base_rule.py).

  python lab/base_exp3d.py free | gouy | slit3d | airy | arago | malus | three

Each prints measured / theory / tolerance / PASS-FAIL / photons / clicks / wave balls / runtime and saves
lab/out/base_<name>_<ver>.json.  Tolerances (TOL) are written before running and never changed.
The 3-D references are independent numpy quadratures of the vector Rayleigh–Sommerfeld integral
(rs3_*), not the ball machinery.
"""
import sys
import time
import numpy as np
from base_rule import *
import optics as OP
from base_exp2d import min_position, bootstrap, ideal_check, report

TOL = dict(
    width_rel=0.02,        # free beam width: max(3σ, 2%)
    gouy_abs=0.03,         # Gouy phase change from −zR to +zR: π/2 ± 0.03 rad
    corr=0.98, rms=0.08, absorbed=0.02,     # long double slit in 3-D: same thresholds as 2-D
    pos_rel=0.02,          # Airy first zero: max(3σ, 2%)
    arago_rel=0.10,        # Arago spot: fraction of clicks in the spot within max(3σ, 10%) of exact
    arago_contrast=3.0,    # spot density / shadow-ring density > 3 (theory ≈ 13)
    malus_abs=0.01,        # P(pass): max(3σ, 0.01)
)


# ─────────────────────────────── independent references (numpy, fine quadrature) ─────────────
def rs3_points(xs, ys, wB, P, z0=0.0, chunk=4000):
    """Vector RS field (x-polarised boundary values wB = Bx·dA at (xs, ys, z0)) at points P (m×3).
    Returns j_z (canonical current, Σ_a Im(E_a* ∂z E_a)) and E_x."""
    out_j, out_e = [], []
    for p in P:
        jz = 0
        dx = p[0] - xs
        dy = p[1] - ys
        z = p[2] - z0
        R = np.sqrt(dx * dx + dy * dy + z * z)
        e = np.exp(1j * K * R) / (4 * np.pi * R)
        h = (1j * K - 1 / R) * e
        h2 = (-K * K - 2j * K / R + 2 / R ** 2) * e
        hr = h / R
        q = (h2 - hr) / R ** 2
        Ex = (-2 * z * hr * wB).sum()
        dzEx = (-2 * (hr + z * z * q) * wB).sum()
        Ez = (2 * dx * hr * wB).sum()
        dzEz = (2 * dx * z * q * wB).sum()
        out_j.append(np.imag(np.conj(Ex) * dzEx + np.conj(Ez) * dzEz))
        out_e.append(Ex)
    return np.array(out_j), np.array(out_e)


def polar_grid(B, r0, r1, nr=500, nphi=256):
    rr = np.linspace(r0, r1, nr + 1)
    rr = 0.5 * (rr[1:] + rr[:-1])
    dr = (r1 - r0) / nr
    ph = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    R_, P_ = np.meshgrid(rr, ph)
    return (R_ * np.cos(P_)).ravel(), (R_ * np.sin(P_)).ravel(), (B(R_) * R_ * dr * 2 * np.pi / nphi).ravel()


def radial_ref(B, r0, r1, rho, L, nr=500, nphi=256):
    """azimuthal average (0°, 45°, 90°) of j_z on the screen at radius rho"""
    xs, ys, w = polar_grid(B, r0, r1, nr, nphi)
    acc = 0
    for ang in (0, np.pi / 4, np.pi / 2):
        P = np.column_stack([rho * np.cos(ang), rho * np.sin(ang), np.full(rho.size, L)])
        acc = acc + rs3_points(xs, ys, w, P)[0]
    return acc / 3


def meta_of(res, sc):
    return dict(photons=res["n"], clicks=res["counts"]["click"], counts=res["counts"], balls=balls_summary(sc),
                runtime=res["time"] + sc.build_time, wall_fraction=res["wall_frac"])


# ─────────────────────────────── 1. free beam (3-D) ─────────────────────────────────────────
def exp_free(n=7000):
    w0 = 3.0
    zr = np.pi * w0 ** 2
    sc = Scene3D()
    src = sc.add(Obj3("source", "source", 0.0, src=dict(w0=w0, zw=0.0)))
    zs = [zr, 2 * zr, 5 * zr]
    for z in zs:
        sc.add(Obj3(f"plane {z:.1f}", "plane", z))
    sc.add(Obj3("plane 1.0", "plane", 1.0))
    sc.add(Obj3("screen", "screen", 5 * zr + 1))
    log("free 3-D:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-150, 150, -150, 150, -1, 5 * zr + 2], seed=12, batch=1000)
    rows = []
    for k, z in enumerate(zs):
        P = np.array([p for i, pl, p in res["rec"] if pl == k + 1])
        wm = 2 * np.sqrt(0.5 * (P[:, 0].var() + P[:, 1].var()))
        wth = w0 * np.sqrt(1 + (z / zr) ** 2)
        sig = 1 / np.sqrt(4 * len(P))
        err = wm / wth - 1
        report(rows, f"3-D free beam width at z={z / zr:.0f} zR", round(wm, 4), round(wth, 4),
               f"max(3σ={3 * sig:.3f}, {TOL['width_rel']})", abs(err) <= max(3 * sig, TOL["width_rel"]),
               f"err {100 * err:+.2f}%  n={len(P)}")
    # the same photons: expansion w(z)/w(z=1) (independent of the statistical fluctuation of the birth sample)
    P1 = np.array([p for i, pl, p in res["rec"] if pl == 4])
    w1 = 2 * np.sqrt(0.5 * (P1[:, 0].var() + P1[:, 1].var()))
    for k, z in enumerate(zs):
        P = np.array([p for i, pl, p in res["rec"] if pl == k + 1])
        wm = 2 * np.sqrt(0.5 * (P[:, 0].var() + P[:, 1].var()))
        rth = np.sqrt(1 + (z / zr) ** 2) / np.sqrt(1 + (1 / zr) ** 2)
        report(rows, f"3-D free beam, same photons: expansion w({z / zr:.0f} zR)/w(1λ)", round(wm / w1, 4), round(rth, 4),
               f"±{TOL['width_rel']} rel.", abs(wm / w1 / rth - 1) <= TOL["width_rel"], f"w(1λ) = {w1:.4f}")
    meta = meta_of(res, sc)
    save_results("free3d", rows, meta)
    return rows, meta


# ─────────────────────────────── 1b. Gouy phase (3-D) ────────────────────────────────────────
def exp_gouy():
    w0 = 2.0
    zr = np.pi * w0 ** 2
    zs = -3 * zr                      # source plane 3 zR before the focus
    sc = Scene3D()
    sc.add(Obj3("source", "source", zs, src=dict(w0=w0, zw=-zs)))
    log("gouy:", balls_summary(sc))
    t0 = time.time()
    z = np.linspace(-2 * zr, 2 * zr, 161)
    P = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
    E, _ = sc.field(P, 0)
    ph = np.unwrap(np.angle(E[:, 0]) - K * (z - zs))
    gouy = -(ph - np.interp(0, z, ph))           # field ∝ e^{ikz − iζ(z)}
    # exact RS on the axis from the same boundary values (fine polar quadrature, independent of the balls)
    src = sc.objs[0]
    B = lambda r: src.src_B(r, np.zeros_like(r))[..., 0]
    xs, ys, w = polar_grid(B, 0, 4 * w0 * np.sqrt(10), 1500, 64)
    _, Er = rs3_points(xs, ys, w, P, z0=zs)
    phr = np.unwrap(np.angle(Er) - K * (z - zs))
    gr = -(phr - np.interp(0, z, phr))
    d_meas = np.interp(zr, z, gouy) - np.interp(-zr, z, gouy)
    d_rs = np.interp(zr, z, gr) - np.interp(-zr, z, gr)
    rows = []
    report(rows, "Gouy phase on axis, ζ(+zR) − ζ(−zR) (read by a photon ball on the axis)", f"{d_meas:.4f} rad",
           f"π/2 = {np.pi / 2:.4f} (paraxial); exact RS {d_rs:.4f}", f"±{TOL['gouy_abs']}",
           abs(d_meas - np.pi / 2) <= TOL["gouy_abs"],
           f"max |ζ − atan(z/zR)| over ±2zR: {np.abs(gouy - np.arctan(z / zr)).max():.4f} rad; vs exact RS {np.abs(gouy - gr).max():.2e}")
    d2 = np.interp(2 * zr, z, gouy) - np.interp(-2 * zr, z, gouy)
    report(rows, "Gouy phase on axis, ζ(+2zR) − ζ(−2zR)", f"{d2:.4f} rad", f"2·atan(2) = {2 * np.arctan(2):.4f}",
           f"±{TOL['gouy_abs']}", abs(d2 - 2 * np.arctan(2)) <= TOL["gouy_abs"])
    meta = dict(photons=0, clicks=0, balls=balls_summary(sc), runtime=time.time() - t0 + sc.build_time,
                note="wave-ball field read on the axis; no photons needed")
    save_results("gouy", rows, meta)
    return rows, meta


# ─────────────────────────────── 2'. long double slit in 3-D (2-D as a special case) ──────────
def exp_slit3d(n=6000):
    wy = 4.0
    rects = [(c - OP.A / 2, c + OP.A / 2, -60.0, 60.0) for c in OP.SLITS]
    sc = Scene3D()
    src = sc.add(Obj3("slits", "source", 0.0, mask=("rects", rects), src=dict(w0=OP.W0, w0y=wy, zw=0.0)))
    sc.add(Obj3("screen", "screen", OP.L))
    log("slit3d:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-200, 200, -200, 200, -1, OP.L + 1], seed=8, batch=700)
    cl = res["clicks"][res["status"] == CLICK]
    xl = cl[:, 0]
    np.save(os.path.join(OUT, "base_slit3d_clicks.npy"), cl)
    ref = OP.screen_reference("flux")
    c = OP.coarse_corr(xl, ref)
    e, rms, mx = OP.fringe_errors(xl, ref)
    base = OP.ideal_baseline(ref, len(xl), 200)
    # absorbed: fraction of |B|² on the wall (analytic, photons that would die there) + photons lost at the coin
    wall = res["wall_frac"]
    absorbed = 1 - (1 - wall) * len(xl) / res["n"]
    rows = []
    report(rows, "3-D long double slit (x-marginal) vs exact 2-D flux: correlation", round(c, 4),
           f"1 (ideal sampler {base[0]:.4f}±{base[1]:.4f}); 2-D rule: see Young", f">= {TOL['corr']}", c >= TOL["corr"])
    report(rows, "3-D long double slit: bright-fringe error rms / max", f"{100 * rms:.1f}% / {100 * mx:.1f}%",
           f"0 (ideal rms {100 * base[2]:.1f}%)", f"rms <= {100 * TOL['rms']:.0f}%", rms <= TOL["rms"],
           f"per fringe {np.round(100 * e).astype(int).tolist()}")
    report(rows, "3-D long double slit: absorbed fraction", round(absorbed, 4), "0.7806 (exact 2-D RS)",
           f"±{TOL['absorbed']}", abs(absorbed - 0.7806) <= TOL["absorbed"], f"wall {wall:.4f}")
    meta = meta_of(res, sc)
    save_results("slit3d", rows, meta)
    return rows, meta


# ─────────────────────────────── 5. Airy ─────────────────────────────────────────────────────
def exp_airy(n=10000):
    D, La = 6.0, 200.0
    sc = Scene3D()
    src = sc.add(Obj3("aperture", "source", 0.0, mask=("circ_in", 0, 0, D / 2), src=dict(w0=np.inf)))
    sc.add(Obj3("screen", "screen", La))
    log("airy:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-400, 400, -400, 400, -1, La + 1], seed=61, batch=1000)
    cl = res["clicks"][res["status"] == CLICK]
    np.save(os.path.join(OUT, "base_airy_clicks.npy"), cl)
    r = np.sqrt(cl[:, 0] ** 2 + cl[:, 1] ** 2)
    th = np.arcsin(1.2197 * LAM / D)
    rf = La * np.tan(th)
    rho = np.linspace(0, 2.2 * rf, 441)
    jr = radial_ref(lambda q: np.ones_like(q), 0, D / 2, rho, La)
    p_r = jr * rho                                   # density of the radius
    sel = np.abs(rho - rf) < 0.2 * rf
    r_rs = rho[sel][np.argmin(jr[sel])]
    win = 0.3 * rf
    est = lambda q: min_position(q, rf, win, norm=lambda x: x)
    rm, sr = est(r), bootstrap(est, r)
    ic = ideal_check(est, rho, p_r, len(r))
    rows = []
    report(rows, f"Airy disk D={D}: first zero radius", f"{rm:.2f} ± {sr:.2f}",
           f"{rf:.2f} = L·tan(asin(1.22λ/D)) [paraxial 1.22λL/D = {1.22 * LAM * La / D:.2f}; exact RS {r_rs:.2f}]",
           f"max(3σ, {100 * TOL['pos_rel']:.0f}%)", abs(rm - rf) <= max(3 * sr, TOL["pos_rel"] * rf),
           f"err {100 * (rm / rf - 1):+.2f}%; same estimator on ideal samples of the exact pattern {ic[0]:.2f} ± {ic[1]:.2f}; "
           f"fraction inside the first zero, of the clicks within 2.2 r1: {np.mean(r[r < 2.2 * rf] < rf):.3f} (exact {np.trapezoid(p_r[rho < r_rs], rho[rho < r_rs]) / np.trapezoid(p_r, rho):.3f} within 2.2 r1)")
    meta = meta_of(res, sc)
    save_results("airy", rows, meta)
    return rows, meta


# ─────────────────────────────── 6. Arago / Poisson spot ─────────────────────────────────────
def exp_arago(n=7000):
    a, w, La = 3.0, 5.0, 12.0
    sc = Scene3D()
    src = sc.add(Obj3("disk", "source", 0.0, mask=("circ_out", 0, 0, a), src=dict(w0=w, zw=0.0)))
    sc.add(Obj3("screen", "screen", La))
    log("arago:", balls_summary(sc))
    res = simulate(sc, src, n, box=[-60, 60, -60, 60, -1, La + 1], seed=71, batch=1000)
    cl = res["clicks"][res["status"] == CLICK]
    np.save(os.path.join(OUT, "base_arago_clicks.npy"), cl)
    r = np.sqrt(cl[:, 0] ** 2 + cl[:, 1] ** 2)
    rs = 2.405 * La / (K * a)                  # first zero of the J0 spot
    B = lambda q: np.exp(-q ** 2 / w ** 2)
    rho = np.linspace(0, 4 * w, 801)
    jr = radial_ref(B, a, 4 * w, rho, La)
    tot = np.trapezoid(jr * rho, rho)
    spot = rho < rs / 2
    p_th = np.trapezoid((jr * rho)[spot], rho[spot]) / tot
    ring = (rho > rs) & (rho < 0.8 * a)
    contrast_th = jr[spot].mean() / jr[ring].mean()
    est = lambda q: np.mean(q < rs / 2)
    p_m, sp = est(r), bootstrap(est, r)
    n_ring = ((r > rs) & (r < 0.8 * a)).sum() / (np.pi * ((0.8 * a) ** 2 - rs ** 2))
    n_spot = (r < rs / 2).sum() / (np.pi * (rs / 2) ** 2)
    contrast = n_spot / max(n_ring, 1e-300)
    jf = radial_ref(B, 0, 4 * w, np.array([0.0]), La)[0]
    rows = []
    report(rows, f"Arago spot (disk a={a}, beam w={w}, L={La}): fraction of clicks within r < {rs / 2:.2f}",
           f"{p_m:.4f} ± {sp:.4f}", f"{p_th:.4f} (exact RS)", f"max(3σ, {100 * TOL['arago_rel']:.0f}%)",
           abs(p_m - p_th) <= max(3 * sp, TOL["arago_rel"] * p_th),
           f"{(r < rs / 2).sum()} clicks in the spot; exact I(0)/I_free(0) = {jr[0] / jf:.3f}")
    report(rows, "Arago spot: density in the spot / density in the shadow ring", f"{contrast:.1f}", f"{contrast_th:.1f} (exact RS)",
           f"> {TOL['arago_contrast']}", contrast > TOL["arago_contrast"])
    meta = meta_of(res, sc)
    save_results("arago", rows, meta)
    return rows, meta


# ─────────────────────────────── 11. Malus and three polarizers (3-D) ────────────────────────
def pol(theta):
    a = np.array([np.cos(theta), np.sin(theta)], complex)
    return np.outer(a, a.conj())


def malus_scene(angles, jones):
    w0 = 2.0
    sc = Scene3D()
    src = sc.add(Obj3("source", "source", 0.0, src=dict(w0=w0, zw=0.0, jones=jones)))
    for k, th in enumerate(angles):
        sc.add(Obj3(f"polarizer {np.degrees(th):.1f}°", "polarizer", 5.0 * (k + 1), T=pol(th)))
    zscr = 5.0 * (len(angles) + 1)
    sc.add(Obj3("screen", "screen", zscr))
    return sc, src, zscr


def run_malus(angles, jones, n, seed):
    sc, src, zscr = malus_scene(angles, jones)
    res = simulate(sc, src, n, box=[-60, 60, -60, 60, -1, zscr + 1], seed=seed, batch=1000, verbose=False)
    passed = res["counts"]["click"]
    born = res["n"] - res["src_coin"]
    return passed, born, res, sc


def exp_malus(n=1500):
    rows = []
    t0 = time.time()
    metas = []
    for deg in (0, 22.5, 45, 67.5, 90):
        th = np.radians(deg)
        k, nb, res, sc = run_malus([th], (1, 0), n, seed=100 + int(deg * 10))
        P = k / nb
        Pth = np.cos(th) ** 2
        sig = np.sqrt(max(Pth * (1 - Pth), 1 / nb) / nb)
        report(rows, f"Malus (absorbing polarizer, 3-D) θ={deg}°: P(pass)", f"{P:.4f}", f"{Pth:.4f} = cos²θ",
               f"max(3σ={3 * sig:.4f}, {TOL['malus_abs']})", abs(P - Pth) <= max(3 * sig, TOL["malus_abs"]),
               f"{k}/{nb} photons")
        metas.append(dict(angle=deg, **meta_of(res, sc)))
    save_results("malus", rows, dict(runs=metas, runtime=time.time() - t0))
    return rows, metas


def exp_three(n=4000):
    rows = []
    t0 = time.time()
    metas = []
    angles = [0.0, np.pi / 4, np.pi / 2]
    # polarized along the first polarizer
    k, nb, res, sc = run_malus(angles, (1, 0), n, seed=301)
    P = k / nb
    sig = np.sqrt(0.25 * 0.75 / nb)
    report(rows, "three polarizers 0/45/90°, light polarized at 0°: P(pass)", f"{P:.4f}", "0.2500",
           f"max(3σ={3 * sig:.4f}, {TOL['malus_abs']})", abs(P - 0.25) <= max(3 * sig, TOL["malus_abs"]), f"{k}/{nb}")
    metas.append(meta_of(res, sc))
    k2, nb2, res2, sc2 = run_malus(angles[::2], (1, 0), n // 2, seed=303)
    report(rows, "two crossed polarizers 0/90° (the middle one removed): P(pass)", f"{k2 / nb2:.4f}", "0",
           "0 exactly", k2 == 0, f"{k2}/{nb2}")
    # unpolarized: a 50/50 mixture of H and V photons (each photon has its own pure polarization)
    rng = np.random.default_rng(305)
    nH = rng.binomial(n, 0.5)
    kH, nbH, rH, scH = run_malus(angles, (1, 0), nH, seed=306)
    kV, nbV, rV, scV = run_malus(angles, (0, 1), n - nH, seed=307)
    P = (kH + kV) / (nbH + nbV)
    sig = np.sqrt(0.125 * 0.875 / (nbH + nbV))
    report(rows, "three polarizers 0/45/90°, unpolarized light: P(pass)", f"{P:.4f}", "0.1250",
           f"max(3σ={3 * sig:.4f}, {TOL['malus_abs']})", abs(P - 0.125) <= max(3 * sig, TOL["malus_abs"]),
           f"{kH + kV}/{nbH + nbV} (H: {kH}/{nbH}, V: {kV}/{nbV})")
    metas += [meta_of(rH, scH), meta_of(rV, scV)]
    save_results("three", rows, dict(runs=metas, runtime=time.time() - t0))
    return rows, metas


if __name__ == "__main__":
    what = sys.argv[1:] or ["free"]
    for w_ in what:
        t = time.time()
        rows, meta = {"free": exp_free, "gouy": exp_gouy, "slit3d": exp_slit3d, "airy": exp_airy,
                      "arago": exp_arago, "malus": exp_malus, "three": exp_three}[w_]()
        log(f"== {w_}: {'PASS' if all(r['ok'] for r in rows) else 'FAIL'}  [{time.time() - t:.0f}s]")

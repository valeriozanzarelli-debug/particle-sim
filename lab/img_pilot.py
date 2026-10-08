"""Imaging with the frozen base rule (level 1/2: wave balls + photon balls), lab/base_rule.py v2.

  python lab/img_pilot.py [n_photons]

Scene (units λ = 1): a point emitter at distance s_o in front of an ideal thin lens (quadratic phase
exp(-i k r²/2f), allowed: it is a phase screen, no material), pupil radius a, sensor at the image distance
s_i (1/s_o + 1/s_i = 1/f, magnification 1).  The source object is the lens plane: B = (spherical wave of the
emitter) × (lens phase) inside the pupil.  The photon balls are born on the pupil with probability |B|²
and guided by the wave balls to the sensor, exactly as in the validated Airy experiment.

What is measured (the PSF = where single photons land):
  * radius of the first dark ring from the clicks vs the exact vector Rayleigh–Sommerfeld reference of the
    same pupil field (independent quadrature, base_exp3d.rs3_points) and vs the textbook 0.61 λ/NA;
  * the click histogram vs the exact wave (level 3): correlation and rms error per radial bin, compared with
    an ideal sampler with the same number of clicks (so the error of level 1/2 vs level 3 is separated
    from shot noise);
  * photons simulated per second on this PC.
Second emitter: placed at the Rayleigh distance 0.61 λ/NA_o; its photons are run in their OWN scene (no
isoplanatism assumption) and the two click sets are merged (incoherent emitters: independent photons).
"""
import json
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base_rule import Obj3, Scene3D, simulate, balls_summary, CLICK, LAM, K, OUT   # noqa: E402
from base_exp3d import rs3_points, polar_grid                                      # noqa: E402
from base_exp2d import min_position, bootstrap, ideal_check, ideal_samples        # noqa: E402

A_PUP, S_O, S_I = 8.0, 40.0, 40.0          # pupil radius, object and image distances (λ)
F = 1 / (1 / S_O + 1 / S_I)
NA = A_PUP / np.hypot(A_PUP, S_I)          # image-side numerical aperture (= object side, M = 1)
TOL = dict(r1_rel=0.03, corr=0.98, sep_rel=0.05)   # written before the runs


class LensSource(Obj3):
    """pupil plane = source: emitter spherical wave × thin-lens quadratic phase (x-polarised)."""

    def __init__(self, xe, ye):
        super().__init__("lens", "source", 0.0, mask=("circ_in", 0, 0, A_PUP), src=dict(w0=np.inf))
        self.xe, self.ye = xe, ye

    def src_B(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        R = np.sqrt((x - self.xe) ** 2 + (y - self.ye) ** 2 + S_O ** 2)
        u = np.exp(1j * K * R) / R * S_O * np.exp(-1j * K * (x * x + y * y) / (2 * F))
        return u[..., None] * np.array([1, 0], complex)


def run_emitter(xe, n, seed):
    sc = Scene3D()
    src = sc.add(LensSource(xe, 0.0))
    sc.add(Obj3("sensor", "screen", S_I))
    print("  balls:", balls_summary(sc), flush=True)
    t0 = time.time()
    res = simulate(sc, src, n, box=[-40, 40, -40, 40, -1, S_I + 1], seed=seed, batch=500, verbose=True)
    dt = time.time() - t0
    cl = res["clicks"][res["status"] == CLICK]
    return cl, dt, sc.build_time


def exact_intensity(xe, X, Y):
    """level-3 reference: |E|-flux j_z on the sensor from the same pupil field (independent quadrature)."""
    src = LensSource(xe, 0.0)
    xs, ys, w = polar_grid(lambda r: np.ones_like(r), 0, A_PUP, nr=160, nphi=256)
    wB = src.src_B(xs, ys)[:, 0] * w
    P = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, S_I)])
    j, _ = rs3_points(xs, ys, wB, P)
    return j.reshape(X.shape)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 4000
    out = {}
    r_ray = 0.61 * LAM / NA
    print(f"lens: a={A_PUP} s_o={S_O} s_i={S_I} f={F:.2f} NA={NA:.4f}  0.61λ/NA={r_ray:.3f}λ", flush=True)
    # ---- exact wave PSF (radial, azimuthal average over 0/45/90 deg) ----
    t0 = time.time()
    rho = np.linspace(0, 2.2 * r_ray, 221)
    acc = 0
    for ang in (0, np.pi / 4, np.pi / 2):
        acc = acc + exact_intensity(0.0, rho * np.cos(ang), rho * np.sin(ang))
    jr = acc / 3
    sel = np.abs(rho - r_ray) < 0.25 * r_ray
    r1_exact = rho[sel][np.argmin(jr[sel])]
    print(f"exact RS first zero r1 = {r1_exact:.3f} (0.61λ/NA = {r_ray:.3f}), {time.time() - t0:.1f}s", flush=True)
    # ---- photons from the pilot rule ----
    ana = "--analyse" in sys.argv
    if ana:   # re-analyse the clicks saved by the last full run (no new photons)
        clA = np.load(os.path.join(OUT, "img_pilot_clicks_A.npy"))
        js = json.load(open(os.path.join(OUT, "img_pilot.json")))
        n, dtA, bA = js["n"], js["t_run"], js["t_build"]
    else:
        clA, dtA, bA = run_emitter(0.0, n, seed=101)
    ra = np.hypot(clA[:, 0], clA[:, 1])
    rate = n / dtA
    print(f"emitter A: {n} photons, {len(clA)} clicks, {dtA:.1f}s run + {bA:.1f}s build -> {rate:.1f} photons/s", flush=True)
    # radial histogram vs exact
    edges = np.linspace(0, 2.2 * r_ray, 23)
    h, _ = np.histogram(ra, edges)
    pr = jr * rho
    cdf = np.concatenate([[0], np.cumsum(0.5 * (pr[1:] + pr[:-1]) * np.diff(rho))])
    pe = np.diff(np.interp(edges, rho, cdf))
    m = (ra < edges[-1]).sum()
    expct = pe / pe.sum() * m
    corr = np.corrcoef(h, expct)[0, 1]
    bright = expct > 0.02 * expct.max()
    err = np.sqrt(np.mean(((h - expct) / expct)[bright] ** 2))
    rng = np.random.default_rng(7)
    ideal = []
    for _ in range(200):
        hi = rng.multinomial(m, pe / pe.sum())
        ideal.append(np.sqrt(np.mean(((hi - expct) / expct)[bright] ** 2)))
    ideal = np.array(ideal)
    chi2 = (((h - expct) ** 2 / np.maximum(expct, 1))[expct > 5]).sum() / max((expct > 5).sum() - 1, 1)
    # first-zero radius from the clicks.  Estimator = the project's validated one (base_exp2d.min_position:
    # weighted cubic on 30 bins in r1 ± 0.3 r1 of the radial density / r), checked on ideal samples of the
    # exact pattern.  (First run used a symmetric parabola: biased on this skew minimum, see parabola_check.)
    est = lambda q: min_position(q, r_ray, 0.3 * r_ray, norm=lambda x: x)
    r1_clicks, sr = est(ra), bootstrap(est, ra)
    ic = ideal_check(est, rho, jr * rho, len(ra))

    def parabola(q):
        rr = q[np.abs(q - r_ray) < 0.3 * r_ray]
        hh, ee = np.histogram(rr, 15)
        cc = 0.5 * (ee[1:] + ee[:-1])
        co = np.polyfit(cc, hh / cc, 2)
        return -co[1] / (2 * co[0])
    par_clicks = parabola(ra)
    par_ideal = np.mean([parabola(ideal_samples(rho, jr * rho, len(ra), 500 + k)) for k in range(40)])
    ok1 = abs(r1_clicks - r1_exact) <= max(3 * sr, TOL["r1_rel"] * r1_exact)
    print(f"TEST pilot PSF first zero: clicks {r1_clicks:.3f} ± {sr:.3f}  exact RS {r1_exact:.3f}  0.61λ/NA {r_ray:.3f}"
          f"  (same estimator on ideal samples of the exact wave {ic[0]:.3f} ± {ic[1]:.3f})  -> {'PASS' if ok1 else 'FAIL'}", flush=True)
    print(f"     parabola_check (estimator of the first run): clicks {par_clicks:.3f}, ideal samples {par_ideal:.3f}"
          f" -> the bias is the estimator's, not the photons'", flush=True)
    ok2 = corr >= TOL["corr"]
    print(f"TEST pilot PSF radial histogram vs exact wave: corr {corr:.4f}, rms err (bright bins) {100 * err:.1f}%"
          f" (ideal sampler {100 * ideal.mean():.1f} ± {100 * ideal.std():.1f}%), chi2/dof {chi2:.2f} -> {'PASS' if ok2 else 'FAIL'}",
          flush=True)
    frac_in = np.mean(ra < r1_exact)
    fe = np.interp(r1_exact, rho, cdf) / cdf[-1]
    print(f"     fraction of clicks inside r1 (of clicks within 2.2 r1): {np.mean(ra[ra < edges[-1]] < r1_exact):.3f}"
          f" exact {fe:.3f}", flush=True)
    out.update(NA=NA, n=n, clicks=len(clA), r1_clicks=r1_clicks, r1_err=sr, r1_exact=r1_exact, r_rayleigh=r_ray,
               corr=corr, rms=err, rms_ideal=ideal.mean(), chi2=chi2, photons_per_s=rate, t_run=dtA, t_build=bA)
    # ---- second emitter at the Rayleigh distance (object side), its own scene ----
    d = 0.61 * LAM / NA
    if ana:
        clB, dtB = np.load(os.path.join(OUT, "img_pilot_clicks_B.npy")), js["t_run_B"]
    else:
        clB, dtB, bB = run_emitter(d, n, seed=202)
    # image of an emitter at +d is at -d (M = -1)
    cx_a, cx_b = np.median(clA[:, 0]), np.median(clB[:, 0])
    sepm = abs(cx_b - cx_a)
    okc = abs(sepm - d) <= TOL["sep_rel"] * d
    print(f"TEST emitter B image centre: separation of the two click clouds {sepm:.3f} (expected {d:.3f}) -> {'PASS' if okc else 'FAIL'}", flush=True)
    both = np.vstack([clA, clB])
    np.save(os.path.join(OUT, "img_pilot_clicks_A.npy"), clA)
    np.save(os.path.join(OUT, "img_pilot_clicks_B.npy"), clB)
    # cut along x through the two images (|y| < 0.5): dip between peaks vs exact incoherent sum
    selb = np.abs(both[:, 1]) < 0.5
    pk_c = ((np.abs(both[:, 0]) < 0.4) & selb).sum() + ((np.abs(both[:, 0] + d) < 0.4) & selb).sum()
    dip_c = 2 * ((np.abs(both[:, 0] + d / 2) < 0.4) & selb).sum() / pk_c
    # exact incoherent sum, same windows (0.8 x 1.0 λ boxes), and pointwise
    gx, gy = np.meshgrid(np.linspace(-0.4, 0.4, 9), np.linspace(-0.5, 0.5, 11))
    def box_mean(x0):
        return (exact_intensity(0.0, gx + x0, gy) + exact_intensity(d, gx + x0, gy)).mean()
    dip_e = 2 * box_mean(-d / 2) / (box_mean(0.0) + box_mean(-d))
    pt = lambda x0: exact_intensity(0.0, np.array([x0]), np.array([0.0]))[0] + exact_intensity(d, np.array([x0]), np.array([0.0]))[0]
    dip_pt = pt(-d / 2) / pt(0.0)
    sd = dip_c * np.sqrt(1 / max(pk_c, 1) + 1 / max(dip_c * pk_c / 2, 1))
    print(f"     two emitters at 0.61λ/NA: dip/peak from clicks {dip_c:.3f} ± {sd:.3f}, exact same windows {dip_e:.3f},"
          f" exact pointwise {dip_pt:.3f} (paraxial Airy theory 0.735)", flush=True)
    out.update(dip_err=sd, dip_point=dip_pt, r1_ideal=ic[0], r1_ideal_sd=ic[1], par_clicks=par_clicks, par_ideal=par_ideal)
    out.update(sep_measured=sepm, sep=d, dip_clicks=dip_c, dip_exact=dip_e, t_run_B=dtB)
    passed = sum([ok1, ok2, okc])
    out["tests"] = f"{passed}/3"
    with open(os.path.join(OUT, "img_pilot.json"), "w") as f:
        json.dump({k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in out.items()}, f, indent=1)
    print(f"RESULT: {passed}/3 PASS", flush=True)


if __name__ == "__main__":
    main()

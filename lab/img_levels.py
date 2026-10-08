"""Task 4 — zoom levels for images: when a group of photons can travel as ONE "giant photon".

  python lab/img_levels.py

Chain used for images (units λ = 1, NA = 0.5 unless stated):
  far   : level 4, rays.  A "giant photon" = one ray carrying w photons, from a point of the source through
          the ideal thin lens (paraxial ray transfer, the lens is a phase screen) to its geometric image point;
  sensor: level 3, the wave.  The ray's landing point is spread by the PSF of the pupil (Airy) — either as
          the expected intensity (deterministic) or by drawing a landing offset per photon;
  clicks: level 2, photon statistics.  The Camera draws Poisson clicks from the expected image.

Tests (tolerances written before the run):
  L1  rays alone (no PSF) put all photons of a point emitter in the pixel of its geometric image.  Error
      = 1 − ensquared energy.  Measured vs pixel pitch; must fall monotonically; report pitch for 10 %, 1 %.
      Check: ensquared energy of a pixel of pitch p → encircled energy of the inscribed/circumscribed
      circles bracket it (EE(p/2) ≤ ESE ≤ EE(p/√2)), exact formula 1 − J0² − J1².
  L2  extended incoherent source (uniform disc of radius R): total-variation distance between the ray
      image and the wave image (disc ⊗ Airy) vs R/(λ/NA); fitted law TV = c·(λ/NA)/R with exponent
      −1 ± 0.1 for R ≥ 4 λ/NA; report R for 10 % and 1 %.
  L3  giant photons of weight w, r rays per pixel (Monte Carlo of the expected image): relative rms
      error in the bright pixels ∝ r^−1/2, fitted exponent −0.5 ± 0.03.
  L4  "one ray = one photon" (w = 1) gives the PHYSICAL shot noise: Fano factor of the pixel counts
      over repeated exposures = 1 within 3σ; giant photons with w = 10 give Fano = 10 (wrong noise)
      unless the Camera redraws Poisson clicks from their mean image: Fano back to 1 within 3σ.
  L5  rays through the ideal thin lens land on the geometric image (M = −s_i/s_o) within 1e-9 λ.
Speed on this PC: rays/s, PSF-sampled photons/s, wave-level image (FFT convolution) time, camera pixels/s,
and the pilot rule's photons/s from lab/out/img_pilot.json.
"""
import json
import os
import sys
import time
import numpy as np
from scipy.special import j0, j1
from scipy.signal import fftconvolve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from img_core import Camera, pixel_psf, airy, sample_airy_offsets   # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
LAM, NA = 1.0, 0.5
U = LAM / NA          # natural unit of the image: λ/NA


def ee(r):
    v = 2 * np.pi * NA * r / LAM
    return 1 - j0(v) ** 2 - j1(v) ** 2


def trace_rays(xe, ye, n, so, si, a, rng):
    """level 4: n rays from a point emitter at (xe, ye, -so) through a thin lens (pupil radius a, focal f
    with 1/so + 1/si = 1/f) to the plane z = si.  Paraxial ray transfer (the lens is a quadratic phase)."""
    f = 1 / (1 / so + 1 / si)
    r = a * np.sqrt(rng.random(n))
    ph = rng.random(n) * 2 * np.pi
    px, py = r * np.cos(ph), r * np.sin(ph)            # hit point on the pupil (uniform = isotropic emitter, small NA)
    ux, uy = (px - xe) / so, (py - ye) / so            # slopes before the lens
    ux, uy = ux - px / f, uy - py / f                  # thin lens
    return px + ux * si, py + uy * si


def main():
    t0 = time.time()
    rng = np.random.default_rng(44)
    passed = total = 0

    # ---------- L1: point emitter, rays vs wave, vs pixel pitch ----------
    pitches = np.array([0.25, 0.5, 1, 2, 4, 8, 16, 32]) * U
    err = []
    ok_br = True
    for p in pitches:
        sub = 400
        s = ((np.arange(sub) + 0.5) / sub - 0.5) * p
        X, Y = np.meshgrid(s, s)
        ese = airy(np.hypot(X, Y), NA, LAM).mean() * p * p / (LAM ** 2 / (np.pi * NA ** 2))
        err.append(1 - ese)
        ok_br &= ee(p / 2) - 1e-3 <= ese <= ee(p / np.sqrt(2)) + 1e-3
    err = np.array(err)
    mono = np.all(np.diff(err) < 0)
    ok = mono and ok_br
    passed += ok; total += 1
    lp = np.log(pitches / U)
    le = np.log(err)
    p10 = np.exp(np.interp(np.log(0.10), le[::-1], lp[::-1]))
    p01 = np.exp(np.interp(np.log(0.01), le[::-1], lp[::-1]))
    print("L1 point emitter, rays only: error (photons outside the geometric pixel) vs pitch: " +
          ", ".join(f"{q / U:g} λ/NA: {100 * e:.2f}%" for q, e in zip(pitches, err)), flush=True)
    print(f"   10 % at pitch ≈ {p10:.2f} λ/NA, 1 % at pitch ≈ {p01:.1f} λ/NA (error ∝ 1/pitch: the Airy tail)"
          f"; bracket check {'ok' if ok_br else 'BROKEN'} -> {'PASS' if ok else 'FAIL'}", flush=True)

    # ---------- L2: extended disc, rays vs wave ----------
    px = 0.1 * U
    Rs = np.array([1, 2, 4, 8, 16, 32]) * U
    tvs = []
    for R in Rs:
        n = int(2 * (R + 12 * U) / px) | 1
        c = (np.arange(n) - n // 2) * px
        X, Y = np.meshgrid(c, c)
        disc = (np.hypot(X, Y) <= R).astype(float)
        disc /= disc.sum()
        m = int(24 * U / px) | 1
        k = (np.arange(m) - m // 2) * px
        KX, KY = np.meshgrid(k, k)
        psf = airy(np.hypot(KX, KY), NA, LAM)
        psf /= psf.sum()
        wave = fftconvolve(disc, psf, mode="same")
        tvs.append(0.5 * np.abs(wave - disc).sum() + 0.5 * (1 - wave.sum()))
    tvs = np.array(tvs)
    sel = Rs >= 4 * U
    slope, icpt = np.polyfit(np.log(Rs[sel] / U), np.log(tvs[sel]), 1)
    ok = abs(slope + 1) <= 0.1
    passed += ok; total += 1
    r10 = np.exp((np.log(0.10) - icpt) / slope)
    r01 = np.exp((np.log(0.01) - icpt) / slope)
    print("L2 uniform disc, rays vs wave (TV distance): " + ", ".join(f"R={R / U:g} λ/NA: {100 * t:.2f}%" for R, t in zip(Rs, tvs)),
          flush=True)
    print(f"   fit TV = {np.exp(icpt):.3f}·(R/(λ/NA))^{slope:.3f} (expect −1) -> {'PASS' if ok else 'FAIL'};"
          f" rays good to 10 % for R > {r10:.1f} λ/NA, to 1 % for R > {r01:.0f} λ/NA (on a {px / U:g} λ/NA grid)", flush=True)

    # ---------- L3: giant photons, Monte Carlo noise vs rays per pixel ----------
    R = 3 * U
    pitch = 0.25 * U
    npx = 48
    edges = (np.arange(npx + 1) - npx / 2) * pitch
    # expected image: very many PSF-sampled photons (reference)
    def mc_image(nrays, w, rng):
        rr = R * np.sqrt(rng.random(nrays))
        ph = rng.random(nrays) * 2 * np.pi
        dx, dy = sample_airy_offsets(nrays, NA, LAM, rng)
        h, _, _ = np.histogram2d(rr * np.sin(ph) + dy, rr * np.cos(ph) + dx, [edges, edges])
        return h * w
    ref = sum(mc_image(2_000_000, 1.0, rng) for _ in range(10)) / 20_000_000
    bright = ref > 0.5 * ref.max()
    npb = bright.sum()
    rpp, rel = [], []
    for nr in (1e3, 1e4, 1e5, 1e6):
        e = []
        for k in range(6):
            img = mc_image(int(nr), 1.0 / nr, rng)
            e.append(np.sqrt(np.mean(((img - ref) / ref)[bright] ** 2)))
        rpp.append(nr * ref[bright].mean())
        rel.append(np.mean(e))
    sl = np.polyfit(np.log(rpp), np.log(rel), 1)[0]
    ok = abs(sl + 0.5) <= 0.03
    passed += ok; total += 1
    print("L3 giant photons: MC rms error in bright pixels vs rays per bright pixel: " +
          ", ".join(f"{r:.0f}: {100 * e:.2f}% (1/√r = {100 / np.sqrt(r):.2f}%)" for r, e in zip(rpp, rel)), flush=True)
    print(f"   exponent {sl:.3f} (−0.5) -> {'PASS' if ok else 'FAIL'}: for an image noise ε one needs ≈ 1/ε² rays per pixel"
          f" (1 %: 10⁴, 10 %: 100), whatever the number of real photons", flush=True)

    # ---------- L4: noise of one ray = one photon vs giant photons ----------
    mean_ph = 50.0                                        # photons per bright pixel per exposure
    ntot = mean_ph / ref[bright].mean()
    reps = 300
    s1, s10, s10c = [], [], []
    cam = Camera()
    for k in range(reps):
        n1 = rng.poisson(ntot)
        s1.append(mc_image(n1, 1.0, rng)[bright])
        n10 = rng.poisson(ntot / 10)
        g = mc_image(n10, 10.0, rng)[bright]
        s10.append(g)
        s10c.append(cam.expose(ref[bright] * ntot, 1.0, rng))   # mean from giant photons (exact here) → Poisson clicks
    def fano(a):
        a = np.array(a)
        f = a.var(0, ddof=1) / a.mean(0)
        return f.mean(), f.std() / np.sqrt(f.size) if f.size > 1 else np.nan
    f1, f10, f10c = fano(s1), fano(s10), fano(s10c)
    sig = np.sqrt(2 / (reps - 1)) / np.sqrt(npb)             # σ of a mean of npb Fano estimates (Poisson, large N)
    ok = abs(f1[0] - 1) <= 3 * max(sig, f1[1]) and abs(f10c[0] - 1) <= 3 * max(sig, f10c[1]) and abs(f10[0] / 10 - 1) < 0.1
    passed += ok; total += 1
    print(f"L4 Fano of pixel counts ({reps} exposures, {npb} bright pixels, {mean_ph:g} photons/pixel):"
          f" 1 ray = 1 photon {f1[0]:.4f}; giant photons w=10 {f10[0]:.3f}; giant photons → mean → Camera Poisson {f10c[0]:.4f}"
          f" (σ ≈ {sig:.4f}) -> {'PASS' if ok else 'FAIL'}", flush=True)

    # ---------- L5: thin-lens ray transfer ----------
    xi, yi = trace_rays(1.3, -0.7, 100000, 40.0, 40.0, 8.0, rng)
    dev = max(np.abs(xi + 1.3).max(), np.abs(yi - 0.7).max())
    xi2, yi2 = trace_rays(1.3, -0.7, 1000, 30.0, 60.0, 8.0, rng)
    dev2 = max(np.abs(xi2 + 2.6).max(), np.abs(yi2 - 1.4).max())
    ok = max(dev, dev2) < 1e-9
    passed += ok; total += 1
    print(f"L5 rays through the ideal thin lens: max distance from the geometric image {max(dev, dev2):.1e} λ (M = −1 and −2)"
          f" -> {'PASS' if ok else 'FAIL'}", flush=True)

    # ---------- speeds ----------
    n = 2_000_000
    t = time.time(); trace_rays(1.0, 0.0, n, 40.0, 40.0, 8.0, rng); tr = n / (time.time() - t)
    t = time.time(); sample_airy_offsets(n, NA, LAM, rng); ts = n / (time.time() - t)
    t = time.time(); mc_image(n, 1.0, rng); tm = n / (time.time() - t)
    big = np.zeros((1024, 1024)); big[400:600, 400:600] = 1
    k = pixel_psf(129, 0.25 * U, NA, LAM)
    t = time.time(); fftconvolve(big, k, mode="same"); tf = time.time() - t
    t = time.time(); Camera(qe=0.8, dark_rate=0.01, read_noise=1.0).expose(np.full((2048, 2048), 30.0), 1.0, rng); tc = 2048 ** 2 / (time.time() - t)
    pil = None
    pj = os.path.join(OUT, "img_pilot.json")
    if os.path.exists(pj):
        pil = json.load(open(pj))["photons_per_s"]
    print(f"SPEED (this PC, numpy, 1 process): rays through the lens {tr:.2e}/s; Airy landing draws {ts:.2e}/s;"
          f" full photon (disc point + PSF + pixel histogram) {tm:.2e}/s; wave-level image 1024² ⊗ 129² PSF (FFT) {tf:.2f} s;"
          f" camera {tc:.2e} pixels/s; pilot rule (level 1/2) {pil if pil is None else f'{pil:.1f}'} photons/s", flush=True)
    if pil:
        print(f"   ratio level 4+3 photons / pilot photons ≈ {tm / pil:.1e}", flush=True)
    print(f"time {time.time() - t0:.1f}s", flush=True)
    print(f"RESULT: {passed}/{total} PASS", flush=True)


if __name__ == "__main__":
    main()

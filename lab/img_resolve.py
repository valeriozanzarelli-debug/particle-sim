"""Task 2 (wave level) and Task 5 — resolving two point emitters; synthetic images.

  python lab/img_resolve.py

Level used: 3 (wave) for the expected image = incoherent sum of two Airy patterns (scalar, paraxial),
then photon clicks drawn by the Camera (Poisson).  Two independent atoms are incoherent emitters: their
photons never interfere with each other, so the image is the SUM of intensities (not of fields).
Error of level 3 vs the exact vector Rayleigh–Sommerfeld wave at NA 0.196: see img_pilot.py (first zero
+1.0 %).  Error of the photon balls (level 1/2) vs level 3: img_pilot.py.

Tests (tolerances written before the run):
  R1  first zero of the Airy pattern at 0.6098 λ/NA (zero of J1 / 2π): |Δ| < 1e-4.
  R2  two emitters at the Rayleigh distance 0.61 λ/NA: dip/peak of the noiseless profile = 0.7346
      (textbook value 0.735, Born & Wolf): |Δ| < 0.002.
  R3  Sparrow distance (dip disappears, d²I/dx² = 0 at the midpoint) = 0.4735 λ/NA (textbook ≈ 0.47):
      |Δ| < 0.005.
  R4  photon statistics: maximum-likelihood estimate of the separation from N = 10⁴ clicks on a pixel
      sensor; its spread equals the Cramér–Rao bound (Poisson Fisher information of the same pixel model)
      within 15 %, at d = 1.0 and 0.5 Rayleigh.
  R5  'resolved' (likelihood-ratio two vs one emitter > 3σ, i.e. 2ΔlnL > 9) at d = Rayleigh with 10⁴
      clicks in ≥ 95 % of the images, and the number of clicks needed for 90 % grows as d drops
      (Rayleigh's curse: ∝ d⁻⁴ for d ≪ λ/NA) — the measured exponent between 0.5 and 0.25 Rayleigh
      is between −5 and −3.
Images: lab/out/img_atom_photons.png, lab/out/img_two_atoms.png
"""
import os
import sys
import time
import numpy as np
from scipy.special import j1
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from img_core import Camera, pixel_psf, airy   # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
LAM, NA = 1.0, 0.5
RAY = 0.61 * LAM / NA                          # Rayleigh distance (textbook 0.61)
RAYX = 3.8317059702 / (2 * np.pi) * LAM / NA   # exact first zero
NPIX, PITCH = 24, 0.25                          # sensor 24×24 pixels of λ/4 (object units, M = 1)
BG = 0.02                                       # dark counts per pixel per image


def profile(x, d):
    return airy(np.abs(x - d / 2), NA, LAM) + airy(np.abs(x + d / 2), NA, LAM)


def two_img(d, cx=0.0):
    return 0.5 * (pixel_psf(NPIX, PITCH, NA, LAM, cx - d / 2, 0, sub=6) + pixel_psf(NPIX, PITCH, NA, LAM, cx + d / 2, 0, sub=6))


def main():
    t0 = time.time()
    passed = total = 0
    rng = np.random.default_rng(31)
    # R1
    v = brentq(lambda q: j1(q), 3, 4.5)
    r1 = v / (2 * np.pi) * LAM / NA
    ok = abs(r1 - RAYX) < 1e-4
    passed += ok; total += 1
    print(f"R1 first zero {r1:.5f} λ/NA·NA... = {r1 * NA:.5f} λ/NA (0.6098) -> {'PASS' if ok else 'FAIL'}", flush=True)
    # R2 (Rayleigh with the exact first zero, as Rayleigh defined it)
    d = RAYX
    dip = profile(np.array([0.0]), d)[0] / profile(np.array([d / 2]), d)[0]
    x = np.linspace(-1.5 * d, 1.5 * d, 30001)
    pk = profile(x, d).max()
    dip = profile(np.array([0.0]), d)[0] / pk
    ok = abs(dip - 0.7346) < 0.002
    passed += ok; total += 1
    print(f"R2 dip/peak at the Rayleigh distance: {dip:.4f} (0.7346) -> {'PASS' if ok else 'FAIL'}", flush=True)
    dip061 = profile(np.array([0.0]), RAY)[0] / profile(x * RAY / d, RAY).max()
    print(f"    (with the rounded 0.61 λ/NA: {dip061:.4f})", flush=True)
    # R3 Sparrow
    h = 1e-4

    def curv(dd):
        return profile(np.array([h]), dd)[0] - 2 * profile(np.array([0.0]), dd)[0] + profile(np.array([-h]), dd)[0]
    ds = brentq(curv, 0.3 / NA, 0.6 / NA)
    ok = abs(ds * NA - 0.4735) < 0.005
    passed += ok; total += 1
    print(f"R3 Sparrow distance {ds * NA:.4f} λ/NA (0.4735; = {ds / RAYX:.3f} Rayleigh) -> {'PASS' if ok else 'FAIL'}", flush=True)

    # R4/R5: maximum likelihood on pixel images (separation along x, centre fitted along x)
    seps = np.linspace(0.0, 1.6 * RAYX, 97)
    cxs = np.linspace(-0.15, 0.15, 13)
    models = np.array([[two_img(s, c) for c in cxs] for s in seps])          # [ns, nc, NPIX, NPIX]
    flat = models.reshape(len(seps), len(cxs), -1)

    def fit(k, N):
        mu = N * flat + BG
        ll = (k[None, None, :] * np.log(mu) - mu).sum(-1)                      # [ns, nc]
        i, j = np.unravel_index(np.argmax(ll), ll.shape)
        # parabolic refinement in separation
        if 0 < i < len(seps) - 1:
            a, b, c = ll[i - 1, j], ll[i, j], ll[i + 1, j]
            den = a - 2 * b + c
            sh = 0.5 * (a - c) / den if den < 0 else 0
            s_hat = seps[i] + sh * (seps[1] - seps[0])
        else:
            s_hat = seps[i]
        return s_hat, 2 * (ll.max() - ll[0].max())

    def fisher_sd(dtrue, N):
        e = 1e-4
        m0 = N * two_img(dtrue) + BG
        dm = N * (two_img(dtrue + e) - two_img(dtrue - e)) / (2 * e)
        dc = N * (two_img(dtrue, e) - two_img(dtrue, -e)) / (2 * e)
        F = np.array([[(dm * dm / m0).sum(), (dm * dc / m0).sum()], [(dm * dc / m0).sum(), (dc * dc / m0).sum()]])
        return np.sqrt(np.linalg.inv(F)[0, 0])

    for frac in (1.0, 0.5):
        dtrue = frac * RAYX
        N = 10000
        mean = N * two_img(dtrue) + BG
        est = np.array([fit(rng.poisson(mean).ravel(), N)[0] for _ in range(300)])
        crb = fisher_sd(dtrue, N)
        sd = est.std(ddof=1)
        ok = abs(sd / crb - 1) <= 0.15
        passed += ok; total += 1
        print(f"R4 d = {frac:.2f} Rayleigh, N = {N}: ML separation {est.mean() / RAYX:.4f} ± {sd / RAYX:.4f} Rayleigh,"
              f" Cramér–Rao {crb / RAYX:.4f}  (ratio {sd / crb:.3f}) -> {'PASS' if ok else 'FAIL'}", flush=True)

    # R5: clicks needed for 90 % 'resolved' (2ΔlnL > 9)
    need = {}
    for frac in (1.0, 0.5, 0.25):
        dtrue = frac * RAYX
        mshape = two_img(dtrue)
        prev = None
        for N in [30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000]:
            res = np.array([fit(rng.poisson(N * mshape + BG).ravel(), N)[1] > 9 for _ in range(100)])
            pr = res.mean()
            if pr >= 0.9:
                # log-interpolate between the previous N and this one
                if prev is not None and prev[1] < 0.9:
                    t = (0.9 - prev[1]) / (pr - prev[1])
                    need[frac] = np.exp(np.log(prev[0]) + t * (np.log(N) - np.log(prev[0])))
                else:
                    need[frac] = N
                print(f"    d = {frac:.2f} Rayleigh: P(resolved) = {pr:.2f} at N = {N}; 90% at N ≈ {need[frac]:.0f} clicks", flush=True)
                break
            prev = (N, pr)
        if frac == 1.0:
            p1e4 = np.array([fit(rng.poisson(10000 * mshape + BG).ravel(), 10000)[1] > 9 for _ in range(200)]).mean()
    ok_a = p1e4 >= 0.95
    expo = np.log(need[0.25] / need[0.5]) / np.log(0.25 / 0.5) if 0.25 in need and 0.5 in need else np.nan
    ok_b = -5 <= expo <= -3
    passed += ok_a + ok_b; total += 2
    print(f"R5 P(resolved) at Rayleigh with 10⁴ clicks = {p1e4:.3f} -> {'PASS' if ok_a else 'FAIL'};"
          f" exponent N_90 ∝ d^{expo:.2f} between 0.5 and 0.25 Rayleigh (expect ≈ −4) -> {'PASS' if ok_b else 'FAIL'}", flush=True)
    # diagnostic (added after the first run, NOT a test): ideal-observer prediction from the expected
    # log-likelihood ratio per photon, 2·KL(two‖one); 90 % detection of a χ²₁ threshold 9 needs N·2KL ≈ 18.3
    from scipy.stats import ncx2
    lam90 = brentq(lambda L: ncx2.sf(9, 1, L) - 0.9, 1, 100)
    kl = {}
    for frac in (1.0, 0.5, 0.25, 0.125, 0.0625):
        m2 = two_img(frac * RAYX) + BG / 1e4
        m1 = two_img(0.0) + BG / 1e4
        kl[frac] = 2 * (m2 * np.log(m2 / m1) - m2 + m1).sum()
    print("    ideal-observer N_90 = " + ", ".join(f"{f:g} Ray: {lam90 / kl[f]:.0f}" for f in kl) +
          f"; local exponents: 0.5→0.25 {np.log(kl[0.25] / kl[0.5]) / np.log(2):.2f}, 0.25→0.125"
          f" {np.log(kl[0.125] / kl[0.25]) / np.log(2):.2f}, 0.125→0.0625 {np.log(kl[0.0625] / kl[0.125]) / np.log(2):.2f}"
          f" (→ −4 only for d ≪ λ/NA)", flush=True)
    print(f"tests done {time.time() - t0:.1f}s", flush=True)

    # ---- images ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cam = Camera(qe=1.0, dark_rate=BG, read_noise=0.0)
    one = pixel_psf(NPIX, PITCH, NA, LAM)
    one /= one.sum()
    Ns = [1, 10, 100, 1000, 10000]
    fig, ax = plt.subplots(1, 6, figsize=(12, 2.4))
    for a, N in zip(ax, Ns):
        img = cam.expose(N * one, 1.0, rng)
        a.imshow(img, cmap="gray", interpolation="nearest")
        a.set_title(f"{N} fotoni rivelati", fontsize=9)
        a.axis("off")
    ax[-1].imshow(one, cmap="gray", interpolation="nearest")
    ax[-1].set_title("media (infiniti)", fontsize=9)
    ax[-1].axis("off")
    fig.suptitle(f"un atomo, NA {NA}, pixel λ/4, buio {BG}/pixel: l'immagine nasce dai clic", fontsize=10)
    fig.savefig(os.path.join(OUT, "img_atom_photons.png"), dpi=70, bbox_inches="tight")
    plt.close(fig)
    pil = None
    pa, pb = os.path.join(OUT, "img_pilot_clicks_A.npy"), os.path.join(OUT, "img_pilot_clicks_B.npy")
    if os.path.exists(pa) and os.path.exists(pb):
        pil = np.vstack([np.load(pa), np.load(pb)])
    fig, ax = plt.subplots(2, 4 if pil is not None else 3, figsize=(10, 5.2))
    for col, frac in enumerate((1.0, 0.75, 0.5)):
        for row, N in enumerate((10000, 300)):
            img = cam.expose(N * two_img(frac * RAYX), 1.0, rng)
            ax[row, col].imshow(img, cmap="gray", interpolation="nearest")
            ax[row, col].set_title(f"d = {frac:.2f}·Rayleigh, {N} clic", fontsize=9)
            ax[row, col].axis("off")
    if pil is not None:
        na_p = 8.0 / np.hypot(8.0, 40.0)           # lens of img_pilot.py (emitter B imaged at x = -d)
        ray_p = 0.61 / na_p
        ex = np.linspace(-ray_p / 2 - 2 * ray_p, -ray_p / 2 + 2 * ray_p, 41)
        ey = np.linspace(-2 * ray_p, 2 * ray_p, 41)
        hh, _, _ = np.histogram2d(pil[:, 1], pil[:, 0], [ey, ex])
        ax[0, 3].imshow(hh, cmap="gray", interpolation="nearest")
        ax[0, 3].set_title(f"regola a palline, {len(pil)} clic\nNA {na_p:.2f}, d = Rayleigh", fontsize=9)
        ax[0, 3].axis("off")
        ax[1, 3].axis("off")
    fig.suptitle("due atomi: al limite di Rayleigh, sotto, e con pochi fotoni", fontsize=10)
    fig.savefig(os.path.join(OUT, "img_two_atoms.png"), dpi=70, bbox_inches="tight")
    plt.close(fig)
    for f in ("img_atom_photons.png", "img_two_atoms.png"):
        print(f"wrote lab/out/{f} ({os.path.getsize(os.path.join(OUT, f)) / 1024:.0f} kB)", flush=True)
    print(f"time {time.time() - t0:.1f}s", flush=True)
    print(f"RESULT: {passed}/{total} PASS", flush=True)


if __name__ == "__main__":
    main()

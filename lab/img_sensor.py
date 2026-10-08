"""Task 1 — the sensor: photon clicks per pixel, and image quality vs photons per pixel.

  python lab/img_sensor.py

Tests (tolerances written before the run):
  T1  ideal camera (QE=1, no dark, no read noise), uniform light of N photons/pixel, N = 1 … 10⁴:
      measured SNR = mean/std over 400 000 pixels equals √N within 3σ (σ of the SNR estimator from
      bootstrap-free formula: rel. error of std ≈ 1/√(2M) for Poisson with N ≥ 10; for N = 1 checked
      directly against the Poisson value 1).
  T2  Fano factor (variance / mean) = 1 for every N (Poisson clicks): |F − 1| < 3·√(2/M)(+kurtosis term).
  T3  real-ish camera (QE 0.8, dark 0.05 e/pixel per exposure, read noise 1.5 e): measured SNR equals
      QE·N / √(QE·N + dark + read²) within 3σ — the read noise dominates below a few photons/pixel.
  T4  clicks from LANDED photons (binomial QE coin, Camera.clicks) give the same SNR as the Poisson chain
      (Camera.expose) when the landed counts are themselves Poisson: |ΔSNR| within 3σ.
  T5  image quality of a picture: relative rms error of a whole picture (smooth test pattern, Airy of a point
      emitter on 64×64 pixels) vs total photons scales as 1/√N_tot: fitted slope −0.5 ± 0.02.
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from img_core import Camera, pixel_psf   # noqa: E402

M = 400_000
NS = [1, 10, 100, 1000, 10000]


def snr_stats(x):
    m, s = x.mean(), x.std(ddof=1)
    n = x.size
    # standard error of SNR = m/s via delta method with the sample's own 3rd/4th moments
    c = x - m
    m3, m4 = (c ** 3).mean(), (c ** 4).mean()
    var_m = s ** 2 / n
    var_s2 = (m4 - s ** 4) / n
    cov = m3 / n
    g = np.array([1 / s, -m / (2 * s ** 3)])
    V = np.array([[var_m, cov], [cov, var_s2]])
    return m / s, float(np.sqrt(g @ V @ g)), s * s / m if m > 0 else np.nan


def main():
    rng = np.random.default_rng(2026)
    t0 = time.time()
    passed, total = 0, 0
    print("N/pixel | ideal SNR (meas ± σ) | √N | Fano | camera SNR (meas ± σ) | theory", flush=True)
    ideal = Camera()
    real = Camera(qe=0.8, dark_rate=0.05, read_noise=1.5)          # dark: per pixel per exposure (exposure=1)
    rows = []
    for N in NS:
        x = ideal.expose(np.full(M, float(N)), 1.0, rng)
        snr, es, fano = snr_stats(x)
        ok1 = abs(snr - np.sqrt(N)) <= 3 * es
        ok2 = abs(fano - 1) <= 3 * np.sqrt(2 / M) * max(1, np.sqrt(1 + 1 / (2 * N)))
        y = real.expose(np.full(M, float(N)), 1.0, rng)
        snr_r, es_r, _ = snr_stats(y - real.dark_rate * 1.0)      # signal = after subtracting the mean dark (calibration)
        th = real.snr_theory(N, 1.0)
        ok3 = abs(snr_r - th) <= 3 * es_r
        landed = rng.poisson(N, M)
        z = real.clicks(landed, 1.0, rng)
        snr_z, es_z, _ = snr_stats(z - real.dark_rate * 1.0)
        ok4 = abs(snr_z - snr_r) <= 3 * np.hypot(es_z, es_r)
        for ok in (ok1, ok2, ok3, ok4):
            passed += ok
            total += 1
        rows.append((N, snr, es, fano, snr_r, es_r, th, snr_z))
        print(f"{N:6d} | {snr:9.4f} ± {es:.4f} | {np.sqrt(N):9.4f} | {fano:.4f} | {snr_r:9.4f} ± {es_r:.4f} | {th:9.4f}"
              f" | landed+QE coin {snr_z:.4f}  [T1 {'PASS' if ok1 else 'FAIL'} T2 {'PASS' if ok2 else 'FAIL'}"
              f" T3 {'PASS' if ok3 else 'FAIL'} T4 {'PASS' if ok4 else 'FAIL'}]", flush=True)
    # crossover where read noise matters: QE N = read² + dark
    sx = (1.5 ** 2 + 0.05) / 3
    print(f"read noise 1.5 e + dark 0.05: SNR is half of the ideal √s at s = QE·N = {sx:.3f} photo-electrons"
          f" (s/√(s+d+r²) = √s/2 -> s = (d+r²)/3); it is within 10% of √s only for s > {(1.5 ** 2 + 0.05) / (1 / 0.9 ** 2 - 1):.1f}", flush=True)
    # T5 picture error vs total photons
    psf = pixel_psf(64, 0.25, 0.5, 1.0)          # pixel = λ/4 at NA 0.5: ~2.4 pixels per Airy radius
    psf /= psf.sum()
    Ntot = np.array([1e2, 1e3, 1e4, 1e5, 1e6])
    errs = []
    for Nt in Ntot:
        e = []
        for k in range(40):
            img = ideal.expose(Nt * psf, 1.0, rng)
            e.append(np.sqrt(((img - Nt * psf) ** 2).sum() / ((Nt * psf) ** 2).sum()))
        errs.append(np.mean(e))
    slope = np.polyfit(np.log(Ntot), np.log(errs), 1)[0]
    ok5 = abs(slope + 0.5) <= 0.02
    passed += ok5
    total += 1
    print("picture rms error vs total photons: " + ", ".join(f"N={int(n)}: {100 * e:.2f}%" for n, e in zip(Ntot, errs)),
          flush=True)
    print(f"T5 slope {slope:.4f} (theory -0.5) -> {'PASS' if ok5 else 'FAIL'}", flush=True)
    print(f"time {time.time() - t0:.1f}s", flush=True)
    print(f"RESULT: {passed}/{total} PASS", flush=True)


if __name__ == "__main__":
    main()

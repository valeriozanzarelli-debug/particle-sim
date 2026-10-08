"""Task 3 — how many photons do we really see from ONE hydrogen atom?

  python lab/img_hydrogen.py

Atom light rules used here (all from EXACT non-relativistic QM of hydrogen, solved numerically in
img_core.hydrogen / polarizability — NOT material tables):
    >>> placeholder: replace with our balls' value (docs/MOLECOLE.md does not exist yet) <<<
  * line energy 1s→2p (Lyman-α), dipole matrix element, Einstein A = Γ;
  * dynamic polarizability α(ω) of 1s → Rayleigh cross-section σ_R(λ) = (8π/3) k⁴ α²;
  * resonant cross-section and saturation intensity of the 1s→2p line treated as a closed two-level
    transition (fine and hyperfine structure ignored: g2/g1 = 3, the spin just doubles both levels);
  * angular pattern: electric dipole, dP/dΩ = 3/(8π) sin²ψ.
Light sources: black body (Planck) with the Sun's nominal temperature and angular size (scene parameters,
IAU 2015 nominal values), and an ideal laser at Lyman-α.  Optics: ideal lens of numerical aperture NA in
vacuum (transmission 1).  Camera: QE and dark rate are DEVICE parameters (scanned).

Tests:
  H1  grid hydrogen vs closed-form exact QM: E1s = −μ/2, E2p = −μ/8, ⟨1s|z|2p0⟩ = 128√2/243 /μ a0,
      α(0) = 4.5/μ³ a0³: relative error < 1e-6 (energies, dipole), < 1e-4 (α).
  H2  dipole collection fraction checked against the closed forms
      perp (4 − 3c − c³)/8, par (2 − 3c + c³)/4, c = cos θ (numerical integral vs formula, |Δ| < 1e-6).
  H3  Planck: total power of the black body at the scene distance equals σT⁴·(R/d)² (Stefan–Boltzmann
      from constants) within 1e-3.
"""
import os
import sys
import time
import numpy as np
from scipy.stats import poisson
from scipy.special import j0, j1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import img_core as IC   # noqa: E402

# scene parameters (geometry of the light source, not material properties)
T_SUN = 5772.0            # K, IAU 2015 nominal solar effective temperature
R_SUN = 6.957e8           # m, IAU 2015 nominal solar radius
D_AU = 1.495978707e11     # m, astronomical unit (exact by definition)


def main():
    t0 = time.time()
    passed = total = 0
    hy = IC.hydrogen()
    mu = hy["mu"]
    e1 = abs(hy["E1s"] / hy["exact_E1s"] - 1)
    e2 = abs(hy["E2p"] / hy["exact_E2"] - 1)
    ez = abs(hy["dip_z"] / hy["exact_z"] - 1)
    a0 = IC.polarizability(hy, [0.0])[0]
    ea = abs(a0 / (4.5 / mu ** 3) - 1)
    ok = e1 < 1e-6 and e2 < 1e-6 and ez < 1e-6 and ea < 1e-4
    passed += ok; total += 1
    lam = hy["lam_lya"]
    G = hy["A21"]
    print(f"H1 hydrogen on the grid (reduced mass μ = {mu:.7f} m_e): rel. errors E1s {e1:.1e}, E2p {e2:.1e},"
          f" <z> {ez:.1e}, α(0) {ea:.1e} -> {'PASS' if ok else 'FAIL'}", flush=True)
    print(f"   Lyman-α: λ = {lam * 1e9:.4f} nm, hν = {IC.H * IC.C / lam / IC.QE_CHARGE:.4f} eV,"
          f" A = Γ = {G:.5e} s⁻¹ (τ = {1e9 / G:.4f} ns, Γ/2π = {G / 2 / np.pi / 1e6:.2f} MHz), α(0) = {a0:.4f} a0³", flush=True)

    # H2 collection fractions vs closed forms
    errs = []
    for na in (0.3, 0.6, 0.9):
        th = np.arcsin(na)
        c = np.cos(th)
        perp_cf = (4 - 3 * c - c ** 3) / 8       # ∫ (3/8π)(1 − sin²θ cos²φ) dΩ over the cone
        par_cf = (2 - 3 * c + c ** 3) / 4
        errs += [abs(IC.collection_fraction(na, "perp") - perp_cf), abs(IC.collection_fraction(na, "par") - par_cf)]
    ok = max(errs) < 1e-6
    passed += ok; total += 1
    print(f"H2 dipole collection integral vs closed forms (4−3c−c³)/8 and (2−3c+c³)/4: max |Δ| {max(errs):.1e}"
          f" -> {'PASS' if ok else 'FAIL'}", flush=True)

    # ---------- (a) white light: black body, Rayleigh scattering ----------
    Om = np.pi * (R_SUN / D_AU) ** 2                        # solid angle of the disc
    lg = np.geomspace(150e-9, 30e-6, 1500)
    Lph = IC.planck_photon_spectral_radiance(lg, T_SUN)     # photons/(s m² sr m)
    flux_ph = Lph * Om                                      # photons/(s m² m) at the atom
    power = np.trapezoid(flux_ph * IC.H * IC.C / lg, lg)
    sb = 2 * np.pi ** 5 * IC.KB ** 4 / (15 * IC.C ** 2 * IC.H ** 3)
    p_sb = sb * T_SUN ** 4 * (R_SUN / D_AU) ** 2
    ok = abs(power / p_sb - 1) < 1e-3
    passed += ok; total += 1
    print(f"H3 black body {T_SUN:.0f} K seen at 1 AU: {power:.1f} W/m² (150 nm–30 µm) vs σT⁴(R/d)² = {p_sb:.1f} W/m²"
          f" -> {'PASS' if ok else 'FAIL'}", flush=True)
    w_au = 2 * np.pi * IC.C_AU * IC.A0 / lg                 # photon energy in Hartree
    alpha = IC.polarizability(hy, w_au)
    sig = IC.rayleigh_sigma(alpha, lg)
    rate_R = np.trapezoid(sig * flux_ph, lg)
    nph = np.trapezoid(flux_ph, lg)
    i500 = np.argmin(abs(lg - 500e-9))
    print(f"(a) sunlight-like black body: {nph:.3e} photons/(s·m²) on the atom;"
          f" σ_R(500 nm) = {sig[i500]:.3e} m² (α = {alpha[i500]:.3f} a0³), σ_R(200 nm) = {sig[np.argmin(abs(lg - 200e-9))]:.3e} m²", flush=True)
    print(f"    Rayleigh-scattered photons per second by ONE H atom: {rate_R:.3e} s⁻¹"
          f"  (one photon every {1 / rate_R / 3.156e7:.1f} years)", flush=True)
    nbar = IC.occupation(lam, T_SUN)
    rate_abs = 3 * G * nbar * Om / (4 * np.pi)               # B12·u: (g2/g1)·A·n̄·(Ω/4π)
    print(f"    resonant Lyman-α absorption of the same black body: n̄(Lyα, {T_SUN:.0f} K) = {nbar:.3e},"
          f" rate = {rate_abs:.3e} s⁻¹", flush=True)

    # ---------- (b) resonant Lyman-α laser ----------
    sig0 = 3 * lam ** 2 / (2 * np.pi)
    hnu = IC.H * IC.C / lam
    isat = hnu * G / (2 * sig0)
    print(f"(b) Lyα laser (closed two-level 1s–2p): σ0 = 3λ²/2π = {sig0:.3e} m², I_sat = ħωΓ/(2σ0) = {isat:.2f} W/m²"
          f" = {isat / 10:.3f} mW/cm²; power for s = 1 in a beam of waist 100 µm: {isat * np.pi * (100e-6) ** 2 / 2 * 1e6:.2f} µW", flush=True)
    for s in (0.1, 1, 10, 1e3):
        print(f"    s = {s:g}: scattering rate Γ/2·s/(1+s) = {G / 2 * s / (1 + s):.3e} photons/s", flush=True)
    # kinematics of one H atom scattering Lyα (exact: momentum of the photon)
    m_h = IC.MP + IC.ME
    k = 2 * np.pi / lam
    vr = IC.HBAR * k / m_h
    er = (IC.HBAR * k) ** 2 / (2 * m_h)
    td = IC.HBAR * G / (2 * IC.KB)
    vrms = np.sqrt(IC.KB * td / m_h)
    print(f"    recoil per photon {vr:.3f} m/s (E_r/k_B = {er / IC.KB * 1e3:.3f} mK); Doppler shift per recoil kv = {k * vr / G:.3f} Γ"
          f" -> a single beam pushes the atom off resonance after ≈ {0.5 / (k * vr / G):.1f} photons; Doppler limit"
          f" T_D = ħΓ/2k_B = {td * 1e3:.2f} mK, v_rms = {vrms:.2f} m/s", flush=True)

    # ---------- (c) collection, detector, exposure ----------
    print("(c) clicks per second = rate × collection(NA, dipole) × QE   [ideal optics, transmission 1]", flush=True)
    rate_sat = G / 2 * 1 / 2          # s = 1
    qe = 0.5
    for na in (0.5, 0.8, 0.95):
        fp, fpar, fi = (IC.collection_fraction(na, o) for o in ("perp", "par", "iso"))
        res = 0.61 * lam / na
        fwhm = 0.5145 * lam / na
        print(f"    NA {na}: collected fraction ⊥ dipole {fp:.4f}, ∥ dipole {fpar:.4f}, isotropic {fi:.4f};"
              f" resolution 0.61λ/NA = {res * 1e9:.1f} nm, FWHM {fwhm * 1e9:.1f} nm;"
              f" clicks/s at s=1, QE {qe}: Lyα {rate_sat * fi * qe:.3e}, sunlight-Rayleigh {rate_R * fi * qe:.3e}", flush=True)
    # time the atom stays inside one resolution element at the Doppler temperature (free, untrapped)
    na = 0.8
    res = 0.61 * lam / na
    t_stay = res / vrms
    print(f"    free atom at the Doppler limit crosses one resolution element ({res * 1e9:.0f} nm) in {t_stay * 1e9:.1f} ns:"
          f" detected photons in that time at s=1, NA 0.8, QE 0.5: {rate_sat * IC.collection_fraction(na, 'iso') * qe * t_stay:.2f}", flush=True)

    # exposure for a recognisable image: detection test on a 3×3 window (pixel = λ/(4NA), Nyquist)
    def frac3x3(na):
        sub = 60
        p = 0.25 * lam / na
        s = (np.arange(3 * sub) + 0.5) / sub * p - 1.5 * p
        X, Y = np.meshgrid(s, s)
        v = 2 * np.pi * na * np.hypot(X, Y) / lam
        I = np.where(v > 1e-12, (2 * j1(np.maximum(v, 1e-12)) / np.maximum(v, 1e-12)) ** 2, 1.0)
        return I.sum() * (p / sub) ** 2 / (lam ** 2 / (np.pi * na ** 2))

    f3 = frac3x3(0.8)
    print(f"    fraction of the photons in a 3×3 window of Nyquist pixels (λ/4NA), NA 0.8: {f3:.4f}", flush=True)

    def exposure_needed(click_rate, dark):
        """shortest t with P(detect) ≥ 0.95 and false alarm ≤ 1e-3 per window (dark: e/pixel/s, 9 pixels)."""
        for t in np.geomspace(1e-9, 1e16, 2501):
            b = 9 * dark * t
            kth = int(poisson.isf(1e-3, b)) + 1 if b > 0 else 1
            if poisson.sf(kth - 1, b + f3 * click_rate * t) >= 0.95:
                return t, click_rate * t, kth
        return np.inf, np.inf, None

    fi = IC.collection_fraction(0.8, "iso")
    for label, r in (("Lyα laser s=1", rate_sat), ("sunlight Rayleigh", rate_R)):
        for dark in (0.0, 1e-3, 1.0):
            t, nd, kth = exposure_needed(r * fi * qe, dark)
            if np.isfinite(t):
                print(f"    exposure for a recognisable single-atom spot [{label}, NA 0.8, QE {qe}, dark {dark:g} e/pix/s]:"
                      f" {t:.3e} s, {nd:.1f} detected photons (threshold {kth} counts in 3×3)", flush=True)
            else:
                sw, bw = f3 * r * fi * qe, 9 * dark
                tg = (3.09 * np.sqrt(bw) + 1.645 * np.sqrt(bw + sw)) ** 2 / sw ** 2     # Gaussian limit, S ≪ B
                print(f"    [{label}, dark {dark:g} e/pix/s]: signal {sw:.2e} clicks/s in the window vs dark {bw:.2e} counts/s"
                      f" -> exposure (z = 3.09 false alarm, 1.645 detection) ≈ {tg:.2e} s = {tg / 3.156e7:.1e} years", flush=True)
    # localisation: σ_PSF of the Airy spot (Gaussian fit width ≈ 0.21 λ/NA) / √N
    nloc = (0.21 / (0.1 * 0.61)) ** 2
    print(f"    to locate the atom to 1/10 of the resolution (σ = 0.061 λ/NA) one needs N ≈ (0.21/0.061)² = {nloc:.0f} photons"
          f" -> Lyα s=1, NA 0.8, QE {qe}: {nloc / (rate_sat * fi * qe) * 1e9:.0f} ns", flush=True)

    # ---------- comparison with real single-atom imaging (Rb-87, 780 nm) ----------
    print("comparison (numbers from the papers; arithmetic here):", flush=True)
    sh = 60e3 * 0.9 * 0.15 * 0.71
    print(f"    Sherson 2010: 60 kHz scattering × 0.9 s × solid angle 15% × transmission 71% = {sh:.0f} photons reaching the"
          f" camera (paper: ~5000 collected), × QE 0.85 = {sh * 0.85:.0f}; NA 0.68 at 780 nm:"
          f" our Airy FWHM {0.5145 * 780 / 0.68:.0f} nm (paper ≈700 nm measured), 0.61λ/NA = {0.61 * 780 / 0.68:.0f} nm", flush=True)
    print(f"    our collection fraction at NA 0.68, isotropic: {IC.collection_fraction(0.68, 'iso'):.4f}"
          f" (paper: 15%)", flush=True)
    print(f"    Bakr 2009: NA 0.8 at 780 nm: our Airy FWHM {0.5145 * 780 / 0.8:.0f} nm (paper: diffraction limit 500 nm,"
          f" measured 570/630 nm); up to 2·10⁴ photons/atom/s collected, exposure 0.2–1 s -> {2e4 * 0.2:.0f}–{2e4 * 1:.0f} per atom", flush=True)
    print(f"time {time.time() - t0:.1f}s", flush=True)
    print(f"RESULT: {passed}/{total} PASS", flush=True)


if __name__ == "__main__":
    main()

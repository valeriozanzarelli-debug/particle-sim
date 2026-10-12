"""h2vib_prereg.py - pre-registration of the H2 vibration anomaly study (B0-B7), written BEFORE the runs.
Writes lab/out/h2vib_prereg.json."""
import json, os, subprocess, datetime

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
P = {
    "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head},
    "anomaly": "docs/MOLECOLE.md 7.2: H2 omega_e = 4339 +- 18 cm-1 (harmonic, from the curvature of the fitted DMC "
               "curve) vs 4401.21 (Huber & Herzberg), -1.4 %; nu(0->1) = 4129 +- 16 vs 4161.17, -0.8 %.",
    "data": "the 18 DMC points (and the per-run values) in assets/molecole/risultati.txt; no new DMC unless declared "
            "later. Same fit (E_inf + sum_{k=1..5} c_k exp(-k beta R), E_inf = 2 E_H(guided), beta grid 0.2..1.36 "
            "step 0.04, weighted) and same finite-difference vibrational solver (0.4..14 bohr, 6000 points) "
            "re-implemented in Python; mu = (m_p + m_e)/2 as in misura.rs.",
    "reference_curve": "exact Born-Oppenheimer potential, Pachucki PRA 82, 032509 (2010) Table IV "
            "(lab/h2vib_pachucki.py), cubic spline in R. Seen before prereg: the table, and by eye that the DMC-minus-"
            "exact residuals are +1e-4..+4e-4 Eh near the minimum.",
    "checklist": {
        "B0_reproduce": "Python pipeline reproduces the Rust numbers (omega_e 4339, nu01 4129, Re 0.7420 A) within "
                        "1 cm-1 / 2e-4 A. Prerequisite.",
        "B1_harmonic_vs_anharmonic": "on the exact BO curve: harmonic omega_e (curvature) and Dunham-like "
            "omega_e, omega_e x_e, omega_e y_e fitted to G(v), v = 0..5; check experiment self-consistency "
            "nu01 = omega_e - 2 omega_e x_e + (13/4) omega_e y_e. CAUSE if harmonic omega_e(exact BO) differs from "
            "4401.21 by > 1%.",
        "B2_reduced_mass": "mu = m_p/2 vs (m_p+m_e)/2. CAUSE if change >= 0.5% (expected ~0.03%).",
        "B3_BO_corrections": "exact-BO nu01 (from the Pachucki curve, computed here) vs experiment 4161.17 = adiabatic + "
            "nonadiabatic + relativistic + QED (literature: Komasa et al. JCTC 7, 3105 (2011): a few cm-1). CAUSE if "
            "|BO - exp| >= 0.5%.",
        "B4_statistics": "synthetic: exact BO energies at our 18 R + Gaussian noise (our error bars x 2.3), 2000 "
            "replicas, identical fit -> distribution of omega_e, nu01 (bias, spread, p-value of our deviation). Plus "
            "residuals DMC - exact, chi2/point with bare error bars. Verdict 'statistical fluke' if deviation within "
            "2.5 sigma of the synthetic spread AND residual chi2/point <= 4.",
        "B5_fit_model_and_range": "same fit on EXACT energies at our 18 R (no noise): model bias. On DMC data: E_inf "
            "free; without R = 0.9 and 5.5; local polynomials deg 2,3,4 over windows 0.2,0.3,0.35,0.45 bohr. CAUSE "
            "(or part) if model bias on exact data >= 0.5% in omega_e.",
        "B6_grid_step": "solver with 1500, 3000, 6000, 12000 points; range 0.4..14 vs 0.6..8. CAUSE if change >= 0.1%.",
        "B7_time_step": "curvature of the residual DMC - exact (weighted quadratic over 1.1..1.9 bohr) converted to an "
            "omega_e shift. CAUSE if |shift| >= 0.5% and > 2 sigma.",
    },
    "corrected_number": "final omega_e, nu01 = DMC curve with the model/range that B5 shows unbiased on exact data, "
            "errors from B4-type bootstrap; compared with experiment after the literature BO->exp correction (B3). "
            "'Resolved' if |corrected - exp| <= 2 sigma for both omega_e and nu01.",
}
json.dump(P, open(os.path.join(OUT, "h2vib_prereg.json"), "w"), indent=1)
print("written", os.path.join(OUT, "h2vib_prereg.json"))

"""glue_prereg.py - pre-registration of the glueball study (G1-G5), written BEFORE the production runs.
Writes lab/out/glue_prereg.json (run plan read by glue_run.py + criteria graded by glue_verdict.py)."""
import json, os, subprocess, datetime

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()

LIT = {
    "AT2017": {"ref": "A. Athenodorou, M. Teper, 'SU(N) gauge theories in 2+1 dimensions: glueball spectra and "
                      "k-string tensions', JHEP 02 (2017) 015, arXiv:1609.03873, Tables A7 and B3 (SU(2))",
               "continuum_m_over_sqrt_sigma": {"0++": [4.7367, 0.0055], "0++*": [6.861, 0.014],
                                               "0++**": [8.382, 0.014], "2++": [7.762, 0.010], "2-+": [7.795, 0.012]},
               "lattice": {"6": {"a_sqrt_sigma": 0.25098, "am0": 1.1908, "lattice_sigma": "25^2x20", "lattice_m": "18^2x20"},
                           "9": {"a_sqrt_sigma": 0.16075, "am0": 0.7624, "lattice_sigma": "36^2x32", "lattice_m": "22^2x32"},
                           "12": {"a_sqrt_sigma": 0.11847, "am0": 0.5594, "lattice_sigma": "50^2x40", "lattice_m": "34^2x40"}}},
    "Teper1998": {"ref": "M. Teper, 'SU(N) gauge theories in 2+1 dimensions', Phys. Rev. D 59 (1998) 014512, "
                         "hep-lat/9804008, Table 28 (SU(2) continuum)",
                  "continuum_m_over_sqrt_sigma": {"0++": [4.718, 0.043], "0++*": [6.83, 0.10], "2++": [7.82, 0.14]},
                  "sqrt_sigma_over_g2": [0.3353, 0.0018],
                  "a_sqrt_sigma": {"6": 0.2538, "9": 0.1622, "12": 0.1179}},
}
R2_LIT = 7.762 / 4.7367      # m(2++)/m(0++)  = 1.6387
R0S_LIT = 6.861 / 4.7367     # m(0++*)/m(0++) = 1.4485

common = dict(n_or=3, ape_alpha=0.5, creutz_Rmax=6)
P = {
    "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head,
             "numba": "numba 0.68.0 was already installed and is used for the update and smearing kernels "
                      "(single thread); results identical to the numpy reference to 1e-15 (tested)."},
    "question": "Valerio: matter made only of light that sticks to itself. Pure SU(2) gauge theory in 2+1 D (no quarks, "
                "no mass term anywhere): do the self-interacting massless gauge fields form massive lumps (glueballs), "
                "and are the parameter-free numbers m(2++)/m(0++), m(0++*)/m(0++), m/sqrt(sigma) the known ones?",
    "rule": "Wilson action S = beta sum_P (1 - 1/2 tr U_P), beta = 4/(a g^2), SU(2) links as unit quaternions, periodic "
            "T x L x L lattice. One compound sweep = 1 Kennedy-Pendleton heat-bath + 3 over-relaxation sweeps "
            "(checkerboard). Cold start, n_therm compound sweeps discarded, then one measurement every `sep` compound sweeps.",
    "literature_checked_before_prereg": LIT,
    "targets": {"R2 = m(2++)/m(0++)": round(R2_LIT, 4), "R0s = m(0++*)/m(0++)": round(R0S_LIT, 4),
                "note": "from AT2017 continuum values; Teper1998 gives 7.82/4.718 = 1.657 and 6.83/4.718 = 1.448"},
    "pilot_seen_before_prereg": "declared: one engine test at beta = 6, 20^2 x 24, 400 measurements (seed 99, files "
                "pilot_b6_*, not used for any graded number). It gave am(0++) = 1.13(13), torelon am = 1.05(10), Creutz "
                "chi(4,4) = 0.080(3); it showed that a plateau rule based only on 1-sigma agreement accepts noise, so the "
                "rule below adds a 20% relative-error requirement. No production data existed when this file was written.",
    "runs": {
        "6":  dict(common, T=24, L=20, seed=6, n_therm=300, sep=2, n_meas=8000, creutz_every=20,
                   smear_levels=[2, 6, 12], sq_sizes=[1, 2, 3, 4], rect_sizes=[1, 2]),
        "9":  dict(common, T=32, L=28, seed=9, n_therm=300, sep=2, n_meas=6000, creutz_every=20,
                   smear_levels=[4, 10, 20], sq_sizes=[1, 2, 3, 4, 6], rect_sizes=[1, 2, 3]),
        "12": dict(common, T=32, L=40, seed=12, n_therm=300, sep=2, n_meas=5000, creutz_every=20,
                   smear_levels=[6, 14, 28], sq_sizes=[2, 3, 4, 6, 8], rect_sizes=[1, 2, 3, 4]),
    },
    "volumes_reasoning": "L sqrt(sigma) ~ 5 at every beta (AT2017 used ~3.5-4 for the mass gap); L chosen so that the "
                "lightest two-torelon state 2 m_P(L) lies above the 2++ mass (beta 6: 2.5 vs 1.96; 9: 1.41 vs 1.25; "
                "12: 1.09 vs 0.92 in lattice units, from AT2017 numbers).",
    "budget_rule": "if wall-clock forces stopping before n_meas, the analysis uses what exists and the actual N is reported.",
    "operators": "spatial Wilson loops built from APE-smeared spatial links (U <- Proj[U + alpha * 2 spatial staples], "
                 "alpha 0.5, levels = cumulative smearing steps), summed over each time slice (zero momentum). "
                 "0++ (A1): squares l x l (sq_sizes) and rectangles 2l x l + l x 2l (rect_sizes), at every level. "
                 "2++ (B1 of the square group, P = + since rectangles are reflection-symmetric): 2l x l - l x 2l. "
                 "Torelon: smeared spatial Polyakov loops sum_y P_x + sum_x P_y. "
                 "Momentum p = 2 pi/L along x and y for the A1 operators (G4).",
    "analysis_rule": {
        "correlators": "C_ij(t) = <O_i(t0+t) O_j(t0)> - <O_i><O_j>, averaged over t0, symmetrised t <-> T-t and i <-> j.",
        "errors": "jackknife over 40 contiguous bins of measurements.",
        "variational": "normalise to C_ii(0)=1; keep eigen-directions of C(0) with eigenvalue >= 1e-3 x max; GEVP "
                       "C(1) v = lambda C(0) v (t0 = 0, td = 1, Teper's choice); vectors from the full sample; state 0 = "
                       "largest lambda (ground state), state 1 = next (0++* in A1).",
        "effective_mass": "cosh form: C(t)/C(t+1) = cosh(m(t-T/2))/cosh(m(t+1-T/2)).",
        "plateau_rule": "t_sig = largest t such that C(t')/dC(t') >= 3 for all 1 <= t' <= t. t_min = first t in "
                        "[1, t_sig-1] such that m_eff(t) and m_eff(t+1) both have relative error <= 20% and agree within "
                        "1 sigma (errors in quadrature). If none: t_min = 1 (flag 'no_plateau_fallback_t1'). Fit range "
                        "[t_min, max(t_sig, t_min+1)], uncorrelated single-cosh fit, jackknife error.",
        "string_tension": "torelon ground-state mass m_P (P channel, same rule); sigma from the Nambu-Goto formula "
                          "m_P = sigma L sqrt(1 - pi/(3 sigma L^2)).",
        "creutz": "chi(R,R) = -ln[W(R,R) W(R-1,R-1)/W(R,R-1)^2] from unsmeared loops (all 6 oriented planes), R = 2..6, "
                  "jackknife over 20 bins.",
        "continuum": "m/sqrt(sigma) and the ratios fitted as c + d (a sqrt(sigma))^2 over the three beta (weighted).",
    },
    "criteria": {
        "G1a_sigma_vs_literature": "a sqrt(sigma) (torelon) within 3% of AT2017 at each beta (the AT2017 statistical "
                                   "errors are < 0.5%). PASS if all three.",
        "G1b_creutz": "at each beta, at the largest R with relative error <= 10%: chi(R,R) - pi/(24 (R-1/2)^2) within 15% "
                      "of sigma a^2 (torelon). PASS if all three.",
        "G2a_mass_gap_vs_literature": "am(0++) within 5% of AT2017 a m_G at each beta. PASS if all three.",
        "G2b_R2_finest": "R2 at beta = 12: |R2 - 1.6387| <= 2 sigma_stat + 0.05 x 1.6387. PRIMARY.",
        "G2c_R0s_finest": "R0s at beta = 12: |R0s - 1.4485| <= 2 sigma_stat + 0.07 x 1.4485. PRIMARY.",
        "G2d_m_over_sqrt_sigma_finest": "at beta = 12: m(0++)/sqrt(sigma) within 2 sigma + 4% of 4.7367; m(2++)/sqrt(sigma) "
                                        "within 2 sigma + 6% of 7.762; m(0++*)/sqrt(sigma) within 2 sigma + 8% of 6.861.",
        "G2e_continuum": "continuum fits (3 beta): R2, R0s, m(0++)/sqrt(sigma), m(2++)/sqrt(sigma) with the same tolerances "
                         "as G2b-d. Secondary (graded, reported separately).",
        "G2f_statistics_sanity": "a G2 PASS needs sigma_stat <= 8% of the value (otherwise 'INCONCLUSIVE', not PASS).",
        "G3_no_mass_parameter": "the only input is beta = 4/(a g^2). beta * a sqrt(sigma) at 6, 9, 12 decreasing towards "
                                "the continuum; linear fit in 1/beta extrapolated to 1/beta = 0 gives sqrt(sigma)/g^2 = "
                                "(beta a sqrt(sigma))/4 within 5% of 0.3353 (Teper1998). PASS if so.",
        "G4_mass_is_energy": "dispersion of the 0++ glueball: c^2 = (E(p)^2 - m^2)/p^2 with p = 2 pi/L, at beta = 9 and 12. "
                             "PASS if |c^2 - 1| <= 2 sigma at both and sigma <= 0.3 (else INCONCLUSIVE).",
        "G5": "figure lab/out/glue_figure.png (not graded).",
    },
    "no_fitted_parameters": "beta, L, T are inputs; there is no mass parameter in the action. The only fits are the "
                            "exponential decays and the continuum extrapolations.",
}
json.dump(P, open(os.path.join(OUT, "glue_prereg.json"), "w"), indent=1)
print("written", os.path.join(OUT, "glue_prereg.json"))

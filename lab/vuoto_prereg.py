"""vuoto_prereg.py - pre-registration of 'Vuoto che polarizza' (prefix vuoto_), written BEFORE the production runs.
Writes lab/out/vuoto_prereg.json."""
import json, os, subprocess, datetime

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
P = {
    "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head,
             "study": "docs/VUOTO_POLARIZZA.md", "engine": "lab/vuoto_core.py (rule = lab/zz_core.py, frozen)"},
    "question": "Does the filled Dirac sea of the zig-zag rule polarize in a static Coulomb potential like the QED "
                "vacuum (Uehling), with the charge-renormalization log 2 alpha/(3 pi), and how much does it shift the "
                "levels (H 2S ~ -27 MHz, Nd59+ 1s)?",
    "method": "exact first-order (in Z alpha = Uehling order) static response of the Floquet sea of the rule. "
              "Sea = all Floquet states with quasi-energy in (-pi, 0) per Floquet period (1-D: one tick; 3-D: two "
              "ticks, the rule alternates xyz/zyx). Potential = the rule's phase e^{-iV} after each coin. dP from the "
              "Riesz formula (occupied-empty pairs only), density measured stroboscopically at the end of the "
              "period (slice 0). k-space integrals (Bloch matrices of zz_core). Units a = 1, m = theta.",
    "development_seen_before_prereg": [
        "Declared honestly. Engine checks: Riesz formula vs finite differences (5e-10); k-space code vs real-space "
        "exact diagonalisation of zz_core walks: 1-D ring N=64 (1e-8 rel) and 3-D 8^3 box (1e-8 rel) with the "
        "species projector; whole one-tick 1-D sea uniform to 2e-15 under V(0) = -0.8.",
        "Species structure seen: 1-D one tick: doubler at k = pi (mass theta, gap at quasi-energy pi). 3-D two "
        "ticks: per pi-cell 4 massive species at (0,0,0), (pi/2,pi/2,0)+perms and 4 MASSLESS cones at quasi-energy "
        "pi at (pi/2,0,0)+perms, (pi/2,pi/2,pi/2) (velocity cos theta).",
        "Development numbers already seen (so thresholds below are NOT blind for these): 1-D physical species vs "
        "continuum err +0.02% (q=0.3..1 m) to -0.46% (q = 10 m) at theta = 0.025, scaling ~theta^2. 3-D: pit(0) "
        "slope between successive theta (0.2..0.025) -11.4%, -4.2%, -1.4% from 1/(6 pi^2). Pi_R along x vs Uehling "
        "at theta 0.1/0.05/0.025: -4.4/-1.8/-0.6% (q = 0.3 m); along (1,1,1) at q = 3 m: -35..-40% (because pit(0) "
        "is ANISOTROPIC: x,y,z = 0.0281, 0.0309, 0.0344 at theta = 0.1, the two time slices swap x and z). With "
        "per-direction subtraction and average over x,y,z: +0.07% (theta 0.1), -0.28% (0.05) at q = 1 m. Whole 3-D "
        "sea: IR slope d pit/d ln q = 0.01692 vs 1/(6 pi^2 cos theta) = 0.01697 at theta = 0.1.",
    ],
    "declared_choices": {
        "species_projector": "physical species = smooth projector Q = sum_k sqrt(w(|k|)) |k><k| around k = 0, "
                             "w = C-infinity step from 1 at |k| = K1 = 0.6 to 0 at K2 = 1.05 (lattice units), the SAME "
                             "for every theta, every q, 1-D and 3-D. This is the declared correction for the doublers "
                             "(like R1 'count only the k = 0 gap'); it is also the regulator. Robustness (reported, "
                             "graded in V2g): K1, K2 = 0.5, 0.9.",
        "renormalization": "on-shell (Thomson) scheme, the same for every case: Pi_R(q n) = Pi(0 n) - Pi(q n) with "
                           "Pi(0 n) = q -> 0 limit along the same direction n (Richardson from q = 0.05 m and 0.025 m). "
                           "Observable = average over n = x, y, z (= average over the two time slices; for s states "
                           "only the angular average enters at first order).",
        "free_parameters": "0 fitted. Inputs: theta (= m a), alpha, Z. K1, K2 are the regulator (declared, varied).",
    },
    "configs": {"theta_1d": [0.2, 0.1, 0.05, 0.025], "theta_3d": [0.1, 0.05, 0.025, 0.0125],
                "q_over_m_exam": [0.3, 1, 3, 10], "q_table_3d": "log grid q_lat from 0.02 theta to 0.6, ~40 points, "
                "directions x, y, z (for r-space and levels)",
                "quadrature_3d": "product sinh-Gauss n = 64 (beta 5, 6 at theta 0.0125); convergence check n = 80",
                "atom": "radial walk of zz_core (RadialWalk), H (Z=1): 1s, 2s, 2p1/2 at theta_atom 0.05, 0.025; "
                        "Nd59+ (Z=60) 1s at theta_atom 0.05, 0.025, 0.0125; electron density form factor from the "
                        "rule's (g, f)"},
    "exams": {
        "V1a_whole_sea_1d_theorem": "1-D one-tick whole sea: S U S = -U predicts NO polarization at all orders. "
                                    "PASS if real-space sea density spread < 1e-12 (N = 128 ring, V(0) = -0.8 per tick) "
                                    "and |chi_whole(q)| < 1e-12 for all q tested.",
        "V1b_1d_physical_species": "chi_lat/chi_cont - 1 (continuum 1+1 D Dirac, per species) for q/m in {0.3,1,3,10}: "
                                   "PASS if <= 1% at theta = 0.025 for all q AND the error falls by >= 3x from "
                                   "theta = 0.05 to 0.025 for q <= 3 m (theta^2 scaling).",
        "V1c_period_independence": "physical-species chi with one-tick vs two-tick Floquet period: PASS if equal to 1e-9.",
        "V2a_whole_sea_3d": "whole two-tick 3-D sea (all species, per massive species): predicted massive minus "
                            "massless-cone response -> IR slope d pit/d ln q = +1/(6 pi^2 cos theta) (anti-screening "
                            "log at LARGE distance). PASS if slope between q = 0.05 m and 0.1 m within 2% at theta = 0.1 "
                            "AND pit_whole(q) < 0 for q/m in {0.05..3}.",
        "V2b_log_coefficient": "pit(0)(theta) (xyz average) fitted as L (-ln theta) + c + d theta^2 on theta in "
                               "{0.1, 0.05, 0.025, 0.0125}: PASS if L within 2% of 1/(6 pi^2) (= 2 alpha/(3 pi) after "
                               "4 pi alpha) AND the finite difference between 0.025 and 0.0125 within 2%.",
        "V2c_uehling_shape": "Pi_R,lat (xyz avg) / Pi_R,Uehling - 1 at q/m in {0.3, 1, 3, 10}: PASS if at theta = 0.0125 "
                             "|err| <= 0.5% for q <= 3 m and <= 1.5% at 10 m.",
        "V2d_q2_coefficient": "pitR(q)/pitR_Uehling(q) at q = 0.2 m (the coefficient that sets the H level shift, "
                              "alpha/(15 pi) q^2/m^2): PASS if within 1% at theta = 0.0125.",
        "V2e_isotropy": "Pi_R along (1,1,1) (per-direction subtraction) vs xyz average at q = 3 m, theta = 0.0125: "
                        "PASS if within 1%.",
        "V2f_anisotropy_of_charge_renormalization": "prediction (from development): pit(0) differs along x, y, z by a "
                        "constant that does NOT vanish as theta -> 0. Reported: (pit_z - pit_x)/(1/6pi^2) and the "
                        "slice-averaged (pit_xz - pit_y)/(1/6pi^2) at each theta. PASS if (pit_z - pit_x) changes by "
                        "< 30% between theta 0.05 and 0.0125 (i.e. it is a cutoff constant, not a theta^2 error).",
        "V2g_regulator_independence": "Pi_R (xyz avg) at q = 1 m and 3 m, theta = 0.025, with (K1,K2) = (0.5,0.9) vs "
                                      "(0.6,1.05): PASS if equal within 0.5%; and L (finite difference 0.05->0.025) equal "
                                      "within 2%.",
        "V3_uehling_potential": "r-space potential V_U,lat(r) from the q-table at theta = 0.0125, nucleus smeared by "
                                "the SAME Gaussian factor exp(-(q/q_c)^2), q_c = 0.3/theta (lattice 0.3), applied to the "
                                "analytic reference: PASS if |V_lat/V_ref - 1| <= 2% at r/lambda_C in {0.05, 0.1, 0.2, "
                                "0.5, 1, 2} and <= 3% at r = 3 (tail ~ e^{-2r}).",
        "V4a_H_2S": "Uehling shift of H 2S1/2 = <V_U,lat> in the rule's 2s state (theta_vac = 0.0125, theta_atom -> 0 by "
                    "Richardson in theta^2 from 0.05, 0.025), infinite nuclear mass. Reference: analytic Uehling "
                    "potential in the Schroedinger 2s density (relativistic correction O((Z alpha)^2) ~ 1e-4, declared); "
                    "leading formula -4 alpha (Z alpha)^4 m c^2/(15 pi n^3) = -27.15 MHz also quoted. PASS if within 1%.",
        "V4b_H_1S": "same for 1S1/2. PASS if within 1%.",
        "V4c_H_2P": "2P1/2 shift: PASS if |dE| < 0.01 MHz (Uehling on p states is (Z alpha)^6-suppressed).",
        "V4d_Nd_1s": "Nd59+ (Z = 60, point nucleus) 1s Uehling shift, theta_vac = 0.0125, theta_atom extrapolated "
                     "(a^(2 gamma) law, gamma = 0.899, from 0.05, 0.025, 0.0125). Reference: analytic Uehling in the "
                     "exact Dirac 1s density (closed form). PASS if within 2%. Only first order in Z alpha (no "
                     "Wichmann-Kroll, declared).",
        "V4e_nonperturbative_consistency": "Nd59+ 1s: E(V + V_U,lat) - E(V) from the radial walk eigenvalue vs <V_U,lat>: "
                                           "PASS if within 1% (theta_atom = 0.025).",
    },
    "what_is_NOT_expected": "self-energy (the +1085 MHz of the 2S Lamb shift), g-2, the 25 MHz of the fine structure "
                            "need quantized light: not addressed; the 1S-2S lattice error (debt #5) is not fixed here.",
}
json.dump(P, open(os.path.join(OUT, "vuoto_prereg.json"), "w"), indent=1, ensure_ascii=False)
print("written", os.path.join(OUT, "vuoto_prereg.json"))

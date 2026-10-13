"""emiss_prereg2.py - PRE-REGISTRATION v4 ("Emissione vera", phase 2).  Written and committed BEFORE the phase-2
production runs.  v3 (lab/out/emiss_prereg.json) is not touched; its thresholds are reused unchanged.
Run: python lab/emiss_prereg2.py -> lab/out/emiss_prereg2.json"""
import json, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = {
    "study": "Emissione vera (prefix emiss_), pre-registration v4 = phase 2",
    "date": "2026-10-13",
    "why": "v3 failed E1v3c, E1v3s, R1_48, S1, K10, SC10 because the 48^3 1s/2p states were pieces of valley multiplets "
           "(docs/BATTIMENTO_1S.md, master 83485bc). Phase 2 redoes them with exact levels of the rule and the learned nucleus "
           "of docs/NUCLEO_1S.md sec. 13, on the FULL cube (sum of the 4 colours; docs/SOTTORETICOLO.md master 178743d: one "
           "colour alone is anisotropic for dipoles), adds the retardation term and judges L1 on the regenerated export.",
    "development_runs_declared": [
        "emiss_v4_dev.py (2400 ticks, Hann, even ticks): spectra of the radial starts. Bare nucleus: 1s+ lines 0.89944 (0.49), "
        "0.8863 (0.24), 0.84038 (0.20); 2p1/2+ one line 0.97558; 2p3/2 one line 0.9788 (the rest is window side lobes). "
        "Learned nucleus: 1s+ ONE line 0.89946 (0.98), 2p1/2+ 0.97596, 2p3/2 0.97882 (single lines at this resolution).",
        "reading of v3 (post-hoc, emiss_posthoc.py): the 2p starts are already single lines with the bare nucleus; the dipole "
        "loss (x0.23-0.36) is mostly from the 1s exact level being a 25/25/25/25 valley mix (BATTIMENTO_1S) - expected to be "
        "cured by the learned nucleus (valleys closed, the smooth valley-0 1s becomes stationary).",
    ],
    "configurations": {
        "atom": "48^3, a = 0.5, Z = 60, full cube (4 colours), nucleus on the vertex of 8 cells; LEARNED NUCLEUS: Coulomb "
                "phase x 0.554593 on the 8 sites at r = sqrt3/2 a and x 1.003549 on the 24 sites at r = sqrt11/2 a "
                "(NUCLEO_1S sec. 13.2, learned with an internal criterion, 2 numbers, not refitted here)",
        "states": "filter diagonalization (Rayleigh-Ritz, as Atomo::uno_s_esatto of BATTIMENTO_1S): from the radial start of "
                  "each state (kappa, m), 1200 ticks, even ticks, Hann window, filters every 0.005 in [0.82, 0.94] (1s) or "
                  "[0.94, 1.01] (2p); Ritz vectors of the two-tick rule in that subspace. BLIND criterion fixed now: the Ritz "
                  "vector with the largest weight in the radial start; its cluster = Ritz vectors within 0.0015 of it (half "
                  "the 2p fine-structure split; valley copies closed by the learned nucleus are within 0.0002-0.0005); "
                  "state = normalised projection of the radial start on the cluster. No Dirac number is used.",
        "states_list": "1s+ (m = +1/2), 2p1/2+, 2p1/2-, 2p3/2+ (m = +1/2), 2p3/2+3 (m = +3/2), 2p3/2- (m = -1/2)",
        "coupling": "coupling_v3 of the v3 pre-registration (energy-consistent field, transverse deposit), light 128^3, "
                    "b = 4 (20 cells/lambda), sponge 12, sigma 1 cell, K = 1",
        "runs": "P4_E1pi, P4_E1sig (semiclassical, pulse as E1v3 at w(2p1/2+) - w(1s+) of the new states, T = 10000); "
                "P4_R1 (qe, pure 2p3/2+3, T = 3000); P4_S1sc (semiclassical control from the pure 2p3/2+3, T = 3000)",
    },
    "tests": {
        "E4a": "pi pulse: strongest line at w(1s->2p3/2) of the new states within 0.5 % (= E1v3a)",
        "E4b": "second line at w(1s->2p1/2) within 0.5 % (= E1v3b)",
        "E4c": "amplitude ratio (3/2 : 1/2) vs prediction from the new lattice dipoles, pulse spectrum and w^2 within 20 % (= E1v3c)",
        "E4s": "sigma- pulse, probe E_y: strongest line at w(1s->2p1/2) within 0.5 % (= E1v3s)",
        "R4": "qe, pure 2p3/2+3 -> 1s+, K = 1: A_norm / (A_rule W2) within 3 %, A_rule from the dipole at t = 0 (= R1_48)",
        "R4b": "A_photons / A_norm within 3 % (= R1b_48)",
        "S4": "1 - norm(T) >= 0.5 A_rule W2 T (qe) AND semiclassical control 1 - p_2p(T) < 1 % of the qe depletion (= S1)",
        "D4": "informative: sum |<1s|r|2p>|^2 of the new states vs Dirac continuum (box 24 lambda_C confines the 2p: < 1 expected); "
              "|x| vs |y| for 2p3/2+3",
        "V4": "informative: the new 1s+ left alone (no light) 3000 time units: max (1 - |<1s|psi(t)>|^2)",
        "N2": "retardation: the rule's own radial states (radial walk, a = 0.025, as P_N1) put into the full golden rule with "
              "e^{-ik.r} (the code of emiss_exact.py with the rule's g, f, w instead of Dirac's): F_ret(rule) = A_EX/A_LW. "
              "Test: (A_N1 / W2) x F_ret(rule) vs Dirac A_EX(2p3/2) = 9.9787e-6 within 3 %. Check inside: the same code with "
              "k -> 0 (velocity form) vs the length form on the rule's states within 2 %.",
        "N2x": "informative: F_ret(rule) vs Dirac 0.91662",
        "L1b": "regenerated Laboratory export (exact 1s, /mnt/project-files/ripresa/dati/atomo_accoppiato/ dated 01:46): the "
               "amplitude of p_exc at the old beat (largest Hann line in [0.010, 0.016], t > 400, amplitude = 2|sum h x e^{iwt}|/sum h) "
               "< 0.01. Reference (informative): the same estimator on the free run of the old Laboratory 1s (emiss_lab1s.json)",
    },
    "not_changed_after": "thresholds and configurations above; changes go to lab/out/emiss_amend2.json, declared",
}
if __name__ == "__main__":
    json.dump(P, open(os.path.join(OUT, "emiss_prereg2.json"), "w"), indent=1)
    print("written", len(P["tests"]), "tests")

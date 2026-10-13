"""emiss_prereg3.py - PRE-REGISTRATION v5 ("Emissione vera", phase 3).  Written and committed BEFORE the phase-3 runs.
v3/v4 are not touched; their thresholds are reused unchanged where a test is repeated.
Run: python lab/emiss_prereg3.py -> lab/out/emiss_prereg3.json"""
import json, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = {
    "study": "Emissione vera (prefix emiss_), pre-registration v5 = phase 3",
    "date": "2026-10-13",
    "inputs_seen_before": [
        "phase 2 (v4): S4 failed on the projection (semiclassical 1 - p_2p = 0.0356 = free drift of the 2p) while the energy the "
        "semiclassical pure 2p gave to the light was 4e-33 vs ~3e-4 for the quantum one: S1/S4 are therefore re-asked here on the "
        "emitted energy, but ONLY on NEW runs (96^3), never re-judged on the phase-1/2 data",
        "phase 2: dipoles 1.37-1.55 x Dirac on 48^3 (box 24 lambda_C; 24 % of the 2p beyond r = 10)",
        "phase 2: E4 spectra have weak lines at +0.0131 from the main lines; 2p3/2 near-degenerate doublet 0.97870/0.97887 "
        "(2400-tick Ritz, development)",
        "docs/SURFISTI_ESATTI.md 12.8: the rule has an exact time reversal with T^2 = -1 (every level is a Kramers pair m <-> -m)",
    ],
    "configurations": {
        "atom96": "96^3, a = 0.5 (box 48 lambda_C), Z = 60, full cube (4 colours), learned nucleus of NUCLEO_1S 13.2 "
                  "(0.554593 / 1.003549, the same 2 numbers: they are per site, independent of the box)",
        "states96": "filter diagonalization as v4 but windows [0.88, 0.92] (1s) and [0.955, 0.995] (2p), filters every 0.005, "
                    "1200 ticks, even ticks, Hann; BLIND criterion of v4 (largest weight in the radial start + Ritz vectors within 0.0015)",
        "runs": "P5_R1_96: qe, pure 2p3/2+3 -> 1s+, K = 1, coupling_v3, light 128^3 b = 4, T = 3000. "
                "P5_S1sc_96: semiclassical control from the same pure 2p, T = 3000",
        "diagnostics48": "48^3 learned nucleus: (a) Ritz with 2400 ticks, filters every 0.0025 in [0.955, 0.995], from the 2p3/2+3 "
                         "and 2p3/2+ radial starts: the two members of the doublet; (b) the final state of P4_E1pi evolved freely "
                         "2400 ticks: spectrum of <psi_end|psi(t)>",
    },
    "tests": {
        "S5": "energy given to the light by the pure 2p (96^3): semiclassical E_sc(T) < 1 % of the quantum E_qe(T)",
        "S5b": "the quantum pure 2p decays by itself: 1 - norm(T) >= 0.5 A_rule W2 T (as S1)",
        "R5": "qe 96^3: A_norm / (A_rule W2) within 3 % (= R1_96 of v3), A_rule from the dipole at t = 0",
        "R5b": "A_photons / A_norm within 3 %",
        "D5": "dipole of the 96^3 states: sum |<1s+|r|2p3/2+3>|^2 vs Dirac continuum R^2/3 = 2.4806 within 5 %",
        "R2v5": "A_norm / W2 vs Dirac A_LW(2p3/2) = 1.08865e-5 within 5 % (= R2_96 of v3; EXPECTED FAIL if w(96^3) stays ~3.7 % "
                "below Dirac: -10.6 % in w^3)",
        "W5": "doublet: valley weights (Fourier of each member, k reduced mod pi per component: valley 0 = all components "
              "within pi/4 of 0 mod pi; valley xy/xz/yz = two components near pi/2, one near 0). Decision rule fixed now: "
              "VALLEYS if each member has >= 0.15 in at least 3 valleys; NOT VALLEYS (cubic/anisotropy of the walk) if each "
              "member has >= 0.80 in valley 0; otherwise UNDECIDED. Also reported: the weight of the 2p3/2+ (m = 1/2) start on the "
              "two members (a cubic/tetragonal split separates |m| = 3/2 from |m| = 1/2)",
        "X5": "weak lines at +0.0131: in the free spectrum of the P4_E1pi final state, H1 = a level near w(1s) - 0.0131 "
              "(1s-type, [0.880, 0.892]) carries more weight than any level near w(2p) + 0.0131 ([0.985, 0.996]); H2 = the "
              "opposite. Reported with the weights; decided by the larger weight (no threshold)",
        "M5": "dynamic multipole coupling: only if time remains; otherwise declared not done",
    },
    "not_changed_after": "thresholds and configurations above; changes go to lab/out/emiss_amend3.json, declared",
}
if __name__ == "__main__":
    json.dump(P, open(os.path.join(OUT, "emiss_prereg3.json"), "w"), indent=1)
    print("written", len(P["tests"]), "tests")

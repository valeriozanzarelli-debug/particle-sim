"""emiss_prereg.py - PRE-REGISTRATION v3 of the light-atom coupling ("Emissione vera", docs/EMISSIONE.md).
Written and committed BEFORE any production run.  The v2 pre-registration (docs/ATOMO.md 3.5, test accoppiato_v2)
is not touched.  Run: python lab/emiss_prereg.py -> lab/out/emiss_prereg.json
"""
import json, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

P = {
    "study": "Emissione vera (prefix emiss_), pre-registration v3",
    "date": "2026-10-13",
    "rule": "zig-zag rule for the electron (3-D Dirac walk, a = 0.5, or its radial form) and for the light (3-D spin-1 walk); "
            "Nd59+ (Z = 60) unless said; charge e = sqrt(4 pi alpha) TRUE (no K boost) unless a test says K",
    "inputs": "alpha, Z, m_e (flip), c, hbar. Exact Dirac/NIST numbers only as exams (lab/out/emiss_exact.json).",
    "development_runs_declared": [
        "levels of the 3-D rule by listening (48^3: 1s 0.89946-0.89957, 2p1/2 0.97555-0.97580, 2p3/2 0.97880; 96^3: 1s 0.89955, 2p1/2 0.97125, 2p3/2 0.97340 from short spectra)",
        "exact Dirac rates (emiss_exact.py): Z=60 A_LW(2p3/2) = 1.08865e-5, A_EX = 9.979e-6 (retardation -8.3 %), A_LW(2p1/2) = 1.0403e-5, A_EX = 1.0705e-5; Z=1 A_EX = 8.0742e-13 = 6.2683e8 s^-1 (x mu/m -> 6.2649e8)",
        "light-only dipole (20 cells/lambda): P / (Larmor W2) = 1.017 (dev_light); Bloch prediction of the light rule 1.006 (20), 1.003 (30), 1.002 (40)",
        "qe 48^3, K=1, coupling as accoppiato.rs (field read before the light tick, local deposit): norm decays 2.2x faster than the photons (dev_qe48A): "
        "a time-lag artifact of the electron's own near field; fixed by the energy-consistent field (1/2 (U F_n + F_n+1)) and by the transverse "
        "(Coulomb-gauge) deposit: dev_qe48B norm/(A_rule W2) 0.994, dev_qe48C 1.057 (photons 1.032) with T = 600 and the S_T = 300 states",
        "radial Nd (a = 0.05) qe: norm/(A_rule W2) = 1.013, 1.006, 0.996 at 20, 30, 40 cells/lambda; photons agree with norm within 0.1 %",
        "toy 1-D (emiss_twa.py, 60 trajectories): hybrid truncated Wigner starts at ~0.58 A and pumps the ground state to ~0.28-0.30; spin TWA ~ A",
        "Laboratory 1s replay (emiss_lab1s.py, running when this file was written): Laboratory 1s = candidate 0.9000 of a 600-tick listening",
    ],
    "configurations": {
        "light": "128^3 cells of b = lambda/20 (3-D tests) or lambda/30 (radial tests), sponge 12 cells (cubic, e^-7.5), Gaussian average sigma = 1 cell",
        "coupling_v3": "energy-consistent field (atom sees 1/2(U F_n + F_n+1), predictor + corrector for the self term) and transverse deposit "
                       "(P_perp of the lattice light walk: axis of its Bloch rotation), unless said 'as v2'",
        "W2": "power factor of the Gaussian average, computed (Light3D.form_factor2), never fitted; rates are reported raw and divided by W2",
        "states": "levels of the 3-D rule listened with a Gaussian window S_T = 1500 (48^3) or 1000 (96^3); radial states by inverse iteration of the radial walk",
        "rate_fit": "A_norm = -slope of ln(norm of the no-photon sector) over t >= 3 periods to the end; A_photons = slope of 2 (energy in box + absorbed)/w over the same window",
        "A_rule": "(4/3) alpha w^3 |<1s|r|2p>|^2 with w and the dipole of the SAME lattice states (listened)",
    },
    "tests": {
        "E1v3a": "lab config 48^3, semiclassical as v2 (K = 1), prescribed pi pulse eE0 = 0.004 at w(1s->2p1/2) of the listened states, tau = one period, "
                 "T = 10000; Hann spectrum of the probe E_z (40 cells behind) after the pulse: strongest line at w(1s->2p3/2) within 0.5 %",
        "E1v3b": "same run: second line at w(1s->2p1/2) within 0.5 %",
        "E1v3c": "same run: amplitude ratio (2p3/2 line / 2p1/2 line) vs prediction from the lattice dipoles, pulse spectrum and w^2 within 20 %",
        "E1v3s": "same with sigma- circular pulse (Delta m = -1), probe E_y: strongest line at w(1s->2p1/2) within 0.5 % (selection rules: 2p1/2 2:1)",
        "E1v3zc": "informative: the v2 estimator (zero crossings, 400-1200 after the pulse) on E1v3a",
        "R0": "light only, prescribed dipole current, 20 cells/lambda: P / (Larmor W2) within 3 % of 1",
        "R1_48": "qe 48^3, pure 2p3/2 m=3/2 -> 1s m=1/2, K = 1, coupling_v3, T = 3000: A_norm / (A_rule W2) within 3 %",
        "R1b_48": "same run: A_photons / A_norm within 3 %",
        "S1": "same run: the pure 2p decays by itself: 1 - norm(T) >= 0.5 A_rule W2 T; control semiclassical from the same pure 2p (coupling_v3, K=1, T=3000): "
              "1 - p_2p(T) < 1 % of the qe depletion",
        "R1_96": "qe 96^3, as R1_48 (T = 3000): A_norm / (A_rule W2) within 3 %",
        "R2_96": "same run: A_norm / W2 vs Dirac A_LW(2p3/2) = 1.08865e-5 within 5 % (EXPECTED FAIL: the 96^3 rule's w(1s->2p3/2) = 0.0739 is 3.7 % "
                 "below Dirac, -10.6 % in w^3, known before this file)",
        "N1": "radial Nd59+, a = 0.025, 30 cells/lambda, qe, K = 1, 25 periods: A_norm/(A_rule W2) within 3 %; A_norm / W2 vs Dirac A_LW(2p3/2) within 3 %",
        "N1x": "informative: A_norm / W2 vs Dirac A_EX (with retardation): expected about -8 % (the coupling is dipole)",
        "H1": "radial hydrogen Z = 1, a = 0.5 and 0.25, 30 cells/lambda, qe, K = 1, 12 periods: each A_norm/(A_rule W2) within 3 %; "
              "Richardson in a^2 of A_norm/W2: vs 6.2683e8 s^-1 (Dirac, infinite nuclear mass) within 3 %",
        "K10": "check, not production: qe 48^3 K = 10 pure 2p, T = 3000: max |norm(t) - exp(-A t)| < 0.02 with A fitted; A(K=10)/100 vs A(K=1) within 5 %",
        "SC10": "semiclassical 48^3 K = 10 from (1s + 2p3/2)/sqrt2, coupling_v3, T = 3000: p(t) follows the logistic of dp/dt = -A p (1-p) with A of K10 "
                "(no fit): max deviation < 0.05",
        "T1": "toy 1-D (2000 trajectories): hybrid truncated Wigner initial rate (p over t = 50..300) / A_QE = 0.5 +- 0.1",
        "T2": "toy: hybrid truncated Wigner from the ground state: p(t = 3000) = 0.293 +- 0.05",
        "T3": "toy: spin truncated Wigner from e: -ln p slope over 50..1500 within 15 % of A_QE",
        "T4": "toy: semiclassical from pure e: p stays 1 within 1e-6",
        "L1": "Laboratory replay: the free beat period of p_exc predicted from the decomposition of the Laboratory 1s on the clean levels "
              "(2 pi / (w_a - w_b), no parameter) vs the beat period of the Laboratory series (serie.csv, t > 400) within 5 %",
    },
    "not_changed_after": "thresholds and configurations above; any change goes to lab/out/emiss_amend.json, declared",
}

if __name__ == "__main__":
    json.dump(P, open(os.path.join(OUT, "emiss_prereg.json"), "w"), indent=1)
    print("written", len(P["tests"]), "tests")

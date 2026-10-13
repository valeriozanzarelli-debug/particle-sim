"""vuoto_f2_prereg_L.py - PRE-REGISTRATION of the Bethe-log exams (phase 2, part L), written and committed BEFORE
lab/vuoto_f2_bethe.py is run on hydrogen. Writes lab/out/vuoto_f2_prereg_L.json."""
import os, json

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = {
    "study": "VUOTO_POLARIZZA phase 2, part L (Bethe logarithm from the resolvent of the radial rule)",
    "date": "2026-10-13",
    "development_seen (declared; hydrogen Z = 1 NEVER computed before this file)": [
        "Z = 20, 40, 5 dense walks and Z = 2, 3 resolvent runs (theta 0.1, 0.05): (1) the unfiltered dipole source r g picks "
        "a tiny content at the band top of the 1-D walk (k a ~ pi, the doubler side) whose dE^3 weight dominates every "
        "sum (D up to 300x the sum rule); FIX: species low-pass of the source + positive-energy projector + smooth band "
        "filter keeping quasi-energies below Omega_n + 0.6 (the same 0.6 = K1 of the species projector). (2) The textbook "
        "ln k0 is a K -> infinity limit of the NON-relativistic sum; the walk is Dirac, and above K ~ 10^3 Ry its spectrum "
        "is (correctly) relativistic, so a K -> infinity extrapolation of the walk is ambiguous (fits from 2.5 to 3.2 "
        "depending on the window). (3) With the fix, Z = 2: walk ln k0(K) = 1.522/1.529 (kappa' +1/-2) at K = 99 Ry and "
        "1.916/1.931 at 309 Ry, vs non-relativistic Schroedinger 1.536 (100 Ry) and ~1.94 (309 Ry).",
        "Exam reference computed before this file: lab/vuoto_f2_nrref.py -> lab/out/vuoto_f2_nrref.json (Schroedinger "
        "finite differences, same cutoff definition, two grids, Richardson h -> 0)."],
    "definition": "ln k0(K) = ln(K/Ry) - int_0^K S(k) dk / (K S(K)), S(k) = sum_m |<m|v|n>|^2 dE/(dE + k) "
                  "(velocity form from the dipole length: v_mn = i dE r_mn); intermediate channels kappa' = +1 (1/3) and -2 "
                  "(2/3); boxes 40 a_B (1S), 60 a_B (2S); theta_atom = 0.05 (graded) and 0.1 (stability)",
    "exams": {
        "L1": {"what": "rule vs non-relativistic QM at the same photon cutoff: |ln k0_rule(K) - ln k0_NR(K)| for 1S and 2S, "
                       "K = 100 and 300 Ry, theta_atom = 0.05", "pass": "<= 0.02 (all four)"},
        "L2": {"what": "lattice stability: |ln k0_rule(300 Ry; theta 0.05) - same at theta 0.1|, 1S and 2S", "pass": "<= 0.01"},
        "L3": {"what": "dipole strength below K: D_K/(2 pi Z alpha |psi(0)|^2) rule vs NR at K = 300 Ry, 1S and 2S",
               "pass": "ratio within 1 %"},
        "L-report (no threshold)": "K = 30, 1000, 3000 Ry (relativistic onset); 'hybrid' ln k0 = rule below K + NR theory "
                                   "above K (= ln k0_rule(K) + 2.984 - ln k0_NR(K)): it is NOT an emergent number, it shows "
                                   "what the rule contributes; MHz per unit of ln k0: 1S 1087, 2S 136 (leading formula)"},
    "honesty": "This tests the low-energy (Bethe) part only, against NR QM with the same cutoff. The textbook number "
               "needs the K -> infinity non-relativistic tail, which the relativistic rule does not have: the matching with the "
               "high-energy part (the k-space loop of part S/G) is NOT done here.",
}
json.dump(P, open(os.path.join(OUT, "vuoto_f2_prereg_L.json"), "w"), indent=1, ensure_ascii=False)
print("written lab/out/vuoto_f2_prereg_L.json")

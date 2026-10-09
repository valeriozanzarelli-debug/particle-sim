"""orb_prereg.py - pre-registration of the ORBIFOLD pre-test (O1-O4). Written and committed BEFORE any orb_* run.
Writes lab/out/orb_prereg.json.

Fact check done on paper before writing (to be confirmed by the runs, corrected in docs/ORBIFOLD.md):
 * T^4/Z2 (x -> -x): 16 fixed points; the singular quotient has b = (1,0,6,0,1), chi = 8; blowing up each of the
   16 A1 points (P^1, chi 1 -> 2) gives K3: b2 = 22, chi = 24. The orbifold ('stringy') Euler number
   (1/|G|) sum_{gh=hg} chi(X^{g,h}) = (0 + 3*16)/2 = 24 already gives the resolved value.
 * T^2/Z3 (hexagonal, 120 deg): 3 fixed points, the quotient is a SPHERE with three cone points (chi = 2).
   It is NOT Calabi-Yau (a CY 1-fold is a T^2): no invariant spinor is expected. It is only a building block.
 * T^6/Z3 (z_i -> w z_i on three hexagonal T^2): 27 fixed points; untwisted h11 = 9, twisted 27 -> h11 = 36,
   h21 = 0, chi = 72; standard embedding: 36 net generations. Singular quotient: b = (1,0,9,2,9,0,1), chi = 18;
   each C^3/Z3 point resolved by O(-3)->P^2 adds 2 -> 18 + 54 = 72.
 * Wilson lines: Ibanez, Nilles, Quevedo PLB 187 (1987) 25 and Ibanez, Kim, Nilles, Quevedo PLB 191 (1987) 282:
   discrete Z3 Wilson lines split the 27 fixed points into classes; 3-family models exist. On T^2/Z3 the Wilson
   lines on e1 and e2 = w e1 must be equal and of order 3.
 * CY <=> Ricci-flat Kahler with holonomy in SU(n) <=> (compact, simply connected case) a parallel spinor.
"""
import json, os

P = {
    "written": "2026-10-09, before any orb_* script was run (only paper fact-check above)",
    "files": "lab/orb_core.py, lab/orb_o1.py, lab/orb_o2.py, lab/orb_o3.py, lab/orb_figure.py",
    "O1_holonomy": {
        "method": "Spinor walkers: (a) 4-D zig-zag split-step walk on an n^4 torus (n = 4, 6), spin matrices = 4 "
                  "anticommuting 4x4 alpha_j, flip = exp(-i theta beta); (b) 2-D/6-D hexagonal Dirac walkers as "
                  "Hamiltonian generator H = sum_k (n_k.sigma) sin + Wilson term (the split-step tick is not exactly "
                  "C3-symmetric: declared). The spinor matrix S of the orbifold operator O = S x P_g is SOLVED from "
                  "the walk (all 4x4 / 2x2 matrices with W O = O W, linear solve). Phase of S fixed by requiring S to "
                  "be the ball's own rotation exp(1/2 sum phi_ab alpha_a alpha_b) (both spin lifts tried). "
                  "Holonomy of a loop around fixed point f = S(g_f); invariant spinors = common +1 eigenspace.",
        "exams": {
            "K3 = T4/Z2 (angles pi,pi)": "walk commutes (residual < 1e-10) and invariant spinors >= 1 (expected 2, "
                                         "same chirality)",
            "T6/Z3 (2pi/3 x3)": "invariant spinors >= 1 (expected 2)",
            "T2/Z3": "expected 0 (not CY) - graded as a prediction",
        },
        "controls": {
            "T4, Z2 on one complex coordinate (angles pi,0)": "0 invariant spinors",
            "T6, Z3 on one complex coordinate (2pi/3,0,0)": "0 invariant spinors",
            "T4, reflection x1 -> -x1 with mass theta > 0": "no S at all (solution space of W O = O W is {0})",
        },
        "pass": "all 6 items as stated",
        "also_reported": "how many invariant spinors appear if an extra U(1) phase on S is allowed (a hidden-circle "
                         "twist): expected to rescue the controls -> this is the loophole, reported not graded",
    },
    "O2_counting": {
        "method": "cell complexes (triangulated hexagonal T^2, L = 6, 9; cubical T^4, n = 4, 6), invariant chain "
                  "complex / Hodge Laplacian zero modes on the quotient -> rational Betti numbers; fixed points "
                  "counted by brute force; T^6/Z3 via Kunneth on the measured Z3 traces of H*(T^2); Dirac zero modes "
                  "of the Wilson-Dirac operator restricted to the orbifold sector.",
        "exams": {
            "T2/Z3": "fixed 3, b = (1,0,1), chi = 2, Dirac zero modes 0",
            "T4/Z2": "fixed 16, b = (1,0,6,0,1), chi = 8, chi_orb = 24 = 8 + 16 (blow-up), b2_resolved = 22, Dirac "
                     "zero modes 2 of one chirality (= A-hat(K3) = 2)",
            "T6/Z3": "fixed 27, b = (1,0,9,2,9,0,1), chi = 18, chi_orb = 72, h11 = 9 + 27 = 36, h21 = 0, "
                     "net generations 36",
        },
        "pass": "every number exact at both lattice sizes",
        "declared_debt": "the blow-up (+1 per A1 point, +2 per C3/Z3 point) is NOT computed by balls; it is the "
                         "textbook local contribution added by hand",
    },
    "O3_wilson_lines": {
        "a": "T2/Z3 lattice: scan AB phases (phi1, phi2) on the two cycles over {0,1/3,2/3}^2 (x 2pi) and "
             "{0,1/2}^2: count which admit a Z3 operator commuting with the hopping. Expected: exactly the 3 Z3 "
             "phases with W(a2) = W(a1)^-1, W^3 = 1 (in the a1, a2 = rotated basis convention measured).",
        "b": "magnetised T2/Z3 (uniform AB flux M = 1..12 quanta, L = 36 triangular lattice): lowest Landau level "
             "has M states (check); split by the Z3 operator into (n0, n1, n2) for each allowed Wilson line. "
             "'Toy reaches 3' = some (M, Wilson line, sector) gives exactly 3. Declared in advance: this is a CHOICE "
             "of M and sector, not a prediction; parameter count reported.",
        "c": "T6/Z3 fixed-point bookkeeping: 27 points -> classes under 0, 1, 2 independent Z3 Wilson lines "
             "(expected 1x27, 3x9, 9x3). Which class carries massless families needs E8xE8 gauge-lattice data "
             "(shift vector, level matching): NOT done, said so.",
    },
    "O4": "optional; pre-registered as SKIPPED unless time remains (then reported as exploratory only)",
}

if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "orb_prereg.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(P, f, indent=1, ensure_ascii=False)
    print("wrote", out)

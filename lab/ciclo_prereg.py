"""ciclo_prereg.py - PRE-REGISTRATION of the study "Il ciclo tra i livelli" (docs/CICLO_LIVELLI.md), prefix ciclo_.
Each phase is ADDED here and committed BEFORE its production runs; entries of earlier phases are never edited
(changes after a commit go in lab/out/ciclo_amend_<n>.json, declared).

Phase 1 "misurare la mappa": re-analysis (read-only) of levels that already exist:
  1-D chains   lab/out/agg_zz_1d.npz  (ladder b = 2: B = 2..32; auto: B = 4, 16, 64, 192, 384; theta = a0/4096)
  3-D chain A  lab/out/agg_zz_3d_A.pkl (free learned generators, B = 16, 256, 512)
  3-D chain B  lab/out/agg_zz_3d_B.pkl (zoom chain, theta = 0.268: B = 16, 32)          [reported, no exam]
  proton chain lab/out/aggfull_proton.pkl (6 free levels, theta_p = 0.030)
  T levels     lab/out/aggfull_relearnT.pkl, aggfull_l1_T.pkl, aggmol_levels.pkl (c2, c4 (,c6) FIXED to the level
               below: NOT used for exponent exams, reported only)
  couplings    lab/out/aggfull_learn.json (B = 256, 512), lab/out/aggmol_learn.json (B = 1024): Darwin gk, field-motion gx,
               spin-orbit gso (FW language)
Declared: some of these numbers are already published (AGGLOMERATO_ZIGZAG A1b/A1c 1-D c4 ratios 0.997..1.028 and 1.497 at
B = 384; Darwin/SO coefficients per level in AGGLOMERATO_COMPLETO s.3 and AGGLOMERATO_DUE_ATOMI s.2). What is NEW here:
the exponents (slopes over levels), the 3-D free c4 and its isotropy, the proton chain c4, the normalized-shape
distances along the chains. Units: hbar = c = m_e = 1; h = cell in Compton lengths of the particle (h = B theta).
"""
import json, os, sys

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

PHASE1 = dict(
    title="Fase 1: misurare la mappa (quanto cambia la regola imparata da un livello all'altro) + registro delle influenze",
    date="2026-10-13",
    method=dict(
        scalar_generator="3-D: g_j = (G_00 + G_11)/2 of the learned generator (spin part 1e-10); 1-D: positive-band "
                         "eigenphase of the learned one-step kernel, G(kappa) = -(phase(kappa) - phase(0))",
        moments="3-D exact stencil moments: c2 = -1/2 sum g_j j_x^2; c4 = 1/24 sum g_j j_x^4; c_xxyy = 1/4 sum g_j j_x^2 j_y^2 "
                "(coefficient of kappa_x^2 kappa_y^2); averaged over the 3 axes / 3 planes. 1-D: least squares of G(kappa) "
                "on kappa^2, kappa^4, kappa^6 for 0 < kappa <= min(0.3 m h, 0.8) (the A1 fit window).",
        continuum="E = k^2/2m - k^4/8m^3 + ..., kappa = k h, step tau = h  =>  c2 = 1/(2 h), c4/c2 = -1/(4 h^2)",
        normalized_shape="S_n(kappa) = G_n(kappa)/c2_n (rest phase removed). Schroedinger fixed point: S* = kappa^2. "
                         "D_n = max_{0<kappa<=pi/2} |S_n - kappa^2| / (pi/2)^2 (1-D: axis; 3-D: axis and body diagonal). "
                         "Continuum prediction D_cont(h) = max |2 h^2 (sqrt(1 + kappa^2/h^2) - 1) - kappa^2| / (pi/2)^2.",
        free_levels="levels whose c2/c4 were NOT constrained: 1-D ladder, 1-D auto, 3-D A, proton chain",
        slope="least squares slope of log|c4/c2| against log h over ALL free levels of a chain (no level dropped)"),
    exams=dict(
        P1a=dict(what="kinetic term marginal (z = 2 fixed point): m_kin/m = 1/(2 h c2) at every free level of every chain "
                      "within 1e-3 of the fine rule's own tan(theta)/theta", tol=1e-3, prediction="PASS"),
        P1b=dict(what="relativistic k^4 term irrelevant with exponent y = -2: per chain (1-D ladder, 1-D auto, 3-D A, proton) "
                      "slope of log|c4/c2| vs log h within -2 +- 0.1", tol=0.1, prediction="PASS (1-D auto may be pulled by B = 384)"),
        P1c=dict(what="the irrelevant k^4 term comes out ISOTROPIC (rotations emerge on the cubic lattice): 3-D A free levels "
                      "c_xxyy / c_xxxx in [1.8, 2.2] (continuum Dirac: 2; a separable cubic artefact: 0)", lo=1.8, hi=2.2,
                 prediction="PASS at B = 16 and 256; B = 512 uncertain (k^4 is 1e-6 of k^2 there)"),
        P1d=dict(what="field-coupling corrections (Darwin gk, field-motion gx, spin-orbit gso) have the SAME coefficient at "
                      "B = 256, 512, 1024 within 2 % (relative spread) => their weight relative to the Coulomb term at the "
                      "cell scale scales as h^-2 (y = -2.00 +- 0.03 by a log fit)", tol=0.02, prediction="PASS (values known)"),
        P1e=dict(what="the map flows to a fixed form: along the 1-D ladder and the 1-D auto chain D_n decreases at every level, "
                      "and for h >= 2 lambda_C D_n is within 25 % of D_cont(h_n)", tol=0.25, prediction="PASS"),
    ),
    reported_not_judged=[
        "successive distance between normalized maps d_{n,n+1} = max |S_n - S_{n+1}| / (pi/2)^2 and its ratio per level",
        "the analytic fixed-point stencil of S* = kappa^2 on the coarse lattice: g*_0 = pi^2/3, g*_j = 2(-1)^j / j^2 per axis "
        "(separable, long range ~ 1/j^2), compared with the learned normalized stencils",
        "T-level, chain-B and aggmol B = 1024 numbers (constrained or too few levels)",
        "influence ledger (table + figure): weight of each term relative to the dominant one per scale, exponent y "
        "(measured where the engine has the term, analytic otherwise, marked)",
        "gravity: weight G m_p^2 / e^2 = 8.1e-37 per proton pair at the atom level; exponent +2 with the size of a "
        "neutral body (energy per particle ~ G M m / R ~ R^2 against a constant cohesive energy); crossover radius "
        "R_x = sqrt(3 eps v / (4 pi G (A m_p)^2)) computed for the engine's own solid H2 (eps, v from aggmol phase 3) and, "
        "for comparison only, with exam-only cohesive energies; the roundness radius needs the yield/cohesion fraction "
        "f_y, NOT derivable in the engine today (declared debt): R_round = sqrt(f_y) R_x",
    ],
)


def build(phases):
    P = dict(study="CICLO_LIVELLI: il ciclo tra i livelli (gruppo di rinormalizzazione automatico)", prefix="ciclo_",
             doc="docs/CICLO_LIVELLI.md", rules="CLAUDE.md, docs/PIANO.md s.2", phases=phases)
    json.dump(P, open(os.path.join(OUT, "ciclo_prereg.json"), "w"), indent=1)


if __name__ == "__main__":
    phases = dict(fase1=PHASE1)
    try:
        import ciclo_prereg_p2 as p2              # later phases append their own dicts (committed before their runs)
        phases.update(p2.PHASES)
    except ImportError:
        pass
    build(phases)
    print("wrote", os.path.join(OUT, "ciclo_prereg.json"), list(phases))

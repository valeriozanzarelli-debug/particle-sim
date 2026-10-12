"""colla_prereg.py - PRE-REGISTRATION of the 'glue' study (COLLA): which pieces are missing so that Skyrme knots bind
nuclei with the right energy. Written and committed BEFORE the production runs. Writes lab/out/colla_prereg.json.

Development seen before writing this (declared): (1) the massless unit fit reproduces ANW (e = 5.4445, F_pi = 129.26 MeV);
(2) lattice inertia of the massless B=1 (48^3, h=0.2) is Lambda = 78.2 vs 106.8 exact (-27 %: the 1/r^2 tail is cut by
the box) -> massless inertia tensors are NOT usable in this box, so the quantisation (C2) is done with massive pions;
(3) timings: one Hessian-vector product on 48^3 = 0.6 s.

Literature checked on the web before pre-registering (2026-10-12):
 * Adkins-Nappi 1984 (physical pion mass 138 MeV, fit to M_N, M_Delta): F_pi = 108 MeV, e = 4.84
   (confirmed via Leese-Manton-Schroers 1995, hep-ph/9502405: F_pi/4e = 5.58 MeV, 2/(e F_pi) = 0.755 fm, hbar = 46.8).
 * Braaten-Carson 1988 (rigid quantisation of the B=2 torus on the 8-dim manifold): deuteron binding ~ 158 MeV,
   charge radius 0.92 fm (as tabulated by LMS 1995).  The brief's "~6 MeV" is NOT Braaten-Carson: it is
   Leese-Manton-Schroers 1995 (quantum mechanics on a 10-dim attractive-channel manifold, massless pions,
   AN units): binding 6.18 MeV, radius 2.18 fm.  We do not reproduce LMS here (out of budget) - declared.
 * Battye-Sutcliffe 2006 (hep-th/0602220) and Battye-Manton-Sutcliffe 2007 (hep-th/0605284): with m = 1 shapes up to
   B ~ 8 unchanged (B=4 cube kept); B=4 cube E/(12 pi^2 B) = 1.307 at m = 1.  Binding not fixed by the pion mass.
 * Gillard-Harland-Speight 2015 (arXiv 1501.05455, NPB 895, 272) 'lightly bound' model:
   E = int (1-alpha)(|d phi|^2 + 2 m^2 (1-sigma)) + Skyrme + alpha (1-sigma)^4, m = 1, alpha = 0.95:
   E/8pi^2 = 1.1082 (B=1), 2.2119 (B=2), 3.3129 (B=3), 4.4111 (B=4); binding per nucleon 1.84, 3.26, 4.56 MeV with the
   calibration E_1 = 938.27 MeV; shapes = point-like lumps on FCC sites (no torus, no cube).
   Their sextic model (no Skyrme term) keeps 3.6-6 % binding: sextic alone does not fix it.
 * Loop (Casimir) corrections: B=1 about -200 MeV (Moussallam review); for the 3D B=2 torus never computed
   (Walliser-Holzwarth 1999). A full one-loop calculation needs counterterms (new parameters): not attempted.
"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
P = dict(
    study="COLLA - the missing glue of nuclear binding in the Skyrme-knot model",
    date="2026-10-12",
    common=dict(
        lattice="c4 stencil (knot2_core), h = 0.2, Dirichlet vacuum; L = 48 for alpha=0 runs, L = 52 for the GHS scan",
        binding="B_bind = B*M_1 - M_B, same lattice and same model for both terms",
        experiment_MeV={"d": 2.2246, "3He": 7.718, "4He": 28.296},
        parameter_count="every route lists: physical constants used, fit inputs, extra free parameters",
    ),
    routes=dict(
        C0="classical massless Skyrme (NUCLEI_DA_NODI): 2 fit inputs (M_N, M_Delta) -> d 72, 3He 178, 4He 310 MeV",
        C1="pion mass m_pi = 138 MeV (constant of nature, 0 extra free parameters) + refit of (e, F_pi) to M_N, M_Delta",
        C2a="rigid-rotor spin/isospin quantisation (0 extra parameters), AN units; d: J=1 I=0 -> hbar^2/V_11; "
            "3He: J=I=1/2 with K+L=0 (assumed, as the hedgehog) -> (3/4)(u+v+2w)/(2(uv-w^2)) hbar^2; 4He: J=I=0 -> 0; "
            "nucleon M_N = E_1 + 3 hbar^2/(8 Lambda)",
        C2b="zero-point energy of the BOUND vibrational modes (omega < m_pi) of B=1 and B=2 (0 extra parameters), "
            "generalised eigenproblem K xi = omega^2 2M xi on the lattice; continuum (Casimir) NOT included",
        C3a="Gillard-Harland-Speight lightly bound potential, 1 extra parameter alpha (m=1 as GHS): scan alpha = 0, 0.5, 0.8, 0.95",
        C3b="sextic term c6 B0^2 added to the AN massive model, 1 extra parameter c6: c6 = 50, 200 (Skyrme units)",
    ),
    exams={
        "C1a": "AN refit reproduces Adkins-Nappi F_pi = 108 MeV and e = 4.84 within 2 %",
        "C1b": "lattice E_1 (AN m) within 1 % of the 1D hedgehog",
        "C1c": "shapes kept with physical pion mass: B=2 torus (asym C_ax 30/45/90 deg < 0.02, one hole), "
               "B=3 tetrahedron (C2x,C2y,C2z,C3 asym < 0.02, inversion asym > 0.3), B=4 cube (Cz90,C3,inv asym < 0.02)",
        "C1d": "PREDICTION 'pion mass alone does not fix binding': classical d binding > 10 x 2.22 MeV AND 4He > 3 x 28.3 MeV",
        "C1e": "literature check m=1 (alpha=0 point of the C3a scan): B=4 cube E/(12 pi^2 B) = 1.307 within 1 %",
        "C2a-i": "lattice B=1 (AN) inertia: U, V, -W diagonal within 3 % of the 1D Lambda",
        "C2a-ii": "torus tensors: U11=U22, V11=V22 within 1 %; |V33|/U33 = 4 and |W33|/U33 = 2 within 3 %; |W11|/U11 < 0.02",
        "C2a-iii": "PREDICTION: rigid quantisation with M_N fitted to 939 MeV makes binding LARGER than classical for d and 4He",
        "C2a-iv": "deuteron rigid-rotor binding within a factor 2 of Braaten-Carson's 158 MeV",
        "C2b-i": "zero modes found: 6 for B=1 and 8 for B=2 with omega^2 < 0.05 m^2 (lattice breaks translations slightly)",
        "C2b-ii": "PREDICTION: zero-point of bound modes REDUCES the d binding (Delta_zp > 0) but d binding after C2a+C2b stays > 10 x 2.22 MeV",
        "C3a-i": "GHS alpha=0.95 1D E_1/8pi^2 = 1.1082 within 1 %; lattice binding per nucleon B=2, B=4 = 1.84, 4.56 MeV within 30 % (E_1 := 938.27 MeV)",
        "C3a-ii": "binding fraction of B=2 and B=4 decreases monotonically with alpha (0, 0.5, 0.8, 0.95)",
        "C3a-iii": "PREDICTION shapes: alpha <= 0.5 torus and cube kept; alpha = 0.95 B=2 and B=4 are separate lumps "
                   "(B=2 not a ring: ring_center_over_max > 0.3 or 2 peaks; B=4 not a cube: inversion asym > 0.1 or 4 peaks)",
        "C3a-iv": "ONE-PARAMETER FIT then PREDICTION: alpha* fixed so that 4He binding = 28.3 MeV (log-linear interpolation in "
                  "the scan, declared fit); predicted d binding at alpha* within a factor 2 of 2.22 MeV and 3He within a factor 2 of 7.72 MeV",
        "C3b-i": "PREDICTION: sextic term lowers the binding fraction monotonically but at c6 = 200 the 4He binding stays > 3 x 28.3 MeV (AN units)",
        "C4": "summary table route x (params, d, 3He, 4He, shapes) and recommendation for the agglomerator",
    },
    notes=[
        "Start fields: relaxed massless N1 fields (rm1..rm4) + 2 % random noise (to let symmetric saddles fall); for alpha=0.95 "
        "also product-ansatz starts (pair d=2.5, tetrahedron edge 2.5, attractive channels); the lowest energy is kept.",
        "AN units: E_u = F_pi/(4e), hbar = 2e^2; GHS scan uses GHS calibration E_1(classical) = 938.27 MeV, alpha=0 point idem.",
        "Each script <= 10 min with checkpoints, 2 threads total.",
    ],
)
if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
    json.dump(P, open(os.path.join(HERE, "out", "colla_prereg.json"), "w"), indent=1)
    print("written", len(P["exams"]), "exams")

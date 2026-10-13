"""nucleo1s_3d_prereg.py -> lab/out/nucleo1s_3d_prereg.json : PRE-REGISTRATION of NUCLEO_1S phase 2 (3-D), written and
committed BEFORE the production runs (lab/nucleo1s_3d_run.py).  The phase-1 pre-registration (nucleo1s_prereg.json) is
not touched.

What was already seen in DEVELOPMENT (declared, lab/nucleo1s_3d_dev.py, lab/out/nucleo1s_3d_dev_*.json):
 * the 3-D rule splits into 4 colours (exact, checked to 3e-15) and, inside a colour, the 1s is a quartet of lines with
   ~1/4 of the Laboratorio start state each: a = 0.5: 0.84037, 0.88627, 0.8993, 0.8996; a = 0.4: 0.86073, 0.89119,
   0.89977, 0.89995; a = 0.3: 0.87795, 0.89473, 0.89989 (x2); a = 0.24: 0.88577, 0.89628, 0.89977 (x2);
   a = 0.2: 0.88986, 0.89704, 0.89963 (x2)  (Dirac 0.899052);
 * the k-space content of the 3 lines at a = 0.5 (4 "tastes": compact momenta 000, 011, 101, 110);
 * the spectrum of the Laboratorio's p_exc (fotogrammi.csv): peaks at 0.0132 and 0.059 mc^2/hbar (951 and 213 ticks),
   i.e. B1 below was SEEN (it is a validation, not a blind prediction);
 * the nearest-site factor s0 at a = 0.5 with s1 = 1: s0 = 0 / 0.5 / 1.5 (lines listed in the document);
 * the hybrid glue: J0 = 2..8 with the phase-1 glue ("row"), and the tau-corrected glue (development, see document).
Everything else below (a = 0.16, the learning, the learned lattices, the excited states, isotropy, cost, the final glue
configuration) has NOT been run.
"""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

pre = dict(
    study="NUCLEO_1S phase 2: the 3-D rule near the nucleus (debts #20, #6 in 3-D) and the hybrid glue",
    date="2026-10-13",
    rule="the cubic 3-D zig-zag rule of zz_core.Walk('dirac') / zz_fast.FastDirac3D / src/atomo/zz.rs, UNCHANGED; one colour "
         "simulated (lab/nucleo1s_3d.py, checked equal to FastDirac3D to 3e-15)",
    inputs="m_e, alpha, Z = 60 (and Z = 1, 2 for the glue); Dirac numbers are exams only",
    config=dict(
        Z=60,
        lab=dict(n=48, a=0.5, nt=2400, start="radial 1s of the same a on the cubic lattice, m=+1/2 (as zz.rs Atomo::new)"),
        lab_data="/mnt/project-files/ripresa/dati/atomo_accoppiato",
        learn=dict(n=48, a=0.5, nt=4000, params="s0: factor on -Za/r of the 2 nearest sites of the colour (r0 = sqrt3/2 a); "
                   "s1: factor on the 6 next sites (r1 = sqrt11/2 a)", criterion="minimise the weighted std of the lines "
                   "(weight >= 0.03, omega in [0.80, 0.95]) of the Laboratorio 1s start state: restore the degeneracy of the "
                   "4 tastes that the continuum has. NO Dirac number is used.", start=[0.6, 1.0],
                   simplex=[[0.6, 1.0], [0.8, 1.0], [0.6, 1.3]], maxfev=40),
        scan=dict(box=24.0, T=500.0, a_bare=[0.5, 0.4, 0.3, 0.24, 0.2, 0.16], a_learned=[0.5, 0.4, 0.3, 0.24, 0.2]),
        excited=dict(box=32.0, T=1500.0, a=[0.5, 0.4], states=["2s", "2p1/2 m=1/2", "2p3/2 m=3/2", "2p3/2 m=1/2"],
                     start="large components of the radial levels of the same a (small component 0), superposed"),
        glue=None,       # filled below
    ),
    exams=dict(
        B1=dict(what="Laboratorio lattice (48^3, a = 0.5), bare rule: the radial-1s start state contains >= 3 lines with "
                     "weight >= 0.1 in [0.80, 0.92]; the 2 strongest peaks of the Laboratorio p_exc(t) coincide with "
                     "differences of these lines", tol=0.03, note="SEEN in development: validation"),
        B2=dict(what="the 4 colours have the same lines (weight >= 0.05)", tol=1e-5),
        B3=dict(what="taste splitting law, NEW lattice a = 0.16 (bare): line distances U - L1 and U - L2 (U = the "
                     "nucleus-blind upper pair) predicted from the dev power laws", pred_U_L1=0.00630, tol_U_L1=0.15,
                pred_U_L2=0.00182, tol_U_L2=0.25),
        N1=dict(what="learning at a = 0.5 succeeds: weighted std of the 1s lines <= max_ratio x bare", max_ratio=0.10),
        N2=dict(what="the SAME 2 numbers at a = 0.4, 0.3, 0.24, 0.2: std <= max_ratio x bare std at the same a",
                max_ratio=0.25),
        N3=dict(what="1s 3-D within 1 % of Dirac (binding of the weighted centroid of the 1s lines) with the learned "
                     "nucleus, at EVERY a in 0.5, 0.4, 0.3, 0.24, 0.2", tol=0.01),
        N4=dict(what="the -11.6 % cured: learned centroid at a = 0.4, 0.3, 0.24 extrapolated as in ZIGZAG A1 "
                     "(quadratic in a^2 through the 3 points) within 1 % of Dirac", tol=0.01),
        N5=dict(what="the beat disappears: Laboratorio lattice, learned nucleus: the largest difference between lines of "
                     "weight >= 0.05 of the start state is <= 2 pi / (6400 ticks x 0.5) = 0.00196 (beat period >= twice "
                     "the Laboratorio run of 3200 ticks)", tol=0.00196),
        E1=dict(what="excited states do not get worse: for 2s, 2p1/2, 2p3/2 (box 32, a = 0.5 and 0.4), the strongest line "
                     "of each start state, learned vs bare: |E_learned - E_Dirac| <= |E_bare - E_Dirac| + 0.005 x binding",
                tol=0.005),
        I1=dict(what="isotropy: 2p3/2 m = 3/2 and m = 1/2 (z) strongest lines differ by <= 2 % of the 2p3/2 binding at "
                     "a = 0.4 (learned), and learned is not more anisotropic than bare by more than 10 %", tol=0.02),
        C1=dict(what="cost: cells (full lattice, box 24) for the 1s within 1 %: learned (the coarsest a passing N3) vs "
                     "bare uniform (a where the bare centroid error, fitted as K a^{2 gamma} + B a^2 on the bare scan, is "
                     "1 %): ratio >= 10", min_ratio=10.0),
        G1=None, G2=None,
    ),
    parameters="2 learned numbers (s0, s1) for Z = 60; 0 fitted to Dirac",
)
glue = dict(method="phase-1 hybrid (zig-zag zoom inside, learned agglomerate B = 256 outside, s wave) with the tau "
                   "(defect) correction: in the overlap the coarse equations get the defect tau = K(w u_in) - T_c(w u_in) of "
                   "the inner zig-zag solution u_in (w: C-infinity window from J0 to J0 + 12 cells; T_c = -u''/(2 m_kin) "
                   "with m_kin of the agglomerate); ghosts j < J0; energy from the coarse equation of row J0 - 1",
            B=256, vmode="point", J0=[3, 4, 5, 6], window=12, states=[[1, 1], [1, 2], [2, 1]])
pre["config"]["glue"] = glue
pre["exams"]["G1"] = dict(what="glue: the Z = 1 1s energy moves by <= 0.02 % over J0 = 3, 4, 5, 6 (max - min)", tol=2e-4)
pre["exams"]["G2"] = dict(what="glue: Z = 1 1s within 0.1 % of Dirac (reduced mass as HY1) for every J0; Z = 1 2s and "
                               "He+ 1s move by <= 0.05 % over J0", tol_1s=1e-3, tol_other=5e-4)
json.dump(pre, open(os.path.join(OUT, "nucleo1s_3d_prereg.json"), "w"), indent=1)
print("written")

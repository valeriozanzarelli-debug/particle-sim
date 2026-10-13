"""nucleo1s_f4_prereg.py -> lab/out/nucleo1s_f4_prereg.json : PRE-REGISTRATION of phase 4 of NUCLEO_1S (written and
committed BEFORE the runs). Thresholds of phases 1-3 are NOT changed (N3, C1, M2, M5 stay FAIL as judged).

(1) TWO LATTICES + TRANSFER LAW. Phase 3 learned (s0, s1) at a = 0.5 with the internal criterion "one line AND at the
nucleus-blind level U of the bare rule" and used the same numbers at every a: -0.95/-1.02/-0.88/-0.72 % at
a = 0.4/0.3/0.24/0.2 (M2 FAIL by 0.017 points). Phase 4: learn (s0, s1) with the SAME criterion at a second lattice,
a = 0.3 (U from scan_bare_0.3), and transfer to the other a with a law fixed in advance, with no new knob:
    s_i(a) = s_i(0.5) + [s_i(0.3) - s_i(0.5)] * (g(a) - g(0.5)) / (g(0.3) - g(0.5)).
Why g(a) = a (MAIN, judged): the Coulomb contact problem has no scale (phase 1: the contact error is K a^{2 gamma} with
K = +1.640 the same at every a, computed without fit); in units of the cell, the only thing that changes with a near the
nucleus is the mass term theta = m a. The nucleus factors cancel the contact phase mismatch at the trained theta, so
what is left is analytic in theta: to first order linear in a. The a^{2 gamma} law then says how the residue reaches the
energy (it is weighted by a^{2 gamma}, hence small at small a).
g(a) = a^{2 - 2 gamma} (transferP): the ratio between bulk (B a^2) and contact (K a^{2 gamma}) error; REPORTED, NOT JUDGED.
Test lattices (never used for training): a = 0.4 (inside), 0.24 and 0.2 (outside). a = 0.16: reported, not judged.
Why 0.3 as second lattice: it brackets 0.4 and is cheap enough to learn (n = 80) within the time budget; the listening
length at 0.3 is that of the scans (T = 500), not that of learn2 (T = 2000): declared.

(2) bcc onion junction (doc section 13.3): only if time remains; declared as a development prototype, no exam here.
"""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = dict(
    study="NUCLEO_1S phase 4", date="2026-10-13", Z=60,
    keeps="all thresholds of nucleo1s_prereg, nucleo1s_3d_prereg, nucleo1s_3d_amend1, nucleo1s_f3_prereg unchanged",
    config=dict(
        learn3=dict(n=80, a=0.3, nt=1667, U="as learn2, from results['scan_bare_0.3']",
                    objective="std^2 + (centroid - U)^2 (as learn2)",
                    start=[0.5546, 1.0035], simplex=[[0.5546, 1.0035], [0.5646, 1.0385], [0.5646, 1.0035]], maxfev=30),
        first_lattice="learn2 of phase 3 (a = 0.5, s = (0.554593, 1.003549)), unchanged",
        law_main="transfer: g(a) = a", law_reported="transferP: g(a) = a^(2 - 2 gamma), gamma = sqrt(1 - (Z alpha)^2)",
        scan=dict(box=24.0, T=500.0, a_test=[0.4, 0.24, 0.2], a_reported=[0.16]),
        excited=dict(cfg="transfer", box=32.0, T=1500.0, a=[0.4]),
    ),
    knobs="2 learned numbers per lattice x 2 lattices = 4 numbers, fixed by the internal criterion; the law has no knob",
    exams=dict(
        Q1=dict(what="training at a = 0.3: spread <= 0.002 mc^2 and |centroid - U(0.3)| <= 5e-4 mc^2", spread=0.002, dU=5e-4),
        Q2=dict(what="MAIN: with the transfer law g = a, at every TEST lattice a = 0.4, 0.24, 0.2 the 1s centroid binding is "
                     "within 1 % of Dirac (omega 0.899052) AND the spread (lines weight >= 0.05 in [0.80, 0.95]) is <= 0.002",
                tol=0.01, spread=0.002),
        Q3=dict(what="excited states not worse at a = 0.4 (box 32), transfer law: mean |binding error| of the 4 excited "
                     "starts <= bare (rule of M4)"),
        Q4=dict(what="cost, only if Q2 passes: cells (box 24) at the coarsest TEST lattice that passes vs 6.1e6 of the bare "
                     "uniform lattice: ratio >= 10", min_ratio=10),
    ),
    reported_not_judged=["transferP at a = 0.4, 0.24, 0.2", "transfer at a = 0.16"],
)
json.dump(P, open(os.path.join(OUT, "nucleo1s_f4_prereg.json"), "w"), indent=1)
print("written")

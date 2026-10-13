"""nucleo1s_f3_prereg.py -> lab/out/nucleo1s_f3_prereg.json : PRE-REGISTRATION of phase 3 of NUCLEO_1S (written and
committed BEFORE the runs). The thresholds of phase 1/2 (nucleo1s_prereg, nucleo1s_3d_prereg, nucleo1s_3d_amend1) are NOT
changed: N3 and C1 stay FAIL as judged in phase 2. Phase 3 adds NEW exams M1..M5 for a NEW learned nucleus.

(1) SECOND INTERNAL CRITERION for the learned nucleus (doc section 12.5). Phase 2 found that the criterion "the 4 lines of
the 1s multiplet become one" leaves a 1-parameter line of solutions (s0, s1), along which the merged level runs from
0.8976 to 0.9049 (Dirac 0.899052). The optimiser stopped at one end (0.9049: -5.8 % of the binding).
Second criterion, without Dirac: the merged level must sit where the bare rule puts the NUCLEUS-BLIND pair U (the two
levels of the bare multiplet that have ~0 density on the 2 sites nearest to the nucleus, so they are not spoiled by the
hand-made point charge between sites). U is a number of the RULE (bare scan, a = 0.5), not of Dirac.
HONEST DECLARATION: we already know (phase 2, V2) that U is within 1 % of Dirac at every a; so this criterion was chosen
AFTER seeing that. It is internal (it does not read Dirac), but it is not blind. What is not known yet, and is the exam:
whether the 2 numbers learned at a = 0.5 keep the merged level within 1 % at the other a (0.4 ... 0.2), where U itself
moves differently (bare U: -0.40, -0.80, -0.83, -0.71, -0.57 % at a = 0.5, 0.4, 0.3, 0.24, 0.2).
Objective = std^2 + (centroid - U)^2 (same lines as phase 2: weight >= 0.03, omega in [0.80, 0.95]); Nelder-Mead,
start on the degenerate line where the phase-2 evaluations crossed U.

(2) 3-D ONION with an inner shell that respects the 4 body diagonals: DESIGN ONLY in this phase (no 3-D refinement
interface of the zig-zag exists yet that keeps the 4 colours separated and is unitary); declared NOT RUN, no exam.
(3) Z = 1 in 3-D: run only if M2 passes AND the cost fits; a0 = 137 lambda_C needs a box of ~1000 lambda_C, i.e.
(2000)^3 cells at a = 0.5: does not fit (<= 2.5 GB). Declared NOT RUN unless a 3-D hybrid exists; no exam.
"""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = dict(
    study="NUCLEO_1S phase 3", date="2026-10-13", Z=60,
    keeps="all thresholds of lab/out/nucleo1s_prereg.json, nucleo1s_3d_prereg.json, nucleo1s_3d_amend1.json unchanged",
    config=dict(
        learn2=dict(n=48, a=0.5, nt=4000,
                    U="weighted mean of the highest lines of results['scan_bare_0.5'] until their weight sum >= 0.4 "
                      "(the nucleus-blind pair)",
                    objective="std^2 + (centroid - U)^2 of the 1s lines of the Laboratorio start state (weight >= 0.03, "
                              "omega in [0.80, 0.95])",
                    start=[0.55, 1.01], simplex=[[0.55, 1.01], [0.58, 1.06], [0.52, 0.95]], maxfev=40),
        scan=dict(cfg="learned2", box=24.0, T=500.0, a=[0.5, 0.4, 0.3, 0.24, 0.2],
                  note="same machinery as the phase-2 scans (part_scan)"),
        excited=dict(cfg="learned2", box=32.0, T=1500.0, a=[0.4]),
    ),
    knobs="2 numbers (s0, s1), same as phase 2; U is not a knob (it is read from the bare rule)",
    exams=dict(
        M1=dict(what="learning: at a = 0.5 the best point has 1s spread (max - min of lines) <= 0.002 mc^2 AND "
                     "|centroid - U| <= 5e-4 mc^2 (both internal criteria met together)", spread=0.002, dU=5e-4),
        M2=dict(what="MAIN EXAM: with the same 2 numbers, at a = 0.4, 0.3, 0.24, 0.2 the 1s centroid binding is within "
                     "1 % of Dirac (omega 0.899052, binding 0.100948) AND the spread is <= 0.002 mc^2 "
                     "(a = 0.5 reported, not judged: it is the training lattice)", tol=0.01, spread=0.002),
        M3=dict(what="the Laboratorio beat is closed: at a = 0.5 the spread of the learned2 1s is <= 1/10 of the bare "
                     "Laboratorio spread (0.0196)", max_spread=0.00196),
        M4=dict(what="excited states not worse: at a = 0.4 (box 32) with learned2, the mean |binding error| of "
                     "2s, 2p1/2, 2p3/2 is <= that of the bare nucleus (same rule as E1)"),
        M5=dict(what="cost: if M2 passes, cells for the 1s within 1 % with learned2 at the coarsest passing a vs "
                     "bare uniform (6.1e6 cells, phase 2): ratio >= 10", min_ratio=10),
    ),
    not_run=dict(onion3d="design only (doc section 13)", Z1_3d="not run: does not fit; see doc section 13"),
)
json.dump(P, open(os.path.join(OUT, "nucleo1s_f3_prereg.json"), "w"), indent=1)
print("written")

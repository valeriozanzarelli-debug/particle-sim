"""arrivo_emendamenti.py - writes lab/out/arrivo_emendamenti.json: changes made AFTER the pre-registration."""
import json
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
em = [
    dict(id="A1", when="after NR production (NR_A, NR_B)",
         what="tau_max of the NR reference: the 50 latest (and the not-arrived) trajectories re-integrated with dt/4 "
              "(checked dt/10: same to 1e-3). With the pre-registered steps 2 of 20000 trajectories (NR_A) and 1 of 4000 "
              "(NR_B) were integrated badly near a node (arrival 43.69 / not arrived / 20.99): with dt/4 they arrive at "
              "8.8, 16.3 and 1.18. The raw sample max (43.69) would PASS E1, the corrected one (41.52) FAILS: we judge "
              "E1 with the corrected value.",
         effect="E1 FAIL instead of PASS; E3 recomputed with the corrected tau_max (still PASS)"),
    dict(id="A2", when="during the 1-D production",
         what="hybrid surfers (Heun) with the pre-registered NR step schedule (2e-4, 2e-3, 1e-2) are NOT converged for "
              "omega = 1000 (mean -2.6 %; test on 1000 samples, D = 40, a = 0.125: 2.7239 / 2.7955 / 2.7982 / 2.7981 for "
              "steps x 1 / 0.2 / 0.1 / 0.05). All 'x' runs redone with steps x 0.1 = (2e-5, 2e-4, 1e-3).",
         effect="E7, E8 judged on the redone runs"),
    dict(id="A3", when="during the 1-D production (bug)",
         what="alpha = 0 hybrid runs: the transverse rotation X, Y at rate omega blew up in the explicit step "
              "(inf * sin 0 = nan): X = Y = 0 for alpha = 0 (it does not touch Z). z runs started before the fix redone.",
         effect="E9 judged on redone runs; the flux records were not affected"),
    dict(id="A4", when="analysis",
         what="E13: besides (1 - v)/(1 + v) with the continuum v we also report it with the lattice group velocity; "
              "the verdict uses the pre-registered continuum formula.",
         effect="none"),
]
json.dump(em, open(os.path.join(OUT, "arrivo_emendamenti.json"), "w"), indent=1)
print("written")

"""surf_amend.py -> lab/out/surf_amend.json : changes declared AFTER the pre-registration (no threshold changed)."""
import json
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A = dict(
    date="2026-10-13",
    changes=[
        dict(id="A1", when="before production", what="run_orb split in checkpointed pieces; old orbital reader with 20000 "
             "surfers instead of 50000 (report-only scheme), to fit the 10-minute rule"),
        dict(id="A2", when="during production", what="run_beat was killed by the 590 s timeout in its first (non "
             "checkpointed) attempt; checkpointing added, same seeds and code path, rerun from tick 0"),
        dict(id="A3", when="after production", what="TV ratio (E5-E7): the code divided by the ideal-sampler TV computed "
             "at t = 0; the pre-registration defines it against an ideal sampler of the same density. surf_tvnull.py "
             "recomputes the reference at the final density (same N, same blocks). Both numbers are reported."),
        dict(id="A4", when="after production", what="E5 v_phi FAILS (0.21 instead of 0.75-1.25). Diagnosis after the run: "
             "the listened 1s of the n = 48, a = 0.5 lattice has <Sigma_z> = 0.51 (exact Dirac 1s: ~0.99) and its own "
             "net current gives v_phi/(Za sin) = 0.23; the Bohm-Dirac formula on the same lattice spinor gives 0.36. "
             "The failure is of the lattice state at theta = 0.5, not of the surfer rule. Threshold NOT changed."),
        dict(id="A5", when="after production", what="explorations not pre-registered, for the section 'Idee geometriche "
             "e visive': surf_hodge.py (Hodge split of the two-tick flux), surf_macro.py (exact surfers on 2^3 blocks)"),
    ],
)
with open(os.path.join(OUT, "surf_amend.json"), "w") as f:
    json.dump(A, f, indent=1)
print("written lab/out/surf_amend.json")

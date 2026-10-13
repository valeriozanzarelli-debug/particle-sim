"""surf_f3_amend.py -> lab/out/surf_f3_amend.json : PHASE 3 changes declared AFTER lab/out/surf_f3_prereg.json."""
import json
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A = dict(date="2026-10-13", changes=[
    dict(id="C1", when="after production", what="the surfer window starting at tick 3990 was never measured: the run loop "
         "stops at t = 4000 before the window closes (bug in surf_f3.py run). F3 and F5 are judged on the 4 measured "
         "windows (500, 1000, 2000, 3000). Thresholds unchanged. The state's own flux v_phi at tick 3900 (0.834) is reported "
         "as the late value."),
])
json.dump(A, open(os.path.join(OUT, "surf_f3_amend.json"), "w"), indent=1)
print("written lab/out/surf_f3_amend.json")

"""surf_f4_amend.py -> lab/out/surf_f4_amend.json : PHASE 4 additions declared AFTER lab/out/surf_f4_prereg.json."""
import json
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A = dict(date="2026-10-13", changes=[
    dict(id="D1", when="after production", what="POST-HOC fit (not an exam, not counted): d(a) = d0 + c a^(2 gamma), gamma = "
         "sqrt(1 - Za^2) = 0.899 (the contact law of NUCLEO_1S): d0 = 0.0058, R^2 = 0.99998. Chosen after seeing that the "
         "pre-registered linear law fails G6 by 0.0008 while d0 + c a^2 gives 0.013."),
])
json.dump(A, open(os.path.join(OUT, "surf_f4_amend.json"), "w"), indent=1)
print("written lab/out/surf_f4_amend.json")

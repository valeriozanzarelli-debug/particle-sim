"""big_r1_audit.py - which old geo_ claims (docs/FORMA_DALLE_PALLINE.md) survive the repaired reader?
Applies the R1 criterion (H1 persistence / diameter > 0.4, full filtration, landmarks) to the stored geodesic
matrices G0 (piece 0) of lab/out/geo_space_{vac,quench}_*.npz, and big_core.read_shape to the stored vacuum MI.
Not graded (pre-registered as an audit in big_prereg_R.json). Writes lab/out/big_r1_audit.json (~20 s)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import big_core as B

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
OLD_B1 = {"vac": {"chain48": 0, "ring48": 1, "torus16": 2, "cyl36x18": 2, "stadium30x16": 1},
          "quench": {"chain48": 0, "ring48": 1, "torus16": 2, "cyl36x18": 2, "stadium30x16": 4}}
DRAWN_B1 = {"chain48": 0, "ring48": 1, "torus16": 2, "cyl36x18": 1, "stadium30x16": 1}
res = {}
for mode in ("vac", "quench"):
    for s in ("chain48", "ring48", "torus16", "cyl36x18", "stadium30x16"):
        z = np.load(os.path.join(OUT, f"geo_space_{mode}_{s}.npz"))
        rel, diam, _ = B._h1_rel(z["G0"])
        e = dict(old_b1=OLD_B1[mode][s], drawn_b1=DRAWN_B1[s], new_holes_G0=int((rel > B.THR).sum()),
                 rel_top3=[float(x) for x in rel[:3]])
        if mode == "vac":
            r = B.read_shape(z["raw"], with_mds=False)
            e["new_holes_from_MI"] = r["holes"]; e["pieces"] = r["pieces"]
        res[f"{mode}_{s}"] = e
        print(mode, s, e, flush=True)
with open(os.path.join(OUT, "big_r1_audit.json"), "w") as f:
    json.dump(res, f, indent=1)

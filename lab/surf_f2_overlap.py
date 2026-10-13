"""surf_f2_overlap.py - POST-HOC (declared in surf_f2_amend.json): is the E5 state (listened at 0.89945, <Sigma_z> = 0.51)
a mixture of the near-degenerate exact levels 0.89930 / 0.89960? Overlaps <level|E5>, and Sigma_z of the reconstruction.
-> lab/out/surf_f2_overlap.json"""
import json
import numpy as np
from surf_core import OUT
import os
from surf_f2_valli import sigma_z, valley_split, CK

files = sorted(f for f in os.listdir(CK) if f.startswith("f2_level_"))
L = {f[9:-4]: np.load(os.path.join(CK, f)) for f in files}
L = {k: v / np.linalg.norm(v) for k, v in L.items()}
e5 = np.load(os.path.join(CK, "state_1s.npy")); e5 = e5 / np.linalg.norm(e5)
ov = {k: complex(np.vdot(v, e5)) for k, v in L.items()}
gram = {f"{a}|{b}": abs(complex(np.vdot(L[a], L[b]))) for a in L for b in L if a < b}
rec = sum(ov[k] * L[k] for k in L)
res = dict(overlap_sq={k: abs(c)**2 for k, c in ov.items()}, captured=float(np.linalg.norm(rec)**2), gram_offdiag=gram)
rn = rec / np.linalg.norm(rec)
res["sigma_z_E5"] = sigma_z(e5); res["sigma_z_reconstruction"] = sigma_z(rn)
pair = ov["0.89930"] * L["0.89930"] + ov["0.89960"] * L["0.89960"]
res["captured_pair"] = float(np.linalg.norm(pair)**2)
pn = pair / np.linalg.norm(pair)
res["sigma_z_pair_only"] = sigma_z(pn)
res["valleys_pair_only"] = {k: dict(weight=v["weight"], sigma_z_lab=v["sigma_z_lab"]) for k, v in valley_split(pn).items()}
print(json.dumps(res, indent=1, default=float))
json.dump(res, open(os.path.join(OUT, "surf_f2_overlap.json"), "w"), indent=1, default=float)

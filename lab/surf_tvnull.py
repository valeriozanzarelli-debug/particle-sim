"""surf_tvnull.py - AMENDMENT (declared after the runs, see lab/out/surf_amend.json): the TV ratio of the 3-D runs was
divided by the TV of an ideal sampler computed at t = 0, while the pre-registration defines it against an ideal sampler
of the SAME density. For the moving/spreading states the reference must be taken at the final density. This script
recomputes the ideal-sampler TV at the final tick (same N, same 2^3 blocks) and writes lab/out/surf_tvnull.json.
No threshold is changed."""
import json
import numpy as np
from surf_core import *
from surf_run import setup_3d

out = {}
rng = np.random.default_rng(SEED + 77)
for state, T in (("free", 300), ("1s", 1000), ("1s2p", 1670)):
    W, psi = setup_3d(state)
    for _ in range(T):
        psi, _ = W.tick(psi)
    rho = density(psi)
    d = json.load(open(os.path.join(OUT, f"surf_a3d_{state}.json")))
    n200 = tv_null(rng, rho, 200000, 2, 10)
    n50 = tv_null(rng, rho, 50000, 2, 10)
    ratio = {k: d["rec"]["tv"][k][-1] / (n200[0] if k in ("flow", "flow2") else n50[0]) for k in d["rec"]["tv"]}
    out[state] = dict(tv_null_final_200k=n200, tv_null_final_50k=n50, tv_ratio_final=ratio)
    print(state, {k: round(v, 3) for k, v in ratio.items()}, flush=True)
json.dump(out, open(os.path.join(OUT, "surf_tvnull.json"), "w"), indent=1)

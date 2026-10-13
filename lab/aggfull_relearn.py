"""aggfull_relearn.py - AMENDMENT run (lab/out/aggfull_amend.json): the free generator re-learned with all samples where the
observable lives (uniform share weight 1e-3 instead of 0.03), w = 4; B = 256 from the fine walk, B = 512 from the re-learned B = 256. Tameness and
the usual checks are evaluated (agg_zz_core.kcheck, state mode). Writes lab/out/aggfull_relearn.pkl/.json"""
import json, os, pickle
import numpy as np
import aggfull_core as ac
import agg_zz_core as azc
from agg_zz_core import Level

TH = ac.TH_A
rng = np.random.default_rng(20261014)
prev = Level("fine3", 3, "dirac4", 4, 1, TH, period=2)
out, js = [], []
# uniform weight per level (amendment 2, aggfull_amend2.json): the largest of 1e-3, 3e-4, 1e-4, ... whose
# state-distributed <T - k^2/2m_kin> agrees with the FINE walk within 5 % (B = 256: 1e-3; B = 512: 3e-4)
WU = {256: 1e-3, 512: 3e-4}
for B, kfac in ((256, 1 / ac.A0), (512, 1 / (2 * ac.A0))):
    b = B // prev.B
    r = azc.learn_h(prev, b, "pauli2", 4, -kfac * TH * B, rng, nk=8000, frac_uniform=0.25, w_uniform=WU[B])
    cand = Level("learned", 3, "pauli2", 2, B, TH, rule=r, b=b, w=4)
    ck = azc.kcheck(prev, cand, b, rng, 1000, kfac * TH, sig=max(1.0, 2.0 / (kfac * TH * B)), k0_cells=0.5 * kfac * TH * B,
                    energy_mode="state")
    ok = ck["packet"] <= 1e-3 and ck["energy_rel"] <= 3e-4 and ck["rest_rel"] <= 3e-4 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01
    print(f"B={B} (b={b} from B={prev.B}): " + ", ".join(f"{k} {v:.2e}" for k, v in ck.items() if isinstance(v, float)) + f"  ok={ok}", flush=True)
    out.append(dict(B=B, b=b, w=4, w_uniform=WU[B], desc="pauli2", C=2, K=r.K, offs=r.offs, phi0=r.phi0, check=ck, ok=bool(ok)))
    js.append(dict(B=B, b=b, w=4, w_uniform=WU[B], check=ck, ok=bool(ok)))
    cand.check = ck
    prev = cand
pickle.dump(out, open(os.path.join(ac.OUT, "aggfull_relearn.pkl"), "wb"))
json.dump(js, open(os.path.join(ac.OUT, "aggfull_relearn.json"), "w"), indent=1, default=float)

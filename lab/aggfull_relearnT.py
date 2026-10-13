"""aggfull_relearnT.py - PHASE 4: the free generator of the a0/16 and a0/8 levels learned with the Taylor-pinned + parity
learner of phase 3 (aggfull_l1.fit mode "T"; constraints c2, c4, c6 along 13 directions read off the level below), and
the Markov version (mode "M") whose infeasibility at these cells is predicted (|c4|/c2 >= 1/12 for positive hops vs
1/(4h^2) needed). B = 256 from the fine walk (b = 256), B = 512 from the learned B = 256 (b = 2, repeat).
Width w: the smallest of (4, 5) passing the usual check (agg_zz_core.kcheck, state mode; budget as aggfull_relearn).
Writes lab/out/aggfull_relearnT.pkl (not committed) and lab/out/aggfull_relearnT.json.  Usage: python aggfull_relearnT.py [dev]"""
import json, os, sys, pickle, time
import numpy as np
import aggfull_core as ac
import aggfull_l1 as L
import agg_zz_core as azc
from agg_zz_core import Level

dev = len(sys.argv) > 1 and sys.argv[1] == "dev"
TH = ac.TH_A
D = L.dirs13()
prev = Level("fine3", 3, "dirac4", 4, 1, TH, period=2)
out, js = [], []
for B, kfac, wu in ((256, 1 / ac.A0, 1e-3), (512, 1 / (2 * ac.A0), 3e-4)):
    b = B // prev.B
    t0 = time.time()
    cons = L.taylor6(prev, b, D, kmax=min(0.224 * B * TH, 1.5))     # same physical k range as level 1 (k <= 0.224 m)
    rng = np.random.default_rng(20261015)
    kap, Gs, sw, phi0 = L.samples(prev, b, kfac * TH * B, rng, 8000, 0.25, wu)
    chosen = None
    rec = dict(B=B, b=b, w_uniform=wu, c6_over_dirac=float(np.mean(cons[:, 2] / (1 / (16 * (B * TH) ** 5)))), tried=[])
    for w in (4, 5):
        r = L.fit(kap, Gs, sw, phi0, w, "T", cons=cons, D=D)
        cand = Level("learned", 3, "pauli2", 2, B, TH, rule=r, b=b, w=w)
        ck = azc.kcheck(prev, cand, b, np.random.default_rng(B + w), 1000, kfac * TH, sig=max(1.0, 2.0 / (kfac * TH * B)),
                        k0_cells=0.5 * kfac * TH * B, energy_mode="state")
        ok = ck["packet"] <= 1e-3 and ck["energy_rel"] <= 3e-4 and ck["rest_rel"] <= 3e-4 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01
        rec["tried"].append(dict(w=w, check=ck, ok=bool(ok)))
        print(f"B={B} T w={w}: " + ", ".join(f"{k} {v:.2e}" for k, v in ck.items() if isinstance(v, float)) + f" ok={ok} ({time.time() - t0:.0f}s)", flush=True)
        if ok:
            chosen = (w, r, cand)
            break
    # Markov at this cell: predicted infeasible
    rM = L.fit(kap, Gs, sw, phi0, 4, "M", cons=cons[:, :2], D=D)
    g, o = L.band_of(rM)
    c2 = np.array([-0.5 * np.sum(g * (o @ d) ** 2) for d in D]); c4 = np.array([np.sum(g * (o @ d) ** 4) / 24 for d in D])
    rec["markov_w4"] = dict(max_c2_dev=float(np.max(np.abs(c2 / cons[:, 0] - 1))), max_c4_dev=float(np.max(np.abs(c4 / cons[:, 1] - 1))),
                            bound_ratio_needed=float(np.mean(np.abs(cons[:, 1]) / cons[:, 0])), bound_min_positive=1 / 12)
    print(f"B={B} Markov w=4: c2 dev {rec['markov_w4']['max_c2_dev']:.2e} c4 dev {rec['markov_w4']['max_c4_dev']:.2e} "
          f"(|c4|/c2 needed {rec['markov_w4']['bound_ratio_needed']:.4f}, positive hops >= 1/12)", flush=True)
    if chosen is None:
        print(f"B={B}: no T candidate passes", flush=True)
        js.append(rec)
        break
    w, r, cand = chosen
    rec.update(w=w)
    out.append(dict(B=B, b=b, w=w, desc="pauli2", C=2, K=r.K, offs=r.offs, phi0=r.phi0, check=rec["tried"][-1]["check"], ok=True))
    js.append(rec)
    prev = cand
tag = "dev_" if dev else ""
pickle.dump(out, open(os.path.join(ac.OUT, tag + "aggfull_relearnT.pkl"), "wb"))
json.dump(js, open(os.path.join(ac.OUT, tag + "aggfull_relearnT.json"), "w"), indent=1, default=float)
print("saved", flush=True)

"""aggfull_p5.py - PHASE 5 (lab/out/aggfull_p5_prereg.json), level-1 part (M2c): Taylor-pinned k^2, k^4, k^6 + parity, Z = 10
observable, width chosen by the agglomerator's own rule (agg_zz_core.auto_level: smallest w passing the WHOLE check, energy
budget 3e-5 included) among w = 4..7, choice made with check sample seed+100; the EXAM is an independent check (sample
seed+200, nk = 10000). Seeds 1..6. Writes lab/out/aggfull_p5.json (checkpoints; rerun until 'P5 L1 DONE')."""
import os, sys, json, time
import numpy as np
import aggfull_l1 as L
import agg_zz_core as azc

OUT = L.OUT
JS = os.path.join(OUT, "aggfull_p5.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
th = L.THETA["A"]; b = B = 16
prev = azc.Level("fine3", 3, "dirac4", 4, 1, th, period=2)
D = L.dirs13()
c246 = np.array(json.load(open(os.path.join(OUT, "aggfull_p4.json")))["c246"])
kf = 10 / L.A0 * th * B


def check(r, w, seed, nk):
    cand = azc.Level("learned", 3, "pauli2", 2, B, th, rule=r, b=b, w=w)
    ck = azc.kcheck(prev, cand, b, np.random.default_rng(seed), 1000, 10 / L.A0 * th, sig=max(1.0, 2.0 / kf), k0_cells=0.5 * kf,
                    nk=nk, energy_mode="state")
    ok = ck["packet"] <= 1e-3 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01 and ck["energy_rel"] <= 3e-5 and ck["rest_rel"] <= 3e-5
    return ck, bool(ok)


t0 = time.time()
for seed in (1, 2, 3, 4, 5, 6):
    key = f"s{seed}"
    if key in res:
        continue
    kap, Gs, sw, phi0 = L.samples(prev, b, kf, np.random.default_rng(seed), 8000, 0.25, 1e-3)
    tried, chosen = [], None
    for w in (4, 5, 6, 7):
        r = L.fit(kap, Gs, sw, phi0, w, "T", cons=c246, D=D)
        ck, ok = check(r, w, seed + 100, 3000)
        tried.append(dict(w=w, energy=ck["energy_rel"], bz=ck["bz_rel90"], neg=ck["neg_gen"], ok=ok))
        print(f"s{seed} w={w}: choice check energy {ck['energy_rel']:.2e} bz {ck['bz_rel90']:.2f} neg {ck['neg_gen']:.1e} ok {ok} ({time.time() - t0:.0f}s)", flush=True)
        if ok:
            chosen = (w, r)
            break
    out = dict(tried=tried, w=None)
    if chosen:
        w, r = chosen
        ck, ok = check(r, w, seed + 200, 10000)
        out.update(w=w, exam_check=ck, exam_ok=ok)
        print(f"s{seed}: chosen w={w}; INDEPENDENT check energy {ck['energy_rel']:.2e} rest {ck['rest_rel']:.1e} bz {ck['bz_rel90']:.2f} ok {ok}", flush=True)
    res[key] = out
    json.dump(res, open(JS, "w"), indent=1, default=float)
    if time.time() - t0 > 420:
        print("checkpoint"); sys.exit(0)
print("P5 L1 DONE", flush=True)

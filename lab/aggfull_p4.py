"""aggfull_p4.py - PHASE 4 (lab/out/aggfull_p4_prereg.json), level-1 part:
  (a) T6: Taylor-pinned up to k^6 + parity, Z = 10 observable, width = the smallest of w = 4, 5, 6 passing the usual check
      (packet <= 1e-3, energy and rest <= 3e-5, bz_rel90 <= 0.2, neg_gen <= 0.01), seeds 1..6 (M2b)
  (b) Markov + k^6 (positive hops, c2 c4 c6 pinned), w = 5, 6, 7, 8 (P3: axis moment bound says w >= 8 needed)
Writes lab/out/aggfull_p4.json (checkpoints; rerun until 'P4 L1 DONE')."""
import os, sys, json, time
import numpy as np
import aggfull_l1 as L
import agg_zz_core as azc

OUT = L.OUT
JS = os.path.join(OUT, "aggfull_p4.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
th = L.THETA["A"]; b = B = 16
prev = azc.Level("fine3", 3, "dirac4", 4, 1, th, period=2)
D = L.dirs13()
if "c246" not in res:
    res["c246"] = L.taylor6(prev, b, D).tolist()
c246 = np.array(res["c246"])
res["c6_over_dirac"] = float(np.mean(c246[:, 2] / (1 / (16 * (B * th) ** 5))))
t0 = time.time()
kf10 = 10 / L.A0 * th * B
for seed in (1, 2, 3, 4, 5, 6):
    key = f"T6_Z10_s{seed}"
    if key in res:
        continue
    kap, Gs, sw, phi0 = L.samples(prev, b, kf10, np.random.default_rng(seed), 8000, 0.25, 1e-3)
    tried = []
    for w in (4, 5, 6):
        r = L.fit(kap, Gs, sw, phi0, w, "T", cons=c246, D=D)
        cand = azc.Level("learned", 3, "pauli2", 2, B, th, rule=r, b=b, w=w)
        ck = azc.kcheck(prev, cand, b, np.random.default_rng(seed + 100), 1000, 10 / L.A0 * th, sig=max(1.0, 2.0 / kf10),
                        k0_cells=0.5 * kf10, energy_mode="state")
        tame = ck["packet"] <= 1e-3 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01
        ok = tame and ck["energy_rel"] <= 3e-5 and ck["rest_rel"] <= 3e-5
        tried.append(dict(w=w, check=ck, tame=bool(tame), ok=bool(ok)))
        print(f"{key} w={w}: energy {ck['energy_rel']:.2e} rest {ck['rest_rel']:.1e} bz {ck['bz_rel90']:.2f} neg {ck['neg_gen']:.1e} ok {ok}", flush=True)
        if tame:
            break
    g, o = L.band_of(r)
    m_rel, m_abs, fr = L.mass_k4(g, o, prev, b, B * th, c246[:, :2])
    res[key] = dict(tried=tried, w=tried[-1]["w"], m_over_mfine_minus_1=m_rel - 1, k4frac=fr)
    json.dump(res, open(JS, "w"), indent=1, default=float)
    if time.time() - t0 > 450:
        print("checkpoint"); sys.exit(0)
kf2 = 2 / L.A0 * th * B
kap, Gs, sw, phi0 = L.samples(prev, b, kf2, np.random.default_rng(1), 8000, 0.25, 1e-3)
for w in (5, 6, 7, 8):
    key = f"M6_w{w}"
    if key in res:
        continue
    t1 = time.time()
    r = L.fit(kap, Gs, sw, phi0, w, "M", cons=c246, D=D)
    g, o = L.band_of(r)
    off = np.any(o != 0, axis=1)
    dev = [float(np.max(np.abs(np.array([s * np.sum(g * (o @ d) ** p) / f for d in D]) / c246[:, i] - 1)))
           for i, (p, f, s) in enumerate(((2, 2, -1), (4, 24, 1), (6, 720, -1)))]
    res[key] = dict(max_dev_c2_c4_c6=dev, wrong=int(np.sum(g[off] > 1e-12 * abs(g[~off][0]))), sec=time.time() - t1)
    json.dump(res, open(JS, "w"), indent=1, default=float)
    print(f"{key}: constraint deviations c2 {dev[0]:.1e} c4 {dev[1]:.1e} c6 {dev[2]:.1e}, wrong taps {res[key]['wrong']} ({time.time() - t1:.0f}s)", flush=True)
    if time.time() - t0 > 450:
        print("checkpoint"); sys.exit(0)
print("P4 L1 DONE", flush=True)

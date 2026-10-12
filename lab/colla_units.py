"""colla_units.py - unit fits and 1D B=1 references for COLLA. Writes lab/out/colla_units.json.
 * ANW (massless) and AN (m_pi = 138 MeV) fits of (e, F_pi) to M_N, M_Delta  [exam C1a]
 * 1D hedgehog E_1 and Lambda for: AN m; m = 1 with alpha = 0, 0.5, 0.8, 0.95 (GHS); AN m with sextic c6 = 50, 200."""
import json, os, math, time
import colla_core as cc
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
t0 = time.time()
anw = cc.fit_units(0.0)
an = cc.fit_units(cc.M_PI)
res = dict(ANW={k: float(v) for k, v in anw.items()}, AN={k: float(v) for k, v in an.items()}, hedgehog={})
print("ANW", res["ANW"]); print("AN", res["AN"])
models = {"an": dict(m=an["m"]), "ghs0": dict(m=1.0), "ghs50": dict(m=1.0, alpha=0.5), "ghs80": dict(m=1.0, alpha=0.8),
          "ghs95": dict(m=1.0, alpha=0.95), "sx50": dict(m=an["m"], c6=50.0), "sx200": dict(m=an["m"], c6=200.0)}
for k, kw in models.items():
    hh = cc.hedgehog_1d(**kw)
    res["hedgehog"][k] = dict(kw=kw, E=hh["E"], E_over_12pi2=hh["E"] / cc.PI2, E_over_8pi2=hh["E"] / (8 * math.pi**2), Lam=hh["Lam"])
    print(k, res["hedgehog"][k], f"{time.time()-t0:.0f}s", flush=True)
res["models"] = models
json.dump(res, open(os.path.join(OUT, "colla_units.json"), "w"), indent=1)

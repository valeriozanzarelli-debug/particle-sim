"""glue_posthoc.py - NOT graded. Diagnostics of the beta = 12 0++ mass (it came out 5.8% above the literature):
chain by chain, and single operators. Writes lab/out/glue_posthoc.json."""
import os, json, glob
import numpy as np
import glue_analysis as A

OUT = A.OUT; RAW = A.RAW
plan = json.load(open(os.path.join(OUT, "glue_prereg.json")))["runs"]["12"]
T = plan["T"]
res = {}
for tag in ("b12", "b12x"):
    files = sorted(glob.glob(os.path.join(RAW, f"{tag}_c*.npz")))
    d = [np.load(f) for f in files]
    names = [str(x) for x in d[0]["names"]]
    ops = np.concatenate([x["ops"] for x in d])
    S, M, n = A.binned_corr(ops, 20)
    C = A.jk_corrs(S, M, n)
    A1 = [i for i, q in enumerate(names) if q.startswith("A1_")]
    r = A.channel(C, A1, T, 1)[0]
    res[tag] = dict(N=len(ops), am0=r["m"], err=r["err"], meff=r["meff"][:6], meff_err=r["meff_err"][:6])
    # single operators: effective mass at t = 2, 3 of each A1 operator alone
    single = {}
    for i in A1:
        Cs = C[:, i, i, :]
        me = [A.jk_err(np.array([A.meff_cosh(Cs[k], t, T) for k in range(len(Cs))])) for t in (1, 2, 3)]
        single[names[i]] = [(float(a), float(b)) for a, b in me]
    best = sorted(single.items(), key=lambda kv: kv[1][1][0] if np.isfinite(kv[1][1][0]) else 9)[:4]
    res[tag]["best_single_ops_meff_t1_t2_t3"] = dict(best)
    print(tag, res[tag]["N"], round(r["m"], 4), round(r["err"], 4), [round(x, 3) for x in r["meff"][:6]])
    for k, v in best:
        print("   ", k, [(round(a, 3), round(b, 3)) for a, b in v])
json.dump(res, open(os.path.join(OUT, "glue_posthoc.json"), "w"), indent=1)

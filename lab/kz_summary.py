"""kz_summary.py - K1 analysis: alpha with bootstrap CI, sensitivity table, pre-registered verdict -> lab/out/kz_k1.json"""
import os, glob, json, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
ALPHA_MF, ALPHA_XY = 0.5, 2 * 0.6717 / (1 + 2 * 0.6717)


def load():
    D = {}
    for fn in glob.glob(os.path.join(OUT, "kz_k1_L*_t*_e*_s*.json")):
        d = json.load(open(fn))
        key = (d["L"], round(d["exponent"], 3))
        for s, v in d["seeds"].items():
            D.setdefault(key, {}).setdefault(d["tauQ"], {})[int(s)] = v
    return D


def alpha_fit(block, ab, taus, nboot=4000, seed=1):
    taus = [t for t in taus if t in block]
    arr = [np.array([block[t][s][ab] for s in sorted(block[t])]) for t in taus]
    x = np.log(taus)
    y = np.log([a.mean() for a in arr])
    a0 = -np.polyfit(x, y, 1)[0]
    rng = np.random.default_rng(seed)
    bs = np.empty(nboot)
    for i in range(nboot):
        yb = np.log([a[rng.integers(0, len(a), len(a))].mean() for a in arr])
        bs[i] = -np.polyfit(x, yb, 1)[0]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return {"alpha": float(a0), "ci95": [float(lo), float(hi)], "boot_sd": float(bs.std()), "tauQ": taus,
            "n_seeds": [len(a) for a in arr], "mean_density": [float(a.mean()) for a in arr],
            "sem_density": [float(a.std(ddof=1) / np.sqrt(len(a))) for a in arr]}


def main():
    D = load()
    res = {"alpha_mean_field": ALPHA_MF, "alpha_3DXY_modelA": ALPHA_XY}
    main_b = D.get((48, 0.573), {})
    allt = [64, 128, 256, 512, 1024, 2048, 4096]
    m = alpha_fit(main_b, "2_0.5", allt)
    res["main_L48_a2_b0.5"] = m
    passed = (m["ci95"][0] > ALPHA_MF or m["ci95"][1] < ALPHA_MF) and m["ci95"][0] <= ALPHA_XY <= m["ci95"][1]
    res["PASS_prediction_beyond_mean_field"] = bool(passed)
    sens = {}
    for a in ("1", "2", "4"):
        for b in ("0.25", "0.5", "1"):
            sens["a%s_b%s" % (a, b)] = alpha_fit(main_b, "%s_%s" % (a, b), allt)
    res["sensitivity_a_b_L48"] = {k: [v["alpha"], v["ci95"]] for k, v in sens.items()}
    res["sensitivity_max_dev"] = float(max(abs(v["alpha"] - m["alpha"]) for v in sens.values()))
    res["window_tau>=256_L48"] = alpha_fit(main_b, "2_0.5", [256, 512, 1024, 2048, 4096])
    if (32, 0.573) in D:
        res["L32"] = alpha_fit(D[(32, 0.573)], "2_0.5", allt)
        res["L32_window_tau>=256"] = alpha_fit(D[(32, 0.573)], "2_0.5", [256, 512, 1024, 2048, 4096])
    if (64, 0.573) in D:
        w = [256, 512, 1024, 2048, 4096]
        res["L64_window_tau>=256"] = alpha_fit(D[(64, 0.573)], "2_0.5", w)
        res["L64_minus_L48_same_window"] = res["L64_window_tau>=256"]["alpha"] - res["window_tau>=256_L48"]["alpha"]
    if (48, 0.5) in D:
        res["meanfield_protocol_L48"] = alpha_fit(D[(48, 0.5)], "2_0.5", allt)
    robust = passed and res["sensitivity_max_dev"] <= 0.05
    if "L64_minus_L48_same_window" in res:
        robust = robust and abs(res["L64_minus_L48_same_window"]) <= 0.05
    if "meanfield_protocol_L48" in res:
        robust = robust and res["meanfield_protocol_L48"]["ci95"][0] > ALPHA_MF
    res["ROBUST"] = bool(robust)
    # local slopes (main)
    md = np.array(m["mean_density"]); tt = np.array(m["tauQ"], float)
    res["local_slopes_main"] = (-np.diff(np.log(md)) / np.diff(np.log(tt))).tolist()
    # mean line spacing in lattice units at the slowest quench (finite-size indicator)
    res["line_spacing_slowest_L48"] = float(md[-1] ** -0.5)
    json.dump(res, open(os.path.join(OUT, "kz_k1.json"), "w"), indent=1)
    print("main alpha = %.3f  CI [%.3f, %.3f]  PASS=%s ROBUST=%s" % (m["alpha"], *m["ci95"], passed, robust))
    for k, v in res["sensitivity_a_b_L48"].items():
        print("  %-10s %.3f [%.3f, %.3f]" % (k, v[0], *v[1]))
    for k in ("window_tau>=256_L48", "L32", "L32_window_tau>=256", "L64_window_tau>=256", "meanfield_protocol_L48"):
        if k in res:
            print("  %-22s %.3f [%.3f, %.3f]" % (k, res[k]["alpha"], *res[k]["ci95"]))
    print("  local slopes", np.round(res["local_slopes_main"], 3), "spacing", res["line_spacing_slowest_L48"])


if __name__ == "__main__":
    main()

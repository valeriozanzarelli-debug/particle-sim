"""assist_val.py - validation exams V1 (Schwinger), V2 (Keldysh crossover), V3 (multiphoton) of the assist_ study.
Pre-registered in lab/out/assist_prereg.json. Writes lab/out/assist_val.json.   usage: python assist_val.py [V1|V2|V3]
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assist_core import Field, walk_pairs, dirac_pairs, sauter_pulse_exact, resonance, floquet_split_walk, \
    floquet_split_dirac, qbar, OUT, PI

PATH = os.path.join(OUT, "assist_val.json")
res = json.load(open(PATH)) if os.path.exists(PATH) else {}
pre = json.load(open(os.path.join(OUT, "assist_prereg.json")))
t0 = time.time()
part = sys.argv[1] if len(sys.argv) > 1 else "all"


def rich(f1, f2):          # f1 at theta, f2 at theta/2
    return (4 * f2 - f1) / 3


if part in ("V1", "all"):
    out = {}
    ok_a, ok_b = True, True
    for c in pre["V1_schwinger"]["configs"]:
        E, gs = c["E"], c["gs"]
        F = Field(E, gs * E)
        p = np.arange(-16, 16.0001, 0.02)
        fe = sauter_pulse_exact(p, E, 1 / (gs * E))
        r = dict(E=E)
        fl = {}
        for th in pre["V1_schwinger"]["thetas"]:
            fl[th] = walk_pairs(F, p, th)
            r[f"N_th{th}"] = float(np.trapezoid(fl[th], p) / (2 * PI))
        th1, th2 = pre["V1_schwinger"]["thetas"][-2:]
        mask = fe > 1e-3 * fe[np.argmin(abs(p))]
        fr = rich(fl[th1], fl[th2])
        r["maxdev_finest"] = float(np.max(np.abs(fl[th2][mask] / fe[mask] - 1)))
        r["maxdev_rich"] = float(np.max(np.abs(fr[mask] / fe[mask] - 1)))
        r["f0_ratio_by_theta"] = {str(th): float(fl[th][np.argmin(abs(p))] / fe[np.argmin(abs(p))]) for th in fl}
        r["N_exact"] = c["N_exact"]; r["N_LCFA"] = c["N_LCFA"]
        r["N_lat_over_LCFA"] = r[f"N_th{th2}"] / c["N_LCFA"]
        r["N_lat_over_exact"] = r[f"N_th{th2}"] / c["N_exact"]
        r["f0_exact_over_schwinger"] = c["f0_exact"] / c["f0_schwinger"]
        r["V1a"] = bool(r["maxdev_finest"] <= 0.03 and r["maxdev_rich"] <= 0.01)
        r["V1b"] = bool(abs(r["N_lat_over_LCFA"] - 1) <= 0.10)
        ok_a &= r["V1a"]; ok_b &= r["V1b"]
        r["curve_p"] = p[::20].tolist(); r["curve_lat"] = fl[th2][::20].tolist(); r["curve_exact"] = fe[::20].tolist()
        out[f"E{E}"] = r
        print("V1", E, {k: v for k, v in r.items() if not k.startswith("curve")}, f"{time.time()-t0:.0f}s", flush=True)
    out["V1a_PASS"] = ok_a; out["V1b_PASS"] = ok_b
    res["V1"] = out
    json.dump(res, open(PATH, "w"), indent=1)

if part in ("V2", "all"):
    out = {}
    ok = True
    for c in pre["V2_keldysh"]["configs"]:
        g, E = c["gamma"], c["E"]
        tau = 1 / (g * E)
        F = Field(E, g * E)
        P = 1 / g + 3
        p = np.arange(-P, P + 1e-9, 0.01)
        fe = sauter_pulse_exact(p, E, tau)
        Ne = float(np.trapezoid(fe, p) / (2 * PI))
        i0 = np.argmin(abs(p))
        r = dict(gamma=g, E=E, f0_exact=float(fe[i0]), N_exact=Ne, f0_schwinger=float(np.exp(-PI / E)))
        fl = {}
        for th in pre["V2_keldysh"]["thetas"]:
            fl[th] = walk_pairs(F, p, th)
            r[f"f0_ratio_th{th}"] = float(fl[th][i0] / fe[i0])
            r[f"N_ratio_th{th}"] = float(np.trapezoid(fl[th], p) / (2 * PI) / Ne)
        r["f0_ratio_rich"] = (4 * r["f0_ratio_th0.025"] - r["f0_ratio_th0.05"]) / 3
        r["N_ratio_rich"] = (4 * r["N_ratio_th0.025"] - r["N_ratio_th0.05"]) / 3
        fd = dirac_pairs(F, [0.0], 0.01)[0]
        r["dirac_over_exact_f0"] = float(fd / fe[i0])
        r["order_check_f0"] = float((r["f0_ratio_th0.1"] - 1) / (r["f0_ratio_th0.05"] - 1))
        r["PASS"] = bool(abs(r["f0_ratio_th0.025"] - 1) <= 0.03 and abs(r["N_ratio_th0.025"] - 1) <= 0.03
                         and abs(r["f0_ratio_rich"] - 1) <= 0.01 and abs(r["N_ratio_rich"] - 1) <= 0.01)
        ok &= r["PASS"]
        r["curve_p"] = p[::10].tolist(); r["curve_lat"] = fl[0.025][::10].tolist(); r["curve_exact"] = fe[::10].tolist()
        out[f"g{g}"] = r
        print("V2", {k: v for k, v in r.items() if not k.startswith("curve")}, f"{time.time()-t0:.0f}s", flush=True)
    out["V2_PASS"] = ok
    res["V2"] = out
    json.dump(res, open(PATH, "w"), indent=1)

if part in ("V3", "all"):
    out = {}
    for xi in pre["V3_multiphoton"]["xis"]:
        for n in range(1, 6):
            r = dict(n=n, xi=xi)
            wd, Dd = resonance(n, xi, lambda w: floquet_split_dirac(w, xi, 3000))
            r["omega_dirac"], r["Delta_dirac"] = wd, Dd
            for a in (0.05, 0.025):
                N = 2 * int(round(PI / (wd * a)))
                wl, Dl = resonance(n, xi, lambda w, N=N: floquet_split_walk(w, xi, N))
                r[f"omega_walk_a{a}"], r[f"Delta_walk_a{a}"], r[f"N_a{a}"] = wl, Dl, N
            r["omega_qbar"] = 2 * qbar(xi) / n
            out[f"xi{xi}_n{n}"] = r
            print("V3", r, f"{time.time()-t0:.0f}s", flush=True)
    ex = {}
    for xi in pre["V3_multiphoton"]["xis"]:
        for n in (2, 4):
            ex[f"V3a_xi{xi}_n{n}"] = out[f"xi{xi}_n{n}"]["Delta_walk_a0.025"] / out[f"xi{xi}_n{n+1}"]["Delta_walk_a0.025"]
    for n in (1, 3, 5):
        ex[f"V3b_slope_n{n}"] = float(np.log2(out[f"xi0.2_n{n}"]["Delta_walk_a0.025"] / out[f"xi0.1_n{n}"]["Delta_walk_a0.025"]))
        ex[f"V3b_slope_dirac_n{n}"] = float(np.log2(out[f"xi0.2_n{n}"]["Delta_dirac"] / out[f"xi0.1_n{n}"]["Delta_dirac"]))
    ex["V3c_maxdev"] = max(abs(out[k]["omega_walk_a0.025"] / out[k]["omega_dirac"] - 1) for k in out)
    ex["qbar_maxdev"] = max(abs(out[k]["omega_qbar"] / out[k]["omega_dirac"] - 1) for k in out)
    ex["V3a_PASS"] = bool(all(v <= 1e-2 for k, v in ex.items() if k.startswith("V3a_")))
    ex["V3b_PASS"] = bool(all(abs(ex[f"V3b_slope_n{n}"] / n - 1) <= 0.05 for n in (1, 3, 5)))
    ex["V3c_PASS"] = bool(ex["V3c_maxdev"] <= 1e-3)
    out["exams"] = ex
    print("V3 exams", ex)
    res["V3"] = out
    json.dump(res, open(PATH, "w"), indent=1)
print("done", f"{time.time()-t0:.0f}s")

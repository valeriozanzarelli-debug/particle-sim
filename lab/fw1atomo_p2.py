"""fw1atomo_p2.py - P2 of the pre-registration: the BATTIMENTO_1S test redone with FW1 next to the frozen rule.
48^3 cells, a = 0,5 lambda_C, Z = 60, bare point nucleus on the corner of 8 cells; start = 1s of the radial form of
the same rule put on the cube (as src/atomo/zz.rs::con_stato, the Laboratorio's radial 1s).
Steps per rule R in {F0, FW1}:  ritz_R  (levels in [0,82; 0,94], 600 periods)
                                surv_R  (survival from the radial 1s, 1600 periods = 3200 ticks)
                                fix_R   (survival from the main exact level, 1600 periods: tool check F3)
Run: OMP_NUM_THREADS=1 timeout 590 python lab/fw1atomo_p2.py <step>  (exit 3 = checkpoint); sum -> lab/out/fw1atomo_p2.json"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import fw1atomo_core as fc

T0 = time.time()
N, A = 48, 0.5
OUT = os.path.join(os.path.dirname(__file__), "out", "fw1atomo_p2.json")
RES = os.path.join(fc.SCR, "p2")
os.makedirs(RES, exist_ok=True)
REF = [0.8993, 0.8996, 0.8863, 0.8404]
H = float(os.environ.get("FW1_H", "0.4"))
TAG = os.environ.get("FW1_TAG", "")


def setup(rule):
    L = fc.Lattice(N, A, "angolo")
    R = L.rule(rule, H)
    psi0, om_rad = L.radial1s()
    return L, R, psi0, om_rad


def surv(R, psi0, nper, ck, tmax=540):
    st = json.load(open(ck + ".json")) if os.path.exists(ck + ".json") else {"p": 0}
    if st["p"]:
        psi = np.load(ck + "_psi.npy")
        c = np.load(ck + "_c.npy")
    else:
        psi, c = psi0.copy(), np.zeros(nper + 1, complex)
        c[0] = 1
    while st["p"] < nper:
        if time.time() - T0 > tmax:
            np.save(ck + "_psi.npy", psi)
            np.save(ck + "_c.npy", c)
            json.dump(st, open(ck + ".json", "w"))
            print("checkpoint", st, flush=True)
            raise SystemExit(3)
        psi = R.period(psi)
        st["p"] += 1
        c[st["p"]] = np.vdot(psi0, psi)
    for f in (".json", "_psi.npy", "_c.npy"):
        if os.path.exists(ck + f):
            os.remove(ck + f)
    return c


def run(step):
    kind, rule = step.split("_")
    fn = os.path.join(RES, step + TAG + ".json")
    if kind == "ritzc":
        fn = os.path.join(RES, f"ritzc{os.environ['FW1_CLS']}_{rule}{TAG}.json")
    if os.path.exists(fn):
        print("done already", step)
        return
    L, R, psi0, om_rad = setup(rule)
    if kind == "ritzc":                       # amendment A1: one period-sublattice at a time
        c = int(os.environ["FW1_CLS"])
        S = psi0 * (L.cls == c)[None]
        S = S / np.linalg.norm(S)
        lv = fc.ritz(R, [S], A, 0.82, 0.94, 0.005, 600, ck=os.path.join(fc.SCR, f"p2c{c}_" + rule + TAG), tmax=540, t0=T0)
        out = []
        for o in lv:
            wfull = float(abs(np.vdot(psi0, o["v"]))**2)
            if o["resid"] < 1e-4 and wfull >= 0.002:
                d = L.diag(o["v"])
                d.update(omega=o["omega"], resid=o["resid"], w=wfull, cls=c)
                out.append(d)
        if c == 0:
            main = max(lv, key=lambda o: abs(np.vdot(psi0, o["v"]))**2 if o["resid"] < 1e-6 else -1)
            np.save(os.path.join(fc.SCR, f"p2_main_{rule}{TAG}.npy"), main["v"])
        fn = os.path.join(RES, f"ritzc{c}_{rule}{TAG}.json")
        json.dump({"levels": out}, open(fn, "w"), indent=1)
        for d in out:
            print(f"c{c} {d['omega']:.5f} res {d['resid']:.1e} w {d['w']:.4f} G {d['val_Gamma']:.3f} r {d['r_mean']:.2f} S {np.round(d['spin'], 2)}", flush=True)
    elif kind == "ritz":
        lv = fc.ritz(R, [psi0], A, 0.82, 0.94, 0.005, 600, ck=os.path.join(fc.SCR, "p2r_" + rule + TAG), tmax=540, t0=T0)
        out = []
        for o in lv:
            if o["resid"] < 1e-4 and o["w_start"][0] >= 0.002:
                d = L.diag(o["v"])
                d.update(omega=o["omega"], resid=o["resid"], w=o["w_start"][0])
                out.append(d)
        main = max(lv, key=lambda o: o["w_start"][0] if o["resid"] < 1e-6 else -1)
        np.save(os.path.join(fc.SCR, f"p2_main_{rule}{TAG}.npy"), main["v"])
        json.dump({"levels": out, "omega_radial_rule": om_rad}, open(fn, "w"), indent=1)
        for d in out:
            print(f"{d['omega']:.5f} res {d['resid']:.1e} w {d['w']:.4f} G {d['val_Gamma']:.3f} N {d['val_N']:.3f} pi {d['val_pi']:.4f} "
                  f"r {d['r_mean']:.2f} S {np.round(d['spin'], 3)}", flush=True)
    elif kind == "surv":
        c = surv(R, psi0, 1600, os.path.join(fc.SCR, "p2s_" + rule + TAG))
        json.dump({"c_re": c.real.tolist(), "c_im": c.imag.tolist()}, open(fn, "w"))
    elif kind == "fix":
        v = np.load(os.path.join(fc.SCR, f"p2_main_{rule}{TAG}.npy"))
        c = surv(R, v, 1600, os.path.join(fc.SCR, "p2f_" + rule + TAG))
        json.dump({"max_p_exc": float(np.max(1 - np.abs(c)**2))}, open(fn, "w"))
    print("DONE", step, time.time() - T0)


def beat_stats(c, T=2 * A):
    pex = 1 - np.abs(c)**2
    t = np.arange(len(c)) * T
    m = (t >= 450) & (t <= 1600)
    win = [(float(t0), float(pex[(t >= t0) & (t < t0 + 200)].max() - pex[(t >= t0) & (t < t0 + 200)].min()))
           for t0 in range(0, 1600, 200)]
    return {"max_p_exc": float(pex.max()), "range_450_1600": float(pex[m].max() - pex[m].min()), "range_per_200": win,
            "lines": sorted(fc.lines(c, T), key=lambda x: -x[1])[:6], "p_exc_every_25": pex[::25].tolist()}


def summary():
    res = {"ref_battimento": REF}
    for rule in ("F0", "FW1"):
        r = {}
        for kind in ("ritz", "surv", "fix"):
            fn = os.path.join(RES, f"{kind}_{rule}{TAG}.json")
            if os.path.exists(fn):
                r[kind] = json.load(open(fn))
        if "surv" in r:
            d = r.pop("surv")
            r["surv"] = beat_stats(np.array(d["c_re"]) + 1j * np.array(d["c_im"]))
        res[rule] = r
    ex = {}
    if "ritz" in res["F0"]:
        lv = [l for l in res["F0"]["ritz"]["levels"] if l["resid"] < 1e-6 and l["w"] >= 0.01]
        dev = [min(abs(l["omega"] - e) for l in lv) if lv else 1 for e in REF]
        ex["V1"] = {"dev": dev, "PASS": bool(max(dev) <= 5e-4)}
    if "surv" in res["F0"]:
        ex["V2"] = {"max_p_exc": res["F0"]["surv"]["max_p_exc"], "PASS": res["F0"]["surv"]["max_p_exc"] >= 0.8}
    import glob
    cfiles = sorted(glob.glob(os.path.join(RES, f"ritzc?_FW1{TAG}.json")))
    if cfiles:
        res["FW1"]["ritzc"] = [l for f in cfiles for l in json.load(open(f))["levels"]]
        res["FW1"]["ritz"] = {"levels": res["FW1"]["ritzc"]}
    if "ritz" in res["FW1"]:
        lv = sorted([l for l in res["FW1"]["ritz"]["levels"] if l["resid"] < 1e-6 and l["w"] >= 0.02], key=lambda l: l["omega"])
        groups = []
        for l in lv:
            if groups and l["omega"] - groups[-1][-1]["omega"] <= 2e-3:
                groups[-1].append(l)
            else:
                groups.append([l])
        ex["F1"] = {"groups": [[round(l["omega"], 5) for l in g] for g in groups],
                    "weights": [round(sum(l["w"] for l in g), 4) for g in groups], "PASS": len(groups) == 1}
    if "surv" in res["FW1"]:
        ex["F2"] = {"range_450_1600": res["FW1"]["surv"]["range_450_1600"], "PASS": res["FW1"]["surv"]["range_450_1600"] < 0.05}
    if "fix" in res["FW1"]:
        ex["F3"] = {"max_p_exc": res["FW1"]["fix"]["max_p_exc"], "PASS": res["FW1"]["fix"]["max_p_exc"] < 1e-6}
    res["esami"] = ex
    out = OUT if not TAG else OUT.replace(".json", TAG + ".json")
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(ex, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "sum":
        summary()
    else:
        run(sys.argv[1])

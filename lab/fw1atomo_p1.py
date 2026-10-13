"""fw1atomo_p1.py - P1 of the pre-registration: the two lines of the FW1 1s at theta = 0,5 (Z = 60).
Each configuration: Ritz (filter diagonalisation) from the given start states; for each level found (residual
< 1e-4, weight >= 0,5 % in some start) the diagnostics of fw1atomo_core.Lattice.diag.
Run:  OMP_NUM_THREADS=1 timeout 590 python lab/fw1atomo_p1.py <config>   (exit 3 = checkpoint: run again)
      python lab/fw1atomo_p1.py list | sum      (sum -> lab/out/fw1atomo_p1.json)"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import fw1atomo_core as fc

T0 = time.time()
OUT = os.path.join(os.path.dirname(__file__), "out", "fw1atomo_p1.json")
RES = os.path.join(fc.SCR, "p1")
os.makedirs(RES, exist_ok=True)


def cfg(rule, h=0.4, n=40, a=0.5, pos="angolo", starts="ud", lo=0.80, hi=0.97, dw=0.005, tlen=1600.0):
    return dict(rule=rule, h=h, n=n, a=a, pos=pos, starts=starts, lo=lo, hi=hi, dw=dw, nper=int(round(tlen / (2 * a))))


CONFIGS = {}
# P1a base
for r in ("FW1", "F0"):
    CONFIGS[f"a_{r}"] = cfg(r)
# P1b one period-sublattice at a time (spin up)
for r in ("FW1", "F0"):
    for c in range(8):
        CONFIGS[f"b_{r}_c{c}"] = cfg(r, starts=f"class{c}")
# P1c h scan (amendment A4: one period-sublattice at a time, 0 and 1)
for h in (0.2, 0.3, 0.6, 0.8):
    for c in (0, 1):
        CONFIGS[f"c_h{h}_c{c}"] = cfg("FW1", h=h, starts=f"class{c}")
# P1d theta scan (box ~20 lambda_C, ~1200 lambda_C/c; A2: dw 0,01; A4: sublattices 0 and 1)
for a, n in ((0.4, 50), (1 / 3, 60), (0.25, 80)):
    for r in ("FW1", "F0"):
        for c in ((0, 1) if r == "FW1" else (0,)):
            CONFIGS[f"d_{r}_a{a:.3f}_c{c}"] = cfg(r, n=n, a=a, starts=f"class{c}", tlen=1200.0, dw=0.01)
# P1e nucleus position (A4: sublattices 0 and 1)
for pos in ("spigolo", "faccia", "sito"):
    for r in ("FW1", "F0"):
        for c in (0, 1):
            CONFIGS[f"e_{r}_{pos}_c{c}"] = cfg(r, pos=pos, starts=f"class{c}")


def starts_of(L, spec):
    if spec == "ud":
        return [L.dirac1s((1, 0)), L.dirac1s((0, 1))]
    if spec == "u":
        return [L.dirac1s((1, 0))]
    if spec.startswith("class"):
        c = int(spec[5:])
        s = L.dirac1s((1, 0)) * (L.cls == c)[None]
        return [s / np.linalg.norm(s)]
    raise ValueError(spec)


def run(name):
    fn = os.path.join(RES, name + ".json")
    if os.path.exists(fn):
        print("already done", name)
        return
    C = CONFIGS[name]
    L = fc.Lattice(C["n"], C["a"], C["pos"])
    R = L.rule(C["rule"], C["h"])
    S = starts_of(L, C["starts"])
    tp = time.time()
    R.period(S[0])
    t_per = time.time() - tp
    lv = fc.ritz(R, S, C["a"], C["lo"], C["hi"], C["dw"], C["nper"], ck=os.path.join(fc.SCR, "p1_" + name), tmax=540, t0=T0)
    out = []
    for o in lv:
        if o["resid"] < 1e-4 and max(o["w_start"]) >= 0.005:
            d = L.diag(o["v"])
            d.update(omega=o["omega"], resid=o["resid"], w_start=o["w_start"])
            out.append(d)
    json.dump({"config": C, "levels": out, "s_per_period": t_per}, open(fn, "w"), indent=1)
    for d in out:
        print(f"{d['omega']:.5f} res {d['resid']:.1e} w {np.round(d['w_start'], 3)} G {d['val_Gamma']:.3f} N {d['val_N']:.4f} "
              f"pi {d['val_pi']:.4f} r {d['r_mean']:.2f} in1 {d['w_r_lt_1']:.3f} S {np.round(d['spin'], 3)} "
              f"col {np.round(d['colours'], 3)} cls {np.round(d['classes'], 3)}", flush=True)
    print("DONE", name, time.time() - T0)


def summary():
    res = {}
    for name in CONFIGS:
        fn = os.path.join(RES, name + ".json")
        if os.path.exists(fn):
            res[name] = json.load(open(fn))
    json.dump(res, open(OUT, "w"), indent=1)
    print("written", OUT, len(res), "configs")


if __name__ == "__main__":
    if sys.argv[1] == "list":
        print(" ".join(CONFIGS))
    elif sys.argv[1] == "sum":
        summary()
    else:
        try:
            run(sys.argv[1])
        except SystemExit as e:
            sys.exit(e.code)

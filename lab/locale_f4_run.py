"""locale_f4_run.py - phase 4 (1): local terms for the pi-valleys, exams Q1-Q4 of lab/out/locale_f4_prereg.json.
Writes lab/out/locale_f4_run.json (checkpointed, exit 3 = run again)."""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_f4_core import walk, measure, E_low

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f4_run.json")
T0 = time.time()
R = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    with open(OUT, "w") as f:
        json.dump(R, f, indent=1)


def verdict(r):
    q1 = r["min_doubler_gap"] >= 0.3
    q2 = r["E0_rel"] <= 1e-12 and r["drift"] <= 1e-9 and r["unitarity"] <= 1e-12
    q3 = max(abs(x) for x in r["mkin_vs_frozen"]) <= 1e-4
    s = r["dE_vs_frozen"]["slope"]
    q4 = (s is not None and s >= 3.5) and r["dE_vs_frozen"]["rel_at_p_eq_m"] <= 1e-4
    return dict(Q1=q1, Q2=q2, Q3=q3, Q4=q4, FOUND=q1 and q2 and q3 and q4)


def mk_x(th, ref, **kw):
    w = walk(th, **kw)
    ps = np.linspace(0.02, 0.1, 9) * th
    d = np.array([1, 0, 0.])
    E = E_low(w, ps[:, None] * d[None]); E0 = E_low(w, np.zeros((1, 3)))[0]
    c = np.linalg.lstsq(np.stack([ps ** 2, ps ** 4], 1), E ** 2 - E0 ** 2, rcond=None)[0]
    Er = E_low(ref, ps[:, None] * d[None]); E0r = E_low(ref, np.zeros((1, 3)))[0]
    cr = np.linalg.lstsq(np.stack([ps ** 2, ps ** 4], 1), Er ** 2 - E0r ** 2, rcond=None)[0]
    return cr[0] / c[0] - 1


CANDS = []
for m in (0.02, 0.05, 0.1, 0.3):
    CANDS.append(("C1", dict(mu_x=m)))
for m in (0.05, 0.1):
    CANDS.append(("C2", dict(mu_y=m)))
for h in (0.02, 0.05, 0.1):
    CANDS.append(("C3", dict(hs=h)))
for m in (0.1, 0.3):
    CANDS.append(("C4", dict(mu_x=m, hu=0.1)))
for m in (0.1, 0.3):
    for hu in (0.0, 0.1):
        CANDS.append(("C5", dict(mu_x=m, hu=hu, hs="tune")))

for th in (0.1,):
    ref = walk(th)
    for fam, kw in CANDS:
        key = f"{fam}|{th}|" + ",".join(f"{k}={v}" for k, v in kw.items())
        if key in R:
            continue
        if time.time() - T0 > 480:
            save(); print("rerun"); sys.exit(3)
        kw = dict(kw)
        extra = {}
        if kw.get("hs") == "tune":
            f = lambda h: mk_x(th, ref, mu_x=kw["mu_x"], hu=kw["hu"], hs=h)
            grid = np.linspace(-0.5, 0.5, 41)
            vals = [f(h) for h in grid]
            roots = []
            for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
                if np.sign(fa) != np.sign(fb):
                    roots.append(brentq(f, a, b, xtol=1e-12))
            extra["roots"] = roots
            if not roots:
                R[key] = dict(no_root=True, scan=[[float(g), float(v)] for g, v in zip(grid, vals)])
                save(); print(key, "no root", flush=True); continue
            kw["hs"] = min(roots, key=abs)
        r = measure(walk(th, **kw), th, ref)
        r.update(extra); r["params"] = kw; r["verdict"] = verdict(r)
        R[key] = r
        print(key, f"gap {r['min_doubler_gap']:.3f} (N {r['N_xy']:.3f} X {r['X_x']:.3f} R {r['R_xyz']:.3f} zpi {r['zone_gpi']:.3f}) "
              f"mk {max(abs(x) for x in r['mkin_vs_frozen']):.1e} slope {r['dE_vs_frozen']['slope']} "
              f"dE(m) {r['dE_vs_frozen']['rel_at_p_eq_m']:.1e} E0 {r['E0_rel']:.0e} -> {r['verdict']}", flush=True)
        save()

# found ones at theta = 0.02 (information)
for key in list(R):
    if key.startswith(("C", )) and "|0.1|" in key and R[key].get("verdict", {}).get("FOUND"):
        k2 = key.replace("|0.1|", "|0.02|")
        if k2 in R:
            continue
        ref = walk(0.02)
        r = measure(walk(0.02, **R[key]["params"]), 0.02, ref)
        r["params"] = R[key]["params"]; r["verdict"] = verdict(r); R[k2] = r; save()
save()
print("ALL DONE", time.time() - T0)

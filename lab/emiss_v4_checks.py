"""emiss_v4_checks.py - informative tests D4 and V4 of lab/out/emiss_prereg2.json (new states, learned nucleus).
Run: python lab/emiss_v4_checks.py -> lab/out/emiss_v4_checks.json"""
import json, os, time
import numpy as np
from emiss_core import OUT, ovl
os.environ["EMISS_CFG"] = "none"
import emiss_configs as EC

EX = json.load(open(os.path.join(OUT, "emiss_exact.json")))
R = EX["Z60"]["2p3/2"]["R_length"]
R12 = EX["Z60"]["2p1/2"]["R_length"]
cont = {"2p1/2+": R12**2 / 9, "2p3/2+": 2 * R**2 / 9, "2p3/2+3": R**2 / 3, "2p1/2-": 2 * R12**2 / 9, "2p3/2-": R**2 / 9}
D = EC.atom4()
s = EC.level4("1s+")
res = {"D4": {}}
for st, c in cont.items():
    d = D.dip_matrix(s, EC.level4(st))
    res["D4"][st] = dict(abs=np.abs(d).tolist(), sum2=float(np.sum(np.abs(d)**2)), continuum=c, ratio=float(np.sum(np.abs(d)**2) / c))
    print(f"D4 {st}: |d| = {np.round(np.abs(d), 4)}, sum |d|^2 = {np.sum(np.abs(d)**2):.4f} vs Dirac {c:.4f} (ratio {np.sum(np.abs(d)**2) / c:.3f})", flush=True)
# V4: free evolution 3000 time units
t0 = time.time()
ps = {k: EC.level4(k).copy() for k in ("1s+", "2p3/2+3")}
ref = {k: EC.level4(k) for k in ps}
drift = {k: [] for k in ps}
nt = int(3000 / D.a)
for t in range(nt + 1):
    if t % 40 == 0:
        for k in ps:
            drift[k].append(1 - abs(ovl(ref[k], ps[k]))**2)
    if t < nt:
        for k in ps:
            D.tick(ps[k], t)
res["V4"] = {k: dict(max=float(np.max(v)), end=float(v[-1]), series=[float(x) for x in v[::5]]) for k, v in drift.items()}
print("V4:", {k: (round(v["max"], 5), round(v["end"], 5)) for k, v in res["V4"].items()}, f"({time.time() - t0:.0f} s)")
json.dump(res, open(os.path.join(OUT, "emiss_v4_checks.json"), "w"), indent=1)

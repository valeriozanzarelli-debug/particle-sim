"""kz_k1.py - K1 production chunk: 3D XY linear quench, vortex-line density at a*t_hat after cooling b*t_hat.
usage: python kz_k1.py L tauQ seed0 nseeds [exponent]     -> lab/out/kz_k1_L{L}_t{tauQ}_e{exp}_s{seed0}.json
"""
import sys, os, json, time, math
import torch
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from kz_core import quench

L = int(sys.argv[1]); tauQ = int(sys.argv[2]); s0 = int(sys.argv[3]); ns = int(sys.argv[4])
ex = float(sys.argv[5]) if len(sys.argv) > 5 else 0.573
A = (1, 2, 4); B = (0.25, 0.5, 1.0)
that = tauQ ** ex
ct = {a * that: [b * that for b in B] for a in A}
out = {"L": L, "tauQ": tauQ, "exponent": ex, "t_hat": that, "a": A, "b": B, "seeds": {}}
fn = os.path.join(HERE, "out", "kz_k1_L%d_t%d_e%.3f_s%d.json" % (L, tauQ, ex, s0))
for s in range(s0, s0 + ns):
    t0 = time.time()
    keep = (L == 48 and tauQ == 1024 and s == 0 and abs(ex - 0.573) < 1e-6)
    r = quench(L, 3, tauQ, list(ct), 1000 * L + 7 * tauQ + s, cool_times=ct, keep_final=keep)
    if keep:
        r, th = r
        torch.save(th, os.path.join(HERE, "out", "kz_final_L48_t1024.pt"))
    out["seeds"][s] = {"%g_%g" % (round(k[0] / that, 3), round(k[1] / that, 3)): v[0] / L ** 3 for k, v in r.items()}
    print(L, tauQ, s, "%.1fs" % (time.time() - t0), out["seeds"][s]["2_0.5"], flush=True)
    with open(fn, "w") as f:
        json.dump(out, f, indent=1)

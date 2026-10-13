"""nucleo1s_f8_posthoc.py - AFTER the verdicts (declared, not judged): does the control's error go to zero when the
zoom chain covers more of the atom?  Exact band, J0 = 32, 48, 64 (R_f = 4.1, 6.1, 8.1 a0/Z), 1s and 2s.
-> lab/out/nucleo1s_f8_posthoc.json (re-run until 'all done')"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_f8_core import CoarseT, Glue8, delta_pair, delta_dirac, A0  # noqa: E402
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
F = os.path.join(OUT, "nucleo1s_f8_posthoc.json")
res = json.load(open(F)) if os.path.exists(F) else {}
T0 = time.time()
cs = {B: CoarseT(B, "exact") for B in (256, 512)}
for J0 in (32, 48, 64):
    for Z, B in ((1, 512), (2, 256)):
        for n in (1, 2):
            k = f"J{J0}_Z{Z}_n{n}"
            if k in res:
                continue
            if time.time() - T0 > 480:
                print("time budget: run again"); sys.exit(0)
            r = Glue8(cs[B], Z, J0, 90 * A0 / Z).solve(n)
            res[k] = r
            json.dump(res, open(F, "w"), indent=1, default=float)
            print(k, f"{r['rel']:+.3e}", flush=True)
for J0 in (32, 48, 64):
    for n in (1, 2):
        res[f"err_J{J0}_{n}s"] = delta_pair(res[f"J{J0}_Z1_n{n}"]["B"], res[f"J{J0}_Z2_n{n}"]["B"]) / delta_dirac(n) - 1
        print(J0, n, f"{res[f'err_J{J0}_{n}s']:+.4f}")
json.dump(res, open(F, "w"), indent=1, default=float)
print("all done")

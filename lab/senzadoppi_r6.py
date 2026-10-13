"""senzadoppi_r6.py - amendment A4: exam E6 redone with a sweep of k_x over ONE period of the copies (pi instead of
2 pi), so that every mode crosses each valley once (no Stueckelberg interference between two Bloch cycles of the same
mode).  Reference: 4 Gamma crossings (transverse positions {0, pi}^2) x one species.  Rules F0, PA, PAW, FW1.
Output lab/out/senzadoppi_r6.json.  Run: OMP_NUM_THREADS=2 timeout 590 python lab/senzadoppi_r6.py  (rerun until 0)"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import senzadoppi_r3 as r3
_ev3, _pm3 = r3.evolve, r3.period_matrix          # original each/once versions
import senzadoppi_r5 as r5            # registers FW1 and patches r3.evolve / r3.period_matrix for 'each_w0'

r3.SWEEP = np.pi
r3.KX0 = -0.75 * np.pi
r3.RAMP_K = 0.3
CK = os.path.join(os.path.dirname(r3.CK), "senzadoppi_r6_ck.json")
OUT = os.path.join(os.path.dirname(__file__), "out", "senzadoppi_r6.json")
T0 = time.time()
_ev5, _pm5 = r3.evolve, r3.period_matrix


def main():
    d = json.load(open(CK)) if os.path.exists(CK) else {}
    nlines = r3.NPERP**2
    for name in ("FW1", "PAW", "PA", "F0"):
        mod = r3
        r3.evolve, r3.period_matrix = (_ev5, _pm5) if name == "FW1" else (_ev3, _pm3)
        for fac in (1.0, 2.0):
            key = f"{name}_E{fac}"
            st = d.setdefault(key, {"sum": 0.0, "done": 0, "Pgrid": []})
            while st["done"] < nlines:
                if time.time() - T0 > 540:
                    json.dump(d, open(CK, "w"))
                    print("checkpoint", key, st["done"], flush=True)
                    sys.exit(3)
                idx = np.arange(st["done"], min(st["done"] + r3.CHUNK, nlines))
                P, Eeff = mod.run_chunk(name, fac * r3.THETA**2, idx)
                st["sum"] += float(np.sum(P))
                st["Pgrid"] += [float(x) for x in P]
                st["Eeff"] = Eeff
                st["done"] = int(idx[-1] + 1)
                print(key, st["done"], round(time.time() - T0), flush=True)
            json.dump(d, open(CK, "w"))
    dA = (2 * np.pi / r3.NPERP)**2
    res = {}
    for key, st in d.items():
        E = st["Eeff"]
        ref = 2 * E * np.exp(-np.pi * r3.THETA**2 / E)
        res[key] = {"E_over_Ec": E / r3.THETA**2, "R": st["sum"] * dA / (4 * ref),
                    "max_P_line": float(max(st["Pgrid"]))}
        print(key, res[key])
    out = {"sweep": "pi", "results": res,
           "E6_PASS_PAW": bool(all(abs(res[f"PAW_E{f}"]["R"] - 1) <= 0.1 for f in (1.0, 2.0))),
           "E6_PASS_FW1": bool(all(abs(res[f"FW1_E{f}"]["R"] - 1) <= 0.1 for f in (1.0, 2.0)))}
    json.dump(out, open(OUT, "w"), indent=1)
    print(out["E6_PASS_PAW"], out["E6_PASS_FW1"])


if __name__ == "__main__":
    main()

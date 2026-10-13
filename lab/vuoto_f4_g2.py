"""phase 4 production: g - 2 with the exact discrete-time vertex (lab/vuoto_f4_core.py), checkpointed.
Output lab/out/vuoto_f4_g2.json.  Exit 3 = not finished (driver reruns)."""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vuoto_f4_core as F4
from zz_core import ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "vuoto_f4_g2.json")
A2PI = ALPHA / (2 * np.pi)
CASES = [("R", "K12", (0.6, 1.05)), ("R", "K12/3", (0.2, 0.35)), ("C", "K12/3", (0.2, 0.35))]
THETAS = [0.05, 0.025, 0.0125, 0.00625]
KMINS = [1e-3, 3e-4]
PAIRS = [(2, 0), (2, 1), (0, 1), (0, 2), (1, 0), (1, 2)]
GRID = dict(nr=64, nct=12, nph=16)
T0 = time.time()


def load():
    if os.path.exists(OUT):
        return json.load(open(OUT))
    return {"grid": GRID, "pairs": PAIRS, "runs": {}}


def save(d):
    json.dump(d, open(OUT + ".tmp", "w"), indent=1)
    os.replace(OUT + ".tmp", OUT)


def main():
    d = load()
    for rule, cn, K12 in CASES:
        for th in THETAS:
            for kf in KMINS:
                key = "%s|%s|%g|%g" % (rule, cn, th, kf)
                run = d["runs"].setdefault(key, {})
                g = dict(GRID, K12=K12, kmin_fac=kf)
                if "dmv" not in run:
                    run["dmv"] = {k: v for k, v in F4.dmv(th, rule, g).items() if k != "rel_P"}
                    save(d)
                dmv = run["dmv"]
                for qa, ja in PAIRS:
                    pk = "%d%d" % (qa, ja)
                    if pk in run:
                        continue
                    if time.time() - T0 > 480:
                        save(d); print("checkpoint", key, pk); sys.exit(3)
                    o = F4.g2_pair(th, rule, qa, ja, g)
                    lin = (o["r"].real - 1) + dmv["dm_over_m"] - dmv["dv"][ja] - dmv["dv"][qa]
                    run[pk] = {"a_over_a2pi": lin / A2PI, "imag_r": o["r"].imag, "resid": o["resid"],
                               "parts_over_a2pi": {k: o[k].real / A2PI for k in ("Vz", "LamT", "LamC", "leg", "rel")}}
                    save(d)
                    print(key, pk, round(lin / A2PI, 5), flush=True)
    print("ALL DONE")


if __name__ == "__main__":
    main()

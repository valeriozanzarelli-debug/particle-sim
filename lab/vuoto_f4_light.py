"""phase 4 production: the dressed light (photon self-energy of the sea, lab/vuoto_f4_core.py), checkpointed.
Output lab/out/vuoto_f4_light.json.  Exit 3 = not finished."""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vuoto_f4_core as F4

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "vuoto_f4_light.json")
CASES = [("R", "K12", (0.6, 1.05)), ("R", "K12/3", (0.2, 0.35)), ("C", "K12/3", (0.2, 0.35))]
THETAS = [0.05, 0.025, 0.0125, 0.00625]
AXES = {"x": 0, "y": 1, "z": 2, "d": [1, 1, 1]}
GRID = dict(nr=64, nct=24, nph=32, wmode="occ_mean")
EXTRA = [("R", "K12", (0.6, 1.05), 0.05, "grid2", dict(nr=96, nct=32, nph=48, wmode="occ_mean"), (0.05, 0.1)),
         ("R", "K12", (0.6, 1.05), 0.05, "q2", GRID, (0.1, 0.2)),
         ("R", "K12", (0.6, 1.05), 0.05, "wsym", dict(GRID, wmode="sym"), (0.05, 0.1))]
T0 = time.time()


def clean(r):
    return {k: v for k, v in r.items() if k not in ("pol_vectors", "offdiag")} | {
        "pol_vectors": [[[x.real, x.imag] for x in v] for v in r["pol_vectors"]],
        "offdiag": [[x.real, x.imag] for x in r["offdiag"]]}


def main():
    d = json.load(open(OUT)) if os.path.exists(OUT) else {"grid": GRID, "runs": {}}
    todo = []
    for rule, cn, K12 in CASES:
        for th in THETAS:
            for an, ax in AXES.items():
                todo.append(("%s|%s|%g|%s" % (rule, cn, th, an), rule, K12, th, ax, GRID, (0.05, 0.1)))
    for rule, cn, K12, th, tag, g, qfs in EXTRA:
        for an in ("x", "y", "z"):
            todo.append(("%s|%s|%g|%s|%s" % (rule, cn, th, an, tag), rule, K12, th, AXES[an], g, qfs))
    for key, rule, K12, th, ax, g, qfs in todo:
        if key in d["runs"]:
            continue
        if time.time() - T0 > 420:
            print("checkpoint", key); sys.exit(3)
        r = F4.photon_dc(th, rule, ax, "species", qfs=qfs, K12=K12, **g)
        d["runs"][key] = clean(r)
        json.dump(d, open(OUT + ".tmp", "w"), indent=1); os.replace(OUT + ".tmp", OUT)
        print(key, np.round(r["dc_over_alpha"], 5), flush=True)
    print("ALL DONE")


if __name__ == "__main__":
    main()

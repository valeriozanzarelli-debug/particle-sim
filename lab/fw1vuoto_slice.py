"""fw1vuoto_slice.py - amendment A1 (lab/out/fw1vuoto_emendamenti.json): stroboscopic density at the OTHER time slice
(period started with tick 1), frozen + projector and FW1 whole sea, so that the slice average tests the x <-> z
symmetry (E3g as registered compared slice-0 numbers, where the frozen rule also has x != z). Output
lab/out/fw1vuoto_slice.json. Run: OMP_NUM_THREADS=1 timeout 590 python lab/fw1vuoto_slice.py (rerun until ALL DONE)"""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fw1vuoto_core as c

FN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "fw1vuoto_slice.json")
T0 = time.time()


def rotated(mv):
    i = [k for k, m in enumerate(mv) if m["phi"]][0] + 1      # end of tick 0
    return mv[i:] + mv[:i]


def main():
    D = json.load(open(FN)) if os.path.exists(FN) else {}
    for rule, mode in (("F0", "ball"), ("FW1", "cellG")):
        for th in (0.025, 0.0125):
            beta = (5.0 if th > 0.02 else 6.0) if mode == "ball" else (5.0 if th > 0.02 else 6.5)
            mv = rotated(c.frozen_moves(th) if rule == "F0" else c.fw1_moves(th))
            for d in "xyz":
                for qm in (0.025, 0.05):
                    k = f"{rule}|{mode}|{th}|{d}|{qm}"
                    if k in D:
                        continue
                    if time.time() - T0 > 540:
                        sys.exit(3)
                    q = qm * th * np.eye(3)[("xyz").index(d)]
                    Ks, WW = c.quad_cell(mode, 64, beta, q)
                    w = c.species_weight() if mode == "ball" else None
                    v = c.integrate(lambda K: c.strobo_F(mv, K, q, w), Ks, WW)
                    D[k] = float(v.real) / (qm * th) ** 2
                    json.dump(D, open(FN, "w"), indent=0)
                    print(k, D[k], flush=True)
    print("ALL DONE")


if __name__ == "__main__":
    main()

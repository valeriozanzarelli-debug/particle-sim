"""vuoto_run3d.py - production (pre-registered, lab/out/vuoto_prereg.json): 3-D static polarization of the zig-zag sea.
Checkpointed: every evaluation is stored in lab/out/vuoto_3d.json; rerun until it prints ALL DONE (exit 3 = rerun me).
Run: timeout 590 python lab/vuoto_run3d.py"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from vuoto_core import chi_3d, chi_full_3d, species_weight

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
W2 = os.environ.get("VUOTO_WORKER2") == "1"          # second worker: reversed order, its own file (merged later)
FN = os.path.join(OUT, "vuoto_3d_b.json" if W2 else "vuoto_3d.json")
FO = os.path.join(OUT, "vuoto_3d.json" if W2 else "vuoto_3d_b.json")
T0 = time.time()
LIMIT = 540
DIRS = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1), "d": (1, 1, 1)}
THETAS = [0.1, 0.05, 0.025, 0.0125]
QEXAM = [0.3, 1, 3, 10]


def qtable(th):
    return list(np.round(np.exp(np.linspace(np.log(0.02), np.log(0.6 / th), 40)), 8))


def tasks():
    T = []
    for th in THETAS:
        beta = 6.0 if th < 0.02 else 5.0
        for dn in ("x", "y", "z", "d"):
            for qm in [0.05, 0.025] + QEXAM:
                T.append(("proj", th, dn, qm, 64, beta, 0.6, 1.05))
        for dn in ("x", "y", "z"):
            for qm in qtable(th):
                T.append(("proj", th, dn, qm, 64, beta, 0.6, 1.05))
    # whole sea (all species), theta = 0.1
    for qm in [0.05, 0.1, 0.2, 0.4, 1, 3]:
        T.append(("full", 0.1, "x", qm, 16, 0, 0, 0))
    # regulator variant
    for th in (0.05, 0.025):
        for dn in ("x", "y", "z"):
            for qm in [0.05, 0.025, 1, 3]:
                T.append(("proj", th, dn, qm, 64, 5.0, 0.5, 0.9))
    # quadrature convergence
    for qm in (0.05, 3):
        T.append(("proj", 0.0125, "x", qm, 80, 6.0, 0.6, 1.05))
    return T


def key(t):
    return "%s|%g|%s|%.8g|%d|%g|%g|%g" % t


def main():
    D = json.load(open(FN)) if os.path.exists(FN) else {}
    other = json.load(open(FO)) if os.path.exists(FO) else {}
    todo = [t for t in tasks() if key(t) not in D and key(t) not in other]
    if W2:
        todo = todo[::-1]
    print(f"{len(D)} done, {len(todo)} to do", flush=True)
    for t in todo:
        if time.time() - T0 > LIMIT:
            json.dump(D, open(FN, "w"), indent=0)
            print("checkpoint, rerun me", flush=True)
            sys.exit(3)
        try:
            other = json.load(open(FO)) if os.path.exists(FO) else {}
        except Exception:
            other = {}
        if key(t) in other:
            continue
        kind, th, dn, qm, n, beta, K1, K2 = t
        d = np.array(DIRS[dn], float); d /= np.linalg.norm(d)
        q = qm * th * d
        if kind == "full":
            c = chi_full_3d(q, th, nseg=n)
        else:
            c = chi_3d(q, th, species_weight(K1, K2), n=n, beta=beta)
        D[key(t)] = [c.real / (qm * th) ** 2, c.imag / (qm * th) ** 2]
        json.dump(D, open(FN, "w"), indent=0)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()

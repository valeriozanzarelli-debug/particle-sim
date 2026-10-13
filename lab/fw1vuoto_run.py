"""fw1vuoto_run.py - production runs pre-registered in lab/out/fw1vuoto_prereg.json (commit 04ef3cb).
Two workers (FW1V_WORKER = 0 / 1) take alternate tasks; each stores every evaluation in its own checkpoint file
lab/out/fw1vuoto_raw_w<i>.json (raw integral per site, complex, and q). Exit 3 = checkpoint, rerun me.
Run: sh lab/fw1vuoto_drive.sh   (two workers, <= 10 min pieces)"""
import os
import sys
import json
import time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fw1vuoto_core as c

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
WID = int(os.environ.get("FW1V_WORKER", "0"))
FN = os.path.join(OUT, f"fw1vuoto_raw_w{WID}.json")
T0 = time.time()
LIMIT = 540
TH = [0.1, 0.05, 0.025, 0.0125]
Q7 = [0.025, 0.05, 0.2, 0.3, 1, 3, 10]
DIRV = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1), "d": (1, 1, 1)}
AX = "xyz"


def beta_proj(th):
    return 5.0 if th > 0.02 else 6.0


def beta_cell(th):
    return 5.0 if th > 0.02 else 6.5


def tasks():
    T = []
    # 1. frozen + projector (validation of Debiti)
    for th in TH:
        for d in "xyz":
            for qm in Q7:
                T.append(("S", "F0", "ball", th, d, qm, "V", "V", 64, beta_proj(th)))
    # 2. frozen whole sea
    for th in (0.1, 0.05, 0.025):
        for qm in (0.05, 0.1, 0.2, 0.4, 1, 3, 10):
            T.append(("S", "F0", "cellS", th, "x", qm, "V", "V", 64, 6.0))
    for qm in (0.05, 3):
        T.append(("S", "F0", "cellS", 0.025, "x", qm, "V", "V", 96, 6.0))
    # 3. FW1 whole sea (main)
    for th in TH:
        for d in "xyz":
            for qm in Q7:
                T.append(("S", "FW1", "cellG", th, d, qm, "V", "V", 64, beta_cell(th)))
    for qm in (0.025, 0.05, 3):
        T.append(("S", "FW1", "cellG", 0.0125, "d", qm, "V", "V", 64, 6.5))
    for th in (0.025, 0.0125):
        for qm in (0.05, 3):
            T.append(("S", "FW1", "cellG", th, "x", qm, "V", "V", 80, beta_cell(th)))
    for qm in np.round(np.exp(np.linspace(np.log(0.02), np.log(48), 14)), 6):
        for d in "xyz":
            T.append(("S", "FW1", "cellG", 0.0125, d, float(qm), "V", "V", 64, 6.5))
    # 4. FW1 + projector (control)
    for th in TH:
        for d in "xyz":
            for qm in Q7:
                T.append(("S", "FW1", "ball", th, d, qm, "V", "V", 64, beta_proj(th)))
    # 5. potential position sensitivity
    for qm in (0.05, 1, 3):
        T.append(("S", "FW1mid", "cellG", 0.025, "x", qm, "V", "V", 64, 5.0))
    # 6. sea energy: chi_E
    for th in TH:
        for d in "xyz":
            for qm in (0.025, 0.05, 0.2, 0.3, 1, 3):
                T.append(("E", "FW1", "cellG", th, d, qm, "V", "V", 64, beta_cell(th)))
    # 7. sea energy: transverse magnetic kernel
    for th in TH:
        for i in range(3):
            for j in range(3):
                if j == i:
                    continue
                for qm in (0.025, 0.05):
                    T.append(("E", "FW1", "cellG", th, AX[i], qm, "A" + AX[j], "A" + AX[j], 64, beta_cell(th)))
    # 8. Ward (longitudinal)
    T.append(("E", "FW1", "cellG", 0.05, "x", 0.05, "Ax", "Ax", 64, 5.0))
    # 9. mixed V-A (odd: CFJ; even: magnetoelectric)
    for th in (0.05, 0.025, 0.0125):
        for i in range(3):
            for m in range(3):
                if m == i:
                    continue
                for qm in (0.05, 0.1):
                    for sg in ("", "-"):
                        T.append(("E", "FW1", "cellG", th, sg + AX[i], qm, "V", "A" + AX[m], 64, beta_cell(th)))
    # 10. antisymmetric A-A (b_0)
    for th in (0.05, 0.0125):
        for i in range(3):
            j, k = [a for a in range(3) if a != i]
            for qm in (0.05, 0.1):
                for sg in ("", "-"):
                    T.append(("E", "FW1", "cellG", th, sg + AX[i], qm, "A" + AX[j], "A" + AX[k], 64, beta_cell(th)))
    # 11. stroboscopic mixed density-to-A: frozen + projector (control, T symmetric) and FW1 whole sea
    for rule, mode in (("F0", "ball"), ("FW1", "cellG")):
        for i in range(3):
            for m in range(3):
                if m == i:
                    continue
                for qm in (0.05, 0.1):
                    T.append(("S", rule, mode, 0.05, AX[i], qm, "A" + AX[m], "A" + AX[m], 64, 5.0))
    return T


def key(t):
    return "%s|%s|%s|%g|%s|%.8g|%s|%s|%d|%g" % t


def pert(p):
    return "V" if p == "V" else ("A", "xyz".index(p[1]))


def evaluate(t):
    meth, rule, mode, th, d, qm, p1, p2, n, beta = t
    sg = -1.0 if d.startswith("-") else 1.0
    dv = np.array(DIRV[d.lstrip("-")], float)
    dv = sg * dv / np.linalg.norm(dv)
    q = qm * th * dv
    if rule == "F0":
        mv = c.frozen_moves(th)
    elif rule == "FW1":
        mv = c.fw1_moves(th)
    else:
        mv = c.fw1_moves(th, phi_mid=True)
    Ks, WW = c.quad_cell(mode, n, beta, q)
    if meth == "S":
        w = c.species_weight() if mode == "ball" else None
        val = c.integrate(lambda K: c.strobo_F(mv, K, q, w, pert(p1)), Ks, WW)
    else:
        val = c.integrate(lambda K: c.energy_F(mv, K, q, pert(p1), pert(p2)), Ks, WW, chunk=12000)
    return [float(val.real), float(val.imag), float(qm * th)]


def main():
    D = json.load(open(FN)) if os.path.exists(FN) else {}
    T = tasks()
    mine = [t for i, t in enumerate(T) if i % 2 == WID and key(t) not in D]
    print(f"worker {WID}: {len(D)} done, {len(mine)} to do (total tasks {len(T)})", flush=True)
    for t in mine:
        if time.time() - T0 > LIMIT:
            print("checkpoint, rerun me", flush=True)
            sys.exit(3)
        t1 = time.time()
        D[key(t)] = evaluate(t)
        json.dump(D, open(FN + ".tmp", "w"), indent=0)
        os.replace(FN + ".tmp", FN)
        print(key(t), D[key(t)][:2], "%.1fs" % (time.time() - t1), flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "count":
        print(len(tasks()))
    else:
        main()

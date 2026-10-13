"""protetta_griglia.py - EMENDAMENTO EM1 (dichiarato DOPO la corsa D1-D3, non cambia soglie): controllo di griglia sui
punti ad alto p/m, dove l'emissione quasi collineare (angolo ~ m/p) richiede piu' angoli.  Griglia fine nr=128,
nct=40, nph=64 contro quella di produzione (96, 20, 32).  -> lab/out/protetta_griglia.json"""
import os, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import protetta_vf3 as F
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
FN = os.path.join(OUT, "protetta_griglia.json")
D = json.load(open(FN)) if os.path.exists(FN) else {}
for cut, K in (("K12/3", (0.2, 0.35)), ("K12", (0.6, 1.05))):
    for th, pd in ((0.00078125, 0.05), (0.00078125, 10), (0.000390625, 0.05), (0.000390625, 76.8)):
        k = "R|%s|%g|%g" % (cut, th, pd)
        if k in D:
            continue
        t = time.time()
        r = F.dv_exact(th, "R", K, pd=pd, axes=(0, 1), nr=128, nct=40, nph=64)
        r["time_s"] = time.time() - t
        D[k] = r
        json.dump(D, open(FN, "w"), indent=1)
        print(k, r["dv_over_alpha"][:2], "%.0f s" % r["time_s"], flush=True)
print("ALL DONE")

"""tavg_f4_drive.py - FASE 4: tutte le corse pre-registrate (lab/out/tavg_prereg_f4.json), una (specie, X) per processo.
Uso: python tavg_f4_drive.py      (fa la prossima corsa non finita; stampa TUTTO_FATTO alla fine)
Ciclo non sorvegliato: bash lab/tavg_f4_loop.sh   (log non tracciato lab/out/tavg_f4.log)
Risultati: lab/out/tavg_f4_runs.json; ampiezze lab/out/tavg_cc_amp/*.npz (non tracciate)."""
import os, sys, json, time
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_cc as C

OUT = os.path.join(HERE, "out"); FN = os.path.join(OUT, "tavg_f4_runs.json")
PRE = json.load(open(os.path.join(OUT, "tavg_prereg_f4.json")))
ORDER = ["controllo_C", "controllo_O", "O", "O+", "N", "N+"]

if __name__ == "__main__":
    res = json.load(open(FN)) if os.path.exists(FN) else {}
    for sp in ORDER:
        Z, cfg = PRE["specie"][sp]
        for X in PRE["base"]["X"]:
            key = f"{sp}_X{X}"
            if res.get(key, {}).get("convergiuto"):
                continue
            t0 = time.time()
            r = C.run(int(Z), cfg, X, name=sp.replace("+", "p"), t_budget=520, log=lambda s: None)
            r["secondi"] = time.time() - t0
            res = json.load(open(FN)) if os.path.exists(FN) else {}
            res[key] = r
            json.dump(res, open(FN, "w"), indent=1)
            print(f"[{key}] convergiuto={r['convergiuto']} ({r['secondi']:.0f} s)", flush=True)
            sys.exit(0)
    print("TUTTO_FATTO", flush=True)

"""tavg_simm.py - la rottura di simmetria dell'ossigeno sulla griglia: e' fisica (HF ristretto a orbitali non
equivalenti) o e' la griglia? Prova decisiva:
  1. energia sulla griglia con gli orbitali del livello radiale (tutti i p con LA STESSA funzione radiale), senza
     rilassare (maxiter = 0): deve ridare l'energia radiale dello pseudo-atomo (errore di griglia puro);
  2. stessa partenza, rilassamento libero: quanto scende e quanto diventano diversi i raggi <r> dei p;
  3. ripetuto a due celle (a0/4, a0/6): se la discesa non dipende dalla cella non e' la griglia.
Uscita: lab/out/tavg_simm.json (checkpoint per caso). Uso: python tavg_simm.py"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_grid as tg
import tavg_run_grid as rg

FN = os.path.join(HERE, "out", "tavg_simm.json")


def mean_r(res):
    mc = res["machine"]
    X, Y, Z = np.meshgrid(mc.x, mc.x, mc.x, indexing="ij", sparse=True)
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    return [float(np.sum(p * p * r) * mc.dV) for p in res["orbitali"]]


def case(Z, h, L, maxiter):
    blk = tg.load_block(Z)
    spin, order, val = rg.occupation(Z, 0, blk)
    r = tg.hf_griglia([(Z, (0.0, 0.0, 0.0))], spin=spin, h=h, L=L, tol=1e-8, maxiter=maxiter, order=order,
                      log=lambda *a: None)
    return dict(E=r["E"], E_val=r["E_val"], dE_vs_radiale_mEh=1e3 * (r["E_val"] - blk["E_ps_ref"]),
                r_medio=mean_r(r), occ=r["occ"], convergito=r["convergito"], iterazioni=len(r["storia"]) - 1,
                secondi=r["secondi"], h=h, L=L)


def main():
    res = json.load(open(FN)) if os.path.exists(FN) else {}
    t0 = time.time()
    for Z, h, L in ((8, 0.25, 16.0), (8, 1 / 6, 12.0)):
        for tag, mi in (("fermo", 0), ("rilassato", 800)):
            key = f"Z{Z}_h{h:.4f}_L{L:.0f}_{tag}"
            if key in res:
                continue
            if time.time() - t0 > 200:
                print("checkpoint: rilanciami", flush=True); sys.exit(3)
            res[key] = case(Z, h, L, mi)
            json.dump(res, open(FN, "w"), indent=1)
            print(key, {k: res[key][k] for k in ("E_val", "dE_vs_radiale_mEh", "r_medio", "iterazioni", "secondi")},
                  flush=True)
    # diagnostica (non contata): l'effetto uova di G6 a cella a0/6, stessa scatola L = 16
    for tag, sh in (("uova_h6_centro", (0, 0, 0)), ("uova_h6_spostato", (0.37, 0.21, 0.11))):
        if tag in res:
            continue
        if time.time() - t0 > 200:
            print("checkpoint: rilanciami", flush=True); sys.exit(3)
        o = rg.run_atom(8, 0, h=1 / 6, L=16.0, shift=sh, log=lambda *a: None)
        res[tag] = dict(E_val=o["E_val"], E=o["E"], convergito=o["convergito"], h=1 / 6, L=16.0, shift_celle=sh)
        json.dump(res, open(FN, "w"), indent=1)
        print(tag, o["E_val"], flush=True)
    print("SIMM DONE", flush=True)


if __name__ == "__main__":
    main()

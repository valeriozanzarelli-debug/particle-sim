"""tavg_blocks_prod.py - produzione dei blocchi-nucleo con la scelta automatica della cella h (pre-registrata):
la cella piu' grossa in [a0/4, a0/6, a0/8] il cui blocco trasferisce il catione entro 0.8 mEh. Esiti in
lab/out/tavg_blocks.json (checkpoint per elemento), blocchi in lab/out/tavg_blocchi/Z??.json.
Uso: python tavg_blocks_prod.py Z1 Z2 ...   (exit 3 = rilanciami)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_radial as tr
import tavg_core as tc
import tavg_blocks as tb

OUT = os.path.join(HERE, "out")
FN = os.path.join(OUT, "tavg_blocks.json")
HS = [0.25, 1 / 6, 0.125]
SEL_BUDGET = 0.8e-3


def ps_energy(b, cfg_full, nc):
    val = tr.config_from_string(cfg_full)[nc:] if cfg_full else []
    val = [v for v in val if v[2] > 0]
    if not val:
        return 0.0
    p = b.pseudo_atom(tc.relabel(val))
    p.scf()
    return p.E


def promotion_config(Z):
    """ns1 np^(n+1) for B..Ne, Na..Ar, Ga..Kr (None otherwise)"""
    if not (5 <= Z <= 10 or 13 <= Z <= 18 or 31 <= Z <= 36):
        return None
    cfg = tb.ground_config(Z)
    toks = cfg.split()
    s = [i for i, t in enumerate(toks) if t[1] == "s"][-1]
    p = [i for i, t in enumerate(toks) if t[1] == "p" and t[0] == toks[s][0]]
    if not p:
        return None
    q = int(toks[p[0]][2:])
    if q >= 6:
        return None
    toks[s] = toks[s][:2] + "1"; toks[p[0]] = toks[p[0]][:2] + str(q + 1)
    return " ".join(toks)


def main():
    Zs = [int(x) for x in sys.argv[1:]]
    res = json.load(open(FN)) if os.path.exists(FN) else {}
    t0 = time.time()
    for Z in Zs:
        if res.get(str(Z), {}).get("done"):
            continue
        if time.time() - t0 > 300:
            print("checkpoint: rilanciami", flush=True); sys.exit(3)
        ts = time.time()
        cfg0 = tb.ground_config(Z); cfg1 = tb.ground_config(Z, 1)
        nc = tb.core_of(Z)[1]
        ae = tb.get_ae(Z, cfg0)
        E1 = None
        if Z > 1:
            a1 = tb.get_ae(Z, cfg1, grid=ae.g); E1 = a1.E
        rec = dict(Z=Z, config=cfg0, config_catione=cfg1, E_fine=ae.E, E_fine_catione=E1, prove=[])
        chosen = None
        for h in HS:
            b, d = tb.learn(Z, h_grid=h, ae=ae, write=False, log=lambda *a: None)
            if Z == 1:
                err = 0.0
            else:
                Ep1 = ps_energy(b, cfg1, nc)
                err = (Ep1 - b.ps_ref["E"]) - (E1 - ae.E)
            rec["prove"].append(dict(h=h, rc={str(k): v for k, v in b.rc.items()}, err_catione_Eh=err))
            print(f"Z={Z} h={h:.4f} rc={ {k: round(v, 3) for k, v in b.rc.items()} } errore catione {1e3 * err:+.3f} mEh",
                  flush=True)
            if abs(err) <= SEL_BUDGET:
                chosen = (h, b, d)
                break
        if chosen is None:
            chosen = (h, b, d); rec["selezione_fallita"] = True
        h, b, d = chosen
        d["h_grid"] = h
        d["err_catione_Eh"] = rec["prove"][-1]["err_catione_Eh"]
        # B3: promotion ns -> np
        pc = promotion_config(Z)
        if pc:
            ap = tb.get_ae(Z, pc, grid=ae.g)
            Epp = ps_energy(b, pc, nc)
            rec["B3"] = dict(config=pc, dE_fine=ap.E - ae.E, dE_blocco=Epp - b.ps_ref["E"],
                             err_Eh=(Epp - b.ps_ref["E"]) - (ap.E - ae.E))
            d["err_promozione_Eh"] = rec["B3"]["err_Eh"]
        os.makedirs(tb.BLOCKDIR, exist_ok=True)
        json.dump(d, open(os.path.join(tb.BLOCKDIR, f"Z{Z:02d}.json"), "w"))
        rec.update(h=h, rc={str(k): v for k, v in b.rc.items()}, E_core=b.E_core, E_ps_ref=b.ps_ref["E"],
                   eps_fine={str(k): v for k, v in b.eps.items()}, done=True, secondi=time.time() - ts)
        res[str(Z)] = rec
        json.dump(res, open(FN, "w"), indent=1)
        print(f"== Z={Z}: h = {h:.4f}, {time.time() - ts:.0f} s", flush=True)
    print("BLOCKS DONE", flush=True)


if __name__ == "__main__":
    main()

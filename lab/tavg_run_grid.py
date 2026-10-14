"""tavg_run_grid.py - corse sulla griglia 3-D (esami G1..G7 della pre-registrazione tavg_).
Uso: python tavg_run_grid.py atom Z [Z ...] | conv Z ... | learned Z ... | shift Z ... | mol H2 N2 ...
Esiti in lab/out/tavg_grid.json (checkpoint per compito). exit 3 = rilanciami."""
import os, sys, json, time, math
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_radial as tr
import tavg_grid as tg
import tavg_blocks as tb

OUT = os.path.join(HERE, "out")
FN = os.path.join(OUT, "tavg_grid.json")
CHAN = {0: [0], 1: [1, 2, 3], 2: [4, 5, 6, 7, 8]}


def load():
    return json.load(open(FN)) if os.path.exists(FN) else {}


def save(res):
    json.dump(res, open(FN, "w"), indent=1)


def valence(Z, charge):
    cfg = tb.ground_config(Z, charge)
    nc = tb.core_of(Z)[1]
    return [v for v in tr.config_from_string(cfg)[nc:] if v[2] > 0] if cfg else []


def occupation(Z, charge, block):
    """spin (2S) and the order of the guess functions: doubly occupied orbitals first, then singly occupied (Hund)"""
    val = valence(Z, charge)
    det = tr.hund_determinant(val)
    spin = sum(s for _, _, s in det)
    chans = sorted(int(l) for l in block["phi"])
    # function index of channel l in atomic_guess (channels in order s, p, d among those present)
    base = {}; i = 0
    for l in (0, 1, 2):
        if l in chans:
            base[l] = list(range(i, i + 2 * l + 1)); i += 2 * l + 1
    dbl, sgl = [], []
    for (n, l, q) in val:
        nd = max(0, q - (2 * l + 1)); ns = q - 2 * nd
        fs = base[l]
        dbl += fs[:nd]; sgl += fs[nd:nd + ns]
    return spin, dbl + sgl, val


def spherical(val):
    """valence density spherical for this determinant? (closed or half-filled shells only)"""
    return all(q in (0, 2 * l + 1, 2 * (2 * l + 1)) or l == 0 for (n, l, q) in val)


def box(Z):
    return 16.0 if Z <= 10 else 20.0


def ps_cation_energy(Z, brec):
    return brec["E_ps_ref"] + (brec["E_fine_catione"] - brec["E_fine"]) + brec["prove"][-1]["err_catione_Eh"]


def run_atom(Z, charge, h=None, L=None, kinetic="regola", shift=(0, 0, 0), log=print):
    blk = tg.load_block(Z)
    spin, order, val = occupation(Z, charge, blk)
    h = h or blk["h_grid"]; L = L or box(Z)
    pos = tuple(s * h for s in shift)
    nel = sum(q for _, _, q in val)
    if nel == 0:
        return dict(E_val=0.0, E=blk["E_core"], secondi=0.0, convergito=True, h=h, L=L)
    r = tg.hf_griglia([(Z, pos)], carica=charge, spin=spin, h=h, L=L, kinetic=kinetic, tol=1e-8, maxiter=600,
                      order=order, log=log)
    return dict(E_val=r["E_val"], E=r["E"], eps=r["eps"], occ=r["occ"], convergito=r["convergito"],
                secondi=r["secondi"], h=h, L=L, N=r["N"], parti=r["parti"], spin=spin,
                sferico=spherical(val), iterazioni=len(r["storia"]))


def task(res, key, fn):
    if res.get(key, {}).get("done"):
        return
    ts = time.time()
    out = fn()
    out["done"] = True
    out["secondi_tot"] = time.time() - ts
    res[key] = out
    save(res)
    print(f"{key}: " + ", ".join(f"{k}={v:.8f}" if isinstance(v, float) else f"{k}={v}" for k, v in out.items()
                                  if k in ("E", "E_val", "convergito", "secondi_tot", "h", "N")), flush=True)


def main():
    mode = sys.argv[1]
    args = sys.argv[2:]
    res = load()
    t0 = time.time()
    quiet = lambda *a: None
    for a in args:
        if time.time() - t0 > 360:
            print("checkpoint: rilanciami", flush=True); sys.exit(3)
        if mode == "atom":
            Z = int(a)
            for q in (0, 1):
                task(res, f"atom_Z{Z}_q{q}", lambda: run_atom(Z, q, log=quiet))
        elif mode == "conv":
            Z = int(a)
            h = tg.load_block(Z)["h_grid"]
            task(res, f"conv_Z{Z}_h1", lambda: run_atom(Z, 0, h=h, L=12.0, log=quiet))
            task(res, f"conv_Z{Z}_h2", lambda: run_atom(Z, 0, h=h * 2 / 3, L=12.0, log=quiet))
        elif mode == "learned":
            Z = int(a)
            task(res, f"learned_Z{Z}", lambda: run_atom(Z, 0, h=0.25, L=16.0, kinetic="imparato", log=quiet))
        elif mode == "shift":
            Z = int(a)
            task(res, f"shift_Z{Z}", lambda: run_atom(Z, 0, shift=(0.37, 0.21, 0.11), log=quiet))
        elif mode == "mol":
            geo = {"H2": [(1, (0, 0, -0.7)), (1, (0, 0, 0.7))],
                   "N2": [(7, (0, 0, -1.034)), (7, (0, 0, 1.034))]}[a]

            def fn():
                h = min(tg.load_block(Z)["h_grid"] for Z, _ in geo)
                r = tg.hf_griglia(geo, h=h, L=16.0, tol=1e-8, maxiter=800, log=quiet)
                return dict(E=r["E"], E_val=r["E_val"], E_nn=r["E_nn"], eps=r["eps"], convergito=r["convergito"],
                            secondi=r["secondi"], h=h, N=r["N"])
            task(res, f"mol_{a}", fn)
    print("GRID DONE", flush=True)


if __name__ == "__main__":
    main()

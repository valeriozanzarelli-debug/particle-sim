"""tavg_f5_drive.py - FASE 5 (lab/out/tavg_prereg_f5.json): una corsa per processo (<= 590 s), ripresa da checkpoint.
Uso: python tavg_f5_drive.py ; ciclo: bash lab/tavg_f5_loop.sh (log non tracciato lab/out/tavg_f5.log)
Risultati lab/out/tavg_f5_runs.json; ampiezze lab/out/tavg_cc_amp/f5_*.npz (non tracciate)."""
import os, sys, json, time
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_cc as C

OUT = os.path.join(HERE, "out"); FN = os.path.join(OUT, "tavg_f5_runs.json")
PRE = json.load(open(os.path.join(OUT, "tavg_prereg_f5.json")))
TASKS = [("He_Z8", "AEc", 3), ("He_Z8", "AEc", 4)]
for sp in ("O", "O+", "N", "N+"):
    TASKS += [(sp, "FCv", 4), (sp, "FCv", 5), (sp, "AEc", 3), (sp, "FCc", 3), (sp, "AEc", 4), (sp, "FCc", 4)]
TASKS += [("He_Z8", "AEc", 5)] + [(sp, k, 5) for sp in ("O", "O+", "N", "N+") for k in ("AEc", "FCc")]   # emendamento 1
KIND = {"FCv": dict(core=False, frozen_core=True), "AEc": dict(core=True, frozen_core=False),
        "FCc": dict(core=True, frozen_core=True)}

if __name__ == "__main__":
    res = json.load(open(FN)) if os.path.exists(FN) else {}
    order = TASKS[::-1] if "--rev" in sys.argv else TASKS
    for sp, kind, X in order:
        key = f"{sp}_{kind}_X{X}"
        if res.get(key, {}).get("convergiuto"):
            continue
        Z, cfg = PRE["specie"][sp]
        t0 = time.time()
        r = C.run(int(Z), cfg, X, name="f5_" + key.replace("+", "p"), t_budget=420, log=lambda s: None, ae=True, **KIND[kind])
        r["secondi"] = time.time() - t0
        import resource; r["RSS_MB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
        res = json.load(open(FN)) if os.path.exists(FN) else {}
        res[key] = r; json.dump(res, open(FN, "w"), indent=1)
        print(f"[{key}] convergiuto={r['convergiuto']} ({r['secondi']:.0f} s, {r['RSS_MB']} MB)", flush=True)
        sys.exit(0)
    print("TUTTO_FATTO", flush=True)

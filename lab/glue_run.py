"""glue_run.py - production Monte Carlo for one beta (resumable, time-boxed).

usage: python glue_run.py BETA [budget_seconds] [planfile]
Reads the run plan from lab/out/glue_prereg.json. Each call runs at most `budget` seconds (default 540),
writes one chunk of measurements to lab/out/glue_raw/b{beta}_c{k:03d}.npz plus a lattice checkpoint, and exits.
Exit code 0 when the planned number of measurements is complete, 3 when more calls are needed.
Single thread (NUMBA_NUM_THREADS=1, OMP_NUM_THREADS=1 set below).
"""
import os, sys, json, time, glob
os.environ.setdefault("NUMBA_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import glue_core as g

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out"); RAW = os.path.join(OUT, "glue_raw")
os.makedirs(RAW, exist_ok=True)


def main():
    beta = float(sys.argv[1]); budget = float(sys.argv[2]) if len(sys.argv) > 2 else 540.0
    t_start = time.time()
    planfile = sys.argv[3] if len(sys.argv) > 3 else "glue_prereg.json"
    plan = json.load(open(os.path.join(OUT, planfile)))["runs"][f"{beta:g}"]
    T, L = plan["T"], plan["L"]
    tag = ("pilot_" if "pilot" in planfile else "") + f"b{beta:g}"
    if len(sys.argv) > 4:          # second independent chain (declared deviation): tag suffix + own seed/length
        tag += sys.argv[4]; plan = dict(plan, seed=plan["seed"] + 1, n_meas=int(sys.argv[5]))
    ck = os.path.join(RAW, f"{tag}_state.npz")
    chunks = sorted(glob.glob(os.path.join(RAW, f"{tag}_c*.npz")))
    nmeas_done = sum(int(np.load(c)["nmeas"]) for c in chunks)
    if nmeas_done >= plan["n_meas"]:
        print(tag, "complete", nmeas_done); return 0
    k = len(chunks)
    lat = g.FastLattice(T, L, beta, seed=plan["seed"] * 1000 + k)   # new RNG stream per chunk
    if os.path.exists(ck):
        st = np.load(ck); lat.U[...] = st["U"]; therm_done = bool(st["therm"])
    else:
        therm_done = False
    if not therm_done:
        for i in range(plan["n_therm"]):
            lat.heatbath()
            for _ in range(plan["n_or"]):
                lat.overrelax()
        lat.reunit()
        np.savez(ck, U=lat.U, therm=True)
        print(tag, "thermalised", lat.plaquette(), time.time() - t_start, flush=True)
    ops_l, pops_l, plaq_l, W_l, act_l = [], [], [], [], []
    names = pnames = None
    n = 0
    tmeas = 0.0
    while nmeas_done + n < plan["n_meas"]:
        t0 = time.time()
        for _ in range(plan["sep"]):
            lat.heatbath()
            for _ in range(plan["n_or"]):
                lat.overrelax()
        lat.reunit()
        ops, pops = g.glue_ops(lat, plan["smear_levels"], plan["ape_alpha"], plan["sq_sizes"], plan["rect_sizes"])
        if names is None:
            names = list(ops.keys()); pnames = list(pops.keys())
        ops_l.append(np.array([ops[q] for q in names], dtype=np.float64))
        pops_l.append(np.array([pops[q] for q in pnames], dtype=np.complex64))
        act = lat.action_density_slices()
        plaq_l.append(1.0 - act.sum() / (3 * T * L * L))
        act_l.append(act.astype(np.float32))
        if (nmeas_done + n) % plan["creutz_every"] == 0:
            W_l.append(g.wilson_loops(lat, plan["creutz_Rmax"]))
        n += 1
        tmeas = max(tmeas, time.time() - t0)
        if time.time() - t_start + 2 * tmeas > budget:
            break
    np.savez(os.path.join(RAW, f"{tag}_c{k:03d}.npz"), nmeas=n, names=np.array(names), pnames=np.array(pnames),
             ops=np.array(ops_l), pops=np.array(pops_l), plaq=np.array(plaq_l), act=np.array(act_l),
             W=np.array(W_l) if W_l else np.zeros((0, plan["creutz_Rmax"] + 1, plan["creutz_Rmax"] + 1)))
    np.savez(ck, U=lat.U, therm=True)
    tot = nmeas_done + n
    print(tag, f"chunk {k}: +{n} meas -> {tot}/{plan['n_meas']}, <P>={np.mean(plaq_l):.6f}, "
               f"{(time.time() - t_start) / max(n, 1):.3f} s/meas", flush=True)
    return 0 if tot >= plan["n_meas"] else 3


if __name__ == "__main__":
    sys.exit(main())

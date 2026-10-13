"""h2v_wf.py - VIDEO H2: rebuild the two-electron wavefunctions of H2 that Debiti computed (lab/aggmol_h2.py), with
checkpoints, so the video can sample them. Same machine (lab/aggmol_ecg.py), same seeds, same sizes, same warm-start
chain as aggmol_h2.py job_curve: the energies must come back equal to lab/out/aggmol_h2.json (exam V0, pre-registered).
Bases are written to lab/out/h2v_basis_R*.npz (big binary: never committed).

    timeout 590 nice -n 10 python3 lab/h2v_wf.py        # run again until it prints ALL DONE
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import aggmol_ecg as E

OUT = E.OUT
JS = os.path.join(OUT, "h2v_wf.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
# warm-start chain of aggmol_h2.py job_curve restricted to what the video needs (R in [1.2, 1.7])
CHAIN = [(1.4, None), (1.35, 1.4), (1.3, 1.35), (1.2, 1.3), (1.45, 1.4), (1.5, 1.45), (1.6, 1.5), (1.7, 1.6)]


def bfile(R):
    return os.path.join(OUT, f"h2v_basis_R{R:.2f}.npz")


def ckfile(R):
    return os.path.join(OUT, f"h2v_ck_R{R:.2f}.npz")


def save():
    json.dump(res, open(JS, "w"), indent=1, default=float)


def save_ck(R, B, rng, stage, nref_done):
    np.savez(ckfile(R), A=B.A, s=B.s, nrm=B.nrm, S=B.S, H=B.H)
    res.setdefault("ck", {})[f"{R:.2f}"] = dict(rng=rng.bit_generator.state, stage=stage, nref_done=nref_done)
    save()


def run(R, Rn, t0, tmax):
    """reproduce aggmol_ecg.build exactly (seed 1), but resumable"""
    key = f"{R:.2f}"
    nuc = [(0, 0, R / 2), (0, 0, -R / 2)]
    B = E.Basis(2, nuc)
    gen = E.make_gen(2, max(R, 1.0))
    rng = np.random.default_rng(1)
    nref = 8 if Rn is None else 4
    ck = res.get("ck", {}).get(key)
    if ck:
        z = np.load(ckfile(R))
        B.A, B.s, B.nrm, B.S, B.H = z["A"], z["s"], z["nrm"], z["S"], z["H"]
        rng.bit_generator.state = ck["rng"]
        nref_done = ck["nref_done"]
    else:
        nref_done = 0
        if Rn is not None:
            st = np.load(bfile(Rn))
            B.A, B.s = st["A"].copy(), st["s"] * (R / Rn)
            B.rebuild()
    e = np.inf
    while len(B.A) < 180:
        if time.time() - t0 > tmax:
            save_ck(R, B, rng, "grow", nref_done); return None
        e = B.add_best(gen, 40, rng)
    while nref_done < nref:
        if time.time() - t0 > tmax:
            save_ck(R, B, rng, "refine", nref_done); return None
        e = B.refine(gen, 40, rng)
        nref_done += 1
        save_ck(R, B, rng, "refine", nref_done)
    e = B.energy(B.H, B.S)
    # ground-state coefficients (normalised basis functions times nrm)
    w, V = __import__("scipy.linalg", fromlist=["eigh"]).eigh(B.H, B.S)
    np.savez(bfile(R), A=B.A, s=B.s, nrm=B.nrm, c=V[:, 0], E_el=w[0])
    return float(e)


def main():
    t0 = time.time()
    ref = json.load(open(os.path.join(OUT, "aggmol_h2.json")))["curve"]
    done = res.setdefault("E", {})
    for R, Rn in CHAIN:
        key = f"{R:.2f}"
        if key in done:
            continue
        print(f"R = {R} ...", flush=True)
        e = run(R, Rn, t0, 400)
        if e is None:
            print("checkpoint: run again", flush=True); return
        Et = e + 1 / R
        done[key] = dict(E_el=e, E_tot=Et, E_debiti=ref[key]["E_tot"], diff=Et - ref[key]["E_tot"], sec_total=time.time() - t0)
        res.get("ck", {}).pop(key, None)
        save()
        print(f"R = {R}: E_tot {Et:.10f}  Debiti {ref[key]['E_tot']:.10f}  diff {Et - ref[key]['E_tot']:.2e}", flush=True)
        if time.time() - t0 > 300:
            print("checkpoint: run again", flush=True); return
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()

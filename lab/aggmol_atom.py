"""aggmol_atom.py - task (2), exams on ONE atom with the new a0/4 level (study aggmol_): 1s, 2s, 2p1/2, 2p3/2 of hydrogen
(box +-16 a0, N = 128 at a0/4), leading order (learned kinetic + gamma0 Galerkin field) plus the first order of the rest of
the learned coupling (Darwin, spin-orbit, ...), exactly as aggfull_run.run_pert. Compared with the fine zig-zag walk at
the same theta (lab/out/aggfull_fineref.json: radial fine rule) and with Dirac.
Jobs: s (1s + 2s, block of two in the even sector), p (2p_z, odd in z), pert_s, pert_p.  Results: lab/out/aggmol_atom.json
Usage: timeout 590 python lab/aggmol_atom.py JOB [JOB ...] [--f 4]"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_core as mc
import aggmol_fz_atom as fa

PFX = "dev_" if mc.DEV else ""
JS = os.path.join(mc.OUT, PFX + "aggmol_atom.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
F = 4
if "--f" in sys.argv:
    F = int(sys.argv[sys.argv.index("--f") + 1])
B = {16: 256, 8: 512, 4: 1024}[F]
H = B * mc.TH
N = int(round(32 * mc.A0 / H))          # box +-16 a0


def save():
    json.dump(res, open(JS, "w"), indent=1, default=float)


def vec(name):
    return os.path.join(mc.OUT, PFX + f"aggmol_atomvec_{name}_a0_{F}.npy")


def setup():
    rec = mc.level_rec(B)
    C, offs, g0 = mc.coupling_rec(B)
    T = mc.kinetic_grid(rec, N, H)
    f = mc.MultiField(N, H, [(0, 0, 0)])
    return T, f, C, offs, g0.real


def job_states(kind):
    T, f, C, offs, g0 = setup()
    sym = (1, 1, 1) if kind == "s" else (1, 1, -1)
    m = mc.Mol(T, f, gamma0=g0, sym=sym)
    if kind == "s":
        X0 = np.stack([mc.guess_mol(N, H, [(0, 0, 0)], "1s"), mc.guess_mol(N, H, [(0, 0, 0)], "2s")], 1)
        if os.path.exists(vec("1s")):
            X0 = np.stack([np.load(vec("1s")), np.load(vec("2s"))], 1)
    else:
        X0 = (np.load(vec("2p")) if os.path.exists(vec("2p")) else mc.guess_mol(N, H, [(0, 0, 0)], "2p"))[:, None]
    X, E, r, conv, its = mc.solve_until(m, X0, 520, rtol=2e-8, tag=f"{kind} a0/{F}")
    names = ["1s", "2s"] if kind == "s" else ["2p"]
    out = {}
    for i, nm in enumerate(names):
        np.save(vec(nm), X[:, i])
        out[nm] = dict(E0_eV=E[i] * mc.MC2, residual=r[i])
    return dict(states=out, converged=conv, its=its, N=N, f=F, gamma0=g0)


def job_pert(kind):
    T, f, C, offs, g0 = setup()
    out = {}
    if kind == "s":
        for nm in ("1s", "2s"):
            v = np.load(vec(nm)).reshape((N,) * 3)
            v = v / np.linalg.norm(v)
            S = fa.coupling_blocks(f, [v], offs)
            M = fa.first_order_matrix(C, offs, S, g0)
            ev = np.linalg.eigvalsh(M)
            out[nm] = dict(dV_eV=(ev * mc.MC2).tolist())
    else:
        pz = np.load(vec("2p")).reshape((N,) * 3)
        pz = pz / np.linalg.norm(pz)
        phis = [np.transpose(pz, (2, 1, 0)), np.transpose(pz, (0, 2, 1)), pz]
        S = fa.coupling_blocks(f, phis, offs)
        M = fa.first_order_matrix(C, offs, S, g0)
        ev = np.linalg.eigvalsh(M)
        out["2p"] = dict(dV_eV=(ev * mc.MC2).tolist(), splitting_eV=float((ev[2:].mean() - ev[:2].mean()) * mc.MC2),
                         interval_ratio=float((ev[2:].mean() - ev.mean()) / (ev.mean() - ev[:2].mean())))
    return out


if __name__ == "__main__":
    for job in [a for a in sys.argv[1:] if not a.startswith("--") and not a.isdigit()]:
        key = f"{job}_a0_{F}"
        if key in res and (res[key].get("converged", True)):
            print(f"{key}: done", flush=True)
            continue
        t0 = time.time()
        out = job_states(job) if job in ("s", "p") else job_pert(job.split("_")[1])
        out["seconds"] = time.time() - t0
        res[key] = out
        save()
        print(key, json.dumps(out, default=float)[:600], flush=True)

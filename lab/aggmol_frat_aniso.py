"""aggmol_frat_aniso.py - study aggmol_frat_ (fractal test with H2), fine level: the FULL anisotropic H2-H2 force.
Phase 3 (aggmol_h2h2.py) computed 4 orientations (H, X, T, L). A rank-2 anisotropic pair force has 5 invariant
coefficients, so one more orientation is needed: S = both molecules parallel, tilted 45 deg from the axis.
Same machine, same monomer (lab/out/aggmol_h2h2_mono.npz, not committed; `python lab/aggmol_h2h2.py mono` remakes it).
Also: the triple-dipole (Axilrod-Teller) coefficient C9 from the same pseudo-states (three-body force, for the trimers).
Jobs: geom (checkpointed) | c9.  Output lab/out/aggmol_frat_aniso.json"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_h2h2 as M

JS = os.path.join(M.OUT, "aggmol_frat_aniso.json")
s2 = 1 / np.sqrt(2)
M.ORIENT["S"] = ((s2, 0, s2), (s2, 0, s2))


def load():
    return json.load(open(JS)) if os.path.exists(JS) else {}


def save(r):
    json.dump(r, open(JS, "w"), indent=1, default=float)


def job_geom(budget=360):
    t0 = time.time()
    mono = dict(np.load(M.NPZ)); spins = M.spin_factors()
    for R in M.RGRID:
        r = load(); cur = r.setdefault("geom", {})
        key = f"S_{R:.1f}"
        if key in cur:
            continue
        if time.time() - t0 > budget:
            print("checkpoint: run again", flush=True); return
        t1 = time.time()
        out = M.geometry("S", R, mono, spins)
        out["sec"] = time.time() - t1
        out["E_int"] = out["E1_HL"] + out["E_ind"] + out["E_disp"]
        cur[key] = out; save(r)
        print(key, f"E_int {out['E_int']:.4e}", flush=True)
    print("ORIENT S DONE", flush=True)


def job_c9():
    """isotropic C6 and C9 from the dipole pseudo-states: alpha(iw) = sum 2 dE F / (dE^2 + w^2) (F = |d|^2/3 summed
    over x, y, z; isotropic average), C6 = 3/pi int alpha^2, C9 = 3/pi int alpha^3 (Gauss-Legendre on w = t/(1-t))"""
    mono = dict(np.load(M.NPZ))
    A, s, c0, D, dE, ng = mono["A"], mono["s"], mono["c0"], mono["D"], mono["dE"], int(mono["nground"])
    nuc = [np.array((0, 0, M.RB / 2)), np.array((0, 0, -M.RB / 2))]
    g = M.transition_gaussians(A, s, c0, np.concatenate([c0[:, None], D], 1), ng, nuc)
    dip = np.einsum("gm,gx->mx", g["T"], g["c"])[1:]
    F = np.sum(dip ** 2, 1) / 3
    x, w = np.polynomial.legendre.leggauss(400)
    t = 0.5 * (x + 1); wt = 0.5 * w
    om = t / (1 - t); jac = 1 / (1 - t) ** 2
    al = np.array([np.sum(2 * dE * F / (dE ** 2 + o ** 2)) for o in om])
    C6 = 3 / np.pi * np.sum(wt * jac * al ** 2)
    C9 = 3 / np.pi * np.sum(wt * jac * al ** 3)
    r = load(); r["C6_quad"] = float(C6); r["C9"] = float(C9); r["alpha0_iso"] = float(al[0]); save(r)
    print("C6", C6, "C9", C9, "alpha(0)", al[0])


if __name__ == "__main__":
    dict(geom=job_geom, c9=job_c9)[sys.argv[1]]()

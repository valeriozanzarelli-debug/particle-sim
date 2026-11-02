"""ciclo_p3_run.py - phase 3 production: python lab/ciclo_p3_run.py <Q|P|P2|F> <mix|coarse|b3|truth>
-> lab/out/ciclo_p3_<scenario>_<kind>.json   (each call well under 10 min; settings = lab/out/ciclo_prereg.json, fase3)"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
import scipy.sparse as sps
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ciclo_p3_mix as M

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
N0 = M.N0
XC = N0 / 2                     # atom centre (base-cell units)
DT, T, TP = 25.0, 30000.0, 15700.0
NCHK, EPS_REF, BUF, MARGIN, PATCH = 16, 1e-2, 27, 1, 3
SIGMA = 2.0
R_ION = 120.0                   # "ionised" = farther than this from the atom (7 rms radii of the ground state)


def base_G(k):
    g = M.p3.hierarchy(2)["g"]
    return g[0] + 2 * sum(g[j] * np.cos(j * k) for j in range(1, len(g)))


def scenario(name, E0):
    if name == "Q":
        return M.Scenario(0.0, XC, SIGMA, 2.4e-3, TP, T, DT)
    if name == "P":
        return M.Scenario(2e-3, XC, SIGMA, 2.4e-3, TP, T, DT)
    if name == "P2":
        return M.Scenario(3e-3, XC + 12, SIGMA, float(base_G(0.7) - E0), TP, T, DT)
    if name == "F":
        return M.Scenario(2e-3, XC + 600, SIGMA, 2.4e-3, TP, T, DT)
    raise ValueError(name)


def main(name, kind):
    pr = M.prolongators()
    P, P1 = pr["P"], pr["P1"]
    Hs = (pr["Hb"] + sps.diags(pr["V0"])).tocsr()
    G = P.conj().T @ (Hs @ P)
    e, U = np.linalg.eigh(G)
    a0 = P @ U[:, 0]                          # start: the atom's ground state ON THE COARSE LEVEL, prolonged
    sc = scenario(name, e[0])
    x = np.arange(N0) + 0.5
    far = np.abs(x - XC) > R_ION
    t0 = time.time()
    res = dict(scenario=name, kind=kind, amp=sc.amp, xd=sc.xd, omega=sc.omega, sigma=SIGMA, dt=DT, T=T, Tp=TP,
               nchk=NCHK, eps_ref=EPS_REF, buf=BUF, margin=MARGIN, patch=PATCH, E0_coarse=float(e[0]))
    if kind == "truth":
        sn, _ = M.truth(sc, Hs, a0, NCHK)
        np.save(os.path.join(OUT, f"ciclo_p3_{name}_truth.npy"), sn)          # not committed (large)
        res.update(cost=float(N0 * sc.nsteps))
    else:
        Pk = P1 if kind == "b3" else P
        r = M.run_mixed(sc, Pk, Hs, a0, NCHK, EPS_REF, BUF, MARGIN, patch=PATCH, check=(kind == "mix"))
        sn = r["snaps"]
        tr = np.load(os.path.join(OUT, f"ciclo_p3_{name}_truth.npy"))
        fid = [float(abs(np.vdot(tr[i], sn[i])) ** 2) for i in range(len(tr))]
        pt = np.abs(tr) ** 2; pm = np.abs(sn) ** 2
        Eb = lambda psi: float(np.real(np.vdot(psi, Hs @ psi)))
        res.update(fid=fid, fid_min=min(fid), P_out=pm[:, far].sum(1).tolist(), P_out_truth=pt[:, far].sum(1).tolist(),
                   E_base_final=Eb(sn[-1]), E_base_final_truth=Eb(tr[-1]), hist=r["hist"], events=r["events"],
                   refined=r["refined"], cost_steps=float(r["cost_steps"]), cost_check=float(r["cost_check"]),
                   cost=float(r["cost_steps"] + r["cost_check"]), norm_loss=r["norm_loss"],
                   dens_final=pm[-1].reshape(-1, 3).sum(1).tolist(), dens_final_truth=pt[-1].reshape(-1, 3).sum(1).tolist())
    res["wall"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, f"ciclo_p3_{name}_{kind}.json"), "w"), default=float)
    print(name, kind, "done", round(res["wall"], 1), "s", {k: res.get(k) for k in ("fid_min", "cost", "norm_loss")})


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

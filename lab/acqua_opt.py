"""acqua_opt.py - ottimizza gli esponenti delle basi even-tempered sui nostri atomi (vedi acqua_basis).
Uso: python acqua_opt.py <fase> <X>   fase in {Hhf, Ohf, Opol, Hpol}. Checkpoint in lab/out/acqua_basi.json."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
import acqua_int as I, acqua_qc as Q, acqua_basis as AB

LOGF = os.path.join(AB.OUT, "acqua_opt_log.txt")


def log(s):
    with open(LOGF, "a") as f:
        f.write(s + "\n")


def atom_E(el, X, par, corr=False, aug=False):
    """spherical ROHF in the s,p part (O) or s part (H), higher l appended as virtuals; optional all-electron CCSD"""
    import acqua_run as AR
    o = AR.atom_calc(el, X, par, aug, corr=corr, triples=False)
    return o["E_hf"] + (o.get("Ec_ccsd", 0.0) if corr else 0.0), o


def h2_E(X, par, R=1.4):
    cen, Z, sh = AB.shells_for(["H", "H"], [[0, 0, 0], [0, 0, R]], X, aug=False, bas={f"H{X}": par})
    B = I.Basis(cen, Z, sh); ints = I.integrals(B); r = Q.scf(B, ints, 1)
    h, G = Q.mo_ints(B, ints, r)
    Ec, t1, t2, c = Q.rccsd(G, r["eps"], 1)
    return r["E"] + Ec


if __name__ == "__main__":
    fase, X = sys.argv[1], int(sys.argv[2])
    bas = AB.load(); t0 = time.time()
    if fase == "Hhf":
        f = lambda x: atom_E("H", X, {"hf": [list(np.exp(x))]})[0]
        x0 = np.log([0.05, 3.5]) if X < 4 else np.log([0.04, 3.2])
        rr = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-10, maxfev=400))
        bas = AB.load(); bas[f"H{X}"] = {"hf": [list(np.exp(rr.x))], "E_hf": rr.fun}
    elif fase == "Ohf":
        def f(x):
            p = np.exp(x); return atom_E("O", X, {"hf": [list(p[:2]), list(p[2:])]})[0]
        x0 = np.log([0.15, 3.0, 0.15, 2.8])
        rr = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-9, maxfev=600))
        p = np.exp(rr.x); bas = AB.load(); bas[f"O{X}"] = {"hf": [list(p[:2]), list(p[2:])], "E_hf": rr.fun}
    elif fase in ("Opol", "Hpol"):
        el = fase[0]; par = dict(bas[f"{el}{X}"]); npol = len(AB.NFUN[X][el]) - AB.NHF[el]
        z0 = {"O": [1.4, 2.4, 2.8], "H": [1.0, 1.6, 2.0]}[el][:npol]

        single = max(AB.NFUN[X][el][AB.NHF[el]:]) == 1     # one function per l: the ratio does not enter

        def unpack(x):
            q = np.exp(x); return {"hf": par["hf"], "pol": [list(q[:npol]), 2.0 if single else float(1.0 + q[npol])]}

        cf = os.path.join(AB.OUT, f"acqua_opt_cache_{fase}{X}.json")
        cache = json.load(open(cf)) if os.path.exists(cf) else {}

        def f(x):
            kk = ",".join(f"{v:.10f}" for v in x)
            if kk in cache:
                return cache[kk]
            e = atom_E("O", X, unpack(x), corr=True)[0] if el == "O" else h2_E(X, unpack(x))
            cache[kk] = e; json.dump(cache, open(cf, "w"))
            log(f"{fase} X={X} {list(np.round(np.exp(x), 5))} E={e:.10f} t={time.time() - t0:.0f}")
            return e
        x0 = np.log(z0 if single else z0 + [1.5])
        simplex = np.vstack([x0] + [x0 + 0.3 * np.eye(len(x0))[k] for k in range(len(x0))])
        rr = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=2e-3, fatol=2e-6, maxfev=60 if el == "O" else 200,
                                                                initial_simplex=simplex))
        p = unpack(rr.x); par.update(p); par["E_corr_opt"] = rr.fun; bas = AB.load(); bas[f"{el}{X}"] = par
    AB.save(bas)
    log(f"DONE {fase} X={X} {json.dumps(bas[f'{fase[0]}{X}'])} t={time.time() - t0:.0f}")
    print("done", fase, X, rr.fun, rr.nfev)

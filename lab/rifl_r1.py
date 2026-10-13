"""rifl_r1.py - esami E0 (NR -> c_bm), E1 (Dirac continuo vs Ashfaque et al. 2019), E4 (nucleo = regola vera).
uso: python lab/rifl_r1.py [E0|E1|E4]   (ogni parte < 10 min; risultati in lab/out/rifl_r1.json)"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rifl_core import *

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rifl_r1.json")
res = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    json.dump(res, open(OUT, "w"), indent=1)


def fit_inv(P, lam, nterm):
    """lam = l_inf - sum_j A_j / P^j, j=1..nterm-1 (minimi quadrati)"""
    P = np.asarray(P, float)
    X = np.column_stack([np.ones_like(P)] + [-(1 / P) ** j for j in range(1, nterm)])
    c = np.linalg.lstsq(X, np.asarray(lam), rcond=None)[0]
    return c[0]


part = sys.argv[1] if len(sys.argv) > 1 else "E0"
t0 = time.time()
if part == "E0":
    Us = [10, 14.142, 20, 28.284, 40, 56.569, 80]
    lam = []
    for U in Us:
        lam.append(nr_backflow(U, cosc=6.0))
        print("U", U, lam[-1], "%.0fs" % (time.time() - t0), flush=True)
    # controllo risoluzione al U piu' grande che costa poco
    chk = nr_backflow(40, cosc=3.0)
    ext = {f"fit{n}_last{m}": fit_inv(Us[-m:], lam[-m:], n) for n in (2, 3, 4) for m in (4, 5, 6) if m > n}
    best = ext["fit4_last6"]
    spread = float(np.std(list(ext.values())))
    res["E0"] = {"U": Us, "lambda": lam, "risoluzione_U40_cosc3_vs_6": chk - lam[4], "estrapolazioni": ext,
                 "c_estrapolato": best, "dispersione_fit": spread, "c_bm_ref": C_BM,
                 "scarto": best - C_BM, "PASS": bool(abs(best - C_BM) <= 2e-6)}
elif part == "E1":
    # uso: E1 <eps>  (una eps per chiamata, < 10 min); poi E1sum
    tab = {0.10: 0.03686, 0.50: 0.03088, 0.80: 0.02722, 1.00: 0.02498, 1.60: 0.01947, 2.00: 0.01660, 2.50: 0.01372}
    eps = float(sys.argv[2])
    ref = tab[round(eps, 2)]
    T = 4 / eps**2
    COSC = 4.0                                # cosc = 12 era sotto-risolto (prova dichiarata nel documento)
    p_, _ = dirac_grid(T, 1000.0, COSC)
    Pmax = min(4000.0, 1000.0 * 6000.0 / len(p_))   # M <= ~6000 (memoria)
    Ps = Pmax / 2.0 ** np.arange(5)[::-1]
    lam = []
    for P in Ps:
        lam.append(dirac_backflow(T, P, cosc=COSC))
        print(eps, P, lam[-1], "%.0fs" % (time.time() - t0), flush=True)
    chk = dirac_backflow(T, Ps[-3], cosc=COSC / 2)
    ext = {f"fit{n}": fit_inv(Ps[-(n + 1):], lam[-(n + 1):], n) for n in (2, 3, 4)}
    best = ext["fit3"]
    unc = max(abs(ext["fit2"] - best), abs(ext["fit4"] - best))
    res.setdefault("E1_eps", {})[f"{eps:.2f}"] = {
        "T": T, "P": list(Ps), "lambda": lam, "cosc": COSC, "risoluzione_P[-3]_cosc/2": chk - lam[-3],
        "estrapolazioni": ext, "c_rbf": best, "incertezza": unc, "tabella": ref,
        "scarto": best - ref, "entro_1e-3": bool(abs(best - ref) <= 1e-3), "entro_5e-5": bool(abs(best - ref) <= 5e-5)}
    print(eps, best, "+-", unc, "tab", ref, flush=True)
    part = f"E1_{eps}"
elif part == "E1sum":
    out = res["E1_eps"]
    res["E1"] = {"per_eps": out, "PASS": all(v["entro_1e-3"] for v in out.values()),
                 "stretto": all(v["entro_5e-5"] for v in out.values())}
    res.pop("E1_vecchio_cosc12", None)
elif part == "E4":
    out = []
    for theta, N, Lr, kmax in ((0.3, 40, 400, np.pi / 2), (0.05, 200, 1200, np.pi / 2), (1.0, 15, 200, np.pi)):
        j = np.arange(1, int(Lr * kmax / (2 * np.pi)))
        k = 2 * np.pi * j / Lr
        k = k[k < kmax]
        A = walk_kernel(k, np.full(len(k), 2 * np.pi / Lr), theta, N)
        lam, v = top_eig(A, True)
        d = lattice_check(theta, N, v, k, Lr)
        out.append({"theta": theta, "N": N, "anello": Lr, "kmax": kmax, "nucleo": lam, "regola": d, "diff": lam - d})
        print(out[-1], flush=True)
    res["E4"] = {"casi": out, "PASS": all(abs(o["diff"]) <= 1e-12 for o in out)}
res.setdefault("tempi", {})[part] = time.time() - t0
save()
print(part, "fatto in %.0f s" % (time.time() - t0))

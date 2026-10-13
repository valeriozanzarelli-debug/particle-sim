"""rifl_r3.py - domande nuove.
 Q1 <eps>  : Dirac continuo c_rbf(eps) per eps grandi (regime ultrarelativistico), taglio estrapolato
 Q1sum     : fit della legge di potenza
 Q2 <theta>: la regola a theta finito: Lambda(theta, N) per N = 1..2000, finestra fisica (k<pi/2) e banda intera (k<pi)
 Q2sum     : massimi
uso: python lab/rifl_r3.py Q1 10   (ogni pezzo < 10 min); risultati in lab/out/rifl_r3.json"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rifl_core import *

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rifl_r3.json")
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
NS = [1, 2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 32, 40, 50, 64, 80, 100, 128, 160, 200, 256, 320, 400, 512, 640,
      800, 1000, 1280, 1600, 2000]


def save():
    json.dump(res, open(OUT, "w"), indent=1)


def fit_inv(P, lam, nterm):
    P = np.asarray(P, float)
    X = np.column_stack([np.ones_like(P)] + [-(1 / P) ** j for j in range(1, nterm)])
    return float(np.linalg.lstsq(X, np.asarray(lam), rcond=None)[0][0])


t0 = time.time()
part = sys.argv[1]
if part == "Q1":
    # griglia geometrica (geo=True) e taglio P = X/T (prereg2); la prima versione aveva P troppo basso
    eps = float(sys.argv[2])
    T = 4 / eps**2
    COSC = 4.0
    Xs = [50, 100, 200, 400, 800, 1600]
    Ps = [X / T for X in Xs]
    lam = []
    for P in Ps:
        lam.append(dirac_backflow(T, P, cosc=COSC, geo=True))
        print("Q1", eps, P, lam[-1], "%.0fs" % (time.time() - t0), flush=True)
    chk = dirac_backflow(T, Ps[-3], cosc=COSC / 2, geo=True) - lam[-3]
    ext = {f"fit{n}": fit_inv(Ps[-(n + 1):], lam[-(n + 1):], n) for n in (2, 3, 4)}
    best = ext["fit3"]
    unc = max(abs(ext["fit2"] - best), abs(ext["fit4"] - best))
    res.setdefault("Q1_eps", {})[f"{eps:g}"] = {"T": T, "P": list(Ps), "lambda": lam, "risoluzione": chk,
                                              "estrapolazioni": ext, "c_rbf": best, "incertezza": unc,
                                              "c_eps2": best * eps**2, "delta_1_su_pi": 1 - np.pi * best * eps**2}
    print("Q1", eps, best, "+-", unc, "c eps^2 =", best * eps**2, flush=True)
    part = f"Q1_{eps:g}"
elif part == "Q1sum":
    d = res["Q1_eps"]
    e = np.array(sorted(float(x) for x in d))
    c = np.array([d[f"{x:g}"]["c_rbf"] for x in e])
    fits = {}
    for lo, hi in ((20, 50), (50, 100), (20, 100)):
        m = (e >= lo) & (e <= hi)
        g, lA = np.polyfit(np.log(e[m]), np.log(c[m]), 1)
        fits[f"{lo}-{hi}"] = {"gamma": -g, "A": float(np.exp(lA))}
    # fit della letteratura (Ashfaque et al. 2019, eq. 19) per confronto
    alpha = 7.2973525693e-3
    lit = {f"{x:g}": C_BM * np.exp(-(4 * x / 9) * (1 - 4 * alpha * x)) for x in e}
    # c * eps^2 (= c * 4/(mT)): se gamma = 2 tende a una costante
    dl = {f"{x:g}": float(1 - np.pi * y * x**2) for x, y in zip(e, c)}
    big = [x for x in e if x >= 200]
    dd = [dl[f"{x:g}"] for x in big]
    q1b = {"delta": dl, "monotono": bool(all(abs(dd[i + 1]) < abs(dd[i]) for i in range(len(dd) - 1))),
           "delta_10000": dl.get("10000"), "PASS": None}
    if q1b["delta_10000"] is not None:
        q1b["PASS"] = bool(q1b["monotono"] and abs(q1b["delta_10000"]) <= 0.01)
    res["Q1b"] = q1b
    res["Q1"] = {"eps": list(e), "c_rbf": list(c), "fit_potenza": fits,
                 "gamma_stabile_0.05": bool(abs(fits["20-50"]["gamma"] - fits["50-100"]["gamma"]) <= 0.05),
                 "c_eps2": list(c * e**2), "fit_letteratura": lit, "per_eps": d}
    print(json.dumps(fits, indent=1))
elif part == "Q2":
    th = float(sys.argv[2])
    rows = res.setdefault("Q2_theta", {}).setdefault(f"{th:g}", {})
    for N in NS:
        if str(N) in rows:
            continue
        if time.time() - t0 > 480:
            print("tempo: riprendere", flush=True)
            break
        lp = walk_backflow_gram(th, N, np.pi / 2)
        lf = walk_backflow_gram(th, N, np.pi)
        rows[str(N)] = {"fisica": lp, "intera": lf, "eps": float(eps_of(th, N))}
        save()
        print("Q2", th, N, lp, lf, "%.0fs" % (time.time() - t0), flush=True)
    part = f"Q2_{th:g}"
elif part == "Q2long":
    # N piu' grandi col metodo del nucleo (matrice M x M, M ~ v_g N): solo theta grandi dove M resta <= ~6000
    th = float(sys.argv[2])
    win = sys.argv[3]                     # "fisica" o "intera"
    kmax = np.pi / 2 if win == "fisica" else np.pi
    rows = res.setdefault("Q2long_theta", {}).setdefault(f"{th:g}", {}).setdefault(win, {})
    for N in [2000, 2560, 3200, 4000, 5000, 6400, 8000, 10000, 12800, 16000, 20000]:
        if str(N) in rows:
            continue
        k, wq = walk_grid(th, N, kmax, 6.0)
        if len(k) > 6000:
            print("M troppo grande", N, len(k), flush=True)
            break
        if time.time() - t0 > 360:
            print("tempo: riprendere", flush=True)
            break
        rows[str(N)] = {"lambda": walk_backflow(th, N, kmax, cosc=6.0), "M": len(k), "eps": float(eps_of(th, N))}
        save()
        print("Q2long", th, win, N, len(k), rows[str(N)]["lambda"], "%.0fs" % (time.time() - t0), flush=True)
    part = f"Q2long_{th:g}_{win}"
elif part == "Q2redo":
    # emendamento A9: rifare N <= 400 con la griglia corretta; controllo a N = 1000 (valore vecchio conservato)
    th = float(sys.argv[2])
    rows = res["Q2_theta"][f"{th:g}"]
    chk = res.setdefault("Q2redo_controllo", {})
    for N in [n for n in NS if n <= 400] + [1000]:
        lp = walk_backflow_gram(th, N, np.pi / 2)
        lf = walk_backflow_gram(th, N, np.pi)
        old = rows.get(str(N), {})
        chk[f"{th:g}_{N}"] = {"fisica_vecchio": old.get("fisica"), "fisica_nuovo": lp,
                              "intera_vecchio": old.get("intera"), "intera_nuovo": lf}
        rows[str(N)] = {"fisica": lp, "intera": lf, "eps": float(eps_of(th, N))}
    save()
    d = [abs(v["fisica_nuovo"] - v["fisica_vecchio"]) for k_, v in chk.items() if k_.startswith(f"{th:g}_")]
    print("Q2redo", th, "max |nuovo - vecchio| =", max(d), flush=True)
    part = f"Q2redo_{th:g}"
elif part == "Q2osc":
    # esplorazione (emendamento A10): banda intera a theta grande oscilla con N (periodo ~2 pi/(pi - 2 theta) tick,
    # interferenza tra i due bordi della banda k~0 e k~pi). Scansione fitta di N per trovare il massimo vero.
    th = float(sys.argv[2])
    N0 = int(sys.argv[3])
    step = int(sys.argv[4])
    per = 2 * np.pi / (np.pi - 2 * th)
    rows = res.setdefault("Q2osc", {}).setdefault(f"{th:g}", {}).setdefault(str(N0), {})
    for N in range(N0, N0 + int(np.ceil(per)) + step, step):
        if str(N) in rows:
            continue
        if time.time() - t0 > 450:
            print("tempo: riprendere", flush=True)
            break
        rows[str(N)] = walk_backflow(th, N, np.pi, cosc=6.0)
        save()
        print("Q2osc", th, N, rows[str(N)], "%.0fs" % (time.time() - t0), flush=True)
    part = f"Q2osc_{th:g}_{N0}"
elif part == "Q2sum":
    d = res["Q2_theta"]
    best_p = max(((float(t), int(N), v["fisica"]) for t, r in d.items() for N, v in r.items()), key=lambda z: z[2])
    best_f = max(((float(t), int(N), v["intera"]) for t, r in d.items() for N, v in r.items()), key=lambda z: z[2])
    res["Q2"] = {"max_fisica": best_p, "max_intera": best_f, "c_bm": C_BM, "tight_binding": 0.0764734,
                 "fisica_supera_c_bm": bool(best_p[2] > C_BM), "intera_supera_c_bm": bool(best_f[2] > C_BM)}
    print(res["Q2"])
res.setdefault("tempi", {})[part] = time.time() - t0
save()
print(part, "fatto in %.0f s" % (time.time() - t0))

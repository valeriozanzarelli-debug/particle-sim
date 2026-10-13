"""rifl_r2.py - esami E2 (limite NR della regola -> c_bm) ed E3 (regola -> Dirac continuo per a -> 0).
uso: python lab/rifl_r2.py E2 <U>  |  E2sum  |  E3 <eps>  |  E3sum     (ogni pezzo < 10 min)
Emendamento A1 (dichiarato in lab/out/rifl_emendamenti.json): in E3 la sequenza di a e' theta = 4/(eps^2 N) con
N intero che raddoppia (theta <= 0.4), perche' per eps >= 2.5 i theta della pre-registrazione danno N non interi."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rifl_core import *

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rifl_r2.json")
res = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    json.dump(res, open(OUT, "w"), indent=1)


def polyfit0(x, y, deg):
    return float(np.polyfit(np.asarray(x, float), np.asarray(y, float), deg)[-1])


def fit_inv(P, lam, nterm):
    P = np.asarray(P, float)
    X = np.column_stack([np.ones_like(P)] + [-(1 / P) ** j for j in range(1, nterm)])
    return float(np.linalg.lstsq(X, np.asarray(lam), rcond=None)[0][0])


t0 = time.time()
part = sys.argv[1]
if part == "E2":
    U = float(sys.argv[2])
    pm = 0.05                       # p_max / m
    T = 2 * U**2 / pm**2            # U = p_max sqrt(T/2m), m = 1
    thetas = [0.2, 0.1, 0.05, 0.025]
    lam = []
    for th in thetas:
        N = int(round(T / th))
        lam.append(walk_backflow(th, N, pm * th, cosc=6.0))
        print("E2", U, th, N, lam[-1], "%.0fs" % (time.time() - t0), flush=True)
    chk = walk_backflow(thetas[0], int(round(T / thetas[0])), pm * thetas[0], cosc=4.0) - lam[0]
    nr = nr_backflow(U, cosc=6.0)
    dc = dirac_backflow(T, pm, cosc=6.0)
    a0 = polyfit0(thetas, lam, 2)
    res.setdefault("E2_U", {})[f"{U:g}"] = {"T": T, "theta": thetas, "lambda": lam, "risoluzione": chk,
                                            "regola_a0": a0, "NR": nr, "Dirac_stesso_taglio": dc,
                                            "scarto_rel_NR": (a0 - nr) / nr, "scarto_rel_Dirac": (a0 - dc) / dc}
    print("E2", U, "a->0", a0, "NR", nr, "Dirac", dc, flush=True)
    part = f"E2_{U:g}"
elif part == "E2sum":
    d = res["E2_U"]
    Us = sorted(float(u) for u in d)
    a0 = [d[f"{u:g}"]["regola_a0"] for u in Us]
    ext = {f"fit{n}": fit_inv(Us[-(n + 1):], a0[-(n + 1):], n) for n in (2, 3)}
    best = ext["fit3"]
    res["E2"] = {"U": Us, "regola_a0": a0, "estrapolazioni_U": ext, "c_regola": best, "c_bm": C_BM,
                 "scarto_rel": (best - C_BM) / C_BM,
                 "PASS": bool(abs(best - C_BM) / C_BM <= 0.01 and all(abs(d[f"{u:g}"]["scarto_rel_NR"]) <= 0.01 for u in Us))}
    print(res["E2"])
elif part == "E3":
    eps = float(sys.argv[2])
    T = 4 / eps**2
    N0 = int(np.ceil(T / 0.4))
    Ns = [N0 * 2**j for j in range(5)]
    rows = []
    for N in Ns:
        th = T / N
        if th > 0.4 + 1e-12:
            continue
        lp = walk_backflow(th, N, np.pi / 2, cosc=6.0)
        lf = walk_backflow(th, N, np.pi, cosc=6.0)
        P = np.pi / (2 * th)
        dP = dirac_backflow(T, P, cosc=4.0) if 3 * T * P < 9000 else None
        rows.append({"N": N, "theta": th, "fisica": lp, "intera": lf, "Dirac_stesso_taglio": dP})
        print("E3", eps, rows[-1], "%.0fs" % (time.time() - t0), flush=True)
        if time.time() - t0 > 420:
            break
    r = rows[-1]
    chk = walk_backflow(r["theta"], r["N"], np.pi / 2, cosc=3.0) - r["fisica"]
    th = [r["theta"] for r in rows]
    res.setdefault("E3_eps", {})[f"{eps:g}"] = {
        "T": T, "righe": rows, "risoluzione_ultima": chk,
        "fisica_a0_lin": polyfit0(th[-3:], [r["fisica"] for r in rows[-3:]], 1),
        "fisica_a0_quad": polyfit0(th[-4:], [r["fisica"] for r in rows[-4:]], 2) if len(rows) >= 4 else None,
        "intera_a0_lin": polyfit0(th[-3:], [r["intera"] for r in rows[-3:]], 1)}
    part = f"E3_{eps:g}"
elif part == "E3sum":
    r1 = json.load(open(os.path.join(os.path.dirname(OUT), "rifl_r1.json")))
    q1 = json.load(open(os.path.join(os.path.dirname(OUT), "rifl_r3.json"))).get("Q1", {}) if os.path.exists(
        os.path.join(os.path.dirname(OUT), "rifl_r3.json")) else {}
    out = {}
    for e, v in res["E3_eps"].items():
        ref = None
        key = f"{float(e):.2f}"
        if key in r1.get("E1_eps", {}):
            ref = r1["E1_eps"][key]["c_rbf"]
        elif "per_eps" in q1 and e in q1["per_eps"]:
            ref = q1["per_eps"][e]["c_rbf"]
        a0 = v["fisica_a0_quad"] if v["fisica_a0_quad"] is not None else v["fisica_a0_lin"]
        out[e] = {"regola_a0": a0, "Dirac": ref, "scarto_rel": None if ref is None else (a0 - ref) / ref,
                  "PASS": None if ref is None else bool(abs(a0 - ref) / ref <= 0.01)}
    res["E3"] = {"per_eps": out, "PASS": all(o["PASS"] for o in out.values() if o["PASS"] is not None)}
    print(res["E3"])
res.setdefault("tempi", {})[part] = time.time() - t0
save()
print(part, "fatto in %.0f s" % (time.time() - t0))

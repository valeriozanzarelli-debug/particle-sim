"""acqua_plan.py - pianificatore della produzione H2O (studio ACQUA).
  python acqua_plan.py hfopt      minimo HF nella base X=2 (centro di partenza della griglia, nessun dato esterno)
  python acqua_plan.py next       scrive in out/acqua_coda.txt i punti mancanti della fase corrente, stampa la fase
  python acqua_plan.py fit        stampa i minimi/le superfici correnti
Superficie: griglia 3x3 simmetrica (r1 = r2 = r, theta), passi 0.01 A e 1.5 gradi, interpolante biquadratica.
Superficie CBS = E_HF(X=3) + E_c(2,3) estrapolata X^-3 punto per punto (E_c = CCSD + (T), tutti gli elettroni).
Ricentra se il minimo dista piu' di 0.004 A o 0.6 gradi dal centro (al massimo 3 volte per base)."""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_run as AR, acqua_basis as AB

DR, DT = 0.01, 1.5
SF = os.path.join(AB.OUT, "acqua_stato.json")
CODA = os.path.join(AB.OUT, "acqua_coda.txt")
FIELDS = [0.002, -0.002, 0.004, -0.004]
DELTAS = [0.01, 0.02]


def st_load():
    return json.load(open(SF)) if os.path.exists(SF) else {}


def st_save(s):
    json.dump(s, open(SF, "w"), indent=1)


def pt(X, r1, r2, th, F=0.0):
    k = AR.key("mol", X, float(r1), float(r2), float(th), float(F), 1)
    f = os.path.join(AR.PD, k + ".json")
    return json.load(open(f)) if os.path.exists(f) else None


def line(X, r1, r2, th, F=0.0):
    return f"mol {X} {r1:.6f} {r2:.6f} {th:.6f} {F}"


def grid(c):
    r0, t0 = c
    return [(round(r0 + i * DR, 6), round(t0 + j * DT, 6)) for i in (-1, 0, 1) for j in (-1, 0, 1)]


def E_of(p, kind):
    if kind == "hf":
        return p["E_hf"]
    return p["E_tot"]


def surface(c, kind):
    """kind: 'hf2', 2, 3 or 'cbs'; returns dict(min=(r,t), E, hess (r [A], t [deg]), coeffs) or None"""
    xs = []; Es = []
    for (r, t) in grid(c):
        if kind == "cbs":
            p2, p3 = pt(2, r, r, t), pt(3, r, r, t)
            if p2 is None or p3 is None:
                return None
            ec = lambda p: p["Ec_ccsd"] + p["Ec_T"]
            E = p3["E_hf"] + (27 * ec(p3) - 8 * ec(p2)) / 19
        else:
            p = pt(kind, r, r, t)
            if p is None:
                return None
            E = p["E_tot"]
        xs.append(((r - c[0]) / DR, (t - c[1]) / DT)); Es.append(E)
    A = np.array([[1, x, y, x * x, x * y, y * y, x * x * y, x * y * y, x * x * y * y] for x, y in xs])
    cf = np.linalg.solve(A, np.array(Es))
    f = lambda x, y: cf @ np.array([1, x, y, x * x, x * y, y * y, x * x * y, x * y * y, x * x * y * y])
    g = lambda x, y: np.array([cf[1] + 2 * cf[3] * x + cf[4] * y + 2 * cf[6] * x * y + cf[7] * y * y + 2 * cf[8] * x * y * y,
                               cf[2] + cf[4] * x + 2 * cf[5] * y + cf[6] * x * x + 2 * cf[7] * x * y + 2 * cf[8] * x * x * y])
    Hh = lambda x, y: np.array([[2 * cf[3] + 2 * cf[6] * y + 2 * cf[8] * y * y, cf[4] + 2 * cf[6] * x + 2 * cf[7] * y + 4 * cf[8] * x * y],
                                [cf[4] + 2 * cf[6] * x + 2 * cf[7] * y + 4 * cf[8] * x * y, 2 * cf[5] + 2 * cf[7] * x + 2 * cf[8] * x * x]])
    z = np.zeros(2)
    for _ in range(50):
        z = z - np.linalg.solve(Hh(*z), g(*z))
    H = Hh(*z) / np.array([[DR * DR, DR * DT], [DR * DT, DT * DT]])
    return dict(min=(c[0] + z[0] * DR, c[1] + z[1] * DT), E=float(f(*z)), hess=H.tolist(), off=z.tolist())


def hfopt():
    from scipy.optimize import minimize
    import acqua_int as I, acqua_qc as Q

    def E(x):
        cen, Z, sh = AB.shells_for(["O", "H", "H"], AR.geom(x[0], x[0], x[1]), 2, aug=True)
        B = I.Basis(cen, Z, sh); return Q.scf(B, I.integrals(B), 5)["E"]
    r = minimize(E, [1.0, 110.0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-9, initial_simplex=[[1.0, 110], [1.05, 110], [1.0, 115]]))
    s = st_load(); s["hf2_min"] = [float(r.x[0]), float(r.x[1]), float(r.fun)]; s["c2"] = [round(r.x[0], 3), round(r.x[1] * 20) / 20]
    s["fase"] = "griglia2"; s["ric2"] = 0; s["ric3"] = 0; st_save(s); print(s)


def far(z):
    return abs(z[0]) * DR > 0.004 or abs(z[1]) * DT > 0.6


def nxt():
    s = st_load(); q = []
    f = s["fase"]
    if f == "griglia2":
        c = s["c2"]; q = [line(2, r, r, t) for r, t in grid(c) if pt(2, r, r, t) is None]
        if not q:
            S = surface(c, 2)
            if far(S["off"]) and s["ric2"] < 3:
                s["c2"] = [round(S["min"][0], 3), round(S["min"][1] * 20) / 20]; s["ric2"] += 1; st_save(s); return nxt()
            s["min2"] = S["min"]; s["c3"] = s["c2"]; s["fase"] = "griglia3"; st_save(s); return nxt()
    elif f == "griglia3":
        c = s["c3"]
        q = [line(X, r, r, t) for X in (3, 2) for r, t in grid(c) if pt(X, r, r, t) is None]
        if not q:
            S = surface(c, "cbs")
            if far(S["off"]) and s["ric3"] < 3:
                s["c3"] = [round(S["min"][0], 3), round(S["min"][1] * 20) / 20]; s["ric3"] += 1; st_save(s); return nxt()
            s["min_cbs"] = [round(S["min"][0], 4), round(S["min"][1], 2)]; s["fase"] = "finale"; st_save(s); return nxt()
    elif f == "finale":
        r, t = s["min_cbs"]
        for X in (3, 2):
            q += [line(X, r, r, t)] if pt(X, r, r, t) is None else []
            q += [line(X, r + d, r - d, t) for d in DELTAS if pt(X, r + d, r - d, t) is None]
            q += [line(X, r, r, t, F) for F in FIELDS if pt(X, r, r, t, F) is None]
        for X in (2, 3):
            for el in ("O", "H"):
                if not os.path.exists(os.path.join(AR.PD, AR.key("atom", el, X, 1) + ".json")):
                    q.append(f"atom {el} {X}")
        if not q:
            s["fase"] = "fatto"; st_save(s)
    open(CODA, "w").write("\n".join(q) + ("\n" if q else ""))
    print(s["fase"], len(q))


if __name__ == "__main__":
    {"hfopt": hfopt, "next": nxt, "fit": lambda: print(st_load())}[sys.argv[1]]()

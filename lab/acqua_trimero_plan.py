"""acqua_trimero_plan.py - pianificatore del trimero (H2O)3 (pre-registrazione lab/out/acqua_trimero_prereg.json).
Fase A  griglia2: E_int^CP MP2 a X=2 su griglia 3x3x3 in (R, td, ps), passi 0.05 A, 5 gradi, 15 gradi, triquadratica,
        ricentra (max 6) come nel dimero. Centro di partenza: R = 2.92 A (il nostro dimero), td = 0, ps = 90.
Fase B  scanR: angoli fermi; E_int^CP MP2/CBS(2,3) in R = R2 - 0.12, -0.08, -0.04, 0 (R2 = minimo X=2); parabola sui
        3 punti attorno al piu' basso; se il piu' basso e' al bordo si aggiunge un punto a 0.04 oltre (max 3 volte).
Fase C  finale: al minimo, tutti e 7 i tipi: MP2 X=2,3 + CCSD(T) X=2 -> energia "focal point"
        E = HF(3) + MP2_cbs(2,3) + [CCSD(T) - MP2](X=2).
  python acqua_trimero_plan.py next | judge"""
import os, sys, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_trimero as T, acqua_dimero as DM, acqua_basis as AB

SF = os.path.join(AB.OUT, "acqua_trimero_stato.json"); CODA = os.path.join(AB.OUT, "acqua_tcoda.txt")
STEP = np.array([0.05, 5.0, 15.0]); TAU = 1e-6; EH_EV = 27.211386245988; KCAL = 627.5094740631
KINDS = ("ABC", "AB", "AC", "BC", "A", "B", "C")


def lm(kind, X, *g):
    f = os.path.join(T.PD, DM.key(kind + "m", X, *g) + ".json")
    return json.load(open(f)) if os.path.exists(f) else None


def lc(kind, X, *g):
    f = os.path.join(T.PD, DM.key(kind, X, *g, f"{TAU:.0e}") + ".json")
    return json.load(open(f)) if os.path.exists(f) else None


def comb(get, w):
    ps = {k: get(k) for k in w}
    return None if any(p is None for p in ps.values()) else sum(w[k] * ps[k] for k in w)


WINT = dict(ABC=1, A=-1, B=-1, C=-1)
W3 = dict(ABC=1, AB=-1, AC=-1, BC=-1, A=1, B=1, C=1)
W2 = {"AB": 1, "AC": 1, "BC": 1, "A": -2, "B": -2, "C": -2}


def mp2cbs(g, w):
    hf3 = comb(lambda k: (lm(k, 3, *g) or {}).get("E_hf"), w)
    c2 = comb(lambda k: (lm(k, 2, *g) or {}).get("E_mp2"), w); c3 = comb(lambda k: (lm(k, 3, *g) or {}).get("E_mp2"), w)
    return None if None in (hf3, c2, c3) else hf3 + (27 * c3 - 8 * c2) / 19


def focal(g, w):
    m = mp2cbs(g, w)
    d = comb(lambda k: None if lc(k, 2, *g) is None else lc(k, 2, *g)["E_tot"] - lc(k, 2, *g)["E_hf"] - lc(k, 2, *g)["E_mp2_full"], w)
    return None if None in (m, d) else m + d


def grid(c):
    return [tuple(round(c[k] + d[k] * STEP[k], 4) for k in range(3)) for d in itertools.product((-1, 0, 1), repeat=3)]


def fit(c):
    pts = grid(c); E = [comb(lambda k: (lm(k, 2, *g) or {}).get("E_tot"), WINT) for g in pts]
    if any(e is None for e in E):
        return None
    xs = [tuple((g[k] - c[k]) / STEP[k] for k in range(3)) for g in pts]
    phi = lambda x: np.array([x[0] ** a * x[1] ** b * x[2] ** cc for a in range(3) for b in range(3) for cc in range(3)])
    cf = np.linalg.solve(np.array([phi(x) for x in xs]), np.array(E)); f = lambda x: cf @ phi(x)
    z = np.zeros(3); h = 1e-4
    for _ in range(60):
        gr = np.array([(f(z + h * e) - f(z - h * e)) / (2 * h) for e in np.eye(3)])
        H = np.array([[(f(z + h * (e1 + e2)) - f(z + h * (e1 - e2)) - f(z - h * (e1 - e2)) + f(z - h * (e1 + e2))) / (4 * h * h)
                       for e2 in np.eye(3)] for e1 in np.eye(3)])
        w, U = np.linalg.eigh(H); st = -U @ ((U.T @ gr) / np.maximum(np.abs(w), 1e-6))
        z = np.clip(z + np.clip(st, -0.5, 0.5), -3, 3)
    return dict(off=z.tolist(), min=(np.array(c) + z * STEP).tolist(), E=float(f(z)), Hpd=bool(np.all(w > 0)))


def nxt():
    s = json.load(open(SF)) if os.path.exists(SF) else dict(fase="griglia2", c=[2.92, 0.0, 90.0], ric=0, storia=[])
    q = []
    if s["fase"] == "griglia2":
        q = [f"mp2 {k} 2 {g[0]} {g[1]} {g[2]}" for g in grid(s["c"]) for k in ("ABC", "A", "B", "C") if lm(k, 2, *g) is None]
        if not q:
            F = fit(s["c"]); s["storia"].append(dict(c=s["c"], fit=F)); z = np.array(F["off"])
            if (np.abs(z) > 0.5).any() and s["ric"] < 6:
                c = np.array(s["c"]) + np.clip(z, -2, 2) * STEP
                s["c"] = [round(c[0], 3), round(c[1], 1), round(c[2], 1)]; s["ric"] += 1
            else:
                m = F["min"]; s["ang"] = [round(m[1], 1), round(m[2], 1)]; s["R2"] = round(m[0], 3)
                s["Rs"] = [round(s["R2"] + d, 3) for d in (-0.12, -0.08, -0.04, 0.0)]; s["ext"] = 0; s["fase"] = "scanR"
            json.dump(s, open(SF, "w"), indent=1); return nxt()
    elif s["fase"] == "scanR":
        gs = [(R, *s["ang"]) for R in s["Rs"]]
        q = [f"mp2 {k} {X} {g[0]} {g[1]} {g[2]}" for g in gs for X in (3, 2) for k in ("ABC", "A", "B", "C") if lm(k, X, *g) is None]
        if not q:
            Rs = sorted(s["Rs"]); E = [mp2cbs((R, *s["ang"]), WINT) for R in Rs]; i = int(np.argmin(E))
            if (i == 0 or i == len(Rs) - 1) and s["ext"] < 3:
                s["Rs"].append(round(Rs[i] + (-0.04 if i == 0 else 0.04), 3)); s["ext"] += 1
            else:
                i = min(max(i, 1), len(Rs) - 2); a, b, c = np.polyfit(Rs[i - 1:i + 2], E[i - 1:i + 2], 2)
                s["Rmin"] = float(-b / (2 * a)); s["fase"] = "finale"
            json.dump(s, open(SF, "w"), indent=1); return nxt()
    elif s["fase"] == "finale":
        g = (round(s["Rmin"], 3), *s["ang"])
        q = [f"mp2 {k} {X} {g[0]} {g[1]} {g[2]}" for X in (3, 2) for k in KINDS if lm(k, X, *g) is None]
        q += [f"calc {k} 2 {g[0]} {g[1]} {g[2]} {TAU}" for k in KINDS if lc(k, 2, *g) is None]
        if not q:
            s["fase"] = "fatto"
    json.dump(s, open(SF, "w"), indent=1)
    open(CODA, "w").write("\n".join(q) + ("\n" if q else "")); print(s["fase"], len(q))


def judge():
    s = json.load(open(SF)); g = (round(s["Rmin"], 3), *s["ang"])
    Ei, E3, E2 = focal(g, WINT), focal(g, W3), focal(g, W2)
    De = -Ei * EH_EV; e3 = E3 * EH_EV; R = s["Rmin"]
    v = lambda val, ref, thr: dict(valore=val, rif=ref, rel=(val - ref) / ref, esito="PASS" if abs(val - ref) / abs(ref) <= thr else "FAIL")
    out = dict(geometria=g, R_min=R, R_X2=s["R2"], angoli=s["ang"], Eint_focal=Ei, E3corpi=E3, E2corpi_somma=E2,
               Eint_mp2cbs=mp2cbs(g, WINT), D_e_eV=De, D_e_kcal=-Ei * KCAL, E3_kcal=E3 * KCAL,
               verdetti=dict(T1_De=v(De, 0.685, 0.10), T2_ROO=v(R, 2.80, 0.02), T3_3corpi=v(e3, -0.100, 0.25)))
    json.dump(out, open(os.path.join(AB.OUT, "acqua_trimero_verdetti.json"), "w"), indent=1); print(json.dumps(out))


if __name__ == "__main__":
    {"next": nxt, "judge": judge}[sys.argv[1]]()

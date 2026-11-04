"""acqua_dimero_plan.py - pianificatore del dimero (H2O)2 (vedi acqua_dimero.py, pre-registrazione
lab/out/acqua_dimero_prereg.json).
Superficie: E_int^CP(R, td, ta) CCSD(T)+FNO a X=2 su griglia 3x3x3 (passi 0.05 A, 5 gradi, 15 gradi), interpolante
triquadratica (27 coefficienti), minimo di Newton; ricentra (salto massimo 2 passi) se il minimo esce di mezzo passo,
al massimo 6 volte. Centro di partenza neutro: R = 3.0 A, td = 0, ta = 0 (nessun dato esterno).
Poi al minimo: AB, A, B a X=2 e X=3 -> D_e CBS(2,3) con counterpoise.
  python acqua_dimero_plan.py next | judge"""
import os, sys, json, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_dimero as DM, acqua_basis as AB

SF = os.path.join(AB.OUT, "acqua_dimero_stato.json"); CODA = os.path.join(AB.OUT, "acqua_dcoda.txt")
STEP = np.array([0.05, 5.0, 15.0]); TAU = 1e-6; EH_EV = 27.211386245988


def load(*a):
    f = os.path.join(DM.PD, DM.key(*a, f"{TAU:.0e}") + ".json")
    return json.load(open(f)) if os.path.exists(f) else None


def eint(X, g, part="tot"):
    ps = [load(k, X, *g) for k in ("AB", "A", "B")]
    if any(p is None for p in ps):
        return None
    f = {"tot": lambda p: p["E_tot"], "hf": lambda p: p["E_hf"], "corr": lambda p: p["E_tot"] - p["E_hf"]}[part]
    return f(ps[0]) - f(ps[1]) - f(ps[2])


def grid(c):
    return [tuple(round(c[k] + d[k] * STEP[k], 4) for k in range(3)) for d in itertools.product((-1, 0, 1), repeat=3)]


def fit(c):
    pts = grid(c); E = [eint(2, g) for g in pts]
    if any(e is None for e in E):
        return None
    xs = [tuple((g[k] - c[k]) / STEP[k] for k in range(3)) for g in pts]
    phi = lambda x: np.array([x[0] ** a * x[1] ** b * x[2] ** cc for a in range(3) for b in range(3) for cc in range(3)])
    cf = np.linalg.solve(np.array([phi(x) for x in xs]), np.array(E))
    f = lambda x: cf @ phi(x)
    z = np.zeros(3); h = 1e-4
    for _ in range(60):
        g = np.array([(f(z + h * e) - f(z - h * e)) / (2 * h) for e in np.eye(3)])
        H = np.array([[(f(z + h * (e1 + e2)) - f(z + h * (e1 - e2)) - f(z - h * (e1 - e2)) + f(z - h * (e1 + e2))) / (4 * h * h)
                       for e2 in np.eye(3)] for e1 in np.eye(3)])
        w, U = np.linalg.eigh(H); step = -U @ ((U.T @ g) / np.maximum(np.abs(w), 1e-6))
        z = np.clip(z + np.clip(step, -0.5, 0.5), -3, 3)
    return dict(off=z.tolist(), min=(np.array(c) + z * STEP).tolist(), E=float(f(z)), Hpd=bool(np.all(w > 0)))


def nxt():
    s = json.load(open(SF)) if os.path.exists(SF) else dict(fase="griglia2", c=[3.0, 0.0, 0.0], ric=0, storia=[])
    q = []
    if s["fase"] == "griglia2":
        q = [f"calc {k} 2 {g[0]} {g[1]} {g[2]} {TAU}" for g in grid(s["c"]) for k in ("AB", "A", "B") if load(k, 2, *g) is None]
        if not q:
            F = fit(s["c"]); s["storia"].append(dict(c=s["c"], fit=F))
            z = np.array(F["off"])
            if (np.abs(z) > 0.5).any() and s["ric"] < 6:
                zc = np.clip(z, -2, 2); c = np.array(s["c"]) + zc * STEP
                s["c"] = [round(c[0], 3), round(c[1], 1), round(c[2], 1)]; s["ric"] += 1
            else:
                m = F["min"]; s["min"] = [round(m[0], 3), round(m[1], 1), round(m[2], 1)]; s["fase"] = "finale"
            json.dump(s, open(SF, "w"), indent=1); return nxt()
    elif s["fase"] == "finale":
        g = s["min"]
        q = [f"calc {k} {X} {g[0]} {g[1]} {g[2]} {TAU}" for X in (3, 2) for k in ("AB", "A", "B") if load(k, X, *g) is None]
        if not q:
            s["fase"] = "fatto"
    json.dump(s, open(SF, "w"), indent=1)
    open(CODA, "w").write("\n".join(q) + ("\n" if q else "")); print(s["fase"], len(q))


def judge():
    s = json.load(open(SF)); g = tuple(s["min"]); F = s["storia"][-1]["fit"]
    hf3 = eint(3, g, "hf"); c2, c3 = eint(2, g, "corr"), eint(3, g, "corr")
    Ecbs = hf3 + (27 * c3 - 8 * c2) / 19
    De = -Ecbs * EH_EV; ROO = F["min"][0]
    out = dict(minimo_X2_CP=F["min"], geometria_usata=g, Eint_X2=eint(2, g), Eint_X3=eint(3, g), Eint_cbs=Ecbs,
               Eint_hf_X3=hf3, D_e_eV=De, D_e_kcal=-Ecbs * 627.5094740631)
    out["verdetti"] = dict(
        D1_De=dict(valore=De, rif=0.217, rel=(De - 0.217) / 0.217, esito="PASS" if abs(De - 0.217) / 0.217 <= 0.10 else "FAIL"),
        D2_ROO=dict(valore=ROO, rif=2.91, rel=(ROO - 2.91) / 2.91, esito="PASS" if abs(ROO - 2.91) / 2.91 <= 0.02 else "FAIL"))
    json.dump(out, open(os.path.join(AB.OUT, "acqua_dimero_verdetti.json"), "w"), indent=1)
    print(json.dumps(out))


if __name__ == "__main__":
    {"next": nxt, "judge": judge}[sys.argv[1]]()

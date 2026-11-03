"""tavg_f3_drive.py - FASE 3: tutte le corse pre-registrate (lab/out/tavg_prereg_f3.json), a spezzoni <= 590 s.
Uso: python tavg_f3_drive.py <g1|g2|g3>   (uno spezzone; stampa TUTTO_FATTO quando il gruppo e' finito)
  g1: controllo di macchina (tavg_pecg2.macchina), poi ECG nel blocco O: 2s2 1S, 2s2p 3P, 2s2p 1P (V1)
  g2: CI di coppia nei blocchi Z = 5..8 (P), poi ioni tipo Li nudi Z = 7, 8 (2S, 2P) per il core-valenza
  g3: verifica a 3 elettroni 2s2 2p 2P nei blocchi C e O (A1)
Risultati: lab/out/tavg_f3_<gruppo>.json (tracciati); basi: lab/out/tavg_f3_base/*.json (non tracciate).
Ciclo non sorvegliato: bash lab/tavg_f3_loop.sh <gruppo>  (log non tracciato in lab/out/tavg_f3_<gruppo>.log)"""
import os, sys, json, time
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tav_ecg as T
import tavg_pecg as P
import tavg_pecg2 as Q

OUT = os.path.join(HERE, "out"); BASE = os.path.join(OUT, "tavg_f3_base")
T_BUDGET = 540.0
T0 = time.time()


def left():
    return T_BUDGET - (time.time() - T0)


def rfile(g):
    return os.path.join(OUT, f"tavg_f3_{g}.json")


def load(g):
    return json.load(open(rfile(g))) if os.path.exists(rfile(g)) else {}


def save(g, d):
    tmp = rfile(g) + ".tmp"
    json.dump(d, open(tmp, "w"), indent=1); os.replace(tmp, rfile(g))


def ecg_task(g, name, mk, Zs, Kmax, tol, nn, Lmax=30.0, seed=7):
    """grow a SafeSVM basis (checkpointed) until K = Kmax or |E(K) - E(K-40)| < tol (K >= 120)"""
    res = load(g)
    r = res.get(name, dict(storia=[], fatto=False))
    if r["fatto"]:
        return True
    os.makedirs(BASE, exist_ok=True)
    bp = os.path.join(BASE, f"{name}.json")
    cls = type("SV", (Q.SafeSVM,), {"NN": nn})
    sv = cls(mk(), np.random.default_rng([seed, 0]), Zs, Lmax=Lmax, trials=30, polish=(0.2, 0.08), npol=8)
    if os.path.exists(bp):
        Q.load_basis(bp, sv)
    sv.rng = np.random.default_rng([seed, len(sv.A)])
    hist = r["storia"]
    t_last = time.time()
    while left() > 0:
        K = len(sv.A)
        Es = {h[0]: h[1] for h in hist}
        if K >= Kmax or (K >= 120 and K - 40 in Es and abs(Es[K] - Es[K - 40]) < tol):
            r["fatto"] = True
            break
        e = sv.add()
        if len(sv.A) % 20 == 0:
            hist.append([len(sv.A), float(e)])
            print(f"[{name}] K={len(sv.A)} E={e:.9f}", flush=True)
        if time.time() - t_last > 120:
            Q.save_basis(bp, sv); t_last = time.time()
    Q.save_basis(bp, sv)
    Es = {h[0]: h[1] for h in hist}
    K = len(sv.A)
    r.update(storia=hist, K=K, E=float(sv.E), NN_MIN=nn,
             barra=(2 * abs(Es[K] - Es[K - 40]) if (K in Es and K - 40 in Es) else None))
    res = load(g); res[name] = r; save(g, res)
    return r["fatto"]


_F = {}


def field(Z):
    if Z not in _F:
        _F[Z] = P.BlockField(Z)
    return _F[Z]


def g1():
    js = os.path.join(OUT, "tavg_pecg2_macchina.json")
    m = json.load(open(js)) if os.path.exists(js) else {}
    if not all(m.get(k, {}).get("K", 0) >= 240 for k in ("O_noee_2s2", "He_nudo")):
        Q.macchina(Kmax=240, t_budget=left())
        return False
    for st, spin, L in (("s2_1S", "s2", 0), ("sp_3P", "t2", 1), ("sp_1P", "s2", 1)):
        if left() <= 0:
            return False
        ok = ecg_task("g1", f"V_O_{st}", lambda: Q.QBlockSystem(field(8), spin, L=L), 6.0, 300, 1e-6, Q.NN_MIN)
        if not ok:
            return False
    return True


CI_STATES = ["s2_1S", "sp_3P", "sp_1P", "pp_3P", "pp_1D", "pp_1S"]


def g2():
    import tavg_ci2 as C
    res = load("g2")
    for Z in (5, 6, 7, 8):
        for st in CI_STATES:
            key = f"P_Z{Z}_{st}"
            if key in res:
                continue
            if left() < 150:
                return False
            a = C.run_state(Z, st, lmax=8, n_rad=24, R_box=30.0, field=field(Z), log=lambda s: None)
            b = C.run_state(Z, st, lmax=8, n_rad=18, R_box=30.0, field=field(Z), log=lambda s: None)
            a["Ec_n18"] = b["Ec"]
            a["barra"] = abs(a["Ec"] - b["Ec"]) + 0.5 * abs(a["E_inf"] - a["E_l"][-1])
            res = load("g2"); res[key] = a; save("g2", res)
            print(f"[{key}] Ec={a['Ec'] * 1e3:.3f} mEh barra={a['barra'] * 1e3:.3f} ({a['secondi']:.0f} s)", flush=True)
    for Z in (7, 8):
        for lab, L in (("2S", 0), ("2P", 1)):
            if left() <= 0:
                return False
            ok = ecg_task("g2", f"Li_Z{Z}_{lab}", lambda: T.System(float(Z), "p3", L=L, m_kin=1.0), float(Z),
                          250, 2e-5, Q.NN_MIN)
            if not ok:
                return False
    return True


def g3():
    for Z in (6, 8):
        if left() <= 0:
            return False
        f = field(Z)
        ok = ecg_task("g3", f"A3_Z{Z}_2s2_2p", lambda: P.BlockSystem(f, "p3", L=1), f.Zv, 250, 2e-5, 1e-5)
        if not ok:
            return False
    return True


if __name__ == "__main__":
    done = {"g1": g1, "g2": g2, "g3": g3}[sys.argv[1]]()
    print("TUTTO_FATTO" if done else "SPEZZONE_FINITO", flush=True)

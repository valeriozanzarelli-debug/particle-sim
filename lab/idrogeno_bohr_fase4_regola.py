"""Idrogeno di Bohr, fase 4 (docs/IDROGENO_BOHR.md §13): la massa del MAESTRO sistemata dall'agglomeratore.

Il maestro (il cammino zig-zag fine da cui l'agglomeratore impara la regola grossa) ha la cella theta = a0/4096 e quindi la
massa cinetica m = tan(theta)/theta = 1 + theta^2/3 = 1.000373 (non 1). La regola grossa la eredita (2p +0.0338 %).
Qui NON si corregge niente a mano: l'agglomeratore riceve un controllo in piu', "il maestro e' abbastanza fine?".
Si rifa' la catena con maestri via via piu' fini, theta_j = theta_A / 2^j (stessa regola zig-zag, stessa cella grossa,
b = 256 * 2^j tick per passo grosso), con ESATTAMENTE il learner della fase 3 (lab/idro_regola_P4.py: Taylor c2, c4, c6
fissati al livello sotto, parita' esatta, w = il piu' piccolo fra 4 e 5 che passa kcheck; a0/8 imparato dall'a0/16 con b = 2),
e ci si ferma al primo j in cui la massa della regola a0/16 cambia meno del budget pre-registrato rispetto al j-1.

Uso: python3 lab/idrogeno_bohr_fase4_regola.py J        (una corsa per maestro, <= 10 min; risultati
         lab/out/idrogeno_bohr_fase4_regola_j<J>.json)
     python3 lab/idrogeno_bohr_fase4_regola.py scegli   (applica il criterio, scrive lab/out/idrogeno_bohr_fase4_regola.json
         nel formato degli asset agglomerato_zz_*.json, letto dal motore Rust con IDRO_ASSET)
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "3")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "3")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_core as ac
import aggfull_l1 as L
import agg_zz_core as azc
from agg_zz_core import Level

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
PRE = json.load(open(os.path.join(OUT, "idrogeno_bohr_fase4_prereg.json"))) if os.path.exists(
    os.path.join(OUT, "idrogeno_bohr_fase4_prereg.json")) else None
EPS_M = 1e-5          # budget sulla variazione della massa fra due maestri successivi (pre-registrato)
JMAX = 4


def corsa(j):
    th = ac.TH_A / 2 ** j
    D = L.dirs13()
    fine = Level("fine3", 3, "dirac4", 4, 1, th, period=2)
    prev = fine
    levels = []
    t00 = time.time()
    for Bc, kfac, wu in ((256, 1 / ac.A0, 1e-3), (512, 1 / (2 * ac.A0), 3e-4)):
        B = Bc * 2 ** j                      # fine cells of the master per coarse cell (same physical cell)
        b = B // prev.B
        t0 = time.time()
        cons = L.taylor6(prev, b, D, kmax=min(0.224 * B * th, 1.5))
        rng = np.random.default_rng(20261015)
        kap, Gs, sw, phi0 = L.samples(prev, b, kfac * th * B, rng, 8000, 0.25, wu)
        chosen = None
        for w in (4, 5):
            r = L.fit(kap, Gs, sw, phi0, w, "T", cons=cons, D=D)
            cand = Level("learned", 3, "pauli2", 2, B, th, rule=r, b=b, w=w)
            # observation time kept PHYSICAL (1000 ticks of the theta_A master = 1000 * 2^j ticks here)
            ck = azc.kcheck(prev, cand, b, np.random.default_rng(Bc + w), 1000 * 2 ** j, kfac * th,
                            sig=max(1.0, 2.0 / (kfac * th * B)), k0_cells=0.5 * kfac * th * B, energy_mode="state")
            ok = ck["packet"] <= 1e-3 and ck["energy_rel"] <= 3e-4 and ck["rest_rel"] <= 3e-4 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01
            print(f"j={j} B={B} w={w}: " + ", ".join(f"{k} {v:.2e}" for k, v in ck.items() if isinstance(v, float)) + f" ok={ok} ({time.time() - t0:.0f}s)", flush=True)
            if ok:
                chosen = (w, r, cand, ck)
                break
        assert chosen is not None, f"j={j} B={B}: nessun candidato passa"
        w, r, cand, ck = chosen
        g, offs = L.band_of(r)
        cf_fine = L.taylor(fine, B, D)
        m_rel, m_abs, fr = L.mass_k4(g, offs, fine, B, B * th, cf_fine)
        kk = np.array([[0.2, 0, 0], [0.4, 0, 0], [0.785, 0, 0], [0.4, 0.4, 0.4], [0.8, 0.8, 0.8]])
        gf = L.fine_band_phase(fine, kk, B)
        gl = L.gband(g, offs, kk)
        h = B * th
        c4d = float(np.mean(cons[:, 1] / (-1 / (8 * h ** 3))))
        c6d = float(np.mean(cons[:, 2] / (1 / (16 * h ** 5))))
        c2m = float(np.mean(1 / (2 * cons[:, 0] * h)))     # mass from the pinned c2 of the level below
        print(f"j={j} B={B}: w {w}, phi0 {r.phi0:.9f}, m_abs(regola) {m_abs:.10f}, m(vincolo) {c2m:.10f}, c4/Dirac {c4d:.6f}, "
              f"c6/Dirac {c6d:.5f}, banda/maestro-1 " + ", ".join(f"{np.linalg.norm(k):.3f}:{gl[i] / gf[i] - 1:+.2e}" for i, k in enumerate(kk))
              + f" ({time.time() - t0:.0f} s)", flush=True)
        n = len(offs)
        levels.append({
            "B_fine_cells": B, "b": b, "w": w, "descriptor": "pauli2", "components": 2, "cell_lambdaC": h,
            "cell_over_a0": h / ac.A0, "rest_phase_phi0": float(r.phi0), "offsets": np.asarray(offs).astype(int).tolist(),
            "generator_G": {"re": [[list(map(float, g)), [0.0] * n], [[0.0] * n, list(map(float, g))]],
                            "im": [[[0.0] * n, [0.0] * n], [[0.0] * n, [0.0] * n]]},
            "spin_part_rel": 0.0, "learner": "T con c6 (AGGLOMERATO_COMPLETO §12) dal maestro theta_A/2^%d" % j,
            "controlli": {"kcheck": {k: float(v) for k, v in ck.items() if isinstance(v, (int, float))},
                          "m_assoluta": float(m_abs), "m_vincolo": c2m, "c4_su_dirac": c4d, "c6_su_dirac": c6d,
                          "k4_frazione": fr, "banda_su_maestro_meno_1": {f"{np.linalg.norm(k):.3f}": float(gl[i] / gf[i] - 1) for i, k in enumerate(kk)}}})
        prev = cand
    res = {"j": j, "theta_maestro": th, "theta_maestro_su_a0": th / ac.A0, "massa_maestro_tan_su_theta": float(np.tan(th) / th),
           "levels": levels, "secondi": time.time() - t00}
    json.dump(res, open(os.path.join(OUT, f"idrogeno_bohr_fase4_regola_j{j}.json"), "w"), default=float)
    print(f"scritto lab/out/idrogeno_bohr_fase4_regola_j{j}.json ({time.time() - t00:.0f} s)")


def scegli():
    rs = []
    for j in range(JMAX + 1):
        fn = os.path.join(OUT, f"idrogeno_bohr_fase4_regola_j{j}.json")
        if os.path.exists(fn):
            rs.append(json.load(open(fn)))
    m = [r["levels"][0]["controlli"]["m_assoluta"] for r in rs]
    scelto = None
    for i in range(1, len(rs)):
        d = abs(m[i] - m[i - 1])
        print(f"j={rs[i]['j']}: m = {m[i]:.10f}, |m_j - m_(j-1)| = {d:.3e} (budget {EPS_M:.0e})")
        if d < EPS_M and scelto is None:
            scelto = rs[i]
    assert scelto is not None, "nessun maestro converge entro JMAX: FAIL dichiarato"
    rich = (4 * m[-1] - m[-2]) / 3 if len(m) >= 2 else None
    out = {"what": "catena A rifatta col learner della fase 3 (P4) dal maestro scelto dall'agglomeratore (convergenza nella finezza "
                   "del maestro, budget |dm| < %.0e) - studio idrogeno_bohr fase 4" % EPS_M,
           "script": "lab/idrogeno_bohr_fase4_regola.py", "maestro_j": scelto["j"],
           "masse_per_j": {str(r["j"]): mm for r, mm in zip(rs, m)},
           "controllo_richardson_ultimi_due": rich,
           "chain_A_3d": {"theta_fine": scelto["theta_maestro"], "fine_cell_lambdaC": scelto["theta_maestro"], "levels": scelto["levels"]}}
    json.dump(out, open(os.path.join(OUT, "idrogeno_bohr_fase4_regola.json"), "w"), default=float)
    print(f"scelto j = {scelto['j']} (theta = a0/{ac.A0 / scelto['theta_maestro']:.0f}); Richardson (ultimi due) m = {rich}")
    print("scritto lab/out/idrogeno_bohr_fase4_regola.json")


if __name__ == "__main__":
    a = sys.argv[1]
    if a == "scegli":
        scegli()
    else:
        corsa(int(a))

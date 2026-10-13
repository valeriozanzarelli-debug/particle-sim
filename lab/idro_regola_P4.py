"""Regola grossa 'P4' per l'idrogeno di Bohr (studio idrogeno_bohr fase 3, esame P3).
Rifà ESATTAMENTE il metodo della fase 4 dell'agglomeratore (docs/AGGLOMERATO_COMPLETO.md §12, lab/aggfull_relearnT.py,
modo "T"): Taylor c2, c4, c6 in 13 direzioni FISSATI al livello sotto (taylor6, kmax = min(0.224 B th, 1.5)), parità
esatta, seme 20261015, w = il più piccolo fra (4, 5) che passa il controllo agg_zz_core.kcheck. B = 256 dal cammino fine
(b = 256), B = 512 dal B = 256 imparato (b = 2). Il pickle della fase 4 non è nel repo (non si committano .pkl): si rifà
qui, senza toccare lab/out/aggfull_relearnT.* (output dell'altro thread). Scrive assets/atomo/agglomerato_zz_P4_v1.json
(nuovo, versionato; stesso formato di agglomerato_zz_T_v1.json). Uso: python3 lab/idro_regola_P4.py"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_core as ac
import aggfull_l1 as L
import agg_zz_core as azc
from agg_zz_core import Level

TH = ac.TH_A
D = L.dirs13()
fine = Level("fine3", 3, "dirac4", 4, 1, TH, period=2)
prev = fine
out = {"what": "catena A rifatta con il learner della fase 4 di AGGLOMERATO_COMPLETO §12 (Taylor c2,c4,c6 fissati + parita', "
               "B=512 dal B=256) — studio idrogeno_bohr fase 3", "script": "lab/idro_regola_P4.py",
       "derived_from": "lab/aggfull_relearnT.py (fit mode T), cammino fine zig-zag a theta = a0/4096",
       "chain_A_3d": {"theta_fine": TH, "fine_cell_lambdaC": TH, "levels": []}}
cf_fine = {B: L.taylor(fine, B, D) for B in (256, 512)}
for B, kfac, wu in ((256, 1 / ac.A0, 1e-3), (512, 1 / (2 * ac.A0), 3e-4)):
    b = B // prev.B
    t0 = time.time()
    cons = L.taylor6(prev, b, D, kmax=min(0.224 * B * TH, 1.5))
    rng = np.random.default_rng(20261015)
    kap, Gs, sw, phi0 = L.samples(prev, b, kfac * TH * B, rng, 8000, 0.25, wu)
    chosen = None
    for w in (4, 5):
        r = L.fit(kap, Gs, sw, phi0, w, "T", cons=cons, D=D)
        cand = Level("learned", 3, "pauli2", 2, B, TH, rule=r, b=b, w=w)
        ck = azc.kcheck(prev, cand, b, np.random.default_rng(B + w), 1000, kfac * TH, sig=max(1.0, 2.0 / (kfac * TH * B)),
                        k0_cells=0.5 * kfac * TH * B, energy_mode="state")
        ok = ck["packet"] <= 1e-3 and ck["energy_rel"] <= 3e-4 and ck["rest_rel"] <= 3e-4 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01
        print(f"B={B} w={w}: " + ", ".join(f"{k} {v:.2e}" for k, v in ck.items() if isinstance(v, float)) + f" ok={ok} ({time.time() - t0:.0f}s)", flush=True)
        if ok:
            chosen = (w, r, cand, ck)
            break
    assert chosen is not None, f"B={B}: nessun candidato passa (come fase 4 non dovrebbe succedere)"
    w, r, cand, ck = chosen
    g, offs = L.band_of(r)
    # controls against the FINE WALK at this cell (as for the T_v1 asset): mass, k4, band at larger kappa
    m_rel, m_abs, fr = L.mass_k4(g, offs, fine, B, B * TH, cf_fine[B])
    kk = np.array([[0.2, 0, 0], [0.4, 0, 0], [0.785, 0, 0], [0.4, 0.4, 0.4], [0.8, 0.8, 0.8]])
    gf = L.fine_band_phase(fine, kk, B)
    gl = L.gband(g, offs, kk)
    c6d = float(np.mean(cons[:, 2] / (1 / (16 * (B * TH) ** 5))))
    print(f"B={B}: w {w}, phi0 {r.phi0:.6f}, m/m_fine-1 {m_rel:+.2e}, m_abs {m_abs}, c6(vincolo)/Dirac {c6d:.5f}, band/fine-1 "
          + ", ".join(f"{np.linalg.norm(k):.3f}:{gl[i] / gf[i] - 1:+.2e}" for i, k in enumerate(kk)) + f" ({time.time() - t0:.0f} s)", flush=True)
    n = len(offs)
    out["chain_A_3d"]["levels"].append({
        "B_fine_cells": B, "b": b, "w": w, "descriptor": "pauli2", "components": 2, "cell_lambdaC": B * TH,
        "cell_over_a0": B * TH / ac.A0, "rest_phase_phi0": float(r.phi0), "offsets": np.asarray(offs).astype(int).tolist(),
        "generator_G": {"re": [[list(map(float, g)), [0.0] * n], [[0.0] * n, list(map(float, g))]],
                        "im": [[[0.0] * n, [0.0] * n], [[0.0] * n, [0.0] * n]]},
        "spin_part_rel": 0.0, "learner": "T con c6 (AGGLOMERATO_COMPLETO §12, fase 4)",
        "controlli": {"kcheck": {k: float(v) for k, v in ck.items() if isinstance(v, (int, float))},
                      "m_su_m_fine_meno_1": float(m_rel), "m_assoluta": m_abs if not isinstance(m_abs, np.ndarray) else m_abs.tolist(),
                      "c6_vincolo_su_dirac": c6d, "k4_frazione": fr if isinstance(fr, dict) else str(fr),
                      "banda_su_fine_meno_1": {f"{np.linalg.norm(k):.3f}": float(gl[i] / gf[i] - 1) for i, k in enumerate(kk)}}})
    prev = cand
json.dump(out, open("assets/atomo/agglomerato_zz_P4_v1.json", "w"), default=float)
print("scritto assets/atomo/agglomerato_zz_P4_v1.json")

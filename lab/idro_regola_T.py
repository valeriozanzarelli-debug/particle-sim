"""Regola grossa 'T' per l'idrogeno di Bohr (studio idrogeno_bohr fase 2, emendamento 2).
Adotta il metodo di docs/AGGLOMERATO_COMPLETO.md §11 (lab/aggfull_l1.py, learner "T": parte senza spin pari
e reale, coefficienti di Taylor di massa e k^4 FISSATI esattamente a quelli del cammino fine in 13 direzioni, G(0)=0)
ai livelli che usa l'idrogeno: catena A, B = 256 (a0/16) e B = 512 (a0/8), imparati direttamente dal cammino fine
(b = B tick), stencil w = 4 (729 posti), osservabile 1s di Z = 1, seme 1. NON tocca agglomerato_zz.json:
scrive assets/atomo/agglomerato_zz_T_v1.json (stesso formato, spin scartato).
Uso: python3 lab/idro_regola_T.py"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_l1 as L
import agg_zz_core as azc

th = L.THETA["A"]
D = L.dirs13()
out = {"what": "catena A rifatta con il learner T di AGGLOMERATO_COMPLETO §11 (Taylor fissato + parita') — studio idrogeno_bohr fase 2",
       "derived_from": "lab/aggfull_l1.py (fit mode T), cammino fine zig-zag a theta = a0/4096", "script": "lab/idro_regola_T.py",
       "chain_A_3d": {"theta_fine": th, "fine_cell_lambdaC": th, "levels": []}}
for B in (256, 512):
    t0 = time.time()
    prev = azc.Level("fine3", 3, "dirac4", 4, 1, th, period=2)
    h_lc = B * th
    cf = L.taylor(prev, B, D)
    kap, Gs, sw, phi0 = L.samples(prev, B, 1.0 / L.A0 * th * B, np.random.default_rng(1), 8000, 0.25, 1e-3)
    r = L.fit(kap, Gs, sw, phi0, 4, "T", cons=cf, D=D)
    g, offs = L.band_of(r)
    m_rel, m_abs, fr = L.mass_k4(g, offs, prev, B, h_lc, cf)
    # parity of the raw rule and a few band checks against the fine walk at larger kappa
    kk = np.array([[0.2, 0, 0], [0.4, 0, 0], [0.8, 0, 0], [0.4, 0.4, 0.4], [0.8, 0.8, 0.8]])
    gf = L.fine_band_phase(prev, kk, B)
    gl = L.gband(g, offs, kk)
    print(f"B={B}: phi0 {phi0:.6f}, m/m_fine-1 {m_rel:+.2e}, m_abs {m_abs}, k4 frac {fr}, band/fine-1 at kappa "
          + ", ".join(f"{np.linalg.norm(k):.2f}:{gl[i] / gf[i] - 1:+.2e}" for i, k in enumerate(kk)) + f" ({time.time() - t0:.0f} s)", flush=True)
    n = len(offs)
    out["chain_A_3d"]["levels"].append({
        "B_fine_cells": B, "b": B, "w": 4, "descriptor": "pauli2", "components": 2, "cell_lambdaC": h_lc,
        "cell_over_a0": h_lc / L.A0, "rest_phase_phi0": float(phi0), "offsets": offs.astype(int).tolist(),
        "generator_G": {"re": [[list(map(float, g)), [0.0] * n], [[0.0] * n, list(map(float, g))]],
                        "im": [[[0.0] * n, [0.0] * n], [[0.0] * n, [0.0] * n]]},
        "spin_part_rel": 0.0, "learner": "T (AGGLOMERATO_COMPLETO §11)",
        "controlli": {"m_su_m_fine_meno_1": float(m_rel), "k4_frazione": fr if isinstance(fr, dict) else str(fr),
                      "banda_su_fine_meno_1": {f"{np.linalg.norm(k):.3f}": float(gl[i] / gf[i] - 1) for i, k in enumerate(kk)}}})
json.dump(out, open("assets/atomo/agglomerato_zz_T_v1.json", "w"), default=float)
print("scritto assets/atomo/agglomerato_zz_T_v1.json")

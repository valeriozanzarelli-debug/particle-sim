"""Regola imparata simmetrizzata (studio idrogeno_bohr fase 2, esame R3).
Legge assets/atomo/agglomerato_zz.json (NON lo modifica) e scrive assets/atomo/agglomerato_zz_simm_v1.json:
per ogni livello della catena A il generatore scalare g = (G00+G11)/2 viene mediato sulle 48 operazioni del
gruppo del cubo O_h (permutazioni degli assi e cambi di segno, parità inclusa). Lo spin imparato (~1e-10) è
scartato (G01 = 0, G00 = G11 = g_simm). Nessun parametro: è la proiezione sulla simmetria del reticolo.
Uso: python3 lab/idro_simm_regola.py"""
import json, itertools
import numpy as np

src = "assets/atomo/agglomerato_zz.json"
dst = "assets/atomo/agglomerato_zz_simm_v1.json"
d = json.load(open(src))
ch = d["chain_A_3d"]
out = {"what": "versione simmetrizzata (parità + cubo O_h) della catena A di agglomerato_zz.json — studio idrogeno_bohr fase 2",
       "derived_from": src, "script": "lab/idro_simm_regola.py",
       "how": "g_simm(o) = media su 48 operazioni di O_h di g(R o), con g = (G00+G11)/2; spin scartato; phi0, cella, offsets invariati",
       "units": d.get("units"), "chain_A_3d": {"theta_fine": ch["theta_fine"], "fine_cell_lambdaC": ch["fine_cell_lambdaC"], "levels": []}}
ops = [(p, s) for p in itertools.permutations(range(3)) for s in itertools.product([1, -1], repeat=3)]
for lv in ch["levels"]:
    offs = [tuple(o) for o in lv["offsets"]]
    idx = {o: i for i, o in enumerate(offs)}
    G = lv["generator_G"]
    g = 0.5 * (np.array(G["re"][0][0]) + 1j * np.array(G["im"][0][0]) + np.array(G["re"][1][1]) + 1j * np.array(G["im"][1][1]))
    gs = np.zeros_like(g)
    for p, s in ops:
        for i, o in enumerate(offs):
            ro = tuple(s[a] * o[p[a]] for a in range(3))
            gs[i] += g[idx[ro]]
    gs /= len(ops)
    # statistics of what changed
    odd = 0.5 * (g - np.array([g[idx[tuple(-x for x in o)]] for o in offs]))
    print(f"B={lv['B_fine_cells']}: max|Im g| {np.abs(g.imag).max():.3e} -> {np.abs(gs.imag).max():.3e}; "
          f"max|g - g_simm| {np.abs(g - gs).max():.3e} (max|g| {np.abs(g).max():.3e}); parity-odd part {np.abs(odd).max():.3e}; "
          f"sum g {g.sum().real:.6e} -> {gs.sum().real:.6e}")
    z = [[0.0] * len(offs)]
    nl = {k: v for k, v in lv.items() if k not in ("generator_G", "check")}
    nl["generator_G"] = {"re": [[list(gs.real), [0.0] * len(offs)], [[0.0] * len(offs), list(gs.real)]],
                         "im": [[[0.0] * len(offs), [0.0] * len(offs)], [[0.0] * len(offs), [0.0] * len(offs)]]}
    nl["spin_part_rel"] = 0.0
    nl["simmetrizzazione"] = {"max_im_prima": float(np.abs(g.imag).max()), "max_cambio": float(np.abs(g - gs).max()),
                              "parte_dispari_max": float(np.abs(odd).max())}
    out["chain_A_3d"]["levels"].append(nl)
json.dump(out, open(dst, "w"))
print("scritto", dst)

"""tavg_verdetti_o.py - verdetti pre-registrati (tavg_prereg.json fase 1, tavg_prereg2.json fase 2) per OSSIGENO
(con N, C, H dove servono). Legge solo i risultati gia' scritti; scrive lab/out/tavg_verdetti_ossigeno.json."""
import os, json
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
EV = 27.211386245988
L = lambda f: json.load(open(os.path.join(OUT, f)))
P1, P2 = L("tavg_prereg.json"), L("tavg_prereg2.json")
fine, grid, blk, simm, corr = L("tavg_fine.json"), L("tavg_grid.json"), L("tavg_blocks.json"), L("tavg_simm.json"), \
    L("tavg_corr.json")
KEY = P1["esami_fase1"]
v = {}


def put(name, val, ok, nota=""):
    v[name] = dict(valore=val, esito="PASS" if ok else "FAIL", nota=nota)


for Z in (7, 8):
    q0, q1 = fine[str(Z)]["q0"], fine[str(Z)]["q1"]
    k0 = KEY["F1"]["chiavi"][str(Z)][0]; k1 = KEY["F2"]["chiavi"][str(Z)]
    put(f"F1_Z{Z}", q0["E"] - k0, abs(q0["E"] - k0) <= 1e-5)
    put(f"F2_Z{Z}", q1["E"] - k1, abs(q1["E"] - k1) <= 1e-5)
    best = q0["candidati"][q0["migliore"]]
    put(f"F3_Z{Z}", dict(massimi=best["maxima"], kato=best["kato"]),
        best["maxima"] == 2 and abs(best["kato"] / Z - 1) < 5e-3)
    put(f"F4_Z{Z}", [q0["migliore"], q1["migliore"]],
        q0["migliore"].endswith(KEY["F4"]["chiavi"]["neutro"][str(Z)])
        and q1["migliore"].endswith(KEY["F4"]["chiavi"]["catione"][str(Z)]))
    b = blk[str(Z)]
    put(f"B3_Z{Z}", b["B3"]["err_Eh"], abs(b["B3"]["err_Eh"]) <= 1.6e-3, "FAIL atteso (pre-registrato)")
    gb = json.load(open(os.path.join(OUT, "tavg_blocchi", f"Z{Z:02d}.json")))
    dE0 = grid[f"atom_Z{Z}_q0"]["E_val"] - b["E_ps_ref"]
    ps1 = b["E_ps_ref"] + (b["E_fine_catione"] - b["E_fine"]) + b["prove"][-1]["err_catione_Eh"]
    dE1 = grid[f"atom_Z{Z}_q1"]["E_val"] - ps1
    sfer0 = Z == 7
    ok0 = abs(dE0) <= 1.6e-3 if sfer0 else (-10e-3 <= dE0 <= 1.6e-3)
    put(f"G1_Z{Z}_neutro", dE0, ok0, "p^4 non sferico: finestra -10..+1.6 mEh" if not sfer0 else "")
    put(f"G1_Z{Z}_catione", dE1, abs(dE1) <= 1.6e-3)
    g3n = grid[f"atom_Z{Z}_q0"]["E"] - k0; g3c = grid[f"atom_Z{Z}_q1"]["E"] - k1
    put(f"G3_Z{Z}_neutro", g3n, abs(g3n) <= 1.6e-3 if sfer0 else (-10e-3 <= g3n <= 1.6e-3))
    put(f"G3_Z{Z}_catione", g3c, abs(g3c) <= 2.5e-3)
    ipg = (grid[f"atom_Z{Z}_q1"]["E"] - grid[f"atom_Z{Z}_q0"]["E"]) * EV
    ipf = (q1["E"] - q0["E"]) * EV
    put(f"G4_Z{Z}", dict(IP_griglia=ipg, IP_fine=ipf, diff_eV=ipg - ipf), abs(ipg - ipf) <= 0.05,
        "O: la differenza e' la rottura di simmetria HF del neutro (-3.09 mEh = -0.084 eV); con orbitali p uguali "
        "(tavg_simm 'fermo') sarebbe +0.003 eV" if Z == 8 else "")
    put(f"REL_Z{Z}", "vedi tavg_fine.json (dIP_rel)", True)
c = grid["conv_Z8_h1"]["E"] - grid["conv_Z8_h2"]["E"]
put("G2_Z8", c, abs(c) <= 1.6e-3)
d = grid["learned_Z8"]["E"] - grid["atom_Z8_q0"]["E"]
put("G5_Z8", d, abs(d) <= 1.6e-3, "FAIL atteso (pre-registrato)")
for Z in (6, 8):
    d = grid[f"shift_Z{Z}"]["E"] - grid[f"atom_Z{Z}_q0"]["E"]
    put(f"G6_Z{Z}", d, abs(d) <= 1e-4)
v["G6_Z8"]["nota"] = ("diagnostica non contata: a cella a0/6 l'effetto uova e' 6e-7 Eh (tavg_simm uova_h6_*): "
                      "la regola di scelta di h dovrebbe guardare anche le uova")
put("G7_H2", grid["mol_H2"]["E"] + 1.133629, abs(grid["mol_H2"]["E"] + 1.133629) <= 3e-3)
put("G7_N2", grid["mol_N2"]["E"] + 108.993826, abs(grid["mol_N2"]["E"] + 108.993826) <= 3e-3,
    "corso dopo la correzione di un bug di hf_griglia (nome Y riusato: rompeva le molecole con proiettori KB)")
put("H_atomo_griglia", grid["atom_Z1_q0"]["E"] + 0.5, abs(grid["atom_Z1_q0"]["E"] + 0.5) <= 1.6e-3, "riportato")
v["simmetria_O"] = dict(esito="DICHIARATO", valore={k: simm[k]["dE_vs_radiale_mEh"] for k in simm if "Z8_h" in k},
                        nota="fermo = tutti i p con la stessa funzione radiale del livello fine: +0.11 mEh (a0/4), "
                             "+0.93 (a0/6, L=12); rilassato -3.09 / -3.05 mEh a entrambe le celle; <r> del p doppio "
                             "1.276 contro 1.204: rottura di simmetria vera dell'HF ristretto (orbitali non "
                             "equivalenti), non la griglia e non UHF (UHF: -10.4 mEh)")
# fase 2
cc = corr["cieco"]
ip = {}
for X, Z in (("O", 8), ("N", 7)):
    ipg = (grid[f"atom_Z{Z}_q1"]["E"] - grid[f"atom_Z{Z}_q0"]["E"]) * EV
    ip[X] = ipg + (cc[X + "+"]["Ec_pred"] - cc[X]["Ec_pred"]) * EV
    ip[X + "_fine"] = (cc[X + "+"]["E_pred"] - cc[X]["E_pred"]) * EV
E2 = P2["esami"]
ecO = cc["O"]["Ec_pred"]; k = E2["C1"]["chiave_Eh"]
put("C1", dict(Ec_O=ecO, chiave=k, rel=ecO / k - 1), abs(ecO / k - 1) <= 0.05)
put("C2", dict(IP_O=ip["O"], IP_O_con_HF_fine=ip["O_fine"], rel=ip["O"] / 13.618 - 1),
    abs(ip["O"] / 13.618 - 1) <= 0.03)
put("C3", dict(IP_N=ip["N"], IP_O=ip["O"]), ip["N"] > ip["O"])
put("C4", ip["N"] - ip["O"], abs(ip["N"] - ip["O"] - 0.916) <= 0.30)
put("C5", dict(IP_N=ip["N"], rel=ip["N"] / 14.534 - 1), abs(ip["N"] / 14.534 - 1) <= 0.03)
v["R1"] = dict(esito="RIPORTATO", valore=dict(Ec_N=cc["N"]["Ec_pred"], chiave=-54.5892 + 54.400934))
v["R2"] = dict(esito="RIPORTATO", valore=dict(dEc_O_mEh=1e3 * (cc["O+"]["Ec_pred"] - cc["O"]["Ec_pred"]), chiave=63.7))
v["modello"] = dict(scelto=corr["scelto"], cv={m: corr["models"][m]["cv_rms_mEh"] for m in corr["models"]},
                    fit=corr["models"][corr["scelto"]]["fit"],
                    nota="d al bordo della griglia pre-registrata (0.05): i dati non fissano d")
v["regge_su_ossigeno"] = v["C1"]["esito"] == "PASS" and v["C2"]["esito"] == "PASS"
json.dump(v, open(os.path.join(OUT, "tavg_verdetti_ossigeno.json"), "w"), indent=1, ensure_ascii=False)
for kk, x in v.items():
    print(kk, x if not isinstance(x, dict) else (x.get("esito"), x.get("valore")))

"""tavg_f5_tavola.py - lit oxygen in the periodic table: lab/out/tavg_dati_ossigeno/ -> assets/esportazioni/O/ (docs/ATOMO.md par. 10.2).

    python lab/tavg_f5_export.py && python lab/tavg_f5_tavola.py

Tier 'approssimato'. The samples are electrons drawn INDEPENDENTLY from the spherical HF radial density: no correlation and no p
shape in the picture (declared in the manifest, shown on screen by the table); the correlation is only in the energies.
"""
import json
import os

import numpy as np

LAB = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(LAB, "out", "tavg_dati_ossigeno")
DST = os.path.join(LAB, "..", "assets", "esportazioni", "O")
N = 4000

S = np.load(os.path.join(SRC, "ossigeno_campioni.npy"))[:N]
liv = json.load(open(os.path.join(SRC, "ossigeno_livelli.json"), encoding="utf8"))
print(liv)
os.makedirs(DST, exist_ok=True)
S.astype("<f4").tofile(os.path.join(DST, "campioni.f32"))
esami = [{"id": k, "cosa": "esame di docs/TAVOLA_GRIGLIA.md, fase 5 (lab/out/tavg_f5_verdetti_em1.json)", "PASS": v == "PASS"} for k, v in liv["esami"].items()]
esami.insert(0, {"id": "E_c", "cosa": f"correlazione {liv['E_c']:.4f} contro {liv['chiavi']['E_c']:.4f} Eh ({100 * (liv['E_c'] / liv['chiavi']['E_c'] - 1):+.2f} %)", "PASS": liv["esami"]["C1"] == "PASS"})
esami.insert(1, {"id": "IP", "cosa": f"ionizzazione {liv['IP_eV']:.3f} contro {liv['chiavi']['IP_eV']:.3f} eV ({100 * (liv['IP_eV'] / liv['chiavi']['IP_eV'] - 1):+.2f} %)", "PASS": liv["esami"]["C2"] == "PASS"})
man = {
    "formato": "particle-sim/esportazione-v1", "nome": "ossigeno", "simbolo": "O", "tipo": "atomo", "Z": 8, "carica": 0, "livello": "approssimato",
    "metodo": "O 3P: ROHF del livello fine a tutti gli elettroni + CCSD(T) a tutti gli elettroni, X=4,5, estrapolazione X^-3, zero manopole; blocco-nucleo non usato",
    "fonte": "lab/tavg_f5_export.py, lab/out/tavg_dati_ossigeno/, docs/TAVOLA_GRIGLIA.md fasi 4-5 (seconda prova dopo il C1 della fase 4: non primo cieco)",
    "nuclei": [{"Z": 8, "pos": [0.0, 0.0, 0.0]}], "n_elettroni": 8, "simmetria": "sferica",
    "energia_Eh": liv.get("E_tot", liv.get("E_totale")), "esami": esami,
    "note": ["i bagliori sono elettroni pescati INDIPENDENTI dalla densita' HF radiale: niente correlazione e niente forma p nel disegno; la correlazione e' solo nelle energie",
             "la buca di Coulomb qui non puo' uscire: i campioni sono indipendenti", "debiti: triple perturbative, base finita + X^-3, nucleo fermo, niente relativita'"],
    "campioni": {"file": "campioni.f32", "n": int(len(S)), "unita": "bohr", "ordine": "n x n_elettroni x 3, float32 little endian"},
}
json.dump(man, open(os.path.join(DST, "manifest.json"), "w", encoding="utf8"), indent=1, ensure_ascii=False)
print("esami", [(e["id"], e["PASS"]) for e in esami], "energia", man["energia_Eh"])

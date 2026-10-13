"""fw1atomo_emendamenti.py - amendments to lab/out/fw1atomo_prereg.json, each written BEFORE the run it changes.
Writes lab/out/fw1atomo_emendamenti.json"""
import json, os
E = [
    {"id": "A1", "quando": "dopo P1a/P1b e dopo la prima Ritz di P2 per FW1, prima delle corse ritzc_ e fix_ di P2",
     "cosa": "P2, esami F1 e F3 per FW1: la Ritz da UNA partenza sull'intero cubo non trova livelli esatti di FW1 (nessun residuo < 1e-4): "
             "P1b ha mostrato che con FW1 gli 8 sottoreticoli di periodo hanno livelli quasi uguali ma non uguali (0,8694-0,8730 per la riga bassa), "
             "16 livelli in 4e-3 che 25 filtri larghi ~0,02 non separano. Siccome la regola di periodo U2 e' esattamente a blocchi sui sottoreticoli, "
             "la Ritz si fa sottoreticolo per sottoreticolo (partenza = 1s radiale ristretto al sottoreticolo c, c = 0..7, stessa finestra e durata); "
             "i pesi sono calcolati con la partenza intera. F3 parte dal livello principale del sottoreticolo 0. Soglie invariate.",
     "motivo": "limite dello strumento, non della soglia"},
    {"id": "A2", "quando": "stesso momento",
     "cosa": "P1d: dw = 0,01 invece di 0,005 (memoria a 80^3). Gia' nel codice prima delle corse P1d.",
     "motivo": "memoria < 2,5 GB"},
    {"id": "A3", "quando": "prima di ogni corsa di P3",
     "cosa": "P3: per lo stesso motivo di A1 la Ritz si fa per sottoreticolo di periodo: partenze 'dirac' ristrette ai sottoreticoli 0 e 1 "
             "(i due che in P1b danno gli estremi della riga bassa) e 'cieca' ristretta al sottoreticolo 0; F0 uguale. L'energia dell'1s di FW1 = "
             "livello con piu' peso nella partenza 'dirac' del sottoreticolo 0 (si riporta anche il sottoreticolo 1). Soglie invariate.",
     "motivo": "limite dello strumento"},
    {"id": "A4", "quando": "prima delle corse P1c, P1d, P1e",
     "cosa": "P1c, P1d, P1e: stessa ragione di A1, la Ritz si fa per sottoreticolo di periodo, partenza spin su ristretta ai sottoreticoli 0 e 1 "
             "(F0 in P1d solo 0: i suoi sottoreticoli sono uguali, P1b). Le regole di decisione di P1 restano quelle scritte; sono gia' state "
             "applicate a P1a/P1b, queste corse misurano solo le dipendenze da h, theta e posizione.",
     "motivo": "limite dello strumento"},
]
json.dump(E, open(os.path.join(os.path.dirname(__file__), "out", "fw1atomo_emendamenti.json"), "w"), indent=1, ensure_ascii=False)
print("ok")

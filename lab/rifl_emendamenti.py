"""rifl_emendamenti.py - cambi DOPO la pre-registrazione (dichiarati). Genera lab/out/rifl_emendamenti.json."""
import json, os
em = [
    {"id": "A1", "esame": "E3", "cambio": ("sequenza di a: theta = 4/(eps^2 N) con N intero che raddoppia e theta <= 0.4, "
                                          "invece di theta in {0.4,0.2,0.1,0.05}: per eps >= 2.5 quei theta danno N non interi"),
     "quando": "prima delle corse E3", "effetto_sulle_soglie": "nessuno"},
    {"id": "A2", "esame": "E0", "cambio": "U fino a 80 (non 160): a U = 160 la matrice non sta nei ~2,5 GB concessi",
     "quando": "prima della corsa E0", "effetto_sulle_soglie": "nessuno (soglia 2e-6 invariata)"},
    {"id": "A3", "esame": "E1", "cambio": ("una prima corsa E1 con griglia cosc = 12 era sotto-risolta (scarti 1e-3 al "
                                          "raddoppio); rifatta con cosc = 4 e controllo a cosc = 2. La prima corsa e' "
                                          "conservata in rifl_r1.json sotto 'E1_vecchio_cosc12_sottorisolto'"),
     "quando": "dopo la prima corsa E1 (difetto numerico, non di fisica)", "effetto_sulle_soglie": "nessuno"},
    {"id": "A4", "esame": "Q2", "cambio": "N fino a 2000 (non 1e4): forma di Gram 2N x 2N, memoria e tempo",
     "quando": "prima delle corse Q2", "effetto_sulle_soglie": "nessuna soglia (domanda aperta)"},
    {"id": "A5", "esame": "E1", "cambio": "taglio massimo in impulso scelto per tenere la matrice <= ~6000 (memoria); "
                                         "un processo con matrice ~10000 e' stato ucciso per memoria",
     "quando": "durante E1", "effetto_sulle_soglie": "nessuno"},
    {"id": "A6", "esame": "tutti", "cambio": ("difetto della griglia adattiva (pannello calcolato col ritmo di fase al bordo "
                                           "sinistro: primo pannello troppo lungo). Visto in E2 a U = 40 (0,0320 invece di "
                                           "0,0373). Corretto (ritmo al bordo destro) e RIFATTE tutte le corse r1, r2"),
     "quando": "dopo la prima corsa E2", "effetto_sulle_soglie": "nessuno"},
]
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rifl_emendamenti.json")
json.dump(em, open(out, "w"), indent=1, ensure_ascii=False)
print("scritto", out)

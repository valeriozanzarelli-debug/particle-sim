"""lv_prereg.py - pre-registrazione dello studio LIMITI_LORENTZ (prefisso lv_).
Scrive lab/out/lv_prereg.json PRIMA delle corse lv_r1.py (controllo numerico) e lv_r2.py (limiti)."""
import json
import os
import numpy as np
import lv_core as lv

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "lv_prereg.json")

dirs = {"asse (1,0,0)": [1, 0, 0], "faccia (1,1,0)": [1, 1, 0], "corpo (1,1,1)": [1, 1, 1],
        "generica (1,2,3)": [1, 2, 3], "generica (3,-1,2)": [3, -1, 2]}
pred_c3 = {}
for name, n in dirs.items():
    A = float(lv.A_of(n))
    pred_c3[name] = {"A": A, "luce: coeff k^3 di omega": -A / 6, "elettrone m=0: coeff k^3": -2 * A / 3}

# limiti previsti (formule analitiche al primo ordine; A = 1/3 caso peggiore, 1/5 media sulle direzioni)
E_QG2_LHAASO = 6e-8 * lv.E_PL
pred_bounds = {}
for A in (1 / 3, 1 / 5):
    a_tof = lv.HBARC_GEVM / (E_QG2_LHAASO * np.sqrt(A / 3))
    Eobs = 1.21e6
    Ea_min = np.sqrt(A) * Eobs ** 2 / (4 * lv.ME_GEV)
    a_dec = lv.HBARC_GEVM / Ea_min
    a_cav = lv.HBARC_GEVM / (1.165e-9 / np.sqrt(1e-17 * 6 / A))
    pred_bounds["A=%.4f" % A] = {
        "tempo di volo GRB 221009A (E_QG2>6e-8 E_Pl), a_max [m]": a_tof, "in l_P": a_tof / lv.L_P,
        "decadimento del fotone 1.21 PeV, a_max [m]": a_dec, "in l_P ": a_dec / lv.L_P,
        "cavita' ottica dc/c<1e-17 a 1.165 eV, a_max [m]": a_cav, "in l_P  ": a_cav / lv.L_P}

prereg = {
    "studio": "LIMITI_LORENTZ (lv_)",
    "data": "2026-10-13",
    "regola": "lab/zz_core.py congelata (bloch_two_ticks / bloch_batch), nessuna modifica",
    "parametri_liberi": 0,
    "ingressi": ["c, hbar, m_e (theta = m c a/hbar)", "passo a = incognita da limitare",
                 "limiti sperimentali (solo come esame)"],
    "derivazione_prima_dei_dati": {
        "luce": "cos(omega dt) = cos(kx a) cos(ky a) cos(kz a)  (composizione di rotazioni: quaternioni, sympy)",
        "elettrone": "cos(2 omega dt) = cos^2 th cos(2kx a)cos(2ky a)cos(2kz a) - sin^2 th "
                     "(trovata con una prova numerica rapida su 4 punti prima della pre-registrazione: dichiarato)",
        "identita'": "omega_e(k; th=0) = omega_gamma(2k)/2",
        "sviluppo luce": "omega = k [1 - (A/6)(ka)^2 + O((ka)^4)], A = nx^2ny^2+ny^2nz^2+nz^2nx^2",
        "sviluppo elettrone": "E^2 = m^2 + p^2 - (4/3) A p^4 a^2 - m^2 p^2 a^2/3",
        "ordine atteso": "nessun termine lineare in E (n=1) grazie al tick palindromico; primo termine n=2 "
                         "(dimensione 6), anisotropo con simmetria cubica, nullo sugli assi, massimo in diagonale; "
                         "NESSUNA birifrangenza (W(k) reale per la luce => |omega| uguale per le due elicita'); "
                         "elettrone 4 volte piu' anisotropo del fotone",
        "un tick solo (diagnosi)": "omega = k [1 -+ (nx ny nz)(ka)/2 + ...] : termine n=1 dispari in k",
        "cinematica": "gamma -> e+e- PERMESSO sopra k_th = 2 sqrt(m E_a) A^(-1/4); e -> e gamma VIETATO; "
                      "gamma -> 3 gamma VIETATO",
    },
    "previsioni_coefficienti": pred_c3,
    "esami": {
        "E1 forma chiusa luce": "max |omega_num - arccos(prod cos)| < 1e-12 su 2000 k casuali in [-pi/2, pi/2]^3",
        "E2 forma chiusa Dirac": "stesso con theta casuale in [0,1], < 1e-12; doppietto di spin degenere < 1e-12",
        "E3 coefficiente k^3": "dalla regola numerica (Richardson su eps = 1e-2, 5e-3, 2.5e-3): |c3_num - c3_analitico| "
                               "< 1e-6 per le 5 direzioni, luce ed elettrone (theta=0)",
        "E4 birifrangenza": "| |omega_+| - |omega_-| | < 1e-12 (luce, 2000 k); elettrone: spin degenere (E2)",
        "E5 niente termine n=1": "omega(k) = omega(-k) < 1e-12 e coefficiente di k^2 in omega/k-1 (fit) < 1e-8; "
                                 "diagnosi un tick: coefficiente = -+ nx ny nz / 2 entro 1e-6",
        "E6 simmetria cubica": "omega invariante sotto le 48 operazioni del cubo < 1e-12 (luce e Dirac theta=0.3)",
        "E7 identita'": "|omega_e(k) - omega_gamma(2k)/2| < 1e-12 per theta=0",
        "E8 soglie": "soglia numerica di gamma->ee con le dispersioni esatte (max su x) entro 1 % da k_th analitica "
                     "per 3 direzioni e 3 valori di E_a; e->e gamma e gamma->3 gamma: max Delta E < 0 su griglia",
        "E9 limiti": "a_max per ogni esperimento con A=1/3 e A=1/5; si dichiara il piu' severo (robusto) e quelli "
                     "condizionati a parte",
    },
    "previsione_limiti": pred_bounds,
    "previsione_piu_severo": "decadimento del fotone (fotoni PeV di LHAASO): a_max ~ 2-4e4 l_P ~ 5e-31 m; "
                             "tempo di volo GRB ~ 5e7 l_P; cavita' ottiche ~ 1e20 l_P",
    "criteri_esclusione": {
        "STRUTTURALE (esclusa per ogni a)": "un termine d=4 (anisotropia che non sparisce a k->0) o d=5 "
                                            "(lineare in E, o birifrangenza) nella dispersione libera",
        "esclusa per a >= l_P": "a_max (esperimento robusto, A = 1/5) < l_P",
        "in tensione": "l_P <= a_max < 10 l_P",
        "compatibile": "a_max >= 10 l_P",
        "fuori esame (dichiarato)": "correzioni radiative (Collins et al. 2004): servono le interazioni; si stima "
                                    "solo l'ordine di grandezza a parole; limiti adronici di Auger e protoni compositi: "
                                    "condizionati al modo in cui l'agglomeratore eredita la dispersione",
    },
    "laboratorio": "a = 0.05-1 lambda_C (1.9e-14 - 3.9e-13 m): previsto escluso da tutto, ~17-18 ordini sopra a_max",
}
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    json.dump(prereg, f, indent=1, ensure_ascii=False)
print(json.dumps(pred_bounds, indent=1))

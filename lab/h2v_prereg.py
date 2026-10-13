"""h2v_prereg.py - pre-registration of VIDEO H2 (writes lab/out/h2v_prereg.json). Committed BEFORE the wavefunctions are
rebuilt and before any frame is rendered."""
import json, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
P = {
    "study": "h2v (video della prima molecola: H2 con due elettroni)",
    "data_source": "Debiti, studio aggmol fase 2 (commit 5109ca9): lab/aggmol_ecg.py + lab/aggmol_h2.py, esiti lab/out/aggmol_h2.json. "
                   "Le basi ECG (.npz) non sono nel repo: si ricostruiscono con la stessa macchina, stessi semi, stessa catena di partenze calde.",
    "what_is_computed": [
        "funzione d'onda a due elettroni Psi(r1, r2) di H2 (ECG, 180 gaussiane, singoletto, gerade) a R = 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7 a0",
        "densita' di un elettrone rho(r) analitica (integrale gaussiano esatto sul secondo elettrone)",
        "densita' condizionata |Psi(r1_fisso, r2)|^2 (dove sta l'altro elettrone se uno e' qui)",
        "livelli vibrazionali chi_v(R) sulla curva di Debiti (18 punti, spline, massa ridotta m_p/2), J = 0",
    ],
    "what_is_chosen_declared": [
        "protoni fermi a R = 1.40 a0 nelle parti A e B (Born-Oppenheimer; R_e di Debiti = 1.4009, il punto di griglia 1.40 dista 0.0009 a0)",
        "il moto delle palline-elettrone e' un campionamento Metropolis di |Psi|^2 (passo piccolo per fotogramma): NON sono traiettorie; "
        "lo stato fondamentale e' reale, la sua corrente e' zero",
        "vibrazione (parte C): stato scelto = (chi_0 + chi_1)/sqrt(2), come una molecola che ha appena assorbito mezzo quanto; "
        "periodo calcolato, rallentato di un fattore dichiarato",
        "asse della molecola fisso nello spazio (la rotazione J = 0 sarebbe una sfera: debito)",
        "nella parte C la densita' degli elettroni e' la miscela sui punti di griglia R pesata con |chi(R,t)|^2 (discretizzazione in R dichiarata)",
        "estetica: colori, dimensione e alone delle palline, guadagno, camera, numero di palline",
    ],
    "exams": {
        "V0": {"what": "energie ECG ricostruite contro lab/out/aggmol_h2.json agli 8 R", "threshold": "|diff| <= 1e-6 Eh a ogni R"},
        "V1": {"what": "teorema del viriale a R = 1.40: <T> contro -E_tot - R dE/dR (dE/dR dalla spline della curva di Debiti)",
               "threshold": "|diff| <= 1e-3 Eh; e normalizzazione analitica int rho = 2 entro 1e-6"},
        "V2": {"what": "<1/r12> a R = 1.40 contro 0.5873 (valore di letteratura citato a memoria, Kolos-Wolniewicz/Rychlewski: dichiarato come da verificare)",
               "threshold": "1 %"},
        "V3": {"what": "palline contro densita' calcolata: densita' di colonna nel piano xz, celle 0.2 a0 su [-4,4]^2, >= 1e6 palline accumulate su fotogrammi distanti",
               "threshold": "correlazione di Pearson >= 0.995 e L1 relativo <= 0.05"},
        "V4": {"what": "correlazione elettrone-elettrone: P_opp = probabilita' che i due elettroni stiano da parti opposte del piano di mezzo",
               "threshold": "palline contro calcolo analitico (quadratura) entro 0.01; e P_opp >= 0.52 (senza correlazione, Hartree-Fock: 0.5 esatto)"},
        "V4b": {"what": "condizionata: elettrone di riferimento a z < -0.35 a0 (dietro il protone sinistro), P(altro a z > 0) dalle palline arancioni contro calcolo",
                "threshold": "entro 0.02"},
        "V5": {"what": "lunghezza di legame letta dall'immagine: distanza dei baricentri dei due bagliori dei protoni in una foto vista di lato", "threshold": "entro 2 % di 1.40 a0"},
        "V6": {"what": "periodo di <R>(t) misurato sui fotogrammi della parte C contro h/(E1 - E0) dai livelli di Debiti (8.006 fs)", "threshold": "1 %"},
        "V7": {"what": "la nuvola segue i protoni: <|z|> degli elettroni (palline) nei fotogrammi a R massimo e minimo contro il calcolo analitico della miscela",
               "threshold": "entro 0.03 a0 ciascuno, e il segno della differenza uguale"},
    },
    "knobs_physics": 0,
    "knobs_visual": ["palette", "SPOT", "GLOW", "GAIN", "numero di camminatori (12000)", "passo Metropolis per fotogramma",
                     "rallentamento del tempo", "camera", "stato vibrazionale scelto (50/50 v0, v1)"],
    "render": "stills 1920x1080; video 1280x720 24 fps (~36 s); CPU moderata (nice, 1 thread)",
}
json.dump(P, open(os.path.join(OUT, "h2v_prereg.json"), "w"), indent=1, ensure_ascii=False)
print("written")

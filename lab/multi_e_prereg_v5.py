"""multi_e_prereg_v5.py - EMENDAMENTO VERSIONATO v5 della pre-registrazione MOLTI ELETTRONI.
SOSTITUISCE la v4 (commit b90b7bb), che NON e' stata eseguita: nessun numero v4 e' stato visto (il driver aspettava la fine
delle DMC di D6' ed e' stato fermato prima di partire). Regola: FASE 5 di docs/AGGLOMERATO_COMPLETO.md (par. 13, aggfull_p5.py,
commit a5506f1): livello 1 con k^2, k^4, k^6 fissati + parita' (T6), larghezza scelta con la regola vera dell'agglomeratore
(la piu' piccola fra w = 4..7 che passa TUTTO il controllo, budget 3e-5 compreso, campione seme+100, 3000 punti), giudicata su
un campione INDIPENDENTE (seme+200, 10000 punti). Camminatori: Markov + k^6 a w = 8 (M6, fase 4, invariato in fase 5).
Il solutore pulito (aggfull_atom.Atom(..., clean=True)) NON si applica al nostro elio: l'ECG non proietta H su un settore
con LOBPCG (le gaussiane correlate sono gia' simmetriche per costruzione); lo diciamo e non cambia nulla.
Scritto e committato PRIMA delle corse v5. v1, v2, v3 restano valide; nessun verdetto vecchio si rigiudica.
Scrive lab/out/multi_e_prereg_v5.json.  Run: python lab/multi_e_prereg_v5.py"""
import json, os, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE).decode().strip()

P = {
 "meta": {
   "study": "MOLTI_ELETTRONI", "version": 5, "amends": "v1, v2, v3 (multi_e_prereg*.json)",
   "written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head,
   "decided_by": "sessione principale (coordinatore), default dichiarato: adottare la fase 5 (M2c 6/6 semi entro 3e-5, w = "
                 "5,5,6,7,5,5) come emendamento pre-registrato, ricontrollare i 2 semi che fallivano il Pauli (3 e 4) e "
                 "riportarli onestamente; il seme 5 (2,93e-5) e' entro il rumore del controllo (~1e-5).",
   "data_already_seen": "aggfull_p4.json e aggfull_p5.json letti: fase 5 w = 5/5/6/7/5/5, esame 0,55/0,98/0,73/0,42/2,93/1,23 e-5; "
                        "massa esatta a 1e-13; M6 fattibile solo a w = 8 (6e-11). Le nostre corse v5 sono RIPRODUZIONI con il "
                        "nostro codice; lo diciamo.",
   "pacing": "corse a pezzi <= 8 min con checkpoint (una sola corsa alla volta finche' girano le DMC di D6')",
 },
 "rule": {
   "learner": "aggfull_l1.samples / fit / band_of / taylor6 e agg_zz_core.kcheck, importati e non modificati; stessa ricetta di "
              "aggfull_p5.py per T6 e aggfull_p4.py per M6 (c246 = taylor6 della regola fine della catena, 13 direzioni). Livelli salvati con i nostri nomi "
              "lab/out/multi_e_v5_*.pkl (non committati), risultati lab/out/multi_e_v5.json (codice lab/multi_e_v5.py)",
   "band_level_T6": "catena A, osservabile Z = 10, learner T con c2 c4 c6 fissati, w con la regola di fase 5, seme 1 (la banda "
                    "dell'elio); semi 3 e 4 rifatti per il ricontrollo del Pauli",
   "walker_level_M6": "catena A, osservabile Z = 2, learner M con c2 c4 c6 fissati, w = 8, seme 1",
   "R3_levels": "T6 sulle catene A2, A, B, osservabile Z = 1, w = 5 fisso, seme 1 (la massa e' fissata: lo spread fra semi e' ~1e-13)",
 },
 "exams": {
   "R3v5": "m_kin(theta -> 0) dal nostro mkin_of sui livelli T6 di A2/A/B, fit m0 + a th^2 + b th^4: |m0 - 1| <= 5e-5 (soglia v1); "
           "riportato anche contro 1e-9",
   "L1c_v5": "controllo Z = 10 a 3e-5 (energia e resto) sul campione indipendente: i nostri semi 1, 3, 4 rifatti (3 e 4 sono "
             "quelli che fallivano) + i 6 semi di aggfull_p5.json letti. PASS se tutti stanno sotto 3e-5. Riportato se i nostri "
             "riproducono i loro (stesso w, energia entro 1e-7). Nota prima delle corse: il seme 5 (2,93e-5) passa con un "
             "margine piu' piccolo del rumore del controllo (~1e-5): un PASS li' non e' robusto, lo diciamo comunque",
   "G0v5": "parita' di T6 e M6: parte dispari / parte pari a k h = 0.01 e 1: < 1e-10",
   "R2v5": "livello M6: peso negativo di K_beta < 1e-10 per beta in {1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 1e-3, 1e-2} (griglie 64 "
           "e 128), |sum K - 1| < 1e-9, D curvatura = D spazio reale a 2e-5, |curtosi in eccesso a beta = 1e-2| < 0.02, "
           "0 salti col segno sbagliato. Riportato lo stesso per T6",
   "use_for_He_if": "R3v5, G0v5, R2v5 PASS (L1c_v4 non blocca: FAIL ad alto k gia' noto, dichiarato)",
 },
 "He_with_v5": {
   "M1v5": "stessa base ECG di v2 (he1S_Z2, K = 200), massa m0_v4: E = E_NR(m_A) x m0_v4 / m_A, tolleranza 0.29 mEh",
   "E1v5": "come E1 v2 con m0_v4 (tolleranza 10 meV)",
   "Q1v5": "Delta_band sull'elio al primo ordine sullo stato ECG NR con le bande T6 e M6 contro la cinetica di Dirac esatta "
           "sullo stesso stato: |band/Dirac - 1| <= 0.01 (T6) e <= 0.02 (M6) (stesse soglie di v3). DEBT-4/9: la parte non "
           "Dirac in mEh per livello vecchio, v3 T/M, v5 T6/M6",
   "D_walkers": "nessuna nuova DMC in v5 (i camminatori restano gaussiani; il pezzo di banda si somma al primo ordine)",
 },
 "kept_declared": "FAIL di L1c_v3 resta com'era; se uno dei nostri semi esce sopra 3e-5 si dichiara FAIL (niente nuovi semi)",
}
with open(os.path.join(OUT, "multi_e_prereg_v5.json"), "w") as f:
    json.dump(P, f, indent=1, ensure_ascii=False)
print("scritto lab/out/multi_e_prereg_v5.json", P["meta"]["written"], head)

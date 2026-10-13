"""senzadoppi_prereg.py - pre-registration of the study "Una regola zig-zag senza doppioni?" (docs/SENZA_DOPPIONI.md).
Written and committed BEFORE the production runs senzadoppi_r1.py, _r2.py, _r3.py. Changes after: separate amendment file."""
import json
import os

out = {
    "study": "senzadoppi_",
    "date": "2026-10-13",
    "question": "Esiste una regola 'corri a c lungo un asse, poi ruota' (locale, unitaria) senza doppioni dell'elettrone, "
                "o con doppioni innocui (massa ~1/a a quasi-energia 0 e pi)? Quale ipotesi del teorema W3=0 va rotta "
                "e quanto costa?",
    "frozen_rule_untouched": "lab/zz_core.py non viene modificato; ogni candidato e' una PROPOSTA.",
    "units": "c = hbar = 1, passo a = 1, theta = m c a / hbar; energie fisiche E*a; quasi-energie di un periodo divise "
             "per la durata T del periodo (tick).",
    "dev_seen_before_prereg": [
        "F0 (congelata): 8 copie, 4+4 nodi di Weyl per copia, E0 = theta, Epi = 0, isotropia 1,54 % a 24 celle/lambda, "
        "massa cinetica +3,37e-3 (theta = 0,1) -- gia' noti da REGOLA_LOCALE",
        "B1 (un tick, ordine fisso, = automa di Weyl bcc di D'Ariano-Perinotti): 2+2 nodi per copia, doppione SENZA massa "
        "a energia 0 nel punto P; isotropia 5,6 %",
        "B2 (un tick, L con ordine rovesciato): E0 = Epi = theta, isotropia 5,6 %; con Wilson ritorto sui collegamenti "
        "bcc (+-1,+-1,+-1) h = 0,1 / 0,3: E0 0,41 / 0,43, Epi 0,41 / 1,20",
        "PA (traslazioni della congelata, UNA moneta 2 theta per periodo): E0 = Epi = theta; con Wilson su +-2e_j "
        "h = 0,1 / 0,3 / 0,6: E0 0,22 / 0,44 / 0,24, Epi 0,14 / 0,31 / 0,38; isotropia 1,27 %; massa cinetica 1,36e-2",
        "Le soglie qui sotto NON sono state aggiustate su questi numeri; h = 0,3 e' scelto come valore unico di "
        "produzione perche' era il centro della griglia vista (dichiarato).",
    ],
    "candidates": {
        "F0": "regola congelata: tick pari ordine xyz per L e R, dispari zyx; moneta theta a ogni tick",
        "B1": "un tick, ordine xyz per L e R, moneta theta (automa di Weyl/Dirac bcc di D'Ariano-Perinotti)",
        "B2": "un tick, R in ordine xyz, L in ordine zyx (U_L = U_R^dagger), moneta theta",
        "B2W": "B2 + Wilson ritorto (asse Y di chiralita') sugli 8 collegamenti bcc (+-1,+-1,+-1), h = 0,3, meta' prima "
               "e meta' dopo la moneta",
        "PA": "traslazioni della congelata (xyz poi zyx) e UNA moneta di angolo 2 theta per periodo di due tick",
        "PAW": "PA + Wilson ritorto sui collegamenti +-2 e_j, h = 0,3, meta' prima e meta' dopo la moneta",
        "S2": "ricerca: tutti i 6^4 = 1296 periodi di due tick con ordini liberi per L e R e moneta a ogni tick",
    },
    "exams": {
        "E1_teorema": {
            "what": "W3 = (1/24 pi^2) int tr(U^-1 dU)^3 per 24 prodotti casuali di fattori a una direzione (2x2 e 4x4, "
                    "passi sugli assi, sulle diagonali bcc, rotazioni sul posto e su collegamenti); controllo positivo: "
                    "mappa non locale (sin k, M - sum cos k) normalizzata, M = 2 (grado atteso +-1) e M = 4 (0)",
            "pass": "|W3| < 0,01 per tutti i 24; controllo positivo |W3| entro 0,01 da +-1 e da 0",
            "prediction": "PASS: il teorema vale per qualsiasi numero di componenti e qualsiasi reticolo, perche' "
                          "tr(U^-1 dU)^3 di un fattore che dipende da una sola combinazione k.v e' zero e W3 e' "
                          "additivo sui prodotti",
        },
        "E2_censimento": {
            "what": "nodi di Weyl (parte R senza massa) per copia esatta, a quasi-energia 0 e pi; valli Dirac a theta = 0,1: "
                    "E0 = min |E| fuori da palle di raggio 0,45 attorno a Gamma e alle sue copie, Epi = distanza minima "
                    "dal bordo pi/T",
            "prediction": {
                "F0": "4 + 4 nodi; E0 = theta, Epi = 0 (validazione di REGOLA_LOCALE)",
                "B1": "2 + 2 nodi; E0 = 0 (doppione senza massa nel punto P: R a energia 0, L a pi)",
                "B2": "2 + 2 nodi; E0 = Epi = theta",
                "PA": "4 + 4 nodi; E0 = Epi = theta (la moneta una volta per periodo non si cancella alle valli pi)",
            },
            "pass": "ogni numero previsto entro 1e-6 (nodi: conteggio esatto)",
        },
        "E3_ricerca_ordini": {
            "what": "S2: per ognuno dei 1296 periodi: E0, Epi a theta = 0,1, isotropia a 24 celle/lambda",
            "prediction": "nessuno ha E0 > theta + 1e-6 (teorema: la parte R ha sempre almeno un gemello a energia 0, "
                          "che al massimo prende la massa theta); nessuno ha insieme Epi >= theta e isotropia <= 2 %",
            "pass": "previsione confermata (se un ordine la smentisce, e' una scoperta e va studiato)",
        },
        "E4_candidati": {
            "what": "PAW e B2W a theta = 0,1 e 0,02",
            "criteria": {
                "C1 doppioni a energia 0": "E0 >= 0,3",
                "C2 doppioni a pi": "Epi >= 0,3",
                "C3 massa a riposo": "|E(0)/m - 1| <= 1e-12",
                "C4 isotropia della luce del reticolo (theta = 0) a 24 celle/lambda": "<= 1,1 x congelata (<= 1,70 %)",
                "C5 massa cinetica": "|m_cin/m - 1| <= 2e-2 a theta = 0,1 e rapporto 0,1 -> 0,02 >= 15 (ordine theta^2)",
                "C6 velocita' massima di gruppo nella zona": "<= 1 (c)",
                "C7 unitarieta' (spazio reale, 50 periodi)": "<= 1e-12",
                "C8 Bloch == spazio reale (onda piana)": "<= 1e-12",
            },
            "prediction": {
                "PAW": "C1-C5, C7, C8 PASS; C6 FAIL (~1,03 c dai collegamenti di raggio 2)",
                "B2W": "C1-C3, C5, C7, C8 PASS; C4 FAIL (errore del primo ordine, ~5,6 %); C6 FAIL",
            },
            "verdict_rule": "PROPOSTA forte se C1-C5 PASS; altrimenti 'limite' con i costi scritti",
        },
        "E5_inversione_tempo": {
            "what": "coppie di Kramers con potenziale a caso (6^3, theta = 0,37): split massimo",
            "prediction": "F0 e PA <= 1e-12 (simmetria esatta); PAW e B2W rotta (split >= 1e-4) per il Wilson ritorto "
                          "(generatore -iY reale, REGOLA_LOCALE fase 5)",
            "pass": "previsione confermata",
        },
        "E6_coppie_senza_proiettore": {
            "what": "campo uniforme E lungo x (gauge del potenziale vettore, k_x spostato a ogni tick), theta = 0,1, "
                    "E = E_c = m^2 e 2 E_c; si spazza k_x su 2 pi per TUTTI i k trasversi della zona; si conta la "
                    "probabilita' che uno stato delle bande E < 0 finisca nelle bande E > 0, sommata su tutta la zona "
                    "(nessun proiettore di valle) e divisa per 8 x la formula di Landau-Zener-Schwinger di UNA specie "
                    "(2 spin x E exp(-pi m^2/E) per attraversamento)",
            "prediction": {"F0": "R ~ 90-100 (4 valli massive + 4 coni pi senza massa, REGOLA_LOCALE §12: 96)",
                           "PA": "R ~ 8 (4 valli a 0 + 4 valli pi, tutte di massa m)",
                           "PAW": "R = 1,0 +- 0,1 (una sola specie)"},
            "pass": "PAW: |R - 1| <= 0,1 a entrambi i campi",
        },
        "E7_battimento_1s": {
            "what": "Z = 60, theta = 0,5 (come il Laboratorio), reticolo 40^3, nucleo sull'angolo di 8 celle; stato "
                    "iniziale = 1s di Dirac del continuo (formula: solo stato di prova, non ingresso della regola); "
                    "evoluzione 1600 periodi; p_exc(t) = 1 - |<psi0|psi(t)>|^2; spettro dell'autocorrelazione",
            "prediction": {"F0": "battimento: max p_exc > 0,5; >= 3 righe con peso >= 2 % fra 0,75 e 0,95",
                           "PAW": "max p_exc <= 0,15; una riga (o un doppietto di spin con split <= 5e-3) con peso "
                                  ">= 80 %"},
            "pass": "PAW: max p_exc <= 0,15 e peso della riga principale (o doppietto) >= 0,8",
        },
    },
    "free_parameters": {
        "PAW": "1 regolatore h = 0,3 (Wilson), nessun parametro fisico nuovo; ingressi c, hbar, m, alpha, Z",
        "B2W": "1 regolatore h = 0,3",
        "PA": "0",
    },
    "scripts": ["lab/senzadoppi_r1.py (E1-E3)", "lab/senzadoppi_r2.py (E4-E5)", "lab/senzadoppi_r3.py (E6)",
                "lab/senzadoppi_r4.py (E7)", "lab/senzadoppi_figure.py"],
}
os.makedirs(os.path.join(os.path.dirname(__file__), "out"), exist_ok=True)
with open(os.path.join(os.path.dirname(__file__), "out", "senzadoppi_prereg.json"), "w") as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
print("written lab/out/senzadoppi_prereg.json")

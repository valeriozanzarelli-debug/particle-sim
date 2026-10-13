"""rifl_prereg.py - pre-registrazione dello studio RIFLUSSO (quantum backflow) sulla regola zig-zag 1+1D.
Scritta e committata PRIMA delle corse di produzione. Genera lab/out/rifl_prereg.json.

Note di sviluppo gia' viste prima di questa pre-registrazione (dichiarate):
 - il nucleo (rifl_core.walk_kernel) coincide con la regola vera (zz_core.Walk 'dirac1') a 1e-16 su 3 anelli;
 - NR continuo con taglio U: lambda(5)=0.029457, (10)=0.033773, (20)=0.036092, (40)=0.037268, (80)=0.037858
   (cioe' c_bm - lambda ~ 0.0475/U): per c_bm serve un'estrapolazione in U;
 - Dirac continuo eps=0.1 con taglio p<=30: 0.036753 (tabella di Ashfaque et al.: 0.03686);
 - una prova della regola: th=0.02, N=20000, k<0.02 -> 0.034493.
Nessun altro numero della regola e' stato guardato.
"""
import json, os

C_BM = 0.0384517
ASHFAQUE_T1 = {"0.10": 0.03686, "0.50": 0.03088, "0.80": 0.02722, "1.00": 0.02498,
               "1.60": 0.01947, "2.00": 0.01660, "2.50": 0.01372}
TB_MAX = 0.0764734   # Torres Arvizu, Ortega, Larralde, arXiv:2511.16867 (tight-binding a tempo continuo)

prereg = {
    "studio": "RIFLUSSO - riflusso quantistico (quantum backflow) sulla regola zig-zag 1+1D",
    "prefisso": "rifl_",
    "data": "2026-10-13",
    "regola": "zz_core.Walk kind='dirac1' (corri a c, poi ruota di theta = m c a/hbar); niente potenziale",
    "osservabile": ("Delta_N = sum_{n<N} (|L_n(1)|^2 - |R_n(0)|^2) = probabilita' che attraversa il taglio (0|1) "
                    "verso sinistra in N tick; massimo su stati della banda positiva con k in una finestra W"),
    "metodo": ("autovalore massimo del nucleo K(k,k') (Gauss-Legendre composta in k), controllato contro "
               "la regola vera su anelli; continuo NR (Bracken-Melloy) e Dirac (Melloy-Bracken) con lo stesso codice"),
    "unita": "hbar = c = m = 1; regola: theta = a, N = T/a, eps = sqrt(4/(theta N))",
    "parametri_liberi": 0,
    "ingressi": ["c", "hbar", "massa (solo tramite theta)"],
    "riferimenti_solo_esame": {
        "c_bm": C_BM, "fonte_c_bm": "Penz, Gruebl, Kreidl, Wagner, J. Phys. A 39, 423 (2006)",
        "Ashfaque_Lynch_Strange_2019_Tab1_eps_to_c": ASHFAQUE_T1,
        "fonte_tab": "Ashfaque, Lynch, Strange, Phys. Scr. 94, 125107 (2019), arXiv:1903.07618 (+-2 sull'ultima cifra)",
        "tight_binding_max": TB_MAX,
    },
    "convergenza_numerica": ("ogni numero: risoluzione in k raddoppiata, variazione relativa < 1e-4; tagli in "
                             "impulso estrapolati (fit lambda_inf - A/P - B/P^2) se servono"),
    "esami": [
        {"id": "E0", "cosa": "risolutore NR continuo: estrapolazione U->inf di lambda_NR(U), U = 10..160",
         "soglia": "|c - 0.0384517| <= 2e-6", "tipo": "validazione del codice"},
        {"id": "E1", "cosa": "Dirac continuo (stesso codice) contro la Tabella 1 di Ashfaque et al. 2019",
         "configurazioni": list(ASHFAQUE_T1.keys()),
         "soglia": "PASS se |nostro - tabella| <= 1e-3 per tutti gli eps; 'stretto' se <= 5e-5 (informativo)",
         "tipo": "validazione (e verifica indipendente della letteratura)"},
        {"id": "E2", "cosa": ("limite non relativistico DELLA REGOLA: p_max/m = 0.05, U = p_max sqrt(T/2m) in "
                              "{10,20,40,80}; theta = a in {0.2,0.1,0.05,0.025}; estrapolazione a->0 poi U->inf"),
         "soglia": "|c_regola - c_bm| / c_bm <= 1%; e a ogni U |regola(a->0) - NR(U)|/NR(U) <= 1%",
         "tipo": "validazione"},
        {"id": "E3", "cosa": ("la regola -> Dirac continuo per a->0 a eps fisso, finestra fisica k in (0, pi/2) "
                              "(cioe' tutti gli impulsi positivi quando a->0)"),
         "configurazioni": {"eps": [0.5, 1.0, 2.5, 5.0, 10.0], "theta": [0.4, 0.2, 0.1, 0.05]},
         "soglia": "|regola(a->0) - Dirac(eps)| / Dirac(eps) <= 1% per ogni eps; scarti che calano con a",
         "tipo": "validazione"},
        {"id": "E4", "cosa": "nucleo = regola vera (zz_core) su anelli, 3 configurazioni",
         "soglia": "|differenza| <= 1e-12", "tipo": "controllo"},
    ],
    "domande_nuove_con_previsioni": [
        {"id": "Q1", "cosa": "legge ultrarelativistica di Dirac: c_rbf(eps) per eps fino a 100",
         "previsione": ("c_rbf -> 0 per eps -> inf (massless 1+1D: banda positiva = solo R, niente riflusso). "
                        "Legge di potenza c ~ A eps^-gamma con gamma = 2 (c ~ m T) o 4 (c ~ (m T)^2); "
                        "il fit della letteratura exp(-4eps/9 (1-4 alpha eps)) non puo' valere per eps > 1/(4 alpha) = 34"),
         "criterio": "gamma da fit su eps in [20,100], stabile entro 0.05 tra [20,50] e [50,100]"},
        {"id": "Q2", "cosa": ("la regola a theta finito (reticolo preso sul serio): Lambda(theta, N) nella finestra "
                              "fisica k<pi/2 e nella banda intera k<pi (con il doppione vicino a k = pi)"),
         "configurazioni": {"theta": [0.01, 0.03, 0.1, 0.3, 0.6, 1.0, 1.3, 1.5], "N": "1 .. 1e4 (log)"},
         "previsioni": ["Lambda_fisica <= c_bm per ogni theta, N (probabilita' soggettiva 60%)",
                        "Lambda_banda_intera > c_bm per qualche theta, N (70%), come il tight-binding (0.0765)",
                        "per un tick (N=1) il riflusso e' finito (nel continuo il flusso istantaneo e' illimitato)"],
         "criterio": "sup trovato su griglia + raffinamento locale; dichiarare dove e' raggiunto"},
        {"id": "Q3", "cosa": "doppione: |Lambda_intera - Lambda_fisica| -> 0 per a -> 0 a eps fisso?",
         "previsione": "si' (i settori k~0 e k~pi si disaccoppiano), probabilita' 50%: domanda aperta per noi"},
    ],
    "cosa_conta_come_scoperta": ("numero non presente nella letteratura trovata, convergente (risoluzione e taglio), "
                                 "e per la regola stabile su >= 3 valori di a. Proprieta' della regola a theta finito = "
                                 "proprieta' del modello (come il tight-binding), NON fisica, a meno che sopravviva per a->0"),
    "bonus": "riflusso radiale in stati legati (RadialWalk): solo se resta tempo; nessuna soglia",
}

os.makedirs(os.path.join(os.path.dirname(__file__), "out"), exist_ok=True)
out = os.path.join(os.path.dirname(__file__), "out", "rifl_prereg.json")
with open(out, "w") as f:
    json.dump(prereg, f, indent=1, ensure_ascii=False)
print("scritto", out)

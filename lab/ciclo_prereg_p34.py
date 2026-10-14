"""ciclo_prereg_p34.py - phase 3 (DRAFT: development only, thresholds NOT yet fixed, no production run) and phase 4
(PLAN only) of the CICLO_LIVELLI pre-registration. Imported by ciclo_prereg.py. Development numbers seen so far are in
lab/out/ciclo_p3_dev.json (declared). Before any phase-3 production run, the thresholds marked TBD are fixed in a
committed amendment of this file, BEFORE the run."""

PHASES = dict(
    fase3=dict(
        status="BOZZA (PARZIALE): solo sviluppo; soglie TBD; nessuna corsa di produzione",
        title="Fase 3: controllo al volo e raffinamento locale",
        date="2026-10-14",
        development_seen=dict(
            block_glue="variabili a blocchi netti incollate con Galerkin alla base: B = 9 diagonale 2.4e-3 contro "
                       "3.9e-4 della regola imparata (6x troppo rigida); 10 celle raffinate -> stati spuri a -6.7e-4 "
                       "(20x il legame dell'atomo, -3.25e-5): FALLISCE",
            smoothed="prolungatore ad aggregazione liscia (blocchi lisciati s volte con la regola del livello, Loewdin): "
                     "s = 3, B = 9: regola di Galerkin / regola imparata = 1.003-1.013 nella finestra lenta, stati "
                     "legati entro 0.73 % della base; s = 1: 37 % e 25 % (non va)"),
        config="base 2187 celle (1.6 lambda_C), atomo Coulomb morbido di fase 2; livello grosso B = 9 (243 celle); "
               "base mista = funzioni lisce mascherate (celle grosse) + celle base (raffinate), ortonormale; "
               "Crank-Nicolson; disturbo = campo locale oscillante V_d = A exp(-(x-x_d)^2/2 sigma^2) sin(Omega t) "
               "inviluppo sin^2, sigma = 2 celle base (sotto la cella grossa), Omega = energia per k = 0.5 + legame; "
               "controllo ogni DT ticks su toppe di 3 celle grosse + tampone 27 celle base; rifinire + ripetere "
               "l'intervallo (rollback); solo raffinamento",
        exams_planned=dict(
            B1="quiete (A = 0): zero raffinamenti; fedelta' finale vs base >= 0.99",
            B2="impulso sull'atomo: celle raffinate tutte dentro il cono |x - x_d| <= 3 sigma + v_max t; frazione "
               "raffinata finale <= 0.5",
            B3="impulso: probabilita' di ionizzazione entro 5 % della verita' (base intera), fedelta' >= 0.99; "
               "il livello grosso da solo sbaglia di > 30 % (previsione)",
            B4="conservazione: |1 - norma| <= 1e-9 a parte le proiezioni di raffinamento (perdita totale <= 1e-4); "
               "energia dopo l'impulso costante fra due raffinamenti entro 1e-6 relativo",
            B5="impulso lontano dall'atomo (niente ampiezza): raffinamenti ~ 0 (rifinisce dove c'e' errore, non dove "
               "c'e' campo)",
            thresholds="EPS_REF, DT_CHECK, A: TBD dallo sviluppo, dichiarati prima della produzione")),
    fase4=dict(
        status="PIANO (non eseguibile ora): docs/MOLECOLE.md s.8.2 dice che l'acqua non c'e' ancora (problema del "
               "segno a 5 elettroni per spin). Si esegue solo quando esiste H2O (e il suo dimero) dalla regola a molti "
               "elettroni",
        title="Fase 4: congelare / irraggiare una regione della goccia (prova nei due sensi)",
        plan="goccia di N molecole al livello 'molecola' (agglomerato delle fasi 1-3); una regione scaldata da un "
             "impulso (irraggiare) e una raffreddata (congelare); il controllo al volo deve rifinire SOLO le regioni "
             "disturbate e tornare grosso quando l'errore rientra (de-raffinamento, nuovo rispetto alla fase 3); "
             "rumore di fluttuazione-dissipazione al livello grosso misurato dalle toppe fini (passo 3 del ciclo)",
        exams_planned=dict(
            C1="conservazione: energia totale (piu' lavoro dell'impulso) entro 1e-3 relativo; numero di molecole esatto",
            C2="FD: varianza del rumore grosso = 2 T gamma (gamma misurato dalle toppe) entro 10 %",
            C3="due sensi: la regione congelata dal livello grosso e quella dalle toppe fini hanno la stessa "
               "temperatura finale entro 5 %; la frazione rifinita torna < 10 % dopo il rilassamento",
            C4="esame (solo esame): densita' e calore latente contro l'acqua vera, scarto dichiarato"),
        debts="potenziale fra molecole dal livello sotto (non ancora esistente); termostato = nessuno (FD emergente)"),
)

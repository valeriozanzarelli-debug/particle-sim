"""ciclo_prereg_p34.py - phase 3 (pre-registered 2026-11-02, before production; development declared) and phase 4
(PLAN only) of the CICLO_LIVELLI pre-registration. Imported by ciclo_prereg.py. Development numbers seen so far are in
lab/out/ciclo_p3_dev.json (declared). Before any phase-3 production run, the thresholds marked TBD are fixed in a
committed amendment of this file, BEFORE the run."""

PHASES = dict(
    fase3=dict(
        status="PRE-REGISTRATA (2026-11-02), committata PRIMA delle corse di produzione",
        title="Fase 3: controllo al volo e raffinamento locale",
        date="2026-11-02",
        development_seen=dict(
            block_glue="variabili a blocchi netti incollate con Galerkin alla base: B = 9 diagonale 2.4e-3 contro "
                       "3.9e-4 della regola imparata (6x troppo rigida); 10 celle raffinate -> stati spuri a -6.7e-4 "
                       "(20x il legame dell'atomo, -3.25e-5): FALLISCE",
            smoothed="aggregazione liscia s = 3: regola di Galerkin / regola imparata = 1.003-1.013, stati legati entro "
                     "0.73 % della base",
            loss_fix="prima versione: le funzioni lisce delle celle raffinate perdevano la coda fuori -> norma persa "
                     "7e-3 al raffinamento; corretto (tutte le funzioni mascherate + ortogonalizzazione canonica, "
                     "spazio vecchio dentro il nuovo): perdita 1e-9",
            scenario_P_seen="lo scenario P (impulso al centro, A = 2e-3) e la quiete Q sono stati GIRATI in sviluppo: "
                            "Q a EPS 0.02: 0 raffinamenti, fondo del controllo 0.0064; P a EPS 0.02: fedelta' min "
                            "0.986, P_out 0.2225 vs 0.2273; P a EPS 0.01: fedelta' 0.9954, P_out 0.2252, 71/243 celle, "
                            "costo 1.28e6 vs 2.62e6; solo-grosso: P_out 1.2e-4 (fedelta' 0.75); B = 3 uniforme: P_out "
                            "0.2223 ma fedelta' 0.715. P2 e F NON sono stati girati (prova cieca)"),
        config=dict(
            base="2187 celle (1.6 lambda_C), stencil della fase 2, atomo Coulomb morbido (g 0.015, s 10 lambda_C)",
            coarse="B = 9 (243 celle): variabili automatiche della catena b = 3 lisciate s = 3; regola = P^dag H P",
            mixed="celle grosse: funzioni lisce mascherate fuori dai siti raffinati (tutte, comprese le code delle "
                  "raffinate), ortogonalizzazione canonica (taglio 1e-8); celle raffinate: i 9 siti base; "
                  "Crank-Nicolson dt = 25 tick; partenza = stato fondamentale del livello grosso, prolungato",
            check="ogni 16 passi (400 tick); toppe di 3 celle grosse; finestra base = toppa +- 27 siti; "
                  "e = ||psi_base - psi_mix||_toppa; si salta la toppa se 2 sqrt(peso finestra) < EPS_REF; "
                  "e > EPS_REF = 0.01 (1.5 x il fondo di quiete 0.0064) -> raffina toppa + 1 cella per lato, rollback",
            drive="V_d = A exp(-(x - x_d)^2 / 2 sigma^2) sin^2(pi t / Tp) sin(Omega t), sigma = 2 siti, Tp = 15700, T = 30000",
            scenarios=dict(Q="A = 0", P="A = 2e-3, x_d = centro, Omega = 2.4e-3 (k ~ 0.5)",
                           P2="A = 3e-3, x_d = centro + 12, Omega = G_base(0.7) - E0 (cieco)",
                           F="A = 2e-3, x_d = centro + 600 (dove l'elettrone non c'e'), Omega = 2.4e-3 (cieco)"),
            kinds="truth = base ovunque; mix = adattivo; coarse = livello grosso senza controllo; b3 = livello B = 3 "
                  "uniforme (729 celle, liscio s = 3)",
            P_out="peso a |x - centro| > 120 siti (7 raggi rms del fondamentale)",
            cost="gradi di liberta' x passi (propagazione) + gradi di liberta' x passi delle finestre di controllo"),
        exams=dict(
            B1=dict(what="Q: zero raffinamenti; fedelta' minima vs base >= 0.999", prediction="PASS"),
            B2=dict(what="P e P2: ogni cella raffinata ha il centro entro |x - x_d| <= 36 + 0.029 t (toppa + margine + "
                         "velocita' massima della banda base, t = istante del raffinamento); frazione raffinata finale "
                         "<= 0.5", prediction="PASS"),
            B3=dict(what="P e P2: |P_out - verita'| / verita' <= 5 % a T e fedelta' minima >= 0.99; il livello grosso "
                         "da solo: errore relativo su P_out > 30 %", prediction="PASS (P2 cieco: incerto)"),
            B4=dict(what="conservazione (tutti gli adattivi): deriva della norma <= 1e-8 compresi i raffinamenti; "
                         "dopo l'impulso energia propria (max - min) / |E| <= 1e-6 compresi i salti ai raffinamenti",
                    prediction="PASS"),
            B5=dict(what="F: zero raffinamenti e fedelta' minima >= 0.999 (si raffina dove c'e' errore, non dove c'e' "
                         "campo)", prediction="PASS"),
            B6=dict(what="COSTO a pari errore (P e P2): risparmio = costo del livello uniforme piu' economico (B = 9, "
                         "B = 3, base) con fedelta' minima >= quella dell'adattivo, diviso il costo dell'adattivo "
                         "(controlli compresi); PASS se >= 1.5", prediction="PASS, ~2 (in 1-D l'atomo occupa una "
                         "frazione grande della scatola: in 3-D il risparmio cresce col volume vuoto)"),
        )),
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

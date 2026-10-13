"""motore_prereg.py - PRE-REGISTRATION of the study "Motore ottimizzato" (physics side of the multilevel engine).
Written and committed BEFORE the production runs (motore_m1.py ... motore_m5.py). Development explorations done
before this file (declared): (i) Fourier-exact light evolution vs real-space ticks (2e-15); (ii) error model vs a few
packets: the first try used a too small periodic box (packet tails wrapped, y-drift off by 45 %); with a box of
+-4.5 sigma the Bloch model matched; (iii) the causal-cone freeze on one 128^2 case (bit-identical); (iv) the
electron pipeline on the 1-D learned levels (one case each B = 16, 192, 384); (v) the transverse curvature of the
light dispersion (diffraction factor) at N = 8..32. The numbers in 'predictions' below come from the Bloch
dispersion of the rule (zz_core) and are the PREDICTIONS the measurements are judged against.
Run: python lab/motore_prereg.py  ->  lab/out/motore_prereg.json"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_core as mc

TIERS = {
    "ultra": dict(hardware="scheda da 16 GB (obiettivo)", light_N=32, electron_B=16, tol=0.02),
    "intel": dict(hardware="Intel UHD del portatile (unica GPU: la GTX 1650 Ti e' rotta)", light_N=12, electron_B=192,
                  tol=0.10),
    "sotto_pavimento": dict(hardware="sonda: un livello sotto il pavimento dinamico", light_N=8, electron_B=384,
                            tol=0.10),
}
LADDER = [8, 12, 16, 24, 32, 48, 64]
ANGLES = [0.0, 22.5, 30.0, 45.0]

pred = {}
for N in LADDER:
    row = {}
    for a in ANGLES:
        dw, dvg, dv = mc.level_errors(N, a)
        k = 2 * np.pi / N
        kh = np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0])
        pe = np.array([-np.sin(np.radians(a)), np.cos(np.radians(a)), 0])
        h = 1e-3 * k
        w = lambda q: mc.omega_plus(k * kh + q * pe)
        c2 = (w(h) - 2 * w(0) + w(-h)) / h ** 2 * k
        row[str(a)] = dict(domega_over_omega=dw, speed_error=dv, walkoff_rad=float(np.arctan2(dvg[1] + np.sin(np.radians(a)), dvg[0] + np.cos(np.radians(a))) - np.radians(a)),
                           diffraction_factor=c2)
    pred[str(N)] = row

# predicted level choice of the scene B1 (the engine's own calculation, made here once)
FOOT = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
rng = np.random.default_rng(20261013)
ANG_B = [float(x) for x in rng.choice([0.0, 22.5, 30.0, 45.0, 60.0, 15.0], size=len(FOOT))]
T_H = 100.0


def choose(p, ang, floor, T=T_H, safety=0.5):
    for N in LADDER:
        if N < floor:
            continue
        dx, dphi = mc.predicted_drift(N, ang, T)
        e, *_ = mc.screen_error_px(dx, dphi, p)
        if e <= safety:
            return N, e
    return None, None


sel = {}
for tier, floor in (("ultra", 32), ("intel", 12)):
    rows = []
    for p, a in zip(FOOT, ANG_B):
        N, e = choose(p, a, floor)
        rows.append(dict(footprint_lambda_per_px=p, angle=a, N=N, predicted_px=e))
    ok = [r for r in rows if r["N"] is not None]
    cost = sum(r["N"] ** 3 for r in ok)
    sel[tier] = dict(rows=rows, flagged=[r["footprint_lambda_per_px"] for r in rows if r["N"] is None],
                     cost_rel_uniform64_unflagged=cost / (len(ok) * 64 ** 3))

prereg = dict(
    study="motore_ottimizzato", date="2026-10-13", branch_note="commit before the production runs",
    perimeter=("Solo LATO FISICO dello schema multilivello (cambio di perimetro del coordinatore): scelta del livello per "
               "regione entro un limite di errore fisico (anche proiettato sullo schermo come criterio), transizioni tra "
               "livelli che conservano carica ed energia, passo di tempo locale. Il disegno (splat, preset, GPU) e' del "
               "thread 'Modello standard a palline'."),
    principle=("Valerio: la mela in 4K e in 480p marciscono tutte e due. Il budget di errore scala con l'hardware; le "
               "DINAMICHE devono restare: leggi di conservazione esatte, frequenze e tempi caratteristici entro la tolleranza "
               "del livello, stessa sequenza di eventi, stesse simmetrie. Possono peggiorare: nitidezza (cella), dettagli "
               "fini (onde piu' corte di ~12 celle, tremolio di Dirac)."),
    tiers=TIERS,
    definitions=dict(
        conservation="norma (carica) e energia del campo sum|F|^2 conservate entro 1e-10 dentro un livello uniforme",
        characteristic="velocita' del centro d'energia, frequenza (fase), apertura per diffrazione, energia cinetica; errore relativo <= tol del livello",
        sequence="due eventi separati nel riferimento da piu' di tol*t restano nello stesso ordine (equivale a errori di velocita' < tol/2)",
        symmetry="purezza di elicita' (luce) e isotropia: max-min della velocita' sugli angoli <= tol",
        may_degrade="nitidezza (cella), contenuto con meno di ~12 celle per lambda, Zitterbewegung (periodo pi hbar/mc^2 < passo grosso)"),
    exams=dict(
        A1="modello d'errore spettrale (Bloch pesato dallo spettro del pacchetto): deriva dell'inviluppo entro 5% + 0.01 lambda e fase entro 0.05 rad (mod 2pi) contro la misura; 16 casi N in {8,12,16,32} x angoli {0,22.5,30,45}, T = 50 lambda, sigma = 3 lambda, scatola +-4.5 sigma. PASS se 16/16",
        A2="modello a un solo k0 (quello economico per il motore): stesse soglie; riportato, PASS se >= 14/16",
        B1="scena 'errore sullo schermo': 8 pacchetti con impronta {0.1..20} lambda/pixel, orizzonte 100 lambda, scala N {8..64}, pavimento del livello (32 Ultra, 12 Intel), margine 0.5 px: errore misurato <= 1 px su tutti i pacchetti non segnalati; costo (celle x tick, 2-D, prop. N^3) sui pacchetti NON segnalati ridotto >= 4x contro gli stessi a N = 64 (Intel) e riportato (Ultra). I pacchetti segnalati (nessun livello della scala basta) sono fuori scala anche per il riferimento: si dichiarano, e la loro misura deve confermare che servirebbe di piu' (errore misurato a N = 64 > 0.5 px). Correzione fatta prima delle misure: la prima stesura contava i segnalati come N = 64 nel costo",
        B2="la memoria dell'errore: camera che si avvicina (impronta 5 -> 0.1 lambda/px in 100 lambda), pacchetto a 30 gradi. Politica ingenua (sceglie il livello col pixel di adesso, ogni 10 lambda, cambio di livello spettrale): errore finale > 1 px (FAIL atteso, e' la dimostrazione). Politica col 'passaporto' (errore accumulato + impronta minima raggiungibile): <= 1 px",
        C1="cono causale: celle fuori dal cono all'indietro della regione guardata congelate; differenza nella regione guardata <= 1e-13 (relativa) per luce 2-D e Dirac 2-D con Coulomb Z = 60 (n = 192, T = 64 tick, finestra 48x48); un cono di 1 cella piu' stretto DEVE sbagliare (> 1e-6); risparmio misurato = conteggio analitico",
        C2="orologio pigro: recupero di m = 1024 tick con W(k)^m in Fourier contro tick per tick: <= 1e-11 relativo (luce 2-D, n = 128)",
        C3="regola agglomerata (1-D, livelli B = 16, 64, 192): coda fuori dal cono dello stencil dopo s passi (descrittivo, senza soglia)",
        D1="transizioni di livello per la luce (N -> N/2 -> N), metodo spettrale contro cubico: |dE/E| al cambio <= 1e-6 (spettrale, N grosso >= 12), salto del centro <= 1e-3 lambda, parte su banda sbagliata del reticolo grosso <= 1e-3 (N grosso >= 12); cubico riportato",
        D2="transizioni per l'elettrone (fine -> livello B -> fine) col blocco dell'agglomeratore e prolungamento spettrale (deconvoluzione del blocco): norma del giro <= 1e-10; prolungamento a gradini riportato",
        EL="dinamica della luce Ultra (32) / Intel (12) / sotto (8) contro Maxwell continuo: EL1 velocita' agli angoli {0,22.5,30,45}; EL2 frequenza; EL3 apertura per diffrazione (fascio w0 = 2 lambda, 0 e 45 gradi); EL4 cucitura tra due livelli (mondo a blocchi di reticolo_core, fine->grosso e grosso->fine, 0 e 30 gradi): T >= 1 - tol/20 (0.999 Ultra, 0.995 Intel), energia persa <= tol/20, elicita' sbagliata <= tol/20, velocita' dopo la cucitura entro tol; EL5 energia nel livello uniforme <= 1e-10; EL6 isotropia (max-min velocita') <= tol e sequenza. PASS per livello se tutte le voci passano",
        EE="dinamica dell'elettrone (livelli 1-D imparati, scala atomica: k0 = 1/a0, sigma = a0/2) Ultra B = 16 / Intel B = 192 / sotto B = 384 contro la regola fine esatta: EE1 velocita' di gruppo; EE2 larghezza al tempo di allargamento 2 m sigma^2; EE3 energie cinetiche a k = 0.5, 1, 2 /a0 (fase dei modi); EE4 norma <= 1e-10; entro tol del livello. EE5 Zitterbewegung (descrittivo: dettaglio che puo' sparire)",
        F="budget (calcolo, non esame): costo e memoria per livello sulla Intel UHD (42 ms per passo con ~33 800 palline misurati nel Laboratorio) e su una scheda da 16 GB; richieste di ponte proposte per misurarli"),
    predictions=dict(
        bloch=pred,
        EL_expected=("Ultra 32: velocita' max -0.49% (45 gradi), frequenza -0.16%, diffrazione 0.987 (0 gradi) -> PASS; "
                     "Intel 12: velocita' -3.6%, frequenza -1.2%, diffrazione 0.907 (-9.3%, al limite) -> PASS; "
                     "sotto 8: velocita' -8.4%, diffrazione 0.785 (-21%) -> FAIL su EL3 (e' la prova che il pavimento esiste)"),
        EE_expected="i livelli imparati hanno errori 1e-4..1e-3 sulla dinamica a scala atomica: PASS a tutti e tre i livelli; ZB visibile solo a B = 16 (passo 0.54 < pi/2), aliasata/mediata a B = 192, 384",
        B1_selection=sel,
        B2_expected="ingenua > 1 px (FAIL atteso); passaporto <= 1 px"),
    configs=dict(ladder=LADDER, angles=ANGLES, B1_footprints=FOOT, B1_angles=ANG_B, horizon_lambda=T_H, safety_px=0.5,
                 electron=dict(theta="a0/4096 (agg_zz_1d)", fine_cells=196608, box_a0="~48", k0="1/a0", sigma="a0/2"),
                 seed=20261013, threads=1),
)
out = os.path.join(mc.OUT, "motore_prereg.json")
json.dump(prereg, open(out, "w"), indent=1)
print("written", out)
print(json.dumps(sel, indent=1))

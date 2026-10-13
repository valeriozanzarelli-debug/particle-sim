"""tunnel_prereg.py - pre-registrazione dello studio TEMPO DI TUNNEL (scritta PRIMA delle corse di produzione).
Le previsioni numeriche vengono dal limite continuo della regola (Dirac 1+1 D) e da Schroedinger (limite non
relativistico, Buttiker 1983): sono calcoli su carta/macchina della teoria, non corse della regola.
Uscita: lab/out/tunnel_prereg.json"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunnel_core import times_dirac, times_schrod, schrod_steps

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

E = 2.0
V0S = [1.2, 1.6, 2.0, 2.4, 2.8]
LS = [1, 2, 4, 8, 12, 16]
A_LIST = [0.1, 0.05, 0.025, 0.0125]           # passi reticolari (theta = m a), estrapolazione a -> 0

pred = {}
for V0 in V0S:
    eps = E - V0
    kap = np.sqrt(1 - eps ** 2)
    rows = []
    for L in LS:
        r = times_dirac(E, 1.0, V0, L)
        rows.append(dict(L=L, T=r["T"], tphi=r["tphi"], tau_y=r["TV"].real, tau_z=-r["TV"].imag,
                         tau_my=r["Tm"].real, tau_mz=-r["Tm"].imag))
    pred[f"V0={V0}"] = dict(eps=eps, kappa=kap, rows=rows,
                            ratio_my_y_L16=rows[-1]["tau_my"] / rows[-1]["tau_y"], eps_over_m=eps,
                            slope_mz=(rows[-1]["tau_mz"] - rows[-3]["tau_mz"]) / 8, m_over_kappa=1 / kap,
                            slope_z=(rows[-1]["tau_z"] - rows[-3]["tau_z"]) / 8, eps_over_kappa=eps / kap)

# pacchetti
def pk(Ei, V0, L, sig):
    r = times_dirac(Ei, 1.0, V0, L)
    v = r["vg"]
    x0 = -7 * sig - 10
    tf = (abs(x0) + L + 4 * sig) / v
    adv = L - v * r["tphi"]
    return dict(E=Ei, V0=V0, L=L, sigma=sig, x0=x0, t_final=tf, v=v, tphi=r["tphi"], T=r["T"],
                peak_advance_pred=adv, margin_vs_light_pred=adv - (1 - v) * tf,
                superluminal_peak_pred=bool(adv - (1 - v) * tf > 0))
packets = dict(PK1=pk(5.0, 5.0, 10.0, 20.0), PK2=pk(2.0, 2.0, 8.0, 20.0))

# limite non relativistico (Buttiker), rettangolo con kappa0 L = 3
nr = []
for x in (0.2, 0.4, 0.6, 0.8):
    V0 = 1.0
    L = 3 / np.sqrt(2 * V0)
    r = times_schrod(x * V0, V0, L)
    nr.append(dict(E_over_V0=x, tau_y_V0=r["TV"].real, tau_z_V0=-r["TV"].imag, tphi_V0=r["tphi"], T=r["T"]))

# Ramos 2020: barriera gaussiana V0 exp(-2 x^2 / w^2), forma tau_y(E/V0) (Schroedinger, dimensionless)
def gauss_tauy(Ek, V0, w, M=1.0, n=600, h=1e-6):
    xs = np.linspace(-3 * w, 3 * w, n + 1)
    xc = 0.5 * (xs[1:] + xs[:-1])
    dx = xs[1] - xs[0]
    prof = np.exp(-2 * xc ** 2 / w ** 2)
    tp = schrod_steps(Ek, V0 * (1 + h) * prof, dx, M)
    tm = schrod_steps(Ek, V0 * (1 - h) * prof, dx, M)
    TV = 1j * np.log(tp / tm) / (2 * h * V0)
    return TV.real, -TV.imag, abs(schrod_steps(Ek, V0 * prof, dx, M)) ** 2
ramos = []
for kw in (1.0, 2.0):
    w = kw / np.sqrt(2.0)       # kappa0 w = kw con V0 = 1, M = 1
    for x in (0.3, 0.5, 0.7, 0.9):
        ty, tz, T = gauss_tauy(x, 1.0, w)
        ramos.append(dict(kappa0_w=kw, E_over_V0=x, tau_y_V0=ty, tau_z_V0=tz, T=T))

prereg = dict(
    studio="TEMPO DI TUNNEL e TEMPO DI ARRIVO sulla regola zig-zag 1+1 D",
    data="2026-10-13", autore="Claude (Opus 5.5) per Valerio",
    regola="lab/zz_core.py dirac1 (corri a c, poi ruota th = m a, fase e^{-iVa}); variante: th e V per sito",
    parametri_liberi=0,
    impostazioni_esperimento="E, V0, L (barriera rettangolare o gaussiana): condizioni dell'esperimento, non manopole",
    passi_reticolari=A_LIST,
    previsioni_continuo_dirac=pred,
    pacchetti=packets,
    limite_NR_buttiker=nr,
    ramos_forma=ramos,
    esami={
        "E1_cono_di_luce": "pacchetto a supporto compatto: dopo n tick |psi|^2 == 0 esattamente oltre x_R0 + n a e prima di x_L0 - n a (con e senza barriera, tutte le copie); l'ultima cella non nulla avanza di 1 cella per tick: v_fronte = c con errore < 1e-10 (atteso 0)",
        "E2_hartman_continuo": "tau_phi della regola (stazionario esatto, a = 0.1, 0.05, 0.025, 0.0125, Richardson in a^2) = Dirac continuo entro 1e-3 relativo per tutti V0 in V0S e L in LS",
        "E3_saturazione": "regola a = 0.025: |tau_phi(16) - tau_phi(8)| / tau_phi(16) < 1e-3 per tutti V0 (Hartman: satura)",
        "E4_velocita_apparente": "L / tau_phi(L=4) > c per tutti V0 (fase superluminale presente nella regola)",
        "E5_picco_pacchetto": "PK1, PK2: spostamento del picco trasmesso rispetto al pacchetto libero = L - v tau_phi entro 0.02 L (a -> 0); segno del margine rispetto al raggio di luce dal picco incidente come previsto (PK1 positivo = picco 'prima della luce', PK2 negativo)",
        "E6_larmor_dwell": "stazionario: tau_y(trasmessi) = tau_y(riflessi) = dwell entro 1e-3 (a -> 0); pacchetto: P_T tau_y^T + P_R tau_y^R = sum_n P_barriera a entro 1e-4 relativo",
        "E7_larmor_pacchetto": "orologio vero (due copie V0 +- delta/2, delta = 1e-4): tau_y = media stazionaria pesata entro 1e-3, tau_z entro 1e-2",
        "E8_larmor_continuo": "tau_y, tau_z della regola (a -> 0) = Dirac continuo entro 1e-3 relativo (V0S, L = 1, 4, 8)",
        "E9_limite_NR": "V0/m = 0.04, 0.02, 0.01, kappa0 L = 3: tau_y V0, tau_z V0 della regola (a -> 0, poi lineare in V0/m -> 0) = Buttiker (Schroedinger) entro 1%",
        "E10_ramos_forma": "barriera gaussiana, V0/m = 0.01: tau_y della regola = Schroedinger entro 1% (E/V0 = 0.3..0.9, kappa0 w = 1 e 2); forma prevista (Schroedinger): per kappa0 w = 2 tau_y cresce con E/V0 (energie piu' basse = tempi piu' corti, come Ramos), per kappa0 w = 1 non monotono (massimo vicino a E/V0 ~ 0.5-0.7); valori in ms per Rb87 con w = 1.3 um riportati senza soglia (altezza della barriera di Ramos non nota)",
        "E11a_rotazioni_libere": "zona senza barriera: orologio delle rotazioni T_m = i dlnt/dm = tempo proprio L m / k entro 1e-3 (a -> 0); identita' f_w (contatore z) = s c dlnt/dth + s^2 n_w entro 1e-6",
        "E11b_rapporto": "barriera spessa (L = 16), a -> 0: tau_my / tau_y = eps/m = (E - V0)/m entro 1% per tutti V0",
        "E11c_pendenze": "tra L = 8 e 16: d tau_mz/dL = m/kappa e d tau_z/dL = eps/kappa entro 1%",
        "E11d_centro_gap": "V0 = E (eps = 0): tau_mz(L) = L/c entro 1% per L >= 4; |tau_my| < 0.01 tau_y",
        "E11e_segno": "segno di tau_my = segno di eps per L >= 8 (l'orologio delle rotazioni va all'indietro quando V0 > E)",
        "E11f_non_causale_sottile": "previsione (non un vincolo): per barriere sottili tau_mz < L/c (Dirac: L=1, eps=0: 0.851); la regola deve ritrovarlo",
        "E12_bohm": "traiettorie dei surfisti (quantili 1-D = Bohm): ogni traiettoria trasmessa passa >= L/c - 2a nella barriera; P_T tau_B^T + P_R tau_B^R = dwell entro 1e-3",
        "E13_arrivo": "rivelatore assorbente ideale della regola (toglie R oltre X_D: condizione 'solo onda uscente', Tumulka): probabilita' assorbita totale = probabilita' passata oltre X_D senza rivelatore entro 1e-6; tempi medi di arrivo flusso vs assorbitore riportati senza soglia",
    },
    note=["Le previsioni sono del continuo di Dirac: se la regola le ritrova e' VALIDAZIONE (Dirac lo sa gia').",
          "La parte nuova e' l'interpretazione: orologio delle rotazioni = orologio di massa = tempo proprio, e le "
          "sue proprieta' nella barriera (rapporto eps/m, pendenza m/kappa, L/c al centro del gap). Da controllare "
          "in letteratura nel documento.",
          "V0 = 5 e V0 = 2..2.8 superano 2 m c^2: in teoria dei campi la barriera creerebbe coppie (Klein). La "
          "lettura a una particella le ignora: DEBITO/limite dichiarato."],
)
with open(os.path.join(OUT, "tunnel_prereg.json"), "w") as f:
    json.dump(prereg, f, indent=1, default=float)
print(json.dumps(packets, indent=1, default=float))
for k, v in pred.items():
    print(k, "ratio", v["ratio_my_y_L16"], v["eps_over_m"], "slopes", v["slope_mz"], v["m_over_kappa"], v["slope_z"], v["eps_over_kappa"])
print(nr)
print(ramos)

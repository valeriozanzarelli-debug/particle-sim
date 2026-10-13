"""dilat_prereg.py - pre-registrazione dello studio DILATAZIONE QUANTISTICA (prefisso dilat_).
Scrive lab/out/dilat_prereg.json PRIMA delle corse della regola.
Le previsioni sono calcolate qui SOLO col continuo (Dirac a energia positiva = orologio di massa m/E) e con le
formule della letteratura (Smith-Ahmadi 2020 eq. 24-25, Paige-Plato-Kim 2020) piu' la nostra previsione
analitica dell'errore reticolare O(a^2). Nessun numero della regola e' usato.
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dilat_core import (smith_ahmadi_gq, gq_offset_nr, lattice_rate_a2, two_packets, overlap_spec,
                        u_plus_dirac, spinor_state, U_dirac)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "dilat_prereg.json")
m = 1.0
A_LIST = [1 / 4, 1 / 8, 1 / 16, 1 / 32]
LBOX = 1024.0


def pgrid(n=2 ** 16, L=LBOX * 4):
    return 2 * np.pi * np.fft.fftfreq(n, d=L / n)


def rate_cont(g, p):
    w = np.abs(g) ** 2
    return float(np.sum(w * m / np.sqrt(p ** 2 + m ** 2)) / np.sum(w))


def dq_cont(al, phi, p1, p2, Dl, d=0.0):
    p = pgrid()
    g, g1, g2 = two_packets(p, np.cos(al), np.sin(al), phi, p1, p2, Dl, 0.0, d)
    sup = rate_cont(g, p)
    mix = np.cos(al) ** 2 * rate_cont(g1, p) + np.sin(al) ** 2 * rate_cont(g2, p)
    return sup, mix, sup - mix


pre = dict(
    studio="DILATAZIONE QUANTISTICA DEL TEMPO con l'orologio interno della regola (dilat_)",
    data="2026-10-13",
    regola="zz_core.Walk('dirac1') congelata, 1+1 D, libera; th = m a; unita' hbar = c = m = 1",
    orologio=("orologio di massa: tau = -arg<psi_{m-d/2}(t)|psi_{m+d/2}(t)>/d (d -> 0: = contatore di flip "
              "Im f_w = -(sin 2th/2a) tau); d finito = qubit di Zych/Smith-Ahmadi/Paige-Plato-Kim"),
    parametri_liberi=0,
    debiti=["la sonda a due masse (qubit con splitting delta) e' messa a mano come l'orologio di Larmor di TEMPO_TUNNEL",
            "pacchetti iniziali gaussiani di sola energia positiva preparati a mano (proiettore di banda dichiarato)",
            "il 'miscuglio' classico si calcola come media pesata delle due corse separate (non e' una dinamica)"],
    passi_reticolari=A_LIST,
    scatola=LBOX,
    convenzione_pacchetti="g_j(p) = exp(-(p-p_j)^2/(2 Dl^2)) e^{-i p x_j}  (Smith-Ahmadi eq. 19; |g|^2 ha varianza Dl^2/2)",
    esami={},
)
E = pre["esami"]

# ---------------------------------------------------------------- E0 identita' del codice
E["E0"] = dict(descr="regola nel tempo (Walk) = motore spettrale; identita' del contatore di flip",
               config=dict(a=1 / 8, nx=4096, ticks=1600, stati=["due pacchetti p=0.5,0.8 Dl=0.3", "chiralita' R pura (ZB)"]),
               soglie={"|O_tempo - O_spettrale|": 1e-10, "|Re f_w - n sin^2 th|/n": 1e-9,
                       "|Im f_w + (sin2th/2a) tau| / |Im f_w|": 1e-6})

# ---------------------------------------------------------------- E1 dilatazione classica
Dl1 = 0.05
pbar = [0.0, 0.5, 1.0, 2.0, 3.0]
pred1 = {}
p = pgrid()
for pb in pbar:
    g, _, _ = two_packets(p, 1.0, 0.0, 0.0, pb, pb, Dl1)
    r = rate_cont(g, p)
    w = np.abs(g) ** 2 / np.sum(np.abs(g) ** 2)
    C = float(np.sum(w * (m / np.sqrt(p ** 2 + 1)) * lattice_rate_a2(p)) / r)
    pred1[str(pb)] = dict(rate_cont=r, inv_gamma_pbar=1 / np.sqrt(1 + pb ** 2), coeff_a2_previsto=C,
                          rate_regola_prevista={str(a): r * (1 + C * a * a) for a in A_LIST})
E["E1"] = dict(descr="pacchetto stretto in impulso: tau/t -> <m/E> (dilatazione classica), convergenza O(a^2)",
               config=dict(Dl=Dl1, pbar=pbar), previsione=pred1,
               soglie={"E1a |R(a=1/32)-R_cont|/R_cont": 1e-3,
                       "E1b rapporto errori a->a/2 (1/8->1/16->1/32, pbar>0)": [3.5, 4.5],
                       "E1c (R-R_cont)/(a^2 R_cont) a=1/32 contro coeff previsto (pbar>0), errore relativo": 0.05})

# ---------------------------------------------------------------- E2 correzione quantistica
# E2a: non relativistico, contro Smith-Ahmadi eq. 25
Dl2 = 0.02
p1 = 0.02
svals = [0.5, 1.0, 2.0, 3.0, 4.0]
alvals = [np.pi / 8, np.pi / 4, 3 * np.pi / 8]
phivals = [0.0, np.pi / 2, np.pi]
pred2a = []
for s in svals:
    for al in alvals:
        for ph in phivals:
            p2 = p1 + s * Dl2
            sup, mix, dq = dq_cont(al, ph, p1, p2, Dl2)
            pred2a.append(dict(s=s, alpha=al, phi=ph, p1=p1, p2=p2, dq_cont=dq,
                               dq_SA=float(smith_ahmadi_gq(al, ph, p1, p2, Dl2)), rate_mix=mix))
E["E2a"] = dict(descr="correzione quantistica Delta_Q = R_sup - R_mix nel regime non relativistico contro S-A eq.25",
                config=dict(Dl=Dl2, p1=p1, s=svals, alpha=alvals, phi=phivals, a="tutti; esame a 1/32"),
                previsione=pred2a,
                soglie={"|dq_regola - dq_SA| / |dq_SA| (dove |dq_SA| > 1e-9); soglia fissata guardando l'errore NR del solo continuo (max 1.8 %)": 0.03,
                        "|dq_regola(1/32) - dq_cont| / |dq_cont| (dove |dq_cont| > 1e-9)": 1e-3,
                        "|dq_regola| per phi = pi/2": 1e-12})

# E2b: relativistico, contro il continuo esatto
Dl3 = 0.3
pr = 0.5
pred2b = []
for s in [0.5, 1.0, 2.0, 3.0]:
    for al in [np.pi / 8, np.pi / 4]:
        for ph in [0.0, np.pi]:
            p2 = pr + s * Dl3
            sup, mix, dq = dq_cont(al, ph, pr, p2, Dl3)
            pred2b.append(dict(s=s, alpha=al, phi=ph, p1=pr, p2=p2, dq_cont=dq,
                               dq_SA_NR=float(smith_ahmadi_gq(al, ph, pr, p2, Dl3)), rate_sup=sup, rate_mix=mix))
E["E2b"] = dict(descr="regime relativistico (p ~ mc): la regola contro il continuo esatto; S-A (NR) riportato",
                config=dict(Dl=Dl3, p1=pr), previsione=pred2b,
                soglie={"|dq_regola(1/32) - dq_cont| / |dq_cont|": 1e-3,
                        "ordine: rapporto errori 1/16->1/32": [3.0, 5.0]})

# E2c: pacchetti separati nello spazio (estensione nostra della stessa algebra)
dvals_D = [0.0, 0.5, 1.0, 2.0, 3.0]
pred2c = []
for dD in dvals_D:
    d = dD / Dl2
    for ph in [0.0, np.pi]:
        p2 = 0.06
        sup, mix, dq = dq_cont(np.pi / 4, ph, p1, p2, Dl2, d)
        pred2c.append(dict(dDl=dD, d=d, phi=ph, dq_cont=dq, dq_formula_NR=float(gq_offset_nr(np.pi / 4, ph, p1, p2, Dl2, d))))
E["E2c"] = dict(descr="pacchetti separati di d nello spazio: Delta_Q(d) contro la formula NR estesa",
                config=dict(Dl=Dl2, p1=p1, p2=0.06, alpha="pi/4"), previsione=pred2c,
                soglie={"|dq_regola - dq_formula| / max(|dq_formula|, 1e-3 |dq(d=0)|)": 0.03})

E["E2d"] = dict(descr=("persistenza: il ritmo dell'orologio della sovrapposizione resta costante nel tempo anche quando "
                       "i due pacchetti si sono separati nello spazio (|psi(k)|^2 conservato)"),
                config=dict(a=1 / 16, p1=0.5, p2=0.8, Dl=0.3, alpha="pi/4", phi=0, t_max=200, nx=8192),
                previsione="R(finestra) costante; Delta_Q come E2b (s = 1)",
                soglie={"max variazione relativa del ritmo tra finestre di 20": 1e-6})

# ---------------------------------------------------------------- E3 Paige-Plato-Kim (delta finito)
pred3 = []
pp = pgrid()
for pb in [0.05, 0.1, 0.3]:
    for dl in [0.2, 0.5]:
        g, _, _ = two_packets(pp, 1.0, 0.0, 0.0, pb, pb, 0.01)
        w = np.abs(g) ** 2 / np.sum(np.abs(g) ** 2)
        Mp, Mm = m + dl / 2, m - dl / 2
        exact = float(np.sum(w * (np.sqrt(pp ** 2 + Mp ** 2) - np.sqrt(pp ** 2 + Mm ** 2))) / dl)
        ppk = float(np.sum(w * (1 - pp ** 2 / (2 * Mp * Mm))))
        meanm = float(np.sum(w * (1 - pp ** 2 / (2 * m * m))))
        pred3.append(dict(pbar=pb, delta=dl, D_cont_esatto=1 - exact, D_PPK_NR=1 - ppk, D_massa_media_NR=1 - meanm))
E["E3"] = dict(descr=("orologio a qubit con splitting finito delta (spinta in impulso): dilatazione D = 1 - ritmo della "
                      "fase relativa, contro Paige-Plato-Kim D = p^2/(2 M+ M-) (non la massa media p^2/2m^2); fase a t = 2"),
               config=dict(Dl=0.01, t=2.0), previsione=pred3,
               soglie={"|D_regola(1/32) - D_cont| / D_cont": 1e-2,
                       "|D_regola - D_PPK| < 0.2 |D_PPK - D_massa_media| per pbar = 0.05, delta = 0.5": True})

# ---------------------------------------------------------------- E4 visibilita'
tvis = [0, 50, 100, 200, 400, 800]
p = pgrid()
u = u_plus_dirac(p, m)
g, g1, g2 = two_packets(p, np.cos(np.pi / 4), np.sin(np.pi / 4), 0.0, 0.5, 1.0, 0.3)
Os = overlap_spec(spinor_state(g, u), p, tvis, m, 0.05, kind="dirac")
O1 = overlap_spec(spinor_state(g1, u), p, tvis, m, 0.05, kind="dirac")
O2 = overlap_spec(spinor_state(g2, u), p, tvis, m, 0.05, kind="dirac")
E["E4"] = dict(descr="visibilita' |O(t)| dell'orologio (delta = 0.05): sovrapposizione contro miscuglio",
               config=dict(p1=0.5, p2=1.0, Dl=0.3, alpha="pi/4", phi=0, delta=0.05, t=tvis),
               previsione=dict(V_sup=[float(abs(x)) for x in Os], V_mix=[float(abs(0.5 * (x + y))) for x, y in zip(O1, O2)]),
               soglie={"max_t |V_regola(1/32) - V_cont|": 1e-3})

# ---------------------------------------------------------------- E5 doppione
pred5 = {}
for qb in [0.0, 0.5, 1.0]:
    g, _, _ = two_packets(p, 1.0, 0.0, 0.0, qb, qb, 0.05)
    pred5[str(qb)] = -rate_cont(g, p)
E["E5"] = dict(descr="pacchetto attorno al doppione k = pi/a - q: l'orologio gira all'indietro, R = -<m/E_q>",
               config=dict(Dl=0.05, q=[0.0, 0.5, 1.0]), previsione=pred5,
               soglie={"|R_regola(1/32) - prev| / |prev|": 2e-3})

# ---------------------------------------------------------------- E6 Zitterbewegung dell'orologio
tz = np.linspace(0, 40, 161)
sx = 1.0
pz = pgrid()
gz = np.exp(-(pz - 1.0) ** 2 * sx ** 2)       # |psi(x)|^2 gaussiana con sigma_x = sx, impulso medio 1
psi0 = np.zeros((len(pz), 2), complex)
psi0[:, 0] = gz / np.sqrt(np.sum(gz ** 2))
Oz = overlap_spec(psi0, pz, tz, m, 1e-5, kind="dirac")
tauz = -np.angle(Oz) / 1e-5
wz = gz ** 2 / np.sum(gz ** 2)
rate_mean = float(np.sum(wz * pz * m / (pz ** 2 + m ** 2)))
E["E6"] = dict(descr=("stato di sola chiralita' R (energie + e - insieme), sigma_x = 1, impulso medio 1: l'orologio "
                      "di massa ha la Zitterbewegung (2E) e il ramo negativo gira all'indietro; ritmo medio <p m/E^2>"),
               config=dict(sigma_x=sx, pbar=1.0, t=list(tz)),
               previsione=dict(tau_cont=list(map(float, tauz)), ritmo_medio=rate_mean),
               soglie={"max_t |tau_regola(1/32) - tau_cont|": 2e-3,
                       "|tau_regola(40)/40 - ritmo_medio| / ritmo_medio": 0.02})

E["E7"] = dict(descr="numeri per uno ione intrappolato (Gerritsma 2010): nessun esame, solo previsione nel documento")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    json.dump(pre, f, indent=1, default=float)
print("scritto", OUT)
for r in pred2a[:9]:
    print(r)
for r in pred2b:
    print({k: r[k] for k in ("s", "alpha", "phi", "dq_cont", "dq_SA_NR")})
for r in pred2c:
    print(r)
print(pred3)
print(E["E4"]["previsione"])
print(pred5)
print("E1", {k: (v["rate_cont"], v["coeff_a2_previsto"]) for k, v in pred1.items()})
print("E6 tau", tauz[::20], rate_mean, tauz[-1]/40)
print("E2a max rel |cont-SA|", max(abs(r["dq_cont"]-r["dq_SA"])/abs(r["dq_SA"]) for r in pred2a if abs(r["dq_SA"])>1e-9))

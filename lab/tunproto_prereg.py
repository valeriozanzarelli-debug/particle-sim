"""tunproto_prereg.py - PRE-REGISTRAZIONE del protocollo per l'orologio di massa nel tunnel.
Scrive lab/out/tunproto_prereg.json: configurazioni, previsioni (dallo stazionario esatto della regola e di Dirac
mediato sullo spettro del pacchetto), numero di ripetizioni per 5 sigma, esami con soglie PASS/FAIL.
Le CORSE (tunproto_rA/rC/rD) fanno gli esperimenti simulati tick per tick (regola), split-step (Dirac) e modi
accoppiati (guide d'onda), con rumore e conteggi finiti, e si giudicano contro questo file."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunproto_core import (OUT, spectral_rule, stationary_rule, stationary_dirac, dirac_spectral)

# ---------------------------------------------------------------- statistica (rumore di proiezione dell'ancella)
def c_tilt(sz, delta):
    """Var(tau_z^est) * N_z  (N_z = trasmessi misurati nella base z)"""
    return 1.0 / (delta ** 2 * (1 - sz ** 2))


def c_plane(sx, sy, delta):
    """Var(tau_y^est) * N_xy  (N_xy trasmessi, meta' in base x e meta' in base y)"""
    r2 = sx ** 2 + sy ** 2
    vsx = 2 * (1 - sx ** 2)
    vsy = 2 * (1 - sy ** 2)
    return (sy ** 2 * vsx + sx ** 2 * vsy) / r2 ** 2 / delta ** 2


def slope_c(Ls, cs):
    """Var(pendenza) * N  per un fit lineare pesato con varianze cs/N"""
    w = 1 / np.asarray(cs)
    L = np.asarray(Ls)
    Lm = np.sum(w * L) / np.sum(w)
    return 1 / np.sum(w * (L - Lm) ** 2)


def lattice_thick(theta, eps_a):
    """formule esatte della regola (barriera spessa): pendenze di tau_mz e tau_z, rapporto tau_my/tau_y"""
    ka = np.arccosh(np.cos(eps_a) / np.cos(theta))
    sm = np.cos(eps_a) * np.sin(theta) / (np.cos(theta) ** 2 * np.sinh(ka))
    sv = np.sin(eps_a) / (np.cos(theta) * np.sinh(ka))
    return dict(slope_m=float(sm), slope_V=float(sv), ratio=float(np.tan(eps_a) / np.tan(theta)))


def main():
    P = dict(studio="PROTOCOLLO SPERIMENTALE per l'orologio di massa nel tunnel (previsione del giro 1)",
             prefisso="tunproto_", data="2026-10-13", autore="Claude (Opus 5.5) per Valerio",
             regola="lab/zz_core.py dirac1 via lab/tunnel_core.step (th e V per sito)",
             parametri_liberi=0,
             misura=("ancella con stati 0/1: massa (o potenziale) dentro la barriera +-delta/2; ancella preparata in "
                     "(|0>+|1>)/sqrt2; si legge il vettore di Bloch dei soli trasmessi. tau_y^est = -arg(s_x+i s_y)/delta, "
                     "tau_z^est = -atanh(s_z)/delta. Previsione stazionaria: s_z = -tanh(delta tau_z), "
                     "s_x + i s_y = e^{-i delta tau_y}/cosh(delta tau_z)."),
             statistica=("N = numero di particelle TRASMESSE e rivelate (con l'ancella misurata). Inclinazione: tutte in base z, "
                         "Var = 1/(N delta^2 (1 - s_z^2)). Nel piano: meta' in base x, meta' in y, propagazione degli errori. "
                         "N_5sigma = 25 * (c1 + c2) / Delta^2 con c = Var*N delle due grandezze confrontate."))

    # ================================================================ A: la misura ridà tau_m (generico, Dirac)
    A = dict(E0=2.0, sigma=10.0, L=4.0, delta=0.05, V0=[1.6, 2.0, 2.4], passi=[0.05, 0.025, 0.0125],
             dirac_split=dict(dx=0.05, dt=0.0125))
    A["previsioni_dirac_spettrale"] = {}
    A["stazionario_dirac_E0"] = {}
    for V0 in A["V0"]:
        A["previsioni_dirac_spettrale"][str(V0)] = {c: dirac_spectral(2.0, V0, 4.0, 10.0, 0.05, c) for c in ("m", "V")}
        A["stazionario_dirac_E0"][str(V0)] = stationary_dirac(2.0, V0, 4.0)
    A["limite_debole_stazionario"] = {}
    for V0 in A["V0"]:
        A["limite_debole_stazionario"][str(V0)] = {
            str(s): {c: dirac_spectral(2.0, V0, 4.0, s, 1e-3, c, nk=801) for c in ("m", "V")} for s in (10, 20, 40, 80)}
    P["A_misura"] = A
    st = A["stazionario_dirac_E0"]
    P["esami_A"] = {
        "A1": "regola tick per tick (a = 0,05) = previsione spettrale della regola stessa: |d tau| <= 1e-5 (tz, ty, due orologi, 3 V0)",
        "A2": "regola tick per tick estrapolata a a -> 0 (lineare nei due passi piu' piccoli) = Dirac spettrale: |d tau| <= 1e-3",
        "A3": "Dirac split-step (dx 0,05, dt 0,0125) = Dirac spettrale: |d tau| <= 5e-3 (tz e ty, orologio di massa, 3 V0); |dP_T|/P_T <= 1e-3",
        "A4": "limite debole e stazionario: Dirac spettrale a delta = 1e-3, sigma = 80 -> stazionario tau(E0) entro 1e-2; scarto decrescente con sigma (10 -> 80)",
        "A5": ("lettura (regola tick per tick, a = 0,0125): V0 = 2,4: tau_my^est < 0 e tau_y^est > 0; V0 = 2: tau_mz^est/L in [0,97; 1,03] "
               "e tau_z^est/L < 0,15; V0 = 1,6 e 2,4: tau_my^est/tau_y^est = (eps/m) entro 0,05 (a L = 4 stazionario: %.4f e %.4f)"
               % (st["1.6"]["tau_my"] / st["1.6"]["tau_y"], st["2.4"]["tau_my"] / st["2.4"]["tau_y"])),
    }

    # ================================================================ B: effetti sistematici generici
    B = {}
    a_b = 0.025
    B["bordi_lisci"] = dict(a=a_b, E0=2.0, V0=2.0, L=[4.0, 6.0], w=[0.0, 0.25, 0.5, 1.0], sigma=10.0, delta=0.05, val={})
    for w in B["bordi_lisci"]["w"]:
        for prof in (False, True):
            key = "w%.2f_%s" % (w, "profilo" if prof else "rettangolo")
            s4 = stationary_rule(a_b, 2.0, 2.0, 4.0, w, clock_profile=prof)
            s6 = stationary_rule(a_b, 2.0, 2.0, 6.0, w, clock_profile=prof)
            B["bordi_lisci"]["val"][key] = dict(L4=s4, L6=s6, slope_m=(s6["tau_mz"] - s4["tau_mz"]) / 2,
                                                slope_V=(s6["tau_z"] - s4["tau_z"]) / 2,
                                                spettrale_L4={c: spectral_rule(a_b, 2.0, 2.0, 4.0, 10.0, 0.05, c, w=w,
                                                                               clock_profile=prof, nk=201) for c in ("m", "V")})
    B["delta_finito"] = dict(E0=2.0, V0=2.0, L=4.0, sigma=10.0,
                            val={str(d): dirac_spectral(2.0, 2.0, 4.0, 10.0, d, "m") for d in (1e-3, 0.01, 0.05, 0.1, 0.2, 0.4, 0.8)})
    B["delta_finito_V2.4"] = {str(d): dirac_spectral(2.0, 2.4, 4.0, 10.0, d, "m") for d in (1e-3, 0.05, 0.2, 0.4, 0.8)}
    B["decoerenza_ancella"] = ("defasamento puro dell'ancella: (s_x, s_y) -> C (s_x, s_y), s_z invariato: tau_y^est non distorto, "
                               "errore x 1/C; tau_z^est invariato. Rilassamento T1 dell'ancella: distorce s_z (va calibrato).")
    P["B_sistematici"] = B
    P["esami_B"] = {
        "B1": "bordi lisci: regola tick per tick (a = 0,025, L = 4) = previsione spettrale con lo stesso profilo: |d tau| <= 1e-4",
        "B1r": "robustezza (dichiarata, calcolata dallo stazionario qui sopra): |pendenza tau_mz - 1| <= 0,05 per w <= 0,5 con orologio a rettangolo",
        "B2": "delta finito: lo scarto relativo di tau_mz^est rispetto a delta = 1e-3 e' <= 1 % finche' delta tau_mz <= 0,8; la regola tick per tick lo riproduce entro 1e-4",
    }

    # ================================================================ C: cammino in fibra = la regola a theta finito
    th = 0.5
    C = dict(piattaforma="cammino quantistico fotonico multiplato nel tempo (moneta programmabile per posto e passo)",
             theta=th, a=th, nota_unita="a = theta: 1 posto = 1 passo = a; tempi in passi = tau/a",
             E0=2.0, omega_a=2.0 * th, V0=[1.6, 2.0, 2.4], N_siti=[2, 4, 6, 8], sigma=2.0, sigma_siti=4,
             taglio="4 sigma (stato preparato solo a sinistra della barriera)", x0=-8.0, passi=40, margin=1.0,
             delta=0.15, delta_moneta_rad=0.15 * th / 2, delta_fase_EOM_rad=0.15 * th / 2,
             formule_regola_spessa={str(V0): lattice_thick(th, (2.0 - V0) * th) for V0 in (1.6, 2.0, 2.4)},
             val={})
    for V0 in C["V0"]:
        for N in C["N_siti"]:
            L = N * th
            C["val"]["V%.1f_N%d" % (V0, N)] = dict(
                staz=stationary_rule(th, 2.0, V0, L),
                spett={c: spectral_rule(th, 2.0, V0, L, 2.0, 0.15, c) for c in ("m", "V")})
    # libero (calibrazione di delta): tempo proprio della regola sulla zona dell'orologio
    C["calibrazione_libera"] = {str(N): spectral_rule(th, 2.0, 0.0, N * th, 2.0, 0.15, "m") for N in (4, 8)}
    # statistica D1: pendenza orologio di massa contro pendenza orologio di potenziale (centro del gap, N = 4,6,8)
    Ls = [N * th for N in (4, 6, 8)]
    vm = [C["val"]["V2.0_N%d" % N]["spett"]["m"] for N in (4, 6, 8)]
    vv = [C["val"]["V2.0_N%d" % N]["spett"]["V"] for N in (4, 6, 8)]
    sl_m = np.polyfit(Ls, [v["tz"] for v in vm], 1)[0]
    sl_v = np.polyfit(Ls, [v["tz"] for v in vv], 1)[0]
    cm = slope_c(Ls, [c_tilt(v["sz"], 0.15) for v in vm])
    cv = slope_c(Ls, [c_tilt(v["sz"], 0.15) for v in vv])
    N5_D1 = 25 * (cm + cv) / (sl_m - sl_v) ** 2
    # D2: segno di tau_my a V0 = 2,4, N = 8
    v = C["val"]["V2.4_N8"]["spett"]["m"]
    N5_D2 = 25 * c_plane(v["sx"], v["sy"], 0.15) / v["ty"] ** 2
    # D1 bis: singolo punto (N = 8): tau_mz contro tau_z
    a8, b8 = C["val"]["V2.0_N8"]["spett"]["m"], C["val"]["V2.0_N8"]["spett"]["V"]
    N5_D1b = 25 * (c_tilt(a8["sz"], 0.15) + c_tilt(b8["sz"], 0.15)) / (a8["tz"] - b8["tz"]) ** 2
    C["discriminazioni"] = dict(
        D1=dict(cosa="pendenza d tau/dL al centro del gap: orologio di massa contro orologio di potenziale (Buttiker-Landauer)",
                pendenza_m=float(sl_m), pendenza_V=float(sl_v), pendenza_m_spessa_regola=1 / np.cos(th),
                N_5sigma_per_punto_per_orologio=float(N5_D1)),
        D1b=dict(cosa="un solo spessore (N = 8 siti): tau_mz contro tau_z", tau_mz=a8["tz"], tau_z=b8["tz"],
                 N_5sigma_per_orologio=float(N5_D1b)),
        D2=dict(cosa="segno di tau_my a V0 = 2,4 (V0 > E), N = 8 siti", tau_my=v["ty"],
                N_5sigma=float(N5_D2), P_T=v["PT"]))
    # delta piu' grande (il bias da delta finito e' piccolo, vedi B2): stesse discriminazioni a delta = 0,6
    d2 = 0.6
    vm2 = [spectral_rule(th, 2.0, 2.0, N * th, 2.0, d2, "m") for N in (4, 6, 8)]
    vv2 = [spectral_rule(th, 2.0, 2.0, N * th, 2.0, d2, "V") for N in (4, 6, 8)]
    sl_m2 = np.polyfit(Ls, [v_["tz"] for v_ in vm2], 1)[0]
    sl_v2 = np.polyfit(Ls, [v_["tz"] for v_ in vv2], 1)[0]
    N5_D1_2 = 25 * (slope_c(Ls, [c_tilt(v_["sz"], d2) for v_ in vm2]) + slope_c(Ls, [c_tilt(v_["sz"], d2) for v_ in vv2])) / (sl_m2 - sl_v2) ** 2
    w2 = spectral_rule(th, 2.0, 2.4, 8 * th, 2.0, d2, "m")
    N5_D2_2 = 25 * c_plane(w2["sx"], w2["sy"], d2) / w2["ty"] ** 2
    C["delta_0.6"] = dict(spett_centro_m=vm2, spett_centro_V=vv2, pendenza_m=float(sl_m2), pendenza_V=float(sl_v2),
                          N_5sigma_D1=float(N5_D1_2), spett_V24_N8=w2, N_5sigma_D2=float(N5_D2_2))
    C["rumore"] = dict(fase_comune=0.05, fase_indipendente=0.05, moneta_comune=0.01, pdl=0.03, realizzazioni=200,
                       nota="rad per sito per passo; pdl = frazione di potenza persa dalla componente R per passo")
    P["C_fibra"] = C
    P["esami_C"] = {
        "C1": "tick per tick (taglio 4 sigma, 40 passi) = previsione spettrale: |d tau| <= max(0,03 |tau|, 0,05), 12 configurazioni x 2 orologi",
        "C2": "pendenza tau_mz^est (N = 4,6,8) al centro entro 10 %% di 1/cos(theta) = %.4f; pendenza tau_z^est < 0,3" % (1 / np.cos(th)),
        "C3": "Monte Carlo (2000 esperimenti, conteggi binomiali) a N = N_5sigma(D1) = %.0f per punto per orologio: significativita' mediana >= 4,5" % N5_D1,
        "C4": "Monte Carlo a N = N_5sigma(D2) = %.0f: significativita' mediana di tau_my < 0 >= 4,5" % N5_D2,
        "C3b_C4b": "come C3 e C4 ma con delta = 0,6 e i suoi N_5sigma (D1: %.0f, D2: %.0f)" % (N5_D1_2, N5_D2_2),
        "C5": "rumore (fase comune 0,05; fase indipendente 0,05; moneta 0,01; pdl 3 %): pendenza di tau_mz entro 10 % del caso pulito, segno di tau_my (V0 = 2,4, N = 8) conservato, contrasto r >= 0,5. Ogni caso PASS/FAIL separato",
    }

    # ================================================================ D: guide d'onda binarie
    D = dict(piattaforma="reticolo binario di guide d'onda (silice, laser fs), orologio = polarizzazione",
             kappa_cm=4.5, sigma_massa_cm=1.5, c_guide_per_cm=9.0, compton_guide=6.0, Z_cm=12.0, E0_cm=3.0,
             V0_cm=[2.4, 3.0, 3.6], L_guide=[12, 18, 24], sx_guide=9, delta_cm=0.3, margin=6,
             birifrangenza_alternata_dn=float(0.3 / 2 * 633e-9 / (2 * np.pi) * 100), eta_kappa=0.01, disordine_cm=0.05,
             realizzazioni_disordine=50,
             previsioni_dirac={})
    for V0 in D["V0_cm"]:
        for Lg in D["L_guide"]:
            D["previsioni_dirac"]["V%.1f_L%d" % (V0, Lg)] = {c: dirac_spectral(2.0, V0 / 1.5, Lg / 6.0, 1.5, 0.2, c) for c in ("m", "V")}
    P["D_guide"] = D
    P["esami_D"] = {
        "D1": "modi accoppiati (apparato) contro Dirac spettrale (stesse unita' Compton): |d tau_mz| <= 10 % a tutti i 9 punti; pendenza al centro (L = 2..4 Compton) entro 15 % di 1",
        "D2": "contaminazione: accoppiamento dipendente dalla polarizzazione eta = 1 %: spostamento di tau_mz^est (centro, L = 3) <= 10 %; dopo sottrazione della corsa a orologio spento <= 1 %",
        "D3": "disordine di fabbricazione 0,05 cm^-1 (50 campioni): media di tau_mz^est entro 5 % del pulito, dispersione tra campioni <= 10 %",
    }
    P["debiti_dichiarati"] = [
        "barriera rettangolare e orologi messi a mano (apparato di misura)",
        "parametri delle piattaforme (theta = 0,5, 40 passi; kappa = 4,5 cm^-1, sigma = 1,5 cm^-1, 12 cm) scelti da me: plausibili, non verificati come fattibili",
        "modello delle guide d'onda (modi accoppiati a primi vicini, Longhi 2010) = modello dell'apparato, non la regola",
        "rumore modellato come gaussiano bianco per sito e per passo",
    ]
    with open(os.path.join(OUT, "tunproto_prereg.json"), "w") as f:
        json.dump(P, f, indent=1, default=float)
    print("delta 0.6: D1", N5_D1_2, "D2", N5_D2_2, sl_m2, sl_v2, w2["ty"]); print("D1 N5 per punto:", N5_D1, " D1b:", N5_D1b, " D2:", N5_D2, " slopes", sl_m, sl_v)


if __name__ == "__main__":
    main()

"""lv_r2.py - esame E9: massimo passo a del reticolo permesso da ogni esperimento (con la dispersione esatta
della regola) + controlli aggiunti (EM2). Scrive lab/out/lv_r2.json. Pochi secondi."""
import json
import os
import numpy as np
import lv_core as lv

HERE = os.path.dirname(os.path.abspath(__file__))
HC = lv.HBARC_GEVM
m = lv.ME_GEV
LP = lv.L_P
out = {"A_casi": {"peggiore (diagonale)": 1 / 3, "media sulle direzioni": 1 / 5}}


def rec(a):
    return {"a_max_m": a, "a_max_lP": a / LP, "E_a_min_GeV": HC / a}


B = {}
for lab, A in out["A_casi"].items():
    b = {}
    # 1. tempo di volo, n=2: v = 1 - (A/2)(E/E_a)^2 = 1 - (3/2)(E/E_QG2)^2  ->  E_a = E_QG2 sqrt(A/3)
    for name, Eqg in (("LHAASO GRB 221009A (PRL 133, 071501, 2024): E_QG2 > 6e-8 E_Pl", 6e-8 * lv.E_PL),
                      ("H.E.S.S. Mrk 501 (ApJ 870, 93, 2019): E_QG2 > 8.5e10 GeV", 8.5e10),
                      ("MAGIC GRB 190114C (PRL 125, 021301, 2020): E_QG2 > 6.3e10 GeV", 6.3e10)):
        b["tempo di volo: " + name] = rec(HC / (Eqg * np.sqrt(A / 3)))
    # 2. decadimento del fotone gamma -> e+e- (soglia k_th = 2 sqrt(m E_a) A^-1/4 > E_osservata)
    for name, E in (("LHAASO J2032+4102, 1.42 PeV (limite 95%: 1.21 PeV; Nature 594, 33, 2021)", 1.21e6),
                    ("Granchio, 1.1 PeV (LHAASO, Science 373, 425, 2021)", 1.1e6)):
        Ea = np.sqrt(A) * E ** 2 / (4 * m)
        b["decadimento del fotone: " + name] = rec(HC / Ea)
    # 3. cavita': anisotropia della velocita' di fase fra due bracci ortogonali, Delta A <= 1/4
    dA = 0.25
    for name, E, d in (("risonatori ottici 1064 nm (Herrmann et al., PRD 80, 105011, 2009): dc/c < 1e-17", 1.165e-9, 1e-17),
                       ("risonatori a microonde 13 GHz (Nagel et al., Nat. Commun. 6, 8174, 2015): dc/c < 1e-18", 5.4e-14, 1e-18)):
        Ea = E * np.sqrt(dA / (6 * d))
        b["cavita': " + name] = rec(HC / Ea)
    # 4. indicativo: Lorentz massimo dell'elettrone (gamma_max = sqrt(E_a / (4 sqrt(A) m))) >= 4.5e9 (2.3 PeV)
    g = 2.3e6 / m
    Ea = 4 * np.sqrt(A) * m * g ** 2
    b["INDICATIVO elettroni del Granchio 2.3 PeV (Lorentz massimo, stile JLM)"] = rec(HC / Ea)
    B[lab] = b
out["limiti"] = B

# variante a un tick (NON la regola): termine n=1  omega = k (1 -+ nxnynz k a /2)
s = 1 / (3 * np.sqrt(3))   # nxnynz massimo (diagonale)
out["variante un tick (diagnosi)"] = {
    "tempo di volo n=1, E_QG1 > 10 E_Pl (LHAASO 2024): a_max": rec(HC / (10 * lv.E_PL) * s),
    "decadimento del fotone n=1 superluminale, E_LV1 >= 1.7e24 GeV (Chen et al. 2021): a_max": rec(HC / (1.7e24 * s)),
    "nota": "nella direzione opposta il fotone e' superluminale: la variante senza palindromo e' esclusa per a > ~1e-5 l_P"}

# ritardo previsto per a = l_P (GRB 221009A, z = 0.151, E = 10 TeV)
z = np.linspace(0, 0.151, 2001)
Om, H0 = 0.315, 67.4e3 / 3.0857e22
K2 = np.trapezoid((1 + z) ** 2 / np.sqrt(Om * (1 + z) ** 3 + 1 - Om), z)
EQG2 = lv.E_PL * np.sqrt(3 / (1 / 3))
dt = 1.5 * (1e4 / EQG2) ** 2 * K2 / H0
out["ritardo previsto a = l_P, GRB 221009A, 10 TeV, diagonale"] = {"s": dt, "E_QG2_GeV": EQG2}
out["soglia decadimento fotone a = l_P"] = {"diagonale_GeV": lv.kth_photon_decay_analytic(1 / 3, lv.E_PL),
                                           "media_GeV": lv.kth_photon_decay_analytic(1 / 5, lv.E_PL)}

# laboratorio
labs = {}
for th in (0.05, 0.1, 0.5, 1.0):
    a = th * lv.LAMBDA_C
    Ea = HC / a
    labs["a = %.2f lambda_C" % th] = {"a_m": a, "a/l_P": a / LP,
                                      "soglia decadimento fotone (diag) MeV": lv.kth_photon_decay_analytic(1 / 3, Ea) * 1e3,
                                      "anisotropia della luce visibile (550 nm)": 0.25 / 6 * (2.25e-9 / Ea) ** 2}
out["laboratorio"] = labs

# controlli aggiunti EM2: stessa specie o stesso coefficiente -> vietato
chk = {}
Ea = 1e4
for nd in ([1, 1, 1], [1, 2, 3]):
    nd = np.array(nd, float) / np.linalg.norm(nd)
    worst3 = -np.inf
    for p in np.linspace(50, 2000, 40):
        Ein = lv.omega_dirac_closed((p * nd / Ea)[None], 1e-9)[0] * Ea
        for x1 in np.linspace(0.02, 0.9, 30):
            x2 = np.linspace(0.01, 1 - x1 - 0.01, 30)
            x3 = 1 - x1 - x2
            Ef = sum(lv.omega_dirac_closed((xx * p)[:, None] * nd[None] / Ea, 1e-9) * Ea for xx in (np.full_like(x2, x1), x2, x3))
            worst3 = max(worst3, np.max(Ein - Ef))
    chk["nu -> 3 nu / e -> 3e (m~0), %s" % np.round(nd, 3).tolist()] = worst3
out["controlli EM2 (max Delta E, deve essere < 0)"] = chk

best = {}
for lab_, b in B.items():
    robust = {k: v for k, v in b.items() if "INDICATIVO" not in k}
    kmin = min(robust, key=lambda k: robust[k]["a_max_m"])
    best[lab_] = {"piu' severo (robusto)": kmin, **robust[kmin]}
out["piu_severo"] = best
crit = best["media sulle direzioni"]["a_max_lP"]
out["verdetto"] = ("compatibile" if crit >= 10 else "in tensione" if crit >= 1 else "esclusa per a >= l_P")


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    return o


with open(os.path.join(HERE, "out", "lv_r2.json"), "w") as f:
    json.dump(clean(out), f, indent=1, ensure_ascii=False)
for lab_, b in B.items():
    print("==", lab_)
    for k, v in b.items():
        print("  %-110s a_max = %.2e m = %.2e l_P" % (k[:110], v["a_max_m"], v["a_max_lP"]))
print(json.dumps(clean({k: out[k] for k in ("variante un tick (diagnosi)", "ritardo previsto a = l_P, GRB 221009A, 10 TeV, diagonale",
                                            "soglia decadimento fotone a = l_P", "laboratorio",
                                            "controlli EM2 (max Delta E, deve essere < 0)", "piu_severo", "verdetto")}), indent=1))

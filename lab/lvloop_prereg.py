"""lvloop_prereg.py - pre-registration of the study lvloop_ (Lorentz violation carried by loops).
Written BEFORE running lvloop_sym.py and BEFORE reading any phase-3 number of the vuoto_ study.
Output: lab/out/lvloop_prereg.json.  Free parameters: 0 (the lattice step a is the unknown).

Contents
 1. symmetry exams S1-S4 on the interacting rule (lvloop_sym.py)
 2. SME coefficients allowed in each symmetry case
 3. bounds used (only as exams), with sources
 4. conditional table: scaling scenario -> a_max -> verdict
 5. the numbers we need from vuoto_ phase 3 (N1-N5) and how the table is applied
"""
import json, os
import numpy as np

ALPHA = 7.2973525693e-3
LAMBDA_C = 3.8615926796e-13      # m, hbar/(m_e c)
L_PL = 1.616255e-35              # m
ME_EV = 0.51099895e6
HBARC_EVM = 1.973269804e-7       # eV m

# ---------------------------------------------------------------- bounds (dimensionless, on |delta v| / c)
# sign convention: delta_iso = (c_e - c_gamma)/c  (electron limiting speed minus photon speed)
BOUNDS = {
  "iso_photon_faster": [   # delta_iso < 0  -> photon decay gamma -> e+ e-
    {"name": "decadimento del fotone negli sciami UHE (Duenkel, Niechciol, Risse, PRD 104, 015010 (2021))",
     "B": 6e-21, "E_eV": 2**0.5 * ME_EV / 6e-21**0.5, "verified": True, "note": "kappa > -6e-21 (98% CL, SIBYLL 2.3d); limite su fotone contro e+e-"},
    {"name": "fotone 1,42 PeV non decade (LHAASO, Nature 594, 33 (2021)); conto nostro 2 m^2/E^2",
     "B": 2 * (ME_EV / 1.42e15)**2, "E_eV": 1.42e15, "verified": True, "note": "soglia k_th = sqrt(2) m / sqrt(delta)"},
    {"name": "Tevatron, fotoni non decadono (Hohensee et al. PRL 102, 170402 (2009))",
     "B": 5.8e-12, "E_eV": 3.4e11, "verified": True, "note": "citato in Altschul PRD 80, 091901 (2009)"},
  ],
  "iso_electron_faster": [  # delta_iso > 0 -> vacuum Cherenkov e -> e gamma
    {"name": "elettroni PeV del Granchio (Li e Ma, PLB 829, 137034 (2022)); conto nostro m^2/(2E^2) a 2,3 PeV",
     "B": (ME_EV / 2.3e15)**2 / 2, "E_eV": 2.3e15, "verified": True, "note": "indicativo: energia dedotta dal modello IC"},
    {"name": "elettroni cosmici ~40 TeV (H.E.S.S. 2024); conto nostro m^2/(2E^2)",
     "B": (ME_EV / 4.0e13)**2 / 2, "E_eV": 4.0e13, "verified": False, "note": "energia massima a memoria"},
    {"name": "LEP, niente Cherenkov a 104,5 GeV (Hohensee et al. 2009)",
     "B": 1.2e-11, "E_eV": 1.045e11, "verified": True, "note": ""},
  ],
  "iso_two_sided": [
    {"name": "LEP, perdite di sincrotrone (Altschul, PRD 80, 091901 (2009))",
     "B": 5e-15, "E_eV": 1.045e11, "verified": True, "note": "|kappa_tr - 4/3 c_00| < 5e-15"},
  ],
  "aniso_electron": [   # rank-2 traceless electron c_JK
    {"name": "Yb+ E3/E2, confronto di orologi (Sanner et al., Nature 567, 204 (2019))",
     "B": 1.6e-20, "verified": True, "note": "c_X-Y 0,9 +- 1,6e-20; c_XY -6,9 +- 8,0e-21; c_XZ 1,3 +- 1,3e-20"},
    {"name": "Yb+ Ramsey a impulsi composti (Dreissen et al., Nat. Commun. 13, 7314 (2022))",
     "B": 1e-21, "verified": True, "note": "'livello 1e-21' dall'abstract; componenti non lette una per una"},
    {"name": "Ca+ (Pruttivarasin et al., Nature 517, 592 (2015))", "B": 2.3e-18, "verified": True, "note": "dalla tabella di Sanner"},
    {"name": "Dy (Hohensee et al., PRL 111, 050401 (2013))", "B": 5.7e-17, "verified": True, "note": "dalla tabella di Sanner"},
  ],
  "aniso_photon_nonbiref": [
    {"name": "risonatori (Nagel et al., Nat. Commun. 6, 8174 (2015)); kappa_e-", "B": 1e-18, "verified": False,
     "note": "ordine di grandezza a memoria"},
  ],
  "photon_biref": [
    {"name": "polarizzazione di sorgenti cosmologiche/GRB (Kostelecky e Mewes, PRL 87, 251304 (2001); PRL 110, 201601 (2013))",
     "B": 1e-32, "verified": False, "note": "d=4 birifrangente <~1e-32 (sorgenti radio/CMB), GRB piu' forte; a memoria"},
  ],
}


def a_max_mass_type(B, p, kappa=1.0):
    """delta = alpha kappa theta^p with theta = m a / hbar (constant at all energies). returns a_max [m] or None"""
    if p == 0:
        return None if ALPHA * kappa > B else np.inf
    th = (B / (ALPHA * kappa))**(1.0 / p)
    return th * LAMBDA_C


def a_max_momentum_type(B, p, E_eV, kappa=1.0):
    """delta(E) = alpha kappa (E a / hbar c)^p, bound B at energy E_eV (threshold energy of the test)"""
    if p == 0:
        return a_max_mass_type(B, 0, kappa)
    return (B / (ALPHA * kappa))**(1.0 / p) * HBARC_EVM / E_eV


def verdict(amax):
    if amax is None:
        return "ESCLUSA per qualunque a (salvo kappa < B/alpha: fine tuning)"
    if amax == np.inf:
        return "salva"
    if amax >= 10 * L_PL:
        return "salva con passo di Planck (a_max = %.1e m = %.1e l_P)" % (amax, amax / L_PL)
    return "salva solo con a < %.1e m = %.1e l_P (sotto Planck: di fatto esclusa)" % (amax, amax / L_PL)


def table():
    rows = []
    tests = [("iso", BOUNDS["iso_photon_faster"][0], "gap isotropo, fotone piu' veloce (decadimento)"),
             ("iso", BOUNDS["iso_photon_faster"][1], "gap isotropo, fotone piu' veloce (LHAASO PeV)"),
             ("iso", BOUNDS["iso_electron_faster"][0], "gap isotropo, elettrone piu' veloce (Cherenkov, Granchio)"),
             ("iso", BOUNDS["iso_two_sided"][0], "gap isotropo, entrambi i segni (LEP sincrotrone)"),
             ("an", BOUNDS["aniso_electron"][1], "anisotropia elettrone rango 2 (Yb+ 2022)"),
             ("an", BOUNDS["aniso_electron"][0], "anisotropia elettrone rango 2 (Yb+ 2019)"),
             ("an", BOUNDS["aniso_photon_nonbiref"][0], "anisotropia luce non birifrangente (risonatori)"),
             ("an", BOUNDS["photon_biref"][0], "birifrangenza della luce (se eps e mu anisotropi diversi)")]
    for kind, b, label in tests:
        for p in (0, 1, 2):
            am = a_max_mass_type(b["B"], p)
            row = {"test": label, "bound": b["B"], "exponent": p, "a_max_mass_type_m": None if am is None else float(am),
                   "verdict_mass_type": verdict(am)}
            if b.get("E_eV"):
                amp = a_max_momentum_type(b["B"], p, b["E_eV"])
                row["a_max_momentum_type_m"] = None if amp is None else float(amp)
                row["verdict_momentum_type"] = verdict(amp)
            rows.append(row)
    return rows


def main():
    pre = {
      "study": "lvloop_ - violazione di Lorentz portata dai loop (Collins-Perez-Sudarsky-Urrutia-Vucetich 2004)",
      "date": "2026-10-13",
      "free_parameters": 0,
      "written_before": ["lvloop_sym.py run", "any vuoto_ phase-3 number"],
      "symmetry_claim": {
        "statement": ("La dispersione LIBERA della regola e' invariante per tutte le 48 operazioni del cubo (lv_ E6), ma "
                      "la DINAMICA (sequenza dei passi x y z C z y x C) e' invariante solo per: inversione P, rotazioni di "
                      "pi attorno agli assi, e lo scambio x<->z combinato con una traslazione di un tick. Gruppo previsto: "
                      "D_4h (16 elementi) con asse distinto y (il passo di mezzo). Le rotazioni di ordine 3 (x->y->z) NON "
                      "sono simmetrie. Con l'ordine a giro (6 permutazioni, variante 'six') la rotazione di ordine 3 diventa "
                      "una traslazione di 2 tick e x<->z diventa l'inversione temporale: gruppo previsto O_h (48)."),
        "consequence": ("Gruppo cubico (anche solo T): ogni tensore simmetrico a traccia nulla di rango 2 invariante e' zero, "
                        "ogni vettore e' zero. Tutti i coefficienti di dimensione 4 della SME minimale (c_JK, c_TJ, kappa_e-, "
                        "kappa_o+, kappa_e+, kappa_o-, b_J, ...) hanno al piu' rango 2 sotto SO(3): sopravvivono solo gli "
                        "scalari pari (c_TT dell'elettrone, kappa_tr della luce: cioe' il gap isotropo c_e - c_gamma). "
                        "D_4h con asse y: sopravvive in piu' UNA componente per ogni tensore di rango 2 pari per P e T: "
                        "c_yy - (c_xx+c_zz)/2 dell'elettrone, e la stessa per kappa_e- (luce non birifrangente) e per "
                        "kappa_e+ (luce BIRIFRANGENTE: permessa dalla simmetria, non forzata). Nessun vettore (P e C2 lo vietano).")
      },
      "exams_symmetry": {
        "S1": "regola congelata, 3-4 impurezze generiche, n = 6 e 8: numero di operazioni con spettro uguale (< 1e-9) = 16, "
              "permutazioni di assi nel gruppo = {identita', x<->z}. PASS se esattamente cosi'.",
        "S2": "variante 'six' (6 permutazioni a giro), n = 6: 48 operazioni. PASS se 48.",
        "S3": "coppia di impurezze n = 8: spettro(sep. x) = spettro(sep. z) entro 1e-9 e spettro(sep. y) diverso da essi "
              "di piu' di 1e-6; idem (1,1,0) = (0,1,1) != (1,0,1). PASS se tutte e quattro.",
        "S4": "variante 'fixed' (xyz a ogni tick): meno di 48 operazioni (solo descrittivo: nessuna soglia).",
        "check": "operatore costruito da lvloop_sym.py contro zz_core.Walk.tick su vettore casuale: < 1e-12"
      },
      "predictions_for_vuoto_phase3": {
        "P1": ("Media sulle due fette di tempo (o definizione indipendente dalla fetta: la quasi-energia del vestito lo e'): "
               "delta_v_x = delta_v_z ESATTAMENTE (simmetria x<->z con traslazione di un tick); la differenza x-z misurata "
               "(oggi -0,005 alpha) e' un errore della costruzione, non fisica."),
        "P2": ("Regola congelata: l'anisotropia delta_v_y - (delta_v_x+delta_v_z)/2 NON e' vietata e ci aspettiamo che resti "
               "FINITA per a -> 0 (q = 0, ordine alpha * O(10^-2)): e' il problema di Collins in forma tetragonale. "
               "NON e' la previsione 'deve sparire' della consegna: quella valeva solo se la dinamica fosse cubica."),
        "P3": ("Variante 'six': anisotropia di rango 2 a ordine a^0 vietata: deve andare a zero almeno come theta^2 "
               "(operatori di dimensione 6, anisotropia cubica n_x^4+n_y^4+n_z^4). PASS della cura se |aniso| cala di "
               ">= x3 per ogni dimezzamento di theta, o e' < 1e-3 alpha a tutti i theta."),
        "P4": ("Gap isotropo c_e - c_gamma: nessuna simmetria della regola lo vieta (tick e cella non sono legati da una "
               "simmetria del reticolo). Previsione: finito per a -> 0 (p = 0), ordine alpha*1e-2, salvo cancellazioni. "
               "ATTENZIONE: deve essere la differenza con la c della LUCE vestita (polarizzazione del vuoto), e il "
               "regolatore (proiettore K1,K2) da solo da' 0,053 alpha anche sul continuo: conta solo la parte che "
               "non dipende dal regolatore.")
      },
      "numbers_needed": {
        "N1": "Delta_iso(theta) = [delta_v_regola - delta_v_continuo, stesso proiettore]_media isotropa, mediato sulle fette, "
              "a >= 3 theta, con la variazione da (K1,K2) = (0,6;1,05) a (0,5;0,9). Fit alpha*kappa_iso*theta^p.",
        "N2": "Delta_an(theta) = delta_v_y - (delta_v_x+delta_v_z)/2 mediato sulle fette, a >= 3 theta. Fit alpha*kappa_an*theta^q.",
        "N3": "delta_v_x - delta_v_z mediato sulle fette (controllo: deve essere 0).",
        "N4": "Stessi N1, N2 per la variante con l'ordine a giro.",
        "N5": "Lato luce: anisotropia y contro (x,z) della risposta ELETTRICA (Pi(0), gia' misurata: 4e-5 dopo media) e di "
              "quella MAGNETICA (Pi trasverso a q -> 0). Se le due anisotropie sono diverse, la luce e' birifrangente (kappa_e+)."
      },
      "verdict_rules": {
        "R1": "a_max = lambda_C * (B/(alpha*kappa))^(1/p) per lo scenario 'tipo massa' (coefficiente costante ~ theta^p); "
              "'tipo impulso' (alpha*kappa*(E a)^p): a_max = (B/(alpha*kappa))^(1/p) * hbar c / E alla soglia del test.",
        "R2": "p = 0 (o q = 0) con kappa ~ 1e-2..1: ESCLUSA per qualunque a, a meno di sintonizzare un parametro nudo a "
              "1 parte su B/(alpha kappa) (manopola + fine tuning: dichiarato come tale, non e' una salvezza).",
        "R3": "a_max >= 10 l_P: 'salva con passo di Planck'. a_max < l_P: 'di fatto esclusa'.",
        "R4": "Il limite giro 1 (fotone PeV, a <~ 4,8e-31 m) resta valido e si combina col minimo."
      },
      "scenario_curable": ("Scenario C (chiesto dal coordinatore): anisotropia di rango 2 FINITA per a -> 0 nella regola "
        "congelata (q = 0) ma che SPARISCE con l'ordine a giro (N4 con q >= 2). Verdetto: la regola congelata e' esclusa "
        "come teoria fondamentale per qualunque a (Yb+: B/alpha ~ 1e-19); la variante 'six' prende il posto del giudizio, "
        "con il suo gap isotropo (N4) da confrontare con le righe 'gap isotropo'. Cambiare l'ordine dei passi e' un cambio "
        "di regola: va deciso con Valerio e coordinato con tutti i thread (non lo faccio io)."),
      "bounds": BOUNDS,
      "conditional_table": table(),
    }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "lvloop_prereg.json")
    json.dump(pre, open(out, "w"), indent=1, ensure_ascii=False, default=float)
    for r in pre["conditional_table"]:
        print(f'{r["test"][:55]:55s} p={r["exponent"]} B={r["bound"]:.1e}  {r["verdict_mass_type"]}'
              + (f' | impulso: {r["verdict_momentum_type"]}' if "verdict_momentum_type" in r else ""))


if __name__ == "__main__":
    main()

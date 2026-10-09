"""trap_eh.py - Parte B2a: la non linearita' VERA del vuoto (Euler-Heisenberg) puo' legare la luce?

Fisica usata (QED a un loop, campi lenti rispetto a m c^2/hbar):
  L_EH = (2 alpha^2 / 45 m^4) [ (E^2-B^2)^2 + 7 (E.B)^2 ]      (unita' naturali, Euler-Heisenberg 1936)
  => una sonda che attraversa un campo forte vede un indice n = 1 + kappa (alpha/45 pi) (E/E_S)^2,
     kappa = 2 / 3.5 (campo B statico, polarizzazioni perpendicolare/parallela, Adler 1971)
     kappa fino a ~ 4 x (11+3) = 56 per un'onda contro-propagante (Heinzl et al. 2006) -> usiamo
     kappa = 56, il caso PIU' FAVOREVOLE, per non sottovalutare la QED.
  Un'onda piana da sola: E^2-B^2 = 0 e E.B = 0  =>  nessuna auto-interazione (zero esatto).
  E_S = m^2 c^3/(e hbar) = 1.32e18 V/m (campo di Schwinger).
  Tasso di creazione di coppie (Schwinger 1951, campo costante):
     w = (e E)^2 / (4 pi^3 hbar^2 c) * sum_n n^-2 exp(-n pi E_S/E)     [coppie / (m^3 s)]

Condizione per intrappolare luce di lunghezza d'onda lam in una regione di raggio a:
  una guida (fibra) con salto d'indice Dn ha un modo legato se V = (2 pi a/lam) sqrt(2 Dn) >= 2.405
  => Dn_needed(a) = 0.5 (2.405 lam / (2 pi a))^2   (in 3-D serve almeno questo in ogni direzione).

uso: python lab/trap_eh.py
"""
import json, os
import numpy as np

c = 299792458.0; hbar = 1.054571817e-34; h = 2 * np.pi * hbar; e = 1.602176634e-19
eps0 = 8.8541878128e-12; me = 9.1093837015e-31; alpha = 7.2973525693e-3
mc2 = me * c**2
ES = me**2 * c**3 / (e * hbar)
lamC = h / (me * c); lambdabar_C = hbar / (me * c)
KAPPA_MAX = 56.0
OUT = os.path.join(os.path.dirname(__file__), "out")


def dn_EH(E, kappa=KAPPA_MAX):
    return kappa * alpha / (45 * np.pi) * (E / ES) ** 2


def E_for_dn(dn, kappa=KAPPA_MAX):
    return ES * np.sqrt(dn * 45 * np.pi / (kappa * alpha))


def dn_needed(a_over_lam):
    return 0.5 * (2.405 / (2 * np.pi * a_over_lam)) ** 2


def schwinger_rate(E):
    nsum = sum(np.exp(-n * np.pi * ES / E) / n**2 for n in range(1, 200))
    return (e * E) ** 2 / (4 * np.pi**3 * hbar**2 * c) * nsum


res = {}
print(f"E_S = {ES:.4e} V/m ; I_S = c eps0 E_S^2/2 = {c*eps0*ES**2/2:.3e} W/m^2 = {c*eps0*ES**2/2/1e4:.3e} W/cm^2")
print(f"lambda_C = {lamC:.4e} m, lambdabar_C = {lambdabar_C:.4e} m")
res["E_S"] = ES

print("\n(1) Indice di rifrazione del vuoto a campi realistici (kappa = 56, il caso piu' favorevole)")
rows = []
for name, E in [("campo nell'atomo di H a r = a0 (5.1e11 V/m)", 5.14e11),
                ("laser record 1.1e23 W/cm^2 (Yoon et al. 2021)", np.sqrt(2 * 1.1e27 / (c * eps0))),
                ("magnetar B = 1e11 T  (cB)", 1e11 * c),
                ("0.1 E_S", 0.1 * ES), ("E_S", ES)]:
    dn = dn_EH(E)
    a_need = 2.405 / (2 * np.pi * np.sqrt(2 * dn))   # raggio minimo di intrappolamento, in lambda
    rows.append(dict(name=name, E=E, E_over_ES=E / ES, dn=dn, a_min_over_lambda=a_need))
    print(f"  {name:48s} E/E_S={E/ES:9.2e}  Dn={dn:9.2e}  raggio minimo per legare = {a_need:9.2e} lambda")
res["realistic"] = rows

print("\n(2) Quanto campo serve per legare la luce in una regione di raggio a = s*lambda?")
rows = []
for s in [0.5, 1, 2, 5, 10]:
    dn = dn_needed(s)
    E = E_for_dn(dn)
    rate = schwinger_rate(E)
    u = 0.5 * eps0 * E**2
    tau = u / (rate * 2 * mc2)       # tempo per convertire tutta l'energia del campo in coppie
    rows.append(dict(s=s, dn_needed=dn, E_over_ES=E / ES, pair_depletion_time_s=tau,
                     tau_in_compton_times=tau / (hbar / mc2)))
    print(f"  a = {s:4.1f} lambda: Dn = {dn:8.2e} -> E = {E/ES:6.2f} E_S ; il campo si svuota in coppie in "
          f"{tau:8.2e} s = {tau/(hbar/mc2):8.2e} tempi Compton (hbar/mc^2 = {hbar/mc2:.2e} s)")
res["needed"] = rows

print("\n(3) L'elettrone come UN fotone di energia m c^2 chiuso nel proprio volume")
rows = []
for label, Lside in [("cubo di lato lambda_C", lamC), ("cubo di lato lambdabar_C", lambdabar_C),
                     ("cubo di lato 1e-18 m (limite sperimentale)", 1e-18)]:
    u = mc2 / Lside**3
    E = np.sqrt(2 * u / eps0)
    dn = dn_EH(E)
    a_over_lam = Lside / 2 / lamC
    need = dn_needed(a_over_lam)
    rows.append(dict(label=label, L=Lside, E_over_ES=E / ES, dn_EH=dn, dn_needed=need, ratio=need / dn))
    print(f"  {label:44s}: E = {E/ES:9.3e} E_S, Dn(EH) = {dn:9.2e}, Dn necessario = {need:9.2e}  "
          f"-> manca un fattore {need/dn:9.2e}")
res["electron_as_photon"] = rows
print("  (a 1e-18 m l'espansione di Euler-Heisenberg non vale piu': E >> E_S e hbar*omega >> m c^2;"
      " il campo produrrebbe coppie, non un oggetto stabile)")

print("\n(4) Un fotone di lunghezza d'onda lam chiuso in un cubo di lato lam: per quale lam si lega da solo?")
# Dn(EH, u = h c/lam / lam^3) = Dn_needed(0.5)  ->  risolvi in lam
lam = np.logspace(-16, -10, 6001)
E = np.sqrt(2 * (h * c / lam) / lam**3 / eps0)
ratio = dn_EH(E) / dn_needed(0.5)
i = np.argmin(abs(np.log(ratio)))
print(f"  servirebbe lam ~ {lam[i]:.2e} m (energia {h*c/lam[i]/e/1e6:.1f} MeV = {h*c/lam[i]/mc2:.0f} m_e c^2),"
      f" campo {E[i]/ES:.1f} E_S: oltre la soglia di coppie (1.022 MeV) e fuori dalla validita' di EH")
res["self_binding_lambda"] = dict(lam=lam[i], energy_MeV=h * c / lam[i] / e / 1e6, E_over_ES=E[i] / ES)

print("\n(5) Gravita' contro elettromagnetismo alla scala dell'elettrone")
G = 6.67430e-11
ratio = G * me**2 / (e**2 / (4 * np.pi * eps0))
print(f"  G m_e^2 / (e^2/4 pi eps0) = {ratio:.3e}   (vale a ogni distanza: entrambi ~ 1/r)")
rs = 2 * G * me / c**2
print(f"  raggio di Schwarzschild dell'elettrone 2Gm/c^2 = {rs:.2e} m  ;  lunghezza di Planck = {np.sqrt(hbar*G/c**3):.2e} m")
res["grav_over_em"] = ratio
res["schwarzschild_e"] = rs
# Einstein-Cartan: l'interazione spin-spin efficace ha densita' ~ G hbar^2 n^2 / c^2 (Hehl 1976).
# diventa confrontabile con la densita' di energia di massa m c^2 n quando n ~ m c^4/(G hbar^2) -> scala di Cartan
lc = (lamC / (2 * np.pi) * (np.sqrt(hbar * G / c**3)) ** 2) ** (1 / 3)
print(f"  scala di Cartan (dove la torsione di Einstein-Cartan conta per un elettrone) "
      f"r_C ~ (lambdabar_C l_P^2)^(1/3) = {lc:.2e} m")
res["cartan_length"] = lc
# accoppiamento che servirebbe: forza 'geometrica' di intensita' alpha_g  =  G_eff m^2/hbar c ~ alpha
Geff = alpha * hbar * c / me**2
print(f"  per avere un effetto geometrico forte quanto l'EM su m_e servirebbe G_eff ~ alpha hbar c/m_e^2 = {Geff:.2e}"
      f" = {Geff/G:.2e} G")
res["G_eff_needed_over_G"] = Geff / G

with open(os.path.join(OUT, "trap_eh.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
print("\nscritto lab/out/trap_eh.json")

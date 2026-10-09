"""trap_atom_energy.py - Parte C: mappa ENERGETICA dell'atomo (fisica stabilita, nessuna speculazione).

Idrogeno (e, se veloce, elio-4): densita' di energia radiale da 1e-18 m a 1e-9 m.
Componenti (tutte in J/m^3, poi "energia per decade di raggio" dE/dln r = 4 pi r^3 u):
  * protone/nucleo: energia di riposo M c^2 distribuita con un profilo di massa
      H:  esponenziale (profilo "dipolo") con raggio quadratico medio di massa 0.755 fm
          (Duran et al., Nature 615, 813 (2023), da J/psi a GlueX; dipende dal modello: il raggio
          scalare dello stesso lavoro e' ~1.0 fm). Carica: esponenziale con r_E = 0.84075 fm (CODATA 2022).
      He: gaussiana con r = 1.67824 fm (Krauth et al., Nature 589, 527 (2021)) per carica E massa
          (approssimazione dichiarata); energia di legame nucleare 28.296 MeV mostrata come "di cui".
  * elettrone: energia di riposo m c^2 distribuita come |psi_1s|^2 (la "nuvola" della MQ; in Bohm
    l'elettrone e' in un punto, ma la sua probabilita' e' la nuvola) ; He: 2 elettroni, Z_eff = 27/16.
  * energia cinetica dell'elettrone (hbar^2/2m)|grad psi|^2  (integrale = 13.6 eV per H)
  * campo elettrico: u = eps0 E^2/2 con E dalla carica del nucleo meno la nuvola elettronica racchiusa.
  * per confronto (tratteggio): campo classico di un elettrone PUNTIFORME fuori da r (diverge come 1/r).
uso: python lab/trap_atom_energy.py
"""
import json, os
import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "out")
c = 299792458.0; e = 1.602176634e-19; eps0 = 8.8541878128e-12; hbar = 1.054571817e-34
me = 9.1093837015e-31; MeV = 1e6 * e; eV = e
alpha = 7.2973525693e-3
a0 = 5.29177210903e-11
mp_c2 = 938.272088 * MeV; me_c2 = 0.51099895 * MeV
r = np.geomspace(1e-19, 1e-8, 40001)
lnr = np.log(r)


def integ(u):              # energia totale da densita' u(r) (simmetria sferica)
    return np.trapezoid(4 * np.pi * r**3 * u, lnr)


def cum(u):
    f = 4 * np.pi * r**3 * u
    return np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(lnr))])


def exp_profile(rrms):     # densita' normalizzata ~ exp(-r/a), <r^2> = 12 a^2
    a = rrms / np.sqrt(12)
    return np.exp(-r / a) / (8 * np.pi * a**3), (lambda x: 1 - np.exp(-x / a) * (1 + x / a + (x / a) ** 2 / 2))


def gauss_profile(rrms):   # <r^2> = 3 s^2
    s = rrms / np.sqrt(3)
    from scipy.special import erf
    dens = np.exp(-r**2 / (2 * s**2)) / ((2 * np.pi) ** 1.5 * s**3)
    enc = lambda x: erf(x / (np.sqrt(2) * s)) - np.sqrt(2 / np.pi) * (x / s) * np.exp(-x**2 / (2 * s**2))
    return dens, enc


def atom(name, Z, M_c2, mass_prof, charge_prof, n_el, a_orb, B_nuc_MeV, E_bind_eV):
    rho_M = M_c2 * mass_prof
    # elettroni 1s: densita' di numero n_el * exp(-2r/a)/(pi a^3)
    ne = n_el * np.exp(-2 * r / a_orb) / (np.pi * a_orb**3)
    rho_e = me_c2 * ne
    enc_e = n_el * (1 - np.exp(-2 * r / a_orb) * (1 + 2 * r / a_orb + 2 * (r / a_orb) ** 2))
    # cinetica: (hbar^2/2m)|dpsi/dr|^2 per elettrone, psi = exp(-r/a)/sqrt(pi a^3)  -> |psi'|^2 = psi^2/a^2
    kin = n_el * hbar**2 / (2 * me) * np.exp(-2 * r / a_orb) / (np.pi * a_orb**3) / a_orb**2
    Qn = Z * charge_prof(r)
    Enet = e * (Qn - enc_e) / (4 * np.pi * eps0 * r**2)
    Enuc = e * Qn / (4 * np.pi * eps0 * r**2)
    u_net = 0.5 * eps0 * Enet**2
    u_nuc = 0.5 * eps0 * Enuc**2
    u_point_e = 0.5 * eps0 * (e / (4 * np.pi * eps0 * r**2)) ** 2
    comps = dict(nucleus_rest=rho_M, electron_rest=rho_e, electron_kinetic=kin, field_net=u_net,
                 field_nucleus_alone=u_nuc, point_electron_field=u_point_e)
    tot = {k: integ(v) for k, v in comps.items()}
    # coda analitica oltre r = 10 nm del campo del nucleo da solo (il campo netto li' e' ~0)
    tot["field_nucleus_alone"] += Z**2 * e**2 / (8 * np.pi * eps0 * r[-1])
    tot["point_electron_field"] += e**2 / (8 * np.pi * eps0 * r[-1])
    # energie di Coulomb (classiche) per controllo: U_net - U_nuc = U_ee(autoenergia nuvola) + U_en
    U_en = -Z * e**2 / (4 * np.pi * eps0) * np.trapezoid(4 * np.pi * r**3 * ne / r * charge_prof_weight(r, charge_prof), lnr)
    regions = [("r < 3 fm (nucleo)", 0, 3e-15), ("3 fm - 1 pm", 3e-15, 1e-12), ("1 pm - 10 nm (atomo)", 1e-12, 1e-8)]
    reg = {}
    for lab, lo, hi in regions:
        msk = (r >= lo) & (r < hi)
        reg[lab] = {k: float(np.trapezoid((4 * np.pi * r**3 * v)[msk], lnr[msk])) / eV for k, v in comps.items()}
    E_total = M_c2 + n_el * me_c2 - E_bind_eV * eV
    print(f"\n=== {name} ===")
    print(f"  energia totale (riposo) = {E_total/MeV:.6f} MeV")
    for k, v in tot.items():
        print(f"  {k:22s}: {v/eV:14.6e} eV  = {v/E_total:10.3e} del totale")
    print(f"  campo netto - campo del nucleo da solo = {(tot['field_net']-tot['field_nucleus_alone'])/eV:.3f} eV "
          f"(= autoenergia classica della nuvola + interazione nucleo-nuvola; la MQ NON conta la prima)")
    if n_el == 1:
        print(f"     atteso: +(5/16) e^2/(4 pi eps0 a) = {5/16*e**2/(4*np.pi*eps0*a_orb)/eV:.3f} eV  "
              f"- e^2/(4 pi eps0 a) = {-e**2/(4*np.pi*eps0*a_orb)/eV:.3f} eV")
    print(f"  energia di Coulomb nucleo-elettroni U_en = {U_en/eV:.3f} eV (H: -27.21 attesi);"
          f" cinetica = {tot['electron_kinetic']/eV:.3f} eV")
    print(f"  campo classico di un elettrone puntiforme fuori da 1e-18 m: "
          f"{np.trapezoid((4*np.pi*r**3*u_point_e)[r >= 1e-18], lnr[r >= 1e-18])/MeV:.1f} MeV (contro m c^2 = 0.511 MeV)")
    print("  per regione (eV):")
    for lab in reg:
        d = reg[lab]
        print(f"    {lab:22s} nucleo {d['nucleus_rest']:10.4e}  e-riposo {d['electron_rest']:10.4e}  "
              f"e-cinetica {d['electron_kinetic']:10.4e}  campo netto {d['field_net']:10.4e}")
    if B_nuc_MeV:
        print(f"  di cui nel nucleo: energia di legame nucleare -{B_nuc_MeV} MeV = {B_nuc_MeV*MeV/E_total:.4%} del totale")
    print(f"  energia di legame elettronica -{E_bind_eV} eV = {E_bind_eV*eV/E_total:.3e} del totale")
    np.savez_compressed(os.path.join(OUT, f"trap_atom_energy_{name}.npz"), r=r[::20],
                        **{k: v[::20] for k, v in comps.items()},
                        totals_eV=np.array([tot[k] / eV for k in comps]), names=np.array(list(comps)))
    return dict(r=r, comps=comps, tot={k: v / eV for k, v in tot.items()}, regions=reg, E_total_MeV=E_total / MeV,
                U_en_eV=U_en / eV)


def charge_prof_weight(rr, enc):
    return enc(rr)


res = {}
mH, encH = exp_profile(0.755e-15)
_, encHc = exp_profile(0.84075e-15)
a_red = a0 * (1 + me / 1.67262192369e-27)
res["H"] = atom("H", 1, mp_c2, mH, encHc, 1, a_red, 0, 13.598434)
mHe, encHe = gauss_profile(1.67824e-15)
res["He4"] = atom("He4", 2, 3727.3794066 * MeV, mHe, encHe, 2, a0 / (27 / 16), 28.29566, 79.00515)

# autoenergia QED dell'elettrone (log-divergente): dm/m = (3 alpha / 2 pi) ln(Lambda/m)
for lab, Lam in (("1 TeV", 1e6), ("scala di Planck 1.22e19 GeV", 1.22e22)):
    print(f"  autoenergia QED dell'elettrone con taglio {lab}: dm/m = {3*alpha/(2*np.pi)*np.log(Lam/0.511):.3f}")

# ------------------------------------------------------------------ figura
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
style = dict(nucleus_rest=("nucleo: energia di riposo", "#1f4e9c", "-"),
             electron_rest=("elettrone/i: energia di riposo (nuvola)", "#c0392b", "-"),
             electron_kinetic=("elettrone: energia cinetica", "#e67e22", "-"),
             field_net=("campo elettrico netto eps0 E^2/2", "#27ae60", "-"),
             field_nucleus_alone=("campo del nucleo da solo", "#27ae60", ":"),
             point_electron_field=("campo di un elettrone PUNTIFORME (classico)", "#7f8c8d", "--"))
for ax, (nm, d) in zip(axs, res.items()):
    rr = d["r"]
    sel = (rr >= 1e-18) & (rr <= 1e-9)
    for k, v in d["comps"].items():
        lab, col, ls = style[k]
        y = (4 * np.pi * rr**3 * v / eV)[sel]
        ax.loglog(rr[sel], np.where(y > 1e-12, y, np.nan), ls, color=col, lw=1.8, label=lab)
    ax.set_xlabel("raggio r [m]")
    ax.set_title(f"{'Idrogeno' if nm=='H' else 'Elio-4'}: energia per decade di raggio")
    ax.axvline(0.84e-15 if nm == "H" else 1.68e-15, color="k", lw=0.6, alpha=0.4)
    ax.axvline(5.29e-11 if nm == "H" else 5.29e-11 / 1.6875, color="k", lw=0.6, alpha=0.4)
    ax.set_ylim(1e-6, 1e10)
    ax.grid(True, which="major", alpha=0.3)
axs[0].axhline(13.6, color="purple", lw=1, ls="-.", label="legame dell'elettrone: 13.6 eV (totale)")
axs[1].axhline(79.0, color="purple", lw=1, ls="-.", label="legame dei 2 elettroni: 79.0 eV (totale)")
axs[1].axhline(28.3e6, color="navy", lw=1, ls="-.", label="legame nucleare: 28.3 MeV (totale)")
axs[1].legend(fontsize=8, loc="upper right", handles=axs[1].get_lines()[-2:])
axs[0].set_ylabel("dE / d(ln r)  [eV]")
axs[0].legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "trap_atom_energy.png"), dpi=110)
print("\nscritte lab/out/trap_atom_energy_H.npz, _He4.npz, trap_atom_energy.png")
with open(os.path.join(OUT, "trap_atom_energy.json"), "w") as f:
    json.dump({k: dict(tot=v["tot"], regions=v["regions"], E_total_MeV=v["E_total_MeV"], U_en_eV=v["U_en_eV"])
               for k, v in res.items()}, f, indent=1)

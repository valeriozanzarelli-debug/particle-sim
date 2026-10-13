"""multi_e_figure.py - figura dello studio MOLTI ELETTRONI (M0 + elio): lab/out/multi_e_figure.png"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from multi_e_core import OUT, THETA, load_json, Band, chain_levels, A0, ALPHA, dirs26

K = load_json("multi_e_kernel.json"); L = load_json("multi_e_lattice.json"); E = load_json("multi_e_ecg.json")
M = load_json("multi_e_mf.json"); Dd = load_json("multi_e_disp.json"); Q = load_json("multi_e_dmc.json")
fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
C1, C2, C3, C4 = "#2a6fdb", "#d9532b", "#2e9d5b", "#7a4fc4"

# (a) R4: stessa regola, due direzioni del tempo
a = ax[0, 0]
E0 = L["det"]["E_lobpcg"]
bs = [2e-4, 1e-4, 5e-5]; Eb = [L["det"][f"beta{b:g}"]["E_growth"] for b in bs]
bs += [2.5e-5, 1.25e-5]; Eb += [L["diag_R4a"][f"beta{b:g}"]["E_growth"] for b in (2.5e-5, 1.25e-5)]
rt = L["diag_R4a"]["realtime"]; ts = np.array(sorted(float(t) for t in rt)); Er = np.array([rt[f"{t:g}"] for t in ts])
a.plot(np.array(bs) ** 2 * 1e8, np.array(Eb) - E0, "o-", c=C1, label="tempo immaginario: $-\\ln\\lambda/\\beta$")
a.plot(ts ** 2 * 1e8, Er - E0, "s-", c=C2, label="tempo reale: $-\\arg\\langle\\psi|W|\\psi\\rangle/\\tau$")
xx = np.linspace(0, 4.2, 50)
a.plot(xx, L["diag_R4a"]["c2_imag"] * xx * 1e-8, ":", c=C1); a.plot(xx, L["diag_R4a"]["c2_real"] * xx * 1e-8, ":", c=C2)
a.axhline(0, c="k", lw=0.8)
a.set_xlabel("$\\beta^2$ o $\\tau^2$  [$10^{-8}$ a.u.$^2$]"); a.set_ylabel("$E - E_{LOBPCG}$  [Eh]")
a.set_title(f"R4 ione Z=16 sul reticolo imparato\nWick: c2 immag. {L['diag_R4a']['c2_imag']:.3g}, reale {L['diag_R4a']['c2_real']:.3g}", fontsize=10)
a.set_xlim(0, 4.2); a.set_ylim(-1.0, 1.0); a.legend(fontsize=8)

# (b) R2: peso negativo del nucleo
a = ax[0, 1]
betas = [1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 1e-3, 1e-2, 3e-2]
for nm, c, lab in (("B16", C1, "livello 1 (a0/256)"), ("B256", C2, "a0/16"), ("B512", C3, "a0/8")):
    y = [max(K["kbeta"][nm][f"n128_b{b:g}"]["neg"], 1e-17) for b in betas]
    a.loglog(betas, y, "o-", c=c, label=lab)
a.axhline(1e-10, c="k", ls="--", lw=0.8); a.axvline(K["R2"]["summary"]["B16"]["beta_min"], c=C1, ls=":")
a.set_xlabel("$\\beta$ [a.u.]"); a.set_ylabel("peso negativo di $K_\\beta$")
a.set_title("R2: dove il tempo immaginario e' una probabilita'", fontsize=10); a.legend(fontsize=8)

# (c) contenuto k^4 della banda imparata
a = ax[0, 2]
b1 = Band(chain_levels("A")[0]["K"], chain_levels("A")[0]["offs"], 16, THETA["A"], sym=True)
m = K["R3"]["m_kin_theta_A_used"]
kk = np.geomspace(2, 60, 40)
for nm, d, ls in (("asse", np.array([1., 0, 0]), "-"), ("faccia", np.array([1., 1, 0]) / 2 ** .5, "--"), ("corpo", np.ones(3) / 3 ** .5, ":")):
    e = b1.eps(kk[:, None] * d[None])
    a.plot(kk, (e - kk ** 2 / (2 * m)) / (-ALPHA ** 2 * kk ** 4 / 8), ls, c=C1, label=f"livello 1 esistente, {nm}")
a.axhline(1, c="k", lw=0.8)
for nm, c in (("relearn_A", C3), ("relearn_A_pauli_rejected", C2)):
    an = Dd["bands"][nm]["anisotropy"]
    a.plot([10, 30], [an["k4frac_mean_k10"], an["k4frac_mean_k30"]], "o", c=c,
               label={"relearn_A": "ri-imparato SCELTO (dirac4, B=2)", "relearn_A_pauli_rejected": "miglior Pauli b=16 (RIFIUTATO)"}[nm])
a.set_ylim(-0.2, 1.3); a.set_xlabel("k [1/a0]"); a.set_ylabel("$(\\varepsilon-k^2/2m)/(-\\alpha^2k^4/8)$")
a.set_title("quanta relativita' (k^4) porta la banda", fontsize=10); a.legend(fontsize=7)

# (d) R3 m_kin(theta)
a = ax[1, 0]
th = np.linspace(0, 0.28, 100)
a.plot(th ** 2, np.tan(th) / th - 1, "k-", lw=0.8, label="legge della regola fine tan(th)/th")
for ch, c in (("A2", C3), ("A", C1), ("B", C2)):
    for sd in (1, 2):
        r = K["mkin"][f"{ch}_s{sd}_sym"]
        a.plot(r["theta"] ** 2, r["m_kin"] - 1, "o", c=c, ms=6 if sd == 1 else 4)
f = np.mean(K["R3"]["fits"], 0)
a.plot(th ** 2, f[0] - 1 + f[1] * th ** 2 + f[2] * th ** 4, "--", c=C4, label=f"fit: m_kin(0)-1 = {K['R3']['m0'] - 1:.1e}")
a.set_xscale("symlog", linthresh=1e-4); a.set_yscale("symlog", linthresh=1e-5)
a.set_xlabel("$\\theta^2$"); a.set_ylabel("$m_{kin}-1$"); a.set_title("R3: la massa imparata va a 1", fontsize=10); a.legend(fontsize=8)

# (e) G1: Hooke e condizionamento
a = ax[1, 1]
rows = E["G1_diag_posthoc"]["Hooke"]["rows"]
a.semilogy([r["K"] for r in rows], [abs(r["E_minus_2"]) for r in rows], "o-", c=C1, label="Hooke |E-2| (post hoc)")
a.semilogy([200], [abs(E["G1"]["Hooke"]["err"])], "x", c=C2, ms=10, label="G1 pre-registrato (K=200)")
a.axhline(1e-8, c="k", ls="--", lw=0.8)
a2 = a.twinx(); a2.semilogy([r["K"] for r in rows], [r["cond_S"] for r in rows], "s:", c=C4); a2.set_ylabel("cond(S)", color=C4)
a.set_xlabel("K funzioni"); a.set_ylabel("errore [Eh]"); a.set_title("G1: atomo di Hooke (FAIL, poi diagnosi)", fontsize=10)
a.legend(fontsize=8, loc="lower left")

# (f) campo medio: ordine della griglia
a = ax[1, 2]
for nm, c in zip(("He", "H-", "Li", "Be", "Ne"), (C1, C2, C3, C4, "#888")):
    r = M.get("G3", {}).get(nm)
    if not r or "E_n6000" not in r:
        continue
    ns = np.array([1500, 3000, 6000]); e = np.array([r[f"E_n{n}"] for n in ns])
    a.loglog(ns, np.abs(e - r["E_richardson"]) / abs(r["target"]), "o-", c=c, label=f"{nm}: p = {r['order_p']:.2f}, Richardson {r['rel_err']:+.1e}")
a.axhline(2e-5, c="k", ls="--", lw=0.8)
a.set_xlabel("punti della griglia"); a.set_ylabel("|E(n) - E_Richardson| / |E|"); a.set_title("G3: campo medio radiale", fontsize=10)
a.legend(fontsize=7)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "multi_e_figure.png"), dpi=110)
print("figura scritta")

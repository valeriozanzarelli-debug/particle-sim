"""assist_figure.py - figure of the assist_ study -> lab/out/assist_figura.png"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assist_core import S_wi, OUT, PI

L = lambda n: json.load(open(os.path.join(OUT, f"assist_{n}.json")))
V, AN, SA, pre = L("val"), L("analysis"), L("specA"), L("prereg")
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

# (a) Keldysh crossover, single Sauter pulse: rule vs exact
a = ax[0, 0]
for key, col in (("g0.1", "C0"), ("g1.0", "C1"), ("g10.0", "C2")):
    r = V["V2"][key]
    p = np.array(r["curve_p"]); fe = np.array(r["curve_exact"]); fl = np.array(r["curve_lat"])
    a.semilogy(p, fe, "-", color=col, lw=1, label=f"esatto, γ = {r['gamma']} (E = {r['E']} E_S)")
    a.semilogy(p[::4], fl[::4], "o", color=col, ms=3)
a.set_xlabel("impulso p / m"); a.set_ylabel("coppie per modo f(p)")
a.set_title("V2: un impulso solo, da tunnel (γ≪1) a multifotone (γ≫1)\nlinee: formula esatta; punti: regola zig-zag (a = 0,025/m)")
a.set_ylim(1e-24, 1e-3); a.legend(fontsize=8)

# (b) exponent vs gamma
a = ax[0, 1]
gg = np.linspace(0.3, 11, 400)
S0 = S_wi(0.1, 0.2)
a.plot(gg, [S_wi(g, 0.2) / S0 for g in gg], "k-", lw=1, label="istantone (ε→0), Schützhold–Gies–Dunne")
gs = pre["A_sauter"]["gammas"]
a.plot(gs, [AN["A2"][str(g)]["S_fit"] / S0 for g in gs], "o", color="C3", label="regola: esponente di Q (impulso di Sauter)")
gc = AN["A3"]["gc_fit"]
a.axvline(pre["A_sauter"]["predictions"]["gamma_crit_gs"], color="k", ls=":", lw=1)
a.axvline(gc, color="C3", ls="--", lw=1, label=f"soglia dal fit γ_c = {gc:.3f} (previsto 1,592)")
gco = pre["O_cos"]["gammas"]
S0o = AN["O2_O3"]["S0_fit"]["S"]
a.plot(gco, [AN["O2_O3"][f"eps0.01_g{g}"]["S_fit"] / S0o for g in gco], "s", color="C0", ms=4,
       label="regola: coseno debole ε = 0,01")
a.plot(gco, [pre["O_cos"]["predictions"][f"eps0.01_g{g}"]["s_ratio_WKB"] for g in gco], "--", color="C0", lw=1,
       label="WKB coseno ε = 0,01")
a.set_xlabel("γ = m ω / (e E1)  (Keldysh combinato)"); a.set_ylabel("esponente / esponente di Schwinger")
a.set_title("A2–A3: la soglia netta in frequenza"); a.legend(fontsize=7); a.set_ylim(0, 1.1)

# (c) spectra at E1 = 0.1, eps = 0.1
a = ax[1, 0]
p = np.array(SA["p"])
rich = lambda d: (4 * np.array(d["th0.025"]) - np.array(d["th0.05"])) / 3
a.semilogy(p, rich(SA["strong_only"]), "k-", lw=1.5, label="solo campo forte E1 = 0,1 E_S")
for g, col in ((1.0, "C0"), (2.0, "C1"), (3.0, "C2"), (6.0, "C3")):
    a.semilogy(p, rich(SA[f"eps0.1_g{g}"]["both"]), "-", color=col, lw=1, label=f"forte + debole (ε = 0,1), γ = {g}")
a.set_xlabel("impulso p / m"); a.set_ylabel("coppie per modo f(p)")
a.set_title("spettri della regola (estrapolati a → 0)"); a.legend(fontsize=8); a.set_ylim(1e-16, 1e-4)

# (d) enhancement
a = ax[1, 1]
for name, eps_list, gl, mk in (("specA", [0.01, 0.1], pre["A_sauter"]["gammas"] + [20.0], "o"),
                               ("specO", [0.001, 0.01], pre["O_cos"]["gammas"], "s")):
    for eps in eps_list:
        y = [AN[name][f"eps{eps}_g{g}"]["enhancement"] for g in gl]
        a.semilogy(gl, y, mk + "-", ms=4, lw=1,
                   label=("Sauter" if name == "specA" else "coseno") + f", ε = {eps}")
a.axhline(1, color="k", lw=0.5)
a.set_xlabel("γ"); a.set_ylabel("N(insieme) / [N(forte) + N(debole)]")
a.set_title("aumento del numero di coppie, E1 = 0,1 E_S (1+1 D)"); a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "assist_figura.png"), dpi=110)
print("ok")

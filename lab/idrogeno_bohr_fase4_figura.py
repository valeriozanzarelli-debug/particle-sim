"""Figura della fase 4 dell'idrogeno di Bohr (docs/IDROGENO_BOHR.md §13): solo bagliori (palline morbide), niente linee.
Legge lab/out/idrogeno_bohr_fase4_esame.json, _grosso_esame.json, _regola.json. Scrive lab/out/idrogeno_bohr_fase4_figura.png."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
E = json.load(open(os.path.join(OUT, "idrogeno_bohr_fase4_esame.json")))["numeri"]
G = json.load(open(os.path.join(OUT, "idrogeno_bohr_fase4_grosso_esame.json")))
R = json.load(open(os.path.join(OUT, "idrogeno_bohr_fase4_regola.json")))
BG, FG = "#07080d", "#d8dbe6"


def glow(ax, x, y, col, s=60):
    for k, a in ((9.0, 0.05), (5.0, 0.09), (2.6, 0.18), (1.4, 0.4), (0.7, 0.95)):
        ax.scatter(x, y, s=s * k, color=col, alpha=a, edgecolors="none", zorder=3)
    ax.scatter(x, y, s=s * 0.25, color="white", alpha=0.85, edgecolors="none", zorder=4)


def stile(ax, t):
    ax.set_facecolor(BG)
    ax.set_title(t, color=FG, fontsize=10.5)
    ax.tick_params(colors=FG, labelsize=8)
    for s in ax.spines.values():
        s.set_color("#3a3f52")
    ax.xaxis.label.set_color(FG); ax.yaxis.label.set_color(FG)


fig, axs = plt.subplots(2, 2, figsize=(11, 8.2), facecolor=BG)
# (a) massa del maestro
ax = axs[0, 0]
js = sorted(int(j) for j in R["masse_per_j"])
m1 = [R["masse_per_j"][str(j)] - 1 for j in js]
glow(ax, js, m1, "#ffb347")
glow(ax, [R["maestro_j"]], [m1[js.index(R["maestro_j"])]], "#7cff9e", s=110)
ax.set_yscale("log"); ax.set_xticks(js); ax.set_xticklabels([f"a0/{4096 * 2 ** j}" for j in js], fontsize=7.5)
ax.set_ylabel("m − 1 della regola a0/16"); ax.set_xlabel("cella del maestro")
stile(ax, "(a) massa del maestro: l'agglomeratore affina finché |Δm| < 1e-5")
# (b) convergenza nel blocco
ax = axs[0, 1]
V = E.get("V", {})
nis = sorted(int(k) for k in V)
glow(ax, nis, [1e4 * V[str(k)]["ibrido_su_grosso"] for k in nis], "#6fd3ff")
vecchio = {3: 0.72, 5: 2.37, 7: 4.88, 9: 3.44, 11: 2.60}      # sviluppo: bordo vecchio, 48 passi, regola T_v1
glow(ax, list(vecchio), list(vecchio.values()), "#ff6f91", s=40)
ax.set_xlabel("blocco fine (celle grosse per lato)"); ax.set_ylabel("ibrido / grosso − 1  [1e-4]")
ax.set_xticks([3, 5, 7, 9, 11]); ax.set_ylim(-0.3, 7.2)
ax.text(0.03, 0.93, "azzurro: fase 4 (bordo giusto, partenza veff, 144 passi)", color="#6fd3ff", transform=ax.transAxes, fontsize=8)
ax.text(0.03, 0.86, "rosa: bordo vecchio (sviluppo, 48 passi)", color="#ff6f91", transform=ax.transAxes, fontsize=8)
stile(ax, "(b) l'energia dell'ibrido al crescere del blocco (N = 128)")
# (c) livelli contro Dirac
ax = axs[1, 0]
nomi = ["1s1/2", "2s1/2", "2p1/2", "2p3/2"]
prima = [0.0283, 0.0347, 0.0377, 0.0377]          # AGGLOMERATO_2S, E2XT (regola T col maestro a0/4096)
dopo = [100 * G["livelli"][k]["rel"] for k in nomi]
x = np.arange(4)
glow(ax, x - 0.12, prima, "#ff6f91")
glow(ax, x + 0.12, dopo, "#7cff9e")
ax.set_xticks(x); ax.set_xticklabels(nomi)
ax.set_ylabel("errore da Dirac [%]")
ax.text(0.03, 0.93, "rosa: prima (maestro a0/4096)", color="#ff6f91", transform=ax.transAxes, fontsize=8)
ax.text(0.03, 0.86, "verde: fase 4 (maestro a0/65536 + bordi ritorti + Born)", color="#7cff9e", transform=ax.transAxes, fontsize=8)
ax.set_ylim(-0.005, 0.056)
stile(ax, "(c) il grosso di default a a0/8: quattro livelli")
# (d) cucitura
ax = axs[1, 1]
cp, ci = E["C_pavimento"], E["C_ibrido"]
lab = ["perdita\nper passo", "stazionarietà\nfuori dal blocco", "corrente\nfuori dal blocco"]
fase3 = [5.1e-6, 3.0e-4, 2.7e-4]
ora = [ci["perdita_mediana_2a_meta"], ci["staz_fuori"], ci["j"][-1][2]]
pav = [np.nan, cp["staz_fuori"], cp["j"][-1][2]]
x = np.arange(3)
glow(ax, x - 0.18, fase3, "#ff6f91")
glow(ax, x, ora, "#6fd3ff")
glow(ax, x[1:] + 0.18, pav[1:], "#9a9fb5", s=40)
ax.set_yscale("log"); ax.set_ylim(1e-8, 5e-3); ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8)
ax.text(0.03, 0.93, "rosa: fase 3   azzurro: fase 4   grigio: sola regola grossa (pavimento)", color=FG, transform=ax.transAxes, fontsize=7.5)
stile(ax, "(d) la cucitura fine/grosso (blocco 7)")
fig.suptitle("Idrogeno di Bohr, fase 4: cucitura, massa del maestro, grosso di default", color=FG, fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(os.path.join(OUT, "idrogeno_bohr_fase4_figura.png"), dpi=130, facecolor=BG)
print("scritto lab/out/idrogeno_bohr_fase4_figura.png")

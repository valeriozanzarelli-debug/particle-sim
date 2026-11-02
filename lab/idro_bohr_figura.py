"""Figura dello studio idrogeno_bohr (docs/IDROGENO_BOHR.md).
Legge i numeri da lab/out/idrogeno_bohr_risultati.json e due fotogrammi dall'export in
/mnt/project-files/ripresa/dati/idrogeno_bohr/ (se manca, i pannelli delle immagini restano vuoti).
Uso: python3 lab/idro_bohr_figura.py"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = json.load(open("lab/out/idrogeno_bohr_risultati.json"))
D = "/mnt/project-files/ripresa/dati/idrogeno_bohr"
fig = plt.figure(figsize=(17, 10.5))
gs = fig.add_gridspec(2, 6, height_ratios=[1.15, 1], hspace=0.45, wspace=0.55)

# (a) errore contro Dirac x mu: stati e livelli
ax = fig.add_subplot(gs[0, 0:2])
e = R["energie_pct"]
nomi = list(e.keys()); val = [e[k] for k in nomi]
col = ["#888" if "grosso" in k else "#2a7" if "ibrido" in k else "#36c" for k in nomi]
ax.barh(range(len(nomi)), val, color=col)
ax.set_yticks(range(len(nomi))); ax.set_yticklabels(nomi, fontsize=8)
ax.axvline(R["massa_rinormalizzata_pct"], color="r", ls="--", lw=1)
ax.text(R["massa_rinormalizzata_pct"], len(nomi) - 0.4, " massa del cammino\n di base (+0,037 %)", color="r", fontsize=7, va="top")
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("errore contro Dirac × massa ridotta (%)")
ax.set_title("(a) energie: chi sbaglia e di quanto", fontsize=10)

# (b) scatola: l'effetto cristallo
ax = fig.add_subplot(gs[0, 2:4])
for lab, pts in R["scatola"].items():
    x, y = zip(*pts); ax.plot(x, y, "o-", label=lab)
ax.set_yscale("symlog", linthresh=0.05); ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("mezzo lato della scatola (a0)"); ax.set_ylabel("errore (%)")
ax.set_title("(b) scatola periodica = cristallo di atomi", fontsize=10); ax.legend(fontsize=7)

# (c) ibrido contro grosso al crescere del blocco fine
ax = fig.add_subplot(gs[0, 4:6])
for lab, pts in R["blocco"].items():
    x, y = zip(*pts); ax.plot(x, np.array(y) * 1e4, "o-", label=lab)
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("lato del blocco fine (celle grosse a0/16)"); ax.set_ylabel("E_ibrido/E_grosso − 1  (×10⁻⁴)")
ax.set_title("(c) quanto pesa il blocco fine (e la colla)", fontsize=10); ax.legend(fontsize=7, loc="upper left")

def leggi(nome, forma):
    p = os.path.join(D, nome)
    if not os.path.exists(p): return None
    return np.fromfile(p, dtype="<f2").astype(np.float32).reshape(forma)

# (d) vista lontana: (1s + 2p_z)/sqrt2, piano y = 0, due istanti
m = R["export"]["lato_lontano"]; F = R["export"]["fotogrammi_lontano"]
rho = leggi("rho_lontano.f16", (F, m, m, m))
for q, (k, tit) in enumerate([(0, "t = 0"), (F // 2, "mezzo periodo")]):
    ax = fig.add_subplot(gs[1, q])
    if rho is not None:
        ext = R["export"]["mezzo_lato_lontano_a0"]
        ax.imshow(np.log10(rho[k, :, m // 2, :].T + 1e-5), origin="lower", extent=[-ext, ext, -ext, ext], cmap="magma", vmin=-4, vmax=-0.5)
    ax.set_title(f"(d) lontano (1s+2p)/√2, {tit}\nlog10 densità, piano y = 0", fontsize=8); ax.set_xlabel("x (a0)", fontsize=8); ax.set_ylabel("z (a0)", fontsize=8)

# (e) zoom: densità e |R|²−|L|² nel blocco fine (vicino al protone)
mf = R["export"]["lato_fine"]
rf = leggi("rho_fine.f16", (32, mf, mf, mf)); cf = leggi("chir_fine.f16", (32, mf, mf, mf))
ef = R["export"]["mezzo_lato_fine_a0"]
ax = fig.add_subplot(gs[1, 2])
if rf is not None: ax.imshow(rf[0, :, mf // 2, :].T, origin="lower", extent=[-ef, ef, -ef, ef], cmap="magma")
ax.set_title("(e) zoom: densità\n(cuspide sul protone)", fontsize=9); ax.set_xlabel("x (a0)", fontsize=8)
ax = fig.add_subplot(gs[1, 3])
if cf is not None:
    v = np.abs(cf[0]).max()
    ax.imshow(cf[0, :, mf // 2, :].T, origin="lower", extent=[-ef, ef, -ef, ef], cmap="RdBu_r", vmin=-v, vmax=v)
ax.set_title("(e) |R|²−|L|², tick 1\n(chi corre su, chi giù, a c)", fontsize=8); ax.set_xlabel("x (a0)", fontsize=8)
ax = fig.add_subplot(gs[1, 4])
if cf is not None:
    d = cf[1, :, mf // 2, :] - cf[0, :, mf // 2, :]; v = np.abs(d).max()
    ax.imshow(d.T, origin="lower", extent=[-ef, ef, -ef, ef], cmap="PuOr", vmin=-v, vmax=v)
ax.set_title("(e) tick 2 − tick 1\n(il tremolio)", fontsize=8); ax.set_xlabel("x (a0)", fontsize=8)

# (f) costi
ax = fig.add_subplot(gs[1, 5])
c = R["costo_s_per_passo"]
ax.bar(range(len(c)), list(c.values()), color="#c63")
ax.set_xticks(range(len(c))); ax.set_xticklabels(list(c.keys()), rotation=60, fontsize=7)
ax.set_ylabel("s per passo grosso (2 thread)", fontsize=8); ax.set_title("(f) costo su CPU", fontsize=9)

fig.suptitle("Idrogeno vero a scala di Bohr: regola imparata lontano, zig-zag vicino al protone (Z = 1)", fontsize=12)
fig.savefig("lab/out/idrogeno_bohr_figure.png", dpi=110)
print("scritto lab/out/idrogeno_bohr_figure.png")

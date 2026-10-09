"""kk_shape_figure.py - one figure for docs/FORMA_NASCOSTA.md -> lab/out/kk_shape.png
(shapes and their mode patterns, mass spectra vs the leptons, level statistics, a scar, the exponential route)"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

md = np.load("out/kk_shape_modes.npz")
st = np.load("out/kk_shape_stats.npz")
lp = np.load("out/kk_shape_leptons.npz")
sp = json.load(open("out/kk_shape_spectra.json"))
le = json.load(open("out/kk_shape_leptons.json"))
sj = json.load(open("out/kk_shape_stats.json"))

shapes = [("disk", "cerchio (disco)"), ("ellipse_2", "ovale 2:1"), ("triangle", "triangolo"), ("stadium_1", "stadio (caotico)")]
fig = plt.figure(figsize=(15, 13))
gs = fig.add_gridspec(4, 6, height_ratios=[1, 1, 1.15, 1.15], hspace=0.45, wspace=0.35)

# rows 1-2: the four shapes, 3 lowest distinct modes each (density on the hidden plane)
for c, (k, title) in enumerate(shapes):
    mask = md[f"{k}__mask"]; dens = md[f"{k}__dens"]; x = md[f"{k}__x"]
    # layout: 4 shapes x 3 modes = 12 panels in rows 0-1 (6 per row): shape c uses cols 3*(c%2)..+2, row c//2
    for r in range(3):
        ax = fig.add_subplot(gs[c // 2, 3 * (c % 2) + r])
        img = np.where(mask, dens[r], np.nan)
        ax.imshow(img.T, origin="lower", cmap="magma", interpolation="nearest")
        ax.contour(mask.T.astype(float), [0.5], colors="c", linewidths=0.8)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{title}\nmodo {r + 1}: m R = {x[r]:.2f}" if r == 0 else f"modo {r + 1}: m R = {x[r]:.2f}", fontsize=9)

# row 3 left: spectra as ratios m_i/m_1 (log) vs the leptons
ax = fig.add_subplot(gs[2, 0:2])
names = ["disk", "ellipse_1.5", "ellipse_2", "ellipse_3", "triangle", "stadium_1"]
for i, k in enumerate(names):
    x = np.array(sp["shapes"][k]["R40"]["x"][:6]); rr = x / x[0]
    ax.scatter([i] * len(rr), rr, marker="_", s=300, color="C0")
ax.axhline(le["targets"]["r2"], color="C3", ls="--"); ax.text(0, le["targets"]["r2"] * 1.15, "muone / elettrone = 206.8", color="C3", fontsize=8)
ax.axhline(le["targets"]["r3"], color="C3", ls=":"); ax.text(0, le["targets"]["r3"] * 1.15, "tau / elettrone = 3477", color="C3", fontsize=8)
ax.axhline(le["E0_theorem"]["mass_ratio_max"], color="gray", lw=0.8)
ax.set_yscale("log"); ax.set_ylim(0.8, 1e4)
ax.set_xticks(range(len(names))); ax.set_xticklabels(["disco", "ov 1.5", "ov 2", "ov 3", "triang.", "stadio"], fontsize=8)
ax.set_ylabel("massa / massa del modo 1"); ax.set_title("masse = modi della forma (cammino, R = 40)\nlinea grigia: limite di Ashbaugh-Benguria", fontsize=9)

# row 3 middle: spacing distributions
ax = fig.add_subplot(gs[2, 2:4])
s = np.linspace(0, 4, 200)
bins = np.linspace(0, 4, 25)
ax.hist(np.concatenate([st["lap_quarter_disk"], st["lap_golden_rectangle"]]), bins=bins, density=True, histtype="step", lw=1.6, color="C0", label="quarto di disco + rettangolo (integrabili)")
ax.hist(st["lap_quarter_stadium"], bins=bins, density=True, histtype="step", lw=1.6, color="C3", label="quarto di stadio (caotico)")
if "walk_w8" in st.files:
    ax.hist(st["walk_w8"], bins=bins, density=True, histtype="step", lw=1.2, ls="--", color="C2", label="cammino, stadio (coppie fuse: NON conclusivo)")
ax.plot(s, np.exp(-s), "k-", lw=0.8, label="Poisson")
ax.plot(s, np.pi * s / 2 * np.exp(-np.pi * s**2 / 4), "k--", lw=0.8, label="Wigner GOE")
ax.plot(s, 32 * s**2 / np.pi**2 * np.exp(-4 * s**2 / np.pi), "k:", lw=0.8, label="Wigner GUE")
ax.set_xlabel("distanza tra livelli vicini (media 1)"); ax.set_title("statistica dei livelli", fontsize=9); ax.legend(fontsize=7)

# row 3 right: scar
ax = fig.add_subplot(gs[2, 4:6])
ax.imshow(st["scar_img"].T, origin="lower", cmap="magma", interpolation="nearest")
ax.contour(st["scar_mask"].T.astype(float), [0.5], colors="c", linewidths=0.8)
n1, n2 = st["scar_mask"].shape
ax.plot([0, n1 - 1], [(n2 - 1) / 2] * 2, "w--", lw=0.6)
ax.set_xticks([]); ax.set_yticks([])
ax.set_title(f"cicatrice (Heller): modo #{sj['scar']['mode_index']} dello stadio,\n{sj['scar']['enhancement_caps_tube']:.1f}x di peso sull'orbita instabile (tratteggio)", fontsize=9)

# row 4 left: E1 curvature mass
ax = fig.add_subplot(gs[3, 0:3])
mu = lp["e1_mu"]
for k, c in (("ellipse", "C0"), ("trefoil", "C1"), ("stadium", "C2"), ("egg", "C4")):
    ax.loglog(mu, lp[f"e1_{k}_r2"], color=c, label=f"{k}: m2/m1")
    ax.loglog(mu, lp[f"e1_{k}_r3"], color=c, ls="--")
ax.axhline(le["targets"]["r2"], color="C3", ls="--", lw=0.8); ax.axhline(le["targets"]["r3"], color="C3", ls=":", lw=0.8)
ax.set_xlabel("mu (massa = mu x curvatura, unita' del reticolo)"); ax.set_ylabel("rapporto di massa")
ax.set_title("'massa = curvatura' sull'anello nascosto (1 parametro): nessun mu da' insieme muone (continuo) e tau (tratteggio)", fontsize=9)
ax.legend(fontsize=7)

# row 4 right: E2 localization
ax = fig.add_subplot(gs[3, 3:6])
ax.semilogy(lp["e2_Md"], lp["e2_r2"], "C0", label="m2/m1 (cammino)")
ax.semilogy(lp["e2_Md"], lp["e2_r3"], "C1", label="m3/m1 (cammino)")
ax.axhline(le["targets"]["r2"], color="C3", ls="--", lw=0.8); ax.axhline(le["targets"]["r3"], color="C3", ls=":", lw=0.8)
ax.axvline(le["E2"]["fit_Md"], color="gray", lw=0.8)
ax.annotate(f"fissato M d = {le['E2']['fit_Md']:.2f} sul muone\n-> prevede tau/e = {le['E2']['pred_r3']:.0f} (vero 3477)",
            (le["E2"]["fit_Md"], le["E2"]["pred_r3"]), xytext=(0.5, 2e4), fontsize=8, arrowprops=dict(arrowstyle="->"))
ax.set_xlabel("M d (massa x lunghezza dei segmenti)"); ax.set_title("pareti di dominio (localizzazione): gerarchia ESPONENZIALE, ma tau sbagliato x12", fontsize=9)
ax.legend(fontsize=7)

fig.suptitle("Forma nascosta: la pallina zig-zag in uno spazio nascosto con una forma", fontsize=13)
fig.savefig("out/kk_shape.png", dpi=110, bbox_inches="tight")
print("saved out/kk_shape.png")

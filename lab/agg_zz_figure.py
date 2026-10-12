"""agg_zz_figure.py - lab/out/agg_zz_figure.png: (a) 1-D dispersion of every learned level vs the fine walk,
(b) 3-D chain A: learned motion energy vs the fine walk, (c) the learned generator stencil of the a0/16 level,
(d) hydrogen levels from the learned rule, (e) the hybrid zoom (fine patch inside the coarse field),
(f) Zitterbewegung spectrum of the fine region."""
import json, os, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 137.035999084
MC2 = 510998.95
fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))

# (a) 1-D
d1 = json.load(open(os.path.join(OUT, "agg_zz_1d.json")))
ar = np.load(os.path.join(OUT, "agg_zz_1d.npz"))
th = d1["theta"]
a = ax[0, 0]
for r in d1["auto"]:
    B = r["B"]
    k, e, ef = ar[f"auto_B{B}_k"], ar[f"auto_B{B}_eps"], ar[f"auto_B{B}_eps_fine"]
    s = k * B <= np.pi / 2
    a.loglog(k[s][1:] / th, (e[s][1:]) / th, "-", lw=2, alpha=0.7, label=f"B={B} (b={r['b']}, w={r['w']})")
kk = np.logspace(-4, 1.5, 200) * th
a.loglog(kk / th, (np.sqrt(th ** 2 + kk ** 2) - th) / th, "k:", lw=1, label="Dirac sqrt(m^2+k^2)-m")
a.loglog(kk / th, kk ** 2 / (2 * th * th), "r--", lw=0.8, label="Schroedinger k^2/2m")
a.set_xlabel("k / m"); a.set_ylabel("energia di moto / m c^2"); a.set_title("(a) 1-D: dispersione dei livelli imparati")
a.legend(fontsize=7)

# (b) 3-D chain A
stA = pickle.load(open(os.path.join(OUT, "agg_zz_3d_A.pkl"), "rb"))
thA = A0 / 16 / 256
a = ax[0, 1]
for l in stA["levels"]:
    dx = l["disp"]["x"]
    kap = np.array(dx["kappa"]); el = np.array(dx["eps_learned"]); ef = np.array(dx["eps_fine"])
    k = kap / l["B"]
    a.plot(k[1:] / thA, el[1:] / thA, "-", lw=2, alpha=0.7, label=f"B={l['B']} ({l['B'] * thA:.2f} lambda_C)")
    a.plot(k[1:] / thA, ef[1:] / thA, "k:", lw=1)
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("k / m (asse x)"); a.set_ylabel("energia di moto / m c^2")
a.set_title("(b) 3-D catena A: livelli imparati (puntini: regola fine)")
a.legend(fontsize=7)

# (c) stencil of the a0/16 level: G_j along the x axis and in the z = 0 plane
L = [l for l in stA["levels"] if l["B"] == 256][0]
offs = L["offs"]; G = np.real(0.5 * (L["K"][0, 0] + L["K"][1, 1]))
w = L["w"]
img = np.zeros((2 * w + 1, 2 * w + 1))
for o, g in zip(offs, G):
    if o[2] == 0:
        img[o[0] + w, o[1] + w] = g
a = ax[0, 2]
v = np.max(np.abs(img))
im = a.imshow(img.T, origin="lower", cmap="RdBu_r", vmin=-v, vmax=v, extent=[-w - 0.5, w + 0.5, -w - 0.5, w + 0.5])
plt.colorbar(im, ax=a, fraction=0.046)
a.set_title(f"(c) stencil imparato a0/16 (B=256, w={w}): G_j, piano z=0\n(rad per passo; centro = {G[len(G) // 2]:.2e})")
a.set_xlabel("j_x"); a.set_ylabel("j_y")

# (d) hydrogen levels
a = ax[1, 0]
if os.path.exists(os.path.join(OUT, "agg_zz_atom.json")):
    at = json.load(open(os.path.join(OUT, "agg_zz_atom.json")))
    names, Es, ex = [], [], []
    for k, r in at.items():
        names.append(f"{r['kind']}\ncella a0/{A0 / r['cell_lambdaC']:.0f}")
        Es.append(r["E_eV"]); ex.append(r["schroedinger_eV"])
    x = np.arange(len(names))
    a.bar(x - 0.2, Es, 0.4, label="regola imparata")
    a.bar(x + 0.2, ex, 0.4, label="Schroedinger esatto")
    for i, (e, t) in enumerate(zip(Es, ex)):
        a.text(i, min(e, t) - 0.6, f"{100 * (e / t - 1):+.3f}%", ha="center", fontsize=8)
    a.set_xticks(x); a.set_xticklabels(names, fontsize=8)
    a.set_ylabel("energia (eV)"); a.legend(fontsize=8)
a.set_title("(d) idrogeno con la regola agglomerata (192^3 celle)")

# (e) hybrid picture and (f) ZB spectrum
if os.path.exists(os.path.join(OUT, "agg_zz_zoom.npz")):
    z = np.load(os.path.join(OUT, "agg_zz_zoom.npz"))
    zj = json.load(open(os.path.join(OUT, "agg_zz_zoom.json")))
    a = ax[1, 1]
    sc = z["slice_coarse"]
    n = sc.shape[0]; c0 = n // 2; hcell = zj["uniform_coarse"]["cell_lambdaC"]
    ext = np.array([-n / 2 - 0.5, n / 2 - 0.5]) * hcell / A0
    a.imshow(np.log10(sc.T + 1e-30), origin="lower", cmap="magma", extent=[ext[0], ext[1], ext[0], ext[1]],
             vmin=np.log10(sc.max()) - 6)
    a.add_patch(plt.Rectangle((-1.5 * hcell / A0, -1.5 * hcell / A0), 3 * hcell / A0, 3 * hcell / A0, fill=False, ec="c", lw=1))
    a.add_patch(plt.Rectangle((-0.5 * hcell / A0, -0.5 * hcell / A0), hcell / A0, hcell / A0, fill=False, ec="w", lw=1))
    a.set_xlabel("x / a0"); a.set_ylabel("y / a0")
    a.set_title("(e) zoom ibrido: 1s grossa (log |phi|^2);\nbianco = cella del protone (regola fine), ciano = alone")
    ins = a.inset_axes([0.62, 0.62, 0.36, 0.36])
    sf = z["slice_fine"].sum(axis=0)
    ins.imshow(np.log10(sf.T + 1e-30), origin="lower", cmap="magma", vmin=np.log10(sf.max()) - 3)
    ins.set_xticks([]); ins.set_yticks([])
    ins.set_title("dentro (fine, 96^2)", fontsize=7, color="w")
    a = ax[1, 2]
    om, sp = z["om"], z["sp"]
    s = om < 4
    a.semilogy(om[s], sp[s] / sp[s].max(), lw=1)
    a.axvline(zj["hybrid"]["zb_expected_2E"], color="r", ls=":", label="2 E (Zitterbewegung)")
    a.set_xlabel("omega / (m c^2 / hbar)"); a.set_ylabel("potenza (norm.)")
    a.set_title("(f) regione fine: spettro di sum 2 Re(L^dag R)")
    a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "agg_zz_figure.png"), dpi=90)
print("figure written")

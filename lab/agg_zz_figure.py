"""agg_zz_figure.py - lab/out/agg_zz_figure.png: (a) 1-D: relative error of every learned level vs the fine walk,
(b) 3-D chain A: the same, (c) the learned generator stencil of the a0/16 level, (d) hydrogen levels from the learned
rule, (e) the hybrid zoom (fine patch inside the coarse field), (f) spectrum of the fine region (A3b) and the post-hoc
Zitterbewegung test."""
import json, os, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 137.035999084
fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))

# (a) 1-D
d1 = json.load(open(os.path.join(OUT, "agg_zz_1d.json")))
ar = np.load(os.path.join(OUT, "agg_zz_1d.npz"))
th = d1["theta"]
a = ax[0, 0]
for r in d1["auto"]:
    B = r["B"]
    k, e, ef = ar[f"auto_B{B}_k"], ar[f"auto_B{B}_eps"], ar[f"auto_B{B}_eps_fine"]
    s = (k * B <= np.pi / 2) & (k > 0)
    a.loglog(k[s] / th, np.abs(e[s] / ef[s] - 1) + 1e-16, "-", lw=1.5, label=f"B={B} (b={r['b']}, w={r['w']})")
a.axvline(1 / 137.036, color="k", ls=":", lw=0.8)
a.text(1 / 137.036, 0.92, " 1/a0", fontsize=7, transform=a.get_xaxis_transform())
a.set_xlabel("k / m"); a.set_ylabel("|E_imparata / E_fine - 1|  (energia di moto)")
a.set_title("(a) 1-D: livelli imparati contro la regola fine"); a.legend(fontsize=7)

# (b) 3-D chain A
stA = pickle.load(open(os.path.join(OUT, "agg_zz_3d_A.pkl"), "rb"))
thA = A0 / 16 / 256
a = ax[0, 1]
for i, l in enumerate(stA["levels"]):
    for nm, ls in (("x", "-"), ("diag", "--")):
        dx = l["disp"][nm]
        kap = np.array(dx["kappa"]); el = np.array(dx["eps_learned"]); ef = np.array(dx["eps_fine"])
        k = kap / l["B"]
        a.loglog(k[1:] / thA, np.abs(el[1:] / ef[1:] - 1) + 1e-16, ls, color=f"C{i}", lw=1.5,
                 label=f"B={l['B']} ({l['B'] * thA:.2f} lambda_C), {nm}")
a.axvline(1 / 137.036, color="k", ls=":", lw=0.8)
a.text(1 / 137.036, 0.92, " 1/a0", fontsize=7, transform=a.get_xaxis_transform())
a.set_xlabel("k / m"); a.set_ylabel("|E_imparata / E_fine - 1|")
a.set_title("(b) 3-D catena A (Pauli): asse x (-) e diagonale (--)")
a.legend(fontsize=6)

# (c) stencil of the a0/16 level
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
a.set_title(f"(c) generatore imparato, cella a0/16 (stencil 9^3): G_j, piano z=0\n"
            f"(rad per passo; centro {G[len(G) // 2]:.3f}, vicini {img[w + 1, w]:.3f})")
a.set_xlabel("j_x"); a.set_ylabel("j_y")

# (d) hydrogen levels: relative error
a = ax[1, 0]
at = json.load(open(os.path.join(OUT, "agg_zz_atom.json")))
names, err, thr = [], [], []
for k, r in at.items():
    names.append(f"{r['kind']}  ({r['E_eV']:.4f} eV)\ncella a0/{A0 / r['cell_lambdaC']:.0f}")
    err.append(100 * r["rel_err"]); thr.append(0.1 if (r["kind"] == "1s" and r["B"] == 256) else 1.0)
x = np.arange(len(names))
a.bar(x, err, 0.5, color=["C2" if abs(e) < t else "C3" for e, t in zip(err, thr)])
for i, t in enumerate(thr):
    a.plot([i - 0.35, i + 0.35], [t, t], "k--", lw=0.8)
    a.plot([i - 0.35, i + 0.35], [-t, -t], "k--", lw=0.8)
a.axhline(0, color="k", lw=0.5)
a.set_xticks(x); a.set_xticklabels(names, fontsize=7)
a.set_ylabel("errore rispetto a Schroedinger (%)")
a.set_title("(d) idrogeno con la regola agglomerata, 192^3 celle\n(tratteggio: soglie pre-registrate)")

# (e) hybrid picture
z = np.load(os.path.join(OUT, "agg_zz_zoom.npz"))
zj = json.load(open(os.path.join(OUT, "agg_zz_zoom.json")))
a = ax[1, 1]
sc = z["slice_coarse"]
n = sc.shape[0]; hcell = zj["uniform_coarse"]["cell_lambdaC"]
ext = np.array([-n / 2 - 0.5, n / 2 - 0.5]) * hcell / A0
a.imshow(np.log10(sc.T + 1e-30), origin="lower", cmap="magma", extent=[ext[0], ext[1], ext[0], ext[1]],
         vmin=np.log10(sc.max()) - 4)
a.add_patch(plt.Rectangle((-1.5 * hcell / A0, -1.5 * hcell / A0), 3 * hcell / A0, 3 * hcell / A0, fill=False, ec="c", lw=1))
a.set_xlabel("x / a0"); a.set_ylabel("y / a0")
a.set_title("(e) zoom ibrido: 1s su celle a0/16 (log |phi|^2);\nriquadro ciano = cella del protone + alone; inserto: regione fine")
ins = a.inset_axes([0.6, 0.6, 0.38, 0.38])
sf = z["slice_fine"].sum(axis=0)
ins.imshow(np.log10(sf.T + 1e-30), origin="lower", cmap="magma", vmin=np.log10(sf.max()) - 1.5)
ins.set_xticks([]); ins.set_yticks([])
ins.set_title("96^2 celle fini (0.27 lambda_C)", fontsize=7, color="w")

# (f) spectrum of the fine region + post-hoc ZB test
a = ax[1, 2]
om, sp = z["om"], z["sp"]
s = om < 4.5
a.semilogy(om[s], sp[s] / sp[(om > 0.5) & s].max(), lw=1, label="A3: 1s ferma (interfaccia giusta)")
TH = A0 / 512
fstep = 2 * np.pi / (32 * TH)
for q in range(1, 7):
    a.axvline(q * fstep, color="0.7", lw=0.6)
a.axvline(zj["hybrid"]["zb_expected_2E"], color="r", ls=":", label="2E (Zitterbewegung)")
if os.path.exists(os.path.join(OUT, "agg_zz_zb.npz")):
    zb = np.load(os.path.join(OUT, "agg_zz_zb.npz"))
    b = zb["naive_interface_chi0"]
    bb = b - np.polyval(np.polyfit(np.arange(len(b)), b, 2), np.arange(len(b)))
    sp2 = np.abs(np.fft.rfft(bb * np.hanning(len(bb)), 8192)) ** 2
    om2 = 2 * np.pi * np.fft.rfftfreq(8192, d=TH)
    s2 = (om2 < 4.5) & (om2 > 0.3)
    a.semilogy(om2[s2], sp2[s2] / sp2[s2].max(), color="C3", lw=1, label="post hoc: + e - mescolati (chi = 0)")
a.set_ylim(1e-10, 3)
a.set_xlabel("omega / (m c^2 / hbar)   (grigio: armoniche del passo grosso)"); a.set_ylabel("potenza (norm.)")
a.set_title("(f) regione fine: spettro di sum psi^dag beta psi")
a.legend(fontsize=7, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "agg_zz_figure.png"), dpi=90)
print("figure written")

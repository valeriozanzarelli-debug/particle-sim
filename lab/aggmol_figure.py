"""aggmol_figure.py - figure of study aggmol_ -> lab/out/aggmol_figure.png"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_core as mc
import aggmol_exact as ex

O = mc.OUT
H = json.load(open(os.path.join(O, "aggmol_h2p.json")))
C = json.load(open(os.path.join(O, "aggmol_costs.json")))["rows"]
V = json.load(open(os.path.join(O, "aggmol_verdict.json")))
fig, ax = plt.subplots(2, 3, figsize=(17, 10.5))

# (a) cost law
a = ax[0, 0]
rows = [r for r in C if r["box_a0"] == 8]
fr = [r["frac"] for r in rows]
a.semilogy(fr, [r["cells"] for r in rows], "o-", label="celle (scatola ±8 a₀)")
a.semilogy(fr, [r["mem_dynamics_GB"] * 1e9 / 1e3 for r in rows], "s-", label="memoria dinamica [kB]")
a.semilogy(fr, [r["two_electron_singlet_GB_c64"] * 1e9 / 1e3 for r in rows], "^-", label="H₂ esatto 2 elettroni (6-D) [kB]")
a.axhline(16e9 / 1e3, color="k", ls="--", lw=0.8); a.text(2.1, 2.2e7, "16 GB", fontsize=8)
a.set_xscale("log", base=2); a.set_xticks(fr); a.set_xticklabels([f"a₀/{f}" for f in fr])
a.set_title("(a) legge dei costi, scatola di due atomi"); a.legend(fontsize=8); a.set_xlabel("cella")

# (b) curve
b = ax[0, 1]
Rs = [1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.0, 4.0]
Rf = np.linspace(0.9, 4.2, 80)
b.plot(Rf, [ex.energy(R)[0] + 1 / R for R in Rf], "k-", lw=1, label="esatto (due centri, nuclei fermi)")
b.plot(Rs, [H[f"curve_4_R{R}"]["E_tot_Eh"] for R in Rs], "o", color="C1", label="regola agglomerata a₀/4 (64³)")
b.axhline(-0.5, color="gray", ls=":", lw=0.8)
b.set_xlabel("R [a₀]"); b.set_ylabel("E [Eh]"); b.set_title("(b) collaudo a due centri: H₂⁺ (UN elettrone)")
b.legend(fontsize=8)
bb = b.inset_axes([0.5, 0.12, 0.45, 0.35])
bb.plot(Rs, [H[f"curve_4_R{R}"]["err_mEh"] for R in Rs], "o-", color="C1")
bb.axhline(1.594, color="r", lw=0.6); bb.axhline(-1.594, color="r", lw=0.6)
bb.set_ylabel("errore [mEh]", fontsize=7); bb.tick_params(labelsize=7); bb.set_ylim(-2, 2)

# (c), (d) states in the x-z plane
for k, (nm, st) in enumerate((("1sg", "1σg: la nube di legame"), ("2sg", "2σg: il nodo è un ellissoide confocale"))):
    c = ax[0, 2] if k == 0 else ax[1, 0]
    N = 96; h = 1024 * mc.TH / mc.A0
    v = np.load(os.path.join(O, f"aggmol_vec_{nm}_a0_4.npy")).astype(float).reshape((N,) * 3)
    sl = v[:, N // 2, :]
    x = (np.arange(N) - N // 2) * h
    m = np.abs(x) <= 5
    S = sl[np.ix_(m, m)]
    if nm == "1sg":
        c.imshow((S ** 2).T, origin="lower", extent=[x[m][0], x[m][-1], x[m][0], x[m][-1]], cmap="magma")
    else:
        lim = np.abs(S).max()
        c.imshow(S.T, origin="lower", extent=[x[m][0], x[m][-1], x[m][0], x[m][-1]], cmap="RdBu_r", vmin=-lim, vmax=lim)
        c.contour(x[m], x[m], S.T, [0], colors="k", linewidths=1.2)
        xi0 = V["C7"]["xi0_exact"]
        t = np.linspace(0, 2 * np.pi, 200)
        c.plot(np.sqrt(xi0 ** 2 - 1) * np.cos(t), xi0 * np.sin(t), "y--", lw=1.5, label=f"ellissoide esatto ξ={xi0:.4f}")
        c.plot(np.cos(t) * 1.7, np.sin(t) * 1.7, ":", color="lime", lw=1, label="una sfera (per confronto)")
        c.legend(fontsize=7, loc="lower left")
    c.plot([0, 0], [-1, 1], "w+", ms=10)
    c.set_title(f"({'c' if k == 0 else 'd'}) {st}, R = 2 a₀, a₀/4"); c.set_xlabel("x [a₀]"); c.set_ylabel("z [a₀]")

# (e) ladder
e = ax[1, 1]
lad = V["C8"]["rel_vs_masskey"]; lex = V["C8"]["rel_vs_exact"]
fs = [4, 8, 16]
e.loglog([1 / f for f in fs], [abs(lad[str(f)]) for f in fs], "o-", label="H₂⁺ R=2 contro esatto con la massa della regola")
e.loglog([1 / f for f in fs], [abs(lex[str(f)]) for f in fs], "s--", label="H₂⁺ R=2 contro esatto (m = 1)")
e.loglog([1 / f for f in fs], [abs(V["C8"]["De_diff_mEh"][str(f)]) / 102.634 for f in fs], "^-", label="|errore D_e| / D_e")
e.set_xlabel("cella [a₀]"); e.set_ylabel("errore relativo"); e.set_title("(e) scala dei livelli (a₀/16 in scatola ±6 a₀: C8 FAIL)")
e.legend(fontsize=7)

# (f) exams
f_ = ax[1, 2]; f_.axis("off")
txt = "Esami (pre-registrati):\n\n" + "\n".join(f"{k}: {'PASS' if v['PASS'] else 'FAIL'}" for k, v in V.items() if k != "summary")
txt += f"\n\n{V['summary']}"
f_.text(0.05, 0.95, txt, va="top", family="monospace", fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(O, "aggmol_figure.png"), dpi=110)
print("saved")

"""colla_figure.py - lab/out/colla_figure.png: (a) binding per route vs experiment, (b) the quantised deuteron budget,
(c) GHS lightly-bound scan, (d) shapes (projections of the baryon density)."""
import json, os, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
R = json.load(open(os.path.join(OUT, "colla_results.json")))["numbers"]
Q = json.load(open(os.path.join(OUT, "colla_quant.json")))
U = json.load(open(os.path.join(OUT, "colla_units.json"))); Eu = U["AN"]["E_unit"]
EXP = {"d": 2.2246, "3He": 7.718, "4He": 28.296}
nuc = ["d", "3He", "4He"]
routes = [("esperimento", EXP, "k"),
          ("C0 Skyrme classico\n(senza massa)", {"d": 72, "3He": 178, "4He": 310}, "#9e9e9e"),
          ("C1 + massa del pione", R["AN_classical_MeV"], "#1f77b4"),
          ("C1+C2a rotore\nquantistico", {k: R["C2a"][k]["quantised"] for k in nuc}, "#d62728")]
if "C2b" in R:
    routes.append(("C1+C2a+C2b\n(+punto zero, solo d)", {"d": R["C2b"]["d_after"]}, "#ff7f0e"))
sc = R["C3a_scan"]
if "0.95" in sc:
    routes.append(("C3a GHS alpha=0.95", {k: sc["0.95"][k]["MeV"] for k in nuc if k in sc["0.95"]}, "#2ca02c"))
st = R.get("C3a_alpha_star")
ds = R.get("C3a_direct_alpha_star")
if ds:
    routes.append((f"C3a GHS alpha*={ds['alpha']:.3f} (4He fit, corsa diretta)", {k: ds[k]["MeV"] for k in nuc if k in ds}, "#98df8a"))
elif st:
    routes.append(("C3a GHS alpha* (4He = fit)", {"d": st.get("d"), "3He": st.get("3He"), "4He": EXP["4He"]}, "#98df8a"))
if "200" in R.get("C3b", {}):
    routes.append(("C3b sestico c6=200", {k: R["C3b"]["200"][k]["MeV"] for k in ("d", "4He")}, "#9467bd"))

fig = plt.figure(figsize=(15, 10.5))
ax = fig.add_subplot(2, 2, 1)
wbar = 0.8 / len(routes)
for i, (lab, d, col) in enumerate(routes):
    xs, ys = [], []
    for j, n in enumerate(nuc):
        v = d.get(n)
        if v is not None and v > 0:
            xs.append(j + (i - len(routes) / 2 + 0.5) * wbar); ys.append(v)
    ax.bar(xs, ys, wbar, color=col, label=lab.replace("\n", " "), edgecolor="k" if lab == "esperimento" else None, lw=0.5)
ax.set_yscale("log"); ax.set_xticks(range(3)); ax.set_xticklabels(["deuterone (B=2)", "3He (B=3)", "4He (B=4)"])
ax.set_ylabel("energia di legame (MeV, scala log)"); ax.set_title("(a) Legame per strada: quale pezzo manca?")
ax.legend(fontsize=7, loc="upper left", ncol=2); ax.grid(alpha=0.3, which="both", axis="y")

# (b) deuteron budget relative to two nucleons at rest (2 M_N)
ax = fig.add_subplot(2, 2, 2)
b = Q["B"]; E1 = b["1"]["E"] * Eu; ErN = b["1"]["Erot"] * Eu; E2 = b["2"]["E"] * Eu; Erd = b["2"]["Erot"] * Eu
zp = R["C2b"]["delta"] if "C2b" in R else 0.0
ref = 2 * (E1 + ErN)
levels = [("2 nucleoni lontani\n(2 M_N = 2 x 939)", 0.0, "k"),
          ("2 nodi classici\n(senza la trottola di spin)", 2 * E1 - ref, "#7f7f7f"),
          ("ciambella classica\n(B=2)", E2 - ref, "#1f77b4"),
          ("+ rotazione J=1, I=0", E2 + Erd - ref, "#d62728"),
          ("+ punto zero vibrazioni", E2 + Erd + zp - ref, "#ff7f0e"),
          ("deuterone vero", -EXP["d"], "g")]
for i, (lab, y, c) in enumerate(levels):
    ax.hlines(y, i - 0.35, i + 0.35, color=c, lw=3)
    ax.text(i, y + 6, f"{y:+.0f}", ha="center", fontsize=9, color=c)
ax.set_xticks(range(len(levels))); ax.set_xticklabels([l[0] for l in levels], fontsize=7)
ax.axhline(0, color="k", lw=0.5, ls=":")
ax.set_ylabel("massa - 2 M_N (MeV)"); ax.set_title("(b) Il deuterone quantizzato (pioni con massa, unita' AN)")
ax.grid(alpha=0.3, axis="y")

# (c) GHS scan
ax = fig.add_subplot(2, 2, 3)
al = sorted(float(a) for a in sc)
for n, c in zip(nuc, ("#1f77b4", "#ff7f0e", "#2ca02c")):
    xs = [a for a in al if n in sc[str(a)]]
    ys = [sc[str(a)][n]["MeV"] for a in xs]
    ax.semilogy(xs, ys, "o-", color=c, label=f"{n} (modello)")
    ax.axhline(EXP[n], color=c, ls="--", lw=1)
if st:
    ax.axvline(st["alpha"], color="k", ls=":", lw=1); ax.text(st["alpha"], 1.0, " alpha*", fontsize=8)
if ds:
    for n, c in zip(nuc, ("#1f77b4", "#ff7f0e", "#2ca02c")):
        if n in ds:
            ax.plot([ds["alpha"]], [ds[n]["MeV"]], "*", ms=13, color=c, mec="k")
ax.set_xlabel("alpha (peso del potenziale GHS, 1 parametro in piu')"); ax.set_ylabel("legame (MeV, calibrazione E1 = 938 MeV)")
ax.set_title("(c) Nodi 'leggermente legati' (Gillard-Harland-Speight)"); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")

# (d) shapes
ax = fig.add_subplot(2, 2, 4); ax.axis("off"); ax.set_title("(d) Forme (densita' barionica proiettata)")
# each panel: preferred job first, then fallbacks (cloud 2026-10-13: the laptop field files are not in the repo, so the
# panels use the fields of the runs made on the cloud: B=3 at alpha = 0 and 0.95, pair/tetrahedron starts at alpha = 0.95)
cands = [[("an_B2", "B=2 pione massivo"), ("ghs0_B3", "B=3 GHS alpha=0 (tetraedro)")],
         [("an_B4", "B=4 pione massivo"), ("ghs50_B3", "B=3 GHS alpha=0.5")],
         [(R["C3a_shapes"].get("0.95_B2", {}).get("job"), "B=2 GHS 0.95"), ("ghs95_B2pair", "B=2 GHS 0.95 (da coppia)")],
         [(R["C3a_shapes"].get("0.95_B4", {}).get("job"), "B=4 GHS 0.95"), ("ghs95_B4tet", "B=4 GHS 0.95 (da tetraedro)")],
         [("ghs95_B3", "B=3 GHS 0.95")]]
jobs = []
for cl in cands:
    ok = [(j, l) for j, l in cl if j is not None and os.path.exists(os.path.join(OUT, f"colla_fields_{j}.npz"))]
    if ok:
        jobs.append(ok[0])
for i, (j, lab) in enumerate(jobs):
    bb = np.clip(np.load(os.path.join(OUT, f"colla_fields_{j}.npz"))["b"], 0, None)
    sub = fig.add_axes([0.52 + 0.093 * i, 0.12, 0.085, 0.25])
    sub.imshow(bb.sum(axis=2).T, origin="lower", cmap="magma"); sub.set_xticks([]); sub.set_yticks([])
    sub.set_title(lab.replace(" (", "\n("), fontsize=7)
    sub2 = fig.add_axes([0.52 + 0.093 * i, 0.06, 0.085, 0.06]); sub2.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "colla_figure.png"), dpi=110)
print("saved")

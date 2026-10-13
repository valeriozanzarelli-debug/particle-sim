"""agg2s_figure.py - figure of "2s all'agglomeratore": lab/out/agg2s_figure.png
(a) the box as a crystal of images: 2s energy vs box size, periodic / anti-periodic / average (development grids, exact
    Schroedinger, cells a0/4) and the learned rule K at a0/8 (L = 32 a0, 4 twists);
(b) the 2s budget term by term (K, a0/8);
(c) the nucleus cell: coarse contact density and <k^4> lost by the cell, before / after the Born tail (1s, 2s, 3s);
(d) alpha^4 shifts of the s states: aggfull style, plain contact, contact rule of the nucleus cell (T and K)."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
ld = lambda f: json.load(open(os.path.join(OUT, f)))
box, V, A = ld("agg2s_dev_box.json"), ld("agg2s_verdict.json"), ld("agg2s_atom.json")
dc = ld("agg2s_dev_contact.json")
ct = ld("agg2s_contact.json")
RY = 13.605693122994
fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))

# (a)
a = ax[0, 0]
for div, mk, ref in ((2, "s", None), (4, "o", None)):
    Ls = sorted({box[k]["L_a0"] for k in box if k.startswith(f"d{div}_")})
    r = lambda L, tw: box[f"d{div}_N{int(L * div)}_tw{tw}"]["2s"] / (-RY / 4) - 1
    E0 = np.array([r(L, "000") for L in Ls])
    E3 = np.array([r(L, "111") for L in Ls])
    if div == 4:            # reference: the contact runs (a0/4, L = 48 a0, pair average)
        Ls = Ls + [48.0]
        E0 = np.append(E0, ct["S2s_000"]["E_rel"]); E3 = np.append(E3, ct["S2s_111"]["E_rel"])
    Lr = np.array(Ls)
    refv = 0.5 * (E0[-1] + E3[-1])
    m = Lr < Lr.max()
    a.plot(Lr[m], 100 * np.abs(E0[m] - refv), mk + "-", color="C3", label=f"scatola periodica (celle a0/{div})")
    a.plot(Lr[m], 100 * np.abs(0.5 * (E0[m] + E3[m]) - refv), mk + "--", color="k", label=f"media periodico + antiperiodico (a0/{div})")
if "box_2s_K" in V["numbers"]:
    a.plot([32], [100 * V["numbers"]["box_2s_K"]], "D", color="C2", ms=9, label="regola imparata K, a0/8: periodico − media")
a.set_yscale("log")
a.set_xlabel("lato della scatola L (a0)")
a.set_ylabel("|errore della scatola| sul 2s (%)")
a.set_title("(a) la scatola = un cristallo di immagini; i bordi ritorti la tolgono")
a.legend(fontsize=7)

# (b)
b = ax[0, 1]
bud = V["numbers"].get("budget_2s_K", {})
names = [("box", "scatola\n(periodica)"), ("mass", "massa della\nregola fine θ²/3"), ("cell", "cella\n(base finita)"),
         ("alpha4_band_shape", "forma della\nbanda appresa"), ("alpha4_contact", "α⁴ contatto\n(punta)")]
vals = [100 * bud.get(k, np.nan) for k, _ in names]
b.bar(range(len(vals)), vals, color=["C3", "C1", "C0", "C4", "C5"])
b.axhline(0, color="k", lw=0.5)
b.set_xticks(range(len(vals)))
b.set_xticklabels([n for _, n in names], fontsize=8)
b.set_ylabel("contributo al 2s (% rispetto a Dirac, + = più legato)")
if bud:
    b.set_title(f"(b) il 2s pezzo per pezzo: periodico {100 * bud['observed_periodic']:+.4f} %, "
                f"con i bordi ritorti {100 * bud['observed_twisted']:+.4f} %", fontsize=9)
for i, v in enumerate(vals):
    if np.isfinite(v):
        b.text(i, v, f"{v:+.4f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=8)

# (c)
c = ax[1, 0]
hs, rows = [], []
for div in (2, 4, 8):
    k = f"d{div}_N{ {2: 32, 4: 64, 8: 128}[div] }_1s"
    if k in dc:
        hs.append(1 / div)
        rows.append((dc[k]["rho_ratio"], dc[k]["rho_born_ratio"], dc[k]["k4_ratio"], dc[k]["k4_born_ratio"]))
rows = np.array(rows)
c.plot(hs, 100 * (1 - rows[:, 0]), "o-", color="C0", label="1s: contatto perso dalla cella |ψ(0)|²")
c.plot(hs, 100 * (1 - rows[:, 2]), "s-", color="C1", label="1s: ⟨p⁴⟩ perso dalla cella")
c.plot(hs, 100 * np.abs(1 - rows[:, 1]), "o--", color="C0", label="|ψ(0)|² dopo la coda di Born (|errore|)")
c.plot(hs, 100 * np.abs(1 - rows[:, 3]), "s--", color="C1", label="⟨p⁴⟩ dopo la coda di Born (|errore|)")
for s, mk in (("1s", "^"), ("2s", "v"), ("3s", "D")):
    nm = f"BORN_a4_{s}"
    if nm in V["numbers"]:
        q = V["numbers"][nm]
        c.plot([0.25], [100 * (1 - q["rho_c_ratio"])], mk, color="C2", ms=8)
        c.plot([0.25], [100 * abs(1 - q["rho_B_ratio"])], mk, color="C2", ms=8, mfc="none", label=f"{s}, a0/4, L = 48 a0" + (" (cieco)" if s == "3s" else ""))
c.set_xscale("log"); c.set_yscale("log")
c.set_xlabel("cella h (a0)")
c.set_ylabel("%")
c.set_title("(c) la punta del nucleo: quello che la cella perde e la coda di Born")
c.legend(fontsize=7)

# (d)
d = ax[1, 1]
lab, pl, ru, ag = [], [], [], []
for X in ("T", "K"):
    for k, v in V["numbers"].get(f"F3S_{X}", {}).items():
        lab.append(f"{X}: {k}")
        pl.append(100 * v["err_plain"]); ru.append(100 * v["err_rule"]); ag.append(100 * v["err_aggfull_style"])
x = np.arange(len(lab))
d.bar(x - 0.27, ag, 0.27, color="C7", label="come aggfull (banda imparata + Darwin imparato)")
d.bar(x, pl, 0.27, color="C1", label="contatto della cella, senza coda")
d.bar(x + 0.27, ru, 0.27, color="C2", label="regola della cella del nucleo (coda di Born)")
d.axhline(5, color="k", ls=":", lw=0.8); d.axhline(-5, color="k", ls=":", lw=0.8)
d.set_xticks(x); d.set_xticklabels(lab, rotation=20, fontsize=8)
d.set_ylabel("errore dello spostamento α⁴ (%)")
d.set_yscale("symlog", linthresh=1)
d.set_title("(d) spostamenti α⁴ degli s (F3): soglia ±5 %")
d.legend(fontsize=7)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "agg2s_figure.png"), dpi=110)
print("saved")

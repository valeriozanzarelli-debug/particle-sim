"""nuclei_figure.py - lab/out/nuclei_figure.png: isosurfaces (in negative: only the knot, no balls) of the baryon density
for B = 1..4 and for the N2 aggregation runs, E/B vs literature, torus vs naive sum, binding vs real nuclei."""
import json, os, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import ndimage
import knot_core as kc
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
H = 0.2
R = json.load(open(os.path.join(OUT, "nuclei_results.json")))


def fld(job):
    p = os.path.join(OUT, f"nuclei_fields_{job}.npz")
    return np.clip(np.load(p)["b"].astype(np.float64), 0, None) if os.path.exists(p) else None


def iso(ax, b, title, col):
    T = kc.isosurface(b, 0.25 * b.max(), H, -(np.array(b.shape) - 1) / 2 * H)
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    light = np.array([0.4, -0.5, 0.75]); light /= np.linalg.norm(light)
    s = 0.35 + 0.65 * np.abs(n @ light)
    c = np.array(matplotlib.colors.to_rgb(col))
    pc = Poly3DCollection(T, facecolors=np.clip(s[:, None] * c[None], 0, 1), edgecolors="none")
    ax.add_collection3d(pc)
    lim = 2.2
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(22, 35)
    ax.set_title(title, fontsize=9)


fig = plt.figure(figsize=(15, 11))
n1 = R["N1"]["per_B"]
names = {1: "B=1 nucleone (sfera)", 2: "B=2 'deuterone' (toro)", 3: "B=3 '3He' (tetraedro)", 4: "B=4 '4He' (cubo)"}
cols = {1: "#7a7a7a", 2: "#c0504d", 3: "#4f81bd", 4: "#9bbb59"}
for i, B in enumerate((1, 2, 3, 4)):
    b = fld(f"rm{B}")
    if b is None: continue
    ax = fig.add_subplot(3, 4, i + 1, projection="3d")
    e = n1.get(str(B), {})
    iso(ax, b, f"{names[B]}\nE/12π²B = {e.get('E_over', float('nan')):.4f} (lett. {e.get('lit', float('nan')):.4f})", cols[B])
lab = {"pair25": "da 2 nodi B=1 (d=2.5)", "pair30": "da 2 nodi B=1 (d=3.0)", "tet4": "da 4 nodi a tetraedro", "sq4": "da 4 nodi a quadrato"}
k = 5
for job in ("pair25", "pair30", "tet4", "sq4"):
    b = fld(job)
    if b is None: continue
    r = R["N2"].get(job, {})
    ax = fig.add_subplot(3, 4, k, projection="3d"); k += 1
    tag = r.get("torus_test") or r.get("cube_test") or ""
    iso(ax, b, f"N2 {lab[job]}\nE/12π²B = {r.get('E_over', float('nan')):.4f}  test forma: {tag}", "#8064a2")
# E/B
ax = fig.add_subplot(3, 4, 9)
Bs = [int(B) for B in n1]; Bs.sort()
ax.plot(Bs, [n1[str(B)]["E_over"] for B in Bs], "o-", label="reticolo (h=0.2, L=48)")
ax.plot([1, 2, 3, 4], [1.2322, 1.1791, 1.1462, 1.1201], "s--", label="Battye–Sutcliffe 2002")
ax.plot([1, 2, 3, 4], [1.232, 1.208, 1.184, 1.137], "^:", label="mappe razionali (ansatz)")
ax.axhline(1, color="k", lw=0.6); ax.set_xlabel("B"); ax.set_ylabel("E / (12π² B)"); ax.legend(fontsize=7)
ax.set_title("energia per nucleone: scende con B", fontsize=9); ax.set_xticks([1, 2, 3, 4])
# torus vs naive: slice of the density in the plane of the torus
ax = fig.add_subplot(3, 4, 10)
b2 = fld("rm2"); b1 = fld("rm1")
if b2 is not None and b1 is not None:
    d = R["N3"]["naive_vs_B2"]["d_naive"]; s = d / 2 / H
    nv = ndimage.shift(b1, (s, 0, 0), order=3) + ndimage.shift(b1, (-s, 0, 0), order=3)
    m = b2.shape[0] // 2
    sl2 = 0.5 * (b2[:, :, m - 1] + b2[:, :, m]); sln = 0.5 * (nv[:, :, m - 1] + nv[:, :, m])
    ext = [-(b2.shape[0] - 1) / 2 * H, (b2.shape[0] - 1) / 2 * H] * 2
    both = np.concatenate([sl2.T, sln.T], 1)
    ax.imshow(both, origin="lower", extent=[ext[0], ext[1] + (ext[1] - ext[0]), ext[2], ext[3]], cmap="magma")
    ax.set_xlim(-3, 3 + (ext[1] - ext[0])); ax.set_ylim(-3, 3)
    ax.set_title("vero B=2 (anello)  |  somma ingenua di due B=1", fontsize=9); ax.set_yticks([]); ax.set_xticks([])
# binding vs nuclei
ax = fig.add_subplot(3, 4, 11)
bd = R["N3"].get("binding", {})
if bd:
    Bs = sorted(int(x) for x in bd)
    ax.bar(np.array(Bs) - 0.2, [bd[str(B)]["E_bind_lattice_MeV"] for B in Bs], 0.4, label="Skyrme classico (unità ANW, fit)")
    ax.bar(np.array(Bs) + 0.2, [bd[str(B)]["exp_MeV"] for B in Bs], 0.4, label="nuclei veri")
    ax.set_yscale("log"); ax.set_xticks(Bs); ax.set_xticklabels(["2H", "3He/3H", "4He"]); ax.set_ylabel("energia di legame (MeV)")
    ax.set_ylim(1, 3000); ax.legend(fontsize=7, loc="upper left"); ax.set_title("legame: troppo forte (FALLIMENTO noto)", fontsize=9)
ax = fig.add_subplot(3, 4, 12); ax.axis("off")
txt = ["Esami pre-registrati:"]
for B in ("2", "3", "4"):
    e = n1.get(B, {})
    txt.append(f"N1 B={B}: forma {e.get('N1a','-')}, E/B {e.get('N1b','-')} ({100*e.get('dev',float('nan')):+.2f}%), B {e.get('N1d','-')}")
txt.append(f"N1c E/B decrescente: {R['N1'].get('N1c','-')}")
for j in ("pair25", "pair30", "tet4", "sq4"):
    if j in R["N2"]:
        r = R["N2"][j]; txt.append(f"N2 {j}: {r.get('torus_test') or r.get('cube_test')} (E/B {r['E_over']:.4f})")
n3 = R["N3"]
txt.append(f"N3a/N3b (non è la somma): {n3.get('naive_vs_B2',{}).get('N3a','-')} / {n3.get('naive_vs_B2',{}).get('N3b','-')}")
txt.append(f"N3c legame deuterone entro x2: {n3.get('N3c','-')}")
txt.append(f"N4 topologia superfici (sfera, toro): {R['N4'].get('exam','-')}")
ax.text(0, 1, "\n".join(txt), va="top", fontsize=8.5, family="monospace")
fig.suptitle("Nuclei dai nodi: nodi B=1 del fluido di palline che si aggregano (modello di Skyrme, pioni senza massa, 0 parametri liberi)", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "nuclei_figure.png"), dpi=110)
print("saved")

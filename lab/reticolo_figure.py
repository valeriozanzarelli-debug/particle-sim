"""reticolo_figure.py - figures of the sparse hierarchical lattice study.
  lab/out/reticolo_zoom_ladder.png      the zoom-out ladder (level 0 -> 20) of the R3 world, one panel per level
  lab/out/reticolo_zoom_L00..L20.png    the same, one frame per level (for an animation)
  lab/out/reticolo_zoom_both.png        R3b: after 8 radio ticks, radio (level 10) and gamma (level 0) in one world
  lab/out/reticolo_figure.png           R2 interface curve, R2a speed, R4 scaling"""
import os, sys, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
H0_NOTE = "se il livello 20 e' un'onda radio di 1 m, h0 = 40 nm e il 'gamma' ha lambda = 0,95 um"


def unit(h):
    return "h0" if h == 1 else ("2^%d h0" % int(round(np.log2(h))))


def panel(ax, Z, k, half, rects, center, big=False):
    v = Z["v%02d" % k]
    re = np.real(v)
    m = np.nanmax(np.abs(re)) or 1.0
    ext = [center[0] - half, center[0] + half, center[1] - half, center[1] + half]
    ax.imshow(re.T, origin="lower", extent=ext, cmap="RdBu_r", vmin=-m, vmax=m, interpolation="nearest")
    for kk, col, lw in ((k + 1, "0.35", 0.5), (k, "k", 0.6), (k - 1, "#e8a000", 0.8)):
        key = "rect%02d_%02d" % (k, kk)
        if key in Z.files:
            for x0, y0, s in Z[key]:
                ax.add_patch(Rectangle((x0, y0), s, s, fill=False, ec=col, lw=lw * (1.6 if big else 1), alpha=0.7))
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("livello %d   cella = %s\nfinestra %d celle" % (k, unit(2.0 ** k), 160), fontsize=8 if not big else 11)


def ladder(tag="ladder_t0"):
    f = os.path.join(OUT, "reticolo_r3_%s.npz" % tag)
    if not os.path.exists(f):
        return
    Z = np.load(f)
    center = Z["center"]
    ks = sorted(int(n[1:]) for n in Z.files if n.startswith("v"))
    fig, axs = plt.subplots(3, 7, figsize=(17, 9.6))
    for ax, k in zip(axs.flat, ks):
        panel(ax, Z, k, float(Z["half%02d" % k]), None, center)
    fig.suptitle("Un solo mondo, 21 livelli: lo stesso campo F = E + icB (parte reale di F_z) visto da vicino (livello 0, "
                 "il 'gamma', lambda = 24 h0) a lontano (livello 20, la 'radio', lambda = 24 x 2^20 h0).\n"
                 "Quadrati neri: blocchi attivi del livello mostrato; arancio: blocchi del livello piu' fine; grigio: piu' "
                 "grosso. Ogni pannello e' 2 volte piu' lontano del precedente.", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(OUT, "reticolo_zoom_ladder.png"), dpi=110)
    plt.close(fig)
    if "wide_v" in Z.files:
        fig, ax = plt.subplots(figsize=(7.5, 7.8))
        wide(ax, Z, "La radio intera al livello 20 (lambda = 24 celle grosse = 2.5e7 h0)\nnero: blocchi del livello 20; "
                    "colori: il camino dei livelli 19, 18, 17 attorno al gamma")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "reticolo_zoom_radio_whole.png"), dpi=100)
        plt.close(fig)
    for k in ks:
        fig, ax = plt.subplots(figsize=(6.4, 6.6))
        panel(ax, Z, k, float(Z["half%02d" % k]), None, center, big=True)
        ax.text(0.01, 0.01, H0_NOTE, transform=ax.transAxes, fontsize=7, color="0.2")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "reticolo_zoom_L%02d.png" % k), dpi=90)
        plt.close(fig)


def wide(ax, Z, title):
    v = np.real(Z["wide_v"])
    c, half, k = Z["wide_c"], float(Z["wide_half"]), int(Z["wide_k"])
    m = np.nanmax(np.abs(v[np.isfinite(v)])) or 1.0
    ext = [c[0] - half, c[0] + half, c[1] - half, c[1] + half]
    ax.imshow(v.T, origin="lower", extent=ext, cmap="RdBu_r", vmin=-m, vmax=m, interpolation="nearest")
    for kk in range(k, k - 4, -1):
        key = "wide_rect_%02d" % kk
        if key in Z.files:
            col = {0: "k", 1: "#e8a000", 2: "#2a9d8f", 3: "#8e44ad"}[k - kk]
            for x0, y0, sz in Z[key]:
                ax.add_patch(Rectangle((x0, y0), sz, sz, fill=False, ec=col, lw=0.4, alpha=0.6))
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=9)


def both():
    f = os.path.join(OUT, "reticolo_r3_both_tend.npz")
    if not os.path.exists(f):
        return
    Z = np.load(f)
    center = Z["center"]
    ks = [0, 2, 4, 6, 8]
    fig, axs = plt.subplots(1, 6, figsize=(19, 3.9))
    if "wide_v" in Z.files:
        wide(axs[0], Z, "livello 10, vista larga: radio a sinistra,\ncamino del gamma a destra (blocchi 10/9/8/7)")
    for ax, k in zip(axs[1:], ks):
        panel(ax, Z, k, float(Z["half%02d" % k]), None, center)
    fig.suptitle("R3b: radio al livello 10 (lambda = 24 x 1024 h0) e gamma al livello 0 nello stesso mondo, dopo 8 tick "
                 "radio = 8192 tick gamma (vista centrata sul gamma)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(os.path.join(OUT, "reticolo_zoom_both.png"), dpi=110)
    plt.close(fig)


def summary():
    cases = [json.load(open(p)) for p in glob.glob(os.path.join(OUT, "reticolo_r2_cases", "*.json"))]
    v1 = [json.load(open(p)) for p in glob.glob(os.path.join(OUT, "reticolo_v1", "*.json"))]
    fig, axs = plt.subplots(2, 2, figsize=(13, 9))
    ax = axs[0, 0]
    for ang, ls in ((0, "-"), (30, "--")):
        for dr, col in (("f2c", "C0"), ("c2f", "C3")):
            c = sorted([r for r in cases if r["d"] == 2 and r["ang"] == ang and r["direction"] == dr], key=lambda r: r["N"])
            if not c:
                continue
            N = [r["N"] for r in c]
            ax.plot(N, [max(abs(1 - r["T"]), 1e-12) for r in c], ls, color=col, marker="o",
                    label="|1 - T| %s %d gradi" % ("fine->grosso" if dr == "f2c" else "grosso->fine", ang))
            ax.plot(N, [max(r["R"], 1e-16) for r in c], ls, color=col, marker="x", alpha=0.5,
                    label="R %s %d gradi" % ("f->g" if dr == "f2c" else "g->f", ang))
    c1 = sorted([r for r in v1 if r["ang"] == 30 and r["direction"] == "c2f"], key=lambda r: r["N"])
    if c1:
        ax.plot([r["N"] for r in c1], [abs(1 - r["T"]) for r in c1], ":", color="0.4", marker="s",
                label="|1 - T| motore v1 (30 gradi, g->f)")
    ax.axvline(16, color="k", lw=0.6)
    ax.axhline(1e-2, color="k", lw=0.6, ls=":")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_ylim(1e-16, 2)
    ax.set_xlabel("celle per lunghezza d'onda sul livello grosso (N)")
    ax.set_title("R2b: energia che non passa (|1-T|) e riflessa (R) all'interfaccia")
    ax.legend(fontsize=6.5, ncol=2)
    ax = axs[0, 1]
    for ang, ls in ((0, "-"), (30, "--")):
        for dr, col in (("f2c", "C0"), ("c2f", "C3")):
            c = sorted([r for r in cases if r["d"] == 2 and r["ang"] == ang and r["direction"] == dr and r["N"] >= 6],
                       key=lambda r: r["N"])
            if c:
                ax.plot([r["N"] for r in c], [max(1 - r["fidelity"], 1e-16) for r in c], ls, color=col, marker="o",
                        label="1 - fedelta' %s %d gradi" % (dr, ang))
                ax.plot([r["N"] for r in c], [max(r["hel_minus"] + r["hel_static"], 1e-16) for r in c], ls, color=col,
                        marker="^", alpha=0.5, label="elicita' sbagliata + modo fermo %s %d" % (dr, ang))
    ax.axhline(1e-3, color="k", lw=0.6, ls=":")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xlabel("N (celle per lambda, livello grosso)")
    ax.set_title("R2c: forma/fase ed elicita' dopo l'interfaccia")
    ax.legend(fontsize=6.5)
    ax = axs[1, 0]
    ra = []
    for p in glob.glob(os.path.join(OUT, "reticolo_r2a_*.json")):
        ra += json.load(open(p))["runs"]
    for ang, col in ((0, "C0"), (30, "C1"), (45, "C2")):
        c = sorted([r for r in ra if r["ang"] == ang], key=lambda r: r["k"])
        if c:
            ax.plot([r["k"] for r in c], [100 * r["err"] for r in c], "o-", color=col, label="%d gradi (misurata)" % ang)
            ax.axhline(100 * (c[0]["v_maxwell_same_packet"] - 1), color=col, ls=":", lw=0.8)
    ax.axhspan(-0.5, 0.5, color="0.9")
    ax.set_xlabel("livello k (cella = 2^k h0)")
    ax.set_ylabel("v/c - 1 (%)")
    ax.set_title("R2a: velocita' del pacchetto (lambda = 24 celle) su ogni livello\n(puntini: Maxwell nel continuo, stesso pacchetto)")
    ax.legend(fontsize=8)
    ax = axs[1, 1]
    try:
        r4 = json.load(open(os.path.join(OUT, "reticolo_r4.json")))
        mm = sorted(int(m) for m in r4["b_engine"])
        ax.plot(mm, [r4["b_engine"][str(m)]["cells"] for m in mm], "o", color="C3", ms=8, label="motore (gamma + radio)")
        if "b_engine_m20_from_R3" in r4:
            ax.plot([20], [r4["b_engine_m20_from_R3"]["cells"]], "o", color="C3", ms=8)
        mod = r4["b_model"]
        M = sorted(int(m) for m in mod)
        rad = r4["b_engine"][str(mm[0])]["blocks"][-1] * 256
        ax.plot(M, [mod[str(m)]["total_blocks_without_radio"] * 256 + rad for m in M], "-", color="0.3",
                label="modello a blocchi (regola di annidamento)")
        a = r4["a_one_pulse"]
        ax.plot([r["k"] for r in a], [r["cells"] for r in a], "s-", color="C0", label="un impulso solo al livello k")
        ax.set_xlabel("livello k dell'impulso radio (o distanza in livelli m)")
        ax.set_ylabel("celle attive (2-D)")
        ax.set_title("R4: il costo non dipende dalla scala assoluta;\nper unire due scale si paga un 'camino' lineare nei livelli")
        ax.legend(fontsize=8)
    except Exception as e:
        ax.text(0.1, 0.5, "R4 mancante: %s" % e)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "reticolo_figure.png"), dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    ladder()
    both()
    summary()
    print("ok")

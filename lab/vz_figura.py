"""vz_figura.py - VIDEO SUPERCRITICO: figure docs/fig/video_supercritico.png (4 scene crops + exams) and the stills
<SCR>/video_supercritico_<nome>.png (full 1080p frames)."""
import json
import os
import shutil
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vz_comune import SCR, TV, T, SC  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FERMI = [(4.4, "si_tuffa"), (30.0, "guscio"), (96.0, "positrone_esce"), (240.0, "vuoto_carico")]


def main():
    E = json.load(open(os.path.join(REPO, "lab", "out", "vz_esami.json")))
    fig = plt.figure(figsize=(16, 9), facecolor="#08090f")
    for i, (t, nome) in enumerate(FERMI):
        iv = int(np.argmin(np.abs(TV - t)))
        f = os.path.join(SCR, "frames_1920", f"{iv:05d}.png")
        im = Image.open(f).convert("RGB")
        shutil.copy(f, os.path.join(SCR, f"video_supercritico_{nome}.png"))
        ax = fig.add_axes([0.01 + 0.245 * i, 0.42, 0.24, 0.54])
        ax.imshow(np.asarray(im)[60:1060, 40:1040])
        ax.set_axis_off()
        ax.set_title(f"t = {TV[iv]:.0f} λ_C/c".replace(".", ","), color="w", fontsize=12)
    ax = fig.add_axes([0.06, 0.07, 0.40, 0.28], facecolor="#08090f")
    ax.plot(T, SC["P_in_rw"], color="#4a7bd6", lw=2, label="dentro r_w: dati")
    ax.plot(T[::3], E["serie_palline"]["P_in_rw"][::3], "o", ms=3, color="#9cc0ff", label="dentro r_w: palline")
    ax.plot(T, SC["passato_r30"], color="#c27a30", lw=2, label="oltre r = 30: dati")
    ax.plot(T[::3], E["serie_palline"]["passato_r30"][::3], "o", ms=3, color="#ffb070", label="oltre r = 30: palline")
    ax.set_xlabel("t [λ_C/c]", color="w")
    ax.tick_params(colors="w")
    ax.legend(fontsize=9, facecolor="#08090f", labelcolor="w", frameon=False)
    ax.set_title(f"V3 carica: scarto massimo {E['V3_carica']['max_P_in_rw']:.1e} (PASS)", color="w", fontsize=11)
    ax2 = fig.add_axes([0.55, 0.07, 0.40, 0.28], facecolor="#08090f")
    ov = np.array([np.nan if x is None else x for x in E["serie_sovrapposizione"]])
    ax2.plot(T, ov, color="#9be29b", lw=1.5)
    ax2.axhline(0.95, color="#ff8c80", lw=1, ls="--")
    ax2.set_ylim(0.94, 1.001)
    ax2.set_xlabel("t [λ_C/c]", color="w")
    ax2.tick_params(colors="w")
    ax2.set_title(f"V1 istogramma palline contro ρ(r): minimo {E['V1_istogramma']['sovrapposizione_min']:.4f} (soglia 0,95, PASS)",
                  color="w", fontsize=11)
    for a in (ax, ax2):
        for sp in a.spines.values():
            sp.set_color("#555")
    out = os.path.join(REPO, "docs", "fig", "video_supercritico.png")
    fig.savefig(out, dpi=110, facecolor=fig.get_facecolor())
    print(out)


if __name__ == "__main__":
    main()

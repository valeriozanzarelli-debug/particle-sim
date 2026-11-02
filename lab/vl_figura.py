"""vl_figura.py - VIDEO ATOMO: the summary figure (docs/fig/video_atomo.png) and the stills.

    python3 lab/vl_figura.py          # needs the rendered 1080p frames in <SCR>/frames_1920
"""
import json
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vl_comune import SCR, OUT, DT, timeline  # noqa: E402

QUI = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(QUI, "..", "docs", "fig")
TV = timeline()
STILLS = {"assorbe": 120.0, "oscilla": 200.0, "riemette": 360.0, "si_svuota": 1400.0}


def iv_of(t):
    return int(np.argmin(np.abs(TV - t)))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fdir = os.path.join(SCR, "frames_1920")
    os.makedirs(FIG, exist_ok=True)
    stills = {}
    for nome, t in STILLS.items():
        iv = iv_of(t)
        im = Image.open(os.path.join(fdir, f"{iv:05d}.png"))
        p = os.path.join(SCR, f"video_atomo_{nome}.png")
        im.save(p, optimize=True)
        stills[nome] = (iv, float(TV[iv]), p)
    se = json.load(open(os.path.join(OUT, "vl_serie.json")))
    ex = json.load(open(os.path.join(OUT, "vl_esami.json")))
    t = np.array(se["k"]) * DT
    fig = plt.figure(figsize=(15, 9.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 1.0], hspace=0.28, wspace=0.28)
    for c, nome in enumerate(["assorbe", "riemette", "si_svuota"]):
        iv, tt, p = stills[nome]
        im = Image.open(p)
        ax = fig.add_subplot(gs[0, c])
        ax.imshow(im.crop((150, 150, 930, 930)))
        ax.set_axis_off()
        ax.set_title(f"elettrone, t = {tt:.0f} λ_C/c ({nome.replace('_', ' ')})", fontsize=10)
    ax = fig.add_subplot(gs[0, 3])
    iv, tt, p = stills["riemette"]
    ax.imshow(Image.open(p).crop((1180, 60, 1920, 560)))
    ax.set_axis_off()
    ax.set_title(f"luce, t = {tt:.0f}: gusci della luce dell'atomo (×1000)", fontsize=10)
    a1 = fig.add_subplot(gs[1, 0])
    a1.plot(t, se["corr"], color="#2e86c1", lw=1.2, label="Pearson")
    a1.axhline(0.97, color="#2e86c1", ls=":", lw=1)
    a1.set_ylim(0.95, 1.001)
    a1.set_xlabel("t (λ_C/c)")
    a1.set_title(f"E1 palline contro |ψ|²: min {ex['E1']['corr_min']:.4f} ({'PASS' if ex['E1']['PASS'] else 'FAIL'})", fontsize=10)
    a1b = a1.twinx()
    a1b.plot(t, se["tv"], color="#c0392b", lw=1, label="TV")
    a1b.set_ylim(0, 0.1)
    a1b.set_ylabel("TV", color="#c0392b")
    a1.set_ylabel("correlazione", color="#2e86c1")
    a2 = fig.add_subplot(gs[1, 1:3])
    a2.plot(t, se["dip_manifest"], color="#7f8c8d", lw=1.2, label="<z> dal calcolo (fotogrammi.csv)")
    a2.plot(t, se["dip_palline"], ".", color="#2e86c1", ms=3, label="<z> dalle palline")
    a2.set_xlabel("t (λ_C/c)")
    a2.set_ylabel("<z> (λ_C)")
    a2.legend(fontsize=8, loc="upper right")
    a2.set_title(f"E2 dipolo: Pearson {ex['E2']['pearson']:.5f}, scarto rms {ex['E2']['rms_lambdaC']:.3f} λ_C", fontsize=10)
    a3 = fig.add_subplot(gs[1, 3])
    a3.plot(t, se["p_manifest"], color="#c0392b", lw=1.4, label="p_exc vero (manifest)")
    a3.plot(t, se["p_reticolo"], color="#7f8c8d", lw=1, label="stimatore a fase dominante, reticolo")
    a3.plot(t, se["p_palline"], ".", color="#e67e22", ms=2.5, label="stesso stimatore, dalle palline")
    a3.set_xlabel("t (λ_C/c)")
    a3.set_ylim(0, 0.2)
    a3.legend(fontsize=7, loc="lower right")
    a3.set_ylim(0, max(0.2, 1.1 * max(se["p_manifest"])))
    a3.set_title(f"E3a {'PASS' if ex['E3a']['PASS'] else 'FAIL'} (palline = reticolo), E3b {'PASS' if ex['E3b']['PASS'] else 'FAIL'} "
                 f"(max {ex['E3b']['max_diff_palline_vs_manifest']:.3f})", fontsize=10)
    fig.suptitle(f"VIDEO ATOMO: l'elettrone (Z = 60) gira, assorbe un impulso e lo rimanda come un'antenna — palline e controlli (dati {ex['impronta_esportazione']})",
                 fontsize=12)
    out = os.path.join(FIG, "video_atomo.png")
    fig.savefig(out, dpi=72, bbox_inches="tight")
    print(out, os.path.getsize(out), stills)


if __name__ == "__main__":
    main()

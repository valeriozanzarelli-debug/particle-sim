"""sc_esporta_gluoni.py - SCENE_CALCOLATE: export of a gluon film for the program (src/atomo/calcolate.rs).

    QG_SCR=<cartella dei dati qg_> python lab/sc_esporta_gluoni.py video_grande [video_piccolo video_tre]

Needs the saved frames <SCR>/qg_ym_<run>_fXXXX.npz of lab/qg_gluoni.py. Writes
    assets/scene_calcolate/gluoni_<run>/manifest.json
    assets/scene_calcolate/gluoni_<run>/palline.bin
palline.bin, little endian, frame after frame (N = manifest "n"):
    q   N x 3 uint16   position / shape * 65535 (balls are stored one after the other: x y z x y z ...)
    rb  N uint8        electric share E^2/(E^2+B^2) * 255 (colour, gauge invariant)
    lv  N uint8        255 * min(eps / eps_max, 1)^0.35, eps_max = 99.5th percentile of frame 0 (brightness 0.1 + 0.9 lv)
N = 20000: the first balls of the run (they are born in random order). Size: 8 N bytes per frame (22 MB for 140 frames).
Exam S8: Pearson r, per 8^3 cell, of the kept balls against all 160000 in every frame (the fields are not in the frames;
the fidelity of all the balls against the field energy is the "fid" of each frame, exam G7, copied into the manifest).
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_gluoni import run_cfg  # noqa: E402
from qg_run import SCR  # noqa: E402

N = 20000
CB = 8
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    for run in sys.argv[1:]:
        c = run_cfg(run)
        shape = np.array(c["shape"])
        fs = sorted(glob.glob(os.path.join(SCR, f"qg_ym_{run}_f[0-9]*.npz")))
        dst = os.path.join(REPO, "assets", "scene_calcolate", f"gluoni_{run}")
        os.makedirs(dst, exist_ok=True)
        emax = float(np.percentile(np.load(fs[0])["eb"].astype(float), 99.5))
        t, fid, r_sub = [], [], []
        ed = [np.arange(0, s + 1, CB) for s in shape]
        with open(os.path.join(dst, "palline.bin"), "wb") as o:
            for f in fs:
                z = np.load(f)
                q = z["q"]
                P = q.astype(float) / 65535 * shape[:, None]
                hf = np.histogramdd(P.T, bins=ed)[0]
                hs = np.histogramdd(P[:, :N].T, bins=ed)[0]
                r_sub.append(float(np.corrcoef(hf.ravel(), hs.ravel())[0, 1]))
                o.write(q[:, :N].T.astype("<u2").tobytes())
                o.write(np.round(np.clip(z["rb"][:N].astype(float), 0, 1) * 255).astype("u1").tobytes())
                o.write(np.round(255 * np.clip(z["eb"][:N].astype(float) / emax, 0, 1) ** 0.35).astype("u1").tobytes())
                t.append(float(z["t"]))
                fid.append(float(z["fid"]))
        man = dict(formato="particle-sim/scena-gluoni-v1", corsa=run, n=N, n_fotogrammi=len(fs), forma=shape.tolist(), ampiezza=c["amp"],
                   dt=c["dt"], t=t, G7_fedelta=fid, G7_min=min(fid), S8_sottoinsieme_r=r_sub, S8_min=min(r_sub), S8_PASS=bool(min(r_sub) >= 0.85),
                   n_palline_corsa=c["n_balls"], completa=len(fs) * c["frame_every"] >= c["steps"])
        json.dump(man, open(os.path.join(dst, "manifest.json"), "w"), indent=1)
        print(run, len(fs), "fotogrammi; G7 min", min(fid), "S8 min", min(r_sub), flush=True)


if __name__ == "__main__":
    main()

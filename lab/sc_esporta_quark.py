"""sc_esporta_quark.py - SCENE_CALCOLATE: export of the three-quark scene for the program (src/atomo/calcolate.rs).

    QG_SCR=<cartella dei dati qg_> python lab/sc_esporta_quark.py [tag=prod]

Needs <SCR>/qg_<tag>_campi.npz and qg_<tag>_analisi.json (lab/qg_run.py, lab/qg_analisi1.py). Writes
    assets/scene_calcolate/quark_Y/manifest.json
    assets/scene_calcolate/quark_Y/palline.f32     n x 5 float32 LE: x, y, z (lattice units from the centroid), quota, livello
quota = electric share of the expelled action (0..1, colour); livello = f/f_max (0..1, brightness 0.1 + 0.9 l^0.35).
Balls come from qg_render1.carica (density = f minus background where f > 2 errors), in random order. Exam S7 (balls per
lattice cell against the sampled density) goes in the manifest; A, B, C in lab/out/qg_verdetto_portatile.json.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qg_render1 as Q  # noqa: E402
from qg_run import SCR  # noqa: E402

N = 60000
D, T = 4, 4
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DST = os.path.join(REPO, "assets", "scene_calcolate", "quark_Y")


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "prod"
    Q.NB = N
    sc = Q.carica(tag, D, T)
    pos, wt = sc["pos"], sc["wt"]
    lev = np.clip((wt - 0.1) / 0.9, 0, 1) ** (1 / 0.35)
    # the hue back as the share (tinta_quota is a palette of it): recomputed here from the fields would repeat carica
    F, fe, an = Q.dati(tag, D, T)
    L = F.shape[1]
    c = np.full(3, D / 3.0)
    pts = pos + c[:, None]
    vt, vs = Q._spline(F[0], pts), Q._spline(F[1], pts)
    g = np.arange(-L / 2, L / 2, 1.0 / Q.SUB) + 0.5 / Q.SUB
    G = np.stack(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1)
    far = np.sqrt((G ** 2).sum(0)) > Q.BG_R
    gp = G + c[:, None]
    vt0 = float(np.median(Q._spline(F[0], gp)[far]))
    vs0 = float(np.median(Q._spline(F[1], gp)[far]))
    quota = np.clip((vt - vt0) / np.maximum((vt - vt0) + (vs - vs0), 1e-12), 0, 1)
    # exam S7: balls per lattice cell against the density they were drawn from
    val = Q._spline(F[2], gp)
    val = val - float(np.median(val[far]))
    err = Q.map_coordinates(fe, gp, order=1, mode="grid-wrap")
    rho = np.where(val > Q.SOGLIA_SIGMA * err, val, 0.0)
    ed = np.arange(-L / 2, L / 2 + 1, 1.0)
    H_rho = np.histogramdd(G.T, bins=[ed, ed, ed], weights=rho)[0]
    H_b = np.histogramdd(pos.T, bins=[ed, ed, ed])[0]
    r = float(np.corrcoef(H_rho.ravel(), H_b.ravel())[0, 1])
    print("S7 Pearson", r, "PASS" if r >= 0.95 else "FAIL", flush=True)
    os.makedirs(DST, exist_ok=True)
    np.stack([pos[0], pos[1], pos[2], quota, lev], 1).astype("<f4").tofile(os.path.join(DST, "palline.f32"))
    fl = an["flux"][f"{D}_{T}"]
    man = dict(formato="particle-sim/scena-quark-v1", n=N, colonne=["x", "y", "z", "quota", "livello"], unita="passi di reticolo dal centro del triangolo",
               d=D, T=T, beta=5.7, quark=(sc["quarks"].tolist()), R_Y=fl["R_Y"], placchetta=an["plaq"],
               a_sqrt_sigma=an["cornell"]["a_sqrt_sigma"], n_misure=an["n_meas"], S7_pearson=r, S7_PASS=bool(r >= 0.95))
    json.dump(man, open(os.path.join(DST, "manifest.json"), "w"), indent=1)
    print("scritto", DST, flush=True)


if __name__ == "__main__":
    main()

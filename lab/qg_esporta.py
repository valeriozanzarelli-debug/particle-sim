"""qg_esporta.py - export of the quark and gluon scenes for the single program (assets/scene/quark_gluoni/).

    python3 lab/qg_esporta.py quark                          -> quark_Y/{manifest.json, palline.f32}
    python3 lab/qg_esporta.py gluoni <run> [<run> ...]      -> gluoni_<run>/{manifest.json, palline.bin}

Same format and byte layout as the reference exporters of the scene thread (lab/sc_esporta_quark.py,
lab/sc_esporta_gluoni.py on branch claude/project-thread-sxhn28: "particle-sim/scena-quark-v1" and
"particle-sim/scena-gluoni-v1"), with two declared differences:
  * quark: quota and livello are the exact per-ball values of qg_render1.carica (field at the sampled grid point,
    background subtracted), not recomputed at the jittered position;
  * gluoni: only one saved frame out of PASSO (plus the last) is written, to keep the folder under ~25 MB; the
    kept frame indices are in the manifest ("fotogrammi_tenuti").
palline.f32 (quark): n x 5 float32 LE, row after row: x, y, z (lattice steps from the triangle centre), quota, livello.
palline.bin (gluoni), little endian, frame after frame, N = manifest "n":
    q   N x 3 uint16   position / forma * 65535, interleaved x y z x y z ...
    rb  N uint8        E^2/(E^2+B^2) * 255
    lv  N uint8        255 * min(eps / eps_max, 1)^0.35, eps_max = 99.5th percentile of eps at the balls in frame 0
"""
import glob
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_gluoni import run_cfg  # noqa: E402
from qg_run import SCR  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DST = os.path.join(REPO, "assets", "scene", "quark_gluoni")
NQ = 60000
NG = 20000
PASSO = 3
CB = 8
D, T = 4, 4
A_FM = 0.17   # lattice spacing at beta 5.7 from r0 = 0.5 fm (r0/a = 2.922): 0.171 fm (Edwards-Heller-Klassen)


def commit():
    try:
        return subprocess.check_output(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "?"


def quark():
    import qg_render1 as Q
    Q.NB = NQ
    sc = Q.carica("prod", D, T)
    pos, quota, lev = sc["pos"], sc["quota"], sc["livello"]
    dst = os.path.join(DST, "quark_Y")
    os.makedirs(dst, exist_ok=True)
    np.stack([pos[0], pos[1], pos[2], quota, lev], 1).astype("<f4").tofile(os.path.join(dst, "palline.f32"))
    an = sc["an"]
    F, fe, _ = Q.dati("prod", D, T)          # measured field on the 12^3 lattice (mid time slice)
    np.stack([F[2], fe, F[0], F[1]]).astype("<f4").tofile(os.path.join(dst, "campo.f32"))
    ver = json.load(open(os.path.join(REPO, "lab", "out", "qg_verdetto.json")))["parte1"]
    fl = an["flux"][f"{D}_{T}"]
    man = dict(formato="particle-sim/scena-quark-v1", n=NQ, colonne=["x", "y", "z", "quota", "livello"],
               dtype="<f4", righe=NQ, byte=NQ * 5 * 4,
               unita="passi di reticolo dal centro del triangolo", passo_reticolo_fm=A_FM,
               d=D, T=T, beta=5.7, L=12, quark=sc["quarks"].tolist(),
               quota="quota elettrica dell'azione espulsa dpt/(dpt+dps), 0..1 (colore)",
               livello="f/f_max, 0..1 (luminosita' 0.1 + 0.9 livello^0.35)",
               R_Y=fl["R_Y"], placchetta=an["plaq"], a_sqrt_sigma=an["cornell"]["a_sqrt_sigma"],
               n_misure=an["n_meas"], voti={k: v["voto"] for k, v in ver.items() if isinstance(v, dict)},
               stato="istantanea in tempo immaginario (Monte Carlo SU(3), quark statici): non e' un film",
               campo=dict(file="campo.f32", dtype="<f4", forma=[4, 12, 12, 12], canali=["f", "f_errore", "dpt", "dps"],
                          assi="canale, x, y, z (indice = sito; il punto di giunzione delle linee e' il sito 0, i quark "
                               "nei siti d*e_x, d*e_y, d*e_z; reticolo periodico)",
                          significato="f = (dpt + dps)/2 = meno la variazione della densita' d'azione attorno ai quark "
                                      "(tubo di flusso dove f > 0); errore jackknife; dpt, dps parti temporale e spaziale; "
                                      "senza sottrazione del fondo"),
               ordine="casuale (si possono accendere le prime k)", seme=7, commit_sorgente=commit())
    json.dump(man, open(os.path.join(dst, "manifest.json"), "w"), indent=1, ensure_ascii=False)
    print("quark", NQ, "palline ->", dst, flush=True)


def gluoni(run):
    c = run_cfg(run)
    shape = np.array(c["shape"])
    fs = sorted(glob.glob(os.path.join(SCR, f"qg_ym_{run}_f[0-9]*.npz")))
    keep = list(range(0, len(fs), PASSO))
    if keep[-1] != len(fs) - 1:
        keep.append(len(fs) - 1)
    dst = os.path.join(DST, f"gluoni_{run}")
    os.makedirs(dst, exist_ok=True)
    emax = float(np.percentile(np.load(fs[0])["eb"].astype(float), 99.5))
    t, fid, r_sub, fid_all = [], [], [], []
    ed = [np.arange(0, s + 1, CB) for s in shape]
    for f in fs:
        fid_all.append(float(np.load(f)["fid"]))
    with open(os.path.join(dst, "palline.bin"), "wb") as o:
        for i in keep:
            z = np.load(fs[i])
            q = z["q"]
            P = q.astype(float) / 65535 * shape[:, None]
            hf = np.histogramdd(P.T, bins=ed)[0]
            hs = np.histogramdd(P[:, :NG].T, bins=ed)[0]
            r_sub.append(float(np.corrcoef(hf.ravel(), hs.ravel())[0, 1]))
            o.write(q[:, :NG].T.astype("<u2").tobytes())
            o.write(np.round(np.clip(z["rb"][:NG].astype(float), 0, 1) * 255).astype("u1").tobytes())
            o.write(np.round(255 * np.clip(z["eb"][:NG].astype(float) / emax, 0, 1) ** 0.35).astype("u1").tobytes())
            t.append(float(z["t"]))
            fid.append(float(z["fid"]))
    man = dict(formato="particle-sim/scena-gluoni-v1", corsa=run, n=NG, n_fotogrammi=len(keep),
               byte_per_fotogramma=8 * NG, forma=shape.tolist(), ampiezza=c["amp"], dt=c["dt"],
               passo_fotogrammi=PASSO, fotogrammi_tenuti=keep, n_fotogrammi_corsa=len(fs), t=t,
               eps_max=emax, G7_fedelta=fid, G7_min=min(fid_all), G7_min_tenuti=min(fid),
               S8_sottoinsieme_r=r_sub, S8_min=min(r_sub), S8_PASS=bool(min(r_sub) >= 0.85),
               n_palline_corsa=c["n_balls"], completa=len(fs) * c["frame_every"] >= c["steps"],
               unita="passi di reticolo (scatola periodica), tempo in unita' di reticolo (c = 1)",
               stato="Yang-Mills SU(3) classico in tempo reale: niente quantistica, niente confinamento",
               pre_registrata=run in ("video_piccolo", "video_grande"), commit_sorgente=commit())
    json.dump(man, open(os.path.join(dst, "manifest.json"), "w"), indent=1, ensure_ascii=False)
    print(run, len(keep), "/", len(fs), "fotogrammi; G7 min", min(fid_all), "S8 min", min(r_sub), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "quark":
        quark()
    else:
        for r in sys.argv[2:]:
            gluoni(r)

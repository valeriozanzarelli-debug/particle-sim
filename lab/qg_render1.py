"""qg_render1.py - PART 1 picture: three static quarks and the flux measured around them, as glows only.

    python3 lab/qg_render1.py still <tag> <d> <T> <view> out.png [--scala 1]   view = alto | obliquo | <angle deg>
    python3 lab/qg_render1.py orbita <tag> <d> <T> <nframes> <i0> <i1> [--scala 1]  -> <SCR>/qg_orbita_<W>/

Balls are sampled ONCE with density proportional to max(f, 0), f = minus the measured change of the action
density (lab/qg_analisi1.py), interpolated from the 12^3 lattice with a periodic cubic spline (declared).
The quarks are three bright glows at the input positions. Nothing moves: it is a snapshot; only the camera turns.
"""
import json
import os
import sys
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import map_coordinates

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_bagliori import bagliori, camera, project  # noqa: E402
from qg_run import SCR  # noqa: E402

LAB = os.path.dirname(os.path.abspath(__file__))
W0, H0 = 1920, 1080
SCALA = 1.0
FONT = "/usr/share/fonts/opentype/inter/Inter-Regular.otf"
FONTB = "/usr/share/fonts/opentype/inter/Inter-SemiBold.otf"
# aesthetics (knobs, counted in the doc)
NB = 220000
SUB = 6              # spline sub-samples per lattice spacing
SPOT = 2.0
GLOW = 16.0
GAIN = 2.0
DEPTH = 18.0
CAM_D = 26.0
FOV = 38.0
SOGLIA_SIGMA = 1.0   # balls only where f exceeds this many jackknife errors (declared: noise is not drawn)
C_LOW = np.array([0.85, 0.20, 0.06])
C_MID = np.array([1.00, 0.62, 0.16])
C_HIGH = np.array([1.00, 0.95, 0.80])
C_QUARK = np.array([0.70, 0.85, 1.00])


def rot_to_z(n):
    n = n / np.linalg.norm(n)
    z = np.array([0, 0, 1.0])
    v = np.cross(n, z)
    s, c = np.linalg.norm(v), n @ z
    if s < 1e-12:
        return np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * ((1 - c) / s ** 2)


def dati(tag, d, T):
    z = np.load(os.path.join(SCR, f"qg_{tag}_campi.npz"))
    F = z[f"F_{d}_{T}"]
    Fe = z[f"Ferr_{d}_{T}"]
    an = json.load(open(os.path.join(SCR, f"qg_{tag}_analisi.json")))
    return F[2], Fe[2], an


def palline(f, fe, d, rng):
    L = f.shape[0]
    c = np.full(3, d / 3.0)
    g = np.arange(-L / 2, L / 2, 1.0 / SUB) + 0.5 / SUB
    G = np.stack(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1)
    pts = G + c[:, None]
    from scipy.ndimage import spline_filter
    val = map_coordinates(spline_filter(f, order=3, mode="grid-wrap"), pts, order=3, mode="grid-wrap", prefilter=False)
    err = map_coordinates(fe, pts, order=1, mode="grid-wrap")
    rho = np.where(val > SOGLIA_SIGMA * err, val, 0.0)
    rho = np.maximum(rho, 0)
    p = rho / rho.sum()
    idx = rng.choice(rho.size, NB, p=p)
    pos = G[:, idx] + (rng.random((3, NB)) - 0.5) / SUB      # relative to the centroid, lattice units
    lev = val[idx] / val.max()
    return pos, np.clip(lev, 0, 1)


def colori(lev):
    a = np.clip(lev / 0.5, 0, 1)[:, None]
    b = np.clip((lev - 0.5) / 0.5, 0, 1)[:, None]
    return (1 - a) * C_LOW + a * (1 - b) * C_MID + a * b * C_HIGH


def font(sz, bold=False):
    return ImageFont.truetype(FONTB if bold else FONT, max(8, int(round(sz * SCALA))))


def num(x, f="{:.3f}"):
    return f.format(x).replace(".", ",")


def scena(pos, lev, d, view, w, h):
    R = rot_to_z(np.ones(3))
    q = np.array([[d, 0, 0], [0, d, 0], [0, 0, d]], float) - d / 3.0
    P = R @ pos
    Q = R @ q.T
    if view == "alto":
        el, az = np.radians(89.0), np.radians(-90)
    elif view == "obliquo":
        el, az = np.radians(28.0), np.radians(-60)
    else:
        el, az = np.radians(26.0 + 6 * np.sin(np.radians(float(view)))), np.radians(float(view))
    eye = CAM_D * np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    up = (0, 1, 0) if view == "alto" else (0, 0, 1)
    cam = camera(eye, (0, 0, 0), FOV, w, h, up=up)
    s = w / 1080
    img, _ = bagliori(P, colori(lev), cam, w, h, SPOT * s, GLOW * s, GAIN, DEPTH, CAM_D,
                      extra=[(Q, np.tile(C_QUARK, (3, 1)), 5.0 * s, 1.0)])
    return img, cam, Q


def frame(tag, d, T, view, an):
    s = SCALA
    W, H = int(W0 * s), int(H0 * s)
    f, fe, _ = dati(tag, d, T)
    rng = np.random.default_rng(7)
    pos, lev = palline(f, fe, d, rng)
    bg = np.array([0.010, 0.011, 0.020])
    canvas = np.zeros((H, W, 3))
    canvas[:] = bg
    wa = H
    sc, cam, Q = scena(pos, lev, d, view, wa, wa)
    sc = np.clip(sc + bg[:, None, None] * (1 - sc), 0, 1) ** (1 / 1.1)
    canvas[:, :wa] = sc.transpose(1, 2, 0)
    im = Image.fromarray((canvas * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im)
    dr.text((22 * s, 18 * s), "tre quark fermi (bagliori azzurri) e il campo misurato attorno a loro", fill=(205, 208, 230), font=font(15))
    dr.text((22 * s, 40 * s), "istantanea in tempo immaginario, non un film; il confinamento qui è misurato, non disegnato",
            fill=(255, 205, 120), font=font(13))
    pannello(dr, wa, W - wa, d, T, an)
    return im


def pannello(dr, x0, w, d, T, an):
    s = SCALA

    def Tx(y, txt, sz=13, col=(215, 218, 235), bold=False):
        dr.text((x0 + 24 * s, y * s), txt, fill=col, font=font(sz, bold))

    a = 0.40 / 1.0
    Tx(22, "Dentro un barione: tre quark e la colla", 24, (245, 245, 252), True)
    Tx(56, "SU(3) su reticolo 12⁴, β = 6/g² = 5,7 (l'unico numero della regola)", 14, (170, 175, 200))
    fl = an["flux"][f"{d}_{T}"]
    cor = an["cornell"]
    pl = an["plaq"]
    ry = fl["R_Y"]
    lato = d * np.sqrt(2)
    asig = cor["a_sqrt_sigma"]
    afm = 0.44 / max(asig[0], 1e-9) if False else None
    lines = [
        ("Cosa vedi", (245, 245, 252), True),
        (f"• i tre bagliori azzurri: tre quark fermi (cariche di colore statiche), ai vertici di un", (205, 208, 230), False),
        (f"  triangolo equilatero di lato {num(lato, '{:.2f}')} passi di reticolo (circa {num(lato * 0.17, '{:.1f}')} fm)", (205, 208, 230), False),
        ("• le palline arancio: dove il campo dei gluoni è diverso dal vuoto perché ci sono i quark", (205, 208, 230), False),
        ("  (fluttuazioni del vuoto espulse = il tubo di flusso); più chiaro = più intenso", (205, 208, 230), False),
        ("• palline solo dove il segnale supera 1 errore: il rumore non si disegna", (170, 175, 200), False),
        ("", None, False),
        ("I numeri (esami pre-registrati)", (245, 245, 252), True),
        (f"placchetta {num(pl[0], '{:.4f}')} ± {num(pl[1], '{:.4f}')} (letteratura 0,5492)", (205, 208, 230), False),
        (f"tensione di corda a√σ = {num(asig[0])} ± {num(asig[1])} (letteratura ≈ 0,40)", (205, 208, 230), False),
        (f"Y o Δ?  R_Y = f(centro)/f(metà dei lati) = {num(ry[0], '{:.2f}')} ± {num(ry[1], '{:.2f}')}", (205, 208, 230), False),
    ]
    ver = an.get("verdetti", {})
    for k in ("A_placchetta", "B_tensione", "C_Y_contro_Delta"):
        if k in ver:
            v = ver[k]
            lines.append((f"{k.split('_')[0]}: {v['voto']}  ({v.get('nota', '')})", (170, 230, 170) if v["voto"].startswith("PASS") else (255, 150, 140), False))
    lines += [
        ("", None, False),
        ("Cosa NON è", (245, 245, 252), True),
        ("• non è un film: è una media su centinaia di 'fotografie' del vuoto in tempo immaginario", (255, 205, 120), False),
        ("• i quark sono fermi e infinitamente pesanti (sorgenti statiche): niente quark che si muovono", (255, 205, 120), False),
        ("• il colore delle palline è solo l'intensità: il 'colore' dei quark non si può disegnare", (255, 205, 120), False),
        ("  (non è una cosa gauge-invariante)", (255, 205, 120), False),
        ("• quark che si muovono e restano chiusi in tempo reale: nessuno sa calcolarlo (problema del segno)", (255, 205, 120), False),
        ("", None, False),
        ("reticolo grossolano (a ≈ 0,17 fm): 12 punti per lato, forme a grana grossa", (150, 155, 185), False),
        ("posizioni dei quark = ingresso; il tubo = misurato (spline cubica fra i punti del reticolo)", (150, 155, 185), False),
    ]
    y = 100
    for txt, col, b in lines:
        if col is not None:
            Tx(y, txt, 15 if b else 12.6, col, b)
        y += 26 if b else 19


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--scala" in args:
        i = args.index("--scala")
        SCALA = float(args[i + 1])
        del args[i:i + 2]
    tag, d, T = args[1], int(args[2]), int(args[3])
    an = json.load(open(os.path.join(SCR, f"qg_{tag}_analisi.json")))
    vp = os.path.join(LAB, "out", "qg_verdetto.json")
    if os.path.exists(vp):
        an["verdetti"] = json.load(open(vp)).get("parte1", {})
    if args[0] == "still":
        t0 = time.time()
        frame(tag, d, T, args[4], an).save(args[5])
        print("still", args[5], f"{time.time() - t0:.1f}s")
    elif args[0] == "orbita":
        n, a, b = int(args[4]), int(args[5]), int(args[6])
        dd = os.path.join(SCR, f"qg_orbita_{int(W0 * SCALA)}")
        os.makedirs(dd, exist_ok=True)
        t0 = time.time()
        for i in range(a, min(b, n)):
            fn = os.path.join(dd, f"{i:05d}.png")
            if os.path.exists(fn):
                continue
            frame(tag, d, T, str(-90 + 360.0 * i / n), an).save(fn + ".tmp.png")
            os.replace(fn + ".tmp.png", fn)
            if time.time() - t0 > 480:
                break
        print("orbita", a, b, f"{time.time() - t0:.0f}s")

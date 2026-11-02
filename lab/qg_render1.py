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
from qg_bagliori import ZZ, bagliori, camera, project, tinta_quota  # noqa: E402
from qg_run import SCR  # noqa: E402

LAB = os.path.dirname(os.path.abspath(__file__))
W0, H0 = 1920, 1080
SCALA = 1.0
FONT = "/usr/share/fonts/opentype/inter/Inter-Regular.otf"
FONTB = "/usr/share/fonts/opentype/inter/Inter-SemiBold.otf"
# aesthetics (knobs, counted in the doc)
NB = 220000
SUB = 6              # spline sub-samples per lattice spacing
SPOT = 1.5
GLOW = 11.0
GAIN = 2.8
DEPTH = 18.0
CAM_D = 26.0
FOV = 38.0
BG_R = 5.0           # lattice units: background region (far from the triangle)
SOGLIA_SIGMA = 2.0   # balls only where f exceeds this many jackknife errors (declared: noise is not drawn)
C_LOW = np.array([0.85, 0.20, 0.06])
C_MID = np.array([1.00, 0.62, 0.16])
C_HIGH = np.array([1.00, 0.95, 0.80])
C_QUARK = np.array([1.00, 1.00, 1.00])


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
    return F, Fe[2], an


def _spline(field, pts):
    from scipy.ndimage import spline_filter
    return map_coordinates(spline_filter(field, order=3, mode="grid-wrap"), pts, order=3, mode="grid-wrap", prefilter=False)


def carica(tag, d, T, seed=7):
    """LOADER (portable: the per-ball data of the scene). Returns dict with
    pos (3,N) lattice units relative to the centroid, col (N,3), wt (N,), quarks (3,3), plus the numbers.
    Balls: density proportional to f - background where f > SOGLIA_SIGMA errors; brightness weight as zzscena.rs
    (0.1 + 0.9 (f/f_max)^0.35); hue = electric share of the expelled action, dpt / (dpt + dps), on the zig-zag
    palette (cyan = magnetic, violet = equal, gold = electric): gauge invariant, measured."""
    F, fe, an = dati(tag, d, T)
    rng = np.random.default_rng(seed)
    L = F.shape[1]
    c = np.full(3, d / 3.0)
    g = np.arange(-L / 2, L / 2, 1.0 / SUB) + 0.5 / SUB
    G = np.stack(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1)
    pts = G + c[:, None]
    far = np.sqrt((G ** 2).sum(0)) > BG_R
    vals = []
    for k in (2, 0, 1):
        v = _spline(F[k], pts)
        # uniform background: on a 2 fm periodic box f does not go back to 0 far away (the whole box feels the
        # quarks + the sum rule of the action). We subtract the median farther than BG_R from the centroid
        # (declared display choice; the exams use f as it is).
        vals.append(v - float(np.median(v[far])))
    val, vt, vs = vals
    err = map_coordinates(fe, pts, order=1, mode="grid-wrap")
    rho = np.where(val > SOGLIA_SIGMA * err, val, 0.0)
    idx = rng.choice(rho.size, NB, p=rho / rho.sum())
    pos = G[:, idx] + (rng.random((3, NB)) - 0.5) / SUB
    lev = np.clip(val[idx] / val.max(), 0, 1)
    quota = np.clip(vt[idx] / np.maximum(vt[idx] + vs[idx], 1e-12), 0, 1)
    q = np.array([[d, 0, 0], [0, d, 0], [0, 0, d]], float) - d / 3.0
    return {"pos": pos, "col": tinta_quota(quota), "wt": 0.1 + 0.9 * lev ** 0.35, "quarks": q, "an": an,
            "quota_media": float(np.average(quota, weights=lev))}


def disegna(sc, view, w, h):
    """DRAW (zig-zag base look, qg_bagliori.ZZ): returns image (3,h,w) in 0..1 with the background."""
    R = rot_to_z(np.ones(3))
    P = R @ sc["pos"]
    Q = R @ sc["quarks"].T
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
    img, _ = bagliori(P, sc["col"], cam, w, h, ZZ["spot"] * s, ZZ["glow"] * s, ZZ["gain"], DEPTH, CAM_D,
                      weight=sc["wt"], keep=ZZ["keep"], glow_mix=ZZ["glow_mix"], gamma=ZZ["gamma"],
                      extra=[(Q, np.tile(C_QUARK, (3, 1)), 6.0 * s, 1.6)])
    bg = np.array(ZZ["bg"])[:, None, None]
    return np.clip(img + bg * (1 - img), 0, 1) ** (1 / ZZ["out_gamma"]), cam


def font(sz, bold=False):
    return ImageFont.truetype(FONTB if bold else FONT, max(8, int(round(sz * SCALA))))


def num(x, f="{:.3f}"):
    return f.format(x).replace(".", ",")


_CACHE = {}


def frame(tag, d, T, view, an):
    s = SCALA
    W, H = int(W0 * s), int(H0 * s)
    key = (tag, d, T)
    if key not in _CACHE:
        _CACHE[key] = carica(tag, d, T)
    sc = _CACHE[key]
    canvas = np.zeros((H, W, 3))
    canvas[:] = np.array(ZZ["bg"])
    wa = H
    img, cam = disegna(sc, view, wa, wa)
    canvas[:, :wa] = img.transpose(1, 2, 0)
    im = Image.fromarray((canvas * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im)
    dr.text((22 * s, 18 * s), "tre quark fermi (tre grandi bagliori bianchi) e il campo misurato attorno a loro", fill=(205, 208, 230), font=font(15))
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
    et = os.environ.get("QG_ETICHETTA", "")
    if et:
        Tx(78, et, 14, (255, 150, 60), True)
    fl = an["flux"][f"{d}_{T}"]
    cor = an["cornell"]
    pl = an["plaq"]
    ry = fl["R_Y"]
    lato = d * np.sqrt(2)
    asig = cor["a_sqrt_sigma"]
    afm = 0.44 / max(asig[0], 1e-9) if False else None
    lines = [
        ("Cosa vedi", (245, 245, 252), True),
        (f"• i tre grandi bagliori bianchi: tre quark fermi (cariche di colore statiche), ai vertici di un", (205, 208, 230), False),
        (f"  triangolo equilatero di lato {num(lato, '{:.2f}')} passi di reticolo (circa {num(lato * 0.17, '{:.1f}')} fm)", (205, 208, 230), False),
        ("• le palline: dove il campo dei gluoni è diverso dal vuoto perché ci sono i quark", (205, 208, 230), False),
        ("  (fluttuazioni del vuoto espulse = il tubo di flusso); più luminose = più intenso", (205, 208, 230), False),
        ("• colore (tavolozza dell'atomo zig-zag): quota elettrica del tubo, ciano = magnetico,", (205, 208, 230), False),
        ("  viola = metà e metà, oro = elettrico (misurata, non dipende dalla gauge)", (205, 208, 230), False),
        ("• palline solo dove il segnale supera 2 errori: il rumore non si disegna; tolto il fondo", (170, 175, 200), False),
        ("  uniforme (f lontano dai quark: la scatola di 2 fm è piccola)", (170, 175, 200), False),
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
        ("• il 'colore' (carica di colore) dei quark non si può disegnare", (255, 205, 120), False),
        ("  (non è una cosa gauge-invariante)", (255, 205, 120), False),
        ("• quark che si muovono e restano chiusi in tempo reale: nessuno sa calcolarlo (problema del segno)", (255, 205, 120), False),
        ("", None, False),
        ("reticolo grossolano (a ≈ 0,17 fm): 12 punti per lato, forme a grana grossa", (150, 155, 185), False),
        ("posizioni dei quark = ingresso; il tubo = misurato (spline cubica fra i punti del reticolo)", (150, 155, 185), False),
    ]
    y = 104
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
            fn = os.path.join(dd, f"{i:05d}.jpg")
            if os.path.exists(fn):
                continue
            frame(tag, d, T, str(-90 + 360.0 * i / n), an).save(fn + ".tmp.jpg", quality=93)
            os.replace(fn + ".tmp.jpg", fn)
            if time.time() - t0 > 480:
                break
        print("orbita", a, b, f"{time.time() - t0:.0f}s")

"""qg_render2.py - PART 2 film: two gluon packets meet (classical SU(3) Yang-Mills), glows only.

    python3 lab/qg_render2.py still <run> <k> out.png [--scala 1]
    python3 lab/qg_render2.py batch <i0> <i1> [--scala 1]        -> <SCR>/qg_film_<W>/%05d.png (whole film list)

Every ball is a fixed share of the field energy, born with density proportional to the energy and carried by
v = S/eps (Yang-Mills Poynting vector over energy density) inside lab/qg_gluoni.py; positions saved every 4 steps
(dt 0.4) and linearly interpolated between saved frames for the slow motion (SUBF sub-frames, declared).
Colour = electric fraction E^2/(E^2+B^2) at the ball (gauge invariant). No lines, no trails.
"""
import glob
import json
import os
import sys
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_bagliori import ZZ, bagliori, camera, tinta_quota  # noqa: E402
from qg_run import SCR  # noqa: E402

LAB = os.path.dirname(os.path.abspath(__file__))
W0, H0 = 1920, 1080
SCALA = 1.0
FONT = "/usr/share/fonts/opentype/inter/Inter-Regular.otf"
FONTB = "/usr/share/fonts/opentype/inter/Inter-SemiBold.otf"
SHAPE = np.array([96, 48, 48])
WSC = 1 / 7.0          # world units per lattice cell
KEEP = 30000          # sharp glows (a fair sample: balls are born in random order); all go in the halo
SUBF = 3               # rendered frames per saved frame
HOLD = 36              # frames held at the start of each run
SPOT = 1.0
GLOW = 5.0
GAIN = 1.6
DEPTH = 22.0
CAM_D = 14.0
FOV = 40.0
C_B = np.array([0.25, 0.55, 1.00])    # magnetic dominated
C_0 = np.array([1.00, 0.78, 0.40])    # balanced (a free travelling wave has E^2 = B^2)
C_E = np.array([1.00, 0.28, 0.50])    # electric dominated
RUNS = ["video_piccolo", "video_grande", "video_tre"]
# short fading glow tail (rule of 13/10 evening): the same glow re-splatted at the positions of the last
# CODA rendered frames (0.25 s at 24 fps), dimmer and smaller; no lines, no chains of separate balls
CODA = 6
CODA_DECAY = 2.2      # weight exp(-k / CODA_DECAY)
CODA_SHRINK = 0.10    # radius x (1 - CODA_SHRINK k)
CODA_FRAC = 2         # one ball out of CODA_FRAC carries a tail (cost)
CODA_TOT = 0.6        # whole tail = 0.6 of the head's light
NORM = sum(np.exp(-k / CODA_DECAY) for k in range(1, CODA + 1))
REFS = {}
EMAX = {}


def files(run):
    return sorted(glob.glob(os.path.join(SCR, f"qg_ym_{run}_f[0-9]*.npz")))


def set_shape(run):
    global SHAPE
    from qg_gluoni import run_cfg
    SHAPE = np.array(run_cfg(run)["shape"])


def load(f):
    z = np.load(f)
    pos = z["q"].astype(float) / 65535 * SHAPE[:, None]
    return pos, z["rb"].astype(float), float(z["t"]), z["eb"].astype(float)


def lista(runs=None):
    """global frame list: (run, saved index a, fraction)"""
    out = []
    for run in (runs or RUNS):
        n = len(files(run))
        out += [(run, 0, 0.0)] * HOLD
        for a in range(n - 1):
            for s in range(SUBF):
                out.append((run, a, s / SUBF))
        out += [(run, n - 1, 0.0)] * HOLD
    return out


def stato(run, a, fr):
    set_shape(run)
    fs = files(run)
    p0, r0, t0, e0 = load(fs[a])
    if fr == 0 or a + 1 >= len(fs):
        return p0, r0, t0, e0
    p1, r1, t1, e1 = load(fs[a + 1])
    dp = p1 - p0
    dp -= SHAPE[:, None] * np.round(dp / SHAPE[:, None])      # periodic: shortest way
    return (p0 + fr * dp) % SHAPE[:, None], (1 - fr) * r0 + fr * r1, (1 - fr) * t0 + fr * t1, (1 - fr) * e0 + fr * e1


def colori(rb):
    """zig-zag palette (zzscena.rs 'tinta') on the electric share E^2/(E^2+B^2), stretched x1.5 around 1/2:
    cyan = magnetic, violet = balanced (free wave), gold = electric (declared)."""
    return tinta_quota(0.5 + 1.5 * (rb - 0.5))


def cam_at(u, w, h):
    az = np.radians(-70 + 50 * u)
    el = np.radians(24 + 6 * np.sin(np.pi * u))
    eye = CAM_D * np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    return camera(eye, (0, 0, 0), FOV, w, h)


def font(sz, bold=False):
    return ImageFont.truetype(FONTB if bold else FONT, max(8, int(round(sz * SCALA))))


def num(x, f="{:.2f}"):
    return f.format(x).replace(".", ",")


def indietro(run, a, fr, k):
    """saved-frame coordinates k rendered frames earlier (clamped at the start of the run)."""
    x = a + fr - k / SUBF
    if x <= 0:
        return 0, 0.0
    ia = int(np.floor(x + 1e-9))
    return ia, x - ia


def carica(run, a, fr):
    """loader: positions (world), colours, weights, sizes of all glows of this frame (tail included)."""
    pos, rb, t, eb = stato(run, a, fr)
    if run not in EMAX:
        EMAX[run] = float(np.percentile(stato(run, 0, 0.0)[3], 99.5))
    lum = 0.1 + 0.9 * np.clip(eb / EMAX[run], 0, 1) ** 0.35        # zzscena.rs: alpha ~ density^0.35
    P = [(pos - SHAPE[:, None] / 2) * WSC]
    C = [colori(rb)]
    Wt = [lum]
    Sz = [np.ones(pos.shape[1])]
    sel = slice(0, None, CODA_FRAC)
    for k in range(1, CODA + 1):
        ia, f2 = indietro(run, a, fr, k)
        if (ia, f2) == (a, fr):
            break
        pk, rk, _, ek = stato(run, ia, f2)
        dp = pk[:, sel] - pos[:, sel]
        dp -= SHAPE[:, None] * np.round(dp / SHAPE[:, None])
        P.append((pos[:, sel] + dp - SHAPE[:, None] / 2) * WSC)
        C.append(colori(rk[sel]))
        n = P[-1].shape[1]
        Wt.append(lum[sel] * CODA_FRAC * CODA_TOT * np.exp(-k / CODA_DECAY) / NORM)
        Sz.append(np.full(n, max(1 - CODA_SHRINK * k, 0.3)))
    return np.concatenate(P, 1), np.concatenate(C), np.concatenate(Wt), np.concatenate(Sz), pos, t


def frame(run, a, fr, u, ver):
    s = SCALA
    W, H = int(W0 * s), int(H0 * s)
    P, col, wt, sz, pos, t = carica(run, a, fr)
    nb = pos.shape[1]
    sharp = np.zeros(P.shape[1], bool)
    KEEP = ZZ["keep"]
    sharp[:min(KEEP, nb)] = True
    off = nb
    while off < P.shape[1]:
        m = min(-(-nb // CODA_FRAC), P.shape[1] - off)
        sharp[off:off + min(-(-KEEP // CODA_FRAC), m)] = True
        off += m
    order = np.r_[np.nonzero(sharp)[0], np.nonzero(~sharp)[0]]
    P, col, wt, sz = P[:, order], col[order], wt[order], sz[order]
    nsharp = int(sharp.sum())
    wa = H
    cam = cam_at(u, wa, wa)
    key = run
    img, refs = bagliori(P, col, cam, wa, wa, ZZ["spot"] * wa / 1080, ZZ["glow"] * wa / 1080, ZZ["gain"], DEPTH, CAM_D,
                         refs=REFS.get(key), weight=wt, size=sz, keep=nsharp, glow_mix=ZZ["glow_mix"], gamma=ZZ["gamma"])
    REFS.setdefault(key, refs)            # exposure fixed by the first frame of each run (declared)
    bg = np.array(ZZ["bg"])
    canvas = np.zeros((H, W, 3))
    canvas[:] = bg
    sc = np.clip(img + bg[:, None, None] * (1 - img), 0, 1) ** (1 / ZZ["out_gamma"])
    canvas[:, :wa] = sc.transpose(1, 2, 0)
    im = Image.fromarray((canvas * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im)
    small = run == "video_piccolo"
    tit = {"video_piccolo": "ampiezza piccola: i due pacchetti si attraversano",
           "video_grande": "ampiezza grande: i due pacchetti si urtano",
           "video_tre": "tre gluoni di colori diversi da tre direzioni si incontrano"}[run]
    dr.text((22 * s, 18 * s), tit, fill=(225, 228, 245), font=font(17, True))
    dr.text((22 * s, 44 * s), "campo classico: niente gluoni quantistici, niente confinamento in tempo reale", fill=(255, 205, 120), font=font(13))
    pannello(dr, wa, W - wa, run, t, pos, ver)
    return im


def pannello(dr, x0, w, run, t, pos, ver):
    s = SCALA

    def Tx(y, txt, sz=13, col=(215, 218, 235), bold=False):
        dr.text((x0 + 24 * s, y * s), txt, fill=col, font=font(sz, bold))

    Tx(22, "Gluoni che si urtano", 24, (245, 245, 252), True)
    Tx(56, f"Yang-Mills SU(3) classico su reticolo {SHAPE[0]} × {SHAPE[1]} × {SHAPE[2]}, tempo reale", 14, (170, 175, 200))
    Tx(80, f"t = {num(t, '{:.1f}')} (in passi di reticolo / c) · rallentato", 15, (235, 235, 245))
    # energy profile along x from the balls (UI chart: dots only)
    cx, cw = x0 + 60 * s, w - 90 * s
    cy, ch = 128 * s, 110 * s
    dr.rectangle([cx, cy, cx + cw, cy + ch], outline=(52, 54, 72), width=1)
    hb = np.histogram(pos[0], np.arange(0, SHAPE[0] + 1, 2))[0] / pos.shape[1]
    for i, v in enumerate(hb):
        X = cx + cw * (i + 0.5) / len(hb)
        Y = cy + ch - ch * min(v / 0.08, 1)
        r = 2.4 * s
        dr.ellipse([X - r, Y - r, X + r, Y + r], fill=(255, 200, 120))
    dr.text((cx, cy - 17 * s), "dove sta l'energia lungo la direzione dell'urto (punti = palline)", fill=(150, 155, 185), font=font(11.5))
    small = run == "video_piccolo"
    p2 = ver.get("parte2", {})
    Dv = p2.get("D", {})
    lines = [
        ("Cosa vedi", (245, 245, 252), True),
        ("• ogni pallina è un pezzetto fisso dell'energia del campo dei gluoni", (205, 208, 230), False),
        ("• si muove con il flusso di energia del campo (vettore di Poynting di Yang-Mills)", (205, 208, 230), False),
        ("• colore (tavolozza dell'atomo zig-zag): viola = elettrico e magnetico in pari (onda libera),", (205, 208, 230), False),
        ("  oro = più elettrico, ciano = più magnetico (E²/(E²+B²): non dipende dalla gauge)", (205, 208, 230), False),
        ("• due pacchetti di colore diverso (T¹ e T²) si vengono incontro lungo x", (205, 208, 230), False),
        ("", None, False),
        ("Il punto", (245, 245, 252), True),
    ]
    if run == "video_tre":
        lines[5] = ("• tre pacchetti (colori T¹, T², T⁴, che non commutano) lungo x, y, z", (205, 208, 230), False)
        lines += [("ampiezza 1: si incontrano al centro e la colla si urta da sola;", (205, 208, 230), False),
                  ("corsa aggiunta dopo la pre-registrazione, solo per il film (dichiarata)", (255, 205, 120), False)]
    elif small:
        lines += [("ampiezza 0,05: si attraversano quasi come la luce (lineare)", (205, 208, 230), False),
                  (f"energia spostata rispetto all'onda lineare: D = {num(100 * Dv.get('a005', float('nan')), '{:.2f}')} %", (170, 230, 170), False)]
    else:
        lines += [("ampiezza 1: la colla porta colore e si attira da sola: si urtano,", (205, 208, 230), False),
                  ("nasce campo nel terzo colore T³ = [T¹, T²] e l'energia si sparpaglia", (205, 208, 230), False),
                  (f"energia spostata rispetto all'onda lineare: D = {num(100 * Dv.get('a1', float('nan')), '{:.0f}')} %", (255, 190, 130), False),
                  (f"stesso colore (T¹ e T¹), stessa ampiezza: D = {num(100 * (p2.get('G5_colore', {}).get('D_stesso_1.0') or float('nan')), '{:.0f}')} %", (205, 208, 230), False)]
    if "G6_scala" in p2:
        lines.append((f"D cresce come ampiezza^{num(p2['G6_scala']['pendenza'], '{:.2f}')} (previsto 2)", (205, 208, 230), False))
    lines += [
        ("", None, False),
        ("Cosa NON è", (245, 245, 252), True),
        ("• un campo classico: niente gluoni uno per uno, niente quantistica", (255, 205, 120), False),
        ("• niente confinamento in tempo reale (problema del segno): qui le onde vanno lontano", (255, 205, 120), False),
        ("• il 'colore' T¹, T² è quello delle condizioni iniziali; l'urto (non commutano) è invariante", (255, 205, 120), False),
        ("• palline fra due fotogrammi salvati: posizione interpolata (rallentatore)", (150, 155, 185), False),
    ]
    y = 262
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
    vp = os.path.join(LAB, "out", "qg_verdetto.json")
    ver = json.load(open(vp)) if os.path.exists(vp) else {}
    rf = os.path.join(SCR, "qg_film_refs.json")
    if os.path.exists(rf):
        REFS.update(json.load(open(rf)))
    if args[0] == "still":
        run, k = args[1], int(args[2])
        frame(run, k, 0.0, 0.5, ver).save(args[3])
    elif args[0] == "refs":
        for run in RUNS:
            REFS.pop(run, None)
            frame(run, 0, 0.0, 0.0, ver)
        json.dump({k: {kk: float(vv) for kk, vv in v.items()} for k, v in REFS.items()}, open(rf, "w"))
        print(REFS)
    elif args[0] in ("batch", "batchrun"):
        if args[0] == "batchrun":
            run = args.pop(1)
            L = lista([run])
            dd = os.path.join(SCR, f"qg_film_{run}_{int(W0 * SCALA)}")
        else:
            L = lista()
            dd = os.path.join(SCR, f"qg_film_{int(W0 * SCALA)}")
        a, b = int(args[1]), int(args[2])
        os.makedirs(dd, exist_ok=True)
        t0 = time.time()
        for i in range(a, min(b, len(L))):
            fn = os.path.join(dd, f"{i:05d}.jpg")
            if os.path.exists(fn):
                continue
            run, k, fr = L[i]
            frame(run, k, fr, i / len(L), ver).save(fn + ".tmp.jpg", quality=93)
            os.replace(fn + ".tmp.jpg", fn)
            if time.time() - t0 > 470:
                break
        print("batch", a, b, len(L), f"{time.time() - t0:.0f}s")

"""vz_render.py - VIDEO SUPERCRITICO: draw the frames. Every pixel of the scene comes from balls (soft additive
gaussian glows, depth weight, wide halo), the same recipe as lab/vl_render.py. No lines, trails or surfaces in the
scene; the right panel is UI (text + two small charts).

    python3 lab/vz_render.py still IV out.png [--scala 1.0]
    python3 lab/vz_render.py refs IV                         # fix the exposure once (from frame IV)
    python3 lab/vz_render.py batch IV0 IV1 [--scala 1.0]     # frames to <SCR>/frames_<W>/, skips existing

Visual knobs (aesthetic only, counted in docs/VIDEO_SUPERCRITICO.md): SPOT, GLOW, GAIN, GLOW_MIX, DEPTH, camera,
palette, colour logistic width LOGI_W, whitening WHITE, cut-away.
"""
import json
import os
import sys
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter, zoom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vz_comune import (SCR, T, DT, NF, RB, RC, RHO, QE, SC, RW, ER, Z, MAN, R0_LOG, RMAX, TV, NV,  # noqa: E402
                       palline, raggi, disegno, Q_at)

W0, H0 = 1920, 1080
SCALA = 1.0
FONT = "/usr/share/fonts/opentype/inter/Inter-Regular.otf"
FONTB = "/usr/share/fonts/opentype/inter/Inter-SemiBold.otf"
LAB = os.path.dirname(os.path.abspath(__file__))
# aesthetics
SPOT = 2.2           # px at 1080p for a ball at the centre distance and reference density
GLOW = 14.0
GAIN = 2.2
GLOW_MIX = 0.5
DEPTH = 14.0
WORLD = 10.0         # drawn radius of r = 40 lambda_C
CAM_D = 31.0
CAM_FOV = 40.0
LOGI_W = 0.25        # width in ln r of the colour change at r_w
WHITE = 0.18         # whitening with |v|/c
KEEP = 70000         # balls drawn as sharp glows (all of them go into the wide halo)
GAMMA = 0.8          # exposure compression: lum = 1 - exp(-gain (D/ref)^GAMMA) (declared)
CUT = os.environ.get("VZ_CUT", "meta")   # "meta": the half of the sphere facing the camera is not drawn
C_VAC = np.array([0.22, 0.48, 1.00])     # weight still around the nucleus: charged vacuum (-e), hole-theory reading
C_POS = np.array([1.00, 0.46, 0.08])     # weight that has left r_w: the positron (+e)
REFS_FILE = os.path.join(SCR, "vz_refs.json")
REFS = json.load(open(REFS_FILE)) if os.path.exists(REFS_FILE) else {}
ESAMI = json.load(open(os.path.join(LAB, "out", "vz_esami.json"))) if os.path.exists(os.path.join(LAB, "out", "vz_esami.json")) else {}
ETICHETTA = os.environ.get("VZ_ETICHETTA", "BOZZA")

M, DIRS = palline()


def camera(eye, centre, fov, w, h):
    fwd = np.asarray(centre, float) - eye
    fwd /= np.linalg.norm(fwd)
    right = np.cross(fwd, np.array([0.0, 0.0, 1.0]))
    right /= np.linalg.norm(right)
    up = np.cross(right, fwd)
    f = 0.5 * h / np.tan(np.radians(fov) / 2)
    return eye, fwd, right, up, f, w, h


def project(p, cam):
    eye, fwd, right, up, f, w, h = cam
    d = p - eye[:, None]
    z = (d * fwd[:, None]).sum(0)
    x = (d * right[:, None]).sum(0)
    y = (d * up[:, None]).sum(0)
    zz = np.maximum(z, 1e-3)
    return w / 2 + f * x / zz, h / 2 - f * y / zz, z


def splat(img, sx, sy, wts, col, w, h):
    ix = np.round(sx).astype(np.int64)
    iy = np.round(sy).astype(np.int64)
    ok = (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h) & (wts > 0)
    lin = iy[ok] * w + ix[ok]
    ww = wts[ok]
    for c in range(3):
        img[c] += np.bincount(lin, ww * col[ok, c], w * h).reshape(h, w)
    img[3] += np.bincount(lin, ww, w * h).reshape(h, w)


def blur(a, sig):
    if sig < 0.3:
        return a
    f = 4 if sig >= 8 else (2 if sig >= 4 else 1)
    if f == 1:
        return gaussian_filter(a, sig, truncate=3.0)
    h, w = a.shape
    H2, W2 = -(-h // f) * f, -(-w // f) * f
    p = np.zeros((H2, W2))
    p[:h, :w] = a
    small = p.reshape(H2 // f, f, W2 // f, f).mean((1, 3))
    small = gaussian_filter(small, sig / f, truncate=3.0)
    return zoom(small, f, order=1)[:h, :w]


def additivo(img, gain, ref=None):
    D = img[3]
    if ref is None:
        ref = np.percentile(D[D > 0], 99.7) if (D > 0).any() else 1.0
    return 1 - np.exp(-gain * np.maximum(img[:3] / ref, 0) ** GAMMA), ref


def cam_scena(iv, w, h):
    u = iv / max(NV - 1, 1)
    az = np.radians(-60 + 80 * u)
    el = np.radians(16 + 8 * np.sin(np.pi * u))
    eye = CAM_D * np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    return camera(eye, (0, 0, 0), CAM_FOV, w, h)


def stato(t):
    """radii, radial speeds, colours of the balls at time t"""
    r = raggi(M, t)
    t2 = min(t + 0.05, T[-1])
    t1 = t2 - 0.05
    v = (raggi(M, t2) - raggi(M, t1)) / 0.05 if t2 > t1 else np.zeros_like(r)
    ok = np.isfinite(r)
    return r, np.where(np.isfinite(v), v, 0.0), ok


def scena(iv, w, h):
    t = TV[iv]
    cam = cam_scena(iv, w, h)
    s = SCALA
    r, v, ok = stato(t)
    r, v, d = r[ok], v[ok], DIRS[:, ok]
    gid = np.nonzero(ok)[0]              # global ball index (fixed sharp subset, no flicker)
    R = disegno(r) * WORLD
    p = d * R
    sp = 1 / (1 + np.exp(-(np.log(r) - np.log(RW)) / LOGI_W))          # 0 = around the nucleus, 1 = left r_w
    col = (1 - sp)[:, None] * C_VAC + sp[:, None] * C_POS
    wv = WHITE * np.clip(np.abs(v), 0, 1)[:, None]
    col = (1 - wv) * col + wv
    # local density in the drawing (balls per drawn volume): sparser -> larger, softer glow
    hb, eb = np.histogram(R, np.linspace(0, WORLD, 161))
    vol = 4 / 3 * np.pi * (eb[1:] ** 3 - eb[:-1] ** 3)
    nd = hb / vol
    nref = len(M) / (4 / 3 * np.pi * (0.5 * WORLD) ** 3)
    dn = np.clip((nd[np.clip(np.digitize(R, eb) - 1, 0, 159)] / nref) ** (-1 / 3), 0.6, 2.2)
    if CUT == "meta":
        far = (p * cam[1][:, None]).sum(0) > 0
        p, col, dn, gid = p[:, far], col[far], dn[far], gid[far]
    sx, sy, z = project(p, cam)
    wt = np.exp(-np.clip(z - CAM_D + 4.0, 0, None) / DEPTH)
    rpx = SPOT * s * (CAM_D / np.maximum(z, 3.0)) ** 1.3 * dn
    b = np.clip(np.round(np.log2(np.maximum(rpx / 0.7, 1.0)) * 2).astype(int), 0, 7)
    layers = np.zeros((8, 4, h, w))
    img2 = np.zeros((4, h, w))
    splat(img2, sx, sy, wt, col, w, h)
    sharp = gid < KEEP                   # mass coordinates are shuffled: the first KEEP balls are a fair sample
    for ib in np.unique(b):
        m = (b == ib) & sharp
        splat(layers[ib], sx[m], sy[m], wt[m] / dn[m] ** 0.5, col[m], w, h)
    img = np.zeros((4, h, w))
    for ib in range(8):
        if layers[ib][3].any():
            img += np.stack([blur(layers[ib][c], 0.7 * 2 ** (ib / 2)) for c in range(4)])
    balls, ref = additivo(img, GAIN, REFS.get("b", 0) / s ** 2 or None)
    glow, ref2 = additivo(np.stack([blur(img2[c], GLOW * s) for c in range(4)]), 1.0, REFS.get("g", 0) / s ** 2 or None)
    REFS.setdefault("b_now", ref * s ** 2)
    REFS.setdefault("g_now", ref2 * s ** 2)
    out = 1 - (1 - np.clip(balls, 0, 1)) * (1 - GLOW_MIX * np.clip(glow, 0, 1))
    return out, cam, r


# ------------------------------------------------------------------ UI
def font(sz, bold=False):
    return ImageFont.truetype(FONTB if bold else FONT, max(8, int(round(sz * SCALA))))


def fase_testo(t):
    if t < 6:
        return "1 · lo stato 1s del nucleo Z = 200 è già sotto −mc²: tuffato nel mare di Dirac"
    if t < 70:
        return "2 · la risonanza si svuota piano (come e^(−Γt), 1/Γ = 33 λ_C/c) e il peso esce come un guscio: è il positrone"
    if t < 150:
        return "3 · il positrone corre via (spugna a r = 40 λ_C); attorno al nucleo resta il vuoto carico"
    return "4 · quasi tutto è uscito: resta lo 0,7 % dentro r_w (in parte livelli legati veri, che non decadono)"


def num(x, f="{:.3f}"):
    return f.format(x).replace(".", ",")


def pannello(dr, x0, y0, w, h, iv, rb):
    s = SCALA
    t = TV[iv]

    def Tx(x, y, txt, sz=13, col=(220, 222, 235), bold=False):
        dr.text((x0 + x * s, y0 + y * s), txt, fill=col, font=font(sz, bold))

    Tx(24, 22, "L'atomo supercritico: l'1s si tuffa", 24, (245, 245, 252), True)
    Tx(24, 54, "e il vuoto emette un positrone", 24, (245, 245, 252), True)
    if ETICHETTA:
        Tx(w / s - 110, 26, ETICHETTA, 18, (255, 150, 60), True)
    Tx(24, 90, f"Z = 200 (nucleo a guscio di 12 fm) · regola zig-zag radiale, passo a = 0,002 λ_C", 13, (170, 175, 200))
    rate = (TV[min(iv + 1, NV - 1)] - TV[max(iv - 1, 0)]) / 2 * 24
    slow = 1.0 / (rate * 1.288e-21)
    ex = int(np.floor(np.log10(slow)))
    sup = str(ex).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))
    Tx(24, 116, f"t = {num(t, '{:.1f}')} λ_C/c = {num(t * 1.288, '{:.0f}')} zeptosecondi · rallentato {num(slow / 10 ** ex, '{:.0f}')}·10{sup} volte", 15, (235, 235, 245))
    Tx(24, 142, fase_testo(t), 13.5, (255, 214, 140))
    lw = max(1, int(round(2 * s)))
    # chart 1: radial profile now (log r): line = data, dots = balls
    cx, cw = x0 + 70 * s, w - 100 * s
    cy, ch = y0 + 192 * s, 150 * s
    lr0, lr1 = np.log10(0.01), np.log10(40)
    def X(r): return cx + cw * (np.log10(r) - lr0) / (lr1 - lr0)
    edges = np.geomspace(0.01, 40, 41)
    qd = np.diff(np.interp(edges, RB, Q_at(t)))
    rr = rb[np.isfinite(rb)]
    hb = np.histogram(rr, edges)[0] / len(M)
    dl = np.log10(edges[1]) - np.log10(edges[0])
    def Y(v): return cy + ch - ch * float(np.clip((np.log10(max(v / dl, 1e-5)) + 4) / 4.3, 0, 1))   # log y: 1e-4 .. 2
    dr.rectangle([cx, cy, cx + cw, cy + ch], outline=(52, 54, 72), width=1)
    xm = np.sqrt(edges[1:] * edges[:-1])
    dr.line([(X(a), Y(b)) for a, b in zip(xm, qd)], fill=(120, 130, 170), width=lw)
    for a, b in zip(xm, hb):
        rr_ = 2.6 * s
        c = (110, 160, 255) if a < RW else (255, 160, 80)
        dr.ellipse([X(a) - rr_, Y(b) - rr_, X(a) + rr_, Y(b) + rr_], fill=c)
    xw = X(RW)
    for yy in np.arange(cy, cy + ch, 6 * s):
        dr.point((xw, yy), fill=(90, 90, 110))
    dr.text((xw + 3 * s, cy + ch - 16 * s), "r_w", fill=(150, 150, 170), font=font(11))
    for rv, lab in ((0.01, "0,01"), (0.1, "0,1"), (1, "1"), (10, "10"), (40, "40 λ_C")):
        dr.text((X(rv) - 6 * s, cy + ch + 3 * s), lab, fill=(120, 120, 140), font=font(11))
    dr.text((x0 + 12 * s, cy + 2 * s), "peso", fill=(200, 200, 225), font=font(12.5, True))
    dr.text((x0 + 12 * s, cy + 19 * s), "per", fill=(150, 150, 175), font=font(11))
    dr.text((x0 + 12 * s, cy + 33 * s), "decade", fill=(150, 150, 175), font=font(11))
    dr.text((cx, cy - 17 * s), "dove sta il peso adesso (asse y logaritmico 10⁻⁴…2): riga = dati, punti = palline", fill=(150, 155, 185), font=font(11.5))
    # chart 2: weight inside r_w and passed r = 30 vs time
    cy2, ch2 = cy + ch + 30 * s, 110 * s
    def Xt(tv): return cx + cw * tv / T[-1]
    def Yt(v): return cy2 + ch2 - ch2 * v
    dr.rectangle([cx, cy2, cx + cw, cy2 + ch2], outline=(52, 54, 72), width=1)
    pin = np.array(SC["P_in_rw"])
    p30 = np.array(SC["passato_r30"])
    dr.line([(Xt(a), Yt(b)) for a, b in zip(T, pin)], fill=(60, 80, 130), width=lw)
    dr.line([(Xt(a), Yt(b)) for a, b in zip(T, p30)], fill=(130, 90, 50), width=lw)
    sp = ESAMI.get("serie_palline", {})
    kk = int(np.floor(t / DT + 1e-9))
    for k in range(0, kk + 1, 2):
        if sp:
            for val, c in ((sp["P_in_rw"][k], (120, 170, 255)), (sp["passato_r30"][k], (255, 165, 90))):
                rr_ = 1.8 * s
                dr.ellipse([Xt(T[k]) - rr_, Yt(val) - rr_, Xt(T[k]) + rr_, Yt(val) + rr_], fill=c)
    pnow = float(np.interp(t, T, pin))
    p30now = float(np.interp(t, T, p30))
    dr.text((x0 + 12 * s, cy2 + 2 * s), "dentro r_w", fill=(120, 170, 255), font=font(11.5, True))
    dr.text((x0 + 12 * s, cy2 + 18 * s), num(pnow), fill=(150, 190, 255), font=font(12))
    dr.text((x0 + 12 * s, cy2 + 44 * s), "oltre 30", fill=(255, 165, 90), font=font(11.5, True))
    dr.text((x0 + 12 * s, cy2 + 60 * s), num(p30now), fill=(255, 190, 130), font=font(12))
    dr.text((cx + cw - 300 * s, cy2 + ch2 / 2), "riga = dati, punti = palline (0 … 1)", fill=(130, 135, 165), font=font(11.5))
    dr.text((cx, cy2 + ch2 + 3 * s), "0", fill=(120, 120, 140), font=font(11))
    dr.text((cx + cw - 50 * s, cy2 + ch2 + 3 * s), "t = 240", fill=(120, 120, 140), font=font(11))
    # captions
    e = ESAMI
    v4 = e.get("V4_velocita_dati", {})
    v5 = e.get("V5_velocita_cinematica", {})
    v2 = e.get("V2_ode", {})
    lines = [
        ("300 000 palline = il peso dello stato tuffato; ognuna è un bagliore (niente scie)", (205, 208, 230)),
        ("blu = peso ancora attorno al nucleo → vuoto carico (−e) · arancio = peso uscito da r_w → positrone (+e)", (205, 208, 230)),
        ("il colore è la lettura della teoria delle buche (ZSTAR §10): i dati hanno un solo segno, non è un segno calcolato", (255, 205, 120)),
        ("i dati sono radiali: direzioni messe a caso, sfera perfetta per costruzione (non calcolata)", (255, 205, 120)),
        ("raggio disegnato in scala logaritmica (compressa): l'interno (r < 0,7 λ_C) occupa metà del raggio", (255, 205, 120)),
        ("ogni pallina tiene la sua quota di carica: r(t) da Q(<r, t), cioè mossa con v = j/ρ, continuità esatta", (205, 208, 230)),
        ("più bianca = più veloce · bagliore più grande = palline più rade · esposizione fissa (il buio finale è vero)", (170, 175, 200)),
    ]
    if v4:
        lines.append((f"velocità del positrone (20–30 λ_C): palline {num(v4['palline'])} c, dati {num(v4['dati_sum_j_su_sum_rho'])} c, "
                      f"cinematica da E_r {num(v5['beta_cinematica_r25'])} c: PASS", (170, 230, 170)))
    if v2:
        lines.append((f"FAIL dichiarato: integrare v = j/ρ a passi devia dentro r < 5 (mediana {num(v2['80']['mediana'] * 100, '{:.1f}')} % a t = 80):", (255, 150, 140)))
        lines.append(("i dati esportati a 4 cifre perdono la continuità a 10⁻⁴ per fotogramma; le palline usano Q (esatto)", (255, 150, 140)))
    lines.append(("una particella sola (stessa equazione di Müller–Rafelski–Greiner), niente calcolo a molti corpi né QED", (170, 175, 200)))
    lines.append((f"dati: lab/out/zstar_film.json (commit fc8a08c) · λ_C = 386 fm, λ_C/c = 1,29·10⁻²¹ s" + (f" · {ETICHETTA.lower()}" if ETICHETTA else ""), (120, 125, 150)))
    y = (cy2 - y0) / s + ch2 / s + 26
    for txt, col in lines:
        Tx(24, y, txt, 12.2, col)
        y += 17.6


def frame(iv):
    s = SCALA
    W, H = int(W0 * s), int(H0 * s)
    bg = np.array([0.010, 0.011, 0.020])
    canvas = np.zeros((H, W, 3))
    canvas[:] = bg
    wa = H
    sc, cam, rb = scena(iv, wa, wa)
    sc = np.clip(sc + bg[:, None, None] * (1 - sc), 0, 1) ** (1 / 1.1)
    canvas[:, :wa] = sc.transpose(1, 2, 0)
    im = Image.fromarray((canvas * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im)
    dr.text((22 * s, 18 * s), "il nucleo al centro (non disegnato) e il suo 1s, Z = 200", fill=(200, 202, 225), font=font(15))
    dr.text((22 * s, 40 * s), "raggio in scala logaritmica · bordo della sfera = 40 λ_C = 15 000 fm", fill=(140, 145, 170), font=font(12.5))
    # radius labels (text only), projected along the camera's right vector
    eye, fwd, right, up, f, w, h = cam
    for rv, lab in ((0.1, "0,1"), (RW, "r_w = 0,7"), (4.0, "4"), (40.0, "40 λ_C")):
        p = (right * disegno(np.array([rv]))[0] * WORLD)[:, None]
        sx, sy, _ = project(p, cam)
        dr.text((sx[0] - 8 * s, sy[0] + 6 * s), lab, fill=(110, 112, 135), font=font(11))
    dr.rectangle([wa, 0, W, H], fill=(8, 9, 16))
    pannello(dr, wa, 0, W - wa, H, iv, rb)
    return im


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--scala" in args:
        i = args.index("--scala")
        SCALA = float(args[i + 1])
        del args[i:i + 2]
    if args[0] == "refs":
        iv = int(args[1])
        REFS.clear()
        scena(iv, int(H0 * SCALA), int(H0 * SCALA))
        out = {"b": REFS["b_now"], "g": REFS["g_now"], "iv": iv}
        json.dump(out, open(REFS_FILE, "w"))
        print(out)
    elif args[0] == "still":
        iv = int(args[1])
        t0 = time.time()
        frame(iv).save(args[2])
        print("still", iv, TV[iv], f"{time.time() - t0:.1f}s")
    elif args[0] == "batch":
        a, b = int(args[1]), int(args[2])
        d = os.path.join(SCR, f"frames_{int(W0 * SCALA)}")
        os.makedirs(d, exist_ok=True)
        t0 = time.time()
        budget = float(os.environ.get("VZ_BUDGET", 480))
        for iv in range(a, min(b, NV)):
            f = os.path.join(d, f"{iv:05d}.png")
            if os.path.exists(f):
                continue
            frame(iv).save(f + ".tmp.png")
            os.replace(f + ".tmp.png", f)
            if time.time() - t0 > budget:
                break
        print("batch", a, b, f"{time.time() - t0:.0f}s")

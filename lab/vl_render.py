"""vl_render.py - VIDEO ATOMO: draw the frames. Every pixel of the two scenes comes from balls (soft gaussian spots,
depth-weighted blending, glow), as in lab/fp_render.py. No lines, tubes or surfaces in the scenes; the right-bottom
panel is a plain chart with captions (UI, not scene).

    python3 lab/vl_render.py still IV out.png [--scala 1.0]
    python3 lab/vl_render.py batch IV0 IV1 [--scala 1.0]     # frames to <SCR>/frames_<W>/, skips existing

Visual knobs (aesthetic only, counted in docs/VIDEO_ATOMO.md): SPOT_E, GLOW_E, WAKE_E, COL_FULL, SPOT_L, GLOW_L,
WAKE_L, STRETCH, camera paths, palette, exposure gains.
"""
import json
import os
import sys
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vl_comune import SCR, DT, A, NA, CHUNK, timeline, idx_light, trilinear, catmull  # noqa: E402

W0, H0 = 1920, 1080
SCALA = 1.0
FONT = "/usr/share/fonts/opentype/inter/Inter-Regular.otf"
FONTB = "/usr/share/fonts/opentype/inter/Inter-SemiBold.otf"
# aesthetics
SPOT_E = 2.6         # px (at 1080p, for a ball at the centre distance)
GLOW_E = 16.0
WAKE_E = (1.0, 0.55, 0.32, 0.18, 0.10, 0.05)
COL_FULL = 0.25      # colour full when |drho| / (0.5 rho_1s(x) + 0.05 rho_1s,max) reaches this (relative change, softened)
SPOT_L = 1.1
GLOW_L = 10.0
WAKE_L = (1.0, 0.5, 0.25)
STRETCH = 7.0
GAIN_E = 0.9
GAIN_L = 1.5
GLOW_MIX_E = 0.7
GLOW_MIX_L = 0.45
DEPTH_E = 9.0        # lambda_C: depth scale of the transparency (near balls cover far ones)
DEPTH_L = 260.0

TV = timeline()
REFS_FILE = os.path.join(SCR, "vl_refs.json")
REFS = json.load(open(REFS_FILE)) if os.path.exists(REFS_FILE) else {}
NV = len(TV)
MAN = json.load(open("/mnt/project-files/ripresa/dati/atomo_accoppiato/manifest.json"))


def _csv(path):
    import csv
    with open(path) as f:
        r = csv.reader(f)
        head = next(r)
        rows = np.array([[float(v) for v in row] for row in r])
    return {h: rows[:, i] for i, h in enumerate(head)}


FR = _csv("/mnt/project-files/ripresa/dati/atomo_accoppiato/fotogrammi.csv")


# ------------------------------------------------------------------ data access
class Dati:
    def __init__(self):
        self.el = {}
        self.lu = {}
        self.rho = np.load(os.path.join(SCR, "rho_s.npy"), mmap_mode="r")
        self.rho0max = float(np.asarray(self.rho[0]).max())
        self.E = np.load(os.path.join(SCR, "E.npy"), mmap_mode="r")
        self.B = np.load(os.path.join(SCR, "B.npy"), mmap_mode="r")
        self.cacheL = {}

    def el_pos(self, iv):
        c = iv // CHUNK
        if c not in self.el:
            if len(self.el) > 2:
                self.el.pop(min(self.el))
            self.el[c] = np.load(os.path.join(SCR, f"el_{c:03d}.npy"), mmap_mode="r")
        return np.asarray(self.el[c][iv - c * CHUNK], np.float64)

    def lu_pos(self, iv):
        c = iv // CHUNK
        if c not in self.lu:
            if len(self.lu) > 2:
                self.lu.pop(min(self.lu))
            self.lu[c] = np.load(os.path.join(SCR, f"lu_{c:03d}.npz"))
        z = self.lu[c]
        i = iv - c * CHUNK
        return np.asarray(z[f"p{i}"], np.float64), np.asarray(z[f"k{i}"])

    def drho(self, p, t):
        k = int(min(np.floor(t / DT + 1e-9), 199))
        s = t / DT - k
        r0 = np.asarray(self.rho[0], np.float64)
        R = (1 - s) * np.asarray(self.rho[k], np.float64) + s * np.asarray(self.rho[k + 1], np.float64) - r0
        R = R / (0.5 * r0 + 0.05 * self.rho0max)
        return trilinear(R[None], p / A + NA / 2 - 0.5)[0]

    def lframe(self, k):
        if k not in self.cacheL:
            if len(self.cacheL) > 6:
                self.cacheL.pop(min(self.cacheL))
            self.cacheL[k] = np.concatenate([np.asarray(self.E[k], np.float64), np.asarray(self.B[k], np.float64)])
        return self.cacheL[k]

    def EB(self, p, t):
        F = catmull(201, t, self.lframe)
        f = trilinear(F, idx_light(p))
        return f[:3], f[3:]


# ------------------------------------------------------------------ camera + splat
def camera(eye, centre, fov, w, h):
    fwd = np.asarray(centre, float) - eye
    fwd /= np.linalg.norm(fwd)
    up0 = np.array([0.0, 0.0, 1.0])
    right = np.cross(fwd, up0)
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


def tone(img, spot, gain, ref=None):
    rgbw = np.stack([gaussian_filter(img[c], spot) for c in range(4)])
    D = rgbw[3]
    col = rgbw[:3] / np.maximum(D, 1e-12)
    if ref is None:
        ref = np.percentile(D[D > 0], 99.7) if (D > 0).any() else 1.0
    lum = 1 - np.exp(-gain * D / ref)
    return col * lum, ref


def palette(hh, sat=0.85):
    a = 2 * np.pi * hh
    rgb = np.stack([0.5 + 0.5 * np.cos(a), 0.5 + 0.5 * np.cos(a - 2.094), 0.5 + 0.5 * np.cos(a + 2.094)], -1)
    return (1 - sat) + sat * rgb


# ------------------------------------------------------------------ electron scene
C_BASE = np.array([0.78, 0.80, 1.00])
C_UP = np.array([1.00, 0.62, 0.10])
C_DN = np.array([0.18, 0.52, 1.00])


CAM_D = 27.0
CAM_FOV = 38.0
KEEP_E = 30000
NBIN = 5


def cam_atomo(iv, w, h):
    u = iv / (NV - 1)
    az = np.radians(-70 + 55 * u)
    el = np.radians(14 + 6 * np.sin(np.pi * u))
    eye = CAM_D * np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    return camera(eye, (0, 0, 0), CAM_FOV, w, h)


def scena_atomo(dat, iv, w, h):
    """balls in NBIN depth slices: nearer slices get bigger, softer spots (perspective size of a ball)"""
    cam = cam_atomo(iv, w, h)
    s = SCALA
    zb = CAM_D + np.array([-12.0, -5.0, -2.0, 2.0, 5.0, 99.0])   # slice edges (distance from the eye)
    layers = np.zeros((NBIN, 4, h, w))
    img2 = np.zeros((4, h, w))
    for kk, wk in enumerate(WAKE_E):
        j = iv - kk
        if j < 0:
            break
        pall = dat.el_pos(j)
        drall = dat.drho(pall, TV[j]) / COL_FULL
        p, dr = pall[:, :KEEP_E], drall[:KEEP_E]
        a = np.clip(np.abs(dr), 0, 1) ** 0.8
        tgt = np.where((dr > 0)[:, None], C_UP, C_DN)
        col = (1 - a)[:, None] * C_BASE + a[:, None] * tgt
        sx, sy, z = project(p, cam)
        wt = wk * np.exp(-np.clip(z - CAM_D + 3.0, 0, None) / DEPTH_E)
        b = np.clip(np.searchsorted(zb, z) - 1, 0, NBIN - 1)
        for ib in range(NBIN):
            m = b == ib
            splat(layers[ib], sx[m], sy[m], wt[m], col[m], w, h)
        if kk == 0:
            a2 = np.clip(np.abs(drall), 0, 1) ** 0.8
            col2 = (1 - a2)[:, None] * C_BASE + a2[:, None] * np.where((drall > 0)[:, None], C_UP, C_DN)
            sx2, sy2, z2 = project(pall, cam)
            splat(img2, sx2, sy2, np.exp(-np.clip(z2 - CAM_D + 3.0, 0, None) / DEPTH_E), col2, w, h)
    img = np.zeros((4, h, w))
    for ib in range(NBIN):
        zc = 0.5 * (zb[ib] + min(zb[ib + 1], CAM_D + 8))
        sp = SPOT_E * s * (CAM_D / max(zc, 3.0)) ** 1.3
        img += np.stack([gaussian_filter(layers[ib][c], sp) for c in range(4)])
    balls, ref = tone(img, 0.01, GAIN_E, REFS.get("eb", 0) / s ** 2 or None)
    glow, ref2 = tone(img2, GLOW_E * s, 1.0, REFS.get("eg", 0) / s ** 2 or None)
    REFS.setdefault("eb_now", ref * s ** 2)
    REFS.setdefault("eg_now", ref2 * s ** 2)
    out = 1 - (1 - np.clip(balls, 0, 1)) * (1 - GLOW_MIX_E * np.clip(glow, 0, 1))
    return out, cam


# ------------------------------------------------------------------ light scene
def cam_luce(iv, w, h):
    u = iv / (NV - 1)
    az = np.radians(-112 + 14 * u)
    el = np.radians(24)
    D = 760.0
    eye = D * np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    return camera(eye, (0, 0, -10), 34, w, h)




def scena_luce(dat, iv, w, h):
    t = TV[iv]
    cam = cam_luce(iv, w, h)
    img = np.zeros((4, h, w))
    img2 = np.zeros((4, h, w))
    s = SCALA
    for kk, wk in enumerate(WAKE_L):
        j = iv - kk
        if j < 0:
            break
        p, kb = dat.lu_pos(j)
        keep = (p[1] > -8.0) & (p[1] < 72.0)   # cut-away: a slab 0 < y < 72 (the half facing the camera is not drawn)
        p, kb = p[:, keep], kb[keep]
        E, B = dat.EB(p, TV[j])
        sx, sy, z = project(p, cam)
        ex, ey, _ = project(p + 30.0 * E / (np.sqrt((E ** 2).sum(0)) + 1e-30), cam)
        dx, dy = ex - sx, ey - sy
        nr = np.hypot(dx, dy) + 1e-9
        ux, uy = dx / nr, dy / nr
        # hue from the direction of E (around the y axis, the view direction): +z rose, -z teal
        hue = (np.arctan2(E[2], E[0]) / (2 * np.pi) - 0.33) % 1.0
        col = palette(hue)
        col = np.where((kb == 1)[:, None], 0.25 + 0.75 * col, col)
        wt = wk * np.exp(-np.clip(z - 760 + 200, 0, None) / DEPTH_L)
        L = STRETCH * s
        for tt, ww in zip(np.linspace(-0.5, 0.5, 5), (0.1, 0.22, 0.36, 0.22, 0.1)):
            splat(img, sx + tt * L * ux, sy + tt * L * uy, ww * wt, col, w, h)
        if kk == 0:
            splat(img2, sx, sy, wt, col, w, h)
    # the atom at its true size inside the light (the same electron balls, a sample)
    pe = dat.el_pos(iv)[:, ::6]
    sx, sy, z = project(pe, cam)
    splat(img, sx, sy, np.full(sx.size, 0.6), np.tile(np.array([1.0, 0.93, 0.85]), (sx.size, 1)), w, h)
    balls, ref = tone(img, SPOT_L * s, GAIN_L, REFS.get("lb", 0) / s ** 2 or None)
    glow, ref2 = tone(img2, GLOW_L * s, 1.0, REFS.get("lg", 0) / s ** 2 or None)
    REFS.setdefault("lb_now", ref * s ** 2)
    REFS.setdefault("lg_now", ref2 * s ** 2)
    out = 1 - (1 - np.clip(balls, 0, 1)) * (1 - GLOW_MIX_L * np.clip(glow, 0, 1))
    return out, cam


# ------------------------------------------------------------------ UI
SERIE = {}


def dip_palline(dat, iv):
    return float(dat.el_pos(iv)[2].mean())


def font(sz, bold=False):
    return ImageFont.truetype(FONTB if bold else FONT, max(8, int(round(sz * SCALA))))


def fase_testo(t):
    if t < 100:
        return "1 · arriva l'impulso di luce (da sinistra, polarizzato lungo z)"
    if t < 240:
        return "2 · l'elettrone viene scosso: la nuvola oscilla su e giù (dipolo)"
    return "3 · l'impulso è uscito: l'atomo oscilla da solo, riemette luce e si svuota piano"


def pannello_info(draw, x0, y0, w, h, iv, dzs):
    s = SCALA
    t = TV[iv]
    zs = MAN["unita"]["tempo"].split("=")[1].strip().split(" ")[0]
    tz = t * 1.28808866819e-21 / 1e-21
    rate = (TV[min(iv + 1, NV - 1)] - TV[max(iv - 1, 0)]) / 2 * 24  # lambda_C/c per second of video
    slow = 1.0 / (rate * 1.28808866819e-21)
    def T(x, y, txt, sz=15, col=(220, 222, 235), bold=False):
        draw.text((x0 + x * s, y0 + y * s), txt, fill=col, font=font(sz, bold))
    T(22, 14, "Un atomo che assorbe luce e la riemette", 23, (245, 245, 252), True)
    T(22, 46, "ione idrogenoide Z = 60 · elettrone e luce dalla regola zig-zag, calcolati insieme", 13.5, (170, 175, 200))
    ex = int(np.floor(np.log10(slow)))
    sup = str(ex).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))
    T(22, 72, f"t = {t:.0f} λ_C/c = {tz:.0f} zeptosecondi  ·  rallentato {slow / 10 ** ex:.0f}·10{sup} volte", 15, (235, 235, 245))
    T(22, 96, fase_testo(t), 14.5, (255, 214, 140))
    # charts: p_exc (from the computation) and <z> (balls vs computation)
    tt = FR["t"]
    pe = FR["p_exc"]
    dz = FR["dipolo_z"]
    cx, cw = x0 + 70 * s, w - 100 * s
    def X(tv): return cx + cw * tv / 1600.0
    lw = max(1, int(round(2 * s)))
    # p_exc
    cy, ch = y0 + 128 * s, 78 * s
    def Yp(v): return cy + ch - ch * v / 0.13
    draw.rectangle([cx, cy, cx + cw, cy + ch], outline=(52, 54, 72), width=1)
    pts = [(X(a), Yp(b)) for a, b in zip(tt, pe)]
    draw.line(pts, fill=(90, 62, 66), width=lw)
    m = tt <= t + 1e-9
    if m.sum() >= 2:
        draw.line([q for q, mm in zip(pts, m) if mm], fill=(255, 120, 110), width=max(1, int(round(2.6 * s))))
    pnow = float(np.interp(t, tt, pe))
    r = 4 * s
    draw.ellipse([X(t) - r, Yp(pnow) - r, X(t) + r, Yp(pnow) + r], fill=(255, 160, 150))
    draw.text((x0 + 14 * s, cy + 2 * s), "p_exc", fill=(255, 140, 130), font=font(13, True))
    draw.text((x0 + 14 * s, cy + 20 * s), f"{pnow:.3f}", fill=(255, 170, 160), font=font(13))
    draw.text((cx + cw - 250 * s, cy + 4 * s), "quota fuori dall'1s (dal calcolo)", fill=(150, 120, 120), font=font(11.5))
    # <z>
    cy2, ch2 = cy + ch + 10 * s, 62 * s
    def Yd(v): return cy2 + ch2 / 2 - (ch2 / 2) * v / 1.3
    draw.rectangle([cx, cy2, cx + cw, cy2 + ch2], outline=(52, 54, 72), width=1)
    draw.line([(X(a), Yd(b)) for a, b in zip(tt, dz)], fill=(46, 58, 82), width=1)
    for (a, b) in dzs:
        rr = 1.7 * s
        draw.ellipse([X(a) - rr, Yd(b) - rr, X(a) + rr, Yd(b) + rr], fill=(110, 170, 255))
    zb = dzs[-1][1] if dzs else 0.0
    draw.text((x0 + 14 * s, cy2 + 2 * s), "<z>", fill=(120, 175, 255), font=font(13, True))
    draw.text((x0 + 14 * s, cy2 + 20 * s), f"{zb:+.2f}", fill=(150, 195, 255), font=font(13))
    draw.text((cx + cw - 330 * s, cy2 + 4 * s), "dipolo λ_C: punti = palline, riga = calcolo", fill=(110, 125, 160), font=font(11.5))
    draw.text((cx, cy2 + ch2 + 3 * s), "0", fill=(120, 120, 140), font=font(11))
    draw.text((cx + cw - 52 * s, cy2 + ch2 + 3 * s), "t = 1600", fill=(120, 120, 140), font=font(11))
    # captions
    lines = [
        ("elettrone: 120 000 palline messe una volta sola dove |ψ|² è grande, mosse dalla corrente", (205, 208, 230)),
        ("oro = qui la densità cresce rispetto all'1s, blu = cala (pieno a circa ±25 % di cambio locale)", (255, 205, 120)),
        ("smussato su ~2 λ_C: tolta la trama a scacchiera del reticolo (4 valli)", (170, 175, 200)),
        ("luce: palline di energia (E²+B²)/2, colore = direzione di E; vista tagliata (fetta 0 < y < 72 λ_C)", (205, 208, 230)),
        ("luce dell'atomo (da t = 240): ogni pallina vale 1/1000 → energia ×1000", (255, 160, 200)),
        ("accoppiamento ×10 (K = 10) per vedere lo svuotamento; λ_C = 3,9·10⁻¹³ m", (170, 175, 200)),
        ("velocità del giro attorno a z: lettura del reticolo, non validata", (170, 175, 200)),
    ]
    y = 128 + 78 + 10 + 62 + 24
    for txt, col in lines:
        T(22, y, txt, 12.6, col)
        y += 19.5


def frame(dat, iv):
    s = SCALA
    W, H = int(W0 * s), int(H0 * s)
    canvas = np.zeros((H, W, 3))
    canvas[:] = np.array([0.012, 0.013, 0.022])
    wa = H
    atom, _ = scena_atomo(dat, iv, wa, wa)
    wl, hl = W - wa, int(round(600 * s))
    luce, _ = scena_luce(dat, iv, wl, hl)
    bg = np.array([0.012, 0.013, 0.022])[:, None, None]
    atom = np.clip(atom + bg * (1 - atom), 0, 1) ** (1 / 1.1)
    luce = np.clip(luce + bg * (1 - luce), 0, 1) ** (1 / 1.1)
    canvas[:, :wa] = atom.transpose(1, 2, 0)
    canvas[:hl, wa:] = luce.transpose(1, 2, 0)
    im = Image.fromarray((canvas * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im)
    # labels on the scenes (text only)
    dr.text((22 * s, 18 * s), "dentro l'atomo  ·  lato 24 λ_C  ·  z in su", fill=(200, 202, 225), font=font(15))
    dr.text((22 * s, 40 * s), "nucleo al centro (non disegnato) · palline = elettrone", fill=(140, 145, 170), font=font(12.5))
    dr.text((wa + 16 * s, 14 * s), "la luce intorno  ·  scatola 416 λ_C ≈ 5 lunghezze d'onda", fill=(200, 202, 225), font=font(14))
    dr.text((wa + 16 * s, 34 * s), "l'atomo è il puntino al centro (vero rapporto di scala)", fill=(140, 145, 170), font=font(12))
    if "dz" not in SERIE:
        f = os.path.join(SCR, "el_dipolo.npy")
        if os.path.exists(f):
            SERIE["dz"] = list(np.load(f))
        else:
            SERIE["dz"] = [float(dat.el_pos(int(np.argmin(np.abs(TV - k * DT))))[2].mean()) for k in range(201)
                           if os.path.exists(os.path.join(SCR, f"el_{int(np.argmin(np.abs(TV - k * DT))) // CHUNK:03d}.npy"))]
    dzs = [(k * DT, z) for k, z in enumerate(SERIE["dz"]) if k * DT <= TV[iv] + 1e-9]
    dr.rectangle([wa, hl, W, H], fill=(8, 9, 16))
    pannello_info(dr, wa, hl, W - wa, H - hl, iv, dzs)
    return im


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--scala" in args:
        i = args.index("--scala")
        SCALA = float(args[i + 1])
        del args[i:i + 2]
    dat = Dati()
    if args[0] == "refs":   # fix the exposure once for the whole video (from frame IV_E for the atom, IV_L for the light)
        ive, ivl = int(args[1]), int(args[2])
        REFS.clear()
        scena_atomo(dat, ive, int(H0 * SCALA), int(H0 * SCALA))
        scena_luce(dat, ivl, int((W0 - H0) * SCALA), int(600 * SCALA))
        out = {"eb": REFS["eb_now"], "eg": REFS["eg_now"], "lb": REFS["lb_now"], "lg": REFS["lg_now"], "iv_e": ive, "iv_l": ivl}
        json.dump(out, open(REFS_FILE, "w"))
        print(out)
    elif args[0] == "still":
        iv = int(args[1])
        t0 = time.time()
        frame(dat, iv).save(args[2])
        print("still", iv, TV[iv], f"{time.time() - t0:.1f}s")
    elif args[0] == "batch":
        a, b = int(args[1]), int(args[2])
        d = os.path.join(SCR, f"frames_{int(W0 * SCALA)}")
        os.makedirs(d, exist_ok=True)
        t0 = time.time()
        budget = float(os.environ.get("VL_BUDGET", 480))
        for iv in range(a, min(b, NV)):
            f = os.path.join(d, f"{iv:05d}.png")
            if os.path.exists(f):
                continue
            if not os.path.exists(os.path.join(SCR, f"lu_{iv // CHUNK:03d}.npz")):
                continue   # light balls not moved this far yet
            frame(dat, iv).save(f + ".tmp.png")
            os.replace(f + ".tmp.png", f)
            if time.time() - t0 > budget:
                break
        print("batch", a, b, f"{time.time() - t0:.0f}s")

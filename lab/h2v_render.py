"""h2v_render.py - VIDEO H2: draw the frames written by lab/h2v_sim.py. Same recipe as lab/vz_render.py / vl_render.py:
every pixel of the scene is a ball (soft additive gaussian glow, depth weight, wide halo). No lines, trails, sticks or
surfaces. Text is UI.

    python3 lab/h2v_render.py still F out.png [--scala 1.0] [--solo-protoni]
    python3 lab/h2v_render.py refs                       # fix the exposure once (frames 0 and 300)
    python3 lab/h2v_render.py batch F0 F1 [--scala 0.6667]   # -> <SCR>/h2v_png_<H>/, skips existing

Visual knobs (aesthetic only, counted in docs/VIDEO_H2.md): SPOT, GLOW, GAIN, GLOW_MIX, DEPTH, palette, camera.
"""
import json, os, sys, time
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter, zoom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h2v_sim import SCR, FR, NA, NB, NC, NF, FPS  # noqa: E402

LAB = os.path.dirname(os.path.abspath(__file__))
W0, H0 = 1920, 1080
SCALA = 1.0
FONT = "/usr/share/fonts/opentype/inter/Inter-Regular.otf"
FONTB = "/usr/share/fonts/opentype/inter/Inter-SemiBold.otf"
SPOT = 1.25          # px at 1080p, sharp ball blur (crisper after Valerio's note of 13 Oct)
NT = 5               # tail: past frames re-splatted (5 frames = 0.21 s at 24 fps)
NS = 10              # sub-steps between consecutive frames (one smear, no chain of dots)
TAIL_W = 0.5         # tail weight at age 0
TAIL_TAU = 1.4       # tail decay (frames)
GLOW = 16.0          # px at 1080p, wide halo
GAIN = float(os.environ.get("H2V_GAIN", 1.9))
GLOW_MIX = 0.42
GAMMA = float(os.environ.get("H2V_GAMMA", 0.65))
DEPTH = 6.0          # a0
CAM_D = 7.5          # a0
FOV = float(np.degrees(0.8))   # camera.rs: 0.8 rad
C_E = np.array([0.36, 0.62, 1.00])        # electron balls
C_E_DEEP = np.array([0.20, 0.38, 1.00])   # sparse electron balls
C_E_HOT = np.array([0.70, 0.95, 1.00])    # dense electron balls
C_O = np.array([1.00, 0.60, 0.16])        # the other electron, given the reference one
C_REF = np.array([0.55, 1.00, 1.00])      # reference electron (moved by us along a chosen loop)
C_P = np.array([1.00, 0.80, 0.62])        # protons
REFS_F = os.path.join(SCR, "h2v_refs.json")
REFS = json.load(open(REFS_F)) if os.path.exists(REFS_F) else {}
ESAMI = json.load(open(os.path.join(LAB, "out", "h2v_esami.json"))) if os.path.exists(os.path.join(LAB, "out", "h2v_esami.json")) else {}
SOLO_P = False
import h2v_zzlook as ZZ  # noqa: E402
HUE_S = 8.0          # s of video per turn of the colour (phase), declared
HUE0 = 3.6           # starting phase: the blue-cyan end of the palette (choice)
P_ALPHA = 2.5        # proton glow brightness (choice)
R_ALPHA = 2.0        # reference electron brightness (choice)
ESPOSIZIONE = float(os.environ.get("H2V_ESPO", 1.0))
G_COL = 0.3          # colour contrast of g (aesthetic): g = 1.3 is fully orange
G_POW = float(os.environ.get("H2V_GPOW", 1.0))   # light of a ball = g ** G_POW (1 = exact conditional density)
_GM = []


def gmap():
    if not _GM:
        import h2v_core as C
        _GM.append(C.Gmap(C.WF(1.4), os.path.join(SCR, "h2v_rho1_R1.40.npz")))
    return _GM[0]


def camera(eye, centre, fov, w, h):
    fwd = np.asarray(centre, float) - eye
    fwd /= np.linalg.norm(fwd)
    right = np.cross(fwd, np.array([0.0, 0.0, 1.0])); right /= np.linalg.norm(right)
    up = np.cross(right, fwd)
    f = 0.5 * h / np.tan(np.radians(fov) / 2)
    return eye, fwd, right, up, f, w, h


def project(p, cam):
    eye, fwd, right, up, f, w, h = cam
    d = p - eye[:, None]
    z = (d * fwd[:, None]).sum(0); x = (d * right[:, None]).sum(0); y = (d * up[:, None]).sum(0)
    zz = np.maximum(z, 1e-3)
    return w / 2 + f * x / zz, h / 2 - f * y / zz, z


def splat(img, sx, sy, wts, col, w, h):
    ix = np.round(sx).astype(np.int64); iy = np.round(sy).astype(np.int64)
    ok = (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h) & (wts > 0)
    lin = iy[ok] * w + ix[ok]
    for c in range(3):
        img[c] += np.bincount(lin, wts[ok] * col[ok, c], w * h).reshape(h, w)
    img[3] += np.bincount(lin, wts[ok], w * h).reshape(h, w)


def blur(a, sig):
    if sig < 0.3:
        return a
    f = 4 if sig >= 8 else (2 if sig >= 4 else 1)
    if f == 1:
        return gaussian_filter(a, sig, truncate=3.0)
    h, w = a.shape
    H2, W2 = -(-h // f) * f, -(-w // f) * f
    p = np.zeros((H2, W2)); p[:h, :w] = a
    small = gaussian_filter(p.reshape(H2 // f, f, W2 // f, f).mean((1, 3)), sig / f, truncate=3.0)
    return zoom(small, f, order=1)[:h, :w]


def expose(img, gain, ref):
    return 1 - np.exp(-gain * np.maximum(img[:3] / ref, 0) ** GAMMA)


def cam_for(f, w, h):
    """camera path (aesthetic, declared): world X = molecular axis"""
    if f < NA:
        u = f / NA
        az, el = 90 - 55 * u, 14 - 4 * u
    elif f < NA + NB:
        u = (f - NA) / NB
        az, el = 35 + 55 * min(1, u * 4), 10 - 4 * min(1, u * 4)
    else:
        u = (f - NA - NB) / NC
        az, el = 90 + 20 * np.sin(np.pi * u), 6 + 6 * np.sin(np.pi * u)
    if az is None:
        pass
    a, e = np.radians(az), np.radians(el)
    eye = CAM_D * np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    return camera(eye, (0, 0, 0), FOV, w, h)


def world(p):
    """physics (x, y, z) with the bond on z -> world (z, x, y)"""
    return np.stack([p[:, 2], p[:, 0], p[:, 1]])


def carica(f):
    """data loader (separable from drawing): one frame of balls + the past NT frames of the same walkers (for the tail)"""
    def one(k):
        z = np.load(os.path.join(FR, f"f_{k:04d}.npz"))
        return z["p"].astype(float), z["c"], z["prot"], json.loads(str(z["meta"]))
    P, cl, prot, meta = one(f)
    past = []
    for k in range(1, NT + 1):
        if f - k < 0:
            break
        Pk, ck, pk, mk = one(f - k)
        if (ck == 0).sum() != (cl == 0).sum():
            break                       # different walker set (part C v2 starts fresh): no tail across it
        past.append((Pk[ck == 0], pk))
    return dict(P=P, cl=cl, prot=prot, meta=meta, past=past, f=f)


def tail_points(cur, past, interp=True):
    """the short fading tail (Valerio's rule of 13 Oct: 0.2-0.4 s): the same glow re-splatted along the last NT frames,
    NS sub-steps between consecutive frames so it reads as one smear; weight decays with age. Returns points, weights"""
    pts, wts, ages = [], [], []
    prev = cur
    ns = NS if interp else 1
    for k, Pk in enumerate(past):
        for j in range(1, ns + 1):
            u = j / ns
            age = k + u
            pts.append(prev + (Pk - prev) * u)
            wts.append(np.full(len(Pk), TAIL_W * np.exp(-age / TAIL_TAU) / ns))
            ages.append(age)
        prev = Pk
    return pts, wts, ages


def disegna(d, w, h, cam=None):
    """draw function: data from carica() -> (3, h, w) image in [0, 1], zig-zag atom look (lab/h2v_zzlook.py)"""
    P, cl, prot, meta = d["P"], d["cl"], d["prot"], d["meta"]
    cam = cam or cam_for(d["f"], w, h)
    segB = meta["seg"] == "B"
    f_px = 0.5 * h / np.tan(ZZ.FOV / 2)
    size = ZZ.size_world(CAM_D)
    e0 = cl == 0
    Pe = P[e0]
    # brightness from the local density (zzscena.rs rule), density from a 3-D histogram of the balls (0.25 a0 cells)
    edg = np.linspace(-5, 5, 41)
    hh, _ = np.histogramdd(Pe, [edg] * 3)
    ii = np.clip(((Pe + 5) / 0.25).astype(int), 0, 39)
    a = ZZ.alpha(hh[ii[:, 0], ii[:, 1], ii[:, 2]] / max(hh.max(), 1)) * (1.5 if meta.get("version") == "C2" else 1.0)
    # colour = phase of the state: nodeless real ground state -> one phase everywhere, turning at its energy frequency
    # (shown slowed: one turn every HUE_S seconds of video, declared)
    ph = 2 * np.pi * d["f"] / (FPS * HUE_S) + HUE0
    col = np.tile(ZZ.tinta(np.array([ph]))[0], (len(Pe), 1))
    if segB:
        g = gmap().g(np.array(meta["ref"]), Pe)
        meta["g_mean"] = float(g.mean())
        kk = np.clip((g - 1) / G_COL, -1, 1)[:, None]
        mid = 0.5 * (ZZ.STOPS[1] + ZZ.STOPS[4])
        col = np.where(kk > 0, (1 - kk) * mid + kk * ZZ.STOPS[4], (1 + kk) * mid + (-kk) * ZZ.STOPS[1])
        a = a * g ** G_POW
    T = ZZ.Tela(w, h)
    if not SOLO_P:
        sx, sy, zc = project(world(Pe), cam)
        T.balls(sx, sy, size * f_px / np.maximum(zc, 0.3), col, a)
        for pt, tw, age in zip(*tail_points(Pe, [p for p, _ in d["past"]], meta.get("version") != "C2")):
            tx, ty, tz = project(world(pt), cam)
            shrink = max(0.4, 1 - 0.12 * age)          # older glow: dimmer and smaller
            T.balls(tx, ty, shrink * size * f_px / np.maximum(tz, 0.3), col, a * tw)
    # protons: one bright white glow each (4x the sprite), with their own short tail when they move
    px, py, pz = project(world(prot), cam)
    T.balls(px, py, 4 * size * f_px / pz, np.ones((2, 3)), np.full(2, P_ALPHA))
    for pt, tw, age in zip(*tail_points(prot, [q for _, q in d["past"]])):
        qx, qy, qz = project(world(pt), cam)
        T.balls(qx, qy, 4 * size * f_px / qz, np.ones((2, 3)), P_ALPHA * tw)
    if segB and not SOLO_P:
        r1 = np.array(meta["ref"])[None]
        rx, ry, rz = project(world(r1), cam)
        T.balls(rx, ry, 4 * size * f_px / rz, ZZ.STOPS[3][None], np.array([R_ALPHA]))
    img = T.hdr()
    out = ZZ.finisci(img, h / H0, ESPOSIZIONE).transpose(2, 0, 1) / 255.0
    raw = dict(b=float(np.percentile(img.sum(0)[img.sum(0) > 0], 99.7)) if (img > 0).any() else 1.0, g=0.0)
    return out, meta, cam, raw


def scena(f, w, h, cam=None):
    return disegna(carica(f), w, h, cam)


def font(sz, bold=False):
    return ImageFont.truetype(FONTB if bold else FONT, max(8, int(round(sz * SCALA))))


def num(x, fmt="{:.2f}"):
    return fmt.format(x).replace(".", ",")


def testi(meta):
    seg = meta["seg"]
    V4 = ESAMI.get("V4", {})
    popp = V4.get("P_opp_analytic", float("nan"))
    if seg == "A":
        return ("H₂, la prima molecola: due protoni, due elettroni",
                ["24 000 palline = posizioni estratte dallo stato a due elettroni |Ψ(r₁, r₂)|² calcolato",
                 "legame calcolato: R_e = 1,4009 a₀ · energia di legame D_e = 4,748 eV (noti: 1,401 a₀ e 4,747 eV)"])
    if seg == "B":
        return ("Dove sta l'altro elettrone? (la correlazione)",
                ["pallina rosa: un elettrone messo da noi qui (giro scelto) · luce delle palline = dove sta allora l'altro (calcolato)",
                 "giallo: l'altro elettrone ci sta più del solito · blu: ci sta meno (scala dei colori scelta)",
                 f"si evitano: lati opposti nel {num(100 * popp, '{:.1f}')} % dei casi (senza correlazione sarebbe il 50 %)"])
    return ("Il legame vibra",
            [f"stato scelto: metà v = 0, metà v = 1 (livelli calcolati) · periodo 8,0 fs mostrato in 4 s",
             f"t = {num(meta['t_fs'], '{:.2f}')} fs   R medio = {num(meta['R'], '{:.3f}')} a₀   la nuvola segue i protoni (Born–Oppenheimer)"])


ONESTA = ["Stato a due elettroni dal risolutore ECG di Debiti (gaussiane correlate esplicite), non ancora palline zig-zag sulla griglia.",
          "Protoni fermi mentre si calcolano gli elettroni (Born–Oppenheimer). Il moto delle palline è campionamento, non traiettoria.",
          "Aspetto dell'atomo zig-zag: bagliori additivi, colore = fase dello stato (senza nodi: una sola, gira rallentata), scia corta, tono approssimato."]


def ui(out, meta, w, h):
    im = Image.fromarray((np.clip(out, 0, 1).transpose(1, 2, 0) * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im)
    s = h / H0
    tit, righe = testi(meta)
    dr.text((40 * s, 30 * s), tit, fill=(245, 245, 252), font=font(34, True))
    for i, r in enumerate(righe):
        dr.text((40 * s, 80 * s + 30 * i * s), r, fill=(200, 205, 225), font=font(19))
    for i, r in enumerate(ONESTA):
        dr.text((40 * s, h - 92 * s + 26 * i * s), r, fill=(140, 145, 170), font=font(15))
    return im


def frame(f, w, h):
    out, meta, cam, raw = scena(f, w, h)
    return ui(out, meta, w, h), meta, raw


def job_refs():
    rs = []
    for f in (0, NA + 100, NA + NB + 50):
        out, meta, cam, raw = scena(f, 1920, 1080)
        rs.append(raw)
    REFS.update(b=float(np.median([r["b"] for r in rs])), g=float(np.median([r["g"] for r in rs])))
    json.dump(REFS, open(REFS_F, "w"))
    print(REFS)


def args():
    global SCALA, SOLO_P
    a = sys.argv
    if "--scala" in a:
        SCALA = float(a[a.index("--scala") + 1])
    SOLO_P = "--solo-protoni" in a


if __name__ == "__main__":
    args()
    W, H = int(round(W0 * SCALA / 2)) * 2, int(round(H0 * SCALA / 2)) * 2
    if sys.argv[1] == "refs":
        job_refs()
    elif sys.argv[1] == "still":
        im, meta, raw = frame(int(sys.argv[2]), W, H)
        im.save(sys.argv[3]); print(meta.get("seg"), raw)
    elif sys.argv[1] == "batch":
        d = os.path.join(SCR, f"h2v_png_{H}"); os.makedirs(d, exist_ok=True)
        t0 = time.time()
        for f in range(int(sys.argv[2]), int(sys.argv[3])):
            fn = os.path.join(d, f"{f:04d}.jpg")
            if os.path.exists(fn):
                continue
            if time.time() - t0 > 540:
                print("checkpoint: run again", flush=True); break
            frame(f, W, H)[0].save(fn, quality=94)
        print("batch done", flush=True)

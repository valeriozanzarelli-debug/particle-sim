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
SPOT = 1.6           # px at 1080p, sharp ball blur
GLOW = 16.0          # px at 1080p, wide halo
GAIN = float(os.environ.get("H2V_GAIN", 1.3))
GLOW_MIX = 0.55
GAMMA = float(os.environ.get("H2V_GAMMA", 0.65))
DEPTH = 6.0          # a0
CAM_D = 7.5          # a0
FOV = 40.0
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


def scena(f, w, h, cam=None):
    s = h / H0
    z = np.load(os.path.join(FR, f"f_{f:04d}.npz"))
    P, cl, prot = z["p"].astype(float), z["c"], z["prot"]
    meta = json.loads(str(z["meta"]))
    cam = cam or cam_for(f, w, h)
    segB = meta["seg"] == "B"
    # colour by local ball density (3-D histogram of the balls themselves, 0.25 a0 cells): dense -> white-cyan (aesthetic)
    e0 = cl == 0
    edg = np.linspace(-5, 5, 41)
    hh, _ = np.histogramdd(P[e0], [edg] * 3)
    ii = np.clip(((P + 5) / 0.25).astype(int), 0, 39)
    loc = hh[ii[:, 0], ii[:, 1], ii[:, 2]] / max(hh.max(), 1)
    k = np.clip(loc ** 0.5, 0, 1)[:, None]
    ce = (1 - k) * C_E_DEEP + k * C_E_HOT
    col = np.where(cl[:, None] == 0, ce, np.where(cl[:, None] == 1, C_O, C_REF))
    amp = np.where(cl == 0, 1.0, 0.0)          # class 1 (orange walkers of the exam V4b) not drawn
    if segB:
        # the same electron balls; light weighted with g(r2 | r_ref): total light = where the other electron is
        g = gmap().g(np.array(meta["ref"]), P[e0])
        meta["g_mean"] = float(g.mean())
        kk = np.clip((g - 1) / G_COL, -1, 1)[:, None]
        cB = np.where(kk > 0, (1 - kk) * C_E + kk * C_O, (1 + kk) * C_E + (-kk) * C_E_DEEP * 0.6)
        col[e0] = cB
        amp[e0] = g ** G_POW
    sx, sy, zc = project(world(P), cam)
    dw = np.exp(-np.clip(zc - CAM_D, 0, None) / DEPTH) * np.exp(np.clip(CAM_D - zc, 0, None) / (2 * DEPTH))
    wt = amp * dw
    sharp = np.zeros((4, h, w)); wide = np.zeros((4, h, w))
    if not SOLO_P:
        splat(sharp, sx, sy, wt, col, w, h)
        splat(wide, sx, sy, wt, col, w, h)
    # protons: one bright glow each (many weight at one point, then the same blur chain) + the reference electron
    big = np.zeros((4, h, w))
    px, py, pz = project(world(prot), cam)
    splat(big, px, py, np.full(2, 900.0) * (CAM_D / pz) ** 2, np.tile(C_P, (2, 1)), w, h)
    if segB and not SOLO_P:
        i = np.where(cl == 2)[0]
        splat(big, sx[i], sy[i], np.full(len(i), 500.0), col[i], w, h)
    img = np.stack([blur(sharp[c], SPOT * s) + blur(sharp[c], 2.6 * SPOT * s) * 0.5 for c in range(4)])
    pimg = np.stack([blur(big[c], 4.5 * s) + 0.35 * blur(big[c], 12 * s) for c in range(4)])
    ref_b = REFS.get("b", 1.0) * s ** -2
    ref_g = REFS.get("g", 1.0) * s ** -2
    balls = expose(img + pimg, GAIN, ref_b)
    glow = expose(np.stack([blur(wide[c] + 0.05 * big[c], GLOW * s) for c in range(4)]), 1.0, ref_g)
    out = 1 - (1 - np.clip(balls, 0, 1)) * (1 - GLOW_MIX * np.clip(glow, 0, 1))
    raw = dict(b=float(np.percentile(img[3][img[3] > 0], 99.7)) * s ** 2 if (img[3] > 0).any() else 1.0,
               g=float(np.percentile(blur(wide[3], GLOW * s)[wide[3] >= 0], 99.9)) * s ** 2)
    return out, meta, cam, raw


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
                ["pallina azzurra: un elettrone messo da noi qui (giro scelto) · luce delle palline = dove sta allora l'altro (calcolato)",
                 "arancio: l'altro elettrone ci sta più del solito · blu scuro: ci sta meno (scala dei colori scelta)",
                 f"si evitano: lati opposti nel {num(100 * popp, '{:.1f}')} % dei casi (senza correlazione sarebbe il 50 %)"])
    return ("Il legame vibra",
            [f"stato scelto: metà v = 0, metà v = 1 (livelli calcolati) · periodo 8,0 fs mostrato in 4 s",
             f"t = {num(meta['t_fs'], '{:.2f}')} fs   ⟨R⟩ = {num(meta['R'], '{:.3f}')} a₀   la nuvola segue i protoni (Born–Oppenheimer)"])


ONESTA = ["Stato a due elettroni dal risolutore ECG di Debiti (gaussiane correlate esplicite), non ancora palline zig-zag sulla griglia.",
          "Protoni fermi mentre si calcolano gli elettroni (Born–Oppenheimer). Il moto delle palline è campionamento, non traiettoria.",
          "Bagliori allargati e luminosità compressa: estetica."]


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
            fn = os.path.join(d, f"{f:04d}.png")
            if os.path.exists(fn):
                continue
            if time.time() - t0 > 540:
                print("checkpoint: run again", flush=True); break
            frame(f, W, H)[0].save(fn)
        print("batch done", flush=True)

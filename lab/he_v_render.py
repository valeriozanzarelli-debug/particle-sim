"""he_v_render.py - VIDEO ELIO: helium ground state as glows. Data from the main thread's export
(/mnt/project-files/ripresa/dati/elio/, see manifest.json): 200 000 electron pairs sampled from |psi(r1, r2)|^2,
the exact one-electron density, the conditional maps (electron 1 held on +z at 1.0 a0).

The state is stationary: every frame draws a FRESH batch of pairs (no trajectories). Tail = in-place persistence of the
last NT frames' glows, fading (Valerio's rule of 13 Oct). Separable: carica(f) -> data, disegna(data, w, h) -> image.

    python3 lab/he_v_render.py still F out.png [--scala 1.0]
    python3 lab/he_v_render.py refs
    python3 lab/he_v_render.py batch F0 F1 [--scala 0.6667]
"""
import json, os, sys, time
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h2v_render import camera, project, splat, blur, expose, FONT, FONTB  # noqa: E402  (same glow recipe)

DATA = os.environ.get("HEV_DATA", "/mnt/project-files/ripresa/dati/elio")
SCR = os.environ.get("HEV_SCR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "he_v_scr"))
os.makedirs(SCR, exist_ok=True)
LAB = os.path.dirname(os.path.abspath(__file__))
W0, H0 = 1920, 1080
SCALA = 1.0
FPS = 24
NA, NB = 192, 240
NF = NA + NB
NPAIR = 12000          # pairs drawn per frame
NT = 5                 # tail: persistence of the last 5 frames (0.21 s)
TAIL_W, TAIL_TAU = 0.5, 1.4
SPOT, GLOW, GAIN, GLOW_MIX, GAMMA = 1.25, 14.0, 1.9, 0.42, 0.65
DEPTH, CAM_D, FOV = 5.0, 6.0, 40.0
R_REF = 1.0            # a0: the reference electron sits on a chosen loop at this radius (the map is for |r1| = 1.0)
G_COL = 0.3
C_E_DEEP = np.array([0.20, 0.38, 1.00]); C_E_HOT = np.array([0.70, 0.95, 1.00]); C_E = np.array([0.36, 0.62, 1.00])
C_O = np.array([1.00, 0.60, 0.16]); C_REF = np.array([0.55, 1.00, 1.00]); C_N = np.array([1.00, 0.78, 0.55])
REFS_F = os.path.join(SCR, "he_v_refs.json")
REFS = json.load(open(REFS_F)) if os.path.exists(REFS_F) else {}
ES_F = os.path.join(LAB, "out", "he_v_esami.json")
ESAMI = json.load(open(ES_F)) if os.path.exists(ES_F) else {}
MAN = json.load(open(os.path.join(DATA, "manifest.json")))
_D = {}


def dati():
    if not _D:
        _D["pairs"] = np.load(os.path.join(DATA, "elio_coppie.npy")).astype(float)
        d = np.load(os.path.join(DATA, "elio_densita.npz"))
        g = d["griglia_xz"]
        # rapporto[iz, ix]: g(r2 | r1 = (0, 0, 1.0)) on the x-z plane (checked: minimum at z = 1.0, x = 0)
        _D["g_ip"] = RegularGridInterpolator((g, g), d["rapporto_r1_1.0"], bounds_error=False, fill_value=None)
        _D["dens"] = d
    return _D


def ref_path(u):
    """chosen loop of the reference electron at r = 1.0 a0 (declared)"""
    a = 2 * np.pi * u
    v = np.array([np.sin(a), 0.35 * np.sin(2 * a), np.cos(a)])
    return R_REF * v / np.linalg.norm(v)


def g_of(ref, P):
    """g(r2 | r1 = ref) from the map, using spherical symmetry: rotate so ref -> +z, map is axially symmetric"""
    u = ref / np.linalg.norm(ref)
    z = P @ u
    rho = np.sqrt(np.maximum(np.sum(P * P, 1) - z * z, 0))
    q = np.stack([np.clip(z, -3, 3), np.clip(rho, 0, 3)], 1)   # outside +-3 a0: edge value (declared)
    return dati()["g_ip"](q)


def batch(f):
    rng = np.random.default_rng(1000 + f)
    i = rng.integers(0, len(dati()["pairs"]), NPAIR)
    return dati()["pairs"][i].reshape(-1, 3)


def carica(f):
    """data loader: fresh batch of pairs for frame f + the batches of the last NT frames (tail)"""
    meta = dict(f=f, seg="A" if f < NA else "B")
    if f >= NA:
        meta["ref"] = ref_path((f - NA) / NB).tolist()
    return dict(f=f, P=batch(f), past=[batch(f - k) for k in range(1, NT + 1) if f - k >= 0], meta=meta)


def cam_for(f, w, h):
    if f < NA:
        u = f / NA; az, el = -30 + 70 * u, 18 - 6 * u
    else:
        u = (f - NA) / NB; az, el = 90, 4      # side view of the loop plane (x-z)
    a, e = np.radians(az), np.radians(el)
    eye = CAM_D * np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    return camera(eye, (0, 0, 0), FOV, w, h)


def colori(P):
    edg = np.linspace(-4, 4, 41)
    hh, _ = np.histogramdd(P, [edg] * 3)
    ii = np.clip(((P + 4) / 0.2).astype(int), 0, 39)
    k = np.clip((hh[ii[:, 0], ii[:, 1], ii[:, 2]] / max(hh.max(), 1)) ** 0.5, 0, 1)[:, None]
    return (1 - k) * C_E_DEEP + k * C_E_HOT


def luce(d, P):
    """colour and light of a set of balls (B: light = g, the exact conditional density over balls that follow rho1)"""
    col = colori(P); amp = np.ones(len(P))
    if d["meta"]["seg"] == "B":
        g = g_of(np.array(d["meta"]["ref"]), P)
        kk = np.clip((g - 1) / G_COL, -1, 1)[:, None]
        col = np.where(kk > 0, (1 - kk) * C_E + kk * C_O, (1 + kk) * C_E + (-kk) * C_E_DEEP * 0.6)
        amp = g
        return col, amp, g
    return col, amp, None


def disegna(d, w, h, cam=None):
    s = h / H0
    cam = cam or cam_for(d["f"], w, h)
    P = d["P"]
    col, amp, g = luce(d, P)
    if g is not None:
        d["meta"]["g_mean"] = float(g.mean())
    sx, sy, zc = project(P.T, cam)
    dw = np.exp(-np.clip(zc - CAM_D, 0, None) / DEPTH) * np.exp(np.clip(CAM_D - zc, 0, None) / (2 * DEPTH))
    sharp = np.zeros((4, h, w)); tail = np.zeros((4, h, w)); big = np.zeros((4, h, w))
    splat(sharp, sx, sy, amp * dw, col, w, h)
    for k, Pk in enumerate(d["past"]):
        ck, ak, _ = luce(d, Pk)
        tx, ty, tz = project(Pk.T, cam)
        tw = TAIL_W * np.exp(-(k + 1) / TAIL_TAU)
        splat(tail, tx, ty, tw * ak * np.exp(-np.clip(tz - CAM_D, 0, None) / DEPTH), ck, w, h)
    nx, ny, nz = project(np.zeros((3, 1)), cam)
    splat(big, nx, ny, np.array([1600.0]), C_N[None], w, h)          # the nucleus (He, Z = 2), fixed
    if d["meta"]["seg"] == "B":
        r = np.array(d["meta"]["ref"])[:, None]
        rx, ry, rz = project(r, cam)
        splat(big, rx, ry, np.array([500.0]), C_REF[None], w, h)
    img = np.stack([blur(sharp[c], SPOT * s) + 0.35 * blur(sharp[c], 2.6 * SPOT * s) + blur(tail[c], 2.2 * SPOT * s)
                    for c in range(4)])
    pimg = np.stack([blur(big[c], 4.0 * s) + 0.3 * blur(big[c], 12 * s) for c in range(4)])
    wide = sharp + tail
    balls = expose(img + pimg, GAIN, REFS.get("b", 1.0) * s ** -2)
    glow = expose(np.stack([blur(wide[c] + 0.05 * big[c], GLOW * s) for c in range(4)]), 1.0, REFS.get("g", 1.0) * s ** -2)
    out = 1 - (1 - np.clip(balls, 0, 1)) * (1 - GLOW_MIX * np.clip(glow, 0, 1))
    raw = dict(b=float(np.percentile(img[3][img[3] > 0], 99.7)) * s ** 2, g=float(np.percentile(blur(wide[3], GLOW * s), 99.9)) * s ** 2)
    return out, raw


def font(sz, bold=False):
    return ImageFont.truetype(FONTB if bold else FONT, max(8, int(round(sz * SCALA))))


def num(x, fmt="{:.3f}"):
    return fmt.format(x).replace(".", ",")


ONESTA = ["Stato dell'elio dal thread principale: gaussiane correlate (ECG, 200 funzioni) con la regola a molti elettroni, non ancora palline zig-zag.",
          "Nucleo fermo. Lo stato è fermo: ogni fotogramma pesca coppie nuove (il tremolio non è moto). Scia = ultimi 5 fotogrammi che si spengono.",
          "Bagliori allargati e luminosità compressa: estetica."]


def testi(meta):
    E2 = ESAMI.get("E2", {}); E3 = ESAMI.get("E3", {})
    if meta["seg"] == "A":
        return ("Elio: un nucleo, due elettroni",
                [f"24 000 palline per fotogramma = 12 000 coppie (elettrone 1, elettrone 2) pescate dallo stato calcolato",
                 f"energia {num(MAN['sorgente']['E_Eh'], '{:.6f}').replace('-', '−')} Eh (esatta −2,903724) · raggio medio {num(E2.get('mean_r', float('nan')))} a₀ (noto 0,9295)"])
    return ("Dove sta un elettrone, l'altro si scansa (buca di Coulomb)",
            ["pallina azzurra: un elettrone messo da noi a 1 a₀ dal nucleo (giro scelto) · luce delle palline = dove sta allora l'altro",
             "arancio: l'altro ci sta più del solito · blu scuro: meno (scala dei colori scelta)",
             f"vicini (meno di 0,5 a₀): {num(100 * E3.get('P_real', float('nan')), '{:.1f}')} % delle coppie vere contro "
             f"{num(100 * E3.get('P_shuffled', float('nan')), '{:.1f}')} % se fossero indipendenti"])


def ui(out, meta, w, h):
    im = Image.fromarray((np.clip(out, 0, 1).transpose(1, 2, 0) * 255).astype(np.uint8))
    dr = ImageDraw.Draw(im); s = h / H0
    tit, righe = testi(meta)
    dr.text((40 * s, 30 * s), tit, fill=(245, 245, 252), font=font(34, True))
    for i, r in enumerate(righe):
        dr.text((40 * s, 80 * s + 30 * i * s), r, fill=(200, 205, 225), font=font(19))
    for i, r in enumerate(ONESTA):
        dr.text((40 * s, h - 92 * s + 26 * i * s), r, fill=(140, 145, 170), font=font(15))
    return im


def frame(f, w, h):
    d = carica(f)
    out, raw = disegna(d, w, h)
    return ui(out, d["meta"], w, h), d["meta"], raw


if __name__ == "__main__":
    a = sys.argv
    if "--scala" in a:
        SCALA = float(a[a.index("--scala") + 1])
    W, H = int(round(W0 * SCALA / 2)) * 2, int(round(H0 * SCALA / 2)) * 2
    if a[1] == "refs":
        rs = [disegna(carica(f), 1920, 1080)[1] for f in (10, NA + 60)]
        REFS.update(b=float(np.median([r["b"] for r in rs])), g=float(np.median([r["g"] for r in rs])))
        json.dump(REFS, open(REFS_F, "w")); print(REFS)
    elif a[1] == "still":
        im, meta, raw = frame(int(a[2]), W, H); im.save(a[3]); print(meta, raw)
    elif a[1] == "batch":
        dd = os.path.join(SCR, f"he_v_png_{H}"); os.makedirs(dd, exist_ok=True)
        t0 = time.time()
        for f in range(int(a[2]), int(a[3])):
            fn = os.path.join(dd, f"{f:04d}.png")
            if os.path.exists(fn):
                continue
            if time.time() - t0 > 540:
                print("checkpoint: run again", flush=True); break
            frame(f, W, H)[0].save(fn)
        print("batch done", flush=True)

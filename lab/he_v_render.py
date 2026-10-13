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
from h2v_render import camera, project, FONT, FONTB  # noqa: E402
import h2v_zzlook as ZZ  # noqa: E402  (the zig-zag atom look, parameters copied from src/atomo)

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
CAM_D, FOV = 6.0, float(np.degrees(0.8))
HUE_S, HUE0 = 8.0, 3.6      # colour turn: one per 8 s of video (slowed, declared); start at the blue-cyan end
N_ALPHA, R_ALPHA = 3.0, 2.0  # nucleus and reference electron brightness (choice)
ESPOSIZIONE = float(os.environ.get("HEV_ESPO", 1.0))
R_REF = 1.0            # a0: the reference electron sits on a chosen loop at this radius (the map is for |r1| = 1.0)
G_COL = 0.15          # colour contrast of g (choice): g = 1.15 fully gold, 0.85 fully blue
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


def luce(d, P):
    """colour and alpha of a set of balls, zig-zag atom rule (h2v_zzlook): alpha from the EXACT one-electron density at the
    ball (rho1 of the export); colour = phase of the state (nodeless ground state: one phase, turning, shown slowed).
    B: colour from g (blue = less, gold = more, as in the main thread's still) and alpha times g"""
    dd = dati()["dens"]
    r = np.linalg.norm(P, axis=1)
    rho = np.interp(r, dd["rho1_esatta_r"], dd["rho1_esatta"], right=0.0)
    a = ZZ.alpha(rho / dd["rho1_esatta"].max())
    ph = 2 * np.pi * d["f"] / (FPS * HUE_S) + HUE0
    col = np.tile(ZZ.tinta(np.array([ph]))[0], (len(P), 1))
    if d["meta"]["seg"] == "B":
        g = g_of(np.array(d["meta"]["ref"]), P)
        kk = np.clip((g - 1) / G_COL, -1, 1)[:, None]
        mid = 0.5 * (ZZ.STOPS[1] + ZZ.STOPS[4])
        col = np.where(kk > 0, (1 - kk) * mid + kk * ZZ.STOPS[4], (1 + kk) * mid + (-kk) * ZZ.STOPS[1])
        return col, a * g, g
    return col, a, None


def disegna(d, w, h, cam=None):
    """draw function: data from carica() -> (3, h, w) in [0, 1]"""
    cam = cam or cam_for(d["f"], w, h)
    f_px = 0.5 * h / np.tan(ZZ.FOV / 2)
    size = ZZ.size_world(CAM_D)
    P = d["P"]
    col, a, g = luce(d, P)
    if g is not None:
        d["meta"]["g_mean"] = float(g.mean())
    T = ZZ.Tela(w, h)
    sx, sy, zc = project(P.T, cam)
    T.balls(sx, sy, size * f_px / np.maximum(zc, 0.3), col, a)
    # tail: the glows of the last NT frames stay where they were and fade/shrink (fresh pairs: no path to smear)
    for k, Pk in enumerate(d["past"]):
        ck, ak, _ = luce(d, Pk)
        tx, ty, tz = project(Pk.T, cam)
        age = k + 1
        T.balls(tx, ty, max(0.4, 1 - 0.12 * age) * size * f_px / np.maximum(tz, 0.3), ck, ak * TAIL_W * np.exp(-age / TAIL_TAU))
    nx, ny, nz = project(np.zeros((3, 1)), cam)
    T.balls(nx, ny, 4 * size * f_px / nz, np.ones((1, 3)), np.array([N_ALPHA]))       # the nucleus, fixed
    if d["meta"]["seg"] == "B":
        rx, ry, rz = project(np.array(d["meta"]["ref"])[:, None], cam)
        T.balls(rx, ry, 4 * size * f_px / rz, ZZ.STOPS[3][None], np.array([R_ALPHA]))
    img = T.hdr()
    out = ZZ.finisci(img, h / H0, ESPOSIZIONE).transpose(2, 0, 1) / 255.0
    return out, dict(b=float(np.percentile(img.sum(0)[img.sum(0) > 0], 99.7)), g=0.0)


def font(sz, bold=False):
    return ImageFont.truetype(FONTB if bold else FONT, max(8, int(round(sz * SCALA))))


def num(x, fmt="{:.3f}"):
    return fmt.format(x).replace(".", ",")


ONESTA = ["Stato dell'elio dal thread principale: gaussiane correlate (ECG, 200 funzioni) con la regola a molti elettroni, non ancora palline zig-zag.",
          "Nucleo fermo. Lo stato è fermo: ogni fotogramma pesca coppie nuove (il tremolio non è moto). Scia = ultimi 5 fotogrammi che si spengono.",
          "Aspetto dell'atomo zig-zag: bagliori additivi, colore = fase dello stato (senza nodi: una sola, gira rallentata), tono approssimato."]


def testi(meta):
    E2 = ESAMI.get("E2", {}); E3 = ESAMI.get("E3", {})
    if meta["seg"] == "A":
        return ("Elio: un nucleo, due elettroni",
                [f"24 000 palline per fotogramma = 12 000 coppie (elettrone 1, elettrone 2) pescate dallo stato calcolato",
                 f"energia {num(MAN['sorgente']['E_Eh'], '{:.6f}').replace('-', '−')} Eh (esatta −2,903724) · raggio medio {num(E2.get('mean_r', float('nan')))} a₀ (noto 0,9295)"])
    return ("Dove sta un elettrone, l'altro si scansa (buca di Coulomb)",
            ["pallina rosa: un elettrone messo da noi a 1 a₀ dal nucleo (giro scelto) · luce delle palline = dove sta allora l'altro",
             "giallo: l'altro ci sta più del solito · blu: meno (scala dei colori scelta)",
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
            fn = os.path.join(dd, f"{f:04d}.jpg")
            if os.path.exists(fn):
                continue
            if time.time() - t0 > 540:
                print("checkpoint: run again", flush=True); break
            frame(f, W, H)[0].save(fn, quality=94)
        print("batch done", flush=True)

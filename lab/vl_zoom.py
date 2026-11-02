"""vl_zoom.py - the close-up Valerio asked to open (13/10 sera): the electron seen from above (looking down z),
zoomed on the nucleus, in slow motion, only soft additive glows (no trails, no streaks, nothing drawn by hand).

Balls: born where |psi|^2 is (inverse sampling of the exported rho, uniform inside the cell), carried by the rule's
current v = j/rho (exported, averaged over each 16-tick frame, linear in time between frames, Heun steps of
0.5 lambda_C/c), each living LIFE lambda_C/c and fading in and out (sin of its age), births staggered. So the cloud
always follows |psi|^2 (no drift of tracers away from the core) and the motion is the current's.
Hue = phase of the wave at the ball (as the Lab's zig-zag scene); glow size from the local density; a gaussian slab
around z = 0 so the turning in the x-y plane is not hidden by the balls above and below.

    python3 lab/vl_zoom.py palline              # ball states for every video frame -> <SCR>/zoom_palline.npy
    python3 lab/vl_zoom.py batch IV0 IV1 [--scala 1.0]
"""
import json
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vl_render as R  # noqa: E402
from vl_comune import SCR, A, NA, trilinear  # noqa: E402
from atomo_export_leggi import carica  # noqa: E402

DIR = "/mnt/project-files/ripresa/dati/atomo_accoppiato"
T0, T1, DTV = 300.0, 1500.0, 2.0     # lambda_C/c: clip window and time per video frame (slow motion)
N = 120000
LIFE = 24.0                           # lambda_C/c (3 data frames)
SUB = 0.5                             # Heun step, lambda_C/c
HALF = 5.0                            # lambda_C: half height of the view
SLAB = 1.6                            # lambda_C: gaussian half-thickness of the slab around z = 0
SPOT, GLOW, GAIN, GLOW_MIX = 2.4, 14.0, 2.2, 0.5
TT = np.arange(T0, T1 + 1e-9, DTV)
PAL = os.path.join(SCR, "zoom_palline.npy")


class Campi:
    def __init__(self):
        self.d = carica(DIR)
        self.c = {}

    def fr(self, k):
        if k not in self.c:
            if len(self.c) > 4:
                self.c.pop(min(self.c))
            r = np.maximum(self.d.rho(k).astype(np.float64), 0)
            self.c[k] = (r, self.d.j(k).astype(np.float64))
        return self.c[k]

    def rj(self, t):
        k = int(min(np.floor(t / 8.0), 199))
        u = t / 8.0 - k
        (r0, j0), (r1, j1) = self.fr(k), self.fr(k + 1)
        return (1 - u) * r0 + u * r1, (1 - u) * j0 + u * j1

    def vel(self, p, t):
        r, j = self.rj(t)
        idx = p / A + NA / 2 - 0.5
        rr = np.maximum(trilinear(r[None], idx)[0], 1e-6 * r.max())
        v = trilinear(j, idx) / rr
        s = np.linalg.norm(v, axis=0)
        return v * np.minimum(1.0, 0.98 / np.maximum(s, 1e-12))   # |v| < c (the rule's own bound)

    def nasci(self, n, t, rng):
        r, _ = self.rj(t)
        c = np.cumsum(r.ravel())
        q = np.searchsorted(c, rng.random(n) * c[-1])
        i, jj, k = np.unravel_index(q, r.shape)
        return np.stack([(i + rng.random(n)) * A - NA * A / 2, (jj + rng.random(n)) * A - NA * A / 2, (k + rng.random(n)) * A - NA * A / 2])


def palline():
    rng = np.random.default_rng(13)
    F = Campi()
    p = F.nasci(N, T0 - LIFE, rng)
    eta = rng.random(N) * LIFE
    t = T0 - LIFE
    out = np.zeros((len(TT), 4, N), np.float32)
    iv = 0
    while iv < len(TT):
        if t >= TT[iv] - 1e-9:
            out[iv, :3] = p
            out[iv, 3] = eta
            iv += 1
            continue
        h = min(SUB, TT[iv] - t)
        v1 = F.vel(p, t)
        v2 = F.vel(p + h * v1, t + h)
        p = p + 0.5 * h * (v1 + v2)
        t += h
        eta += h
        old = eta >= LIFE
        if old.any():
            p[:, old] = F.nasci(int(old.sum()), t, rng)
            eta[old] -= LIFE
    np.save(PAL, out)
    # check: the balls follow |psi|^2 (correlation of the ball histogram with rho, frames every 100 lambda_C/c)
    cs = []
    for i in range(0, len(TT), 50):
        r, _ = F.rj(TT[i])
        idx = np.floor(out[i, :3] / A + NA / 2).astype(int)
        ok = ((idx >= 0) & (idx < NA)).all(0)
        hst = np.zeros(r.shape)
        np.add.at(hst, tuple(idx[:, ok]), 1.0)
        cs.append(float(np.corrcoef(R.gaussian_filter(hst, 1).ravel(), R.gaussian_filter(r, 1).ravel())[0, 1]))
    json.dump({"correlazione_palline_rho": cs, "min": min(cs)}, open(os.path.join(SCR, "zoom_controllo.json"), "w"))
    print("correlazione palline-rho, minima:", min(cs))


def scena(F, P, iv, w, h, refs):
    s = w / 1920.0
    p = P[iv, :3].astype(np.float64)
    eta = P[iv, 3].astype(np.float64)
    t = TT[iv]
    fade = np.sin(np.pi * np.clip(eta / LIFE, 0, 1))
    wt = fade * np.exp(-0.5 * (p[2] / SLAB) ** 2)
    m = (wt > 1e-3) & (np.abs(p[0]) < 1.2 * HALF * w / h) & (np.abs(p[1]) < 1.2 * HALF)
    p, wt = p[:, m], wt[m]
    hue = (R.Dati.fase_hue(F, p, t))
    col = R.palette(hue, 0.9)
    r, _ = F.rj(t)
    rr = np.maximum(trilinear(r[None], p / A + NA / 2 - 0.5)[0], 1e-12)
    dn = np.clip((rr / r.max()) ** (-1.0 / 3.0), 0.7, 2.6)
    k = 0.5 * h / HALF
    sx, sy = w / 2 + k * p[0], h / 2 - k * p[1]
    b = R.bins_by_size(SPOT * s * dn, 8)
    img = np.zeros((4, h, w))
    for ib in np.unique(b):
        q = b == ib
        lay = np.zeros((4, h, w))
        R.splat(lay, sx[q], sy[q], wt[q] / dn[q] ** 0.5, col[q], w, h)
        img += np.stack([R.blur(lay[c], 0.7 * 2 ** (ib / 2)) for c in range(4)])
    balls, r1 = R.additivo(img, GAIN, refs.get("b"))
    glow, r2 = R.additivo(np.stack([R.blur(img[c], GLOW * s) for c in range(4)]), 1.0, refs.get("g"))
    refs.setdefault("b", r1)
    refs.setdefault("g", r2)
    return 1 - (1 - np.clip(balls, 0, 1)) * (1 - GLOW_MIX * np.clip(glow, 0, 1))


def scritte(im, iv, s):
    d = ImageDraw.Draw(im)
    f1, f2 = R.font(int(30 * s), True), R.font(int(17 * s))
    d.text((28 * s, 22 * s), "L'elettrone visto da sopra, vicino al nucleo", fill=(235, 236, 245), font=f1)
    d.text((28 * s, 62 * s), f"t = {TT[iv]:.0f} λ_C/c  ·  lato {2 * HALF:.0f} λ_C (~{2 * HALF * 0.386:.1f} pm)  ·  rallentato ~10²⁰ volte",
           fill=(210, 212, 225), font=f2)
    righe = [
        "palline = bagliori nati dove |ψ|² è grande e mossi dalla corrente della regola zig-zag; ognuna vive 24 λ_C/c e sfuma · niente scie",
        f"colore = fase dell'onda · grandezza dalla densità · fetta |z| < {SLAB:.1f} λ_C intorno al nucleo",
        "il giro attorno a z: ~0,83 della velocità di Dirac (passo del reticolo 0,5 λ_C: da solo 0,867) — errore del reticolo, niente ritocchi",
        "1s di partenza in valle k = 0 con nucleo imparato: deriva del 3,5 % in 4000 tick (non è fisica)",
        f"dati: esportazione {R.IMPRONTA[0][:10]} (atomo Z = 60 accoppiato alla luce)",
    ]
    for i, r in enumerate(righe):
        d.text((28 * s, (1080 - 30 - 24 * (len(righe) - i)) * s), r, fill=(200, 202, 215), font=f2)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "palline":
        palline()
        sys.exit()
    sc = float(a[a.index("--scala") + 1]) if "--scala" in a else 1.0
    w, h = int(1920 * sc), int(1080 * sc)
    out = os.path.join(SCR, f"zoom_{w}")
    os.makedirs(out, exist_ok=True)
    F = Campi()
    F.fa = np.memmap(os.path.join(DIR, "fase.f16"), "<f2", "r", shape=(201, NA, NA, NA))
    P = np.load(PAL, mmap_mode="r")
    rf = os.path.join(SCR, f"zoom_refs_{w}.json")
    refs = json.load(open(rf)) if os.path.exists(rf) else {}
    if not refs:   # one exposure for the whole clip, from the middle (no flicker)
        scena(F, P, len(TT) // 2, w, h, refs)
        json.dump(refs, open(rf, "w"))
    for iv in range(int(a[1]), min(int(a[2]), len(TT))):
        f = os.path.join(out, f"{iv:05d}.png")
        if os.path.exists(f):
            continue
        im = Image.fromarray((np.clip(scena(F, P, iv, w, h, refs), 0, 1).transpose(1, 2, 0) * 255).astype(np.uint8))
        scritte(im, iv, sc)
        im.save(f + ".tmp.png")
        os.replace(f + ".tmp.png", f)

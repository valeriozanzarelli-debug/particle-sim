"""Studio "bagliori": esami I1-I3 su una foto del Laboratorio (pre-registrazione lab/out/gl_prereg.json).

uso: python3 lab/gl_immagini.py foto.png [etichetta]   -> stampa un JSON con I1, I2, I3

I pannelli dell'interfaccia (testo in alto e in basso) sono esclusi: si guarda solo la fascia centrale
(righe 12 %-82 % dell'altezza).
"""
import json
import sys

import numpy as np
from PIL import Image


def esami(path):
    im = np.asarray(Image.open(path).convert("RGB")).astype(float)
    h, w, _ = im.shape
    im = im[int(0.12 * h):int(0.82 * h)]
    lum = 0.2126 * im[..., 0] + 0.7152 * im[..., 1] + 0.0722 * im[..., 2]
    fondo = np.median(lum)
    lum_n = np.clip(lum - fondo, 0, None)
    acceso = lum_n > 0.1 * lum_n.max()
    # I1 saturazione
    sat_px = (im.max(axis=2) >= 250) & acceso
    i1 = float(sat_px.sum() / max(acceso.sum(), 1))
    # I2 nucleo contro bordo
    yy, xx = np.mgrid[0:lum_n.shape[0], 0:lum_n.shape[1]]
    tot = lum_n.sum()
    cy, cx = (lum_n * yy).sum() / tot, (lum_n * xx).sum() / tot
    r = np.hypot(yy - cy, xx - cx)
    ordine = np.argsort(r.ravel())
    cum = np.cumsum(lum_n.ravel()[ordine])
    r50 = float(r.ravel()[ordine][np.searchsorted(cum, 0.5 * tot)])
    nucleo = lum_n[r < 0.5 * r50].mean()
    corona = lum_n[(r > 1.2 * r50) & (r < 1.6 * r50)].mean()
    i2 = float(nucleo / max(corona, 1e-9))
    # I3 saturazione del colore (HSV) dei pixel accesi
    mx, mn = im.max(axis=2), im.min(axis=2)
    s = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1e-9), 0)
    i3 = float(s[acceso].mean())
    return dict(foto=path, I1_frazione_saturi=i1, I1=i1 <= 0.05, I2_rapporto=i2, I2=i2 >= 2.0, R50_px=r50,
                centro_px=[float(cx), float(cy + int(0.12 * h))], I3_saturazione=i3, I3=i3 >= 0.35,
                pixel_accesi=int(acceso.sum()))


if __name__ == "__main__":
    print(json.dumps(esami(sys.argv[1]), indent=1))

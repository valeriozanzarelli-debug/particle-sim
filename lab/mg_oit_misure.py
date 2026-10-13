"""mg_oit_misure.py - exams V1 (saturated pixels in the back view) and V3 (PSNR of the side view) of the study
mg_oit (docs/MOTORE_GRAFICO_OIT.md, pre-registration lab/out/mg_oit_prereg.json).

usage: python3 mg_oit_misure.py v1 img.png [img2.png ...]
       python3 mg_oit_misure.py v3 nuovo.png vecchio.png
Pixel lit: max(R,G,B) >= 40; saturated: max >= 250 or min >= 215 (washed-out white). Fraction over lit pixels.
"""
import sys, json
import numpy as np
from PIL import Image


def carica(p):
    return np.asarray(Image.open(p).convert('RGB')).astype(np.float64)


def v1(p):
    a = carica(p)
    mx, mn = a.max(-1), a.min(-1)
    acceso = mx >= 40
    sat = (mx >= 250) | (mn >= 215)
    n = int(acceso.sum())
    lum = (0.2126 * a[..., 0] + 0.7152 * a[..., 1] + 0.0722 * a[..., 2])
    return {'file': p, 'accesi': n, 'saturi': int((sat & acceso).sum()),
            'frazione_satura': float((sat & acceso).sum() / max(n, 1)),
            'frazione_max250': float(((mx >= 250) & acceso).sum() / max(n, 1)),
            'frazione_bianco215': float(((mn >= 215) & acceso).sum() / max(n, 1)),
            'contrasto_lum_accesi': float(lum[acceso].std()) if n else 0.0,
            'lum_media_accesi': float(lum[acceso].mean()) if n else 0.0}


def psnr(a, b):
    mse = ((carica(a) - carica(b)) ** 2).mean()
    return float(10 * np.log10(255.0 ** 2 / max(mse, 1e-12)))


if __name__ == '__main__':
    if sys.argv[1] == 'v1':
        print(json.dumps([v1(p) for p in sys.argv[2:]], indent=1))
    else:
        print(json.dumps({'nuovo': sys.argv[2], 'vecchio': sys.argv[3], 'psnr_db': psnr(sys.argv[2], sys.argv[3])}))

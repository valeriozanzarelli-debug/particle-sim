"""qg_bagliori.py - the glow renderer shared by both parts (same recipe as lab/vz_render.py / lab/vl_render.py:
every pixel of the scene comes from balls drawn as soft additive gaussian glows + a wide halo, depth weight;
NO lines, trails, tubes or surfaces). Visual knobs (aesthetic only) are the module constants, counted in
docs/QUARK_GLUONI.md.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, zoom

GLOW_MIX = 0.5
GAMMA = 0.8


def camera(eye, centre, fov, w, h, up=(0.0, 0.0, 1.0)):
    fwd = np.asarray(centre, float) - eye
    fwd /= np.linalg.norm(fwd)
    right = np.cross(fwd, np.asarray(up, float))
    right /= np.linalg.norm(right)
    upv = np.cross(right, fwd)
    f = 0.5 * h / np.tan(np.radians(fov) / 2)
    return eye, fwd, right, upv, f, w, h


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


def bagliori(p, col, cam, w, h, spot, glow, gain, depth, cam_d, refs=None, keep=None, weight=None, extra=None):
    """p (3,N) world positions, col (N,3) in 0..1. Returns image (3,h,w) in 0..1 and the refs used.
    spot: glow radius in px at the reference distance; glow: halo sigma px; depth: depth fading length.
    weight: optional per-ball weight (default 1). extra: list of (p, col, radius_px, weight) bright sources."""
    sx, sy, z = project(p, cam)
    wt = np.exp(-np.clip(z - cam_d + 4.0, 0, None) / depth)
    if weight is not None:
        wt = wt * weight
    rpx = spot * (cam_d / np.maximum(z, 1e-3)) ** 1.3
    b = np.clip(np.round(np.log2(np.maximum(rpx / 0.7, 1.0)) * 2).astype(int), 0, 7)
    layers = np.zeros((8, 4, h, w))
    img2 = np.zeros((4, h, w))
    splat(img2, sx, sy, wt, col, w, h)
    sharp = np.ones(p.shape[1], bool) if keep is None else (np.arange(p.shape[1]) < keep)
    for ib in np.unique(b):
        m = (b == ib) & sharp
        splat(layers[ib], sx[m], sy[m], wt[m], col[m], w, h)
    img = np.zeros((4, h, w))
    for ib in range(8):
        if layers[ib][3].any():
            img += np.stack([blur(layers[ib][c], 0.7 * 2 ** (ib / 2)) for c in range(4)])
    refs = refs or {}
    balls, rb = additivo(img, gain, refs.get("b"))
    halo, rg = additivo(np.stack([blur(img2[c], glow) for c in range(4)]), 1.0, refs.get("g"))
    out = 1 - (1 - np.clip(balls, 0, 1)) * (1 - GLOW_MIX * np.clip(halo, 0, 1))
    if extra:
        for (pe, ce, rad, we) in extra:
            ex, ey, ez = project(pe, cam)
            src = np.zeros((4, h, w))
            splat(src, ex, ey, np.full(ex.size, we), ce, w, h)
            core = np.stack([blur(src[c], rad) for c in range(3)]) * (2 * np.pi * rad ** 2)
            wide = np.stack([blur(src[c], 4 * rad) for c in range(3)]) * (2 * np.pi * rad ** 2) * 0.35
            add = 1 - np.exp(-(core + wide))
            out = 1 - (1 - out) * (1 - np.clip(add, 0, 1))
    return out, {"b": rb, "g": rg}

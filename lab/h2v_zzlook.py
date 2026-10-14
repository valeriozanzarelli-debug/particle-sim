"""h2v_zzlook.py - the zig-zag atom look (src/atomo/zzscena.rs + render_zz.rs + balls.wgsl + main.rs) for the offline
H2 / helium clips, so the frames look like the live scene (Valerio's rule of 13 Oct: the zig-zag graphics are the base
look). Parameters COPIED from the Rust code, not invented:

  sprite (balls.wgsl, additive branch): f(q) = exp(-soft q^2) (1 - q^2) + core exp(-14 q^2), q = distance / s, q < 1,
      soft = 3.0, core = 0.3, gain = lab.luce = 1.0; colour rgb * alpha * gain * f, additive blending
  size (render_zz.rs, Scena::Tavola, bohr units): s = min(0.004 r_cam, 0.08 max(raggio_cam, 1)) * grandezza (= 1)
  alpha of a ball (zzscena.rs palline): a = 0.03 + 0.30 (rho / rho_max)^0.35
  colour (zzscena.rs tinta): cyclic 5-stop palette of the photon laboratory, angle = phase of the wave
  camera (camera.rs): vertical fov 0.8 rad;  background (main.rs ClearColor): sRGB (0.008, 0.01, 0.02)
  HDR + Bloom NATURAL intensity 0.18 (energy conserving) + TonyMcMapface tone mapping (main.rs)

APPROXIMATED here (declared): the bloom is three gaussian blurs (4, 12, 32 px at 1080p) mixed at 0.18; TonyMcMapface
(a 3-D LUT in Bevy) is replaced by a hue-preserving x/(1+x) on luminance with a gentle desaturation of highlights.
Interface for the main thread's renderer: a scene = a list of balls (position in bohr, size factor, rgba) exactly as
render.rs `Ball` = ([f32; 3], f32, [f32; 4]); the loader + draw split in h2v_render.py / he_v_render.py produces those.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, convolve
from scipy.signal import oaconvolve

BIN = 1.4                      # sprite radii grouped in bins of x1.4 (one convolution per bin)

FOV = 0.8                      # rad, camera.rs
SOFT, CORE, GAIN = 3.0, 0.3, 1.0
BG_SRGB = np.array([0.008, 0.01, 0.02])
BLOOM = 0.18
STOPS = np.array([[0.20, 0.95, 1.00], [0.25, 0.45, 1.00], [0.70, 0.30, 1.00], [1.00, 0.35, 0.70], [1.00, 0.85, 0.30]])


def tinta(angle):
    """zzscena.rs tinta(): cyclic palette for a phase (array of angles) -> (n, 3)"""
    t = np.mod(angle, 2 * np.pi) / (2 * np.pi)
    x = t * 5.0
    i = np.floor(x).astype(int) % 5
    s = (x - np.floor(x))[..., None]
    return STOPS[i] + (STOPS[(i + 1) % 5] - STOPS[i]) * s


def alpha(rho_rel):
    """zzscena.rs: brightness from the local density"""
    return 0.03 + 0.30 * np.clip(rho_rel, 0, None) ** 0.35


def size_world(r_cam, raggio_cam=1.0):
    return min(0.004 * r_cam, 0.08 * max(raggio_cam, 1.0))


def kernel(s_px):
    n = int(np.ceil(s_px))
    y, x = np.mgrid[-n:n + 1, -n:n + 1]
    q2 = (x * x + y * y) / max(s_px, 1e-6) ** 2
    k = np.where(q2 < 1, np.exp(-SOFT * q2) * (1 - q2) + CORE * np.exp(-14 * q2), 0.0)
    return k


class Tela:
    """an HDR canvas: balls are binned by on-screen sprite radius and convolved with the sprite of that radius"""
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.bins = {}

    def balls(self, sx, sy, s_px, rgb, a):
        """sx, sy: pixel centres; s_px: sprite radius in px; rgb (n, 3); a: alpha (brightness) per ball"""
        ok = (sx > -20) & (sx < self.w + 20) & (sy > -20) & (sy < self.h + 20) & (a > 0) & (s_px > 0.3)
        sx, sy, s_px, rgb, a = sx[ok], sy[ok], s_px[ok], rgb[ok], a[ok]
        key = np.clip(np.round(np.log(s_px) / np.log(BIN)).astype(int), -8, 30)
        ix = np.clip(np.round(sx).astype(int), 0, self.w - 1); iy = np.clip(np.round(sy).astype(int), 0, self.h - 1)
        lin = iy * self.w + ix
        for k in np.unique(key):
            m = key == k
            acc = self.bins.setdefault(int(k), np.zeros((3, self.h, self.w)))
            for c in range(3):
                acc[c] += np.bincount(lin[m], rgb[m, c] * a[m] * GAIN, self.w * self.h).reshape(self.h, self.w)

    def hdr(self):
        img = np.zeros((3, self.h, self.w))
        for k, acc in self.bins.items():
            ker = kernel(BIN ** k)
            if ker.shape[0] <= 5:
                for c in range(3):
                    img[c] += convolve(acc[c], ker, mode="constant")
            else:
                img += oaconvolve(acc, ker[None], mode="same", axes=(1, 2))
        return img


def finisci(img, scala_px=1.0, esposizione=1.0):
    """bloom (approx.), background, tone map (approx.), sRGB encode -> (h, w, 3) uint8"""
    bg = np.where(BG_SRGB <= 0.04045, BG_SRGB / 12.92, ((BG_SRGB + 0.055) / 1.055) ** 2.4)
    x = img * esposizione + bg[:, None, None]
    bl = np.stack([sum(gaussian_filter(x[c], s * scala_px, truncate=2.5) for s in (4.0, 12.0, 32.0)) / 3 for c in range(3)])
    x = (1 - BLOOM) * x + BLOOM * bl
    L = 0.2126 * x[0] + 0.7152 * x[1] + 0.0722 * x[2]
    Lt = L / (1 + L)
    y = x * (Lt / np.maximum(L, 1e-9))[None]
    wmix = np.clip((L - 1.0) / 6.0, 0, 0.6)[None]
    y = (1 - wmix) * y + wmix * Lt[None]
    y = np.clip(y, 0, 1)
    srgb = np.where(y <= 0.0031308, 12.92 * y, 1.055 * y ** (1 / 2.4) - 0.055)
    return (srgb.transpose(1, 2, 0) * 255 + 0.5).astype(np.uint8)

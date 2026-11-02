"""fp_render.py - drawing of FOTONE A PALLINE. Every pixel comes from balls (exam F4): each ball is a soft
gaussian spot, stretched along the local E, coloured by the direction of E around the flight axis, brighter
where the energy is denser; a short wake of the same balls at earlier ticks. No lines, tubes or surfaces.

usage: python3 fp_render.py still POL TICK VIEW out.png
       python3 fp_render.py movie POL VIEW out.mp4
VIEW: tre_quarti | lato | dentro
Visual knobs (aesthetic only, counted in the doc): SPOT, STRETCH, EXPO, GLOW, WAKE, SAT.
"""
import sys, os, subprocess
import numpy as np
from scipy.ndimage import gaussian_filter
from PIL import Image

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
W, H = 1920, 1080
SPOT = 1.4          # gaussian width of a ball (px, at 1080p)
STRETCH = 9.0       # max extra length along E (px)
EXPO = 1.6
GLOW = 0.45
WAKE = (1.0, 0.45, 0.22, 0.10)    # this tick and the previous ones
WAKE_DENTRO = tuple(0.8 ** k for k in range(14))
SAT = 0.9
CUT = False
DEPTH = 0.9        # depth of the transparency (lambda)
LAMPX = 24.0
KEEP = 25000       # balls drawn sharp (the others only glow)
SAMPLES = 7         # samples along the stretched axis
AXC = (96.0, 96.0)       # beam axis (site coords), box 192^3 (amendment 1)


def load(pol, t):
    d = np.load(os.path.join(OUT, f'fp_frames_{pol}', f't{t:03d}.npz'))
    return d['pos'].astype(np.float64), d['E'].astype(np.float64), d['B'].astype(np.float64)


def palette(h):
    """cyclic soft palette (pastel, like the old Photon Lab): h in [0,1)"""
    a = 2 * np.pi * h
    r = 0.5 + 0.5 * np.cos(a)
    g = 0.5 + 0.5 * np.cos(a - 2.094)
    b = 0.5 + 0.5 * np.cos(a + 2.094)
    rgb = np.stack([r, g, b], -1)
    return (1 - SAT) + SAT * rgb


def camera(view, centre, lam):
    c = np.array(centre)
    if view == 'tre_quarti':
        eye = c + np.array([-6.5, 3.0, -5.0]) * lam
        fov = 34
    elif view == 'lato':
        eye = c + np.array([-9.0, 0.5, 0.0]) * lam
        fov = 36
    elif view == 'dentro':
        eye = c + np.array([0.0, 0.0, -4.5]) * lam
        fov = 46
    else:
        raise ValueError(view)
    fwd = c - eye
    fwd /= np.linalg.norm(fwd)
    up0 = np.array([0.0, 1.0, 0.0])
    right = np.cross(fwd, up0); right /= np.linalg.norm(right)
    up = np.cross(right, fwd)
    f = 0.5 * H / np.tan(np.radians(fov) / 2)
    return eye, fwd, right, up, f


def project(p, cam):
    eye, fwd, right, up, f = cam
    d = p - eye[:, None]
    z = (d * fwd[:, None]).sum(0)
    x = (d * right[:, None]).sum(0)
    y = (d * up[:, None]).sum(0)
    zz = np.maximum(z, 1e-3)
    return W / 2 + f * x / zz, H / 2 - f * y / zz, z


def splat(img, pos, E, B, cam, weight, emax, umax, axis_xy, keep=None):
    """img: (4,H,W) = weighted rgb sums + weight sum (colour = average hue, brightness = density)"""
    if CUT:   # cut-away view: the half of the photon facing the camera is not drawn (a view, not geometry)
        m = pos[0] > AXC[0]
        pos, E, B = pos[:, m], E[:, m], B[:, m]
    if keep is not None:
        pos, E, B = pos[:, :keep], E[:, :keep], B[:, :keep]
    sx, sy, z = project(pos, cam)
    ex, ey, _ = project(pos + 3.0 * E / emax, cam)
    dx, dy = ex - sx, ey - sy
    nrm = np.hypot(dx, dy) + 1e-9
    ux, uy = dx / nrm, dy / nrm
    L = np.minimum(nrm, 3 * STRETCH) * 2.0
    u = 0.5 * ((E**2).sum(0) + (B**2).sum(0)) / umax
    hue = (np.arctan2(E[1], E[0]) / (2 * np.pi)) % 1.0
    col = palette(hue)
    wball = weight * (0.35 + 0.65 * np.sqrt(np.clip(u, 0, 1)))
    if DEPTH > 0:   # weighted blended transparency: near balls cover far ones (McGuire-Bavoil 2013)
        zref = np.percentile(z, 2)
        wball = wball * np.exp(-np.clip(z - zref, 0, None) / (DEPTH * LAMPX))
    near = z > 1
    ts = np.linspace(-0.5, 0.5, SAMPLES)
    wts = np.exp(-0.5 * (ts / 0.3)**2); wts /= wts.sum()
    for t, wt in zip(ts, wts):
        ix = np.round(sx + t * L * ux).astype(int); iy = np.round(sy + t * L * uy).astype(int)
        ok = near & (ix >= 0) & (ix < W) & (iy >= 0) & (iy < H)
        idx = (iy[ok], ix[ok])
        ww = wt * wball[ok]
        for c in range(3):
            np.add.at(img[c], idx, ww * col[ok, c])
        np.add.at(img[3], idx, ww)


def tone(img, spot, gain):
    rgbw = np.stack([gaussian_filter(img[c], spot) for c in range(4)])
    D = rgbw[3]
    col = rgbw[:3] / np.maximum(D, 1e-12)
    ref = np.percentile(D[D > 0], 99.7) if (D > 0).any() else 1
    lum = 1 - np.exp(-gain * D / ref)
    return col * lum, D / ref


def frame(pol, t, view, lam=24, centre=None, nwake=None):
    wake = WAKE_DENTRO if view == 'dentro' else WAKE
    nwake = len(wake) if nwake is None else nwake
    pos, E, B = load(pol, t)
    if centre is None:
        centre = (AXC[0], AXC[1], float(np.median(pos[2])))
    cam = camera(view, centre, lam)
    emax = np.percentile(np.sqrt((E**2).sum(0)), 99.5)
    umax = np.percentile(0.5 * ((E**2).sum(0) + (B**2).sum(0)), 99.5)
    # layer 1: a subset of balls drawn sharp, with a wake made of the same balls at earlier ticks
    img = np.zeros((4, H, W))
    for k in range(nwake):
        if t - k < 0:
            break
        p, e, b = (pos, E, B) if k == 0 else load(pol, t - k)
        splat(img, p, e, b, cam, wake[k], emax, umax, AXC, keep=KEEP)
    balls, _ = tone(img, SPOT, 1.6)
    # layer 2: all the balls, blurred wide = glow where the energy is dense
    img2 = np.zeros((4, H, W))
    splat(img2, pos, E, B, cam, 1.0, emax, umax, AXC)
    glow, _ = tone(img2, 14.0, 1.2)
    out = 1 - (1 - np.clip(balls, 0, 1)) * (1 - GLOW * np.clip(glow, 0, 1))
    bg = np.array([0.010, 0.012, 0.020])[:, None, None]
    out = np.clip(out + bg * (1 - out), 0, 1) ** (1 / 1.1)
    return (out.transpose(1, 2, 0) * 255).astype(np.uint8)


if __name__ == '__main__':
    mode = sys.argv[1]
    if '--taglio' in sys.argv:
        sys.argv.remove('--taglio'); CUT = True
    if mode == 'still':
        pol, t, view, dst = sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
        Image.fromarray(frame(pol, t, view)).save(dst)
    elif mode == 'movie':
        pol, view, dst = sys.argv[2], sys.argv[3], sys.argv[4]
        tmp = os.path.join(OUT, f'fp_movie_{pol}_{view}')
        os.makedirs(tmp, exist_ok=True)
        t0, t1 = int(sys.argv[5]) if len(sys.argv) > 5 else 0, int(sys.argv[6]) if len(sys.argv) > 6 else 96
        for t in range(t0, t1 + 1):
            f = os.path.join(tmp, f'{t:03d}.png')
            if not os.path.exists(f):
                Image.fromarray(frame(pol, t, view)).save(f)
        subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', '24', '-i', os.path.join(tmp, '%03d.png'),
                        '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '18', dst], check=True)

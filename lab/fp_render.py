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
SPOT = 1.1          # gaussian width of a ball (px, at 1080p)
STRETCH = 7.0       # max extra length along E (px)
EXPO = 1.6
GLOW = 0.55
WAKE = (1.0, 0.45, 0.22, 0.10)    # this tick and the previous ones
SAT = 0.78
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
        eye = c + np.array([-5.0, 3.2, -4.2]) * lam
        fov = 34
    elif view == 'lato':
        eye = c + np.array([-7.5, 0.6, 0.0]) * lam
        fov = 36
    elif view == 'dentro':
        eye = c + np.array([0.35, 0.25, -3.6]) * lam
        fov = 70
    else:
        raise ValueError(view)
    fwd = c - eye
    fwd /= np.linalg.norm(fwd)
    up0 = np.array([0.0, 1.0, 0.0]) if view != 'dentro' else np.array([0.0, 1.0, 0.0])
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


def splat(img, pos, E, B, cam, weight, emax, umax, axis_xy):
    sx, sy, z = project(pos, cam)
    ex, ey, _ = project(pos + 2.0 * E / emax, cam)
    dx, dy = ex - sx, ey - sy
    nrm = np.hypot(dx, dy) + 1e-9
    ux, uy = dx / nrm, dy / nrm
    emag = np.sqrt((E**2).sum(0)) / emax
    u = 0.5 * ((E**2).sum(0) + (B**2).sum(0)) / umax
    L = STRETCH * np.clip(emag, 0, 1) * (cam[4] / np.maximum(z, 1)) / (cam[4] / 300.0)
    L = np.minimum(L, 3 * STRETCH)
    hue = (np.arctan2(E[1], E[0]) / (2 * np.pi)) % 1.0
    col = palette(hue) * (weight * (0.25 + 0.75 * np.sqrt(np.clip(u, 0, 1))))[:, None]
    near = z > 1
    ts = np.linspace(-0.5, 0.5, SAMPLES)
    wts = np.exp(-0.5 * (ts / 0.3)**2); wts /= wts.sum()
    for t, wt in zip(ts, wts):
        px = sx + t * L * ux
        py = sy + t * L * uy
        ix = np.round(px).astype(int); iy = np.round(py).astype(int)
        ok = near & (ix >= 0) & (ix < W) & (iy >= 0) & (iy < H)
        for c in range(3):
            np.add.at(img[c], (iy[ok], ix[ok]), wt * col[ok, c])


def frame(pol, t, view, lam=24, centre=None, nwake=len(WAKE)):
    pos, E, B = load(pol, t)
    if centre is None:
        centre = (AXC[0], AXC[1], float(np.median(pos[2])))
    cam = camera(view, centre, lam)
    emax = np.percentile(np.sqrt((E**2).sum(0)), 99.5)
    umax = np.percentile(0.5 * ((E**2).sum(0) + (B**2).sum(0)), 99.5)
    img = np.zeros((3, H, W))
    for k in range(nwake):
        if t - k < 0:
            break
        p, e, b = (pos, E, B) if k == 0 else load(pol, t - k)
        splat(img, p, e, b, cam, WAKE[k], emax, umax, AXC)
    # soft balls: small blur = the ball itself; wide blur = glow where energy is dense
    scale = img.mean() * 1e3 + 1e-12
    img /= scale
    ball = np.stack([gaussian_filter(img[c], SPOT) for c in range(3)])
    glow = np.stack([gaussian_filter(img[c], 9.0) + 0.5 * gaussian_filter(img[c], 28.0) for c in range(3)])
    hdr = ball * 22.0 + GLOW * glow * 60.0
    out = 1 - np.exp(-EXPO * hdr)
    bg = np.array([0.012, 0.014, 0.022])[:, None, None]
    out = np.clip(out + bg * (1 - out), 0, 1) ** (1 / 1.15)
    return (out.transpose(1, 2, 0) * 255).astype(np.uint8)


if __name__ == '__main__':
    mode = sys.argv[1]
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

"""vl_comune.py - shared bits of VIDEO ATOMO: scratch paths, video timeline, trilinear gather."""
import os

import numpy as np

SCR = os.environ.get("VL_SCR", "/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vl")
QUI = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(QUI, "out")
DT = 8.0          # time between data frames (lambda_C/c)
T_END = 1600.0
A, NA = 0.5, 48   # atom grid
H, NL = 8.0, 64   # light grid (downsampled)
BOX_L = 208.0     # light: inside the sponge
CHUNK = 60        # video frames per stored chunk


def rate(t):
    """video frames per data frame: 12 (slow) until t=280, 4 after t=360, linear blend in between"""
    if t <= 280:
        return 12.0
    if t >= 360:
        return 4.0
    return 12.0 + (4.0 - 12.0) * (t - 280) / 80.0


def timeline():
    """video times: every data interval [8k, 8k+8] is split into round(rate) equal steps (data frames are hit)"""
    ts = [0.0]
    for k in range(int(T_END / DT)):
        r = int(round(rate(k * DT + DT / 2)))
        ts += list(k * DT + DT * np.arange(1, r + 1) / r)
    return np.array(ts)


def idx_atom(p):
    return p / A + (NA / 2 - 0.5)


def idx_light(p):
    return (p / 4.0 + NL - 0.5) / 2.0


def trilinear(F, ic):
    """F: (c, n, n, n) fields; ic: (3, N) fractional indices -> (c, N); clamps at the border"""
    n = F.shape[1]
    ic = np.clip(ic, 0, n - 1 - 1e-6)
    i0 = np.floor(ic).astype(np.int64)
    f = ic - i0
    out = np.zeros((F.shape[0], ic.shape[1]), np.float64)
    Ff = F.reshape(F.shape[0], -1)
    for dx in (0, 1):
        wx = f[0] if dx else 1 - f[0]
        for dy in (0, 1):
            wy = f[1] if dy else 1 - f[1]
            for dz in (0, 1):
                wz = f[2] if dz else 1 - f[2]
                lin = ((i0[0] + dx) * n + (i0[1] + dy)) * n + (i0[2] + dz)
                out += Ff[:, lin] * (wx * wy * wz)
    return out


def catmull(frames, t, get):
    """Catmull-Rom in time of grid frames: get(k) -> array; frames every DT"""
    nf = frames
    s = t / DT
    k = int(min(max(np.floor(s), 0), nf - 2))
    u = s - k
    km, kp, kpp = max(k - 1, 0), k + 1, min(k + 2, nf - 1)
    p0, p1, p2, p3 = get(km), get(k), get(kp), get(kpp)
    return 0.5 * ((2 * p1) + (-p0 + p2) * u + (2 * p0 - 5 * p1 + 4 * p2 - p3) * u ** 2 + (-p0 + 3 * p1 - 3 * p2 + p3) * u ** 3)

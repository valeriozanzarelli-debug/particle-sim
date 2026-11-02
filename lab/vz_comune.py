"""vz_comune.py - VIDEO SUPERCRITICO: data access, balls (mass coordinate), timeline. Shared by vz_esami / vz_render.

The data (lab/out/zstar_film.json at commit fc8a08c, branch claude/problemi-aperti-56mnej) are radial: rho(r), j(r),
Q(<r) for 301 frames. Each ball has a fixed mass coordinate m_i and a fixed isotropic direction; its radius is
r_i(t) = Q^{-1}(m_i, t), the exact integral of dr/dt = j/rho in 1-D (Q linear in r inside a cell and linear in t
between frames, i.e. j constant over the interval as in the data). See lab/out/vz_prereg.json.
"""
import json
import os
import subprocess

import numpy as np

SCR = os.environ.get("VZ_SCR", "/tmp/claude-0/-home-claude-particle-sim/a13a264e-04a0-55c0-8c14-0e293e6cf6fd/scratchpad/vz")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMMIT = "fc8a08c"
N = 300000
SEED = 20261013
R0_LOG = 0.02      # draw scale: R = ln(1 + r/R0_LOG) / ln(1 + RMAX/R0_LOG)
RMAX = 40.0


def carica_dati():
    os.makedirs(SCR, exist_ok=True)
    f = os.path.join(SCR, "zstar_film.json")
    if not os.path.exists(f):
        subprocess.run(["git", "-C", REPO, "fetch", "-q", "origin", "claude/problemi-aperti-56mnej"], check=True)
        with open(f, "w") as o:
            o.write(subprocess.run(["git", "-C", REPO, "show", f"{COMMIT}:lab/out/zstar_film.json"],
                                   check=True, capture_output=True, text=True).stdout)
    return json.load(open(f))


D = carica_dati()
MAN = D["manifest"]
RB = np.array(D["r_bordi"])                 # 385 edges
RC = np.array(D["r_centri"])
T = np.array(D["t"])                        # 301 times
DT = float(MAN["dt_fotogramma"])
RHO = np.array(D["rho"])
J = np.array(D["j"])
Q0 = np.array(D["Q"])
QE = np.concatenate([np.zeros((len(T), 1)), Q0], 1)
QE = np.maximum.accumulate(QE, 1)            # 4-digit rounding can break monotonicity: cumulative max (declared)
N_NONMONO = int((np.diff(np.concatenate([np.zeros((len(T), 1)), Q0], 1), axis=1) < 0).sum())
SC = D["fotogrammi_scalari"]
RW = float(MAN["risonanza"]["r_w"])
ER = float(MAN["risonanza"]["E_r"])
Z = int(MAN["Z"])
ZA = Z / 137.035999
NF = len(T)


def palline(n=N, seed=SEED):
    """mass coordinates (stratified) and isotropic unit directions, fixed for the whole film"""
    rng = np.random.default_rng(seed)
    m = (np.arange(n) + rng.random(n)) / n
    rng.shuffle(m)
    u = rng.uniform(-1, 1, n)
    ph = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - u * u)
    d = np.stack([s * np.cos(ph), s * np.sin(ph), u])
    return m, d


def Q_at(t):
    k = int(min(np.floor(t / DT + 1e-9), NF - 2))
    u = t / DT - k
    return (1 - u) * QE[k] + u * QE[k + 1]


def raggi(m, t):
    """r_i(t) = Q^{-1}(m_i, t); nan if the ball has left through r = 40 (the computation's sponge)"""
    q = Q_at(t)
    e = np.searchsorted(q, m, side="left")           # q[e-1] < m <= q[e]
    out = e >= len(q)
    e = np.clip(e, 1, len(q) - 1)
    lo, hi = q[e - 1], q[e]
    f = np.where(hi > lo, (m - lo) / np.maximum(hi - lo, 1e-300), 0.0)
    r = RB[e - 1] + f * (RB[e] - RB[e - 1])
    r[out] = np.nan
    return r


def disegno(r):
    """radial draw scale (logarithmic, declared on screen)"""
    return np.log1p(r / R0_LOG) / np.log1p(RMAX / R0_LOG)


# ------------------------------------------------------------------ timeline: n(k) video frames per data frame
def timeline():
    n = 2.5 + 9.0 * np.exp(-np.arange(NF - 1) / 30.0)
    tv = [0.0]
    for k in range(NF - 1):
        nk = int(round(n[k]))
        for i in range(1, nk + 1):
            tv.append(T[k] + DT * i / nk)
    return np.array(tv)


TV = timeline()
NV = len(TV)

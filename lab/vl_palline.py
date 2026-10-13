"""vl_palline.py - VIDEO ATOMO: move the balls (electron and light) along the whole video timeline.

    python3 lab/vl_palline.py elettrone     # resumable: runs <= ~8 min, then checkpoints; call again until DONE
    python3 lab/vl_palline.py luce

Electron balls: seeded once (t=0) with density rho_s (|psi|^2 smoothed sigma = 1 lambda_C), moved with
v = J/rho_s (J = the rule's current, smoothed, minimally corrected to obey continuity; face fluxes, vl_campi.py).
Light balls: pulse balls seeded once from u = (E^2+B^2)/2; atom-light balls born on the sphere R0 = 40 from
t = 240 with probability prop. to the outgoing flux, each worth 1/1000 of a pulse ball (declared x1000).
Both moved with the energy velocity 2 E x B / (E^2 + B^2).
Output (scratch, not in git): <SCR>/el_XXX.npy / lu_XXX.npz chunks of CHUNK video frames, el_dati.npy (positions at
the 201 data frames, for the exams), lu_dati_kKK.npy.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vl_comune import (SCR, OUT, DT, A, NA, H, NL, BOX_L, CHUNK, timeline, idx_light, trilinear,  # noqa: E402
                       catmull)

BUDGET = float(os.environ.get("VL_BUDGET", 470))
N_E = 120000
N_A = 150000
SEED = 20261013
GAIN = 1000.0
R0 = 40.0
T_ON = 240.0
VMAX = 0.95


def eventi(tv):
    """union of video times and data-frame times (sorted), with flags"""
    td = np.arange(0, 201) * DT
    allt = np.unique(np.round(np.concatenate([tv, td]), 9))
    return allt


# ---------------------------------------------------------------- electron
class Elettrone:
    """RT0 flow inside each cell: face fluxes J (constant over a frame interval), rho of the cell linear in time"""
    def __init__(self):
        self.rho = np.load(os.path.join(SCR, "rho_s.npy"), mmap_mode="r")
        self.F = [np.load(os.path.join(SCR, f"F{c}.npy"), mmap_mode="r") for c in "xyz"]
        self.cache_k = None
        self.clipped = 0
        self.outside = 0

    def load(self, k):
        if self.cache_k != k:
            self.cache = ([np.asarray(f[k], np.float64) for f in self.F],
                          np.asarray(self.rho[k], np.float64).ravel(), np.asarray(self.rho[k + 1], np.float64).ravel())
            self.cache_k = k
        return self.cache

    def v(self, p, t):
        k = int(min(np.floor(t / DT + 1e-12), 199))
        s = t / DT - k
        F, r0, r1 = self.load(k)
        x = np.clip(p / A + NA / 2, 0, NA - 1e-9)
        c = np.floor(x).astype(np.int64)
        f = x - c
        lin = (c[0] * NA + c[1]) * NA + c[2]
        rho = np.maximum((1 - s) * r0[lin] + s * r1[lin], 1e-14)
        v = np.empty_like(p)
        for a in range(3):
            lo = [c[0], c[1], c[2]]
            hi = list(lo); hi[a] = hi[a] + 1
            v[a] = ((1 - f[a]) * F[a][lo[0], lo[1], lo[2]] + f[a] * F[a][hi[0], hi[1], hi[2]]) / rho
        sp = np.sqrt((v ** 2).sum(0))
        m = sp > VMAX
        if m.any():
            v[:, m] *= VMAX / sp[m]
        return v, int(m.sum())

    def step(self, p, t, dt):
        v1, _ = self.v(p, t)
        pm = p + 0.5 * dt * v1
        v2, c = self.v(pm, t + 0.5 * dt)
        self.clipped += c
        q = p + dt * v2
        lim = (NA / 2) * A - 1e-6
        out = np.abs(q) > lim
        self.outside += int(out.any(0).sum())
        return np.clip(q, -lim, lim)


def semina_elettrone():
    rng = np.random.default_rng(SEED)
    rho0 = np.load(os.path.join(SCR, "rho_s.npy"), mmap_mode="r")[0].astype(np.float64)
    p = rho0.ravel() / rho0.sum()
    cells = rng.choice(p.size, size=N_E, p=p)
    i, j, k = np.unravel_index(cells, rho0.shape)
    pos = np.stack([i, j, k]).astype(np.float64) + rng.random((3, N_E))
    return (pos - NA / 2) * A   # cell i spans [(i - 24) a, (i - 23) a]


# ---------------------------------------------------------------- light
class Luce:
    def __init__(self):
        self.E = np.load(os.path.join(SCR, "E.npy"), mmap_mode="r")
        self.B = np.load(os.path.join(SCR, "B.npy"), mmap_mode="r")
        self.cache = {}
        M = 3000
        i = np.arange(M) + 0.5
        th = np.arccos(1 - 2 * i / M)
        ph = np.pi * (1 + 5 ** 0.5) * i
        self.nh = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])

    def frame(self, k):
        if k not in self.cache:
            if len(self.cache) > 6:
                self.cache.pop(min(self.cache))
            self.cache[k] = np.concatenate([np.asarray(self.E[k], np.float64), np.asarray(self.B[k], np.float64)])
        return self.cache[k]

    def fields(self, t):
        return catmull(201, t, self.frame)

    def v(self, p, t, F=None):
        F = self.fields(t) if F is None else F
        f = trilinear(F, idx_light(p))
        E, B = f[:3], f[3:]
        S = np.cross(E, B, axis=0)
        u2 = (E ** 2).sum(0) + (B ** 2).sum(0)
        v = 2 * S / np.maximum(u2, 1e-30)
        sp = np.sqrt((v ** 2).sum(0))
        m = sp > 1.0
        v[:, m] /= sp[m]
        return v

    def step(self, p, t, dt):
        v1 = self.v(p, t)
        v2 = self.v(p + 0.5 * dt * v1, t + 0.5 * dt)
        return p + dt * v2

    def flux_points(self, t):
        F = self.fields(t)
        f = trilinear(F, idx_light(R0 * self.nh))
        S = np.cross(f[:3], f[3:], axis=0)
        return (S * self.nh).sum(0)   # S.n at each point (E-units^2)


def semina_luce():
    rng = np.random.default_rng(SEED + 1)
    E = np.load(os.path.join(SCR, "E.npy"), mmap_mode="r")[0].astype(np.float64)
    B = np.load(os.path.join(SCR, "B.npy"), mmap_mode="r")[0].astype(np.float64)
    u = 0.5 * ((E ** 2).sum(0) + (B ** 2).sum(0))
    X = (2 * np.arange(NL) + 0.5 - NL) * 4.0
    inside = (np.abs(X) + H / 2 < BOX_L)
    m = inside[:, None, None] & inside[None, :, None] & inside[None, None, :]
    u = np.where(m, u, 0.0)
    U = u.sum() * H ** 3
    p = u.ravel() / u.sum()
    cells = rng.choice(p.size, size=N_A, p=p)
    i, j, k = np.unravel_index(cells, u.shape)
    pos = np.stack([X[i], X[j], X[k]]) + (rng.random((3, N_A)) - 0.5) * H
    return pos, U / N_A


# ---------------------------------------------------------------- driver
def run(kind):
    tv = timeline()
    nv = len(tv)
    ck = os.path.join(SCR, f"{kind}_ckpt.npz")
    t0 = time.time()
    if kind == "elettrone":
        sim = Elettrone()
        if os.path.exists(ck):
            z = np.load(ck)
            pos, iv, dati = z["pos"], int(z["iv"]), z["dati"]
            sim.clipped, sim.outside = int(z["clipped"]), int(z["outside"])
        else:
            pos = semina_elettrone()
            iv, dati = 0, np.full((201, 3, N_E), np.nan, np.float32)
            dati[0] = pos
    else:
        sim = Luce()
        if os.path.exists(ck):
            z = np.load(ck, allow_pickle=False)
            pos, kindb, iv = z["pos"], z["kindb"], int(z["iv"])
            q_a, born, dead = float(z["q_a"]), int(z["born"]), int(z["dead"])
        else:
            pos, q_a = semina_luce()
            kindb = np.zeros(pos.shape[1], np.int8)
            iv, born, dead = 0, 0, 0
        q_b = q_a / GAIN
    buf = []
    chunk_start = (iv // CHUNK) * CHUNK
    # reload partially filled chunk (if any)
    part = os.path.join(SCR, f"{kind}_part.npz")
    if os.path.exists(part) and iv > chunk_start:
        zp = np.load(part, allow_pickle=True)
        buf = list(zp["buf"])
    if iv == 0:
        buf = [(pos.astype(np.float16), kindb.copy()) if kind == "luce" else pos.astype(np.float16)]
        iv = 1
    td = np.arange(201) * DT
    while iv < nv:
        ta, tb = tv[iv - 1], tv[iv]
        # sub-steps, landing on data-frame times
        cuts = [ta] + [x for x in td if ta < x < tb] + [tb]
        for a_, b_ in zip(cuts[:-1], cuts[1:]):
            nsub = int(np.ceil((b_ - a_) / (0.5 if kind == "elettrone" else 2.0) - 1e-9))
            hdt = (b_ - a_) / nsub
            for s in range(nsub):
                tt = a_ + s * hdt
                if kind == "elettrone":
                    pos = sim.step(pos, tt, hdt)
                else:
                    pos = sim.step(pos, tt, hdt)
                    # deaths in the sponge
                    alive = (np.abs(pos) < BOX_L).all(0)
                    if not alive.all():
                        dead += int((~alive).sum())
                        pos, kindb = pos[:, alive], kindb[alive]
                    # births of the atom's light
                    tm = tt + 0.5 * hdt
                    if tm >= T_ON:
                        Sn = sim.flux_points(tm)
                        Sp = np.clip(Sn, 0, None)
                        Pout = Sp.mean() * 4 * np.pi * R0 ** 2
                        rng = np.random.default_rng([SEED + 2, int(round(tm * 1000))])
                        nb = rng.poisson(Pout * hdt / q_b)
                        if nb > 0 and Sp.sum() > 0:
                            ii = rng.choice(Sp.size, size=nb, p=Sp / Sp.sum())
                            d = sim.nh[:, ii] + 0.04 * rng.standard_normal((3, nb))
                            d /= np.linalg.norm(d, axis=0)
                            # born somewhere along the last sub-step's path (no shells from the step size)
                            rr = R0 + rng.random(nb) * hdt
                            pos = np.concatenate([pos, d * rr], axis=1)
                            kindb = np.concatenate([kindb, np.ones(nb, np.int8)])
                            born += nb
            if kind == "elettrone":
                kd = int(round(b_ / DT))
                if abs(b_ - kd * DT) < 1e-6:
                    dati[kd] = pos
            else:
                kd = int(round(b_ / DT))
                if abs(b_ - kd * DT) < 1e-6 and kd in (10, 20, 30, 60, 100, 200):
                    np.save(os.path.join(SCR, f"lu_dati_k{kd:03d}.npy"), pos[:, kindb == 0].astype(np.float32))
        buf.append((pos.astype(np.float16), kindb.copy()) if kind == "luce" else pos.astype(np.float16))
        iv += 1
        if iv % CHUNK == 0 or iv == nv:
            c = (iv - 1) // CHUNK
            if kind == "elettrone":
                np.save(os.path.join(SCR, f"el_{c:03d}.npy"), np.stack(buf))
            else:
                np.savez(os.path.join(SCR, f"lu_{c:03d}.npz"),
                         **{f"p{i}": b[0] for i, b in enumerate(buf)}, **{f"k{i}": b[1] for i, b in enumerate(buf)})
            buf = []
            if os.path.exists(part):
                os.remove(part)
        if time.time() - t0 > BUDGET or iv == nv:
            if buf:
                arr = np.empty(len(buf), object)
                arr[:] = buf
                np.savez(part, buf=arr)
            if kind == "elettrone":
                np.savez(ck, pos=pos, iv=iv, dati=dati, clipped=sim.clipped, outside=sim.outside)
                if iv == nv:
                    np.save(os.path.join(SCR, "el_dati.npy"), dati)
            else:
                np.savez(ck, pos=pos, kindb=kindb, iv=iv, q_a=q_a, born=born, dead=dead)
            info = {"kind": kind, "iv": iv, "nv": nv, "t": float(tv[iv - 1]), "secondi": time.time() - t0}
            if kind == "elettrone":
                info.update(clipped=sim.clipped, outside=sim.outside)
            else:
                info.update(alive=int(pos.shape[1]), atom_balls=int((kindb == 1).sum()), born=born, dead=dead, q_a=q_a)
            print(json.dumps(info))
            if iv == nv:
                with open(os.path.join(SCR, f"{kind}_DONE.json"), "w") as f:
                    json.dump(info, f)
            return


if __name__ == "__main__":
    run(sys.argv[1])

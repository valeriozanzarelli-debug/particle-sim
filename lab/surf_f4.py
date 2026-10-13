"""surf_f4.py - PHASE 4 of SURFISTI_ESATTI: the 1s on valley k = 0 with the TRANSFERRED learned nucleus (NUCLEO_1S sec. 14,
straight line s(a) through a = 0.5 and 0.3) at a = 0.5, 0.4, 0.3, 0.25 (box 24 lambda_C): does the missing 13 % of the
Dirac rotation go down with a, and with what law? Exams G1-G6 in lab/out/surf_f4_prereg.json.
usage: python surf_f4.py prep <a>   -> spectrum pass + 2 filtering passes, projection on valley 0   [checkpointed]
       python surf_f4.py run <a>    -> physical time 2000 lambda_C/c, 2e5 exact surfers, 5 windows   [checkpointed]
       python surf_f4.py judge      -> lab/out/surf_f4_results.json, lab/out/surf_f4_figure.png
"""
import json
import os
import pickle
import sys
import time
import numpy as np
from surf_core import Walk3D, net_flux, SEED, OUT, born_sample, chi2_z, flow_tick, density, radial_state, spinor3d_from_radial, ALPHA
from zz_core import PAULI

ZA = 60 * ALPHA
L_BOX = 24.0
CK = os.path.join(OUT, "surf_ckpt")
SA, SB = np.array([0.554593, 1.003549]), np.array([0.57482, 1.002538])     # learn2 (a = 0.5), learn3 (a = 0.3)
T_PHYS = 2000.0
W_PHYS = (250.0, 500.0, 1000.0, 1500.0)                                   # + an end window closing at T_PHYS
NS = 200000
VALLEYS = {"0": (0, 0, 0), "xy": (np.pi / 2, np.pi / 2, 0), "xz": (np.pi / 2, 0, np.pi / 2), "yz": (0, np.pi / 2, np.pi / 2)}


def ev(x):
    return int(2 * round(x / 2))


class Grid:
    def __init__(self, a, t_phys=T_PHYS):
        self.a = a
        self.n = int(round(L_BOX / a))
        x = (np.arange(self.n) - self.n / 2 + 0.5) * a
        self.X, self.Y, self.Z = np.meshgrid(x, x, x, indexing="ij")
        self.R = np.sqrt(self.X**2 + self.Y**2 + self.Z**2)
        self.s = SA + (SB - SA) * (a - 0.5) / (0.3 - 0.5)
        V = -ZA / self.R
        V = np.where(np.abs(self.R - np.sqrt(3) / 2 * a) < 1e-6 * a, self.s[0] * V, V)
        self.V = np.where(np.abs(self.R - np.sqrt(11) / 2 * a) < 1e-6 * a, self.s[1] * V, V)
        self.n0 = int(np.sum(np.abs(self.R - np.sqrt(3) / 2 * a) < 1e-6 * a))
        self.n1 = int(np.sum(np.abs(self.R - np.sqrt(11) / 2 * a) < 1e-6 * a))
        self.dur = ev(5.0 / a)                                            # window length ~ 5 lambda_C / c
        self.T = ev(t_phys / a)
        self.wins = [ev(t / a) for t in W_PHYS if t < t_phys - 10] + [self.T - self.dur]
        self.listen = ev(1024.0 / a)                                      # 2048 ticks at a = 0.5

    def walk(self):
        return Walk3D(self.n, self.a, V=self.V * self.a)

    def seed(self):
        om, chi, rw = radial_state(-1, 1, self.a, ZA, 40.0)
        p = spinor3d_from_radial(rw, chi, -1, self.X, self.Y, self.Z)
        return p / np.linalg.norm(p)

    def valley(self, psi, name="0"):
        k = 2 * np.pi * np.fft.fftfreq(self.n)
        KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
        best, lab = None, None
        for i, (nm, q) in enumerate(VALLEYS.items()):
            dmin = None
            for m in np.ndindex(2, 2, 2):
                qq = np.array(q) + np.pi * np.array(m)
                d = sum(np.abs((Kc - qc + np.pi) % (2 * np.pi) - np.pi)**2 for Kc, qc in zip((KX, KY, KZ), qq))
                dmin = d if dmin is None else np.minimum(dmin, d)
            if best is None:
                best, lab = dmin, np.zeros(dmin.shape, np.int8)
            else:
                better = dmin < best
                lab[better] = i; best = np.minimum(best, dmin)
        del KX, KY, KZ, best
        mask = (lab == list(VALLEYS).index(name))
        return np.fft.ifftn(np.fft.fftn(psi, axes=(1, 2, 3)) * mask[None], axes=(1, 2, 3))

    def vphi_flux(self, F, rho, per_tick):
        Fc = [0.5 * (F[j] + np.roll(F[j], 1, axis=j)) / per_tick for j in range(3)]
        rc = np.sqrt(self.X**2 + self.Y**2)
        m = (self.R > 1) & (self.R < 6) & (np.abs(self.Z) < 0.5 * self.R)
        J = (-self.Y * Fc[0] + self.X * Fc[1]) / rc
        return float(np.sum((J * self.R / rc)[m]) / np.sum(rho[m]) / ZA)


def sigma_z(psi):
    return sum(float(np.real(np.sum(np.einsum("a...,ab,b...->...", psi[q:q + 2].conj(), PAULI[2], psi[q:q + 2])))) for q in (0, 2)) / float(np.sum(np.abs(psi)**2))


def tag(a):
    return f"{a:.2f}"


def prep(a, budget=540.0):
    g = Grid(a)
    f = os.path.join(CK, f"f4_prep_{tag(a)}.pkl")
    t_start = time.time()
    TL = g.listen
    if os.path.exists(f):
        S = pickle.load(open(f, "rb"))
    else:
        s0 = g.seed()
        S = dict(pas=0, t=0, psi0=s0, psi=s0.copy(), ov=[], om=None, acc=None)
    W = g.walk()
    while S["pas"] <= 2 and time.time() - t_start < budget:
        W.tick_count = S["t"]
        psi = S["psi"]
        while S["t"] < TL and time.time() - t_start < budget:
            t = S["t"]
            if t % 2 == 0:
                w = np.sin(np.pi * (t + 1.0) / TL)**2
                if S["pas"] == 0:
                    S["ov"].append(np.vdot(S["psi0"], psi))
                else:
                    S["acc"] += w * np.exp(1j * S["om"] * t * a) * psi
            psi, _ = W.tick(psi)
            S["t"] += 1
        S["psi"] = psi
        if S["t"] < TL:
            break
        if S["pas"] == 0:
            ov = np.array(S["ov"]); tt = np.arange(0, TL, 2) * a
            hann = np.sin(np.pi * (np.arange(0, TL, 2) + 1.0) / TL)**2
            fine = np.arange(0.86, 0.94, 0.0001)
            spec = np.abs((ov * hann)[None, :] @ np.exp(1j * np.outer(tt, fine)))[0]**2
            S["om"] = float(fine[int(np.argmax(spec))])
            print(a, "spectrum peak", S["om"], flush=True)
            S.update(pas=1, t=0, psi=S["psi0"].copy(), acc=np.zeros_like(S["psi0"]), ov=None)
        elif S["pas"] == 1:
            v = S["acc"] / np.linalg.norm(S["acc"])
            S.update(pas=2, t=0, psi=v, acc=np.zeros_like(v))
            print(a, "filter pass 1 done", flush=True)
        else:
            v = S["acc"] / np.linalg.norm(S["acc"])
            v0 = g.valley(v)
            wv0 = float(np.sum(np.abs(v0)**2))
            p = v0 / np.sqrt(wv0)
            info = dict(a=a, n=g.n, s=g.s.tolist(), cells_r0=g.n0, cells_r1=g.n1, peak=S["om"],
                        filtered=dict(sigma_z=sigma_z(v), valley0=wv0),
                        prepared=dict(sigma_z=sigma_z(p), valley0=float(np.sum(np.abs(g.valley(p))**2))))
            np.save(os.path.join(CK, f"f4_psi0_{tag(a)}.npy"), p)
            json.dump(info, open(os.path.join(OUT, f"surf_f4_prep_{tag(a)}.json"), "w"), indent=1, default=float)
            os.remove(f) if os.path.exists(f) else None
            print(a, "PREP DONE", info, flush=True)
            return
    pickle.dump(S, open(f, "wb")); print(a, "prep checkpoint pass", S["pas"], "tick", S["t"], flush=True); sys.exit(3)


def run(a, t_phys=T_PHYS, budget=540.0):
    g = Grid(a, t_phys)
    f = os.path.join(CK, f"f4_run_{tag(a)}.pkl")
    t_start = time.time()
    p0 = np.load(os.path.join(CK, f"f4_psi0_{tag(a)}.npy"))
    if os.path.exists(f):
        S = pickle.load(open(f, "rb"))
    else:
        rng = np.random.default_rng(SEED + 400 + int(100 * a))
        c, x = born_sample(rng, density(p0), NS)
        S = dict(t=0, psi=p0.copy(), c=c, x=x, start={}, win={},
                 rec=dict(t=[], fid=[], valley0=[], sigma_z=[], flux_vphi=[], z=[]))
    W = g.walk()
    every = ev(t_phys / a / 40)
    pos = lambda: (S["c"] + S["x"] - g.n / 2) * a
    while time.time() - t_start < budget:
        t = S["t"]
        for w0 in g.wins:
            if t == w0 + g.dur and w0 in S["start"]:
                q0 = S["start"].pop(w0); q1 = pos()
                rr = np.sqrt(np.sum(q0**2, 1)); rc = np.sqrt(q0[:, 0]**2 + q0[:, 1]**2)
                sel = (rr > 1) & (rr < 6) & (np.abs(q0[:, 2]) < 0.5 * rr)
                dphi = (np.arctan2(q1[:, 1], q1[:, 0]) - np.arctan2(q0[:, 1], q0[:, 0]) + np.pi) % (2 * np.pi) - np.pi
                vphi = rc * dphi / (g.dur * a)
                r1 = np.sqrt(np.sum(q1**2, 1))
                S["win"][w0] = dict(t_phys=w0 * a, vphi_over_Za_sin=float(np.mean(vphi[sel] / (rc[sel] / rr[sel])) / ZA),
                                    vr_rms_over_Za=float(np.sqrt(np.mean(((r1 - rr) / (g.dur * a))[sel]**2)) / ZA))
                print(a, "window", w0, S["win"][w0], flush=True)
            if t == w0:
                S["start"][w0] = pos().copy()
        if t >= g.T:
            break
        psi = S["psi"]
        ra = density(psi)
        W.tick_count = t
        F2 = [0, 0, 0]
        for _ in range(2):
            psi, subs = W.tick(psi)
            for (j, A, B) in subs:
                F2[j] = F2[j] + net_flux(j, A, B)
        if t % every == 0:
            S["rec"]["t"].append(t)
            S["rec"]["fid"].append(abs(complex(np.vdot(p0, S["psi"])))**2)
            S["rec"]["valley0"].append(float(np.sum(np.abs(g.valley(S["psi"]))**2)))
            S["rec"]["sigma_z"].append(sigma_z(S["psi"]))
            S["rec"]["flux_vphi"].append(g.vphi_flux(F2, ra, 2.0))
            S["rec"]["z"].append(float(chi2_z(S["c"], ra, 2)[0]))
            print(a, t, {k: round(v[-1], 4) for k, v in S["rec"].items() if k != "t"}, round(time.time() - t_start), flush=True)
        flow_tick(S["c"], S["x"], ra, density(psi), F2)
        S["psi"] = psi
        S["t"] = t + 2
    if S["t"] < g.T:
        pickle.dump(S, open(f, "wb")); print(a, "run checkpoint", S["t"], flush=True); sys.exit(3)
    out = dict(a=a, n=g.n, ticks=g.T, dur=g.dur, rec=S["rec"], windows=S["win"],
               final_fid=abs(complex(np.vdot(p0, S["psi"])))**2,
               final_valley0=float(np.sum(np.abs(g.valley(S["psi"]))**2)),
               final_z=float(chi2_z(S["c"], density(S["psi"]), 2)[0]))
    json.dump(out, open(os.path.join(OUT, f"surf_f4_run_{tag(a)}.json"), "w"), indent=1, default=float)
    os.remove(f) if os.path.exists(f) else None
    print(a, "RUN DONE", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "prep":
        prep(float(sys.argv[2]))
    else:
        run(float(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else T_PHYS)   # phase 5: a = 0.2, t_phys = 600

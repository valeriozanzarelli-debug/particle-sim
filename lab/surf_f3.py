"""surf_f3.py - PHASE 3 of SURFISTI_ESATTI: do the doubler valleys kill the current of a polarised 1s?
Exams F1-F6 in lab/out/surf_f3_prereg.json.
usage: python surf_f3.py check              -> identity check on a 4^3 random-potential lattice (development, declared)
       python surf_f3.py prep learned|bare  -> filtered 1s (2 listening passes of 2048 ticks), projected on valley k = 0
       python surf_f3.py run learned|bare   -> 4000 ticks, 2e5 exact surfers (two-tick flow), checkpointed
       python surf_f3.py judge              -> lab/out/surf_f3_results.json, lab/out/surf_f3_figure.png
"""
import json
import os
import pickle
import sys
import time
import numpy as np
from surf_core import Walk3D, net_flux, SEED, OUT, born_sample, chi2_z, flow_tick, density
from surf_f2_valli import seed, valley_split, sigma_z, CK, N3, A3, X, Y, Z, R, ZA
from surf_f2_kramers import Theta, W2 as W2flux, vphi_from_flux

S_LEARNED = (0.554593, 1.003549)          # NUCLEO_1S sec. 13.2 (learn2, Z = 60, a = 0.5)
R0, R1 = np.sqrt(3) / 2 * A3, np.sqrt(11) / 2 * A3
T_RUN = 4000
WIN = (500, 1000, 2000, 3000, 3990)
NS = 200000


def potential(kind):
    V = -ZA / R
    if kind == "learned":
        V = np.where(np.abs(R - R0) < 1e-6 * A3, S_LEARNED[0] * V, V)
        V = np.where(np.abs(R - R1) < 1e-6 * A3, S_LEARNED[1] * V, V)
    return V


def walk(kind):
    return Walk3D(N3, A3, V=potential(kind) * A3)


def P0(psi):
    return valley_split(psi)["0"]["comp"]


def check():
    """F_tick2(psi) = -F_tick1(Theta W2 psi) for ANY psi; and F1(psi) = F1(Theta psi) only inside one Kramers doublet"""
    n, th = 4, 0.37
    rng = np.random.default_rng(SEED + 11)
    V = rng.normal(size=(n, n, n))
    W = Walk3D(n, th, V=V)
    err = 0.0
    for _ in range(10):
        x = rng.normal(size=(4, n, n, n)) + 1j * rng.normal(size=(4, n, n, n)); x /= np.linalg.norm(x)
        y, _, F1, F2 = W2flux(W, x)
        _, _, G1, _ = W2flux(W, Theta(y, th, V))
        err = max(err, max(float(np.abs(F2[j] + G1[j]).max()) for j in range(3)))
    D = 4 * n**3
    M = np.zeros((D, D), complex)
    for k in range(D):
        e = np.zeros(D, complex); e[k] = 1
        M[:, k] = W2flux(W, e.reshape(4, n, n, n))[0].ravel()
    lam, U = np.linalg.eig(M)
    o = np.argsort(np.angle(lam))
    ii, jj, kk = np.meshgrid(*(np.arange(n),) * 3, indexing="ij")
    col = [((ii + jj) % 2 == u) & ((jj + kk) % 2 == v) for u in (0, 1) for v in (0, 1)]
    cw = lambda vec: np.array([np.sum(np.abs(vec.reshape(4, n, n, n) * m[None])**2) for m in col])
    a = U[:, o[10]]
    ca = int(np.argmax(cw(a)))
    # a second eigenvector on the SAME colour (colours never exchange flux) with a different eigenphase
    b = next(U[:, q] for q in o if abs(np.angle(lam[q]) - np.angle(lam[o[10]])) > 1e-6 and int(np.argmax(cw(U[:, q]))) == ca)
    x = ((a + b) / np.sqrt(2)).reshape(4, n, n, n)
    _, _, F1, F2 = W2flux(W, x)
    cross = max(float(np.abs(F1[j] + F2[j]).max()) for j in range(3))
    _, _, F1a, _ = W2flux(W, x)
    _, _, F1t, _ = W2flux(W, Theta(x, th, V))
    res = dict(identity_err=err, net_flux_two_doublet_combo=cross,
               F1_psi_minus_F1_Theta_psi=max(float(np.abs(F1a[j] - F1t[j]).max()) for j in range(3)))
    json.dump(res, open(os.path.join(OUT, "surf_f3_check.json"), "w"), indent=1)
    print(res)


def prep(kind, budget=540.0):
    TL = 2048
    f = os.path.join(CK, f"f3_prep_{kind}.pkl")
    t_start = time.time()
    if os.path.exists(f):
        S = pickle.load(open(f, "rb"))
    else:
        s0 = seed()
        S = dict(pas=1, t=0, psi0=s0, psi=s0.copy(), cand=np.arange(0.86, 0.94 + 1e-9, 0.0025), ov=[], peaks=[])
        S["acc"] = np.zeros((len(S["cand"]),) + s0.shape, complex)
    while S["pas"] <= 2 and time.time() - t_start < budget:
        W = walk(kind); W.tick_count = S["t"]
        psi = S["psi"]
        while S["t"] < TL and time.time() - t_start < budget:
            t = S["t"]
            if t % 2 == 0:
                w = np.sin(np.pi * (t + 1.0) / TL)**2
                ph = w * np.exp(1j * S["cand"] * t * A3)
                for k in range(len(S["cand"])):
                    S["acc"][k] += ph[k] * psi
                S["ov"].append(np.vdot(S["psi0"], psi))
            psi, _ = W.tick(psi)
            S["t"] += 1
        S["psi"] = psi
        if S["t"] < TL:
            break
        ov = np.array(S["ov"]); tt = np.arange(0, TL, 2) * A3
        hann = np.sin(np.pi * (np.arange(0, TL, 2) + 1.0) / TL)**2
        fine = np.arange(0.86, 0.94, 0.0001)
        spec = np.abs((ov * hann)[None, :] @ np.exp(1j * np.outer(tt, fine)))[0]**2
        om = float(fine[int(np.argmax(spec))])
        kc = int(np.argmin(np.abs(S["cand"] - om)))
        v = S["acc"][kc] / np.linalg.norm(S["acc"][kc])
        S["peaks"].append(om)
        print(kind, "pass", S["pas"], "peak", om, flush=True)
        if S["pas"] == 1:
            S.update(pas=2, t=0, psi0=v, psi=v.copy(), cand=np.array([om]), ov=[])
            S["acc"] = np.zeros((1,) + v.shape, complex)
        else:
            vs = valley_split(v)
            p = P0(v); p = p / np.linalg.norm(p)
            info = dict(peaks=S["peaks"], filtered=dict(sigma_z=sigma_z(v), valleys={k: x["weight"] for k, x in vs.items()}),
                        prepared=dict(sigma_z=sigma_z(p), valley0=valley_split(p)["0"]["weight"]))
            np.save(os.path.join(CK, f"f3_psi0_{kind}.npy"), p)
            json.dump(info, open(os.path.join(OUT, f"surf_f3_prep_{kind}.json"), "w"), indent=1, default=float)
            print(kind, "PREP DONE", info, flush=True)
            os.remove(f) if os.path.exists(f) else None
            return
    pickle.dump(S, open(f, "wb")); print(kind, "prep checkpoint pass", S["pas"], "tick", S["t"], flush=True); sys.exit(3)


def pos(c, x):
    return (c + x - N3 / 2) * A3


def run(kind, budget=540.0):
    f = os.path.join(CK, f"f3_run_{kind}.pkl")
    t_start = time.time()
    if os.path.exists(f):
        S = pickle.load(open(f, "rb"))
    else:
        p0 = np.load(os.path.join(CK, f"f3_psi0_{kind}.npy"))
        rng = np.random.default_rng(SEED + 300 + (kind == "bare"))
        c, x = born_sample(rng, density(p0), NS)
        S = dict(t=0, psi=p0.copy(), c=c, x=x, start={}, rec=dict(t=[], fid=[], valley0=[], sigma_z=[], flux_vphi=[], z=[]),
                 win={})
    p0 = np.load(os.path.join(CK, f"f3_psi0_{kind}.npy"))
    W = walk(kind)
    rho_reg = None
    while S["t"] < T_RUN and time.time() - t_start < budget:
        t = S["t"]; psi = S["psi"]
        if t in WIN:
            S["start"][t] = pos(S["c"], S["x"]).copy()
        if t - 10 in WIN:
            p0p = S["start"].pop(t - 10); p1 = pos(S["c"], S["x"])
            rr = np.sqrt(np.sum(p0p**2, 1)); rc = np.sqrt(p0p[:, 0]**2 + p0p[:, 1]**2)
            sel = (rr > 1) & (rr < 6) & (np.abs(p0p[:, 2]) < 0.5 * rr)
            dphi = (np.arctan2(p1[:, 1], p1[:, 0]) - np.arctan2(p0p[:, 1], p0p[:, 0]) + np.pi) % (2 * np.pi) - np.pi
            vphi = rc * dphi / (10 * A3)
            r1 = np.sqrt(np.sum(p1**2, 1))
            S["win"][t - 10] = dict(vphi_over_Za_sin=float(np.mean(vphi[sel] / (rc[sel] / rr[sel])) / ZA),
                                    vr_rms_over_Za=float(np.sqrt(np.mean(((r1 - rr) / (10 * A3))[sel]**2)) / ZA))
            print(kind, "window", t - 10, S["win"][t - 10], flush=True)
        ra = density(psi)
        W.tick_count = t
        F2 = [0, 0, 0]
        for _ in range(2):
            psi, subs = W.tick(psi)
            for (j, A, B) in subs:
                F2[j] = F2[j] + net_flux(j, A, B)
        if t % 100 == 0:
            vs = valley_split(S["psi"])
            S["rec"]["t"].append(t)
            S["rec"]["fid"].append(abs(complex(np.vdot(p0, S["psi"])))**2)
            S["rec"]["valley0"].append(vs["0"]["weight"])
            S["rec"]["sigma_z"].append(sigma_z(S["psi"]))
            S["rec"]["flux_vphi"].append(vphi_from_flux(F2, ra, 2.0))
            S["rec"]["z"].append(float(chi2_z(S["c"], ra, 2)[0]))
            print(kind, t, {k: round(v[-1], 4) for k, v in S["rec"].items() if k != "t"}, round(time.time() - t_start), flush=True)
        flow_tick(S["c"], S["x"], ra, density(psi), F2)
        S["psi"] = psi
        S["t"] = t + 2
    if S["t"] < T_RUN:
        pickle.dump(S, open(f, "wb")); print(kind, "run checkpoint", S["t"], flush=True); sys.exit(3)
    out = dict(rec=S["rec"], windows=S["win"], final_fid=abs(complex(np.vdot(p0, S["psi"])))**2,
               final_valley0=valley_split(S["psi"])["0"]["weight"], final_z=float(chi2_z(S["c"], density(S["psi"]), 2)[0]))
    json.dump(out, open(os.path.join(OUT, f"surf_f3_run_{kind}.json"), "w"), indent=1, default=float)
    os.remove(f) if os.path.exists(f) else None
    print(kind, "RUN DONE", flush=True)


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "check":
        check()
    elif m == "prep":
        prep(sys.argv[2])
    elif m == "run":
        run(sys.argv[2])

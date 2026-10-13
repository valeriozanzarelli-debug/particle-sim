"""surf_f2_valli.py - PHASE 2, part 2: the 4 valleys of the 3-D zig-zag and the spin / circulation of the 1s (exams S1-S5
of lab/out/surf_f2_prereg.json). Atom of the Laboratory: Z = 60, n = 48, a = 0.5.
usage: python surf_f2_valli.py ritz      -> exact levels by filter diagonalisation (S1)       [checkpointed]
       python surf_f2_valli.py smooth    -> 1s listened with the smeared nucleus (S4 states)  [checkpointed]
       python surf_f2_valli.py analyze   -> valley weights / spins (S2, S3), surfers (S4, S5) -> lab/out/surf_f2_valli.json
"""
import json
import pickle
import sys
import time
import numpy as np
from scipy.special import erf
from surf_core import *
from zz_core import PAULI

N3, A3 = 48, 0.5
ZA = 60 * ALPHA
CK = os.path.join(OUT, "surf_ckpt")
x1 = (np.arange(N3) - N3 / 2 + 0.5) * A3
X, Y, Z = np.meshgrid(x1, x1, x1, indexing="ij")
R = np.sqrt(X**2 + Y**2 + Z**2)
SIG = PAULI


def walk(smooth=False):
    V = -ZA * (erf(R / 1.0) if smooth else 1.0) / R
    return Walk3D(N3, A3, V=V * A3)


def seed():
    om, chi, rw = radial_state(-1, 1, A3, ZA, 40.0)
    p = spinor3d_from_radial(rw, chi, -1, X, Y, Z)
    return p / np.linalg.norm(p)


# ---------------------------------------------------------------- S1: Ritz
def ritz(budget=540.0):
    t_start = time.time()
    f = os.path.join(CK, "f2_ritz.pkl")
    T = 1200
    cand = np.arange(0.82, 0.94 + 1e-9, 0.005)
    W = walk()
    if os.path.exists(f):
        S = pickle.load(open(f, "rb")); W.tick_count = S["t"]
    else:
        p = seed()
        S = dict(t=0, psi=p, acc=np.zeros((len(cand),) + p.shape, complex))
    psi, acc = S["psi"], S["acc"]
    while S["t"] < T and time.time() - t_start < budget:
        t = S["t"]
        if t % 2 == 0:
            w = np.sin(np.pi * (t + 1.0) / T)**2
            ph = w * np.exp(1j * cand * t * A3)
            for k in range(len(cand)):
                acc[k] += ph[k] * psi
        psi, _ = W.tick(psi)
        S["t"] += 1
    S["psi"] = psi
    if S["t"] < T:
        pickle.dump(S, open(f, "wb")); print("ritz checkpoint", S["t"], flush=True); sys.exit(3)
    # orthonormal basis, U2 restricted
    B = acc.reshape(len(cand), -1).T
    Q, _ = np.linalg.qr(B)
    UQ = np.zeros_like(Q)
    for k in range(Q.shape[1]):
        Wk = walk()
        v = Q[:, k].reshape(4, N3, N3, N3)
        v, _ = Wk.tick(v); v, _ = Wk.tick(v)
        UQ[:, k] = v.ravel()
    Hm = Q.conj().T @ UQ
    lam, C = np.linalg.eig(Hm)
    om = -np.angle(lam) / (2 * A3)
    o = np.argsort(om)
    levels = []
    s0 = seed().ravel()
    for i in o:
        v = Q @ C[:, i]; v /= np.linalg.norm(v)
        Wk = walk(); u, _ = Wk.tick(v.reshape(4, N3, N3, N3)); u, _ = Wk.tick(u)
        l2 = np.vdot(v, u.ravel())
        resid = float(np.linalg.norm(u.ravel() - l2 * v))
        levels.append(dict(omega=float(om[i]), resid=resid, w_seed=float(abs(np.vdot(v, s0))**2)))
        if 0.83 < om[i] < 0.91 and resid < 1e-4:
            np.save(os.path.join(CK, f"f2_level_{om[i]:.5f}.npy"), v.reshape(4, N3, N3, N3))
    json.dump(levels, open(os.path.join(CK, "f2_ritz_levels.json"), "w"), indent=1)
    os.remove(f) if os.path.exists(f) else None
    for L in levels:
        print(L)


# ---------------------------------------------------------------- S4 state: listening with the smeared nucleus
def smooth(budget=540.0):
    t_start = time.time()
    TL = 4096
    f = os.path.join(CK, "f2_smooth.pkl")
    if os.path.exists(f):
        S = pickle.load(open(f, "rb"))
    else:
        S = dict(pas=1, t=0, psi0=seed(), psi=seed(), acc=None, cand=np.arange(0.86, 0.94 + 1e-9, 0.0025), ov=[])
        S["acc"] = np.zeros((len(S["cand"]),) + S["psi0"].shape, complex)
    while S["pas"] <= 2 and time.time() - t_start < budget:
        W = walk(True); W.tick_count = S["t"]
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
        om = fine[int(np.argmax(spec))]
        kc = int(np.argmin(np.abs(S["cand"] - om)))
        v = S["acc"][kc] / np.linalg.norm(S["acc"][kc])
        print("smooth pass", S["pas"], "peak", om, flush=True)
        if S["pas"] == 1:
            S.update(pas=2, t=0, psi0=v, psi=v.copy(), cand=np.array([om]), ov=[])
            S["acc"] = np.zeros((1,) + v.shape, complex)
        else:
            np.save(os.path.join(CK, "f2_smooth_1s.npy"), v)
            S["pas"] = 3
            os.remove(f) if os.path.exists(f) else None
            print("smooth DONE", flush=True)
            return
    pickle.dump(S, open(f, "wb")); print("smooth checkpoint pass", S["pas"], "tick", S["t"], flush=True); sys.exit(3)


# ---------------------------------------------------------------- analysis
VALLEYS = {"0": (0, 0, 0), "xy": (np.pi / 2, np.pi / 2, 0), "xz": (np.pi / 2, 0, np.pi / 2), "yz": (0, np.pi / 2, np.pi / 2)}


def valley_split(psi):
    k = 2 * np.pi * np.fft.fftfreq(N3)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    best, lab = None, None
    for name, q in VALLEYS.items():
        dmin = None
        for m in np.ndindex(2, 2, 2):
            qq = np.array(q) + np.pi * np.array(m)
            d = sum(np.abs((Kc - qc + np.pi) % (2 * np.pi) - np.pi)**2 for Kc, qc in zip((KX, KY, KZ), qq))
            dmin = d if dmin is None else np.minimum(dmin, d)
        if best is None:
            best, lab = dmin, np.zeros(dmin.shape, int)
        else:
            better = dmin < best
            lab[better] = list(VALLEYS).index(name); best = np.minimum(best, dmin)
    pk = np.fft.fftn(psi, axes=(1, 2, 3))
    out = {}
    for i, name in enumerate(VALLEYS):
        comp = np.fft.ifftn(pk * (lab == i)[None], axes=(1, 2, 3))
        w = float(np.sum(np.abs(comp)**2))
        sz = sum(float(np.real(np.sum(np.einsum("a...,ab,b...->...", comp[q:q + 2].conj(), SIG[2], comp[q:q + 2])))) for q in (0, 2))
        out[name] = dict(weight=w, sigma_z_lab=sz / max(w, 1e-300), comp=comp)
    return out


def sigma_z(psi):
    return sum(float(np.real(np.sum(np.einsum("a...,ab,b...->...", psi[q:q + 2].conj(), SIG[2], psi[q:q + 2])))) for q in (0, 2)) / float(np.sum(np.abs(psi)**2))


def surfers(W, psi, T, rng, Ns=200000):
    """two-tick exact flow, like E5: z checks, v_phi / v_r over ticks 500-510, r drift growth 100 -> T"""
    rho = density(psi)
    c, x = born_sample(rng, rho, Ns)
    r0 = np.sqrt(np.sum(((c + x - N3 / 2) * A3)**2, 1))
    zs, drift = [], {}
    p500 = None
    out = {}
    for t in range(0, T + 1, 2):
        if t % 50 == 0:
            zs.append(chi2_z(c, density(psi), 2)[0])
        if t in (100, T):
            r = np.sqrt(np.sum(((c + x - N3 / 2) * A3)**2, 1)); drift[t] = float(np.sqrt(np.mean((r - r0)**2)))
        if t == 500:
            p500 = (c + x - N3 / 2) * A3
        if t == 510:
            p1 = (c + x - N3 / 2) * A3; p0 = p500
            rr = np.sqrt(np.sum(p0**2, 1)); rc = np.sqrt(p0[:, 0]**2 + p0[:, 1]**2)
            sel = (rr > 1) & (rr < 6) & (np.abs(p0[:, 2]) < 0.5 * rr)
            dphi = (np.arctan2(p1[:, 1], p1[:, 0]) - np.arctan2(p0[:, 1], p0[:, 0]) + np.pi) % (2 * np.pi) - np.pi
            vphi = rc * dphi / (10 * A3)
            r1 = np.sqrt(np.sum(p1**2, 1))
            out["vphi_over_Za_sin"] = float(np.mean(vphi[sel] / (rc[sel] / rr[sel])) / ZA)
            out["vr_rms_over_Za"] = float(np.sqrt(np.mean(((r1 - rr) / (10 * A3))[sel]**2)) / ZA)
        if t == T:
            break
        ra = density(psi); F2 = [0, 0, 0]
        for _ in range(2):
            psi, subs = W.tick(psi)
            for (j, A, B) in subs:
                F2[j] = F2[j] + net_flux(j, A, B)
        flow_tick(c, x, ra, density(psi), F2)
    out.update(max_abs_z=float(np.max(np.abs(zs))), drift=drift,
               drift_growth=drift[T] / max(drift[100], 1e-300) if 100 in drift and T in drift else None)
    return out


def flux_vphi(W, psi):
    """azimuthal velocity of the two-tick flux field itself (cell-centred), 1 < r < 6, |cos| < 0.5"""
    r0 = density(psi); F2 = [0, 0, 0]
    for _ in range(2):
        psi, subs = W.tick(psi)
        for (j, A, B) in subs:
            F2[j] = F2[j] + net_flux(j, A, B)
    Fc = [0.5 * (F2[j] + np.roll(F2[j], 1, axis=j)) / 2 for j in range(3)]
    rc = np.sqrt(X**2 + Y**2)
    vphi = (-Y * Fc[0] + X * Fc[1]) / np.maximum(rc, 1e-9) / r0
    m = (R > 1) & (R < 6) & (np.abs(Z) < 0.5 * R)
    return float(np.sum(r0[m] * vphi[m] / (rc[m] / R[m])) / np.sum(r0[m]) / ZA)


def analyze():
    rng = np.random.default_rng(SEED + 200)
    levels = json.load(open(os.path.join(CK, "f2_ritz_levels.json")))
    res = dict(ritz_levels=[L for L in levels if 0.82 < L["omega"] < 0.94])
    files = sorted(f for f in os.listdir(CK) if f.startswith("f2_level_"))
    states = {f[9:-4]: np.load(os.path.join(CK, f)) for f in files}
    states["E5_listened"] = np.load(os.path.join(CK, "state_1s.npy"))
    states["radial_seed"] = seed()
    sm = os.path.join(CK, "f2_smooth_1s.npy")
    if os.path.exists(sm):
        states["smooth_nucleus_1s"] = np.load(sm)
    res["states"] = {}
    for name, p in states.items():
        p = p / np.linalg.norm(p)
        vs = valley_split(p)
        d = dict(sigma_z=sigma_z(p), valleys={k: dict(weight=v["weight"], sigma_z_lab=v["sigma_z_lab"]) for k, v in vs.items()})
        d["sum_w_sigma"] = sum(v["weight"] * v["sigma_z_lab"] for v in vs.values())
        d["flux_vphi_over_Za_sin"] = flux_vphi(walk(name == "smooth_nucleus_1s"), p)
        # the k = 0 valley alone, instantaneous flux (not a level: only indicative)
        p0 = vs["0"]["comp"] / np.linalg.norm(vs["0"]["comp"])
        d["valley0_only_flux_vphi"] = flux_vphi(walk(name == "smooth_nucleus_1s"), p0)
        res["states"][name] = d
        print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in d.items() if k != "valleys"},
              {k: (round(v["weight"], 3), round(v["sigma_z_lab"], 3)) for k, v in d["valleys"].items()}, flush=True)
    # S5: E5 on the exact level nearest 0.8996; S4: surfers with the smooth nucleus
    ex = min((k for k in states if k[0].isdigit()), key=lambda k: abs(float(k) - 0.8996))
    res["S5_level"] = ex
    res["S5"] = surfers(walk(), states[ex] / np.linalg.norm(states[ex]), 1000, rng)
    print("S5", res["S5"], flush=True)
    if "smooth_nucleus_1s" in states:
        res["S4"] = surfers(walk(True), states["smooth_nucleus_1s"], 520, rng)
        print("S4", res["S4"], flush=True)
    json.dump(res, open(os.path.join(OUT, "surf_f2_valli.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    {"ritz": ritz, "smooth": smooth, "analyze": analyze}[sys.argv[1]]()

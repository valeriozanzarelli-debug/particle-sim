"""surf_states.py - the 3-D states of the zig-zag atom used by the exams (Z = 60, n = 48, a = 0.5: the LIVE atom of
src/atomo/zz.rs). Each lattice level is obtained by LISTENING to the rule (Hann-windowed lock-in, as Rust ascolta),
twice: pass 1 from the radial state put on the cubic lattice, on a grid of candidate frequencies (step 0.0025) +
the overlap spectrum; pass 2 from the pass-1 wave at the refined frequency. Checkpointed (run again until it prints
DONE). Writes lab/out/surf_ckpt/state_<name>.npy (not in git) and lab/out/surf_states.json.
usage: python surf_states.py 1s|2p
"""
import json
import sys
import time
import numpy as np
from surf_core import *

N3, A3 = 48, 0.5
ZA = 60 * ALPHA
TL = 4096                       # ticks per listening pass (2048 hbar/mc^2)
CK = os.path.join(OUT, "surf_ckpt")
os.makedirs(CK, exist_ok=True)
SPEC = {"1s": dict(kappa=-1, n=1, lo=0.885, hi=0.915), "2p": dict(kappa=1, n=2, lo=0.93, hi=0.995)}


def atom_walk():
    x = (np.arange(N3) - N3 / 2 + 0.5) * A3
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    return Walk3D(N3, A3, V=-ZA / R * A3), (X, Y, Z)


def seed_state(name):
    sp = SPEC[name]
    W, (X, Y, Z) = atom_walk()
    om, chi, rw = radial_state(sp["kappa"], sp["n"], A3, ZA, 40.0)
    psi = spinor3d_from_radial(rw, chi, sp["kappa"], X, Y, Z)
    return psi / np.linalg.norm(psi), om


def run_pass(name, pas, budget=540.0):
    t_start = time.time()
    sp = SPEC[name]
    f = os.path.join(CK, f"listen_{name}_{pas}.npz")
    W, _ = atom_walk()
    if os.path.exists(f):
        d = np.load(f)
        t0, psi, acc, cand, ov, psi0 = int(d["t"]), d["psi"], d["acc"], d["cand"], list(d["ov"]), d["psi0"]
    else:
        if pas == 1:
            psi0, _ = seed_state(name)
            cand = np.arange(sp["lo"], sp["hi"] + 1e-9, 0.0025)
        else:
            p1 = np.load(os.path.join(CK, f"state_{name}_p1.npy"))
            psi0 = p1
            info = json.load(open(os.path.join(CK, f"info_{name}_p1.json")))
            cand = np.array([info["omega"]])
        psi = psi0.copy()
        acc = np.zeros((len(cand),) + psi0.shape, complex)
        t0, ov = 0, []
    W.tick_count = t0
    t = t0
    while t < TL and time.time() - t_start < budget:
        # only EVEN ticks: the rule's period is two ticks (axis order xyz, then zyx); an eigenstate of the two-tick
        # operator U2 = W_zyx W_xyz is recovered by summing psi(2m) e^{i omega 2m a}
        if t % 2 == 0:
            w = np.sin(np.pi * (t + 1.0) / TL)**2
            ph = np.exp(1j * cand * t * A3) * w
            for k in range(len(cand)):
                acc[k] += ph[k] * psi
            ov.append(np.vdot(psi0, psi))
        psi, _ = W.tick(psi)
        t += 1
    if t < TL:
        np.savez(f, t=t, psi=psi, acc=acc, cand=cand, ov=np.array(ov), psi0=psi0)
        print(f"{name} pass {pas}: checkpoint at tick {t}/{TL}", flush=True)
        return False
    ov = np.array(ov)
    tt = np.arange(0, TL, 2) * A3
    hann = np.sin(np.pi * (np.arange(0, TL, 2) + 1.0) / TL)**2
    fine = np.arange(sp["lo"], sp["hi"], 0.0001)
    spec = np.abs((ov * hann)[None, :] @ np.exp(1j * np.outer(tt, fine)))[0]**2
    kb = int(np.argmax(spec[1:-1])) + 1
    y0, y1, y2 = spec[kb - 1], spec[kb], spec[kb + 1]
    om = fine[kb] + 0.0001 * 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    kc = int(np.argmin(np.abs(cand - om)))
    v = acc[kc]
    v /= np.linalg.norm(v)
    # purity: residual of the eigen-equation of the TWO-tick operator (starting at an even tick)
    Wt, _ = atom_walk()
    v1, _ = Wt.tick(v)
    v1, _ = Wt.tick(v1)
    lam = np.vdot(v, v1)
    resid = float(np.linalg.norm(v1 - lam * v))
    om_ray = float(-np.angle(lam) / (2 * A3))
    np.save(os.path.join(CK, f"state_{name}_p{pas}.npy"), v)
    info = dict(omega=om_ray, omega_peak=float(om), cand=float(cand[kc]), resid_two_ticks=resid,
                overlap_with_start=float(np.abs(np.vdot(psi0, v))**2), ticks=TL)
    json.dump(info, open(os.path.join(CK, f"info_{name}_p{pas}.json"), "w"), indent=1)
    os.remove(f) if os.path.exists(f) else None
    print(f"{name} pass {pas} DONE:", info, flush=True)
    return True


if __name__ == "__main__":
    name = sys.argv[1]
    for pas in (1, 2):
        if os.path.exists(os.path.join(CK, f"state_{name}_p{pas}.npy")):
            continue
        if not run_pass(name, pas):
            sys.exit(3)
    v = np.load(os.path.join(CK, f"state_{name}_p2.npy"))
    np.save(os.path.join(CK, f"state_{name}.npy"), v)
    allinfo = {}
    jf = os.path.join(OUT, "surf_states.json")
    if os.path.exists(jf):
        allinfo = json.load(open(jf))
    allinfo[name] = {p: json.load(open(os.path.join(CK, f"info_{name}_{p}.json"))) for p in ("p1", "p2")}
    json.dump(allinfo, open(jf, "w"), indent=1)
    print("DONE", name)

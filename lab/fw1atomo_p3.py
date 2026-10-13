"""fw1atomo_p3.py - P3 of the pre-registration: the 3-D atom (equivalent of lab/zz_atom3d.py, ZIGZAG A1) with FW1 and
with the frozen rule.  Z = 60, nucleus on the corner of 8 cells, box ~16 lambda_C, a = 0,4 / 0,3 / 0,24 / 0,2.
Steps:  ritz_<rule>_<start>_<a>   Ritz in [0,78; 0,96] (dw 0,01, 600 lambda_C/c) from start 'dirac' (continuum
                                   Dirac 1s, spin up) or 'cieca' (random spinors on the 4^3 sites around the nucleus)
        v0                         validation: F0 with the zz_atom3d protocol at a = 0,4, box 24 (projected blob,
                                   Gaussian lock-in std 27, Rayleigh quotient) -> 0,130043
Run: OMP_NUM_THREADS=1 timeout 590 python lab/fw1atomo_p3.py <step>  (exit 3 = checkpoint); list; sum"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import fw1atomo_core as fc
import zz_core as zz

T0 = time.time()
OUT = os.path.join(os.path.dirname(__file__), "out", "fw1atomo_p3.json")
RES = os.path.join(fc.SCR, "p3")
os.makedirs(RES, exist_ok=True)
AS = {0.4: 40, 0.3: 54, 0.24: 66, 0.2: 80}
ZA = fc.ZA60
B_DIRAC = 1 - zz.dirac_formula(1, -1, ZA)
STEPS = ["v0"] + [f"ritz_{r}_{s}_{a}_c{c}" for a in AS for r in ("FW1", "F0") for s in ("dirac", "cieca") for c in ((0, 1) if s == "dirac" else (0,))]


def run_ritz(step):
    _, rule, start, a, c = step.split("_")
    a, c = float(a), int(c[1:])
    L = fc.Lattice(AS[a], a, "angolo")
    R = L.rule(rule, 0.4)
    S = L.dirac1s((1, 0)) if start == "dirac" else L.blind()
    S = S * (L.cls == c)[None]                    # amendment A3: one period-sublattice (U2 is block diagonal)
    S = [S / np.linalg.norm(S)]
    tp = time.time()
    for _ in range(3):
        R.period(S[0])
    t_per = (time.time() - tp) / 3
    lv = fc.ritz(R, S, a, 0.78, 0.96, 0.01, int(round(600 / (2 * a))), ck=os.path.join(fc.SCR, "p3_" + step), tmax=540, t0=T0)
    out = []
    for o in lv:
        if o["resid"] < 1e-4 and o["w_start"][0] >= 1e-4:
            d = L.diag(o["v"])
            d.update(omega=o["omega"], resid=o["resid"], w=o["w_start"][0])
            out.append(d)
    return {"levels": out, "s_per_period": t_per, "n": AS[a], "a": a}


def run_v0():
    """zz_atom3d.run(0.4) with the rules of fw1atomo_core (F0): same blob, projector, window, Rayleigh"""
    a, half, S_T = 0.4, 12.0, 27.0
    n = int(round(2 * half / a))
    L = fc.Lattice(n, a, "angolo")
    R = L.rule("F0")
    r = L.r
    psi = np.exp(-r * ZA)[None] * (np.array([1, 0, 1, 0], complex) / np.sqrt(2))[:, None, None, None]
    pk = np.fft.fftn(psi, axes=(1, 2, 3)).reshape(4, -1)
    kk = 2 * np.pi * np.fft.fftfreq(n)
    K = np.stack([q.ravel() for q in np.meshgrid(kk, kk, kk, indexing="ij")], -1)
    out = np.zeros_like(pk)
    for b0 in range(0, K.shape[0], 200000):
        sl = slice(b0, b0 + 200000)
        w, v = np.linalg.eig(zz.bloch_batch(K[sl], a, 0.5))
        om = -np.angle(w)
        o = np.argsort(-om, axis=1)[:, :2]
        Q, _ = np.linalg.qr(np.take_along_axis(v, o[:, None, :], axis=2))
        c = np.einsum("nab,na->nb", Q.conj(), pk[:, sl].T)
        out[:, sl] = np.einsum("nab,nb->na", Q, c).T
    psi = np.fft.ifftn(out.reshape(4, n, n, n), axes=(1, 2, 3))
    psi /= np.sqrt(np.sum(np.abs(psi)**2))

    def rayleigh(p):
        q = R.tick(R.tick(p.copy(), 0), 1)
        return -np.angle(np.vdot(p, q) / np.vdot(p, p)) / (2 * a)
    om0 = rayleigh(psi)
    nt = int(round(4.5 * S_T / a))
    nt += nt % 2
    c0, s = nt / 2, S_T / a
    acc = np.zeros_like(psi)
    for t in range(nt + 1):
        acc += np.exp(-0.5 * ((t - c0) / s)**2) * np.exp(1j * om0 * a * t) * psi
        psi = R.tick(psi, t)
    acc /= np.linalg.norm(acc)
    om1 = rayleigh(acc)
    return {"omega_start": float(om0), "omega": float(om1), "binding": float(1 - om1), "ref_zz_atom3d": 0.13004250114978633}


def run(step):
    fn = os.path.join(RES, step + ".json")
    if os.path.exists(fn):
        print("done already", step)
        return
    res = run_v0() if step == "v0" else run_ritz(step)
    json.dump(res, open(fn, "w"), indent=1)
    if step != "v0":
        for d in res["levels"]:
            if d["resid"] < 1e-6 and d["w"] > 2e-3:
                print(f"{d['omega']:.5f} res {d['resid']:.1e} w {d['w']:.4f} G {d['val_Gamma']:.3f} N {d['val_N']:.3f} "
                      f"pi {d['val_pi']:.3f} r {d['r_mean']:.2f} S {np.round(d['spin'], 2)}", flush=True)
    else:
        print(res)
    print("DONE", step, time.time() - T0)


if __name__ == "__main__":
    if sys.argv[1] == "list":
        print(" ".join(STEPS))
    else:
        run(sys.argv[1])

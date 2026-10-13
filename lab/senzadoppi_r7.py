"""senzadoppi_r7.py - amendment A3 diagnostic of E7 at theta = 0,25: Z = 60, 64^3 cells of a = 0,25 lambda_C (box 16
lambda_C ~ 7 Bohr radii), nucleus on the corner of 8 cells, continuum Dirac 1s (spin up) as TEST initial state,
3200 periods.  Pass A: survival c(t) -> lines.  Pass B: same evolution, states filtered at the 4 strongest lines ->
valley content (Fourier weight within 0,3 of Gamma + {0,pi}^3, of the N valleys (pi/2,pi/2,0)-type and of the pi
valleys (pi/2,0,0)/(pi/2,pi/2,pi/2)-type, all modulo the copies) and spin <sigma_z>.
Run: OMP_NUM_THREADS=2 timeout 590 python lab/senzadoppi_r7.py <F0|PAW|FW1>   (rerun until it prints DONE)
     python lab/senzadoppi_r7.py sum   -> lab/out/senzadoppi_r7.json"""
import os
import sys
import json
import time
import itertools
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import zz_core as zz
import senzadoppi_real as sr
import senzadoppi_r5  # noqa: F401  registers FW1
from senzadoppi_r2 import REAL
from senzadoppi_r4 import lines

SCR = "/tmp/claude-0/-home-claude-particle-sim/9de435bc-26c2-5d88-b70d-07efb8033e5b/scratchpad"
OUT = os.path.join(os.path.dirname(__file__), "out", "senzadoppi_r7.json")
Z, A, N, NPER = 60, 0.25, 64, 3200
ZA = Z * zz.ALPHA
T0 = time.time()


def initial():
    x = (np.arange(N) - N / 2 + 0.5) * A
    X, Y, Zc = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Zc**2)
    gam = np.sqrt(1 - ZA**2)
    g = r**(gam - 1) * np.exp(-ZA * r)
    f = -np.sqrt((1 - gam) / (1 + gam)) * g
    nx, ny, nz = X / r, Y / r, Zc / r
    u = np.array([g, 0 * g])
    dd = 1j * f * np.array([nz, nx + 1j * ny])
    psi = np.concatenate([(u - dd) / np.sqrt(2), (u + dd) / np.sqrt(2)]).astype(complex)
    return psi / np.sqrt(np.sum(np.abs(psi)**2)), r


def valley_content(phi):
    F = np.fft.fftn(phi, axes=(1, 2, 3))
    P = np.sum(np.abs(F)**2, axis=0)
    P /= P.sum()
    k = np.fft.fftfreq(N) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")

    def ball(c):
        m = np.zeros_like(P, bool)
        for G in itertools.product((0, np.pi), repeat=3):
            d2 = 0
            for Kc, cc, gg in ((KX, c[0], G[0]), (KY, c[1], G[1]), (KZ, c[2], G[2])):
                d2 = d2 + ((Kc - cc - gg + np.pi) % (2 * np.pi) - np.pi)**2
            m |= d2 < 0.3**2
        return m
    h = np.pi / 2
    groups = {"Gamma": [(0, 0, 0)], "N": [(h, h, 0), (h, 0, h), (0, h, h)],
              "pi": [(h, 0, 0), (0, h, 0), (0, 0, h), (h, h, h)]}
    out = {}
    for gname, cs in groups.items():
        m = np.zeros_like(P, bool)
        for c in cs:
            m |= ball(c)
        out[gname] = float(P[m].sum())
    L, R = phi[:2], phi[2:]
    sz = float(np.real(np.sum(np.conj(L[0]) * L[0] - np.conj(L[1]) * L[1] + np.conj(R[0]) * R[0] - np.conj(R[1]) * R[1])))
    out["sigma_z"] = sz / float(np.sum(np.abs(phi)**2))
    return out


def run(name):
    ck = os.path.join(SCR, f"senzadoppi_r7_{name}_ck.json")
    st = json.load(open(ck)) if os.path.exists(ck) else {"pass": "A", "p": 0}
    ticks, mode, h, vecs = REAL[name]
    psi0, r = initial()
    RR = sr.RealRule(ticks, A, mode, h, vecs, V=-ZA / r * A)
    T = len(ticks) * A
    fpsi = os.path.join(SCR, f"senzadoppi_r7_{name}_psi.npy")
    fc = os.path.join(SCR, f"senzadoppi_r7_{name}_c.npy")
    facc = os.path.join(SCR, f"senzadoppi_r7_{name}_acc.npy")
    if st["pass"] == "done":
        print("DONE")
        return
    psi = np.load(fpsi) if st["p"] > 0 else psi0.copy()
    c = np.load(fc) if (st["p"] > 0 or st["pass"] == "B") else np.zeros(NPER + 1, complex)
    w = np.hanning(NPER + 1)
    if st["pass"] == "B":
        Es = np.array(st["E_lines"])
        acc = np.load(facc) if st["p"] > 0 else np.zeros((len(Es),) + psi0.shape, complex)
        if st["p"] == 0:
            acc += w[0] * psi0[None]
    if st["pass"] == "A" and st["p"] == 0:
        c[0] = 1
    while st["p"] < NPER:
        if time.time() - T0 > 520:
            np.save(fpsi, psi)
            np.save(fc, c)
            if st["pass"] == "B":
                np.save(facc, acc)
            json.dump(st, open(ck, "w"))
            print("checkpoint", name, st["pass"], st["p"], flush=True)
            sys.exit(3)
        psi = RR.period(psi)
        st["p"] += 1
        p = st["p"]
        if st["pass"] == "A":
            c[p] = np.vdot(psi0, psi)
        else:
            ph = np.exp(1j * Es * p * T) * w[p]
            acc += ph[:, None, None, None, None] * psi[None]
    if st["pass"] == "A":
        np.save(fc, c)
        _, _, pk = lines(c, T)
        pk = sorted(pk, key=lambda x: -x[1])[:4]
        st = {"pass": "B", "p": 0, "E_lines": [e for e, _ in pk], "W_lines": [wt for _, wt in pk],
              "norm_A": float(np.sum(np.abs(psi)**2))}
        json.dump(st, open(ck, "w"))
        print("pass A done", name, pk, flush=True)
        sys.exit(3)
    vc = [valley_content(acc[i]) for i in range(len(Es))]
    st["valleys"] = vc
    st["pass"] = "done"
    json.dump(st, open(ck, "w"))
    print("DONE", name, json.dumps(st), flush=True)


def summary():
    res = {}
    for name in ("F0", "PAW", "FW1"):
        ck = os.path.join(SCR, f"senzadoppi_r7_{name}_ck.json")
        if not os.path.exists(ck):
            continue
        st = json.load(open(ck))
        c = np.load(os.path.join(SCR, f"senzadoppi_r7_{name}_c.npy"))
        pex = 1 - np.abs(c)**2
        _, _, pk = lines(c, 2 * A)
        pk = sorted(pk, key=lambda x: -x[1])
        main = pk[0]
        wmain = main[1] + sum(wt for e, wt in pk[1:] if abs(e - main[0]) <= 5e-3)
        res[name] = {"max_p_exc": float(pex.max()), "lines_ge_2pc": pk, "main_weight_with_doublet": float(wmain),
                     "filtered_lines": st.get("E_lines"), "valleys": st.get("valleys"),
                     "dirac_1s": float(zz.dirac_formula(1, -1, ZA))}
        print(name, json.dumps(res[name]))
    json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "sum":
        summary()
    else:
        run(sys.argv[1])

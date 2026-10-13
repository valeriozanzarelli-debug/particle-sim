"""emiss_v4_dev.py - DEVELOPMENT (declared in the phase-2 pre-registration): spectra of the radial starts of 1s and 2p
on the 48^3 lattice with the bare nucleus and with the learned nucleus of NUCLEO_1S sec. 13 (s0 = 0.554593 on the 8
sites at r = sqrt3/2 a, s1 = 1.003549 on the 24 sites at r = sqrt11/2 a; Z = 60, a = 0.5; full cube = 4 colours).
Question: does the learned nucleus close the valley multiplets of 1s AND of the 2p?
Run: python lab/emiss_v4_dev.py [bare|learned] [nt]"""
import sys, json, os, time
import numpy as np
from emiss_core import Dirac3D, radial_state, dirac_3d_from_radial, ovl, OUT, ALPHA

S_LEARNED = (0.554593, 1.003549)


def make_atom(n=48, a=0.5, s=(1.0, 1.0)):
    D = Dirac3D(n, a)
    r0, r1 = np.sqrt(3) / 2 * a, np.sqrt(11) / 2 * a
    f = np.ones_like(D.r)
    f[np.abs(D.r - r0) < 1e-6 * a] = s[0]
    f[np.abs(D.r - r1) < 1e-6 * a] = s[1]
    D.pc = np.exp(1j * D.Z * ALPHA / D.r * a * f)
    D.sfac = s
    return D


STARTS = {"1s+": (-1, 1, 0.5), "2p1/2+": (1, 2, 0.5), "2p3/2+": (-2, 2, 0.5), "2p3/2+3": (-2, 2, 1.5),
          "2p1/2-": (1, 2, -0.5), "2p3/2-": (-2, 2, -0.5)}


def start_state(D, st):
    kap, nr, mj = STARTS[st]
    ar = 0.05 if kap == -1 else 0.5
    om, rg, g, f = radial_state(kap, nr, ar, 60 * ALPHA)
    return dirac_3d_from_radial(rg, g, f, kap, mj, D.x)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "learned"
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 2400
    D = make_atom(s=S_LEARNED if mode == "learned" else (1.0, 1.0))
    names = ["1s+", "2p1/2+", "2p3/2+", "2p3/2+3"]
    ps = [start_state(D, s) for s in names]
    p0 = [p.copy() for p in ps]
    C = {s: [] for s in names}
    t0 = time.time()
    for t in range(nt):
        if t % 2 == 0:
            for s, p, q in zip(names, ps, p0):
                C[s].append(ovl(q, p))
        for p in ps:
            D.tick(p, t)
    res = {}
    for s in names:
        c = np.array(C[s])
        tt = np.arange(len(c)) * 2 * D.a
        h = np.sin(np.pi * (np.arange(len(c)) + 0.5) / len(c))**2
        lo, hi = (0.82, 0.92) if s == "1s+" else (0.94, 1.0)
        om = np.arange(lo, hi, 2e-5)
        S = np.abs(np.exp(1j * np.outer(om, tt)) @ (h * c)) / np.sum(h)
        pk = [(float(om[i]), float(S[i])) for i in range(1, len(om) - 1) if S[i] > S[i - 1] and S[i] > S[i + 1] and S[i] > 0.01 * S.max()]
        pk.sort(key=lambda x: -x[1])
        res[s] = pk[:6]
        print(mode, s, [(round(w, 5), round(v, 3)) for w, v in pk[:6]], flush=True)
    json.dump(res, open(os.path.join(OUT, f"emiss_v4_dev_{mode}.json"), "w"), indent=1)
    print(f"({time.time() - t0:.0f} s)")

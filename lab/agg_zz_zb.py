"""agg_zz_zb.py - POST-HOC (written after A3b failed): where is the Zitterbewegung in the zoom?
The stationary 1s has no trembling (it is an energy eigenstate: no + / - energy superposition). The trembling appears
when the fine region holds a superposition of positive and negative energies. Test: the coarse 1s (chain B, a0/16)
prolonged into a larger fine patch (5 coarse cells = 160 fine cells per side; the central cell is untouched by the
patch edge for 2 B = 64 ticks) in two ways: (i) with the small component chi of the walk's own positive-energy
relation (as in A3), (ii) with chi = 0 (a 'naive' interface: the large component alone = half + and half - energy
near the nucleus). Observable: sum over the central cell of 2 Re(L^dag R) = psi^dag beta psi, per fine tick.
Writes lab/out/agg_zz_zb.json / .npz."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import json, pickle, time
import numpy as np
from agg_zz_core import *
from agg_zz_core import GenRule
from agg_zz_atom import Atom, solve

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
TH = A0 / 16 / 32
NC = 128
st = pickle.load(open(os.path.join(OUT, "agg_zz_3d_B.pkl"), "rb"))
L = max([l for l in st["levels"] if l["B"] <= 32], key=lambda l: l["B"])
B = L["B"]
c0 = NC // 2
t0 = time.time()
at = Atom(L, TH, NC, None, centre="centre")
E1, resid, v, nit = solve(at, "1s")
c = v.reshape((NC,) * 3)
print(f"coarse 1s {E1 * zz.MC2_EV:.5f} eV ({time.time() - t0:.0f} s)", flush=True)


def fit_G():
    G = []
    for j in range(3):
        k = np.zeros(3); k[j] = 1e-4
        lam, V = np.linalg.eig(zz.bloch_two_ticks(k, TH))
        U = V[:, np.argsort(np.angle(lam))[:2]]
        G.append(((U[2:] - U[:2]) / np.sqrt(2)) @ np.linalg.inv((U[:2] + U[2:]) / np.sqrt(2)) / 1e-4)
    return G


G = fit_G()
NP = 5 * B
u = (np.arange(NP) + 0.5) / B - 2.5
nodes = np.arange(-4, 5)


def lagrange_matrix(u, nodes):
    P = np.zeros((len(u), len(nodes)))
    for i, x in enumerate(u):
        i0 = int(np.floor(x)) - 1
        idx = [i0 + q for q in range(4)]
        for a in idx:
            w = 1.0
            for bb in idx:
                if bb != a:
                    w *= (x - bb) / (a - bb)
            P[i, list(nodes).index(a)] = w
    return P


PL = lagrange_matrix(u, nodes)
win = c[c0 - 5:c0 + 6, c0 - 5:c0 + 6, c0 - 5:c0 + 6]
lap = sum(np.roll(win, s, a) for a in range(3) for s in (1, -1)) - 6 * win
p = (win - lap / 24.0)[1:-1, 1:-1, 1:-1] / B ** 1.5
f = np.einsum("ia,jb,kc,abc->ijk", PL, PL, PL, p, optimize=True)
phi = np.zeros((2, NP, NP, NP), complex); phi[0] = f
chi = np.zeros_like(phi)
for j in range(3):
    chi += np.einsum("ab,bxyz->axyz", G[j], -1j * np.gradient(phi, axis=1 + j))
xf = np.arange(NP) - NP / 2 + 0.5
X, Y, Z = np.meshgrid(xf, xf, xf, indexing="ij", sparse=True)
PH = np.exp(1j * zz.ALPHA / (np.sqrt(X ** 2 + Y ** 2 + Z ** 2) * TH) * TH)
I0, I1 = 2 * B, 3 * B
res, ser = {}, {}
for name, ch in (("positive_energy_interface", chi), ("naive_interface_chi0", 0 * chi)):
    psi = np.concatenate([(phi - ch) / np.sqrt(2), (phi + ch) / np.sqrt(2)])[None]
    beta = []
    for t in range(2 * B):
        psi = fine3_tick_fast(psi, TH, t) * PH[None, None]
        Li, Ri = psi[0, 0:2, I0:I1, I0:I1, I0:I1], psi[0, 2:4, I0:I1, I0:I1, I0:I1]
        beta.append(float(np.sum(2 * np.real(np.conj(Li) * Ri))))
    beta = np.array(beta)
    bb = beta - np.polyval(np.polyfit(np.arange(len(beta)), beta, 2), np.arange(len(beta)))
    nz = 8192
    sp = np.abs(np.fft.rfft(bb * np.hanning(len(bb)), nz)) ** 2
    om = 2 * np.pi * np.fft.rfftfreq(nz, d=TH)
    band = (om > 0.5) & (om < 6)
    ipk = np.argmax(np.where(band, sp, 0))
    res[name] = dict(peak_omega=float(om[ipk]), rel_osc_amplitude=float(np.std(bb) / np.mean(beta)),
                     ticks=2 * B)
    ser[name] = beta
    print(f"{name}: peak at {om[ipk]:.3f} m c^2/hbar, relative oscillation {np.std(bb) / np.mean(beta):.2e} "
          f"({time.time() - t0:.0f} s)", flush=True)
res["expected_2E"] = 2 * (1 - 0.5 * zz.ALPHA ** 2)
res["note"] = "post hoc, after A3b failed; not a pre-registered exam"
json.dump(res, open(os.path.join(OUT, "agg_zz_zb.json"), "w"), indent=1)
np.savez_compressed(os.path.join(OUT, "agg_zz_zb.npz"), **ser)

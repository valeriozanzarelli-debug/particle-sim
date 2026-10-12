"""agg_zz_zoom.py - A3: the hybrid ZOOM. Learned coarse rule (chain B top level, cell a0/16) everywhere, the fine zig-zag
(4 components, point Coulomb phase every tick, theta_B = a0/512 = 0.268) in the coarse cell that holds the proton.
Interface (one coarse step = B fine ticks):
  prolongation: the 26 coarse cells around the proton cell (the halo, one coarse cell thick = the light cone of B ticks)
                are refilled every coarse step from the coarse field: block means deconvolved to point values
                (p = c - h^2/24 Lap c), tricubic Lagrange interpolation to the fine cells, small component from the fine
                walk's own positive-energy relation chi = sum_j G_j (-i d_j) phi (G_j fitted on the walk's Bloch eigenvectors
                at small k), (L, R) = ((phi - chi), (phi + chi))/sqrt 2;
  restriction:  the proton cell's coarse value = block mean of phi = (L + R)/sqrt 2 over its B^3 fine cells (x B^{3/2}).
The fine interior keeps its own state across coarse steps (only the halo is refilled).
Energy: phase slope of <psi_0|psi(n)> on the coarse field, hybrid vs uniform coarse (same box, same start = the uniform
coarse 1s). Zitterbewegung: spectrum of the fine interior's scalar density sum 2 Re(L^dag R) over fine ticks.
Run: python lab/agg_zz_zoom.py  (checkpointed: lab/out/agg_zz_zoom.npz / .json)"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import json, pickle, sys, time
import numpy as np
from agg_zz_core import *
from agg_zz_atom import Atom, solve, cell_avg_inv_r

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
MC2 = zz.MC2_EV
TH = A0 / 16 / 32
NC = 128
NSTEP = int(sys.argv[1]) if len(sys.argv) > 1 else 48
CK = os.path.join(OUT, "agg_zz_zoom_state.pkl")
t_start = time.time()
BUDGET = 500

st_B = pickle.load(open(os.path.join(OUT, "agg_zz_3d_B.pkl"), "rb"))
L = max([l for l in st_B["levels"] if l["B"] <= 32], key=lambda l: l["B"])
B = L["B"]
h = B * TH
c0 = NC // 2


# ---------- the positive-energy relation of the fine walk: chi = sum_j G_j k_j phi (linear order, fitted) ----------
def fit_G():
    G = []
    for j in range(3):
        k = np.zeros(3); k[j] = 1e-4
        W = zz.bloch_two_ticks(k, TH)
        lam, V = np.linalg.eig(W)
        U = V[:, np.argsort(np.angle(lam))[:2]]          # positive band: omega = -angle/2 > 0 -> angle < 0
        Pl = (U[:2] + U[2:]) / np.sqrt(2)
        Ps = (U[2:] - U[:2]) / np.sqrt(2)
        G.append(Ps @ np.linalg.inv(Pl) / 1e-4)
    return G


G = fit_G()


def lagrange_matrix(u, nodes):
    """cubic Lagrange interpolation matrix from values at integer nodes to points u (4 nearest nodes)"""
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


NP = 3 * B
u = (np.arange(NP) + 0.5) / B - 1.5                      # fine cell centres in coarse-cell units, relative to c0 centre
nodes = np.arange(-3, 4)
PL = lagrange_matrix(u, nodes)


def prolong(c):
    """coarse scalar spin-up field c (NC^3) -> fine 4-spinor on the patch (4, NP, NP, NP)"""
    win = c[c0 - 4:c0 + 5, c0 - 4:c0 + 5, c0 - 4:c0 + 5]
    lap = (np.roll(win, 1, 0) + np.roll(win, -1, 0) + np.roll(win, 1, 1) + np.roll(win, -1, 1) + np.roll(win, 1, 2)
           + np.roll(win, -1, 2) - 6 * win)
    p = (win - lap / 24.0)[1:-1, 1:-1, 1:-1] / B ** 1.5   # point values at nodes -3..3
    f = np.einsum("ia,jb,kc,abc->ijk", PL, PL, PL, p, optimize=True)
    phi = np.zeros((2, NP, NP, NP), complex)
    phi[0] = f
    chi = np.zeros_like(phi)
    for j in range(3):
        dphi = np.gradient(phi, axis=1 + j)               # per fine cell (a_f = 1 in walk units)
        chi += np.einsum("ab,bxyz->axyz", G[j], -1j * dphi)
    return np.concatenate([(phi - chi) / np.sqrt(2), (phi + chi) / np.sqrt(2)])[None]


xf = (np.arange(NP) - NP / 2 + 0.5)                       # fine cells, nucleus at the central corner (NP even)
Xf, Yf, Zf = np.meshgrid(xf, xf, xf, indexing="ij", sparse=True)
Vf = -zz.ALPHA / (np.sqrt(Xf ** 2 + Yf ** 2 + Zf ** 2) * TH)  # m c^2 units, r in lambda_C
PH = np.exp(-1j * Vf * TH)                                 # one tick = a_f = TH (hbar/mc^2)
I0, I1 = B, 2 * B


def restrict(psi):
    phi = (psi[0, 0:2] + psi[0, 2:4]) / np.sqrt(2)
    blk = phi[:, I0:I1, I0:I1, I0:I1]
    return blk.mean(axis=(1, 2, 3)) * B ** 1.5            # (up, down)


# ---------- uniform coarse atom (nucleus at the centre of cell c0) ----------
res = json.load(open(os.path.join(OUT, "agg_zz_zoom.json"))) if os.path.exists(os.path.join(OUT, "agg_zz_zoom.json")) else {}
if os.path.exists(CK):
    S = pickle.load(open(CK, "rb"))
    at = Atom(L, TH, NC, None, centre="centre")
else:
    t0 = time.time()
    at = Atom(L, TH, NC, None, centre="centre")
    E1, resid, v, nit = solve(at, "1s")
    c_init = v.reshape((NC,) * 3)
    c_init = c_init / np.linalg.norm(c_init)
    res["uniform_coarse"] = dict(E_eV=E1 * MC2, resid=resid, iterations=nit, rest_phase=float(at.rest), B=B,
                                 cell_lambdaC=h, box_a0=NC * h / A0, seconds=time.time() - t0)
    print(f"uniform coarse 1s (nucleus at a cell centre): {E1 * MC2:.6f} eV, residual {resid:.1e} ({time.time() - t0:.0f} s)", flush=True)
    psi_f = prolong(c_init)
    S = dict(c0vec=c_init, c_h=c_init.copy(), c_c=c_init.copy(), psi=psi_f, n=0, C_h=[1.0 + 0j], C_c=[1.0 + 0j],
             beta=[], pint=[], beta_halo=[], down=[], t_fine=0)

Ms = at.Mr * np.exp(1j * at.rest)                         # the learned kernel WITH its rest phase (fine region has it too)


def coarse_step(c):
    return at.D * ifftn(Ms * fftn(at.D * c, 3), 3)


while S["n"] < NSTEP and time.time() - t_start < BUDGET:
    t0 = time.time()
    psi = S["psi"]
    # halo refill from the coarse field (interior kept)
    pro = prolong(S["c_h"])
    keep = psi[:, :, I0:I1, I0:I1, I0:I1].copy()
    psi = pro
    psi[:, :, I0:I1, I0:I1, I0:I1] = keep
    for t in range(B):
        psi = fine3_tick(psi, TH, t) * PH[None, None]
        Li, Ri = psi[0, 0:2, I0:I1, I0:I1, I0:I1], psi[0, 2:4, I0:I1, I0:I1, I0:I1]
        S["beta"].append(float(np.sum(2 * np.real(np.conj(Li) * Ri))))
        S["pint"].append(float(np.sum(np.abs(Li) ** 2 + np.abs(Ri) ** 2)))
        Lh, Rh = psi[0, 0:2, I0 - B // 2:I0, I0:I1, I0:I1], psi[0, 2:4, I0 - B // 2:I0, I0:I1, I0:I1]
        S["beta_halo"].append(float(np.sum(2 * np.real(np.conj(Lh) * Rh))))
    S["psi"] = psi
    cn = coarse_step(S["c_h"])
    up, dn = restrict(psi)
    cn[c0, c0, c0] = up
    S["down"].append(abs(dn) / max(abs(up), 1e-300))
    S["c_h"] = cn
    S["c_c"] = coarse_step(S["c_c"])
    S["C_h"].append(complex(np.vdot(S["c0vec"], S["c_h"])))
    S["C_c"].append(complex(np.vdot(S["c0vec"], S["c_c"])))
    S["n"] += 1
    print(f"step {S['n']}: |C_h| {abs(S['C_h'][-1]):.8f}, arg(C_h/C_c) {np.angle(S['C_h'][-1] / S['C_c'][-1]):+.3e}, "
          f"norm {np.linalg.norm(S['c_h']):.8f}, down/up {S['down'][-1]:.1e} ({time.time() - t0:.1f} s)", flush=True)
pickle.dump(S, open(CK, "wb"))
json.dump(res, open(os.path.join(OUT, "agg_zz_zoom.json"), "w"), indent=1, default=float)

if S["n"] >= NSTEP:
    n = np.arange(len(S["C_h"]))
    tau_phys = B * TH
    sel = n >= 8
    ph_h = np.unwrap(np.angle(np.array(S["C_h"])))
    ph_c = np.unwrap(np.angle(np.array(S["C_c"])))
    sh = np.polyfit(n[sel], ph_h[sel], 1)[0]
    sc = np.polyfit(n[sel], ph_c[sel], 1)[0]
    rest = at.rest
    Eh = -np.angle(np.exp(1j * (sh - rest))) / tau_phys
    Ec = -np.angle(np.exp(1j * (sc - rest))) / tau_phys
    # Zitterbewegung spectrum of the fine interior
    beta = np.array(S["beta"])
    nb = len(beta)
    skip = 4 * B
    bb = beta[skip:] - np.mean(beta[skip:])
    win = np.hanning(len(bb))
    sp = np.abs(np.fft.rfft(bb * win)) ** 2
    om = 2 * np.pi * np.fft.rfftfreq(len(bb), d=TH)        # m c^2 / hbar units
    band = om > 0.5
    ipk = np.argmax(np.where(band, sp, 0))
    frac = float(sp[ipk - 3:ipk + 4].sum() / sp[band].sum())
    res["hybrid"] = dict(E_eV=Eh * MC2, E_coarse_run_eV=Ec * MC2, rel_diff=Eh / Ec - 1, coarse_steps=int(S["n"]),
                         fine_ticks=int(S["n"] * B), fit_from_step=8, C_h_abs_end=abs(S["C_h"][-1]),
                         down_over_up_max=float(np.max(S["down"])), zb_peak_omega=float(om[ipk]), zb_peak_power_frac=frac,
                         zb_expected_2E=float(2 * (1 - 0.5 * zz.ALPHA ** 2)),
                         interior_prob_mean=float(np.mean(S["pint"][skip:])),
                         cost=dict(uniform_fine_cells=float((NC * B) ** 3), uniform_fine_cell_ticks=float((NC * B) ** 3 * S["n"] * B),
                                   hybrid_cell_ticks=float(NC ** 3 * S["n"] + NP ** 3 * S["n"] * B),
                                   ratio=float((NC * B) ** 3 * S["n"] * B / (NC ** 3 * S["n"] + NP ** 3 * S["n"] * B))))
    np.savez_compressed(os.path.join(OUT, "agg_zz_zoom.npz"), beta=beta, beta_halo=np.array(S["beta_halo"]),
                        pint=np.array(S["pint"]), C_h=np.array(S["C_h"]), C_c=np.array(S["C_c"]), om=om, sp=sp,
                        slice_fine=np.abs(S["psi"][0, :, :, :, NP // 2]) ** 2,
                        slice_coarse=np.abs(S["c_h"][:, :, c0]) ** 2)
    print(f"hybrid 1s: {Eh * MC2:.6f} eV vs uniform coarse {Ec * MC2:.6f} eV ({100 * (Eh / Ec - 1):+.4f} %); "
          f"ZB peak at {om[ipk]:.4f} m c^2/hbar (power fraction {frac:.2f}); cost ratio {res['hybrid']['cost']['ratio']:.2e}", flush=True)
    json.dump(res, open(os.path.join(OUT, "agg_zz_zoom.json"), "w"), indent=1, default=float)

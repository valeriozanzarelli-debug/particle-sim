"""multi_e_core.py - pezzi comuni dello studio MOLTI ELETTRONI (docs/MOLTI_ELETTRONI_PIANO.md).
Costanti, lettura dei livelli imparati della regola zig-zag agglomerata (lab/out/agg_zz_3d_{A,B}.pkl, gli stessi stencil
di assets/atomo/agglomerato_zz.json), banda simmetrizzata in Eh, direzioni sferiche, scrittura a checkpoint dei JSON."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("NUMBA_NUM_THREADS", "2")
import json, pickle, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import zz_core as zz                                   # noqa: E402

ALPHA = zz.ALPHA
A0 = 1 / ALPHA                                         # raggio di Bohr in lambda_C
EH_EV = 27.211386245988
THETA = {"A": A0 / 16 / 256, "B": A0 / 16 / 32, "A2": A0 / 16 / 512}
NUC_MASS = {1: 1836.15267343, 2: 7294.29954142, 3: 12786.393, 4: 16424.20, 10: 36433.98}


def load_json(name, default=None):
    p = os.path.join(OUT, name)
    if os.path.exists(p):
        return json.load(open(p))
    return {} if default is None else default


def save_json(name, obj):
    p = os.path.join(OUT, name)
    tmp = p + ".tmp"
    json.dump(obj, open(tmp, "w"), indent=1, default=float)
    os.replace(tmp, p)


def chain_levels(chain):
    """livelli imparati (dict con K = generatore G_j (2,2,ntap), offs, B, phi0, ...)"""
    st = pickle.load(open(os.path.join(OUT, f"agg_zz_3d_{chain}.pkl"), "rb"))
    return st["levels"]


class Band:
    """banda di moto (media di spin) di un livello pauli2 in Eh, k in a.u. (1/a0).
    eps(k) = [G_s(k) - G_s(0)] / tau / alpha^2, G_s = (G00 + G11)/2; sym=True -> stencil Re(G_j) (parte pari, DEBT-8)."""

    def __init__(self, K, offs, B, theta, sym=True):
        g = 0.5 * (K[0, 0] + K[1, 1])
        self.g = np.real(g) if sym else g
        self.offs = np.asarray(offs, float)
        self.B, self.theta = B, theta
        self.h_lc = B * theta                          # cella in lambda_C (= tick in hbar/mc^2)
        self.h = self.h_lc / A0                        # cella in a0
        self.tau_au = self.h_lc * ALPHA ** 2           # tick grosso in unita' atomiche
        self.sym = sym

    def gen(self, kap):
        """G_s(kappa), kappa in 1/cella (..., 3); reale se sym, altrimenti la sua parte reale (G Hermitiana -> reale)"""
        ph = np.exp(1j * np.tensordot(kap, self.offs.T, axes=(-1, 0)))
        return np.real(ph @ self.g)

    def eps(self, k_au):
        """energia di moto in Eh a impulsi k (a.u.) (..., 3)"""
        kap = np.asarray(k_au) * self.h
        z = np.zeros(3)
        return (self.gen(kap) - self.gen(z)) / self.h_lc * A0 ** 2   # (fase per passo)/tau [mc^2] -> Eh (x 1/alpha^2)

    def hj_eh(self):
        """stencil in Eh: H_j = g_j / tau / alpha^2"""
        return np.real(self.g) / self.h_lc * A0 ** 2


def fib_dirs(n=400):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def dirs26():
    out = []
    for v in np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij")).reshape(3, -1).T:
        if np.any(v):
            out.append(v / np.linalg.norm(v))
    return np.array(out)

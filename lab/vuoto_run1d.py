"""vuoto_run1d.py - production (pre-registered): 1+1 D exams V1a, V1b, V1c. Writes lab/out/vuoto_1d.json (~2 min)."""
import os, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from vuoto_core import chi_1d, chi1d_cont, chi1d_cont_direct, species_weight, ueig, PI, lr_terms


def walk1d_matrix(N, th, V):
    """one tick of zz_core 'dirac1' on a ring (same as vuoto_dev_check.py, checked there against k-space)"""
    S = np.zeros((2 * N, 2 * N))
    x = np.arange(N)
    S[(x + 1) % N, x] = 1
    S[N + (x - 1) % N, N + x] = 1
    c, s = np.cos(th), np.sin(th)
    C = np.block([[c * np.eye(N), -1j * s * np.eye(N)], [-1j * s * np.eye(N), c * np.eye(N)]])
    return np.diag(np.exp(-1j * np.concatenate([V, V]))) @ C @ S


def sea_density(W, ncomp, N):
    Om, lam, V = ueig(W[None])
    Vo = V[0][:, Om[0] < 0]
    return np.real(np.einsum("ij,ij->i", Vo, Vo.conj())).reshape(ncomp, -1).sum(0)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
t0 = time.time()
R = {}
# V1a theorem: whole one-tick sea, strong potential, N = 128
N = 128
x = np.arange(N)
Vs = -0.8 / np.sqrt((x - N / 2) ** 2 + 1.0)
for th in (0.1, 0.3):
    n = sea_density(walk1d_matrix(N, th, Vs), 2, N)
    R[f"V1a_spread_theta{th}"] = float(n.max() - n.min())
wc = [abs(chi_1d(qm * 0.05, 0.05, None, n=4096)) for qm in (0.3, 1, 3, 10)]
R["V1a_max_abs_chi_whole_kspace"] = float(max(wc))
R["V1a_PASS"] = bool(max(R["V1a_spread_theta0.1"], R["V1a_spread_theta0.3"], R["V1a_max_abs_chi_whole_kspace"]) < 1e-12)
# V1b physical species vs continuum
w = species_weight()
tab = {}
for qm in (0.3, 1, 3, 10):
    ref = chi1d_cont(qm)
    tab[qm] = {"cont": ref, "cont_direct_check": chi1d_cont_direct(qm) / ref - 1}
    for th in (0.2, 0.1, 0.05, 0.025):
        c = chi_1d(qm * th, th, w, n=3000, beta=4.0)
        tab[qm][th] = float(c.real / ref - 1)
R["V1b_table_relerr"] = {str(k): {str(a): b for a, b in v.items()} for k, v in tab.items()}
okA = all(abs(tab[qm][0.025]) <= 0.01 for qm in tab)
okB = all(abs(tab[qm][0.05]) >= 3 * abs(tab[qm][0.025]) for qm in (0.3, 1, 3))
R["V1b_PASS"] = bool(okA and okB)
# V1c period independence
d = []
for qm in (0.3, 3):
    for th in (0.1, 0.025):
        c1 = chi_1d(qm * th, th, w, n=3000, beta=4.0, period=1)
        c2 = chi_1d(qm * th, th, w, n=3000, beta=4.0, period=2)
        d.append(abs(c1 - c2) / abs(c1))
R["V1c_max_rel_diff"] = float(max(d))
R["V1c_PASS"] = bool(max(d) < 1e-9)
R["time_s"] = time.time() - t0
json.dump(R, open(os.path.join(OUT, "vuoto_1d.json"), "w"), indent=1)
print(json.dumps(R, indent=1))

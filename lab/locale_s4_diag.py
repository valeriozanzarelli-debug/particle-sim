"""locale_s4_diag.py - POST-HOC diagnostics of phase 2 (NOT pre-registered, they do not change any verdict):
why S4 failed P3 (order for s states) and P4 (Schwinger at theta = 0.1). Writes lab/out/locale_s4_diag.json.
  D1 hydrogen S4 vs frozen on a finer theta scan (s and p states).
  D2 the 'staggered twin': P = (-1)^i commutes with every coin and anticommutes with the shift, so P W P = (-1)^(sum|m|) W:
     for S4 (10 shifts per period) every radial state has an exactly degenerate staggered twin; check it numerically.
  D3 Schwinger at smaller theta with the same sweep in units of m (k from -4m to +4m): does S4 converge?"""
import os, sys, json
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_s4_core import SEQ_S4, SEQ_FROZEN, dirac1_bloch, RadialS4
from zz_core import dirac_formula, ALPHA

PI = np.pi
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_s4_diag.json")
D = {}
SEQS = {"S4": SEQ_S4, "frozen": SEQ_FROZEN}

# D1
D["D1_hydrogen"] = {}
for sname, seq in SEQS.items():
    for th in (0.2, 0.15, 0.1, 0.07, 0.06, 0.05, 0.04, 0.03, 0.025):
        row = {}
        for lab, n, kap in (("1s", 1, -1), ("2s", 2, -1), ("2p1/2", 2, 1), ("2p3/2", 2, -2)):
            E = dirac_formula(n, kap)
            om, x, res = RadialS4(kap, th, 30 / ALPHA, seq).eigen(E, iters=12)
            row[lab] = float(((1 - om) - (1 - E)) / (1 - E))
        D["D1_hydrogen"][f"{sname}|{th}"] = row
        print(sname, th, {k: f"{v:.2e}" for k, v in row.items()}, flush=True)

# D2
rs = RadialS4(-1, 0.1, 30 / ALPHA, SEQ_S4)
W = rs.matrix()
Pd = sp.diags(np.repeat((-1.0) ** np.arange(rs.N), 2).astype(complex))
D["D2_twin"] = dict(norm_PWP_minus_W=float(abs(Pd @ W @ Pd - W).max()),
                    norm_PWP_plus_W_frozen=float(abs(Pd @ RadialS4(-1, 0.1, 30 / ALPHA, SEQ_FROZEN).matrix() @ Pd
                                                     + RadialS4(-1, 0.1, 30 / ALPHA, SEQ_FROZEN).matrix()).max()))
print(D["D2_twin"], flush=True)


# D3
def lz(seq, th, EEc, kmax_m=4.0):
    E = EEc * th ** 2
    T = float(sum(seq))
    kmax = kmax_m * th
    nper = int(round(2 * kmax / (E * T)))
    def bands(Wm):
        w, v = np.linalg.eig(Wm); i = np.argsort(-np.angle(w)); return v[:, i]
    psi = bands(dirac1_bloch(-kmax, th, seq)[0])[:, 0].copy()
    for p in range(nper):
        psi = dirac1_bloch(-kmax, th, seq, A=lambda t: -E * t, t0=p * T)[0] @ psi
    kend = -kmax + E * nper * T
    v1 = bands(dirac1_bloch(kend, th, seq)[0])
    # continuum, same sweep
    n = int(np.ceil((kend + kmax) / E / (0.01 * th / 0.1)))
    dt = (kend + kmax) / E / n
    sz = np.diag([1.0, -1.0]); sx = np.array([[0, 1.0], [1.0, 0]])
    H = lambda k: k * sz + th * sx
    psic = np.linalg.eigh(H(-kmax))[1][:, 0].astype(complex)
    for j in range(n):
        kk = -kmax + E * (j + 0.5) * dt; e = np.hypot(kk, th)
        psic = (np.cos(e * dt) * np.eye(2) - 1j * np.sin(e * dt) * H(kk) / e) @ psic
    Pc = abs(np.vdot(np.linalg.eigh(H(kend))[1][:, 1], psic)) ** 2
    return float(abs(np.vdot(v1[:, 1], psi)) ** 2), float(Pc), nper


D["D3_schwinger"] = {}
for sname, seq in SEQS.items():
    for th in (0.1, 0.05, 0.025, 0.0125):
        for EEc in (1, 2, 4):
            P, Pc, nper = lz(seq, th, EEc)
            Ps = float(np.exp(-PI / EEc))
            D["D3_schwinger"][f"{sname}|{th}|{EEc}"] = dict(P=P, P_cont_same_sweep=Pc, P_schwinger=Ps,
                                                            vs_schwinger=P / Ps - 1, vs_continuum=P / Pc - 1,
                                                            periods=nper, kick_per_period_over_gap=EEc * th ** 2 *
                                                            sum(seq) / (2 * th))
            print(sname, th, EEc, f"P/Ps-1 {P / Ps - 1:+.3f}  P/Pc-1 {P / Pc - 1:+.3f}", flush=True)

with open(OUT, "w") as f:
    json.dump(D, f, indent=1)
print("done")

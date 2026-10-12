"""agg_zz_1d.py - A1: the agglomerator on the 1-D zig-zag walk (chain A, theta = a0/4096).
(1) the fixed b = 2 ladder B = 2..32 (w chosen by the automatic rule), (2) the automatic hierarchy (b in {4,3,2}) up to
B <= 512. Each level is learned from the level below (true repetition). Exams A1a-c (lab/out/agg_zz_prereg.json).
Writes lab/out/agg_zz_1d.json and lab/out/agg_zz_1d.npz (kernels, dispersions).
Run: python lab/agg_zz_1d.py"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import json, time
import numpy as np
from agg_zz_core import *

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
TH = A0 / 16 / 256
rng = np.random.default_rng(20261012)
t_start = time.time()
NC = 256          # coarse cells of the learning box
EPS_LEVEL = 1e-3
res = dict(theta=TH, ladder=[], auto=[])


def fine_exact_eps(k):
    """positive branch of the fine 1-D walk per fine time unit, minus the rest energy: cos w = cos th cos k"""
    return np.arccos(np.cos(TH) * np.cos(k)) - TH


def exam_true_fine(lv, S=4, Nt=64):
    """A1a: coarse rule vs the TRUE fine walk over ~1000 fine ticks"""
    B = lv.B
    n = Nt * B
    P = gauss_packets(rng, S, 2, n, 1, 4 * B, 0.5 / B)
    nst = max(1, 1000 // B)
    F = P.copy()
    for _ in range(nst * B):
        F = fine1_tick(F, TH)
    e = rel_err(lv.rule.step(block(P, B, 1), nst), block(F, B, 1))
    return float(e.max()), nst, nst * B


def dispersion(lv):
    """positive branch of the learned kernel (eigenvector with the larger large-component weight), unwrapped from k = 0"""
    B = lv.B
    kap = np.linspace(0, np.pi, 4001)
    M = lv.rule.mult_at(kap[:, None])                       # (n, 2, 2) in (R, L)
    lam, V = np.linalg.eig(M)
    wphi = np.abs(V[:, 0, :] + V[:, 1, :]) ** 2 / 2
    ip = np.argmax(wphi, axis=1)
    lp = lam[np.arange(len(kap)), ip]
    ph = np.unwrap(np.angle(lp))
    eps = -(ph - ph[0]) / B                                  # per fine time unit
    k = kap / B
    rest_phase = float(np.angle(lp[0]))
    mod = np.abs(lam).max(axis=1)
    return k, eps, rest_phase, mod, lam


def fitc(k, eps, kmax):
    s = (k > 0) & (k <= kmax)
    A = np.stack([k[s] ** 2, k[s] ** 4, k[s] ** 6], 1)
    c = np.linalg.lstsq(A, eps[s], rcond=None)[0]
    return c


def analyse(lv, tag):
    e, nst, nt = exam_true_fine(lv)
    k, eps, rp, mod, lam = dispersion(lv)
    B, m = lv.B, TH
    kmax = min(0.3 * m, 0.8 / B)
    c = fitc(k, eps, kmax)
    cf = fitc(k, fine_exact_eps(k), kmax)
    mkin = 1 / (2 * c[0])
    c4n = -c[1] * 8 * m ** 3
    c4f = -cf[1] * 8 * m ** 3
    sel = k <= 0.5 * np.pi / B
    dev = float(np.max(np.abs(eps[sel] - fine_exact_eps(k[sel]))) * B)          # phase error per coarse step, kappa<pi/2
    taps = len(lv.rule.offs)
    rec = dict(tag=tag, B=B, b=lv.b, w=lv.w, check=lv.check, taps=int(len(lv.rule.offs)), exam_err=e, exam_coarse_steps=nst, exam_fine_ticks=nt,
               cost_reduction_cells_ticks=B * B, cost_reduction_with_stencil=B * B * 4 / (taps * 4),
               rest_phase=rp, rest_phase_expected=float(np.angle(np.exp(-1j * TH * B))),
               m_kin_over_m=mkin / m, m_kin_fine_over_m=(1 / (2 * cf[0])) / m, c4_over_1_8m3=c4n, c4_fine_over_1_8m3=c4f,
               fit_kmax_over_m=kmax / m, fit_kappa_max=kmax * B, max_phase_err_per_step_band=dev,
               max_abs_eig=float(mod.max()), seconds=round(time.time() - t_start, 1))
    print(f"  [{tag}] B={B:4d} b={lv.b} w={lv.w}: exam |dpsi|^2 = {e:.2e} over {nt} fine ticks ({nst} steps); "
          f"m_kin/m = {mkin / m:.6f} (fine {rec['m_kin_fine_over_m']:.6f}); c4*8m^3 = {c4n:.4f} (fine {c4f:.4f}); "
          f"rest phase {rp:+.5f} vs {rec['rest_phase_expected']:+.5f}; max|K| {mod.max():.6f}", flush=True)
    return rec, dict(k=k[::20], eps=eps[::20], eps_fine=fine_exact_eps(k[::20]), K=lv.rule.K, offs=lv.rule.offs)


fine = Level("fine1", 1, "RL", 2, 1, TH)
arrays = {}
for chain, bs, Bmax, key in (("ladder", (2,), 32, "ladder"), ("auto", (4, 3, 2), 512, "auto")):
    print(f"=== chain {chain} ===", flush=True)
    prev = fine
    while True:
        lv, tried = auto_level(prev, rng, NC, eps_packet=1e-3, eps_energy=1e-4, T_obs=1000, k_obs=0.3 * TH, bs=bs,
                               ws=(1, 2, 3, 4, 5, 6), bmax_total=Bmax, ks_rule=lambda B: float(np.clip(0.3 * TH * B, 0.05, 0.6)))
        if lv is None:
            print(f"  stop at B = {prev.B} (no admissible candidate below B <= {Bmax})", flush=True)
            break
        rec, arr = analyse(lv, chain)
        rec["candidates"] = tried
        res[key].append(rec)
        for kk, vv in arr.items():
            arrays[f"{chain}_B{lv.B}_{kk}"] = vv
        prev = lv
    json.dump(res, open(os.path.join(OUT, "agg_zz_1d.json"), "w"), indent=1, default=float)
np.savez_compressed(os.path.join(OUT, "agg_zz_1d.npz"), **arrays)

# verdicts
ver = []
for key in ("ladder", "auto"):
    for r in res[key]:
        ver.append((f"A1a {key} B={r['B']}", r["exam_err"] < 0.01))
        if r["B"] >= 32:
            ver.append((f"A1b {key} B={r['B']} m_kin", abs(r["m_kin_over_m"] - 1) < 0.01))
            ver.append((f"A1c {key} B={r['B']} c4", abs(r["c4_over_1_8m3"] - 1) < 0.10))
res["verdict"] = ver
for v in ver:
    print(f"[{'PASS' if v[1] else 'FAIL'}] {v[0]}")
res["seconds"] = time.time() - t_start
json.dump(res, open(os.path.join(OUT, "agg_zz_1d.json"), "w"), indent=1, default=float)
print(f"done {time.time() - t_start:.0f} s")

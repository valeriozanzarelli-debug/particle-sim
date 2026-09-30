"""TASK 5 — polarization with identical balls (each ball carries a Jones vector inside its ψ).

M2 (polarizing beam splitter, no coin): a free Gaussian beam (w0 = 3), polarization H, hits a PBS at angle θ
at z = 0: every ball's ψ becomes  a⟨a|H⟩ e^{+ikαx} + b⟨b|H⟩ e^{−ikαx}  (a = (cosθ, sinθ), b ⟂ a, α = 0.25).
The balls then follow rule C (ψ carried, current summed over the two polarization components). The two
branches separate; the fraction of balls that end up in branch a (x > 0 at z = 40) is compared with cos²θ.
Nothing random happens at the PBS: which branch a ball joins is decided by where it sits in the beam.

M1 (absorbing polarizer): each ball projects its ψ on a; to keep 'ball density = |ψ|²' a ball must vanish
with probability 1 − |⟨a|ψ_i⟩|²/|ψ_i|² (= sin²θ for every ball here) — a local coin per ball.

  python lab/identical_polar.py [N=1000]
"""
import sys, time, numpy as np
import identical_core as ic
from identical_core import K, run_psi, free_start

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
ALPHA, ZEND = 0.25, 25.0
H = np.array([1.0, 0.0], complex)
print(f"M2: PBS, N={N} balls (quantile start), α={ALPHA}, balls counted at z={ZEND}", flush=True)
for deg in (0.0, 22.5, 45.0, 67.5, 90.0):
    th = np.deg2rad(deg); t = time.time()
    a = np.array([np.cos(th), np.sin(th)], complex); b = np.array([-np.sin(th), np.cos(th)], complex)
    x, g = free_start(N)
    psi = (g[:, None] * (np.vdot(a, H) * a[None, :] * np.exp(1j * K * np.sin(ALPHA) * x)[:, None]
                         + np.vdot(b, H) * b[None, :] * np.exp(-1j * K * np.sin(ALPHA) * x)[:, None]))
    sn, pp, info = run_psi(x, psi, (ZEND,), 3, 2, lambda x, z: 0.05, scheme="cn2", tmax=500)
    X = sn[0]
    if X is None:
        print(f"  θ={deg:5.1f}°  FAILED at z={info['z']:.1f}", flush=True); continue
    Pa = (X > 0).mean()
    # polarization carried by the balls of branch a: |<a|ψ>|²/|ψ|²
    P = pp[0]; sel = X > 0
    pol = (np.abs(P[sel] @ a.conj()) ** 2 / (np.abs(P[sel]) ** 2).sum(-1)).mean() if sel.any() else np.nan
    print(f"  θ={deg:5.1f}°  fraction of balls in branch a {Pa:.4f}   cos²θ {np.cos(th)**2:.4f}   "
          f"(granularity 1/N = {1/N:.4f})   mean |<a|ψ>|²/|ψ|² in branch a {pol:.4f}   swaps {info['swaps']} [{time.time()-t:.0f}s]",
          flush=True)

print(f"M1: absorbing polarizer, local coin per ball, N={N}", flush=True)
rng = np.random.default_rng(3)
for deg in (0.0, 22.5, 45.0, 67.5, 90.0):
    th = np.deg2rad(deg); a = np.array([np.cos(th), np.sin(th)], complex)
    x, g = free_start(N); psi = g[:, None] * H[None, :]
    p = np.abs(psi @ a.conj()) ** 2 / (np.abs(psi) ** 2).sum(-1)
    alive = rng.random(N) < p
    print(f"  θ={deg:5.1f}°  surviving balls {alive.mean():.4f}   cos²θ {np.cos(th)**2:.4f}", flush=True)

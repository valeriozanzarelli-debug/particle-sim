"""kk_ab.py - test 3b: Aharonov-Bohm ring (1+1 D ring of L sites x hidden circle of N sites).
Two packets with +p and -p leave x = 0, go round the two halves of the ring and meet at x = L/2.
The hidden circle is twisted by a total Phi (in radians of the hidden circle) around the ring: on the
ring itself there is NO curvature (no 'field'), only the hole has it. Measured: phase of the
interference term at the meeting point vs n * Phi; same physics if the twist sits on one link or is
spread over all links (gauge invariance; the initial state is gauge-transformed accordingly).
Writes lab/out/kk_ab.json/.npz"""
import json, time
import numpy as np
import kk_core as K
from zz_core import bloch_batch

t0 = time.time()
out = {}
# ------------------------- (b) Aharonov-Bohm ring -------------------------
L = 400; N = 32; th = 0.2; p = 0.5; sig = 10.0
R = N / (2 * np.pi)
def vg(n):
    f = lambda k: np.max(-np.angle(np.linalg.eigvals(bloch_batch(np.array([[k, 0, 2 * np.pi * n / N]]), th)[0])) / 2)
    return (f(p + 1e-4) - f(p - 1e-4)) / 2e-4
ab = []; prof = {}
for n in (+1, -1, 0, 2):
    ticks = int(round((L / 2) / abs(vg(n)) / 2)) * 2
    ref = None
    for Phi in (0, np.pi / 4, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi):
        Dx = np.zeros((L, 1)); Dx[L // 4] = R * Phi          # all the twist on ONE link (x = L/4)
        rho, po, pf = K.ring_ab(L, N, th, n, Dx, p=p, sigma=sig, ticks=ticks)
        if ref is None:
            ref = (po, pf)
        d_o = float(np.angle(np.exp(1j * (po - ref[0]))))
        d_f = float(np.angle(np.exp(1j * (pf - ref[1]))))
        pred = float(np.angle(np.exp(1j * n * Phi)))
        gauge = None
        if Phi == np.pi / 2:                                  # same total twist spread on all L links
            Ds = np.full((L, 1), R * Phi / L)
            xx = np.arange(L)
            lam = n * Phi * (xx / L - (xx >= L // 4))          # gauge function: one-link -> spread
            rho_s, _, _ = K.ring_ab(L, N, th, n, Ds, p=p, sigma=sig, ticks=ticks, init_phase=np.exp(-1j * lam))
            gauge = float(np.max(np.abs(rho_s - rho)) / rho.max())
        ab.append(dict(n=n, Phi=float(Phi), phase_shift_interference_term=d_o, phase_shift_fringe_fit=d_f,
                       n_Phi_mod_2pi=pred, ticks=ticks, rel_diff_one_link_vs_spread=gauge))
        prof[(n, round(Phi, 3))] = rho
        print(f"(b) n={n:+d} Phi={Phi:5.3f}: AB phase (interference term) {d_o:+.5f}, fringe fit {d_f:+.4f}, "
              f"n*Phi mod 2pi = {pred:+.5f}" + (f"; one-link vs spread twist: max|drho|/max rho = {gauge:.1e}" if gauge is not None else ""))
out["aharonov_bohm"] = ab
out["ab_params"] = dict(L=L, N=N, theta=th, p=p, sigma=sig, twist_link=L // 4, window=60, note="packets start at x=0 and meet at x=L/2")
out["time_s"] = time.time() - t0
json.dump(out, open("out/kk_ab.json", "w"), indent=1)
np.savez_compressed("out/kk_ab_tmp.npz", r1_0=prof[(1, 0.0)], r1_h=prof[(1, round(np.pi / 2, 3))],
                    rm_h=prof[(-1, round(np.pi / 2, 3))], r0_h=prof[(0, round(np.pi / 2, 3))])
print("time %.1f s" % out["time_s"])

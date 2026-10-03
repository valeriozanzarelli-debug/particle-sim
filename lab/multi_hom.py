"""Hong–Ou–Mandel with the multi-photon pilot rule.

Two photons enter a 50:50 beam splitter from ports a and b of the SAME arm (they share its wave balls).
Identical bosons: Ψ = [ψ_a(x1)ψ_b(x2) + ψ_b(x1)ψ_a(x2)]/√2 (with polarization and spectral tags).
Distinguishability: a delay τ is a phase e^{iωτ} on photon 2's spectral tag (Gaussian spectrum, σ_ω = 1);
or a rotation of photon 2's polarization by θ.  A coincidence = the two photon balls land in
different output ports (c and d).  Theory: P_cd = ½(1 − |<τ1|τ2>|² |<p1|p2>|²).
"""
import time
import numpy as np
from multi_core import System, bs_arm, lin, H

T0 = time.time()
W = np.linspace(-5, 5, 41)
F = np.exp(-W ** 2 / 4); F /= np.linalg.norm(F)


def run(tau, theta, N, seed):
    arm = bs_arm("BS")
    t1, t2 = F, F * np.exp(1j * W * tau)
    p1, p2 = H, lin(np.deg2rad(theta))
    s = 1 / np.sqrt(2)
    sysm = System([arm], ["BS", "BS"], [(s, [("a", p1, t1), ("b", p2, t2)]),
                                        (s, [("b", p2, t2), ("a", p1, t1)])])
    X0 = sysm.sample(N, np.random.default_rng(seed))
    pa, pb = sysm.ports(sysm.run(X0))
    pcd = (pa != pb).mean()
    th = 0.5 * (1 - abs(np.vdot(t1, t2)) ** 2 * abs(np.vdot(p1, p2)) ** 2)
    return pcd, th, ((pa == "c") & (pb == "c")).mean(), ((pa == "d") & (pb == "d")).mean()


def main():
    N = 3000
    res = []
    print(f"Hong–Ou–Mandel, {N} pairs per point, BS arm (w0 = 5, deflection ±0.25 rad)", flush=True)
    print(" delay τ·σω   P(coinc.)  theory   both-c  both-d   verdict", flush=True)
    pts = {}
    for i, tau in enumerate((0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0)):
        t = time.time()
        p, th, cc, dd = run(tau, 0.0, N, 10 + i)
        sd = np.sqrt(max(th * (1 - th), 1 / N) / N)
        ok = abs(p - th) < max(3 * sd, 3 / N)
        res.append(ok); pts[tau] = p
        print(f"   {tau:4.1f}       {p:.4f}    {th:.4f}   {cc:.3f}   {dd:.3f}   {'PASS' if ok else 'FAIL'}  [{time.time()-t:.0f}s]", flush=True)
    V = 1 - pts[0.0] / pts[5.0]
    ok = V >= 0.95
    res.append(ok)
    print(f"dip visibility V = 1 − P(0)/P(5/σω) = {V:.4f}  (required ≥ 0.95; theory {1 - 0/0.5:.3f})  {'PASS' if ok else 'FAIL'}", flush=True)
    print(" polarization angle θ between the photons (τ = 0):", flush=True)
    for i, th_deg in enumerate((0, 30, 45, 60, 90)):
        p, th, cc, dd = run(0.0, th_deg, N, 50 + i)
        sd = np.sqrt(max(th * (1 - th), 1 / N) / N)
        ok = abs(p - th) < max(3 * sd, 3 / N)
        res.append(ok)
        print(f"   θ = {th_deg:2d}°   P(coinc.) {p:.4f}   theory {th:.4f}   {'PASS' if ok else 'FAIL'}", flush=True)
    print(f"RESULT: {sum(res)}/{len(res)} PASS   total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

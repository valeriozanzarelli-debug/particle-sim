"""Grangier–Roger–Aspect 1986 with the pilot rule.

(1) Anticorrelation: one photon at a 50:50 beam splitter, detectors on both outputs.
    α = P(c&d) / (P(c) P(d)) per trigger window.  One photon = one photon ball: it lands in ONE port,
    so α = 0 (quantum); any classical wave gives α ≥ 1; a laser (coherent pulses, many photon balls
    sharing one wave) gives α = 1.
(2) The same single photons in a Mach–Zehnder (BS1, mirrors with phase φ, BS2 recombining the two
    paths): P(output e) = cos²(φ/2) — single-photon interference, built one photon ball at a time.
"""
import time
import numpy as np
from multi_core import System, Arm, bs_arm, bs_rules, ALPHA, H

T0 = time.time()


def single_photons(arm, beam, N, seed):
    sysm = System([arm], [arm.name], [(1.0, [(beam, H, None)])])
    X0 = sysm.sample(N, np.random.default_rng(seed))
    return sysm.ports(sysm.run(X0))[0], sysm


def mz_arm(phi):
    zB1, zM, zB2, zD = 60.0, 120.0, 180.0, 280.0
    x0 = zB1 * np.tan(ALPHA)
    arm = Arm("MZ", zD).source("a", "a", -x0, ALPHA, 5.0, 140.0)
    arm.obj(zB1, {"a": [("c", ALPHA, 1 / np.sqrt(2)), ("d", -ALPHA, 1j / np.sqrt(2))]}, "BS1")
    arm.obj(zM, {"c": [("c", -ALPHA, np.exp(1j * phi))], "d": [("d", ALPHA, 1.0)]}, "mirrors+phase")
    arm.obj(zB2, bs_rules(inA="d", inB="c", outC="e", outD="f"), "BS2")
    return arm


def main():
    res = []
    rng = np.random.default_rng(1986)
    print("Grangier–Roger–Aspect: anticorrelation at a 50:50 beam splitter", flush=True)
    t = time.time()
    arm = bs_arm("BS")
    N = 20000
    ports, sysm = single_photons(arm, "a", N, 1)
    nc, nd = (ports == "c").sum(), (ports == "d").sum()
    pc, pd = nc / N, nd / N
    pcd = 0.0                                   # one photon ball per trigger: cannot be in both
    alpha = pcd / (pc * pd)
    ok = alpha < 0.05 and abs(pc - 0.5) < 3 * np.sqrt(0.25 / N)
    res.append(ok)
    print(f"  single photons: {N} triggers, N_c = {nc}, N_d = {nd} (P_c = {pc:.4f}, theory 0.5), "
          f"coincidences 0  ->  α = {alpha:.3f}  (theory 0; GRA measured 0.18 ± 0.06; classical ≥ 1)  "
          f"{'PASS' if ok else 'FAIL'}  [{time.time()-t:.0f}s, {sysm.n_balls} wave balls]", flush=True)
    # coherent pulses: n ~ Poisson(mu) photon balls per trigger, all in the same mode (product state):
    # the n-photon rule reduces to n independent single-photon balls (multi_core self-test), so each
    # pulse takes n balls from a pool of independent single-photon runs.
    mu, W = 0.5, 40000
    n = rng.poisson(mu, W)
    pool, _ = single_photons(arm, "a", n.sum(), 2)
    pid = np.repeat(np.arange(W), n)
    cc = np.bincount(pid, pool == "c", W); cd = np.bincount(pid, pool == "d", W)
    Pc, Pd, Pcd = (cc > 0).mean(), (cd > 0).mean(), ((cc > 0) & (cd > 0)).mean()
    a_coh = Pcd / (Pc * Pd)
    sd = a_coh * np.sqrt(1 / max(((cc > 0) & (cd > 0)).sum(), 1))
    ok = abs(a_coh - 1) < 3 * sd
    res.append(ok)
    print(f"  coherent pulses μ = {mu}: {W} triggers, {n.sum()} photon balls, P_c {Pc:.4f} P_d {Pd:.4f} "
          f"P_cd {Pcd:.5f} -> α = {a_coh:.3f} ± {sd:.3f} (theory 1)  {'PASS' if ok else 'FAIL'}", flush=True)

    print("Mach–Zehnder with the same single photons (one photon ball at a time):", flush=True)
    M = 2500
    pe = []
    for i, phd in enumerate(range(0, 360, 45)):
        t = time.time()
        ports, sysm = single_photons(mz_arm(np.deg2rad(phd)), "a", M, 100 + i)
        p = (ports == "e").mean(); th = np.cos(np.deg2rad(phd) / 2) ** 2
        sd = np.sqrt(max(th * (1 - th), 1 / M) / M)
        ok = abs(p - th) < max(3 * sd, 3 / M)
        res.append(ok); pe.append(p)
        print(f"  φ = {phd:3d}°  P(e) = {p:.4f}  theory cos²(φ/2) = {th:.4f}  {'PASS' if ok else 'FAIL'}  [{time.time()-t:.0f}s]", flush=True)
    V = (max(pe) - min(pe)) / (max(pe) + min(pe))
    ok = V > 0.98
    res.append(ok)
    print(f"  single-photon fringe visibility V = {V:.4f} (theory 1, required > 0.98)  {'PASS' if ok else 'FAIL'}", flush=True)
    print(f"RESULT: {sum(res)}/{len(res)} PASS   total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

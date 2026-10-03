"""Hanbury Brown–Twiss with the pilot rule: g2 = <n_c n_d> / (<n_c><n_d>) behind a 50:50 beam splitter.

Same arm, same rule for the three kinds of light; only the SOURCE differs (how many photon balls it puts
into the one shared wave in each time window):
  * single photon : exactly 1 ball per window                         -> g2(0) = 0
  * laser         : n ~ Poisson(μ)                                     -> g2(0) = 1
  * thermal       : M emitters with random phases and random frequencies add up to the mode amplitude
                    A(t) = Σ_m e^{i(φ_m + ω_m t)}/√M;  n ~ Poisson(μ|A(t)|²)   -> g2(0) = 2 − 1/M,
                    g2(τ) = 1 + |g1(τ)|² (1 − 1/M) with g1(τ) = e^{−τ²/2}  (σ_ω = 1)
All photons of a window are in the same mode (product state): the n-photon rule reduces to n
independent balls sharing one set of wave balls (multi_core self-test), so the balls of all windows
are run as one batch.  The bunching therefore comes from the wave intensity fluctuations, not from the
beam splitter.
"""
import time
import numpy as np
from multi_core import System, bs_arm, H

T0 = time.time()


def ball_pool(arm, N, seed):
    sysm = System([arm], [arm.name], [(1.0, [("a", H, None)])])
    X0 = sysm.sample(N, np.random.default_rng(seed))
    return sysm.ports(sysm.run(X0))[0]


def counts(n, pool):
    W = n.size
    pid = np.repeat(np.arange(W), n)
    p = pool[:n.sum()]
    return np.bincount(pid, p == "c", W), np.bincount(pid, p == "d", W)


def g2(nc, nd, lag=0):
    a, b = (nc, nd) if lag == 0 else (nc[:-lag], nd[lag:])
    return np.mean(a * b) / (np.mean(a) * np.mean(b))


def g2_err(nc, nd, lag=0, blocks=40):
    W = nc.size // blocks
    v = [g2(nc[i * W:(i + 1) * W], nd[i * W:(i + 1) * W], lag) for i in range(blocks)]
    return np.std(v) / np.sqrt(blocks)


def main():
    rng = np.random.default_rng(1956)
    arm = bs_arm("BS")
    W, mu, M, dt = 60000, 0.5, 50, 0.1
    res = []
    print(f"Hanbury Brown–Twiss behind a 50:50 BS; {W} windows per source", flush=True)
    t = time.time()
    # thermal source
    phi = rng.uniform(0, 2 * np.pi, M); om = rng.normal(0, 1, M)
    tt = np.arange(W) * dt
    A = np.zeros(W, complex)
    for m in range(M):
        A += np.exp(1j * (phi[m] + om[m] * tt))
    A /= np.sqrt(M)
    I = np.abs(A) ** 2
    n_th = rng.poisson(mu * I / I.mean())
    n_coh = rng.poisson(mu, W)
    n_one = np.ones(W // 4, int)
    need = n_th.sum() + n_coh.sum() + n_one.sum()
    pool = ball_pool(arm, need, 7)
    print(f"  {need} photon balls guided through the BS in one batch, sharing one set of wave balls "
          f"({arm.n_balls}) [{time.time()-t:.0f}s]", flush=True)
    i0 = 0
    out = {}
    for name, n in (("single photon", n_one), ("laser (coherent)", n_coh), ("thermal", n_th)):
        nc, nd = counts(n, pool[i0:i0 + n.sum()]); i0 += n.sum()
        out[name] = (nc, nd)
    # theory for the finite-M thermal field actually generated: g2 = <I^2>/<I>^2
    th_thermal = np.mean(I ** 2) / np.mean(I) ** 2
    for name, th in (("single photon", 0.0), ("laser (coherent)", 1.0), ("thermal", th_thermal)):
        nc, nd = out[name]
        g = g2(nc, nd); e = g2_err(nc, nd)
        ok = abs(g - th) < max(3 * e, 1e-9)
        res.append(ok)
        print(f"  {name:17s}: windows {nc.size}, clicks c {nc.sum()} d {nd.sum()}, g2(0) = {g:.3f} ± {e:.3f}  "
              f"theory {th:.3f}{' (2 − 1/M = %.3f for M = %d emitters)' % (2 - 1 / M, M) if name == 'thermal' else ''}  "
              f"{'PASS' if ok else 'FAIL'}", flush=True)
    nc, nd = out["thermal"]
    print("  thermal g2(τ) (bunching disappears after the coherence time 1/σω):", flush=True)
    for tau in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):
        lag = int(round(tau / dt))
        g = g2(nc, nd, lag); e = g2_err(nc, nd, lag)
        Il = I[:-lag] if lag else I
        Ir = I[lag:] if lag else I
        th = np.mean(Il * Ir) / (np.mean(Il) * np.mean(Ir))
        ideal = 1 + np.exp(-tau ** 2) * (1 - 1 / M)
        ok = abs(g - th) < 3 * e
        res.append(ok)
        print(f"    τ = {tau:3.1f}  g2 = {g:.3f} ± {e:.3f}   theory (this field) {th:.3f}   ideal 1+e^(−τ²) = {1+np.exp(-tau**2):.3f}  "
              f"{'PASS' if ok else 'FAIL'}", flush=True)
    print(f"RESULT: {sum(res)}/{len(res)} PASS   total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

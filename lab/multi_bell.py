"""Bell / CHSH with polarization-entangled pairs |Φ+> = (HH + VV)/√2 — the multi-photon pilot rule.

Each photon flies in its own arm to a polarizing beam splitter (PBS) set at angle a (b); the outcome is
the output port where its photon ball lands (+ = transmitted e1, − = reflected e2).  No collapse step:
the balls follow the joint wave; the velocity of A depends on where B is (the one nonlocal ingredient).

Also run: (1) the Born-rule collapse model of src/photon/entangle.rs (re-implemented here for the same
numbers of pairs), (2) the same pilot rule with A's PBS before B's, and (3) the dependence of a single
pair's outcome on the distant setting (same initial ball positions, other setting).
"""
import time, sys
import numpy as np
from multi_core import System, pbs_arm, lin, H, V

DEG = np.pi / 180
A_SET, B_SET = (0.0, 45.0), (22.5, 67.5)
T0 = time.time()


def pilot(a, b, N, seed, zA=60.0, zB=60.0, X0=None):
    armA = pbs_arm("A", lin(a * DEG), zP=zA, zD=zA + 100, zw=zA)
    armB = pbs_arm("B", lin(b * DEG), zP=zB, zD=zB + 100, zw=zB)
    s = 1 / np.sqrt(2)
    sysm = System([armA, armB], ["A", "B"], [(s, [("s", H, None), ("s", H, None)]),
                                              (s, [("s", V, None), ("s", V, None)])])
    if X0 is None:
        X0 = sysm.sample(N, np.random.default_rng(seed))
    Xf = sysm.run(X0)
    pa, pb = sysm.ports(Xf)
    return np.where(pa == "+", 1, -1), np.where(pb == "+", 1, -1), X0, sysm


def collapse(a, b, N, rng):
    """entangle.rs: A measured first with Born probability 1/2, B collapses to A's result."""
    plusA = rng.random(N) < 0.5
    ang = np.where(plusA, a * DEG, a * DEG + np.pi / 2)
    plusB = rng.random(N) < np.cos(ang - b * DEG) ** 2
    return np.where(plusA, 1, -1), np.where(plusB, 1, -1)


def chsh(fn, N):
    E, sig, marg = {}, {}, {}
    for a in A_SET:
        for b in B_SET:
            ra, rb = fn(a, b, N)
            E[a, b] = np.mean(ra * rb); sig[a, b] = np.sqrt((1 - E[a, b] ** 2) / N)
            marg[a, b] = ((ra > 0).mean(), (rb > 0).mean())
    a, a2 = A_SET; b, b2 = B_SET
    S = E[a, b] - E[a, b2] + E[a2, b] + E[a2, b2]
    sS = np.sqrt(sum(v ** 2 for v in sig.values()))
    return S, sS, E, sig, marg


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    print(f"Bell/CHSH, |Φ+>, {4 * N} pairs per model (N = {N} per setting pair)", flush=True)
    res = []

    rng = np.random.default_rng(17)
    S, sS, E, sig, marg = chsh(lambda a, b, n: collapse(a, b, n, rng), N)
    print(f"[collapse model, entangle.rs logic] S = {S:.3f} ± {sS:.3f}  (Rust validation.rs, 40,000 pairs: 2.840)", flush=True)

    t = time.time()
    seeds = iter(range(100, 200))
    store = {}

    def fn(a, b, n):
        ra, rb, X0, sm = pilot(a, b, n, next(seeds))
        store[a, b] = (ra, rb, X0)
        return ra, rb
    S, sS, E, sig, marg = chsh(fn, N)
    dt = time.time() - t
    print(f"[PILOT, PBSs at the same z] S = {S:.4f} ± {sS:.4f}   theory 2√2 = {2*np.sqrt(2):.4f}   [{dt:.0f}s]", flush=True)
    for (a, b), e in E.items():
        th = np.cos(2 * (a - b) * DEG)
        ok = abs(e - th) < 3 * sig[a, b] + 1e-9
        print(f"   E({a:4.1f}°,{b:4.1f}°) = {e:+.4f} ± {sig[a, b]:.4f}   theory {th:+.4f}   {'PASS' if ok else 'FAIL'}"
              f"   P(A+) {marg[a, b][0]:.4f}  P(B+) {marg[a, b][1]:.4f}", flush=True)
        res.append(ok)
    okS = abs(S - 2 * np.sqrt(2)) < 3 * sS
    print(f"   CHSH: |S − 2√2| = {abs(S - 2*np.sqrt(2)):.4f} vs 3σ = {3*sS:.4f}  {'PASS' if okS else 'FAIL'}", flush=True)
    res.append(okS)
    # no signalling: A's marginal must not depend on b, B's must not depend on a
    ns = []
    for a in A_SET:
        p1, p2 = marg[a, B_SET[0]][0], marg[a, B_SET[1]][0]
        ns.append((f"P(A+|a={a:g}) b={B_SET[0]:g} vs {B_SET[1]:g}", p1, p2))
    for b in B_SET:
        p1, p2 = marg[A_SET[0], b][1], marg[A_SET[1], b][1]
        ns.append((f"P(B+|b={b:g}) a={A_SET[0]:g} vs {A_SET[1]:g}", p1, p2))
    for name, p1, p2 in ns:
        sd = np.sqrt(0.5 / N)
        ok = abs(p1 - p2) < 3 * sd and abs(p1 - 0.5) < 3 * sd / np.sqrt(0.5) and abs(p2 - 0.5) < 3 * sd / np.sqrt(0.5)
        print(f"   no-signalling {name}: {p1:.4f} vs {p2:.4f} (diff {p1-p2:+.4f}, 3σ {3*sd:.4f})  {'PASS' if ok else 'FAIL'}", flush=True)
        res.append(ok)

    # E(0, b) curve
    M = max(N // 5, 1000)
    print(f"E(0°, b) curve, {M} pairs per point:", flush=True)
    for b in (0, 15, 30, 45, 60, 75, 90):
        ra, rb, _, _ = pilot(0.0, b, M, 300 + b)
        e = np.mean(ra * rb); sd = np.sqrt(max(1 - e ** 2, 1 / M) / M)
        th = np.cos(2 * b * DEG)
        ok = abs(e - th) < max(3 * sd, 3 / M)
        res.append(ok)
        print(f"   b = {b:2d}°  E = {e:+.4f} ± {sd:.4f}  theory {th:+.4f}  {'PASS' if ok else 'FAIL'}", flush=True)

    # A's PBS before B's (A at z=60, B at z=100)
    t = time.time()
    seeds2 = iter(range(400, 500))
    S2, sS2, E2, sig2, _ = chsh(lambda a, b, n: pilot(a, b, n, next(seeds2), 60.0, 100.0)[:2], M)
    ok = abs(S2 - 2 * np.sqrt(2)) < 3 * sS2
    res.append(ok)
    print(f"[PILOT, A's PBS at z=60, B's at z=100] S = {S2:.4f} ± {sS2:.4f}  {'PASS' if ok else 'FAIL'}  [{time.time()-t:.0f}s]", flush=True)

    # dependence of one pair's outcome on the distant setting (same initial ball positions)
    print("outcome of A for the SAME initial positions when only Bob's setting changes (22.5° -> 67.5°):", flush=True)
    for zB, label in ((60.0, "simultaneous PBSs"), (100.0, "A measured first ")):
        X0 = None
        ra1, rb1, X0, _ = pilot(0.0, 22.5, M, 900, 60.0, zB)
        ra2, rb2, _, _ = pilot(0.0, 67.5, M, 900, 60.0, zB, X0=X0)
        fa, fb = (ra1 != ra2).mean(), (rb1 != rb2).mean()
        print(f"   {label}: A's outcome changes in {100*fa:.2f}% of pairs, B's in {100*fb:.1f}%", flush=True)
    print(f"RESULT: {sum(res)}/{len(res)} PASS   total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

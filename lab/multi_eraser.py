"""Delayed-choice quantum eraser (Kim, Yu, Kulik, Shih, Scully, PRL 84, 1, 2000) with the pilot rule.

A pair is born either in region A or in region B of the crystal: (|A>_s|a>_i + |B>_s|b>_i)/√2.
Signal photon: from A (x = −3) or B (x = +3), both narrow beams, reaches the screen at z = 50 and clicks.
Idler photon (its own arm, longer path): beam a or b;
  z = 30  "which-path taps": a -> ½ continues, ½ deflected to D3;  b -> ½ continues, ½ to D4
  z = 60  50:50 BS mixing a and b -> D1, D2 (erases which-path)
  z = 160 detectors.
The idler's "choice" (which detector) is made AFTER the signal photon has already clicked.
Expected: signal clicks show no fringes overall; conditioned on D1 fringes, on D2 anti-fringes,
on D3/D4 no fringes.  In the pilot picture the signal ball lands first; the idler ball's later route is
correlated with where the signal ball landed (nonlocal guidance), never the other way round.
"""
import time
import numpy as np
from multi_core import System, Arm, bs_rules, ALPHA, H, DU

T0 = time.time()
ZS, D, WS = 50.0, 6.0, 1.5
TAP = 0.75


def build():
    sig = Arm("S", ZS).source("A", "s", -D / 2, 0.0, WS, 0.0).source("B", "s", D / 2, 0.0, WS, 0.0)
    x0 = 60 * np.tan(ALPHA)
    idl = Arm("I", 160.0).source("a", "a", -x0, ALPHA, 5.0, 60.0).source("b", "b", x0, -ALPHA, 5.0, 60.0)
    idl.obj(30.0, {"a": [("a", ALPHA, 1 / np.sqrt(2)), ("D3", -TAP, 1j / np.sqrt(2))],
                   "b": [("b", -ALPHA, 1 / np.sqrt(2)), ("D4", TAP, 1j / np.sqrt(2))]}, "taps")
    idl.obj(60.0, bs_rules(inA="a", inB="b", outC="D1", outD="D2"), "BS_E")
    s = 1 / np.sqrt(2)
    return System([sig, idl], ["S", "I"], [(s, [("A", H, None), ("a", H, None)]), (s, [("B", H, None), ("b", H, None)])])


def fringe(x, q):
    """2<e^{-iqx}>: amplitude and phase of the fringe component (q = fringe spatial frequency)."""
    z = 2 * np.mean(np.exp(-1j * q * x))
    return np.abs(z), np.angle(z)


def main():
    N = 8000
    t = time.time()
    sysm = build()
    print(f"delayed-choice quantum eraser: {N} pairs; {sysm.n_balls} wave balls; tables {sysm.build_time:.0f}s", flush=True)
    X0 = sysm.sample(N, np.random.default_rng(2000))
    Xf = sysm.run(X0)
    xs = Xf[:, 0]
    port = sysm.ports(Xf)[1]
    print(f"  {N} pairs guided [{time.time()-t:.0f}s]", flush=True)
    # theory: Born rule on the same waves, P(x_s, D) = Re Σ c_k* c_l Q^s_kl(x_s) ∫_D Q^i_kl
    xg = np.linspace(-40, 40, 1601); dx = xg[1] - xg[0]
    Qs, _ = sysm.overlaps(0, xg, ZS, False)
    arm = sysm.arms["I"]
    xi = np.concatenate([np.linspace(ch.center(160) - ch.Uh, ch.center(160) + ch.Uh, 2000) for ch in arm.out])
    xi = np.unique(xi)
    Qi, _ = sysm.overlaps(1, xi, 160.0, False)
    pi = arm.port(xi)
    w = np.gradient(xi)
    th = {}
    for Dn in ("D1", "D2", "D3", "D4"):
        I_D = (Qi[pi == Dn] * w[pi == Dn, None, None]).sum(0)
        th[Dn] = np.real((sysm.P[None] * Qs * I_D[None]).sum((1, 2)))
    tot = sum(th.values())
    norm = tot.sum() * dx
    q = 2 * np.pi * D / (ZS * 1.0)      # nominal fringe frequency 2π d/(λ z)
    # refine q on the theory D1 pattern (near field): maximize fringe amplitude
    qs = np.linspace(0.7 * q, 1.3 * q, 241)
    amp = [abs((th["D1"] * np.exp(-1j * qq * xg)).sum()) for qq in qs]
    q = qs[int(np.argmax(amp))]

    def th_fringe(p):
        z = 2 * (p * np.exp(-1j * q * xg)).sum() / p.sum()
        return np.abs(z), np.angle(z)

    res = []
    rows = [("all signal clicks", np.ones(N, bool), tot)] + [(f"signal | idler at {Dn}", port == Dn, th[Dn]) for Dn in ("D1", "D2", "D3", "D4")]
    phases = {}
    print(f"  fringe frequency q = {q:.4f} rad/λ (period {2*np.pi/q:.2f} λ); fringe amplitude V = 2|<e^(−iqx)>|", flush=True)
    for name, sel, p in rows:
        n = sel.sum()
        V, ph = fringe(xs[sel], q)
        Vt, pht = th_fringe(p)
        sd = np.sqrt(2.0 / max(n, 1))
        frac, frac_t = n / N, p.sum() * dx / norm
        okV = abs(V - Vt) < 3 * sd
        okF = abs(frac - frac_t) < 3 * np.sqrt(frac_t * (1 - frac_t) / N) + 1e-9
        h, e = np.histogram(xs[sel], np.linspace(-30, 30, 25))
        pt = np.interp(0.5 * (e[1:] + e[:-1]), xg, p)
        corr = np.corrcoef(h, pt)[0, 1]
        res += [okV, okF]
        phases[name] = (ph, pht)
        print(f"  {name:22s}: {n:5d} clicks ({frac:.4f}, theory {frac_t:.4f}) {'PASS' if okF else 'FAIL'}  "
              f"V = {V:.3f} ± {sd:.3f}  theory {Vt:.3f}  {'PASS' if okV else 'FAIL'}   phase {ph:+.2f} (theory {pht:+.2f})   "
              f"histogram vs theory corr {corr:.3f}", flush=True)
    d = (phases["signal | idler at D1"][0] - phases["signal | idler at D2"][0]) % (2 * np.pi)
    ok = abs(d - np.pi) < 0.3
    res.append(ok)
    print(f"  D1 fringes vs D2 anti-fringes: phase difference {d:.2f} rad (theory π)  {'PASS' if ok else 'FAIL'}", flush=True)
    print(f"RESULT: {sum(res)}/{len(res)} PASS   total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

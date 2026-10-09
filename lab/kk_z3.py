"""kk_z3.py - test 4: the user's 'three': a hidden circle of N = 3 sites (120 degrees apart).
Hidden momenta live in Z3 = {0, +1, -1}. Each ball's Aharonov-Bohm phase is MEASURED on a ring with the
zig-zag rule (fringe shift, as in kk_force.py) for a twist of exactly one hidden site (a 120 deg turn,
a pure permutation of the 3 sites) and for a generic fractional twist. Balls do not interact here, so
a composite's phase is the sum of its balls' phases: neutral = phase 0 mod 2 pi.
Writes lab/out/kk_z3.json."""
import json, time, itertools
import numpy as np
import kk_core as K
from zz_core import bloch_batch

t0 = time.time()
L, N, th, p, sig = 400, 3, 0.2, 0.5, 10.0
xa = np.arange(L); win = np.abs(xa - L // 2) < 30


def E(n, k):
    return np.max(-np.angle(np.linalg.eigvals(bloch_batch(np.array([[k, 0, 2 * np.pi * n / N]]), th)[0])) / 2)


def fringe_phase(n, Dsites):
    """AB phase of one ball: phase of the interference term at the meeting point (kk_core.ring_ab)"""
    v = (E(n, p + 1e-4) - E(n, p - 1e-4)) / 2e-4
    ticks = int(round((L / 2) / abs(v) / 2)) * 2
    Dx = np.zeros((L, 1)); Dx[L // 4] = Dsites          # slide by Dsites hidden SITES on one link
    _, po, pf = K.ring_ab(L, N, th, n, Dx, p=p, sigma=sig, ticks=ticks)
    return po, float(v), pf


out = {"rest_energy": {}, "phases": {}}
for n in (0, 1, -1):
    out["rest_energy"][n] = dict(E_walk=float(E(n, 0)), E_continuum=float(np.sqrt(th**2 + (2 * np.pi * n / N)**2)))
    print(f"N=3, n={n:+d}: rest energy on the walk {E(n,0):.4f} (continuum sqrt(m^2+(n/R)^2) = "
          f"{np.sqrt(th**2+(2*np.pi*n/N)**2):.4f}) -> a 3-site circle is far from the continuum")
ph = {}
for twist in (1.0, 0.37):
    for n in (0, 1, -1):
        p0, v, _ = fringe_phase(n, 0.0)
        p1, _, _ = fringe_phase(n, twist)
        d = float(np.angle(np.exp(1j * (p1 - p0))))
        ph[(twist, n)] = d
        print(f"twist {twist:.2f} site(s) = {360*twist/3:6.1f} deg, n={n:+d}: measured AB phase {d:+.4f} rad "
              f"(|2 pi n D/N| = {abs(2*np.pi*n*twist/N):.4f}); group velocity {v:+.3f}")
        out["phases"][f"twist={twist} n={n}"] = d


def label(c):
    return "(" + ",".join(f"{x:+d}" for x in c) + ")"


comp = {}
for twist in (1.0, 0.37):
    tab = {}
    for size in (1, 2, 3):
        neutral = []
        for c in itertools.combinations_with_replacement((0, 1, -1), size):
            tot = sum(ph[(twist, n)] for n in c)
            if abs(np.exp(1j * tot) - 1) < 1e-3:
                neutral.append(label(c))
        allc = list(itertools.combinations_with_replacement((0, 1, -1), size))
        tab[size] = dict(n_combinations=len(allc), neutral=neutral)
        print(f"twist {twist:.2f}: {size} ball(s): {len(neutral)}/{len(allc)} neutral: {', '.join(neutral)}")
    comp[str(twist)] = tab
out["composites"] = comp
# ordered states (balls distinguishable): count of the 3^k states with sum = 0 mod 3
out["ordered_fraction_neutral"] = {k: sum(1 for c in itertools.product((0, 1, -1), repeat=k) if sum(c) % 3 == 0) / 3**k
                                   for k in (1, 2, 3)}
print("ordered (distinguishable) fraction neutral:", out["ordered_fraction_neutral"])
# the user's statement: remove one ball from the neutral (0,+1,-1)
for miss in (0, 1, -1):
    rest = [n for n in (0, 1, -1) if n != miss]
    tot = sum(ph[(1.0, n)] for n in rest)
    print(f"(0,+1,-1) without {miss:+d}: left {rest}, AB phase {np.angle(np.exp(1j*tot)):+.4f} rad")
    out[f"missing {miss}"] = float(np.angle(np.exp(1j * tot)))
out["time_s"] = time.time() - t0
json.dump(out, open("out/kk_z3.json", "w"), indent=1, default=str)
print("time %.1f s" % out["time_s"])

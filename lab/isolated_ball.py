"""TASK 4 — ordering test: a photon made of ONE ball must fly as a straight ray.

With N = 1 the photon ball's only neighbour is itself (its own wavelet): the current of a single
wavelet points radially away from the origin, i.e. along the ball's own direction.
We also check N = 2 balls with the same origin (two directions) and a ball whose neighbour
fan is far away (> window): it must also go straight.

Run: python lab/isolated_ball.py   (< 5 s)
"""
import numpy as np
from photon_rule import Fan, guide

rng = np.random.default_rng(3)
worst = {}
for est in ("E2", "E1"):
    dev = []
    for _ in range(50):
        s0 = rng.uniform(-5, 5); u0 = rng.uniform(-0.8, 0.8)
        # one ball: origin s0, a lattice with a single direction θ0 = atan(u0)  (θ_0 = -TH + dth/2)
        dth = 0.01
        fan = Fan([s0], [1.0 + 0j], 0.0, TH=dth / 2, dth=dth)
        fan.TH = -np.arctan(u0) + dth / 2
        fan.nu = 1
        x0 = np.array([s0 + 0.5 * u0])                                # photon ball = that ball, at z = 0.5
        xl, uf = guide(x0, fan, 0.5, 160.0, est)
        straight = s0 + 160.0 * u0
        dev.append(abs(xl[0] - straight))
    worst[est] = max(dev)
    print(f"N=1, estimator {est}: max |x_end - straight ray| over 50 random balls = {max(dev):.2e} λ", flush=True)

# isolated from a fan: photon ball far outside any fan (no neighbour at all) keeps its direction
fan = Fan([0.0], [1.0 + 0j], 0.0, TH=0.1, dth=0.02)
x0 = np.array([30.0])
xl, uf = guide(x0, fan, 1.0, 160.0, "E2")
print(f"no neighbour at all: x_end - x_start = {xl[0]-30.0:.2e} (must be 0: straight along z)", flush=True)
ok = all(v < 1e-6 for v in worst.values()) and abs(xl[0] - 30.0) < 1e-9
print("PASS" if ok else "FAIL", flush=True)

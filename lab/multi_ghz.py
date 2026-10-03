"""GHZ with three photons (Pan, Bouwmeester, Daniell, Weinfurter, Zeilinger, Nature 403, 515, 2000)
with the multi-photon pilot rule.

State (HHH + VVV)/√2, one photon per arm, each arm a PBS measuring X (±45°) or Y (R/L circular):
the outcome is the port where the photon ball lands.  Quantum: E(XYY) = E(YXY) = E(YYX) = −1 and
E(XXX) = +1.  Local realism: X1Y2Y3 · Y1X2Y3 · Y1Y2X3 = X1X2X3 (each Y appears twice), so it must
predict E(XXX) = (−1)³ = −1: an all-or-nothing contradiction.
"""
import time
import numpy as np
from multi_core import System, pbs_arm, H, V, basis

T0 = time.time()
X45, Y = basis("45"), basis("R")


def theory(bases, state):
    ops = []
    for e1 in bases:
        e2 = np.array([-np.conj(e1[1]), np.conj(e1[0])])
        ops.append(np.outer(e1, np.conj(e1)) - np.outer(e2, np.conj(e2)))
    O = np.kron(np.kron(ops[0], ops[1]), ops[2])
    return np.real(np.conj(state) @ O @ state)


def main():
    N = 4000
    s = 1 / np.sqrt(2)
    state = s * (np.kron(np.kron(H, H), H) + np.kron(np.kron(V, V), V))
    res, E = [], {}
    print(f"GHZ (HHH+VVV)/√2, {N} triples per setting", flush=True)
    for i, name in enumerate(("XYY", "YXY", "YYX", "XXX")):
        t = time.time()
        bases = [X45 if c == "X" else Y for c in name]
        arms = [pbs_arm(f"A{j}", bases[j]) for j in range(3)]
        sysm = System(arms, ["A0", "A1", "A2"], [(s, [("s", H, None)] * 3), (s, [("s", V, None)] * 3)])
        X0 = sysm.sample(N, np.random.default_rng(20 + i))
        p = sysm.ports(sysm.run(X0))
        r = np.prod([np.where(pp == "+", 1, -1) for pp in p], 0)
        E[name] = r.mean()
        th = theory(bases, state)
        sd = np.sqrt(max(1 - E[name] ** 2, 1 / N) / N)
        ok = abs(E[name] - th) < max(3 * sd, 3 / N)
        res.append(ok)
        marg = [(pp == "+").mean() for pp in p]
        print(f"  E({name}) = {E[name]:+.4f}  theory {th:+.0f}  (wrong-sign triples {np.mean(r != np.sign(th)):.4f}; "
              f"single-photon marginals {marg[0]:.3f} {marg[1]:.3f} {marg[2]:.3f})  {'PASS' if ok else 'FAIL'}  [{time.time()-t:.0f}s]", flush=True)
    lr = E["XYY"] * E["YXY"] * E["YYX"]
    print(f"  local realism, from the three measured XYY-type values: E(XXX) = {np.sign(lr):+.0f}; "
          f"measured E(XXX) = {E['XXX']:+.4f}  -> contradiction {'YES' if E['XXX'] * np.sign(lr) < 0 else 'NO'}", flush=True)
    res.append(E["XXX"] * np.sign(lr) < 0)
    print(f"RESULT: {sum(res)}/{len(res)} PASS   total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

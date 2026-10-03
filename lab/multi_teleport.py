"""Teleportation (Bouwmeester et al. 1997) and entanglement swapping (Pan et al. 1998) with linear optics
and the multi-photon pilot rule.

Bell-state measurement (BSM) = two photons enter a 50:50 BS from ports a and b of the same arm; success =
the two photon balls land in DIFFERENT output ports (c and d): that projects the pair on |Ψ−> (the only
antisymmetric polarization state), probability 1/4 with linear optics.
Teleportation: photon 1 in state p, photons 2-3 in |Ψ−>; after a BSM success photon 3 is in state p.
Measured by a PBS on photon 3 in the basis {p, p⊥}: fidelity F = P(port + | success). Theory 1; the
best classical (measure-and-resend) strategy gives 2/3.
Swapping: |Ψ−>_12 |Ψ−>_34, BSM on 2-3; photons 1 and 4 never met, yet after success they are in |Ψ−>:
F = (1 − E_ZZ − E_XX − E_YY)/4 from the three correlations (theory 1, separable states ≤ 1/2).
"""
import time
import numpy as np
from multi_core import System, pbs_arm, bs_arm, H, V, basis, lin

T0 = time.time()
STATES = {"H": H, "V": V, "+45": basis("45"), "-45": basis("-45"), "R": basis("R"),
          "L": np.array([1, -1j]) / np.sqrt(2)}


def sym(terms, i, j):
    """Bosonic symmetrization of photons i and j (same arm): add each term with their specs swapped."""
    out = []
    for c, sp in terms:
        out.append((c / np.sqrt(2), sp))
        sp2 = list(sp); sp2[i], sp2[j] = sp[j], sp[i]
        out.append((c / np.sqrt(2), sp2))
    return out


def teleport(p, N, seed):
    s = 1 / np.sqrt(2)
    base = [(s, [("a", p, None), ("b", H, None), ("s", V, None)]),
            (-s, [("a", p, None), ("b", V, None), ("s", H, None)])]
    sysm = System([bs_arm("BS"), pbs_arm("B", p)], ["BS", "BS", "B"], sym(base, 0, 1))
    X0 = sysm.sample(N, np.random.default_rng(seed))
    p1, p2, p3 = sysm.ports(sysm.run(X0))
    succ = p1 != p2
    return succ, p3 == "+"


def swap(e1, e4, N, seed):
    h = 0.5
    base = []
    for c, (q1, q2), (q3, q4) in ((+h, (H, V), (H, V)), (-h, (H, V), (V, H)), (-h, (V, H), (H, V)), (+h, (V, H), (V, H))):
        base.append((c, [("s", q1, None), ("a", q2, None), ("b", q3, None), ("s", q4, None)]))
    sysm = System([pbs_arm("A", e1), bs_arm("BS"), pbs_arm("D", e4)], ["A", "BS", "BS", "D"], sym(base, 1, 2))
    X0 = sysm.sample(N, np.random.default_rng(seed))
    p1, p2, p3, p4 = sysm.ports(sysm.run(X0))
    succ = p2 != p3
    r = np.where(p1 == "+", 1, -1) * np.where(p4 == "+", 1, -1)
    return succ, r, len(base) * 2


def main():
    res = []
    N = 3000
    print(f"Teleportation, {N} events per input state (3 photon balls, 4 symmetrized terms)", flush=True)
    Fs, Ps = [], []
    for i, (name, p) in enumerate(STATES.items()):
        t = time.time()
        succ, plus = teleport(p, N, 30 + i)
        ps = succ.mean(); F = plus[succ].mean(); Ff = plus[~succ].mean(); Fa = plus.mean()
        ns = succ.sum()
        sdF = np.sqrt(max(F * (1 - F), 1 / ns) / ns)
        okF = abs(F - 1) < max(3 * sdF, 3 / ns)
        okP = abs(ps - 0.25) < 3 * np.sqrt(0.25 * 0.75 / N)
        res += [okF, okP]
        Fs.append(F); Ps.append(ps)
        print(f"  input {name:>3s}: BSM success {ps:.4f} (theory 0.25) {'PASS' if okP else 'FAIL'};  fidelity | success "
              f"F = {F:.4f} ({ns} events, theory 1, classical bound 2/3) {'PASS' if okF else 'FAIL'};  "
              f"F | BSM failed {Ff:.3f} (theory 1/3);  F without BSM {Fa:.3f} (theory 1/2)  [{time.time()-t:.0f}s]", flush=True)
    print(f"  average: success {np.mean(Ps):.4f}, fidelity {np.mean(Fs):.4f}", flush=True)

    M = 6000
    print(f"Entanglement swapping, {M} events per basis (4 photon balls, 8 symmetrized terms)", flush=True)
    E = {}
    for i, (nm, e) in enumerate((("ZZ", H), ("XX", basis("45")), ("YY", basis("R")))):
        t = time.time()
        succ, r, K = swap(e, e, M, 60 + i)
        E[nm] = r[succ].mean(); ns = succ.sum()
        sd = np.sqrt(max(1 - E[nm] ** 2, 1 / ns) / ns)
        ok = abs(E[nm] + 1) < max(3 * sd, 3 / ns)
        okP = abs(succ.mean() - 0.25) < 3 * np.sqrt(0.25 * 0.75 / M)
        res += [ok, okP]
        print(f"  E_{nm}(1,4) | success = {E[nm]:+.4f} ({ns} events; theory −1) {'PASS' if ok else 'FAIL'};  "
              f"without BSM {r.mean():+.4f} (theory 0);  success {succ.mean():.4f} (0.25) {'PASS' if okP else 'FAIL'}  [{time.time()-t:.0f}s]", flush=True)
    F = (1 - E["ZZ"] - E["XX"] - E["YY"]) / 4
    ok = F > 0.97
    res.append(ok)
    print(f"  fidelity of photons 1-4 with |Ψ−> after swapping: F = {F:.4f} (theory 1; > 1/2 = entangled) {'PASS' if ok else 'FAIL'}", flush=True)
    print(f"RESULT: {sum(res)}/{len(res)} PASS   total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

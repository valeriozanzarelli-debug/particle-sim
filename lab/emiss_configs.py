"""emiss_configs.py - configurations of the runs (development ones are marked dev_; production ones are listed in
lab/out/emiss_prereg.json and must not change after it is committed)."""
import os
import numpy as np
from emiss_core import Dirac3D, CK

LEVELS = {}


def level(n, st):
    k = (n, st)
    if k not in LEVELS:
        LEVELS[k] = np.load(os.path.join(CK, f"level_{n}_{st.replace('/', '_')}.npy"))
    return LEVELS[k]


def atom(n, a=0.5):
    return lambda: Dirac3D(n, a)


def base(name, mode, n_atom, T, K=1.0, nL=128, b=4.0, sigma=1.0, **kw):
    c = dict(name=name, mode=mode, T=T, K=K, nL=nL, b=b, sub=int(round(b / 0.5)), sigma=sigma, sponge=12, probe=40)
    if mode != "light":
        c["atom_factory"] = atom(n_atom)
    c.update(kw)
    return c


CONFIGS = {}


def add(c):
    CONFIGS[c["name"]] = c


# --------------------------------------------------------------- development (declared)
add(base("dev_qe48", "qe", 48, 600,
         psi0_factory=lambda at: level(48, "2p3/2+3").copy(),
         excited=lambda at: [level(48, "2p3/2+3")],
         ground=lambda at: level(48, "1s+"),
         channels=lambda at: [level(48, "1s+")]))
add(base("dev_sc48", "sc", 48, 600,
         psi0_factory=lambda at: level(48, "2p3/2+3").copy(),
         excited=lambda at: [level(48, "2p3/2+3")],
         ground=lambda at: level(48, "1s+")))
add(base("dev_qe48_K10", "qe", 48, 600, K=10.0,
         psi0_factory=lambda at: level(48, "2p3/2+3").copy(),
         excited=lambda at: [level(48, "2p3/2+3")],
         ground=lambda at: level(48, "1s+"),
         channels=lambda at: [level(48, "1s+")]))
add(base("dev_light", "light", None, 1200, I0=1.0, omega=0.0770, ramp=100.0))

for tag, kw in (("A", {}), ("B", dict(centered=True)), ("C", dict(centered=True, transverse=True)), ("D", dict(transverse=True))):
    add(base("dev_qe48" + tag, "qe", 48, 600, K=1.0,
             psi0_factory=lambda at: level(48, "2p3/2+3").copy(),
             excited=lambda at: [level(48, "2p3/2+3")],
             ground=lambda at: level(48, "1s+"),
             channels=lambda at: [level(48, "1s+")], **kw))

for tag, kw in (("B", dict(centered=True)), ("C", dict(centered=True, transverse=True))):
    add(base("dev_long48" + tag, "qe", 48, 2000, K=1.0,
             psi0_factory=lambda at: level(48, "2p3/2+3").copy(),
             excited=lambda at: [level(48, "2p3/2+3")],
             ground=lambda at: level(48, "1s+"),
             channels=lambda at: [level(48, "1s+")], **kw))

# --------------------------------------------------------------- radial form of the rule (any Z, small cells)
from emiss_radial import RadialAtom
RADIAL = {}


def radial_atom(Z, a, rmax):
    k = (Z, a, rmax)
    if k not in RADIAL:
        at = RadialAtom(Z, a, rmax)
        oe, ce, re_ = at.eigen("e", 2)
        og, cg, rg = at.eigen("g", 1)
        RADIAL[k] = (at, ce, cg, oe, og)
    return RADIAL[k]


def radial_cfg(name, Z, a, rmax, T_periods, cells_per_lambda=20, K=1.0, **kw):
    at, ce, cg, oe, og = radial_atom(Z, a, rmax)
    w = oe - og
    lam = 2 * np.pi / w
    sub = int(round(lam / cells_per_lambda / a))
    b = sub * a
    c = dict(name=name, mode="qe", T=T_periods * lam, K=K, nL=128, b=b, sub=sub, sigma=1.0, sponge=12, probe=40,
             centered=True, transverse=True,
             atom_factory=lambda: radial_atom(Z, a, rmax)[0],
             psi0_factory=lambda at_: radial_atom(Z, a, rmax)[1].copy(),
             excited=lambda at_: [radial_atom(Z, a, rmax)[1]],
             ground=lambda at_: radial_atom(Z, a, rmax)[2],
             channels=lambda at_: [radial_atom(Z, a, rmax)[2]], omega=w)
    c.update(kw)
    return c


add(radial_cfg("dev_rad60", 60, 0.05, 40.0, 25))
add(radial_cfg("dev_rad60_c30", 60, 0.05, 40.0, 25, cells_per_lambda=30))
add(radial_cfg("dev_rad60_c40", 60, 0.05, 40.0, 25, cells_per_lambda=40))

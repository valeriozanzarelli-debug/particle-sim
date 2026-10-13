"""fw1atomo_check.py - tool check (not a production run): the rules of fw1atomo_core equal those of
lab/senzadoppi_real.py bit for bit (F0, FW1 with a Coulomb potential), and the period map is unitary.
Run: python lab/fw1atomo_check.py"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import fw1atomo_core as fc
import senzadoppi_real as sr

L = fc.Lattice(12, 0.5, "angolo")
rng = np.random.default_rng(1)
psi = rng.normal(size=(4, 12, 12, 12)) + 1j * rng.normal(size=(4, 12, 12, 12))
for name, mode, h in (("F0", "each", 0.0), ("FW1", "each_w0", 0.4)):
    a = L.rule(name).period(psi)
    b = sr.RealRule(fc.TICKS, 0.5, mode, h, fc.VEC, V=L.V).period(psi)
    print(name, "max |diff| vs senzadoppi_real:", float(np.max(np.abs(a - b))),
          "norm change:", float(abs(np.linalg.norm(a) / np.linalg.norm(psi) - 1)))
L2 = fc.Lattice(12, 0.5, "sito")
print("site-centred V at nucleus (lattice units):", float(L2.V.min()), "expected", -L2.ZA * fc.CUBE_INV_R)

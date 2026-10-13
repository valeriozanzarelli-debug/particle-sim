"""aggmol_posthoc16.py - POST HOC (after the C8 FAIL, not pre-registered): where does the extra binding of the a0/16
level come from? 1s with the same Galerkin field and box (+-6 a0, N = 192) but the exact parabola k^2/(2 m_kin) instead of
the learned T kinetic; plus <T_learned - k^2/2m_kin> on the learned 1s. Writes lab/out/aggmol_posthoc16.json"""
import json, os, numpy as np, aggmol_core as mc
N, h = 192, 256 * mc.TH
mk = np.tan(mc.TH) / mc.TH
f = mc.MultiField(N, h, [(0, 0, 0)], dtype=np.float32)
out = {}
Tl = mc.kinetic_grid(mc.level_rec(256), N, h)
Te = mc.kinetic_exact(N, h) / mk
for nm, T in (("exact_parabola_mkin", Te),):
    m = mc.Mol(T, f, sym=(1, 1, 1))
    X, E, r, conv, its = mc.solve_until(m, mc.guess_mol(N, h, [(0, 0, 0)], "1s")[:, None], 400)
    out[nm] = dict(E_Eh=E[0] / mc.ALPHA ** 2, rel_vs_masskey=E[0] / mc.ALPHA ** 2 / (-0.5 * mk) - 1)
    v = X[:, 0].reshape((N,) * 3); F = np.abs(mc.fftn(v)) ** 2; F /= F.sum()
    d = (Tl - Te) / mc.ALPHA ** 2
    out["<T_learned - k^2/2m>_Eh"] = float(np.sum(F * d))
    k = np.sqrt(2 * mk * Te) * mc.A0                       # |k| in 1/a0
    for lo, hi in ((0, 4), (4, 8), (8, 16), (16, 99)):
        msk = (k >= lo) & (k < hi)
        out[f"shell_{lo}_{hi}_Eh"] = float(np.sum((F * d)[msk]))
print(out)
json.dump(out, open(os.path.join(mc.OUT, "aggmol_posthoc16.json"), "w"), indent=1)
# box check at a0/8 with the exact parabola: +-6 a0 (N = 96) vs +-8 a0 (N = 128)
for N8 in (96, 128):
    h8 = 512 * mc.TH
    f8 = mc.MultiField(N8, h8, [(0, 0, 0)])
    m = mc.Mol(mc.kinetic_exact(N8, h8) / mk, f8, sym=(1, 1, 1))
    X, E, r, conv, its = mc.solve_until(m, mc.guess_mol(N8, h8, [(0, 0, 0)], "1s")[:, None], 200)
    out[f"a0_8_box{N8 // 16}_rel_vs_masskey"] = float(E[0] / mc.ALPHA ** 2 / (-0.5 * mk) - 1)
print(out)
json.dump(out, open(os.path.join(mc.OUT, "aggmol_posthoc16.json"), "w"), indent=1)

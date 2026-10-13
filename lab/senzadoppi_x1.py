"""senzadoppi_x1.py - explorations AFTER the pre-registration (declared, they change no verdict):
 (a) kinetic mass measured without the k^4 bias: fit of E(k) - E(0) to sqrt(m*^2 + k^2) - m* along 26 directions;
 (b) FW1 with the Wilson along the mass axis X (time-reversal friendly?) and h scan of FW1 (0,2 ... 0,8);
 (c) limiting speed of the electron band (group velocity at k = 10 m, theta = 0,02) for F0, PA, FW1.
Output lab/out/senzadoppi_x1.json"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import senzadoppi_core as s
import senzadoppi_r1 as r1
import senzadoppi_r5 as r5

out = {}


def mstar(R, theta, ndir=26):
    D = s.fib_dirs(ndir)
    E0 = s.quasi(R.bloch(np.zeros((1, 3)), theta), R.T)[0]
    e0 = np.min(E0[E0 > 0])
    res = []
    for kk in (theta / 10, theta / 20):
        Ek = s.quasi(R.bloch(kk * D, theta), R.T)
        pos = np.array([np.min(e[e > 0]) for e in Ek])
        # E(k) = sqrt(m^2 + k^2/(m*/m)) ... solve for m* with E(0) = e0: dE = sqrt(e0^2 + e0 k^2 / m*) - e0
        mk = e0 * kk**2 / ((pos - e0 + e0)**2 - e0**2)
        res.append(mk / theta - 1)
    return float(np.max(np.abs(res[1]))), float(np.mean(res[1])), float(np.max(np.abs(res[0] - res[1])))


rules = {k: r1.RULES[k] for k in ("F0", "PA", "PAW")}
rules["FW1"] = r5.FW1(0.4)
out["a_mkin"] = {}
for nm, R in rules.items():
    for th in (0.1, 0.05, 0.02):
        out["a_mkin"][f"{nm}_{th}"] = mstar(R, th)
    v = [out["a_mkin"][f"{nm}_{th}"][1] for th in (0.1, 0.05, 0.02)]
    out["a_mkin"][nm + "_order"] = [float(np.log(v[0] / v[1]) / np.log(2)), float(np.log(v[1] / v[2]) / np.log(2.5))]
    print(nm, [out["a_mkin"][f"{nm}_{th}"] for th in (0.1, 0.05, 0.02)], out["a_mkin"][nm + "_order"], flush=True)


class FW1X(r5.FW1):
    def bloch(self, K, theta):
        K = np.atleast_2d(np.asarray(K, float))
        C = s.coin(theta)[None]
        W = self.tick_trans(K, 0)
        W = s.wilson_links(K, self.hh / 2, gen=s.XCHI, vecs=self.vecs[::-1]) @ C @ \
            s.wilson_links(K, self.hh / 2, gen=s.XCHI, vecs=self.vecs) @ W
        return C @ self.tick_trans(K, 1) @ W


out["b"] = {}
for nm, R in (("FW1X_h0.4", FW1X(0.4)), ("FW1X_h0.1", FW1X(0.1))) + tuple((f"FW1_h{h}", r5.FW1(h)) for h in (0.2, 0.4, 0.6, 0.8)):
    vs = s.valley_scan(R, 0.1, n=24)
    out["b"][nm] = {"E0": vs["E0"], "Epi": vs["Epi"], "iso24": s.isotropy(R, 24, 300)[0], "mkin": mstar(R, 0.1)[1]}
    print(nm, out["b"][nm], flush=True)


def vlim(R, theta=0.02, k=0.2, n=12):
    D = s.fib_dirs(n)
    dk = 1e-4
    Ep = s.quasi(R.bloch((k + dk) * D, theta), R.T)
    Em = s.quasi(R.bloch((k - dk) * D, theta), R.T)
    v = (np.array([np.min(e[e > 0]) for e in Ep]) - np.array([np.min(e[e > 0]) for e in Em])) / (2 * dk)
    vc = k / np.sqrt(k**2 + theta**2)
    return float(np.mean(v / vc - 1)), float(np.min(v / vc - 1)), float(np.max(v / vc - 1))


out["c_vlim_k10m_theta0.02"] = {nm: vlim(R) for nm, R in rules.items()}
print(out["c_vlim_k10m_theta0.02"])
json.dump(out, open(os.path.join(os.path.dirname(__file__), "out", "senzadoppi_x1.json"), "w"), indent=1)

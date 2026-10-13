"""dilat_segno.py - ESPLORATIVO (non pre-registrato, dichiarato): la correzione quantistica cambia segno dove la
curva m/E(p) ha il flesso, p = m/sqrt(2).  Pesi uguali (alpha = pi/4), phi = 0, d = 0:
    Delta_Q ~ -[O/(1+O)] f''(pbar) (p2 - p1)^2 / 8,   f(p) = m/E,  f'' = m (2p^2 - m^2)/E^5,  O = e^{-(dp)^2/4Dl^2}
(per f = 1 - p^2/2m^2 torna Smith-Ahmadi eq. 25 a theta = pi/4).  Scansione in pbar col continuo e con la regola.
Scrive lab/out/dilat_r_segno.json."""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dilat_core import two_packets, mean_rate, smith_ahmadi_gq

HERE = os.path.dirname(os.path.abspath(__file__))
m = 1.0
Dl, dp = 0.05, 0.05
pg = 2 * np.pi * np.fft.fftfreq(2 ** 16, d=4096 / 2 ** 16)
rows = []
for a in [1 / 16, 1 / 32, None]:
    kk = pg if a is None else 2 * np.pi * np.fft.fftfreq(int(1024 / a), d=a)
    kind = "cont" if a is None else "lattice"
    for pb in np.arange(0.0, 2.01, 0.05):
        p1, p2 = pb - dp / 2, pb + dp / 2
        g, g1, g2 = two_packets(kk, np.cos(np.pi / 4), np.sin(np.pi / 4), 0.0, p1, p2, Dl)
        dq = mean_rate(g, kk, a, m, kind) - 0.5 * (mean_rate(g1, kk, a, m, kind) + mean_rate(g2, kk, a, m, kind))
        O = np.exp(-dp ** 2 / (4 * Dl ** 2))
        fpp = m * (2 * pb ** 2 - m ** 2) / (pb ** 2 + m ** 2) ** 2.5
        rows.append(dict(a=a, pbar=pb, dq=dq, formula=-(O / (1 + O)) * fpp * dp ** 2 / 8,
                         SA_NR=float(smith_ahmadi_gq(np.pi / 4, 0.0, p1, p2, Dl))))
cont = [r for r in rows if r["a"] is None]
zero = None
for r0, r1 in zip(cont[:-1], cont[1:]):
    if r0["dq"] > 0 >= r1["dq"] or r0["dq"] < 0 <= r1["dq"]:
        zero = r0["pbar"] + (r1["pbar"] - r0["pbar"]) * r0["dq"] / (r0["dq"] - r1["dq"])
lat = [r for r in rows if r["a"] == 1 / 32]
zl = None
for r0, r1 in zip(lat[:-1], lat[1:]):
    if r0["dq"] * r1["dq"] <= 0 and r0["dq"] != r1["dq"]:
        zl = r0["pbar"] + (r1["pbar"] - r0["pbar"]) * r0["dq"] / (r0["dq"] - r1["dq"])
print("zero del continuo a pbar = %.4f, regola a=1/32 %.4f, flesso m/sqrt2 = %.4f" % (zero, zl, 1 / np.sqrt(2)))
for r in cont[::4]:
    print("pbar %.2f dq %.3e formula %.3e SA %.3e" % (r["pbar"], r["dq"], r["formula"], r["SA_NR"]))
json.dump(dict(Dl=Dl, dp=dp, zero_cont=zero, zero_regola_32=zl, flesso=1 / np.sqrt(2), righe=rows),
          open(os.path.join(HERE, "out", "dilat_r_segno.json"), "w"), indent=1, default=float)

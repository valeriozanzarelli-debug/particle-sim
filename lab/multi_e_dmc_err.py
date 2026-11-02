"""multi_e_dmc_err.py - DIAGNOSI POST HOC delle barre d'errore della DMC (dichiarata, non cambia la lettura pre-registrata).
Perche': la regola pre-registrata (errore = il piu' grande fra 32/16/8 blocchi) usa poche medie, quindi la stima di sigma
e' rumorosa (+-13..25 %), e le corse si fermano appena sigma <= obiettivo: si fermano di preferenza quando la stima e'
uscita BASSA (bias di arresto). Qui: sigma da autocorrelazione integrata (finestra 5 a.u., >> tau_int ~ 0.4-0.9 a.u.) e
da blocchi con 64..256 medie di almeno 2 a.u., sulle serie gia' salvate (nessuna nuova corsa)."""
import os
import numpy as np
from multi_e_core import OUT, load_json, save_json

res = load_json("multi_e_dmc.json")
out = {}
for key in [k for k in res if (k.startswith("He_dmc") or k.startswith("He3S_dmc")) and "_NW" not in k]:
    t = res[key]["tau"]
    z = np.load(os.path.join(OUT, f"multi_e_dmc_{key}.npz"), allow_pickle=True)
    E = np.array(z["Et"])[int(round(20 / t)):]; n = len(E)
    x = E - E.mean(); f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n]; ac /= ac[0]
    Wn = int(round(5 / t)); tint = 2 * np.sum(ac[:Wn + 1]) - 1
    s_ac = float(np.sqrt(x.var() * tint / n))
    sb = []
    for nb in (64, 128, 256):
        L = n // nb
        if L * t >= 2.0:
            sb.append(float(E[:nb * L].reshape(nb, L).mean(1).std(ddof=1) / np.sqrt(nb)))
    out[key] = dict(tau=t, sigma_prereg=res[key]["sigma"], sigma_autocorr=s_ac, tau_int_au=float(tint * t),
                    sigma_blocks_64_256=sb, ratio=s_ac / res[key]["sigma"])
    print(key, {k: (round(v, 7) if isinstance(v, float) else v) for k, v in out[key].items()})
res["err_posthoc"] = out
save_json("multi_e_dmc.json", res)

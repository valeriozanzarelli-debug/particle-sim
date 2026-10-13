"""tunproto_rC.py - corse C: il protocollo sul cammino quantistico in fibra (= la regola a theta = 0,5).
Uso: python3 tunproto_rC.py pulito | mc | rumore      (checkpoint: un JSON per pezzo)."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunproto_core import OUT, rule_run, estimators

PRE = json.load(open(os.path.join(OUT, "tunproto_prereg.json")))
C = PRE["C_fibra"]
th = C["theta"]


def run(V0, N, delta, **kw):
    return rule_run(th, C["E0"], V0, N * th, C["sigma"], delta, "mV", x0=C["x0"], cut=4.0, nt=C["passi"],
                    margin=C["margin"], **kw)


def piece_clean():
    res = {}
    for d in (C["delta"], 0.6):
        for V0 in C["V0"]:
            for N in C["N_siti"]:
                r = run(V0, N, d)
                res["d%.2f_V%.1f_N%d" % (d, V0, N)] = r
                print(d, V0, N, "m tz %.4f ty %.4f | V tz %.4f ty %.4f PT %.2e" % (
                    r["m"]["tz"], r["m"]["ty"], r["V"]["tz"], r["V"]["ty"], r["m"]["PT"]), flush=True)
    # calibrazione: zona dell'orologio senza barriera
    for N in (4, 8):
        res["libero_N%d" % N] = run(0.0, N, C["delta"])
    json.dump(res, open(os.path.join(OUT, "tunproto_rC_pulito.json"), "w"), indent=1)


def _sample_tilt(rng, sz, N, M):
    n0 = rng.binomial(N, (1 + sz) / 2, size=M)
    s = np.clip(2 * n0 / N - 1, -1 + 1e-9, 1 - 1e-9)
    return s


def piece_mc():
    """esperimenti simulati con conteggi finiti (rumore di proiezione dell'ancella), 2000 ripetizioni"""
    cl = json.load(open(os.path.join(OUT, "tunproto_rC_pulito.json")))
    rng = np.random.default_rng(12345)
    M = 2000
    out = {}
    for d, key5 in ((C["delta"], None), (0.6, "delta_0.6")):
        N1 = C["discriminazioni"]["D1"]["N_5sigma_per_punto_per_orologio"] if key5 is None else C[key5]["N_5sigma_D1"]
        N2 = C["discriminazioni"]["D2"]["N_5sigma"] if key5 is None else C[key5]["N_5sigma_D2"]
        N1, N2 = int(np.ceil(N1)), int(np.ceil(N2))
        Ls = np.array([4, 6, 8]) * th
        # D1: pendenza orologio di massa - pendenza orologio di potenziale
        slopes = {}
        errs = {}
        for c in ("m", "V"):
            T = np.zeros((M, 3))
            E = np.zeros((M, 3))
            for j, N in enumerate((4, 6, 8)):
                sz = cl["d%.2f_V2.0_N%d" % (d, N)][c]["sz"]
                s = _sample_tilt(rng, sz, N1, M)
                T[:, j] = -np.arctanh(s) / d
                E[:, j] = 1 / (d * np.sqrt(N1 * (1 - s ** 2)))
            w = 1 / E ** 2
            Lm = np.sum(w * Ls, 1) / np.sum(w, 1)
            sl = np.sum(w * (Ls - Lm[:, None]) * T, 1) / np.sum(w * (Ls - Lm[:, None]) ** 2, 1)
            se = 1 / np.sqrt(np.sum(w * (Ls - Lm[:, None]) ** 2, 1))
            slopes[c], errs[c] = sl, se
        z1 = (slopes["m"] - slopes["V"]) / np.sqrt(errs["m"] ** 2 + errs["V"] ** 2)
        # D2: segno di tau_my (V0 = 2,4, N = 8), N2/2 in base x e N2/2 in base y
        r = cl["d%.2f_V2.4_N8" % d]["m"]
        sxs = 2 * rng.binomial(N2 // 2, (1 + r["sx"]) / 2, size=M) / (N2 // 2) - 1
        sys_ = 2 * rng.binomial(N2 // 2, (1 + r["sy"]) / 2, size=M) / (N2 // 2) - 1
        phi = np.angle(sxs + 1j * sys_)
        ty = -phi / d
        r2 = sxs ** 2 + sys_ ** 2
        vphi = (sys_ ** 2 * 2 * (1 - sxs ** 2) / N2 + sxs ** 2 * 2 * (1 - sys_ ** 2) / N2) / r2 ** 2
        z2 = -ty / (np.sqrt(vphi) / d)
        out["delta%.2f" % d] = dict(N_D1=N1, N_D2=N2, z_D1_mediana=float(np.median(z1)), z_D1_q16=float(np.quantile(z1, 0.16)),
                                    pendenza_m_media=float(slopes["m"].mean()), pendenza_V_media=float(slopes["V"].mean()),
                                    z_D2_mediana=float(np.median(z2)), z_D2_q16=float(np.quantile(z2, 0.16)),
                                    frazione_D2_negativo=float(np.mean(ty < 0)), ty_D2_media=float(ty.mean()))
        print(d, out["delta%.2f" % d], flush=True)
    json.dump(out, open(os.path.join(OUT, "tunproto_rC_mc.json"), "w"), indent=1)


def piece_noise():
    R = C["rumore"]
    cases = {"fase_comune": dict(noise=dict(phase=R["fase_comune"], common=True)),
             "fase_indipendente": dict(noise=dict(phase=R["fase_indipendente"], common=False)),
             "moneta_comune": dict(noise=dict(coin=R["moneta_comune"], common=True)),
             "pdl": dict(pdl=R["pdl"])}
    res = {}
    for name, kw in cases.items():
        t0 = time.time()
        nr = R["realizzazioni"] if "noise" in kw else 1
        rr = {}
        for V0, N in ((2.0, 4), (2.0, 6), (2.0, 8), (2.4, 8)):
            rr["V%.1f_N%d" % (V0, N)] = run(V0, N, C["delta"], nreal=nr, seed=7 + N, **kw)
        Ls = np.array([4, 6, 8]) * th
        rr["pendenza_m"] = float(np.polyfit(Ls, [rr["V2.0_N%d" % N]["m"]["tz"] for N in (4, 6, 8)], 1)[0])
        rr["pendenza_V"] = float(np.polyfit(Ls, [rr["V2.0_N%d" % N]["V"]["tz"] for N in (4, 6, 8)], 1)[0])
        res[name] = rr
        print(name, rr["pendenza_m"], rr["pendenza_V"], rr["V2.4_N8"]["m"]["ty"], rr["V2.4_N8"]["m"]["r"],
              "%.0fs" % (time.time() - t0), flush=True)
    json.dump(res, open(os.path.join(OUT, "tunproto_rC_rumore.json"), "w"), indent=1)


if __name__ == "__main__" and sys.argv[1] != "scan":
    {"pulito": piece_clean, "mc": piece_mc, "rumore": piece_noise}[sys.argv[1]]()


def piece_scan():
    """POST-HOC (non pre-registrato, dichiarato): soglia di rumore tollerabile. Scansione dell'ampiezza del rumore
    di fase comune, di fase indipendente e di moneta; per ogni livello: pendenza di tau_mz, tau_my (V0 = 2,4, N = 8),
    contrasto, probabilita' trasmessa (il rumore che dipende dal tempo non conserva l'energia: 'scalda' il pacchetto
    e lo fa passare sopra la barriera)."""
    res = {}
    Ls = np.array([4, 6, 8]) * th
    for kind, mk in (("fase_comune", lambda s: dict(noise=dict(phase=s, common=True))),
                     ("fase_indipendente", lambda s: dict(noise=dict(phase=s, common=False))),
                     ("moneta_comune", lambda s: dict(noise=dict(coin=s, common=True)))):
        for s in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05):
            rr = {}
            for V0, N in ((2.0, 4), (2.0, 6), (2.0, 8), (2.4, 8)):
                rr["V%.1f_N%d" % (V0, N)] = run(V0, N, C["delta"], nreal=200, seed=11 + N, **mk(s))
            out = dict(pendenza_m=float(np.polyfit(Ls, [rr["V2.0_N%d" % N]["m"]["tz"] for N in (4, 6, 8)], 1)[0]),
                       pendenza_V=float(np.polyfit(Ls, [rr["V2.0_N%d" % N]["V"]["tz"] for N in (4, 6, 8)], 1)[0]),
                       tmy_V24=rr["V2.4_N8"]["m"]["ty"], r_V24=rr["V2.4_N8"]["m"]["r"],
                       PT_N8=rr["V2.0_N8"]["m"]["PT"])
            res["%s_%g" % (kind, s)] = out
            print(kind, s, out, flush=True)
    json.dump(res, open(os.path.join(OUT, "tunproto_rC_scan.json"), "w"), indent=1)


if __name__ == "__main__" and sys.argv[1] == "scan":
    piece_scan()

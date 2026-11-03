"""tavg_f4_judge.py - FASE 4: verdetti come pre-registrati in lab/out/tavg_prereg_f4.json.
Uso: python tavg_f4_judge.py   -> lab/out/tavg_f4_verdetti.json, lab/out/tavg_f4_figura.png"""
import os, sys, json
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_blocks as tb
import tavg_corr as tc

OUT = os.path.join(HERE, "out"); HA_EV = 27.211386245988
PRE = json.load(open(os.path.join(OUT, "tavg_prereg_f4.json")))
R = json.load(open(os.path.join(OUT, "tavg_f4_runs.json")))
F3 = json.load(open(os.path.join(OUT, "tavg_f3_verdetti.json")))
FULL = {"O": (8, "1s2 2s2 2p4", 2, 4), "O+": (8, "1s2 2s2 2p3", 2, 3), "N": (7, "1s2 2s2 2p3", 2, 3), "N+": (7, "1s2 2s2 2p2", 2, 2)}


def cbs(sp):
    e = {X: R[f"{sp}_X{X}"]["Ec_val"] for X in (3, 4, 5)}
    return (125 * e[5] - 64 * e[4]) / 61, e


def main():
    V = {}
    k1 = {}
    for sp, Z in (("controllo_C", "6"), ("controllo_O", "8")):
        inf, e = cbs(sp)
        a = F3["A1"]["valore"][Z]
        ex = a["E_ECG"] - R[f"{sp}_X5"]["E_rohf_blocco"]
        k1[sp] = dict(Ec_X=e, Ec_inf=inf, Ec_esatta=ex, diff=inf - ex)
    V["K1"] = dict(valore=k1, esito="PASS" if all(abs(v["diff"]) <= 1e-3 for v in k1.values()) else "FAIL")
    pr = {}
    for sp, (Z, cfg, ns, npp) in FULL.items():
        inf, e = cbs(sp)
        core = tc._ecg(f"he1S_Z{Z}")[0] - tb.get_ae(Z, "1s2").E
        cv = F3["core_valenza_per_elettrone"][str(Z)]
        ecv = ns * cv["2S"] + npp * cv["2P"]
        pr[sp] = dict(E_HF=tb.get_ae(Z, cfg).E, Ec_val=inf, Ec_val_X=e, Ec_core=core, Ec_cv=ecv, Ec=core + ecv + inf)
    V["previsioni"] = pr
    ip = {}
    for a, b in (("O", "O+"), ("N", "N+")):
        ip[a] = (pr[b]["E_HF"] - pr[a]["E_HF"] + pr[b]["Ec"] - pr[a]["Ec"]) * HA_EV
    c1 = (pr["O"]["Ec"] + 0.257902) / -0.257902
    V["C1"] = dict(valore=pr["O"]["Ec"], rel=c1, esito="PASS" if abs(c1) <= 0.05 else "FAIL")
    c2 = (ip["O"] - 13.618) / 13.618
    V["C2"] = dict(valore=ip["O"], rel=c2, esito="PASS" if abs(c2) <= 0.03 else "FAIL")
    V["C3"] = dict(valore=[ip["N"], ip["O"]], esito="PASS" if ip["N"] > ip["O"] else "FAIL")
    V["C4"] = dict(valore=ip["N"] - ip["O"], esito="PASS" if abs(ip["N"] - ip["O"] - 0.916) <= 0.3 else "FAIL")
    c5 = (ip["N"] - 14.534) / 14.534
    V["C5"] = dict(valore=ip["N"], rel=c5, esito="PASS" if abs(c5) <= 0.03 else "FAIL")
    V["R1"] = dict(valore=pr["N"]["Ec"], chiave=-0.1888)
    # R2: how much the declared cv debt weighs (reported, not counted)
    V["R2"] = dict(Ec_O_senza_cv=pr["O"]["Ec"] - pr["O"]["Ec_cv"], rel_senza_cv=(pr["O"]["Ec"] - pr["O"]["Ec_cv"] + 0.257902) / -0.257902,
                   IP_O_senza_cv=ip["O"] - (pr["O+"]["Ec_cv"] - pr["O"]["Ec_cv"]) * HA_EV)
    json.dump(V, open(os.path.join(OUT, "tavg_f4_verdetti.json"), "w"), indent=1)
    for k in ("K1", "C1", "C2", "C3", "C4", "C5"):
        print(k, V[k]["esito"], json.dumps(V[k].get("rel", V[k]["valore"]))[:300])
    print("R2", V["R2"])
    for sp in pr:
        print(sp, {k: (round(v * 1e3, 2) if isinstance(v, float) else {x: round(y * 1e3, 2) for x, y in v.items()}) for k, v in pr[sp].items() if k != "E_HF"})
    figure(V, k1)


def figure(V, k1):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    pr = V["previsioni"]
    for sp, d in list(pr.items()) + [(k, dict(Ec_val_X=v["Ec_X"], Ec_val=v["Ec_inf"])) for k, v in k1.items()]:
        Xs = [3, 4, 5]; ys = [d["Ec_val_X"][str(x)] if str(x) in d["Ec_val_X"] else d["Ec_val_X"][x] for x in Xs]
        l, = ax[0].plot([x ** -3 for x in Xs], np.array(ys) * 1e3, "o-", label=sp)
        ax[0].plot([0], [d["Ec_val"] * 1e3], "*", color=l.get_color(), ms=10)
    for k, v in k1.items():
        ax[0].plot([0], [v["Ec_esatta"] * 1e3], "kx", ms=9)
    ax[0].set_xlabel("X^-3 (X = l massimo)"); ax[0].set_ylabel("E_c valenza [mEh]")
    ax[0].set_title("a) CCSD(T) nel blocco, base crescente (x = esatti ECG a 3 e-)"); ax[0].legend(fontsize=7)
    names = ["O", "O+", "N", "N+"]
    b0 = np.array([pr[n]["Ec_core"] for n in names]) * 1e3; b1 = np.array([pr[n]["Ec_cv"] for n in names]) * 1e3
    b2 = np.array([pr[n]["Ec_val"] for n in names]) * 1e3
    ax[1].bar(names, b0, label="1s2"); ax[1].bar(names, b1, bottom=b0, label="core-valenza (debito)")
    ax[1].bar(names, b2, bottom=b0 + b1, label="valenza CCSD(T)")
    ax[1].plot(["O"], [-257.9], "k*", ms=12, label="O chiave"); ax[1].plot(["N"], [-188.8], "kx", ms=10, label="N riportato")
    ax[1].set_ylabel("E_c [mEh]"); ax[1].legend(fontsize=7)
    ax[1].set_title(f"b) cieco: E_c(O) {pr['O']['Ec'] * 1e3:.1f} mEh, IP(O) {V['C2']['valore']:.3f} eV")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "tavg_f4_figura.png"), dpi=110)


if __name__ == "__main__":
    main()

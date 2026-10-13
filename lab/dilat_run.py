"""dilat_run.py - corse dello studio dilat_ (regola zig-zag libera, orologio di massa). Uso:
    python3 dilat_run.py E0|E1|E2|E2d|E3|E4|E5|E6|R1|all
Ogni pezzo scrive lab/out/dilat_r_<pezzo>.json (checkpoint) e dura meno di 10 minuti."""
import json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dilat_core import (two_packets, u_plus_lattice, spinor_state, overlap_spec, tau_spec, mean_rate,
                        run_time_domain, rate_lattice, vgroup_lattice, smith_ahmadi_gq, gq_offset_nr)

HERE = os.path.dirname(os.path.abspath(__file__))
PRE = json.load(open(os.path.join(HERE, "out", "dilat_prereg.json")))
EX = PRE["esami"]
A_LIST = PRE["passi_reticolari"]
LBOX = PRE["scatola"]
m = 1.0


def kgrid(a, L=LBOX):
    nx = int(round(L / a))
    return 2 * np.pi * np.fft.fftfreq(nx, d=a), nx


def save(name, obj):
    with open(os.path.join(HERE, "out", "dilat_r_%s.json" % name), "w") as f:
        json.dump(obj, f, indent=1, default=lambda x: float(np.real(x)) if np.isscalar(x) else str(x))


def rate_rule(g, k, a, T=10.0):
    """ritmo letto dall'orologio a qubit (delta -> 0) sulla regola, stato di energia positiva"""
    u = u_plus_lattice(k, a, m)
    n = int(round(T / a))
    tau, V = tau_spec(spinor_state(g, u), k, n, m, 1e-6, "lattice", a)
    return float(tau[0] / (n * a))


# ======================================================================== E0
def run_E0():
    a, nx = 1 / 8, 4096
    k = 2 * np.pi * np.fft.fftfreq(nx, d=a)
    out = {}
    g, _, _ = two_packets(k, np.cos(np.pi / 4), np.sin(np.pi / 4), 0.0, 0.5, 0.8, 0.3, 256, 256)
    states = {"due_pacchetti": spinor_state(g, u_plus_lattice(k, a, m))}
    gz = np.exp(-(k - 1.0) ** 2) * np.exp(-1j * k * 256)
    pz = np.zeros((nx, 2), complex)
    pz[:, 0] = gz / np.linalg.norm(gz)
    states["chiralita_R"] = pz
    th = m * a
    for name, psi in states.items():
        rec = run_time_domain(psi, nx, a, m, 1e-5, 1600, 200)
        rows = []
        for (n, O, fw, nrm, rho) in rec:
            Os = overlap_spec(psi, k, n, m, 1e-5, "lattice", a)[0]
            tau = -np.angle(O) / 1e-5
            rows.append(dict(n=n, dO=abs(O - Os), dRe=abs(fw.real - n * np.sin(th) ** 2) / n,
                             dIm=abs(fw.imag + np.sin(2 * th) / (2 * a) * tau) / max(abs(fw.imag), 1e-30),
                             tau=tau, Imfw=fw.imag, norma_m1=nrm - 1))
        out[name] = rows
    worst = dict(dO=max(r["dO"] for v in out.values() for r in v),
                 dRe=max(r["dRe"] for v in out.values() for r in v),
                 dIm=max(r["dIm"] for v in out.values() for r in v))
    out["peggiori"] = worst
    s = EX["E0"]["soglie"]
    out["PASS"] = bool(worst["dO"] <= 1e-10 and worst["dRe"] <= 1e-9 and worst["dIm"] <= 1e-6)
    save("E0", out)
    print("E0", worst, out["PASS"])


# ======================================================================== E1
def run_E1():
    cfg = EX["E1"]
    out = {"righe": []}
    for pb in cfg["config"]["pbar"]:
        prev = cfg["previsione"][str(pb)]
        row = dict(pbar=pb, R_cont=prev["rate_cont"], coeff_prev=prev["coeff_a2_previsto"], R={}, err={})
        for a in A_LIST:
            k, nx = kgrid(a)
            g, _, _ = two_packets(k, 1.0, 0.0, 0.0, pb, pb, cfg["config"]["Dl"])
            R = rate_rule(g, k, a)
            row["R"][str(a)] = R
            row["err"][str(a)] = (R - prev["rate_cont"]) / prev["rate_cont"]
            row.setdefault("check_mean_rate", {})[str(a)] = abs(R - mean_rate(g, k, a, m))
        e = [row["err"][str(a)] for a in A_LIST]
        row["rapporti"] = [e[i] / e[i + 1] for i in range(len(e) - 1)]
        row["coeff_misurato"] = e[-1] / A_LIST[-1] ** 2
        out["righe"].append(row)
    okA = all(abs(r["err"][str(A_LIST[-1])]) <= 1e-3 for r in out["righe"])
    okB = all(3.5 <= q <= 4.5 for r in out["righe"] if r["pbar"] > 0 for q in r["rapporti"][1:])
    okC = all(abs(r["coeff_misurato"] / r["coeff_prev"] - 1) <= 0.05 for r in out["righe"] if r["pbar"] > 0)
    out["PASS"] = dict(E1a=okA, E1b=okB, E1c=okC)
    save("E1", out)
    for r in out["righe"]:
        print("E1 p=%.1f R(1/32)=%.8f Rc=%.8f err=%s rapp=%s coeff %.4f/%.4f" % (
            r["pbar"], r["R"][str(A_LIST[-1])], r["R_cont"], ["%.2e" % x for x in r["err"].values()],
            ["%.2f" % x for x in r["rapporti"]], r["coeff_misurato"], r["coeff_prev"]))
    print(out["PASS"])


# ======================================================================== E2 (a, b, c)
def dq_rule(a, al, phi, p1, p2, Dl, d=0.0):
    k, nx = kgrid(a)
    g, g1, g2 = two_packets(k, np.cos(al), np.sin(al), phi, p1, p2, Dl, 0.0, d)
    sup = rate_rule(g, k, a)
    mix = np.cos(al) ** 2 * rate_rule(g1, k, a) + np.sin(al) ** 2 * rate_rule(g2, k, a)
    return sup, mix, sup - mix


def run_E2():
    out = {}
    # E2a
    rows = []
    cfg = EX["E2a"]
    for pr in cfg["previsione"]:
        r = dict(pr)
        r["dq"] = {}
        for a in A_LIST:
            r["dq"][str(a)] = dq_rule(a, pr["alpha"], pr["phi"], pr["p1"], pr["p2"], cfg["config"]["Dl"])[2]
        rows.append(r)
    af = str(A_LIST[-1])
    big = [r for r in rows if abs(r["dq_SA"]) > 1e-9]
    relSA = max(abs(r["dq"][af] - r["dq_SA"]) / abs(r["dq_SA"]) for r in big)
    relC = max(abs(r["dq"][af] - r["dq_cont"]) / abs(r["dq_cont"]) for r in rows if abs(r["dq_cont"]) > 1e-9)
    zero = max(abs(r["dq"][af]) for r in rows if abs(r["phi"] - np.pi / 2) < 1e-9)
    out["E2a"] = dict(righe=rows, max_rel_SA=relSA, max_rel_cont=relC, max_phi_pi2=zero,
                      PASS=bool(relSA <= 0.03 and relC <= 1e-3 and zero <= 1e-12))
    print("E2a rel SA %.3e  rel cont %.3e  phi=pi/2 %.2e" % (relSA, relC, zero), out["E2a"]["PASS"])
    # E2b
    rows = []
    cfg = EX["E2b"]
    for pr in cfg["previsione"]:
        r = dict(pr)
        r["dq"] = {}
        for a in A_LIST:
            r["dq"][str(a)] = dq_rule(a, pr["alpha"], pr["phi"], pr["p1"], pr["p2"], cfg["config"]["Dl"])[2]
        e = [r["dq"][str(a)] - r["dq_cont"] for a in A_LIST]
        r["rapporto_16_32"] = e[-2] / e[-1] if e[-1] != 0 else np.nan
        r["rel_cont"] = abs(e[-1]) / abs(r["dq_cont"])
        rows.append(r)
    relC = max(r["rel_cont"] for r in rows)
    ratios = [r["rapporto_16_32"] for r in rows]
    out["E2b"] = dict(righe=rows, max_rel_cont=relC, rapporti=ratios,
                      PASS=dict(precisione=bool(relC <= 1e-3), ordine=bool(all(3.0 <= q <= 5.0 for q in ratios))))
    print("E2b rel cont %.3e ratios %s" % (relC, ["%.2f" % q for q in ratios]), out["E2b"]["PASS"])
    # E2c
    rows = []
    cfg = EX["E2c"]
    dq0 = max(abs(pr["dq_formula_NR"]) for pr in cfg["previsione"] if pr["dDl"] == 0)
    for pr in cfg["previsione"]:
        r = dict(pr)
        r["dq"] = dq_rule(A_LIST[-1], np.pi / 4, pr["phi"], 0.02, 0.06, cfg["config"]["Dl"], pr["d"])[2]
        r["rel"] = abs(r["dq"] - pr["dq_formula_NR"]) / max(abs(pr["dq_formula_NR"]), 1e-3 * dq0)
        r["rel_cont"] = abs(r["dq"] - pr["dq_cont"]) / max(abs(pr["dq_cont"]), 1e-3 * dq0)
        rows.append(r)
    worst = max(r["rel"] for r in rows)
    out["E2c"] = dict(righe=rows, max_rel=worst, PASS=bool(worst <= 0.03))
    print("E2c worst %.3e" % worst, out["E2c"]["PASS"], [("%.1f" % r["dDl"], "%.3e" % r["dq"]) for r in rows])
    save("E2", out)


# ======================================================================== E2d (nel tempo, regola vera)
def run_E2d():
    a, nx = 1 / 16, 8192
    k = 2 * np.pi * np.fft.fftfreq(nx, d=a)
    al = np.pi / 4
    res = {}
    for name, (ca, sa) in {"sup": (np.cos(al), np.sin(al)), "p1": (1.0, 0.0), "p2": (0.0, 1.0)}.items():
        g, _, _ = two_packets(k, ca, sa, 0.0, 0.5, 0.8, 0.3, 128.0, 128.0)
        psi = spinor_state(g, u_plus_lattice(k, a, m))
        t0 = time.time()
        rec = run_time_domain(psi, nx, a, m, 1e-5, 3200, 320)
        taus = [0.0] + [-np.angle(O) / 1e-5 for (n, O, fw, nrm, rho) in rec]
        ts = [0.0] + [n * a for (n, *_r) in rec]
        rhos = [rho[::4].tolist() for (n, O, fw, nrm, rho) in rec] if name == "sup" else None
        res[name] = dict(t=ts, tau=taus, rho=rhos, sec=time.time() - t0)
        print(name, "%.1f s" % res[name]["sec"])
    sup = np.array(res["sup"]["tau"])
    ts = np.array(res["sup"]["t"])
    rates = np.diff(sup) / np.diff(ts)
    mix = 0.5 * (np.array(res["p1"]["tau"]) + np.array(res["p2"]["tau"]))
    dq_t = rates - np.diff(mix) / np.diff(ts)
    var = (rates.max() - rates.min()) / rates.mean()
    # separazione spaziale dei pacchetti: centro di massa delle due corse singole non serve, si calcola da v
    v1 = vgroup_lattice(np.array([0.5]), a, m)[0]
    v2 = vgroup_lattice(np.array([0.8]), a, m)[0]
    out = dict(t=ts.tolist(), ritmi_finestra=rates.tolist(), dq_finestra=dq_t.tolist(), variazione_rel=var,
               separazione_finale=(v2 - v1) * ts[-1], sigma_x=1 / (np.sqrt(2) * 0.3),
               rho_sup=res["sup"]["rho"], dx_rho=4 * a, PASS=bool(var <= 1e-6))
    save("E2d", out)
    print("E2d variazione %.2e  dq %s  separazione %.1f" % (var, ["%.6f" % x for x in dq_t[:3]], out["separazione_finale"]),
          out["PASS"])


# ======================================================================== E3
def run_E3():
    cfg = EX["E3"]
    rows = []
    for pr in cfg["previsione"]:
        r = dict(pr)
        r["D"] = {}
        for a in A_LIST:
            k, nx = kgrid(a)
            g, _, _ = two_packets(k, 1.0, 0.0, 0.0, pr["pbar"], pr["pbar"], cfg["config"]["Dl"])
            u = u_plus_lattice(k, a, m)
            n = int(round(cfg["config"]["t"] / a))
            O = overlap_spec(spinor_state(g, u), k, n, m, pr["delta"], "lattice", a)[0]
            # NB: u e' l'autovettore della massa m: per delta finito lo stato non e' autostato esatto di m +- delta/2
            # (vale come in Paige-Plato-Kim: stato iniziale preparato una volta, poi le due masse)
            r["D"][str(a)] = 1 + np.angle(O) / (pr["delta"] * n * a)
            r.setdefault("V", {})[str(a)] = abs(O)
        r["rel_cont"] = abs(r["D"][str(A_LIST[-1])] - pr["D_cont_esatto"]) / pr["D_cont_esatto"]
        rows.append(r)
    sel = [r for r in rows if r["pbar"] == 0.05 and r["delta"] == 0.5][0]
    dPPK = abs(sel["D"][str(A_LIST[-1])] - sel["D_PPK_NR"])
    gap = abs(sel["D_PPK_NR"] - sel["D_massa_media_NR"])
    out = dict(righe=rows, max_rel_cont=max(r["rel_cont"] for r in rows), dPPK=dPPK, gap=gap,
               PASS=dict(cont=bool(max(r["rel_cont"] for r in rows) <= 1e-2), PPK=bool(dPPK < 0.2 * gap)))
    save("E3", out)
    for r in rows:
        print("E3 p=%.2f d=%.1f D(1/32)=%.6e cont=%.6e PPK=%.6e media=%.6e V=%.6f" % (
            r["pbar"], r["delta"], r["D"][str(A_LIST[-1])], r["D_cont_esatto"], r["D_PPK_NR"], r["D_massa_media_NR"],
            r["V"][str(A_LIST[-1])]))
    print(out["PASS"], dPPK, gap)


# ======================================================================== E4
def run_E4():
    cfg = EX["E4"]
    tv = cfg["config"]["t"]
    out = {"V_sup": {}, "V_mix": {}}
    for a in A_LIST:
        k, nx = kgrid(a)
        u = u_plus_lattice(k, a, m)
        g, g1, g2 = two_packets(k, np.cos(np.pi / 4), np.sin(np.pi / 4), 0.0, 0.5, 1.0, 0.3)
        ns = [int(round(t / a)) for t in tv]
        Os = overlap_spec(spinor_state(g, u), k, ns, m, 0.05, "lattice", a)
        O1 = overlap_spec(spinor_state(g1, u), k, ns, m, 0.05, "lattice", a)
        O2 = overlap_spec(spinor_state(g2, u), k, ns, m, 0.05, "lattice", a)
        out["V_sup"][str(a)] = np.abs(Os).tolist()
        out["V_mix"][str(a)] = np.abs(0.5 * (O1 + O2)).tolist()
    prev = cfg["previsione"]
    af = str(A_LIST[-1])
    d = max(max(abs(x - y) for x, y in zip(out["V_sup"][af], prev["V_sup"])),
            max(abs(x - y) for x, y in zip(out["V_mix"][af], prev["V_mix"])))
    out["max_diff"] = d
    out["err_per_a"] = {str(a): max(abs(x - y) for x, y in zip(out["V_sup"][str(a)], prev["V_sup"])) for a in A_LIST}
    out["PASS"] = bool(d <= 1e-3)
    # curva fine per la figura (a = 1/16)
    a = 1 / 16
    k, nx = kgrid(a)
    u = u_plus_lattice(k, a, m)
    g, g1, g2 = two_packets(k, np.cos(np.pi / 4), np.sin(np.pi / 4), 0.0, 0.5, 1.0, 0.3)
    tf = np.linspace(0, 800, 161)
    ns = [int(round(t / a)) for t in tf]
    Os = overlap_spec(spinor_state(g, u), k, ns, m, 0.05, "lattice", a)
    O1 = overlap_spec(spinor_state(g1, u), k, ns, m, 0.05, "lattice", a)
    O2 = overlap_spec(spinor_state(g2, u), k, ns, m, 0.05, "lattice", a)
    out["curva"] = dict(t=tf.tolist(), V_sup=np.abs(Os).tolist(), V_mix=np.abs(0.5 * (O1 + O2)).tolist())
    save("E4", out)
    print("E4 max diff %.2e" % d, out["err_per_a"], out["PASS"])
    print("V_sup", ["%.4f" % x for x in out["V_sup"][af]])
    print("V_mix", ["%.4f" % x for x in out["V_mix"][af]])


# ======================================================================== E5
def run_E5():
    cfg = EX["E5"]
    out = {"righe": []}
    for q in cfg["config"]["q"]:
        row = dict(q=q, prev=cfg["previsione"][str(q)], R={})
        for a in A_LIST:
            k, nx = kgrid(a)
            kc = np.pi / a - q
            dk = (k - kc + np.pi / a) % (2 * np.pi / a) - np.pi / a       # distanza periodica dal centro
            g = np.exp(-dk ** 2 / (2 * cfg["config"]["Dl"] ** 2))
            g = g / np.linalg.norm(g)
            row["R"][str(a)] = rate_rule(g, k, a)
        row["rel"] = abs(row["R"][str(A_LIST[-1])] - row["prev"]) / abs(row["prev"])
        out["righe"].append(row)
    # elettrone + doppione in sovrapposizione (q = 0.5): l'orologio quasi si ferma
    a = A_LIST[-1]
    k, nx = kgrid(a)
    ge = np.exp(-(k - 0.5) ** 2 / (2 * 0.05 ** 2))
    dk = (k - (np.pi / a - 0.5) + np.pi / a) % (2 * np.pi / a) - np.pi / a
    gd = np.exp(-dk ** 2 / (2 * 0.05 ** 2))
    g = ge / np.linalg.norm(ge) + gd / np.linalg.norm(gd)
    g = g / np.linalg.norm(g)
    out["elettrone_piu_doppione"] = rate_rule(g, k, a)
    worst = max(r["rel"] for r in out["righe"])
    out["PASS"] = bool(worst <= 2e-3)
    save("E5", out)
    for r in out["righe"]:
        print("E5 q=%.1f R=%s prev=%.6f rel=%.2e" % (r["q"], ["%.6f" % x for x in r["R"].values()], r["prev"], r["rel"]))
    print("elettrone+doppione", out["elettrone_piu_doppione"], out["PASS"])


# ======================================================================== E6
def run_E6():
    cfg = EX["E6"]
    tz = np.array(cfg["config"]["t"])
    prev = np.array(cfg["previsione"]["tau_cont"])
    out = {"tau": {}, "maxdiff": {}}
    for a in A_LIST:
        k, nx = kgrid(a)
        gz = np.exp(-(k - 1.0) ** 2 * cfg["config"]["sigma_x"] ** 2)
        psi0 = np.zeros((nx, 2), complex)
        psi0[:, 0] = gz / np.linalg.norm(gz)
        ns = [int(round(t / a)) for t in tz]
        tau, V = tau_spec(psi0, k, ns, m, 1e-5, "lattice", a)
        out["tau"][str(a)] = tau.tolist()
        out["maxdiff"][str(a)] = float(np.max(np.abs(tau - prev)))
    af = str(A_LIST[-1])
    rm = out["tau"][af][-1] / tz[-1]
    out["ritmo_medio"] = rm
    out["rel_ritmo"] = abs(rm - cfg["previsione"]["ritmo_medio"]) / cfg["previsione"]["ritmo_medio"]
    out["PASS"] = dict(traccia=bool(out["maxdiff"][af] <= 2e-3), ritmo=bool(out["rel_ritmo"] <= 0.02))
    save("E6", out)
    print("E6 maxdiff", out["maxdiff"], "ritmo", rm, cfg["previsione"]["ritmo_medio"], out["PASS"])


# ======================================================================== R1 (rapporto, non esame)
def run_R1():
    """sul reticolo: r^2 + v^2 = 1 - sin^2 ka sin^2 th / sin^2 wa  (tempo proprio != sqrt(1 - v^2))"""
    rows = []
    for a in A_LIST:
        for p in [0.5, 1.0, 2.0, 3.0]:
            r = rate_lattice(np.array([p]), a, m)[0]
            v = vgroup_lattice(np.array([p]), a, m)[0]
            wa = np.arccos(np.cos(m * a) * np.cos(p * a))
            rows.append(dict(a=a, p=p, r=r, v=v, sqrt1mv2=np.sqrt(1 - v * v), r_su_sqrt=r / np.sqrt(1 - v * v),
                             identita=r * r + v * v - (1 - np.sin(p * a) ** 2 * np.sin(m * a) ** 2 / np.sin(wa) ** 2),
                             prev_a2=-0.5 * (v * m * a) ** 2 / (1 - v * v)))
    save("R1", dict(righe=rows))
    for r in rows:
        print("R1 a=%.4f p=%.1f r/sqrt(1-v2)-1=%.3e prev %.3e id %.1e" % (r["a"], r["p"], r["r_su_sqrt"] - 1, r["prev_a2"], r["identita"]))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    fns = dict(E0=run_E0, E1=run_E1, E2=run_E2, E2d=run_E2d, E3=run_E3, E4=run_E4, E5=run_E5, E6=run_E6, R1=run_R1)
    for nme, fn in fns.items():
        if which in (nme, "all"):
            t0 = time.time()
            fn()
            print("[%s] %.1f s" % (nme, time.time() - t0))

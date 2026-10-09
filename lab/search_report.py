"""search_report.py - figure, finalisti e test di predizione (docs/RICERCA_REGOLE.md).

  python lab/search_report.py finalists   sceglie i finalisti (dopo la ricerca), li rivaluta con i campi, figure
  python lab/search_report.py predict     test di predizione (rapporti) + look-elsewhere
  python lab/search_report.py figs        figure: storia della fitness, sacchi di luce
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

TARGETS = {"mu/e": 206.768283, "tau/e": 3477.23, "p/e": 1836.15267, "Delta/N": 1.3121,
           "r_p * m_pi+ (raggio protone x massa pione)": 0.8409 * 139.57039 / 197.3269804}


def load(kind):
    return json.load(open(os.path.join(OUT, f"search_rules_{kind}.json")))


def finalists():
    import search_rules as S
    c, r = load("cma"), load("rand")
    allev = [("cma", e) for e in c["evals"]] + [("rand", e) for e in r["evals"]]
    feas = [(k, e) for k, e in allev if e["n_pass"] == 5]
    best = max(allev, key=lambda ke: ke[1]["F"])
    bestr = max(r["evals"], key=lambda e: e["F"])
    light = [(k, e) for k, e in feas if e["Q1"]["fEM"] >= 0.05]
    fin = {"migliore_assoluto": dict(src=best[0], **best[1])}
    if light:
        bl = max(light, key=lambda ke: ke[1]["F"])
        fin["migliore_con_luce"] = dict(src=bl[0], post_hoc=False, **bl[1])
    else:
        # ANALISI POST-HOC (dichiarata): nessun candidato fattibile ha f_EM >= 0.05; si prende quello con luce
        # (f_EM >= 0.05) che conserva la carica (C3) e supera piu' vincoli, e lo si rilassa piu' a lungo.
        cand = [(k, e) for k, e in allev if e.get("Q1") and e["Q1"].get("fEM", 0) >= 0.05 and e["passed"].get("C3")]
        bl = max(cand, key=lambda ke: ke[1]["rank"])
        fin["migliore_con_luce"] = dict(src=bl[0], post_hoc=True, **bl[1])
    names = list(fin)
    res = S.evaluate_population([fin[k]["u"] for k in names], seed=4242, keep_fields=True, steps=4000)
    summary = {}
    for k, rr in zip(names, res):
        rr["src"] = fin[k]["src"]
        rr["F_in_ricerca"] = fin[k]["F"]
        rr["post_hoc"] = fin[k].get("post_hoc", False)
        rr["passed_in_ricerca"] = fin[k]["passed"]
        with open(os.path.join(OUT, f"search_finalist_{'luce' if 'luce' in k else 'assoluto'}.json"), "w") as f:
            json.dump(rr, f, default=float)
        summary[k] = {kk: rr.get(kk) for kk in ("src", "post_hoc", "P", "passed_in_ricerca", "passed", "F", "F_in_ricerca",
                                                 "S", "vac", "E_star", "Q1_pert")}
        summary[k]["Q1"] = {kk: (rr.get("Q1") or {}).get(kk) for kk in
                            ("E", "Q", "rE", "rq", "fEM", "flux", "Bmax", "Ep", "Epp", "drift", "fmax", "frac_in")}
        summary[k]["spectrum"] = {q: {kk: v for kk, v in s.items() if kk != "dens"} for q, s in (rr.get("spectrum") or {}).items()}
    summary["migliore_casuale"] = {kk: bestr.get(kk) for kk in ("P", "passed", "F", "S")}
    # statistiche della ricerca
    st = {}
    for kind, d in (("cma", c), ("rand", r)):
        ev = d["evals"]
        st[kind] = dict(n=len(ev), n_feasible=sum(e["n_pass"] == 5 for e in ev),
                        n_feasible_light=sum(e["n_pass"] == 5 and e["Q1"]["fEM"] >= 0.05 for e in ev),
                        best_F=max(e["F"] for e in ev),
                        fail_counts={cname: sum(not e["passed"][cname] for e in ev) for cname in ("C0", "C1", "C2a", "C2b", "C3")},
                        prereg_sha256=d.get("prereg_sha256"))
        fe = [e for e in ev if e["n_pass"] == 5]
        if fe:
            st[kind]["g_dei_fattibili"] = [round(e["P"]["g"], 3) for e in fe]
            st[kind]["fEM_dei_fattibili"] = [round(e["Q1"]["fEM"], 4) for e in fe]
            st[kind]["F_dei_fattibili"] = [round(e["F"], 3) for e in fe]
    summary["statistiche"] = st
    with open(os.path.join(OUT, "search_finalists_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))
    fig_objects()


def fig_objects():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import search_rules as S
    files = [f for f in sorted(os.listdir(OUT)) if f.startswith("search_finalist_") and f.endswith(".json")]
    rows = []
    for fn in files:
        d = json.load(open(os.path.join(OUT, fn)))
        if "fields" not in d:
            continue
        panels = [("Q=1", np.array(d["fields"]["dens"]))]
        for q in ("2", "3"):
            s = (d.get("spectrum") or {}).get(q, {})
            if "dens" in s:
                panels.append((f"Q={q}", np.array(s["dens"])))
        panels.append(("campo magnetico B (Q=1)", np.array(d["fields"]["B"])))
        rows.append((fn.replace("search_finalist_", "").replace(".json", ""), d, panels))
    if not rows:
        return
    ncol = max(len(p) for _, _, p in rows)
    fig, axs = plt.subplots(len(rows), ncol, figsize=(3.0 * ncol, 3.0 * len(rows)), squeeze=False)
    ext = [-S.L_BOX / 2, S.L_BOX / 2, -S.L_BOX / 2, S.L_BOX / 2]
    for i, (name, d, panels) in enumerate(rows):
        for j in range(ncol):
            ax = axs[i, j]
            if j >= len(panels):
                ax.axis("off")
                continue
            t, a = panels[j]
            cmap = "RdBu_r" if t.startswith("campo") else "magma"
            vm = np.abs(a).max() + 1e-12
            ax.imshow(a.T, origin="lower", extent=ext, cmap=cmap,
                      vmin=-vm if t.startswith("campo") else 0, vmax=vm)
            ax.set_title(f"{name}\n{t}", fontsize=8)
            ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "search_objects.png"), dpi=80)
    print("figura search_objects.png")


def predict():
    rows = []
    files = [f for f in sorted(os.listdir(OUT)) if f.startswith("search_finalist_") and f.endswith(".json")]
    known = json.load(open(os.path.join(OUT, "search_rules_known.json")))
    cands = [(fn.replace("search_finalist_", "").replace(".json", ""), json.load(open(os.path.join(OUT, fn)))) for fn in files]
    cands += [("nota:" + k, v) for k, v in known.items() if v.get("n_pass") == 5]
    preds = []
    for name, d in cands:
        sp = d.get("spectrum") or {}
        pr = {}
        for q in ("2", "3"):
            if sp.get(q, {}).get("ratio"):
                pr[f"E{q}/E1"] = sp[q]["ratio"]
        m2 = d["vac"].get("m2_charged")
        if m2 and m2 > 0 and d.get("Q1"):
            pr["r_E * m_q (raggio x massa del quanto carico)"] = d["Q1"]["rE"] * np.sqrt(m2)
            pr["E1 / (energia minima topologica 4 pi)"] = d["Q1"]["E"] / (4 * np.pi)
        comp = {}
        for pk, pv in pr.items():
            best = min(TARGETS.items(), key=lambda kv: abs(np.log(pv / kv[1])))
            comp[pk] = dict(predetto=pv, bersaglio_piu_vicino=best[0], valore=best[1], scarto_pct=100 * (pv / best[1] - 1))
        preds.append(dict(regola=name, P=d.get("P"), predizioni=comp))
        rows.extend(comp.values())
    npred, ntar = len(rows), len(TARGETS)
    out = dict(predizioni=preds, n_predizioni=npred, n_bersagli=ntar)
    for x in (0.01, 0.05, 0.10):
        p1 = 2 * np.log(1 + x) / np.log(4000)
        out[f"P_caso_almeno_una_entro_{int(x*100)}pct"] = 1 - (1 - p1) ** (npred * ntar)
    out["entro_5pct"] = [r for r in rows if abs(r["scarto_pct"]) < 5]
    out["entro_10pct"] = [r for r in rows if abs(r["scarto_pct"]) < 10]
    with open(os.path.join(OUT, "search_predictions.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


def bag_row(d):
    h = d["hist"]
    t = np.array([r["t"] for r in h])
    EL = np.array([r["EL_disk"] for r in h])
    Ed = np.array([r["E_disk"] for r in h])
    sg = np.array([r["sig_disk_mean"] for r in h])
    row = dict(nome=d["name"], E0=d["E0"])
    if d["cfg"]["init"] == "collide":
        EL0 = d["E0"]  # all'inizio tutta l'energia e' luce (sigma ~ 1e-3)
        tf = 2 * d["cfg"].get("x0", 6.0) + 6.0 + 2 * d["tc"]
        i = int(np.argmin(np.abs(t - tf)))
        row.update(e0=d["cfg"]["e0"], luce_nel_disco_a_t_form=float(EL[i] / EL0),
                   max_sigma_medio_disco=float(sg.max()), t_max_sigma=float(t[sg.argmax()]),
                   t_fine=float(t[-1]), luce_nel_disco_fine=float(EL[-1] / EL0))
        above = t[sg > 0.5]
        row["bolla_sigma>0.5_da_a"] = [float(above.min()), float(above.max())] if len(above) else None
    else:
        tc = d["tc"]
        for k in (1, 5, 10, 25, 50, 100):
            i = int(np.argmin(np.abs(t - k * tc)))
            if abs(t[i] - k * tc) < tc:
                row[f"luce_frac_{k}tc"] = float(EL[i] / max(EL[0], 1e-12))
        row["sigma_medio_disco_inizio_fine"] = [float(sg[0]), float(sg[-1])]
        row["E_disco_inizio_fine"] = [float(Ed[0]), float(Ed[-1])]
        row["t_fine"] = float(t[-1])
    return row


def bag():
    rows = []
    for fn in sorted(os.listdir(OUT)):
        if fn.startswith("search_bag_") and fn.endswith(".json") and fn != "search_bag_mit.json":
            d = json.load(open(os.path.join(OUT, fn)))
            if "hist" in d:
                rows.append(bag_row(d))
    with open(os.path.join(OUT, "search_bag_summary.json"), "w") as f:
        json.dump(rows, f, indent=1)
    for r in rows:
        print(json.dumps(r))


def figs():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.6))
    ax = axs[0]
    for kind, col in (("cma", "C0"), ("rand", "C1")):
        try:
            d = load(kind)
        except FileNotFoundError:
            continue
        F = np.array([e["F"] for e in d["evals"]])
        ax.plot(np.arange(1, len(F) + 1), np.maximum.accumulate(F), color=col, lw=2,
                label=f"{'CMA-ES' if kind == 'cma' else 'casuale'}: migliore finora")
        ax.plot(np.arange(1, len(F) + 1), F, ".", color=col, ms=3, alpha=0.5)
        fe = np.array([e["n_pass"] == 5 and e["Q1"]["fEM"] >= 0.05 for e in d["evals"]])
        ax.plot(np.arange(1, len(F) + 1)[fe], F[fe], "*", color=col, ms=8, mec="k", mew=0.4)
    ax.set_xlabel("valutazioni")
    ax.set_ylabel("fitness F (pre-registrata)")
    ax.set_title("ricerca: CMA-ES vs casuale (stella: fEM>=5%)", fontsize=9)
    ax.legend(fontsize=7)
    # sacchi di luce: bolle pre-formate
    ax = axs[1]
    for n, lab in [("seed_D20_bubble", "dielettrico eps_in/eps_out=20"), ("seed_D400_bubble", "dielettrico 400"),
                   ("seed_U0.6_wgm", "accoppiamento a u (indice 2)"), ("seed_U0.6_wgm_amp0.6", "idem, 4x luce"),
                   ("seed_M8_shell", "guscio con fotone massivo"),
                   ("seed_M8_bagout_dUneg", "Friedberg-Lee: fuori fotone massivo"),
                   ("seed_M8_bagout_dUneg_amp0.6", "idem, 4x luce")]:
        p = os.path.join(OUT, f"search_bag_{n}.json")
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        t = np.array([h["t"] for h in d["hist"]])
        EL = np.array([h["EL_disk"] for h in d["hist"]])
        ax.semilogy(t / d["tc"], np.maximum(EL / EL[0], 1e-5), label=lab, lw=1.2)
    ax.axhline(0.9, color="k", ls=":", lw=0.8)
    ax.set_xlabel("tempo / tempo di attraversamento")
    ax.set_ylabel("energia della LUCE nel disco / iniziale")
    ax.set_title("bolle pre-formate (caso piu' favorevole)", fontsize=9)
    ax.legend(fontsize=6)
    ax = axs[2]
    for kind, mk in (("D20", "o"), ("D400", "s"), ("U0.6", "^"), ("M8", "v")):
        xs, ys = [], []
        for fn in os.listdir(OUT):
            if fn.startswith(f"search_bag_col_{kind}_e") and fn.endswith(".json"):
                d = json.load(open(os.path.join(OUT, fn)))
                xs.append(d["cfg"]["e0"])
                ys.append(max(bag_row(d)["luce_nel_disco_a_t_form"], 1e-9))
        if xs:
            o = np.argsort(xs)
            ax.loglog(np.array(xs)[o], np.array(ys)[o], mk + "-", label=kind)
    ax.set_xlabel("campo dei pacchetti e0")
    ax.set_ylabel("energia di luce nel disco / iniziale (t_uscita + 2 t_c)")
    ax.set_title("urto gamma-gamma nel vuoto normale", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "search_history_bag.png"), dpi=80)
    print("figura search_history_bag.png")


def babycheck(kappa=1.0, mu2=0.25):
    """controllo del reticolo: baby-skyrmione radiale (1D, molto fine) per kappa=1, V = 0.25 (1-n3)."""
    import torch
    torch.set_default_dtype(torch.float64)
    torch.set_num_threads(1)
    R, n = 20.0, 4000
    r = torch.linspace(0, R, n + 1)
    dr = float(r[1])
    z = (np.pi * torch.exp(-r[1:-1] / 1.5)).clone().requires_grad_(True)

    def E(z):
        f = torch.cat([torch.tensor([np.pi]), z, torch.zeros(1)])
        fm, rm = 0.5 * (f[1:] + f[:-1]), 0.5 * (r[1:] + r[:-1])
        fp = (f[1:] - f[:-1]) / dr
        e = 0.5 * (fp ** 2 + torch.sin(fm) ** 2 / rm ** 2) + 0.5 * kappa ** 2 * torch.sin(fm) ** 2 * fp ** 2 / rm ** 2             + mu2 * (1 - torch.cos(fm))
        return 2 * np.pi * (e * rm).sum() * dr
    opt = torch.optim.LBFGS([z], lr=1, max_iter=5000, tolerance_grad=1e-12, tolerance_change=1e-15,
                            line_search_fn="strong_wolfe")

    def cl():
        opt.zero_grad()
        e = E(z)
        e.backward()
        return e
    opt.step(cl)
    val = float(E(z))
    print(f"baby-skyrmione radiale kappa={kappa}, mu^2={mu2}: E = {val:.4f}")
    return val


if __name__ == "__main__":
    if sys.argv[1] == "babycheck":
        babycheck()
        sys.exit()
    {"finalists": finalists, "predict": predict, "figs": figs, "objfig": fig_objects, "bag": bag}[sys.argv[1]]()

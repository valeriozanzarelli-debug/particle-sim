"""ciclo_p1_map.py - Phase 1 of CICLO_LIVELLI: measure the map between levels on the levels that ALREADY exist
(read-only: lab/out/agg_zz_1d.npz, agg_zz_3d_{A,B}.pkl, aggfull_proton.pkl, aggfull_relearnT.pkl, aggfull_l1_T.pkl,
aggmol_levels.pkl, aggfull_learn.json, aggmol_learn.json). Exams P1a-P1e of lab/out/ciclo_prereg.json (fase1).
Writes lab/out/ciclo_p1_map.json.   Run: python lab/ciclo_p1_map.py   (~1 min)"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import agg_zz_core as ac
import zz_core as zz
ac.WK = 1

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
TH_A = A0 / 16 / 256
TH_B = A0 / 16 / 32
KW = np.pi / 2                       # window of the normalized-shape distance


def D_cont(h, kw=KW):
    k = np.linspace(1e-4, kw, 2001)
    return float(np.max(np.abs(2 * h * h * (np.sqrt(1 + k * k / (h * h)) - 1) - k * k)) / kw ** 2)


# ------------------------------------------------------------------ 3-D generators
def scalar_stencil(K):
    return 0.5 * (K[0, 0] + K[1, 1])


def moments3(g, offs):
    j = offs.astype(float)
    c2 = np.array([-0.5 * np.sum(g * j[:, a] ** 2) for a in range(3)])
    c4 = np.array([np.sum(g * j[:, a] ** 4) / 24 for a in range(3)])
    c6 = np.array([-np.sum(g * j[:, a] ** 6) / 720 for a in range(3)])
    cxy = np.array([np.sum(g * j[:, a] ** 2 * j[:, b] ** 2) / 4 for a, b in ((0, 1), (1, 2), (0, 2))])
    return dict(c2=float(np.real(c2.mean())), c4=float(np.real(c4.mean())), c6=float(np.real(c6.mean())),
                cxxyy=float(np.real(cxy.mean())), c2_aniso=float(np.real(c2.max() - c2.min()) / abs(np.real(c2.mean()))),
                imag_c2=float(np.abs(np.imag(c2)).max()))


def G3(g, offs, kv):
    return np.real(np.exp(1j * kv @ offs.T.astype(float)) @ g)


def shape3(g, offs, c2):
    k = np.linspace(1e-3, KW, 400)
    out = {}
    for name, d in (("axis", np.array([1.0, 0, 0])), ("diag", np.ones(3) / np.sqrt(3))):
        kv = k[:, None] * d[None, :]
        S = (G3(g, offs, kv) - G3(g, offs, np.zeros((1, 3)))[0]) / c2
        out[name] = S
    return k, out


def level3(rec, theta, chain, B, kind, mass=1.0):
    g = scalar_stencil(rec["K"])
    offs = rec["offs"]
    h = B * theta                              # cell in Compton lengths of THIS particle (theta already m a)
    mo = moments3(g, offs)
    k, S = shape3(g, offs, mo["c2"])
    Dn = max(float(np.max(np.abs(S[n] - k * k))) / KW ** 2 for n in S)
    spin = float(np.abs(rec["K"][0, 1]).max() + np.abs(rec["K"][0, 0] - rec["K"][1, 1]).max()) / float(np.abs(g).max())
    return dict(chain=chain, kind=kind, B=B, h=h, w=int(np.abs(offs).max()), taps=int(len(offs)),
                m_kin=1 / (2 * h * mo["c2"]), c4_over_c2=mo["c4"] / mo["c2"], c4_ratio_dirac=-(mo["c4"] / mo["c2"]) * 4 * h * h,
                c6_over_c2=mo["c6"] / mo["c2"], c6_ratio_dirac=(mo["c6"] / mo["c2"]) * 8 * h ** 4,
                iso=mo["cxxyy"] / mo["c4"] if mo["c4"] != 0 else float("nan"), c2_aniso=mo["c2_aniso"], spin_rel=spin,
                D=Dn, D_cont=D_cont(h), _k=k, _S=S, _g=g, _offs=offs, _c2=mo["c2"])


# ------------------------------------------------------------------ 1-D kernels
def level1(K, offs, B, chain):
    th = TH_A
    h = B * th
    r = ac.Rule(K, offs, 1, 2, True)
    kap = np.linspace(0, np.pi, 8001)
    M = r.mult_at(kap[:, None])
    lam, V = np.linalg.eig(M)
    wphi = np.abs(V[:, 0, :] + V[:, 1, :]) ** 2 / 2
    ip = np.argmax(wphi, axis=1)
    lp = lam[np.arange(len(kap)), ip]
    ph = np.unwrap(np.angle(lp))
    G = -(ph - ph[0])
    kf = min(0.3 * h, 0.8)
    s = (kap > 0) & (kap <= kf)
    A = np.stack([kap[s] ** 2, kap[s] ** 4, kap[s] ** 6], 1)
    c = np.linalg.lstsq(A, G[s], rcond=None)[0]
    c2, c4, c6 = c
    w = kap <= KW
    S = G[w] / c2
    Dn = float(np.max(np.abs(S - kap[w] ** 2))) / KW ** 2
    # exact positive branch of the fine 1-D walk seen at this cell (cos w = cos th cos k), per coarse step
    Gf = (np.arccos(np.cos(th) * np.cos(kap / B)) - th) * B
    Dlearn = float(np.max(np.abs(G[w] - Gf[w]))) / (c2 * KW ** 2)
    # generator stencil of the positive band (Fourier series of G over the BZ, even part)
    kk = 2 * np.pi * np.arange(512) / 512
    kk = np.where(kk > np.pi, kk - 2 * np.pi, kk)
    Gb = np.interp(np.abs(kk), kap, G)
    gst = np.real(np.fft.fft(Gb)) / 512
    return dict(chain=chain, kind="free", B=B, h=h, w=int(np.abs(offs).max()), taps=int(len(offs)),
                m_kin=1 / (2 * h * c2), c4_over_c2=c4 / c2, c4_ratio_dirac=-(c4 / c2) * 4 * h * h,
                c6_ratio_dirac=(c6 / c2) * 8 * h ** 4, D=Dn, D_cont=D_cont(h), D_learn_vs_fine=Dlearn,
                gen_stencil_norm=(gst[:13] / c2).tolist(), _k=kap[w], _S={"axis": S}, _c2=c2)


def slope(recs):
    h = np.array([r["h"] for r in recs])
    y = np.abs(np.array([r["c4_over_c2"] for r in recs]))
    p = np.polyfit(np.log(h), np.log(y), 1)
    return float(p[0]), float(p[1])


def main():
    res = dict(chains={}, exams={})
    chains = {}
    # 1-D
    z = np.load(os.path.join(OUT, "agg_zz_1d.npz"))
    j1 = json.load(open(os.path.join(OUT, "agg_zz_1d.json")))
    for key in ("ladder", "auto"):
        chains["1D_" + key] = [level1(z[f"{key}_B{r['B']}_K"], z[f"{key}_B{r['B']}_offs"], r["B"], "1D_" + key) for r in j1[key]]
    # 3-D chain A (free), chain B (zoom), proton
    for name, fn, th in (("3D_A", "agg_zz_3d_A.pkl", TH_A), ("3D_B", "agg_zz_3d_B.pkl", TH_B)):
        st = pickle.load(open(os.path.join(OUT, fn), "rb"))
        chains[name] = [level3(lv, th, name, lv["B"], "free") for lv in st["levels"]]
    st = pickle.load(open(os.path.join(OUT, "aggfull_proton.pkl"), "rb"))
    pj = json.load(open(os.path.join(OUT, "aggfull_proton.json")))
    thp = pj["theta_p"]
    chains["proton"] = [level3(lv, thp, "proton", lv["B"], "free") for lv in st["levels"]]
    # constrained (T) levels, reported
    T = []
    l1 = pickle.load(open(os.path.join(OUT, "aggfull_l1_T.pkl"), "rb"))
    T.append(level3(l1, TH_A, "T", l1["B"], "T"))
    for lv in pickle.load(open(os.path.join(OUT, "aggfull_relearnT.pkl"), "rb")):
        T.append(level3(lv, TH_A, "T", lv["B"], "T"))
    for lv in pickle.load(open(os.path.join(OUT, "aggmol_levels.pkl"), "rb"))["free"]:
        T.append(level3(lv, TH_A, "T", lv["B"], "T(aggmol)"))
    chains["T"] = T

    tant = {"1D_ladder": np.tan(TH_A) / TH_A, "1D_auto": np.tan(TH_A) / TH_A, "3D_A": None, "3D_B": None, "proton": None}
    # 3-D fine rule's own kinetic mass: from its Bloch phase at small k (same as aggfull: tan th / th to 5e-13)
    tant["3D_A"] = np.tan(TH_A) / TH_A
    tant["3D_B"] = np.tan(TH_B) / TH_B
    tant["proton"] = np.tan(thp) / thp
    free = ["1D_ladder", "1D_auto", "3D_A", "proton"]
    # P1a
    worst = 0.0
    for c in free:
        for r in chains[c]:
            r["m_kin_over_fine"] = r["m_kin"] / tant[c]
            worst = max(worst, abs(r["m_kin_over_fine"] - 1))
    res["exams"]["P1a"] = dict(value=worst, tol=1e-3, ok=worst <= 1e-3)
    # P1b
    sl = {}
    for c in free + ["3D_B"]:
        if len(chains[c]) >= 2:
            sl[c] = slope(chains[c])[0]
    res["slopes_c4"] = sl
    res["exams"]["P1b"] = dict(value={c: sl[c] for c in free}, tol=0.1, ok=all(abs(sl[c] + 2) <= 0.1 for c in free))
    # P1c
    iso = [r["iso"] for r in chains["3D_A"]]
    res["exams"]["P1c"] = dict(value=iso, lo=1.8, hi=2.2, ok=all(1.8 <= x <= 2.2 for x in iso))
    # P1d
    lj = json.load(open(os.path.join(OUT, "aggfull_learn.json")))
    mj = json.load(open(os.path.join(OUT, "aggmol_learn.json")))
    phys = [(lv["B"], lv["physics"]) for lv in lj["levels"] if lv.get("physics")]
    phys.append((1024, mj["B1024_from512"]["coupling"]["physics"]))
    coup = {}
    okd = True
    for key, ref in (("gk", -0.125), ("gx", 0.25), ("gso", 0.25)):
        vals = np.array([p[key][0] for _, p in phys])
        hB = np.array([B * TH_A for B, _ in phys])
        spread = float((vals.max() - vals.min()) / abs(vals.mean()))
        # weight relative to Coulomb at the cell scale: |coef| (kappa/h)^2 at kappa = 1 -> exponent
        yfit = float(np.polyfit(np.log(hB), np.log(np.abs(vals) / hB ** 2), 1)[0])
        coup[key] = dict(B=[B for B, _ in phys], values=vals.tolist(), dirac=ref, spread=spread, y=yfit)
        okd = okd and spread <= 0.02 and abs(yfit + 2) <= 0.03
    res["couplings"] = coup
    res["exams"]["P1d"] = dict(value={k: (v["spread"], v["y"]) for k, v in coup.items()}, tol=0.02, ok=okd)
    # P1e
    oke = True
    det = {}
    for c in ("1D_ladder", "1D_auto"):
        D = [r["D"] for r in chains[c]]
        mono = all(D[i + 1] < D[i] for i in range(len(D) - 1))
        rel = [(r["B"], r["D"] / r["D_cont"] - 1) for r in chains[c] if r["h"] >= 2]
        okc = mono and all(abs(x) <= 0.25 for _, x in rel)
        det[c] = dict(D=D, monotone=mono, rel_to_cont_h_ge_2=rel, ok=okc)
        oke = oke and okc
    res["exams"]["P1e"] = dict(value=det, tol=0.25, ok=oke)
    # successive distances (reported)
    for c, L in chains.items():
        for i in range(len(L) - 1):
            a, b = L[i], L[i + 1]
            ka = a["_k"]; Sa = a["_S"]["axis"]; Sb = np.interp(ka, b["_k"], b["_S"]["axis"])
            a["d_next"] = float(np.max(np.abs(Sa - Sb))) / KW ** 2
    # analytic fixed-point stencil (1-D, separable in 3-D): g*_j = 2(-1)^j / j^2, g*_0 = pi^2/3
    res["fixed_point_stencil"] = [np.pi ** 2 / 3] + [2 * (-1) ** j / j ** 2 for j in range(1, 13)]
    res["chains"] = {c: [{k: v for k, v in r.items() if not k.startswith("_")} for r in L] for c, L in chains.items()}
    # arrays for the figure
    arr = {}
    for c, L in chains.items():
        for r in L:
            arr[f"{c}_B{r['B']}_k"] = r["_k"]
            for n, S in r["_S"].items():
                arr[f"{c}_B{r['B']}_S_{n}"] = S
    np.savez_compressed(os.path.join(OUT, "ciclo_p1_shapes.npz"), **arr)
    json.dump(res, open(os.path.join(OUT, "ciclo_p1_map.json"), "w"), indent=1, default=float)
    for c, L in chains.items():
        print(f"== {c}")
        for r in L:
            print(f"  B={r['B']:5d} h={r['h']:9.4f} w={r['w']} m_kin={r['m_kin']:.6f} c4/c2*(-4h^2)={r['c4_ratio_dirac']:+.4f} "
                  f"iso={r.get('iso', float('nan')):+.3f} D={r['D']:.3e} Dcont={r['D_cont']:.3e} d_next={r.get('d_next', float('nan')):.2e}")
    print("slopes", sl)
    for k, v in res["exams"].items():
        print(f"[{'PASS' if v['ok'] else 'FAIL'}] {k}: {v['value']}")
    print("couplings", json.dumps(coup))


if __name__ == "__main__":
    main()

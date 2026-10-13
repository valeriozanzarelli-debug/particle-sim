"""multi_e_disp.py - M0: livello 1 RI-IMPARATO con l'osservabile Z = 10 (1s di taglia a0/10), stesso algoritmo e stessi
iper-parametri di AGGLOMERATO par. 9 (solo l'osservabile cambia, budget energia 3e-5), a theta_A e theta_A/2.
Poi: banda simmetrizzata mediata sulle direzioni eps_bar_L(k), tabella di anisotropia, cucitura con la banda fine oltre
la zona (|k| > pi/h), F(s) = E_{p~N(0,sI3)}[eps_bar_L(|p|)] e F4(s) = 15 s^2 (tabelle per l'ECG), anche per il livello
1 esistente (catena A). Uscite: lab/out/multi_e_disp.json (riassunto) e lab/out/multi_e_Ftab.npz (tabelle, non committate).
Run: python lab/multi_e_disp.py learn ; python lab/multi_e_disp.py tables"""
import os, sys, time, pickle
from multi_e_core import zz, OUT, A0, ALPHA, THETA, Band, chain_levels, load_json, save_json, fib_dirs, dirs26
import numpy as np
import agg_zz_core as agg

JS = "multi_e_disp.json"
PKL = os.path.join(OUT, "multi_e_relearn.pkl")
NPZ = os.path.join(OUT, "multi_e_Ftab.npz")
res = load_json(JS)
sec = sys.argv[1] if len(sys.argv) > 1 else "learn"
ZOBS = 10


def fine_band_eh(k_au, th):
    """banda positiva della regola fine (media di spin), energia di moto in Eh, k in a.u. lungo le direzioni date (..., 3)"""
    kap = np.asarray(k_au) * th / A0                     # 1/cella fine (cella fine = th lambda_C = th/A0 a0)
    sh = kap.shape[:-1]
    W = zz.bloch_batch(kap.reshape(-1, 3), th)
    T = agg.conv_k(W, "dirac4", "pauli2")
    ph = agg.band_phase(T, "pauli2")                    # fase in due tick
    W0 = zz.bloch_batch(np.zeros((1, 3)), th)
    ph0 = agg.band_phase(agg.conv_k(W0, "dirac4", "pauli2"), "pauli2")[0]
    om = -np.angle(np.exp(1j * (ph - ph0))) / (2 * th)  # mc^2
    return (om / ALPHA ** 2).reshape(sh)


if sec == "learn":
    st = pickle.load(open(PKL, "rb")) if os.path.exists(PKL) else {}
    for chain in ("A", "A2"):
        if chain in st:
            continue
        TH = THETA[chain]
        kfac = ZOBS / A0
        rng = np.random.default_rng(20261013)
        prev = agg.Level("fine3", 3, "dirac4", 4, 1, TH, period=2)
        t1 = time.time()
        logl = []
        lv, tried = agg.auto_level(prev, rng, 20, eps_packet=1e-3, eps_energy=3e-5, T_obs=1000, k_obs=kfac * TH,
                                   bs=(16, 8, 4, 3, 2), ws=(1, 2, 3, 4), log=lambda s: (print(s, flush=True), logl.append(s)),
                                   ks_rule=lambda B: -kfac * TH * B, learner="auto", packets="phys", energy_mode="state")
        if lv is None:
            st[chain] = None
            res.setdefault("relearn", {})[chain] = dict(chosen=None, tried=tried, log=logl)
        else:
            st[chain] = dict(K=lv.rule.K, offs=lv.rule.offs, B=lv.B, b=lv.b, w=lv.w, phi0=getattr(lv.rule, "phi0", None),
                             desc=lv.desc, C=lv.C, check=lv.check)
            # diagnostica: il miglior candidato Pauli (b = 16, w = 4), RIFIUTATO dal controllo, imparato di nuovo uguale
            rng2 = np.random.default_rng(20261014)
            rp = agg.learn_h(prev, 16, "pauli2", 4, -kfac * TH * 16, rng2)
            cp = agg.Level("learned", 3, "pauli2", 2, 16, TH, rule=rp, b=16, w=4)
            ckp = agg.kcheck(prev, cp, 16, rng2, 1000, kfac * TH, sig=max(1.0, 2.0 / (kfac * TH * 16)),
                             k0_cells=0.5 * kfac * TH * 16, energy_mode="state")
            st[chain + "_pauli_rejected"] = dict(K=rp.K, offs=rp.offs, B=16, b=16, w=4, phi0=rp.phi0, desc="pauli2", C=2,
                                                 check=ckp)
            res.setdefault("relearn", {})[chain] = dict(chosen=dict(B=lv.B, b=lv.b, w=lv.w, desc=lv.desc, check=lv.check,
                                                                     cell_a0=lv.B * TH / A0), n_tried=len(tried),
                                                         best_pauli_rejected_check=ckp,
                                                         sec=time.time() - t1, log=logl)
        pickle.dump(st, open(PKL, "wb"))
        save_json(JS, res)
        print(chain, "fatto", flush=True)

class D4Band:
    """banda positiva (vista dal descrittore pauli2, media di spin, simmetrizzata +-k) di un livello a 4 componenti"""

    def __init__(self, L, TH):
        self.rule = agg.Rule(L["K"], L["offs"], 3, 4, unitary=True)
        self.h_lc = L["B"] * TH
        self.h = self.h_lc / A0
        z = np.zeros((1, 3))
        self.ph0 = agg.band_phase(agg.conv_k(self.rule.mult_at(z), "dirac4", "pauli2"), "pauli2")[0]

    def _e(self, kap):
        ph = agg.band_phase(agg.conv_k(self.rule.mult_at(kap), "dirac4", "pauli2"), "pauli2")
        return -np.angle(np.exp(1j * (ph - self.ph0))) / self.h_lc * A0 ** 2

    def eps(self, k_au):
        k = np.asarray(k_au, float)
        sh = k.shape[:-1]
        kap = k.reshape(-1, 3) * self.h
        return (0.5 * (self._e(kap) + self._e(-kap))).reshape(sh)


def make_band(L, TH):
    if L.get("desc", "pauli2") == "pauli2":
        return Band(L["K"], L["offs"], L["B"], TH, sym=True)
    return D4Band(L, TH)


if sec == "tables":
    st = pickle.load(open(PKL, "rb"))
    bands = {"existing_A": (chain_levels("A")[0], THETA["A"])}
    for ch in ("A", "A2"):
        if st.get(ch):
            bands[f"relearn_{ch}"] = (st[ch], THETA[ch])
        if st.get(ch + "_pauli_rejected"):
            bands[f"relearn_{ch}_pauli_rejected"] = (st[ch + "_pauli_rejected"], THETA[ch])
    D = fib_dirs(400)
    tabs = {}
    summ = {}
    for name, (L, TH) in bands.items():
        t1 = time.time()
        b = make_band(L, TH)
        kz = np.pi / b.h                                   # bordo della zona lungo gli assi (a.u.)
        kf = np.pi * A0 / TH * 0.999                       # bordo della zona fine
        # curvatura (m_kin del livello)
        kap = np.linspace(1e-3, 4e-3, 7)
        rows, ys = [], []
        for d in dirs26():
            k_au = kap[:, None] * d[None, :] / b.h
            k2 = np.sum(k_au ** 2, 1); rows.append(np.stack([k2, k2 ** 2], 1)); ys.append(b.eps(k_au))
        c = np.linalg.lstsq(np.concatenate(rows), np.concatenate(ys), rcond=None)[0]
        m_kin = 1 / (2 * c[0])
        # tabella radiale dentro la zona (media sferica)
        kin = np.unique(np.concatenate([[0.0], np.geomspace(1e-3, min(200.0, kz), 2001), np.geomspace(min(200.0, kz), kz, 600)]))
        eb = np.concatenate([[0.0]] + [np.mean(b.eps(kin[i:i + 20, None, None] * D[None]), axis=1)
                                        for i in range(1, len(kin), 20)])
        # anisotropia: asse / faccia / corpo a k = 1, 10, 30, 100 a.u.
        an = {}
        for kk in (1.0, 10.0, 30.0, 100.0):
            for nm, d in (("axis", np.array([1., 0, 0])), ("face", np.array([1., 1, 0]) / 2 ** .5), ("body", np.ones(3) / 3 ** .5)):
                an[f"{nm}_k{kk:g}"] = float(b.eps(kk * d))
            an[f"mean_k{kk:g}"] = float(np.mean(b.eps(kk * D)))          # media diretta (non interpolata)
            # contenuto k^4: (eps - k^2/2m_kin) / (-alpha^2 k^4/8)
            for nm in ("axis", "face", "body", "mean"):
                an[f"k4frac_{nm}_k{kk:g}"] = float((an[f"{nm}_k{kk:g}"] - kk ** 2 / (2 * m_kin)) / (-ALPHA ** 2 * kk ** 4 / 8))
        # oltre la zona: banda fine (media sferica), cucita
        kout = np.geomspace(kz, kf, 1500)[1:]
        Dsm = fib_dirs(100)
        efine = np.array([np.mean(fine_band_eh(k * Dsm, TH)) for k in kout])
        efine_at_z = float(np.mean(fine_band_eh(kz * Dsm, TH)))
        jump = float(eb[-1] - efine_at_z)
        ktab = np.concatenate([kin, kout]); etab = np.concatenate([eb, efine])
        # F(s): media su p ~ N(0, s I3) di eps_bar(|p|) = int chi_3(u) eps(u sqrt s) du
        u = np.linspace(0, 14, 7001)[1:]
        wu = np.sqrt(2 / np.pi) * u ** 2 * np.exp(-u ** 2 / 2) * (u[1] - u[0])
        svals = np.geomspace(1e-6, 1e8, 1401)
        F = np.empty_like(svals); Fb = np.empty_like(svals)
        # interpolazione del RAPPORTO rho(k) = eps/(k^2/2) (liscio, -> 1/m_kin a k -> 0): interpolare eps stesso
        # linearmente sovrastima k^2 di ~D k^2/6 (errore del 5% a s = 0.01 nella prima versione)
        rtab = np.concatenate([[1 / m_kin], etab[1:] / (ktab[1:] ** 2 / 2)])
        # sotto k = 1 a.u. la differenza G(k) - G(0) e' ~1e-9 rad e perde ~1e-5 relativo per cancellazione: li' si usa
        # la serie della curvatura (stesso fit di m_kin), eps = c2 k^2 + c4 k^4
        lo = ktab < 1.0
        rtab[lo] = 2 * c[0] + 2 * c[1] * ktab[lo] ** 2
        for i, s in enumerate(svals):
            p = u * np.sqrt(s)
            e = np.interp(p, ktab, rtab) * p ** 2 / 2
            F[i] = np.sum(wu * e)
            Fb[i] = np.sum(wu * np.where(p > kz, e, 0.0))
        R = F / (1.5 * svals)                              # -> 1/m_kin per s -> 0
        tabs[name] = dict(s=svals, R=R, Fb=Fb, ktab=ktab, etab=etab, m_kin=m_kin, kz=kz)
        summ[name] = dict(B=int(L["B"]), w=int(L.get("w", 0)), desc=L.get("desc", "pauli2"), check=L.get("check"), theta=TH, cell_a0=b.h, k_zone_au=kz, m_kin=m_kin,
                          c4_over_dirac=float(c[1] / (-ALPHA ** 2 / 8)), anisotropy=an, splice_jump_Eh=jump,
                          splice_jump_rel=jump / efine_at_z, eps_zone_Eh=float(eb[-1]),
                          R_at_s={f"{s:g}": float(np.interp(np.log(s), np.log(svals), R)) for s in (1e-2, 1, 1e2, 1e4, 1e5, 1e6)})
        print(name, "m_kin", m_kin, "k4 frac mean k10", an["k4frac_mean_k10"], "axis/face/body", f"{time.time() - t1:.0f}s",
              an["k4frac_axis_k10"], an["k4frac_face_k10"], an["k4frac_body_k10"], "jump", jump, flush=True)
    np.savez(NPZ, **{f"{n}__{k}": v for n, t in tabs.items() for k, v in t.items()})
    res["bands"] = summ
    save_json(JS, res)
    print("tabelle scritte")

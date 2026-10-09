"""search_bag.py - il meccanismo "a sacco" di Valerio: luce chiusa da un vuoto esterno
che spinge (docs/RICERCA_REGOLE.md, parte A).

Tre pezzi:
  A1  mit      : modello a sacco MIT per fotoni e per gluoni/quark (analitico + radici di Bessel),
                 termine di Casimir (Boyer), predizioni e look-elsewhere.
  A2  run CFG  : simulazione 2D (modo TE, psi = A_z) di luce + campo di "struttura" sigma con
                 doppia buca (sigma=0 vuoto normale, sigma=1 fase "sacco"), tre accoppiamenti:
                   D  dielettrico Lorentz-invariante: L = 1/2 eps(sigma) (E^2 - B^2)
                   U  accoppiamento alla densita' di energia u (NON Lorentz-invariante):
                      L = 1/2 (1+g s) E^2 - 1/2 (1-g s) B^2   (=> + g u s(sigma) nella lagrangiana)
                   M  massa del fotone nella fase sigma=1:  - 1/2 m0^2 s(sigma) A^2
                 Configurazioni: urto di due pacchetti contrapposti nel vuoto normale (sigma=0),
                 oppure bolla "pre-formata" (caso piu' favorevole), oppure bolle senza luce (Coleman).
  A3  table    : raccoglie i risultati.

Unita': c = 1; lambda_luce = 1 per gli urti.  Uso:
  python lab/search_bag.py mit
  timeout 590 python lab/search_bag.py run <nome>      (vedi CONFIGS)
  python lab/search_bag.py table
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq
from scipy.special import spherical_jn

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

HBARC = 197.3269804  # MeV fm
M_E, M_P, M_MU, M_TAU, M_DELTA, M_N = 0.51099895, 938.27208816, 105.6583755, 1776.86, 1232.0, 938.92
TARGETS = {"mu/e": 206.768283, "tau/e": 3477.23, "p/e": 1836.15267, "Delta/N": 1.3121}


# ============================================================================ A1: MIT bag
def bessel_roots():
    f_te = lambda x: spherical_jn(1, x) + x * spherical_jn(1, x, derivative=True)  # (x j1)'=0
    f_tm = lambda x: spherical_jn(1, x)
    f_q = lambda x: spherical_jn(0, x) - spherical_jn(1, x)  # quark senza massa, MIT
    x_te = brentq(f_te, 2.0, 3.5)
    x_tm = brentq(f_tm, 4.0, 5.0)
    x_q = brentq(f_q, 1.5, 2.5)
    # modi l=2 (per la torre)
    f_te2 = lambda x: spherical_jn(2, x) + x * spherical_jn(2, x, derivative=True)
    x_te2 = brentq(f_te2, 3.0, 4.5)
    return x_te, x_tm, x_q, x_te2


def bag_min(S, B):
    """E(R) = S/R + 4pi/3 B R^3  -> R* = (S/4piB)^(1/4), E* = (4/3) S / R*.  (unita' hbar=c=1)"""
    R = (S / (4 * np.pi * B)) ** 0.25
    return R, 4.0 / 3.0 * S / R


def mit():
    x_te, x_tm, x_q, x_te2 = bessel_roots()
    res = {"x_TE(1+,gluone/fotone)": x_te, "x_TM(1-)": x_tm, "x_quark(massless)": x_q, "x_TE(l=2)": x_te2}
    Z_casimir = -0.04618  # Boyer 1968: E_C = +0.04618 hbar c / a  (repulsivo) => Z0 = -0.046
    res["Z0_Casimir_Boyer"] = Z_casimir
    lines = []
    # --- adroni (controllo del metodo): B fissato dal nucleone
    for Z0, lab in [(0.0, "Z0=0"), (1.84, "Z0=1.84 (DeGrand et al. 1975)")]:
        S_N = 3 * x_q - Z0
        # E* = (4/3) S^(3/4) (4 pi B)^(1/4)  => B
        B14 = (M_N / HBARC / (4 / 3 * S_N ** 0.75 * (4 * np.pi) ** 0.25)) ** 1.0  # fm^-1 units (B^(1/4))
        B = B14 ** 4
        R_N, E_N = bag_min(S_N, B)
        S_G = 2 * x_te - Z0
        R_G, E_G = bag_min(S_G, B)
        S_D = S_N  # senza termine colore-magnetico N e Delta sono degeneri
        lines.append({"caso": lab, "B^1/4 [MeV]": B14 * HBARC, "R_N [fm]": R_N, "M_N [MeV]": E_N * HBARC,
                      "glueball 2 gluoni TE [MeV]": E_G * HBARC, "R_glueball [fm]": R_G,
                      "Delta/N senza termine colore-magnetico": 1.0})
    res["adroni"] = lines
    # --- fotoni: un fotone nel modo TE; scala fissata dall'elettrone o dal protone
    phot = []
    for mass, name in [(M_E, "elettrone"), (M_P, "protone")]:
        for Z0, lab in [(0.0, "senza Casimir"), (Z_casimir, "con Casimir di Boyer")]:
            S1 = x_te - Z0
            B14 = mass / HBARC / (4 / 3 * S1 ** 0.75 * (4 * np.pi) ** 0.25)
            R1, E1 = bag_min(S1, B14 ** 4)
            phot.append({"scala": name, "Casimir": lab, "B^1/4 [MeV]": B14 * HBARC, "R* [fm]": R1,
                         "R* [m]": R1 * 1e-15, "E* [MeV]": E1 * HBARC})
    res["fotoni"] = phot
    # --- predizioni (rapporti) di un sacco di fotoni: E ~ S^(3/4), S = sum n_i x_i
    ratios = {"2 fotoni TE / 1": 2 ** 0.75, "3 fotoni TE / 1": 3 ** 0.75,
              "1 fotone TM / 1 TE": (x_tm / x_te) ** 0.75, "1 fotone TE(l=2) / TE(l=1)": (x_te2 / x_te) ** 0.75,
              "TE+TM / TE": ((x_te + x_tm) / x_te) ** 0.75}
    res["rapporti_predetti"] = ratios
    comp = []
    for tn, tv in TARGETS.items():
        best = min(ratios.items(), key=lambda kv: abs(np.log(kv[1] / tv)))
        comp.append({"bersaglio": tn, "valore": tv, "predizione piu' vicina": best[0], "rapporto": best[1],
                     "scarto %": 100 * (best[1] / tv - 1)})
    res["confronto"] = comp
    # --- look-elsewhere: se il numero N di fotoni e' libero, E_N/E_1 = N^(3/4) copre tutto
    lee = []
    for tn, tv in TARGETS.items():
        N = tv ** (4 / 3)
        Nr = max(1, round(N))
        lee.append({"bersaglio": tn, "N necessario": N, "N intero": Nr, "N^(3/4)": Nr ** 0.75,
                    "scarto %": 100 * (Nr ** 0.75 / tv - 1),
                    "spaziatura relativa vicina (%)": 100 * 0.75 / Nr})
    res["look_elsewhere_N_libero"] = lee
    # probabilita' di "colpire" a caso un bersaglio entro 5% con 5 rapporti tra 1 e 4000 (log-uniforme)
    p_one = 2 * np.log(1.05) / np.log(4000)
    res["P(caso) un rapporto entro 5% di un dato bersaglio (log-uniforme 1..4000)"] = p_one
    res["P(caso) almeno una coincidenza 5 rapporti x 4 bersagli"] = 1 - (1 - p_one) ** 20
    with open(os.path.join(OUT, "search_bag_mit.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


# ============================================================================ A2: 2D dynamics
def smooth(s):
    return 3 * s ** 2 - 2 * s ** 3


def dsmooth(s):
    return 6 * s * (1 - s)


class Sim:
    def __init__(s, cfg):
        s.cfg = cfg
        s.N = cfg.get("N", 192)
        s.dx = cfg.get("dx", 1 / 8)
        s.L = s.N * s.dx
        x = (np.arange(s.N) - s.N // 2) * s.dx
        s.X, s.Y = np.meshgrid(x, x, indexing="ij")
        s.R = np.hypot(s.X, s.Y)
        W = cfg.get("sponge", 3.0)
        d = np.maximum(0, np.maximum(np.abs(s.X), np.abs(s.Y)) - (s.L / 2 - W)) / W
        s.gam = cfg.get("gamma0", 4.0) * d ** 2
        s.inner = d == 0
        s.kind = cfg["kind"]
        s.g = cfg.get("g", 0.0)
        s.eps1 = cfg.get("eps1", 1.0)
        s.m0 = cfg.get("m0", 0.0)
        s.ms = cfg.get("msig", 2.0)
        s.dU = cfg.get("dU", 0.05)

    # coefficienti a(sigma), b(sigma), m2(sigma) e derivate
    def coeffs(s, sig):
        sm, ds = smooth(np.clip(sig, -0.5, 1.5)), dsmooth(np.clip(sig, -0.5, 1.5))
        if s.kind == "D":
            a = 1 + (s.eps1 - 1) * sm
            da = (s.eps1 - 1) * ds
            return a, da, a, da, 0.0, 0.0
        if s.kind == "U":
            return 1 + s.g * sm, s.g * ds, 1 - s.g * sm, -s.g * ds, 0.0, 0.0
        if s.kind == "M":
            one = np.ones_like(sig)
            return one, 0.0, one, 0.0, s.m0 ** 2 * sm, s.m0 ** 2 * ds
        if s.kind == "N":  # nessun accoppiamento (controllo / Coleman)
            one = np.ones_like(sig)
            return one, 0.0, one, 0.0, 0.0, 0.0
        raise ValueError

    def U(s, sig):
        return 0.5 * s.ms ** 2 * sig ** 2 * (1 - sig) ** 2 + s.dU * smooth(sig)

    def dU_(s, sig):
        return s.ms ** 2 * sig * (1 - sig) * (1 - 2 * sig) + s.dU * dsmooth(sig)

    def rhs(s, st):
        psi, p, sig, pi = st
        a, da, b, db, m2, dm2 = s.coeffs(sig)
        h = s.dx
        gx = (np.roll(psi, -1, 0) - psi) / h
        gy = (np.roll(psi, -1, 1) - psi) / h
        bx = 0.5 * (b + np.roll(b, -1, 0)) if np.ndim(b) else b
        by = 0.5 * (b + np.roll(b, -1, 1)) if np.ndim(b) else b
        fx, fy = bx * gx, by * gy
        div = (fx - np.roll(fx, 1, 0)) / h + (fy - np.roll(fy, 1, 1)) / h
        psit = p / a
        dp = div - m2 * psi - s.gam * p
        # forza su sigma = -dH/dsigma
        lap = (np.roll(sig, 1, 0) + np.roll(sig, -1, 0) + np.roll(sig, 1, 1) + np.roll(sig, -1, 1) - 4 * sig) / h ** 2
        F = lap - s.dU_(sig) + 0.5 * da * psit ** 2 - dm2 * 0.5 * psi ** 2
        if np.ndim(db):
            g2x, g2y = gx ** 2, gy ** 2
            F = F - 0.25 * db * (g2x + np.roll(g2x, 1, 0) + g2y + np.roll(g2y, 1, 1))
        dpi = F - s.gam * pi
        return np.array([psit, dp, pi, dpi])

    def energy_density(s, st):
        psi, p, sig, pi = st
        a, _, b, _, m2, _ = s.coeffs(sig)
        h = s.dx
        gx = (np.roll(psi, -1, 0) - psi) / h
        gy = (np.roll(psi, -1, 1) - psi) / h
        eL = 0.5 * p ** 2 / a + 0.5 * b * 0.5 * (gx ** 2 + np.roll(gx, 1, 0) ** 2 + gy ** 2 + np.roll(gy, 1, 1) ** 2) \
            + 0.5 * m2 * psi ** 2
        sx = (np.roll(sig, -1, 0) - sig) / h
        sy = (np.roll(sig, -1, 1) - sig) / h
        eS = 0.5 * pi ** 2 + 0.5 * (sx ** 2 + sy ** 2) + s.U(sig)
        return eL, eS

    def init(s):
        c = s.cfg
        rng = np.random.default_rng(c.get("seed", 1))
        psi = np.zeros((s.N, s.N)); p = np.zeros_like(psi); sig = np.zeros_like(psi); pi = np.zeros_like(psi)
        mode = c["init"]
        if mode == "collide":
            e0, k = c["e0"], 2 * np.pi
            for sgn in (+1, -1):
                x0 = -sgn * c.get("x0", 6.0)
                env = np.exp(-((s.X - x0) ** 2) / (2 * 1.0 ** 2) - s.Y ** 2 / (2 * 2.0 ** 2))
                ph = k * (s.X - x0)
                psi += (e0 / k) * env * np.cos(ph)
                # onda f(x - sgn t): psi_t = -sgn dpsi/dx
                dpsidx = (e0 / k) * env * (-k * np.sin(ph) - (s.X - x0) / 1.0 ** 2 * np.cos(ph))
                p += -sgn * dpsidx  # a = 1 nel vuoto normale
            sig += c.get("noise", 1e-3) * rng.standard_normal(sig.shape)
        elif mode in ("bubble", "shell", "bag_outside", "empty_bubble"):
            R0, w = c.get("R0", 3.0), 0.35
            if mode in ("bubble", "empty_bubble"):
                sig = 0.5 * (1 - np.tanh((s.R - R0) / w))
            elif mode == "shell":
                d = c.get("shell", 1.0)
                sig = 0.5 * (np.tanh((s.R - R0) / w) - np.tanh((s.R - R0 - d) / w))
            else:
                sig = 0.5 * (1 + np.tanh((s.R - R0) / w))
            if mode != "empty_bubble":
                a, _, _, _, _, _ = s.coeffs(sig)
                om, ell, rm = c["omega"], c["ell"], c.get("rm", 0.75 * R0)
                th = np.arctan2(s.Y, s.X)
                g = np.exp(-((s.R - rm) ** 2) / (2 * c.get("wr", 0.5) ** 2))
                psi = c["amp"] * g * np.cos(ell * th)
                p = a * c["amp"] * g * om * np.sin(ell * th)  # modo che ruota: cos(l th - om t)
        st = np.array([psi, p, sig, pi])
        if c.get("perturb", 0) > 0:
            for i in range(4):
                fld = st[i]
                sc = np.abs(fld).max()
                st[i] = fld * (1 + c["perturb"] * rng.standard_normal(fld.shape)) if sc > 0 else fld
        return st

    def run(s, tmax, dt, rdisk, every=1.0, t_leave=None, tc=None, verbose=True):
        st = s.init()
        disk = s.R < rdisk
        eL, eS = s.energy_density(st)
        E0 = (eL + eS)[s.inner].sum() * s.dx ** 2
        hist = []
        t = 0.0
        nsteps = int(round(tmax / dt))
        nev = max(1, int(round(every / dt)))
        t0 = time.time()
        dead = 0
        for n in range(nsteps + 1):
            if n % nev == 0:
                eL, eS = s.energy_density(st)
                e = eL + eS
                rec = dict(t=round(t, 3), E_in=float(e[s.inner].sum() * s.dx ** 2),
                           E_disk=float(e[disk].sum() * s.dx ** 2), EL_disk=float(eL[disk].sum() * s.dx ** 2),
                           sig_max=float(st[2][s.inner].max()), sig_disk_mean=float(st[2][disk].mean()))
                hist.append(rec)
                if not np.isfinite(rec["E_in"]):
                    print("NaN, stop"); break
                if verbose and n % (nev * 20) == 0:
                    print(f"t={t:7.2f} E_in={rec['E_in']:.4f} E_disk={rec['E_disk']:.4f} EL_disk={rec['EL_disk']:.4f}"
                          f" sigmax={rec['sig_max']:.3f}  [{time.time()-t0:.0f}s]", flush=True)
                # arresto anticipato: niente nel disco dopo che i pacchetti sono usciti
                if t_leave is not None and t > t_leave + 3 * tc:
                    dead = dead + 1 if rec["E_disk"] < 0.01 * E0 else 0
                    if dead * every > 3 * tc:
                        print("early stop: disco vuoto"); break
            k1 = s.rhs(st)
            k2 = s.rhs(st + 0.5 * dt * k1)
            k3 = s.rhs(st + 0.5 * dt * k2)
            k4 = s.rhs(st + dt * k3)
            st = st + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
        return E0, hist, st


BASE = dict(N=160, dx=1 / 6, msig=2.0, dU=0.05)
CONFIGS = {}
# urti nel vuoto normale (sigma = 0 ovunque, luce libera) -- il test chiave
for kind, extra, amps in [("D", dict(eps1=20.0), (0.1, 0.5, 1.5)), ("D", dict(eps1=400.0), (0.1, 0.5, 1.5)),
                          ("U", dict(g=0.6), (0.5, 1.5, 4.0)), ("M", dict(m0=8.0), (0.5, 1.5, 4.0))]:
    for e0 in amps:
        tag = f"col_{kind}{list(extra.values())[0]:g}_e{e0:g}"
        CONFIGS[tag] = dict(BASE, kind=kind, init="collide", e0=e0, **extra)
        if extra.get("eps1", 0) > 100:
            CONFIGS[tag]["dt"] = 0.015  # forza su sigma ~ eps1 (E^2-B^2): passo piu' piccolo
# bolle pre-formate (caso piu' favorevole), luce dentro.  omega = pi (lambda = 2) per risolvere bene l'interno
CONFIGS["seed_D400_bubble"] = dict(BASE, kind="D", eps1=400.0, init="bubble", R0=3.0, omega=np.pi, ell=4,
                                   amp=0.3 / np.sqrt(20.0), rm=1.8, wr=0.6, dt=0.015)
CONFIGS["seed_D20_bubble"] = dict(BASE, kind="D", eps1=20.0, init="bubble", R0=3.0, omega=np.pi, ell=4, amp=0.3,
                                  rm=1.8, wr=0.6)
CONFIGS["seed_U0.6_wgm"] = dict(BASE, kind="U", g=0.6, init="bubble", R0=3.0, omega=np.pi, ell=15, amp=0.3, rm=2.4,
                                wr=0.4)
CONFIGS["seed_M8_shell"] = dict(BASE, kind="M", m0=8.0, init="shell", R0=3.0, shell=1.5, omega=np.pi, ell=4,
                                amp=0.3, rm=1.8, wr=0.6)
CONFIGS["seed_M8_bagout_dUneg"] = dict(BASE, kind="M", m0=8.0, dU=-0.05, init="bag_outside", R0=3.0, omega=np.pi,
                                       ell=4, amp=0.3, rm=1.8, wr=0.6)
# stesse bolle con piu' luce dentro (la pressione della luce deve bilanciare tensione + Delta U)
for amp in (0.45, 0.6):
    CONFIGS[f"seed_M8_bagout_dUneg_amp{amp:g}"] = dict(CONFIGS["seed_M8_bagout_dUneg"], amp=amp)
    CONFIGS[f"seed_U0.6_wgm_amp{amp:g}"] = dict(CONFIGS["seed_U0.6_wgm"], amp=amp)
# Coleman: bolle senza luce
CONFIGS["coleman_false_R3"] = dict(BASE, kind="N", dU=0.05, init="empty_bubble", R0=3.0)
CONFIGS["coleman_true_R4"] = dict(BASE, kind="N", dU=-0.05, init="empty_bubble", R0=4.0)
CONFIGS["coleman_true_R8"] = dict(BASE, kind="N", dU=-0.05, init="empty_bubble", R0=8.0)


def run_cfg(name, perturb=0.0):
    cfg = dict(CONFIGS[name])
    if perturb:
        cfg["perturb"] = perturb
        cfg["seed"] = 7
        name = name + f"_pert{perturb:g}"
    sim = Sim(cfg)
    dt = cfg.get("dt", 0.06)
    if cfg["init"] == "collide":
        rdisk, tc = 3.0, 6.0
        t_leave = 2 * cfg.get("x0", 6.0) + 6.0   # i pacchetti escono dal disco
        tmax = t_leave + 100 * tc
        E0, hist, st = sim.run(tmax, dt, rdisk, every=0.5, t_leave=t_leave, tc=tc)
    else:
        rdisk = cfg.get("R0", 3.0) + 2.0
        tc = 2 * cfg.get("R0", 3.0)
        t_leave = None
        tmax = 100 * tc
        E0, hist, st = sim.run(tmax, dt, rdisk, every=0.5)
    # sintesi
    T = np.array([h["t"] for h in hist]); Ed = np.array([h["E_disk"] for h in hist])
    summ = dict(name=name, cfg={k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in cfg.items()},
                E0=E0, t_end=float(T[-1]), rdisk=rdisk, tc=tc)
    if cfg["init"] == "collide":
        tf = t_leave + 2 * tc
        i0 = np.searchsorted(T, tf)
        Ef = Ed[i0] if i0 < len(Ed) else 0.0
        summ["yield_E_disk_over_E0_at_t_form"] = float(Ef / E0)
        later = Ed[i0:]
        summ["C2_keep90_for_100tc"] = bool(Ef > 0.01 * E0 and T[-1] >= tf + 100 * tc - 1e-6 and later.min() >= 0.9 * Ef)
        summ["sig_max_ever"] = float(max(h["sig_max"] for h in hist))
    else:
        summ["E_disk_frac"] = {str(tt): float(Ed[np.searchsorted(T, tt)] / Ed[0]) for tt in
                               [tc, 5 * tc, 10 * tc, 25 * tc, 50 * tc, 100 * tc] if tt <= T[-1]}
        summ["sig_disk_mean_end"] = hist[-1]["sig_disk_mean"]
        summ["sig_disk_mean_start"] = hist[0]["sig_disk_mean"]
        # raggio della bolla: dove sigma attraversa 1/2 lungo l'asse x
        sig = st[2]
        mid = sig[sim.N // 2:, sim.N // 2]
        summ["sigma_profile_x_end"] = [float(v) for v in mid[::4]]
        summ["C2_keep90_for_100tc"] = bool(T[-1] >= 100 * tc - 1e-6 and Ed.min() >= 0.9 * Ed[0])
    summ["hist"] = hist[::4]
    np.save(os.path.join(OUT, f"search_bag_{name}_final.npy"), st[[0, 2]].astype(np.float32))
    with open(os.path.join(OUT, f"search_bag_{name}.json"), "w") as f:
        json.dump(summ, f)
    print(json.dumps({k: v for k, v in summ.items() if k != "hist"}, indent=1))


def table():
    rows = []
    for fn in sorted(os.listdir(OUT)):
        if fn.startswith("search_bag_") and fn.endswith(".json") and fn != "search_bag_mit.json":
            d = json.load(open(os.path.join(OUT, fn)))
            r = {k: v for k, v in d.items() if k not in ("hist", "cfg", "sigma_profile_x_end")}
            rows.append(r)
    for r in rows:
        print(json.dumps(r))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "mit":
        mit()
    elif cmd == "run":
        run_cfg(sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 0.0)
    elif cmd == "list":
        print("\n".join(CONFIGS))
    elif cmd == "table":
        table()

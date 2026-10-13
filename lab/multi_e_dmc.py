"""multi_e_dmc.py - camminatori in tempo immaginario per l'elio (DV6: Python + numba; il bin Rust 'elettroni' viene dopo).
Regola: il passo libero e' la rotazione di Wick del livello 1 (R2: nucleo positivo per beta >= 5e-5; limite continuo =
diffusione gaussiana con D = 1/(2 m_kin) letto da lab/out/multi_e_kernel.json, DEBT-9); Coulomb istantaneo (DEBT-1).
Algoritmo (brief 3.7): passi elettrone per elettrone con drift UNR limitato, accetta/rifiuta di Reynolds, un passo che
attraversa il nodo della guida e' rifiutato, branching con la media di E_L prima e dopo, T_pop = 1 a.u., errori a blocchi.
Guida (DV7): orbitali del campo medio (multi_e_meanfield.py carrier) x exp(u(r12)), u = c r/(1 + b r), c = cuspide di Kato
derivata (m_kin/2 spin opposti, m_kin/4 stesso spin), b = 1 fisso; c = 0 per D6.
Run: python lab/multi_e_dmc.py <stato> <tau> [modo] ;  stato in {H, He, He3S, He_u0}, modo in {dmc, vmc}
Checkpoint: lab/out/multi_e_dmc_<stato>_<modo>_tau<tau>.npz (non committato); risultati in lab/out/multi_e_dmc.json"""
import os, sys, time
from multi_e_core import OUT, load_json, save_json
import numpy as np
import numba as nb
from scipy.interpolate import CubicSpline

JS = "multi_e_dmc.json"


def orbital_tables(state):
    """cubic spline of R(x) = P(r)/r, x = ln r, as piecewise-polynomial coefficients (for numba)"""
    if state in ("H", "Hb"):
        m = load_json("multi_e_kernel.json")["R3"]["m_kin_theta_A_used"]
        x = np.linspace(np.log(1e-9), np.log(60), 6000); r = np.exp(x)
        a = m if state == "H" else 0.8 * m              # Hb: guida volutamente sbagliata (test della macchina DMC)
        P = np.array([r * np.exp(-a * r)])
        return x, P, m
    nm = {"He": "He", "He_u0": "He", "He3S": "He3S"}[state]
    z = np.load(os.path.join(OUT, f"multi_e_orb_{nm}.npz"))
    return z["x"], z["P"], float(z["m_kin"])


def spline_coeffs(x, P, r):
    out = []
    for Pi in P:
        R = Pi / r
        cs = CubicSpline(x, R)
        out.append(cs.c)                       # (4, n-1)
    return np.array(out)


@nb.njit(cache=True)
def orb_eval(C, x0, hx, nseg, k, rr):
    """value, d/dr, laplacian of the radial orbital k at distance rr"""
    if rr < 1e-12:
        rr = 1e-12
    x = np.log(rr)
    i = int((x - x0) / hx)
    if i < 0:
        i = 0
    if i >= nseg:
        i = nseg - 1
    t = x - (x0 + i * hx)
    c0, c1, c2, c3 = C[k, 0, i], C[k, 1, i], C[k, 2, i], C[k, 3, i]
    R = ((c0 * t + c1) * t + c2) * t + c3
    Rx = (3 * c0 * t + 2 * c1) * t + c2
    Rxx = 6 * c0 * t + 2 * c1
    dR = Rx / rr
    lap = (Rxx + Rx) / (rr * rr)
    return R, dR, lap


@nb.njit(cache=True)
def guide(pos, C, x0, hx, nseg, kind, Z, cu, bu, D):
    """pos (ne, 3). kind 0: 1 electron (orbital 0); 1: singlet phi0(1)phi0(2) e^u; 2: triplet det[phi0, phi1] e^u.
    returns psi (signed), grad ln|psi| (ne, 3), local energy"""
    ne = pos.shape[0]
    g = np.zeros((ne, 3))
    V = 0.0
    r = np.zeros(ne)
    for i in range(ne):
        r[i] = np.sqrt(pos[i, 0] ** 2 + pos[i, 1] ** 2 + pos[i, 2] ** 2)
        V -= Z / r[i]
    if kind == 0:
        R, dR, lap = orb_eval(C, x0, hx, nseg, 0, r[0])
        for d in range(3):
            g[0, d] = dR / R * pos[0, d] / r[0]
        return R, g, -D * lap / R + V
    # orbitals of each electron
    va = np.zeros((2, 2)); da = np.zeros((2, 2)); la = np.zeros((2, 2))      # [orbital, electron]
    for k in range(2 if kind == 2 else 1):
        for i in range(2):
            a, b_, c = orb_eval(C, x0, hx, nseg, k, r[i])
            va[k, i] = a; da[k, i] = b_; la[k, i] = c
    # Phi and its derivatives (radial orbitals: grad = d/dr * rhat)
    gradP = np.zeros((2, 3)); lapP = np.zeros(2)
    if kind == 1:
        Phi = va[0, 0] * va[0, 1]
        for d in range(3):
            gradP[0, d] = da[0, 0] * pos[0, d] / r[0] * va[0, 1]
            gradP[1, d] = va[0, 0] * da[0, 1] * pos[1, d] / r[1]
        lapP[0] = la[0, 0] * va[0, 1]; lapP[1] = va[0, 0] * la[0, 1]
    else:
        Phi = va[0, 0] * va[1, 1] - va[1, 0] * va[0, 1]
        for d in range(3):
            gradP[0, d] = (da[0, 0] * va[1, 1] - da[1, 0] * va[0, 1]) * pos[0, d] / r[0]
            gradP[1, d] = (va[0, 0] * da[1, 1] - va[1, 0] * da[0, 1]) * pos[1, d] / r[1]
        lapP[0] = la[0, 0] * va[1, 1] - la[1, 0] * va[0, 1]
        lapP[1] = va[0, 0] * la[1, 1] - va[1, 0] * la[0, 1]
    # Jastrow u(r12)
    r12v = np.zeros(3)
    for d in range(3):
        r12v[d] = pos[0, d] - pos[1, d]
    r12 = np.sqrt(r12v[0] ** 2 + r12v[1] ** 2 + r12v[2] ** 2)
    V += 1.0 / r12
    den = 1 + bu * r12
    u = cu * r12 / den
    up = cu / (den * den)
    upp = -2 * cu * bu / (den * den * den)
    lapU = upp + 2 * up / r12
    psi = Phi * np.exp(u)
    lapsum = 0.0
    for i in range(2):
        sgn = 1.0 if i == 0 else -1.0
        dot = 0.0; gu2 = 0.0
        for d in range(3):
            gU = sgn * up * r12v[d] / r12
            gp = gradP[i, d] / Phi
            g[i, d] = gp + gU
            dot += gp * gU; gu2 += gU * gU
        lapsum += lapP[i] / Phi + 2 * dot + lapU + gu2
    EL = -D * lapsum + V
    return psi, g, EL


@nb.njit(cache=True)
def run_block(P, Wt, nsteps, tau, D, C, x0, hx, nseg, kind, Z, cu, bu, Eref, NT, Tpop, Ebest, dmc, seed, hist, hbin):
    """P (nw, ne, 3) walkers; returns new walkers, sums. hist: intracule histogram (opposite/any spin pair) weighted"""
    np.random.seed(seed)
    ne = P.shape[1]
    te = 2 * D * tau                          # variance per component of the free move
    sumE = 0.0; sumW = 0.0; acc = 0.0; prop = 0.0; accr2 = 0.0; propr2 = 0.0
    Etrace = np.zeros(nsteps); Wtrace = np.zeros(nsteps); Ntrace = np.zeros(nsteps)
    cnum = 0.0; cden = 0.0
    for st in range(nsteps):
        nw = P.shape[0]
        newW = np.zeros(nw)
        Eacc = 0.0; Wacc = 0.0
        for w in range(nw):
            x = P[w].copy()
            psi, g, EL0 = guide(x, C, x0, hx, nseg, kind, Z, cu, bu, D)
            for i in range(ne):
                # UNR-capped drift for electron i
                v2 = g[i, 0] ** 2 + g[i, 1] ** 2 + g[i, 2] ** 2
                f = 1.0
                if v2 * te > 1e-12:
                    f = (-1 + np.sqrt(1 + 2 * v2 * te)) / (v2 * te)
                xn = x.copy()
                dr2 = 0.0
                for d in range(3):
                    xn[i, d] = x[i, d] + te * f * g[i, d] + np.sqrt(te) * np.random.normal()
                    dr2 += (xn[i, d] - x[i, d]) ** 2
                psin, gn, ELn = guide(xn, C, x0, hx, nseg, kind, Z, cu, bu, D)
                prop += 1; propr2 += dr2
                if psin * psi <= 0:
                    continue                              # nodo attraversato: rifiutato
                v2n = gn[i, 0] ** 2 + gn[i, 1] ** 2 + gn[i, 2] ** 2
                fn = 1.0
                if v2n * te > 1e-12:
                    fn = (-1 + np.sqrt(1 + 2 * v2n * te)) / (v2n * te)
                lf = 0.0; lb = 0.0
                for d in range(3):
                    lf += (xn[i, d] - x[i, d] - te * f * g[i, d]) ** 2
                    lb += (x[i, d] - xn[i, d] - te * fn * gn[i, d]) ** 2
                ratio = (psin / psi) ** 2 * np.exp(-(lb - lf) / (2 * te))
                if np.random.random() < ratio:
                    x = xn; psi = psin; g = gn
                    acc += 1; accr2 += dr2
            psi, g, EL1 = guide(x, C, x0, hx, nseg, kind, Z, cu, bu, D)
            P[w] = x
            if dmc:
                tef = tau * (accr2 / propr2 if propr2 > 0 else 1.0)
                ww = Wt[w] * np.exp(-tef * (0.5 * (EL0 + EL1) - Eref))
            else:
                ww = 1.0
            newW[w] = ww
            Eacc += ww * EL1; Wacc += ww
            if hbin > 0 and ne == 2:
                r12 = np.sqrt((x[0, 0] - x[1, 0]) ** 2 + (x[0, 1] - x[1, 1]) ** 2 + (x[0, 2] - x[1, 2]) ** 2)
                b = int(r12 / hbin)
                if b < hist.shape[0]:
                    hist[b] += ww
        Etrace[st] = Eacc / Wacc; Wtrace[st] = Wacc / nw; Ntrace[st] = nw
        sumE += Eacc; sumW += Wacc
        cnum = 0.995 * cnum + Eacc / nw; cden = 0.995 * cden + Wacc / nw
        Ebest = cnum / cden
        if dmc:
            # branching: copie intere floor(w + u)
            m = 0
            cnt = np.zeros(nw, np.int64)
            for w in range(nw):
                c = int(newW[w] + np.random.random())
                cnt[w] = c; m += c
            if m < 1:
                m = 1; cnt[0] = 1
            Pn = np.empty((m, ne, 3))
            k = 0
            for w in range(nw):
                for c in range(cnt[w]):
                    Pn[k] = P[w]; k += 1
            P = Pn
            Wt = np.ones(m)
            Eref = Ebest - np.log(m / NT) / Tpop
    return P, Wt, Etrace, Wtrace, Ntrace, Eref, Ebest, acc / max(prop, 1.0), accr2 / max(propr2, 1e-300)


def blocking(E, W, nb_list=(32, 16, 8)):
    errs = []
    for nb_ in nb_list:
        L = len(E) // nb_
        Eb = np.array([np.sum(E[i * L:(i + 1) * L] * W[i * L:(i + 1) * L]) / np.sum(W[i * L:(i + 1) * L]) for i in range(nb_)])
        errs.append(Eb.std(ddof=1) / np.sqrt(nb_))
    return float(np.sum(E * W) / np.sum(W)), float(max(errs))


if __name__ == "__main__":
    state = sys.argv[1]; tau = float(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else "dmc"
    T_total = float(sys.argv[4]) if len(sys.argv) > 4 else 400.0     # a.u. di produzione (tetto)
    sig_goal = float(sys.argv[5]) if len(sys.argv) > 5 else 6e-5
    budget = float(os.environ.get("DMC_BUDGET", "540"))
    kern = load_json("multi_e_kernel.json")
    D = kern["R2"]["D"]                                         # D dalla curvatura del livello 1 (non 1/2 scritto a mano)
    x, P, m_orb = orbital_tables(state)
    r = np.exp(x)
    C = spline_coeffs(x, P, r)
    x0, hx, nseg = float(x[0]), float(x[1] - x[0]), C.shape[2]
    kind = {"Hb": 0, "H": 0, "He": 1, "He_u0": 1, "He3S": 2}[state]
    Z = {"Hb": 1.0, "H": 1.0, "He": 2.0, "He_u0": 2.0, "He3S": 2.0}[state]
    m_kin = 1 / (2 * D)
    cu = {"Hb": 0.0, "H": 0.0, "He": m_kin / 2, "He_u0": 0.0, "He3S": m_kin / 4}[state]
    bu = 1.0
    ne = 1 if kind == 0 else 2
    NT = 4000
    key = f"{state}_{mode}_tau{tau:g}"
    ck = os.path.join(OUT, f"multi_e_dmc_{key}.npz")
    hbin = 0.01 if state == "He_u0" else 0.0
    nbins = 200
    if os.path.exists(ck):
        z = np.load(ck)
        Pw, Et, Wt_, Nt = z["P"], list(z["Et"]), list(z["Wt"]), list(z["Nt"]); Eref, Ebest = float(z["Eref"]), float(z["Ebest"])
        hist = z["hist"]; nstep_done = int(z["nstep"]); accs = list(z["accs"])
    else:
        rng = np.random.default_rng(11)
        Pw = rng.normal(size=(NT, ne, 3)) * (1.0 / Z)
        Et, Wt_, Nt, accs = [], [], [], []
        Eref = Ebest = {"Hb": -0.5, "H": -0.5, "He": -2.9, "He_u0": -2.9, "He3S": -2.17}[state] * m_kin
        hist = np.zeros(nbins); nstep_done = 0
    neq = int(round(20.0 / tau))                               # 20 a.u. di equilibratura scartati
    t0 = time.time()
    Wts = np.ones(len(Pw))
    while time.time() - t0 < budget:
        nchunk = int(round(2.0 / tau))
        hb = np.zeros(nbins)
        Pw, Wts, E_, W_, N_, Eref, Ebest, a, ar = run_block(Pw, Wts, nchunk, tau, D, C, x0, hx, nseg, kind, Z, cu, bu, Eref,
                                                            NT, 1.0, Ebest, mode == "dmc", 1000 + nstep_done // nchunk, hb,
                                                            hbin)
        if nstep_done >= neq:
            hist += hb
        Et += list(E_); Wt_ += list(W_); Nt += list(N_); accs.append(a)
        nstep_done += nchunk
        np.savez(ck, P=Pw, Et=Et, Wt=Wt_, Nt=Nt, Eref=Eref, Ebest=Ebest, hist=hist, nstep=nstep_done, accs=accs)
        Ea, Wa = np.array(Et[neq:]), np.array(Wt_[neq:]) if mode == "dmc" else np.ones(len(Et[neq:]))
        if len(Ea) > 32 * 50:
            E, s = blocking(Ea, Wa)
            if s <= sig_goal or (nstep_done - neq) * tau >= T_total:
                break
    Ea = np.array(Et[neq:]); Wa = np.array(Wt_[neq:]) if mode == "dmc" else np.ones(len(Ea))
    res = load_json(JS)
    if len(Ea) > 32 * 10:
        E, s = blocking(Ea, Wa)
        res[key] = dict(E=E, sigma=s, prod_time_au=len(Ea) * tau, walkers_mean=float(np.mean(Nt[neq:])), acceptance=float(np.mean(accs)),
                        tau=tau, D=D, m_kin=m_kin, cusp_c=cu, b=bu, mode=mode, done=bool(s <= sig_goal or len(Ea) * tau >= T_total))
        if hbin:
            res[key]["intracule_hist"] = hist.tolist(); res[key]["hbin"] = hbin
        save_json(JS, res)
        print(key, res[key] if not hbin else {k: v for k, v in res[key].items() if k != "intracule_hist"}, flush=True)
    else:
        print(key, "equilibratura", nstep_done, flush=True)

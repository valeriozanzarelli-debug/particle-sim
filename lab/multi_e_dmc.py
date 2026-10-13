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


@nb.njit(cache=True, fastmath=False)
def guide(p, C, x0, hx, nseg, kind, Z, cu, bu, D):
    """p: 6 numbers (x1,y1,z1,x2,y2,z2). kind 0: 1 electron (orbital 0, electron 2 ignored); 1: singlet phi0(1) phi0(2) e^u;
    2: triplet det[phi0, phi1] e^u.  Returns psi (signed), g1x..g2z (grad ln|psi|), local energy. No allocations."""
    x1, y1, z1, x2, y2, z2 = p[0], p[1], p[2], p[3], p[4], p[5]
    r1 = np.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
    if kind == 0:
        R, dR, lap = orb_eval(C, x0, hx, nseg, 0, r1)
        f = dR / R / r1
        return R, f * x1, f * y1, f * z1, 0.0, 0.0, 0.0, -D * lap / R - Z / r1
    r2 = np.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    a1, d1, l1 = orb_eval(C, x0, hx, nseg, 0, r1)
    a2, d2, l2 = orb_eval(C, x0, hx, nseg, 0, r2)
    if kind == 1:
        Phi = a1 * a2
        gp1 = d1 / a1 / r1; gp2 = d2 / a2 / r2           # grad Phi / Phi = gp * r_vec
        lp1 = l1 / a1; lp2 = l2 / a2
        G1x, G1y, G1z = gp1 * x1, gp1 * y1, gp1 * z1
        G2x, G2y, G2z = gp2 * x2, gp2 * y2, gp2 * z2
    else:
        b1, e1, m1 = orb_eval(C, x0, hx, nseg, 1, r1)
        b2, e2, m2 = orb_eval(C, x0, hx, nseg, 1, r2)
        Phi = a1 * b2 - b1 * a2
        f1 = (d1 * b2 - e1 * a2) / Phi / r1
        f2 = (a1 * e2 - b1 * d2) / Phi / r2
        lp1 = (l1 * b2 - m1 * a2) / Phi
        lp2 = (a1 * m2 - b1 * l2) / Phi
        G1x, G1y, G1z = f1 * x1, f1 * y1, f1 * z1
        G2x, G2y, G2z = f2 * x2, f2 * y2, f2 * z2
    dx, dy, dz = x1 - x2, y1 - y2, z1 - z2
    r12 = np.sqrt(dx * dx + dy * dy + dz * dz)
    den = 1 + bu * r12
    u = cu * r12 / den
    up = cu / (den * den)
    upp = -2 * cu * bu / (den * den * den)
    lapU = upp + 2 * up / r12
    q = up / r12
    U1x, U1y, U1z = q * dx, q * dy, q * dz                  # grad_1 u; grad_2 u = -grad_1 u
    gu2 = up * up
    lap = lp1 + lp2 + 2 * ((G1x * U1x + G1y * U1y + G1z * U1z) - (G2x * U1x + G2y * U1y + G2z * U1z)) + 2 * (lapU + gu2)
    EL = -D * lap - Z / r1 - Z / r2 + 1.0 / r12
    return (Phi * np.exp(u), G1x + U1x, G1y + U1y, G1z + U1z, G2x - U1x, G2y - U1y, G2z - U1z, EL)


@nb.njit(cache=True)
def walker_step(xw, tau, te, D, C, x0, hx, nseg, kind, Z, cu, bu, ne):
    """moves one walker electron by electron (in place); returns EL before, EL after, accepted dr2, proposed dr2, r12"""
    o = guide(xw, C, x0, hx, nseg, kind, Z, cu, bu, D)
    psi = o[0]; EL0 = o[7]
    g = np.empty(6)
    for k in range(6):
        g[k] = o[1 + k]
    xn = np.empty(6)
    acc2 = 0.0; prop2 = 0.0
    for i in range(ne):
        v2 = g[3 * i] ** 2 + g[3 * i + 1] ** 2 + g[3 * i + 2] ** 2
        f = 1.0
        if v2 * te > 1e-12:
            f = (-1 + np.sqrt(1 + 2 * v2 * te)) / (v2 * te)
        for k in range(6):
            xn[k] = xw[k]
        dr2 = 0.0
        lf = 0.0
        for d in range(3):
            eta = np.random.normal()
            xn[3 * i + d] = xw[3 * i + d] + te * f * g[3 * i + d] + np.sqrt(te) * eta
            dr2 += (xn[3 * i + d] - xw[3 * i + d]) ** 2
            lf += te * eta * eta
        on = guide(xn, C, x0, hx, nseg, kind, Z, cu, bu, D)
        prop2 += dr2
        if on[0] * psi <= 0:
            continue                                     # nodo attraversato: rifiutato
        v2n = on[1 + 3 * i] ** 2 + on[2 + 3 * i] ** 2 + on[3 + 3 * i] ** 2
        fn = 1.0
        if v2n * te > 1e-12:
            fn = (-1 + np.sqrt(1 + 2 * v2n * te)) / (v2n * te)
        lb = 0.0
        for d in range(3):
            lb += (xw[3 * i + d] - xn[3 * i + d] - te * fn * on[1 + 3 * i + d]) ** 2
        ratio = (on[0] / psi) ** 2 * np.exp(-(lb - lf) / (2 * te))
        if np.random.random() < ratio:
            for d in range(3):
                xw[3 * i + d] = xn[3 * i + d]
            psi = on[0]
            for k in range(6):
                g[k] = on[1 + k]
            acc2 += dr2
    o = guide(xw, C, x0, hx, nseg, kind, Z, cu, bu, D)
    r12 = np.sqrt((xw[0] - xw[3]) ** 2 + (xw[1] - xw[4]) ** 2 + (xw[2] - xw[5]) ** 2)
    return EL0, o[7], acc2, prop2, r12


@nb.njit(cache=True, parallel=True)
def sweep(P, tau, te, D, C, x0, hx, nseg, kind, Z, cu, bu, ne):
    nw = P.shape[0]
    E0 = np.empty(nw); E1 = np.empty(nw); A2 = np.empty(nw); P2 = np.empty(nw); R12 = np.empty(nw)
    for w in nb.prange(nw):
        E0[w], E1[w], A2[w], P2[w], R12[w] = walker_step(P[w], tau, te, D, C, x0, hx, nseg, kind, Z, cu, bu, ne)
    return E0, E1, A2, P2, R12


@nb.njit(cache=True)
def seed_all(s):
    np.random.seed(s)


def run_block(P, nsteps, tau, D, C, x0, hx, nseg, kind, Z, cu, bu, Eref, NT, Tpop, Ebest, dmc, rng, hist, hbin, state):
    """P (nw, 6). Branching in numpy (integer copies floor(w + u)); mixed estimator recorded at every step."""
    ne = 1 if kind == 0 else 2
    te = 2 * D * tau
    Et = np.zeros(nsteps); Wt = np.zeros(nsteps); Nt = np.zeros(nsteps)
    for st in range(nsteps):
        E0, E1, A2, P2, R12 = sweep(P, tau, te, D, C, x0, hx, nseg, kind, Z, cu, bu, ne)
        state["acc2"] += A2.sum(); state["prop2"] += P2.sum()
        if dmc:
            tef = tau * state["acc2"] / state["prop2"]
            w = np.exp(-tef * (0.5 * (E0 + E1) - Eref))
        else:
            w = np.ones(len(E1))
        Et[st] = np.sum(w * E1) / np.sum(w); Wt[st] = np.mean(w); Nt[st] = len(P)
        if hbin > 0:
            b = (R12 / hbin).astype(int); ok = b < len(hist)
            np.add.at(hist, b[ok], w[ok])
        state["cn"] = 0.995 * state["cn"] + np.mean(w * E1); state["cd"] = 0.995 * state["cd"] + np.mean(w)
        Ebest = state["cn"] / state["cd"]
        if dmc:
            c = (w + rng.random(len(w))).astype(int)
            if c.sum() < 1:
                c[0] = 1
            P = np.repeat(P, c, axis=0)
            Eref = Ebest - np.log(len(P) / NT) / Tpop
    return P, Et, Wt, Nt, Eref, Ebest


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
    S0 = int(os.environ.get("DMC_SEED", "0"))                  # 0 = semi della v1; D6' (v2) usa 1
    NT = int(os.environ.get("DMC_NW", "4000"))                 # 4000 = pre-registrato; altro solo per la diagnosi del bias di popolazione
    key = f"{state}_{mode}_tau{tau:g}" + ("" if NT == 4000 else f"_NW{NT}") + os.environ.get("DMC_TAG", "")   # TAG: corse diagnostiche post hoc
    ck = os.path.join(OUT, f"multi_e_dmc_{key}.npz")
    hbin = 0.01 if state == "He_u0" else 0.0
    nbins = 200
    if os.path.exists(ck):
        z = np.load(ck, allow_pickle=True)
        Pw, Et, Wt_, Nt = z["P"], list(z["Et"]), list(z["Wt"]), list(z["Nt"]); Eref, Ebest = float(z["Eref"]), float(z["Ebest"])
        hist = z["hist"]; nstep_done = int(z["nstep"]); stt = z["state"].item()
        HB = [list(r) for r in z["HB"]] if "HB" in z.files else []          # istogrammi per pezzo (errori di D6)
    else:
        rng0 = np.random.default_rng(11 + S0)
        Pw = rng0.normal(size=(NT, 6)) * (1.0 / Z)
        Et, Wt_, Nt = [], [], []
        Eref = Ebest = {"Hb": -0.5, "H": -0.5, "He": -2.9, "He_u0": -2.9, "He3S": -2.17}[state] * m_kin
        hist = np.zeros(nbins); nstep_done = 0; HB = []
        stt = dict(acc2=0.0, prop2=0.0, cn=0.0, cd=0.0)
    neq = int(round(20.0 / tau))                               # 20 a.u. di equilibratura scartati
    t0 = time.time()
    rng = np.random.default_rng(1000 + nstep_done + 7919 * S0)
    seed_all(1000 + nstep_done + 7919 * S0)
    while time.time() - t0 < budget:
        nchunk = int(round(2.0 / tau))
        hb = np.zeros(nbins)
        Pw, E_, W_, N_, Eref, Ebest = run_block(Pw, nchunk, tau, D, C, x0, hx, nseg, kind, Z, cu, bu, Eref, NT, 1.0, Ebest,
                                                mode == "dmc", rng, hb, hbin, stt)
        if nstep_done >= neq:
            hist += hb
            if hbin:
                HB.append(list(hb) + [float(np.sum(W_ * N_))])          # ultimo elemento: peso totale del pezzo
        Et += list(E_); Wt_ += list(W_); Nt += list(N_)
        nstep_done += nchunk
        np.savez(ck, P=Pw, Et=Et, Wt=Wt_, Nt=Nt, Eref=Eref, Ebest=Ebest, hist=hist, nstep=nstep_done, state=np.array(stt, dtype=object),
                 HB=np.array(HB if HB else np.zeros((0, nbins + 1))))
        Ea = np.array(Et[neq:]); Wa = np.array(Wt_[neq:]) if mode == "dmc" else np.ones(len(Et[neq:]))
        if len(Ea) > 32 * 50:
            E, s = blocking(Ea, Wa)
            if s <= sig_goal or (nstep_done - neq) * tau >= T_total:
                break
    Ea = np.array(Et[neq:]); Wa = np.array(Wt_[neq:]) if mode == "dmc" else np.ones(len(Ea))
    res = load_json(JS)
    if len(Ea) > 32 * 10:
        E, s = blocking(Ea, Wa)
        res[key] = dict(E=E, sigma=s, prod_time_au=len(Ea) * tau, walkers_mean=float(np.mean(Nt[neq:])),
                        acc_ratio_r2=float(stt["acc2"] / stt["prop2"]),
                        tau=tau, D=D, m_kin=m_kin, cusp_c=cu, b=bu, mode=mode, done=bool(s <= sig_goal or len(Ea) * tau >= T_total))
        if hbin:
            res[key]["intracule_hist"] = hist.tolist(); res[key]["hbin"] = hbin
        save_json(JS, res)
        print(key, res[key] if not hbin else {k: v for k, v in res[key].items() if k != "intracule_hist"}, flush=True)
    else:
        print(key, "equilibratura", nstep_done, flush=True)

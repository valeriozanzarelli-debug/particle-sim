"""multi_e_lattice.py - R4: la stessa regola (livello 1 simmetrizzato, catena A) in due direzioni del tempo, ione Z = 16
su 192^3 celle (cella a0/256, scatola +-6 a_Z), nucleo su uno spigolo di cella, V = -Z <1/r>_cella.
  det   (a) LOBPCG di H = T + V (T = banda imparata H(k) - H(0) in Eh);
        (b) tempo immaginario: autovalore massimo di T_beta = e^{-beta V/2} K_beta e^{-beta V/2}, K_beta = IFFT e^{-beta T(k)},
            a beta = 2e-4, 1e-4, 5e-5 (iterazione di potenza accelerata con Lanczos sullo stesso operatore), E(beta) = -ln
            lambda/beta, estrapolazione quadratica a beta -> 0; stimatore misto deterministico <H psi_T, phi_beta>/<psi_T, phi_beta>;
        (c) tempo reale: -arg<psi|W|psi>/tau, W = e^{-iV tau/2} e^{-iT tau} e^{-iV tau/2}, tau = tick grosso 2.85e-5 a.u.
  walk <beta>  camminatori sul reticolo (posizioni intere, salti dalla tabella alias di K_beta, branching di Strang),
        stimatore misto con guida NON esatta psi_T = exp(-Z m_kin r).
Risultati: lab/out/multi_e_lattice.json. Run: python lab/multi_e_lattice.py det ; python lab/multi_e_lattice.py walk 1e-4"""
import os, sys, time
from multi_e_core import OUT, A0, ALPHA, THETA, Band, chain_levels, load_json, save_json
import numpy as np
import scipy.fft as sfft
import scipy.sparse.linalg as spl
from agg_zz_atom import cell_avg_inv_r

JS = "multi_e_lattice.json"
Z = 16
N = 192
WK = 2
res = load_json(JS)


def setup():
    L = chain_levels("A")[0]
    b = Band(L["K"], L["offs"], L["B"], THETA["A"], sym=True)
    hj = b.hj_eh()
    offs = np.asarray(L["offs"]).astype(int)
    A = np.zeros((N, N, N))
    np.add.at(A, tuple(offs[:, i] % N for i in range(3)), hj)
    Hk = np.real(sfft.fftn(A, workers=WK).conj())
    Tk = Hk - hj.sum()                                   # energia di moto in Eh, T(0) = 0
    inv, x = cell_avg_inv_r(N, b.h, "corner")
    V = -Z * inv
    return b, Tk, V, x


def fft(a):
    return sfft.fftn(a, workers=WK)


def ifft(a):
    return sfft.ifftn(a, workers=WK)


sec = sys.argv[1] if len(sys.argv) > 1 else "det"

if sec == "det":
    t0 = time.time()
    b, Tk, V, x = setup()
    X, Y, Zc = np.meshgrid(x, x, x, indexing="ij", sparse=True)
    r = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
    m_kin = 1.0003832623372946
    out = res.get("det", {})
    H = lambda f: np.real(ifft(Tk * fft(f))) + V * f
    N3 = N ** 3
    pre = 1 / (Tk + 0.5 * Z ** 2)

    def Aop(Xm):
        Xm = np.asarray(Xm).reshape(N3, -1)
        return np.stack([H(Xm[:, c].reshape(N, N, N)).ravel() for c in range(Xm.shape[1])], 1)

    def Pop(Xm):
        Xm = np.asarray(Xm).reshape(N3, -1)
        return np.stack([np.real(ifft(pre * fft(Xm[:, c].reshape(N, N, N)))).ravel() for c in range(Xm.shape[1])], 1)
    vf = os.path.join(OUT, "multi_e_lattice_psi.npy")
    if "E_lobpcg" not in out:
        g0 = np.exp(-Z * r).ravel()[:, None]
        w, v, hist = spl.lobpcg(spl.LinearOperator((N3, N3), matvec=Aop, matmat=Aop, dtype=float), g0,
                                M=spl.LinearOperator((N3, N3), matvec=Pop, matmat=Pop, dtype=float), tol=1e-11, maxiter=200,
                                largest=False, retResidualNormsHistory=True)
        psi = v[:, 0].reshape(N, N, N); psi /= np.linalg.norm(psi)
        Hp = H(psi)
        E = float(np.sum(psi * Hp)); resid = float(np.linalg.norm(Hp - E * psi))
        np.save(vf, psi)
        out.update(E_lobpcg=E, lobpcg_resid=resid, lobpcg_iters=len(hist), hydrogenic_minus_mZ2_2=-m_kin * Z ** 2 / 2,
                   rel_to_hydrogenic=E / (-m_kin * Z ** 2 / 2) - 1, sec_lobpcg=time.time() - t0)
        res["det"] = out; save_json(JS, res)
        print(f"LOBPCG: E = {E:.10f} Eh (idrogenoide -m Z^2/2 = {-m_kin * Z ** 2 / 2:.6f}), residuo {resid:.1e}, "
              f"{len(hist)} it, {time.time() - t0:.0f} s", flush=True)
    psi = np.load(vf)
    E0 = out["E_lobpcg"]
    if "decomp" not in out:
        # scomposizione dello scarto dall'idrogenoide: parte non quadratica della banda (k^4 imparato) e resto
        fp = fft(psi); w2 = np.abs(fp) ** 2 / np.sum(np.abs(fp) ** 2)
        k1 = 2 * np.pi * np.fft.fftfreq(N, b.h)
        KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij", sparse=True)
        k2 = KX ** 2 + KY ** 2 + KZ ** 2
        Tb = float(np.sum(w2 * Tk)); Tq = float(np.sum(w2 * k2 / (2 * m_kin))); Vv = float(np.sum(psi ** 2 * V))
        k4 = float(np.sum(w2 * k2 ** 2))
        out["decomp"] = dict(T_band=Tb, T_quadratic=Tq, band_nonquadratic=Tb - Tq, V=Vv,
                             dirac_mass_velocity_on_this_state=-(1 / 137.035999084) ** 2 * k4 / 8,
                             dirac_mv_hydrogenic=-5 / 8 * Z ** 4 / 137.035999084 ** 2,
                             deviation_from_hydrogenic=E0 + m_kin * Z ** 2 / 2)
        res["det"] = out; save_json(JS, res)
        print("decomp", out["decomp"], flush=True)
    psiT = np.exp(-Z * m_kin * r)
    HpsiT = H(psiT)
    for beta in (2e-4, 1e-4, 5e-5):
        key = f"beta{beta:g}"
        if key in out:
            continue
        t1 = time.time()
        Kk = np.exp(-beta * Tk)
        hv = np.exp(-beta * V / 2)
        Top = lambda f: (hv * np.real(ifft(Kk * fft(hv * f.reshape(N, N, N))))).ravel()
        lam, vv = spl.eigsh(spl.LinearOperator((N3, N3), matvec=Top, dtype=float), k=1, which="LA", v0=psi.ravel(), tol=1e-15)
        phi = vv[:, 0].reshape(N, N, N); phi *= np.sign(phi.sum())
        Eb = float(-np.log(lam[0]) / beta)
        Emix = float(np.sum(HpsiT * phi) / np.sum(psiT * phi))
        out[key] = dict(E_growth=Eb, E_mixed_det=Emix, lam=float(lam[0]), sec=time.time() - t1)
        np.save(os.path.join(OUT, f"multi_e_lattice_phi_{key}.npy"), phi)
        res["det"] = out; save_json(JS, res)
        print(f"beta {beta:g}: -ln(lambda)/beta = {Eb:.10f}, misto det. {Emix:.8f} ({time.time() - t1:.0f} s)", flush=True)
    bs = np.array([2e-4, 1e-4, 5e-5]); Es = np.array([out[f"beta{x:g}"]["E_growth"] for x in bs])
    c = np.polyfit(bs, Es, 2)
    Eext = float(c[-1])
    # (c) tempo reale
    tau = b.tau_au
    Wf = np.exp(-0.5j * V * tau) * ifft(np.exp(-1j * Tk * tau) * fft(np.exp(-0.5j * V * tau) * psi))
    Ert = float(-np.angle(np.vdot(psi, Wf)) / tau)
    out.update(E_imag_beta0=Eext, quad_fit=c.tolist(), rel_imag_vs_lobpcg=abs(Eext - E0) / abs(E0), E_realtime=Ert, tau_au=tau,
               rel_realtime_vs_lobpcg=abs(Ert - E0) / abs(E0))
    out["PASS_a"] = bool(out["rel_imag_vs_lobpcg"] < 1e-6 and out["rel_realtime_vs_lobpcg"] < 3e-5)
    res["det"] = out; save_json(JS, res)
    print(f"R4a: imag (beta->0) {Eext:.10f}, LOBPCG {E0:.10f} (rel {out['rel_imag_vs_lobpcg']:.1e}); tempo reale {Ert:.10f} "
          f"(rel {out['rel_realtime_vs_lobpcg']:.1e}) -> {'PASS' if out['PASS_a'] else 'FAIL'}", flush=True)

# ============================================================================================ walkers (R4 b)
if sec == "walk":
    import numba as nb

    beta = float(sys.argv[2])
    nsteps_target = int(float(sys.argv[3])) if len(sys.argv) > 3 else 40000
    NW = 20000
    key = f"walk_beta{beta:g}"
    b, Tk, V, x = setup()
    m_kin = 1.0003832623372946
    X, Y, Zc = np.meshgrid(x, x, x, indexing="ij", sparse=True)
    r = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
    psiT = np.exp(-Z * m_kin * r)
    EL = (np.real(ifft(Tk * fft(psiT))) + V * psiT) / psiT          # energia locale della guida sul reticolo
    # nucleo K_beta nello spazio reale, supporto +-6 sigma, pesi negativi (<1e-13) azzerati e rinormalizzato
    Kx = np.real(ifft(np.exp(-beta * Tk)))
    sig = np.sqrt(2 * 0.4998084422 * beta) / b.h
    R = int(np.ceil(6 * sig))
    idx = np.arange(-R, R + 1)
    sub = Kx[np.ix_(idx % N, idx % N, idx % N)]
    negw = float(-sub[sub < 0].sum() / sub.sum()); lost = float(1 - sub.sum())
    p = np.clip(sub, 0, None).ravel(); p /= p.sum()
    J = np.stack(np.meshgrid(idx, idx, idx, indexing="ij"), -1).reshape(-1, 3)
    # tabella alias (Walker/Vose)
    n = len(p); prob = np.zeros(n); alias = np.zeros(n, np.int64)
    sc = p * n
    small = [i for i in range(n) if sc[i] < 1]; large = [i for i in range(n) if sc[i] >= 1]
    while small and large:
        s_, l_ = small.pop(), large.pop()
        prob[s_] = sc[s_]; alias[s_] = l_
        sc[l_] = sc[l_] - (1 - sc[s_])
        (small if sc[l_] < 1 else large).append(l_)
    for i in large + small:
        prob[i] = 1.0
    Vf = V.ravel().copy(); ELf = EL.ravel().copy(); PTf = psiT.ravel().copy()

    @nb.njit(cache=True)
    def run(pos, nsteps, beta, prob, alias, J, Vf, ELf, PTf, Eref0, NW, Tpop, seed, rec_every):
        np.random.seed(seed)
        n = len(prob)
        Eref = Eref0
        nrec = nsteps // rec_every
        num = np.zeros(nrec); den = np.zeros(nrec); pop = np.zeros(nrec); eref = np.zeros(nrec)
        Nn = 192
        Ebest = Eref0; cn = 0.0; cd = 0.0
        newpos = np.empty((4 * NW, 3), np.int64)
        for st in range(nsteps):
            nw = len(pos)
            w = np.empty(nw)                                  # (prima versione: w di taglia fissa -> scrittura fuori bordo)
            for i in range(nw):
                ix = pos[i, 0] * Nn * Nn + pos[i, 1] * Nn + pos[i, 2]
                wi = np.exp(-beta * (0.5 * Vf[ix] - Eref))
                u = np.random.random() * n
                k = int(u)
                if u - k >= prob[k]:
                    k = alias[k]
                for d in range(3):
                    pos[i, d] = (pos[i, d] + J[k, d]) % Nn
                ix = pos[i, 0] * Nn * Nn + pos[i, 1] * Nn + pos[i, 2]
                w[i] = wi * np.exp(-beta * 0.5 * Vf[ix])
            # stimatore misto (pesi prima del branching)
            sn = 0.0; sd = 0.0; sw = 0.0
            for i in range(nw):
                ix = pos[i, 0] * Nn * Nn + pos[i, 1] * Nn + pos[i, 2]
                sn += w[i] * PTf[ix] * ELf[ix]; sd += w[i] * PTf[ix]; sw += w[i]
            cn = 0.999 * cn + sn / nw; cd = 0.999 * cd + sd / nw          # media mobile (~1000 passi)
            Ebest = cn / cd
            rr = st // rec_every
            if rr < nrec:
                num[rr] += sn / nw; den[rr] += sd / nw; pop[rr] += nw; eref[rr] += Eref
            # branching (copie intere = floor(w + u))
            m = 0
            for i in range(nw):
                c = int(w[i] + np.random.random())
                for _ in range(c):
                    if m < 4 * NW:
                        newpos[m, 0] = pos[i, 0]; newpos[m, 1] = pos[i, 1]; newpos[m, 2] = pos[i, 2]; m += 1
            pos = newpos[:m].copy()
            Eref = Ebest - np.log(m / NW) / Tpop            # controllo proporzionale della popolazione, T_pop in a.u.
        return pos, num, den, pop, eref, Eref

    st = res.get(key, {})
    ck = os.path.join(OUT, f"multi_e_lattice_walk_{key}.npz")
    if os.path.exists(ck):
        z = np.load(ck); pos = z["pos"]; num_all = list(z["num"]); den_all = list(z["den"]); Eref = float(z["Eref"])
    else:
        rng = np.random.default_rng(7)
        pp = (psiT ** 2).ravel(); pp /= pp.sum()
        ii = rng.choice(N ** 3, NW, p=pp)
        pos = np.stack(np.unravel_index(ii, (N, N, N)), 1).astype(np.int64)
        num_all, den_all = [], []
        Eref = -Z ** 2 / 2 * m_kin
    t0 = time.time()
    rec_every = int(round(0.002 / beta))                  # un record ogni 0.002 a.u.
    done_steps = len(num_all) * rec_every
    while done_steps < nsteps_target and time.time() - t0 < 480:
        chunk = max(rec_every, (min(2000 * rec_every, nsteps_target - done_steps) // rec_every) * rec_every)
        pos, num, den, pop, eref, Eref = run(pos, chunk, beta, prob, alias, J, Vf, ELf, PTf, Eref, NW, 1.0,
                                             1000 + done_steps, rec_every)
        num_all += list(num / rec_every); den_all += list(den / rec_every)
        done_steps += chunk
        np.savez(ck, pos=pos, num=np.array(num_all), den=np.array(den_all), Eref=Eref)
        print(f"  {done_steps} passi ({done_steps * beta:.2f} a.u.), pop {len(pos)}, Eref {Eref:.4f}, {time.time() - t0:.0f} s", flush=True)
    num_a, den_a = np.array(num_all), np.array(den_all)
    neq = int(0.05 / 0.002)                               # 0.05 a.u. di equilibratura scartati
    nn, dd = num_a[neq:], den_a[neq:]
    E = nn.sum() / dd.sum()
    errs = []
    for nb_ in (32, 16, 8):
        L_ = len(nn) // nb_
        if L_ < 1:
            continue
        Eb = np.array([nn[i * L_:(i + 1) * L_].sum() / dd[i * L_:(i + 1) * L_].sum() for i in range(nb_)])
        errs.append(Eb.std(ddof=1) / np.sqrt(nb_))
    sigma = float(max(errs)) if errs else float("nan")
    det = res.get("det", {}).get(f"beta{beta:g}", {})
    Edet = det.get("E_mixed_det")
    st = dict(E_mixed_walkers=float(E), sigma=sigma, sigma_rel=sigma / abs(E), steps=done_steps, time_au=done_steps * beta,
              walkers=NW, kernel_support_cells=R, kernel_neg_weight=negw, kernel_mass_outside=lost, E_mixed_det=Edet,
              z=(float(E) - Edet) / sigma if Edet is not None else None)
    if Edet is not None:
        st["PASS_b"] = bool(abs(E - Edet) <= 3 * sigma and sigma / abs(E) <= 1e-4)
    res[key] = st; save_json(JS, res)
    print(st, flush=True)

# ============================================================================================ diagnostica post hoc R4a
if sec == "diag":
    # (dopo il FAIL di R4a) l'errore di Strang e' pari in beta (nessun termine lineare) e grande a Z = 16: si scende in
    # beta (2.5e-5, 1.25e-5; il nucleo non e' positivo li', ma l'iterazione deterministica non ne ha bisogno) e in tau
    # per il tempo reale (tau/2, tau/4, tau/8), e si estrapola con la forma pari E0 + c2 x^2 + c4 x^4. Non rigiudica R4a.
    b, Tk, V, x = setup()
    out = res.get("diag_R4a", {})
    psi = np.load(os.path.join(OUT, "multi_e_lattice_psi.npy"))
    N3 = N ** 3
    for beta in (2.5e-5, 1.25e-5):
        key = f"beta{beta:g}"
        if key in out:
            continue
        t1 = time.time()
        Kk = np.exp(-beta * Tk); hv = np.exp(-beta * V / 2)
        Top = lambda f: (hv * np.real(ifft(Kk * fft(hv * f.reshape(N, N, N))))).ravel()
        lam, vv = spl.eigsh(spl.LinearOperator((N3, N3), matvec=Top, dtype=float), k=1, which="LA", v0=psi.ravel(), tol=1e-15)
        out[key] = dict(E_growth=float(-np.log(lam[0]) / beta), sec=time.time() - t1)
        res["diag_R4a"] = out; save_json(JS, res)
        print(key, out[key], flush=True)
    tau0 = b.tau_au
    rt = {}
    for f_ in (1, 2, 4, 8):
        tau = tau0 / f_
        Wf = np.exp(-0.5j * V * tau) * ifft(np.exp(-1j * Tk * tau) * fft(np.exp(-0.5j * V * tau) * psi))
        rt[f"{tau:g}"] = float(-np.angle(np.vdot(psi, Wf)) / tau)
    out["realtime"] = rt
    E0 = res["det"]["E_lobpcg"]
    bs = [2e-4, 1e-4, 5e-5, 2.5e-5, 1.25e-5]
    Eb = [res["det"][f"beta{x:g}"]["E_growth"] if f"beta{x:g}" in res["det"] else out[f"beta{x:g}"]["E_growth"] for x in bs]
    bb = np.array(bs[-3:]); ee = np.array(Eb[-3:])
    ci = np.linalg.solve(np.stack([np.ones(3), bb ** 2, bb ** 4], 1), ee)
    tt = np.array([tau0 / 2, tau0 / 4, tau0 / 8]); er = np.array([rt[f"{x:g}"] for x in tt])
    cr = np.linalg.solve(np.stack([np.ones(3), tt ** 2, tt ** 4], 1), er)
    out.update(imag_even_fit_last3=ci.tolist(), imag_rel=float(abs(ci[0] - E0) / abs(E0)),
               real_even_fit_last3=cr.tolist(), real_rel=float(abs(cr[0] - E0) / abs(E0)),
               c2_imag=float(ci[1]), c2_real=float(cr[1]),
               note="Wick: il coefficiente di x^2 nel tempo reale deve essere l'opposto di quello immaginario")
    res["diag_R4a"] = out; save_json(JS, res)
    print({k: v for k, v in out.items()}, flush=True)

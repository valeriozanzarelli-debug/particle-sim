"""aggmol_frat_run.py - study aggmol_frat_: production jobs (pre-registration lab/out/aggmol_frat_prereg.json).
  step1               dimers / trimers / B2 (aggmol_frat_pot.job_step1)
  clusters            Mackay 13/55/147: minimum + MD at 5, 10, 15, 20 K
  learn               blob level from the 147 cluster at 20 K (force matching + friction from the blob VACF)
  drop T              1000 molecules from a cube at T (checkpointed); per-frame positions saved
  coarse T            the blob rule on 77 blobs from the same cube (and on the 147 cluster)
  analyze             verdicts -> lab/out/aggmol_frat.json
Big arrays (.npy/.npz) are not committed; small JSON frames for the glow engine are."""
import os, sys, json, time
import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_frat_pot as P
import aggmol_frat_md as D

OUT = P.OUT
JS = os.path.join(OUT, "aggmol_frat.json")
M = P.M_H2_AMU
DT = 0.004
GAMMA = 0.2
BLOB = 13
TPS = D.TU_PS


def load():
    return json.load(open(JS)) if os.path.exists(JS) else {}


def save(r):
    json.dump(r, open(JS, "w"), indent=1, default=float)


def table():
    r, v, dv = P.table_K_A()
    return D.Table(r, v, dv)


def well_A():
    r, v, dv = P.table_K_A()
    i = int(np.argmin(v))
    return float(r[i]), float(v[i])


def lindemann(traj):
    """traj (nf, N, 3): sqrt(<r^2>-<r>^2)/<r> averaged over pairs"""
    nf, N, _ = traj.shape
    iu = np.triu_indices(N, 1)
    s = np.zeros(len(iu[0])); s2 = np.zeros_like(s)
    for x in traj:
        d = np.linalg.norm(x[iu[0]] - x[iu[1]], axis=1)
        s += d; s2 += d * d
    s /= nf; s2 /= nf
    return float(np.mean(np.sqrt(np.maximum(s2 - s * s, 0)) / s))


def largest_cluster(x, link=5.2):
    pr = cKDTree(x).query_pairs(link, output_type="ndarray")
    A = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(x), len(x)))
    n, lab = connected_components(A, directed=False)
    big = np.bincount(lab).argmax()
    return lab == big


def job_step1():
    P.job_step1()


def job_clusters():
    res = load(); out = res.setdefault("clusters", {})
    tab = table(); rm, vm = well_A()
    for ns in (1, 2, 3):
        x0 = D.mackay(ns, rm)
        x, e = D.minimise(x0, tab)
        N = len(x)
        row = out.setdefault(str(N), {}); row["Emin_per"] = e / N
        np.save(os.path.join(OUT, f"aggmol_frat_mackay{N}.npy"), x)
        for T in (5.0, 10.0, 15.0, 20.0):
            if f"{T:.0f}K" in row:
                continue
            S = D.System(x.copy(), M, tab); rng = np.random.default_rng(int(10 * T) + N)
            S.v = rng.normal(size=S.x.shape) * np.sqrt(T / M)
            nst = int(round(110 / TPS / DT)); neq = nst // 5
            S.baoab(DT, T, GAMMA, rng, neq)
            E, traj = [], []
            for k in range((nst - neq) // 25):
                S.baoab(DT, T, GAMMA, rng, 25)
                E.append((S.epot + S.ekin()) / N)
                if k % 4 == 0:
                    traj.append(S.x - S.x.mean(0))
            big = largest_cluster(S.x)
            row[f"{T:.0f}K"] = dict(E_per=float(np.mean(E)), Epot_per=float(S.epot / N), lindemann=lindemann(np.array(traj)),
                                     bound_fraction=float(big.mean()))
            print(N, T, row[f"{T:.0f}K"], flush=True)
            save(res)
    save(res)


# ------------------------------------------------------------------------------------------- blob level
def blobs(x, k, cen=None, rng=None):
    if cen is None:
        cen, lab = kmeans2(x, k, minit="++", seed=rng)
    else:
        cen, lab = kmeans2(x, cen, minit="matrix", iter=20)
    return cen, lab


def blob_data(x, f, v, lab, k):
    m = np.bincount(lab, minlength=k).astype(float)
    C = np.array([np.bincount(lab, x[:, a], k) for a in range(3)]).T / np.maximum(m, 1)[:, None]
    F = np.array([np.bincount(lab, f[:, a], k) for a in range(3)]).T
    V = np.array([np.bincount(lab, v[:, a], k) for a in range(3)]).T / np.maximum(m, 1)[:, None]
    return C, F, V, m


HATS = np.arange(2.0, 30.0 + 1e-9, 1.0)
LAMBDAS = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
WINDOW = int(os.environ.get("AGGMOL_FRAT_W", 10))     # frames per window (0.22 ps each)
NRUNS = 3           # smoothness penalty, chosen by 2-fold cross-validation


def hat_design(C, k):
    """design matrix for F_I = sum_J phi(R_IJ) R_IJ_hat, phi = sum_b w_b hat_b(R)"""
    i, j = np.triu_indices(k, 1)
    d = C[j] - C[i]; r = np.linalg.norm(d, axis=1); u = d / r[:, None]
    h = HATS[1] - HATS[0]
    nb = len(HATS)
    X = np.zeros((k, 3, nb))
    pos = (r - HATS[0]) / h
    b0 = np.floor(pos).astype(int); fr = pos - b0
    for bb, ww in ((b0, 1 - fr), (b0 + 1, fr)):
        ok = (bb >= 0) & (bb < nb)
        for a in range(3):
            # phi > 0 = repulsive: force on I along -u, on J along +u
            np.add.at(X[:, a, :], (i[ok], bb[ok]), -ww[ok] * u[ok, a])
            np.add.at(X[:, a, :], (j[ok], bb[ok]), ww[ok] * u[ok, a])
    return X.reshape(3 * k, nb)


def job_learn():
    res = load()
    tab = table()
    x = np.load(os.path.join(OUT, "aggmol_frat_mackay147.npy"))
    N = len(x); k = int(round(N / BLOB)); T = 20.0
    nfr = 2000; every = int(round(0.22 / TPS / DT))
    W = WINDOW                                            # time window: membership fixed, forces and centres averaged
    nb = len(HATS)
    XtX = [np.zeros((nb, nb)), np.zeros((nb, nb))]; Xty = [np.zeros(nb), np.zeros(nb)]; yy = [0.0, 0.0]
    Vtr = []; sizes = []; dists = []
    for run in range(NRUNS):
        rng = np.random.default_rng(147 + run)
        S = D.System(x.copy(), M, tab); S.v = rng.normal(size=S.x.shape) * np.sqrt(T / M)
        S.baoab(DT, T, GAMMA, rng, int(round(22 / TPS / DT)))
        cen, lab = blobs(S.x, k, rng=rng)
        for t in range(nfr // W):
            cen, lab = blobs(S.x, k, cen=cen)
            Cs = 0; Fs = 0
            for q in range(W):
                S.baoab(DT, T, GAMMA, rng, every)
                C, F, V, m = blob_data(S.x, S.f, S.v, lab, k)
                Cs = Cs + C / W; Fs = Fs + F / W
                if run == 0:
                    Vtr.append(V)
            X = hat_design(Cs, k)
            h = (t // max(1, 200 // W)) % 2                 # two folds of interleaved 44 ps blocks
            XtX[h] += X.T @ X; Xty[h] += X.T @ Fs.ravel(); yy[h] += float(Fs.ravel() @ Fs.ravel())
            sizes.append(m)
            dists.append(np.linalg.norm(Cs[:, None] - Cs[None], axis=-1)[np.triu_indices(k, 1)])
    D2 = np.diff(np.eye(nb), 2, axis=0); Pen = D2.T @ D2 + 1e-2 * np.eye(nb)
    sc = np.trace(XtX[0] + XtX[1]) / nb
    cv = []
    for lam in LAMBDAS:
        err = 0.0
        for h in (0, 1):
            w_ = np.linalg.solve(XtX[h] + lam * sc * Pen + 1e-9 * sc * np.eye(nb), Xty[h])
            o = 1 - h
            err += yy[o] - 2 * w_ @ Xty[o] + w_ @ XtX[o] @ w_
        cv.append(err / (yy[0] + yy[1]))
    lam = LAMBDAS[int(np.argmin(cv))]
    XtXa = XtX[0] + XtX[1]; Xtya = Xty[0] + Xty[1]
    w = np.linalg.solve(XtXa + lam * sc * Pen + 1e-9 * sc * np.eye(nb), Xtya)
    used = np.diag(XtXa) > 1e-6 * np.max(np.diag(XtXa))
    resid = 1 - (2 * w @ Xtya - w @ XtXa @ w) / (yy[0] + yy[1])
    # blob VACF -> friction (exp decay time); identities carried by the k-means centres
    Vtr = np.array(Vtr)
    v0 = Vtr - Vtr.mean(1, keepdims=True)
    c = np.array([np.mean(np.sum(v0[:-l or None] * v0[l:], -1)) for l in range(40)]); c /= c[0]
    lag = np.arange(40) * every * DT
    tau = float(lag[np.argmax(c < np.exp(-1))]) if np.any(c < np.exp(-1)) else float(lag[-1])
    res["learn"] = dict(k=k, frames=nfr, window=W, runs=NRUNS, lambda_cv=lam, cv=cv, force_residual_fraction=float(resid), w=w.tolist(), used=used.tolist(), hats=HATS.tolist(), vacf=c.tolist(),
                        vacf_lag_tu=lag.tolist(), tau_tu=tau, gamma_B=1 / tau, mean_size=float(np.mean(sizes)),
                        size_spread=float(np.std(sizes)), dmin=float(np.min(np.concatenate(dists))))
    save(res)
    print("learned", dict(k=k, tau_ps=tau * TPS, mean_size=np.mean(sizes), dmin=res["learn"]["dmin"]), flush=True)


def blob_table(res):
    L = res["learn"]
    h = np.array(L["hats"]); w = np.array(L["w"]); used = np.array(L["used"])
    # inside the sampled range: phi; below it: continue the steepest repulsion linearly (declared extrapolation)
    lo = np.argmax(used)
    r = np.linspace(0.5, h[-1], 6000)
    phi = np.interp(r, h, w)
    sl = max(abs((w[lo + 1] - w[lo]) / (h[1] - h[0])), 10.0)
    phi = np.where(r < h[lo], max(w[lo], 0.0) + sl * (h[lo] - r), phi)     # rising wall below the sampled range
    # potential U(r) = int_r^rmax phi (phi = -dU/dr ... here phi>0 repulsive means dU/dr = -phi)
    U = np.cumsum(phi[::-1])[::-1] * (r[1] - r[0])
    U -= U[-1]
    return D.Table(r, U, -phi)


def cube1000():
    a = (1 / 0.0211) ** (1 / 3)
    g = np.arange(10) * a
    return np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T


def job_drop(T, budget=430):
    t0 = time.time()
    res = load(); key = f"drop_{T:.0f}K"
    ck = os.path.join(OUT, f"aggmol_frat_drop{T:.0f}_ck.npz")
    tab = table()
    nst = int(round(165 / TPS / DT)); nav0 = int(round(55 / TPS / DT)); blk = 100
    if os.path.exists(ck):
        z = np.load(ck); x, v, step = z["x"], z["v"], int(z["step"]); acc = json.loads(str(z["acc"]))
        rng = np.random.default_rng(int(z["seed"]) + step)
    else:
        x = cube1000(); rng = np.random.default_rng(int(T) * 7); v = rng.normal(size=x.shape) * np.sqrt(T / M); step = 0
        acc = dict(t=[], rg=[], asph=[], big=[], epot=[], frames=[])
    S = D.System(x, M, tab); S.v = v
    frames_dir = os.path.join(OUT, f"aggmol_frat_drop{T:.0f}_traj.npy")
    traj = list(np.load(frames_dir)) if os.path.exists(frames_dir) and step > 0 else []
    while step < nst:
        if time.time() - t0 > budget:
            np.savez(ck, x=S.x, v=S.v, step=step, acc=json.dumps(acc), seed=int(T) * 7)
            np.save(frames_dir, np.array(traj))
            print(f"{key}: checkpoint at {step * DT * TPS:.1f} ps", flush=True); return
        S.baoab(DT, T, GAMMA, rng, blk); step += blk
        big = largest_cluster(S.x)
        rg, asp = D.gyration(S.x[big])
        acc["t"].append(step * DT * TPS); acc["rg"].append(rg); acc["asph"].append(asp); acc["big"].append(float(big.mean()))
        acc["epot"].append(S.epot / S.N)
        if step % (blk * 10) == 0:
            traj.append(S.x.astype(np.float32))
    np.save(frames_dir, np.array(traj))
    # final analysis on the stored frames of the last 110 ps
    tr = np.array(traj); tt = np.arange(1, len(tr) + 1) * blk * 10 * DT * TPS
    use = tr[tt > 55]
    out = analyse_drop(use, tab, T)
    sel = np.array(acc["t"]) > 55
    out.update(rg=float(np.mean(np.array(acc["rg"])[sel])), asph=float(np.mean(np.array(acc["asph"])[sel])),
               big=float(np.mean(np.array(acc["big"])[sel])), epot_per=float(np.mean(np.array(acc["epot"])[sel])),
               series=dict(t=acc["t"][::5], rg=acc["rg"][::5], asph=acc["asph"][::5], big=acc["big"][::5]))
    res = load(); res[key] = out; save(res)
    if os.path.exists(ck):
        os.remove(ck)
    # small JSON frames for the glow engine (every ~11 ps, positions in A to 0.01)
    fr = tr[::max(1, len(tr) // 15)]
    json.dump(dict(units="angstrom", T_K=T, dt_ps=float(blk * 10 * DT * TPS * max(1, len(tr) // 15)),
                   note="classical H2 molecule-balls, liquid H2 standing in for water; regenerate: python lab/aggmol_frat_run.py drop %g" % T,
                   frames=[np.round(f - f.mean(0), 2).tolist() for f in fr]),
              open(os.path.join(OUT, f"aggmol_frat_drop{T:.0f}_frames.json"), "w"))
    print(key, {k: v for k, v in out.items() if k != "series"}, flush=True)


def analyse_drop(tr, tab, T):
    rhos, peaks, ucs, Us, Res, atms = [], [], [], [], [], []
    for n, x in enumerate(tr):
        x = x.astype(float)
        big = largest_cluster(x); xb = x[big]
        S = D.System(xb, M, tab)
        e = S.epot_per()
        c = xb.mean(0); r = np.linalg.norm(xb - c, axis=1)
        Rg = np.sqrt(np.mean(r ** 2)); Rs = np.sqrt(5 / 3) * Rg
        core = r < 0.5 * Rs
        rho = core.sum() / (4 / 3 * np.pi * (0.5 * Rs) ** 3)
        Re = (3 * len(xb) / (4 * np.pi * rho)) ** (1 / 3)
        rr, g = D.gofr_core(xb, 0.5 * Rs, 8.0, 160)
        peaks.append(float(rr[np.argmax(g)]))
        rhos.append(rho); ucs.append(float(np.mean(e[core]))); Us.append(float(S.epot)); Res.append(Re)
        if n % max(1, len(tr) // 5) == 0 and len(atms) < 5:
            atms.append(atm_energy(xb) / len(xb))
    rho = float(np.mean(rhos)); uc = float(np.mean(ucs))
    gcm3 = rho * M / 0.602214076
    # tail beyond the cutoff (bulk, analytic, -C6/r^6): per molecule 2 pi rho int_rc^inf V r^2 dr
    C6K = P.Pair().C6 * P.EH_K * P.A0_A ** 6
    rc = tab.rmax
    tail = -2 * np.pi * rho * C6K / (3 * rc ** 3)
    # the shift of the table at rc also lowers u by v(rc)/2 per neighbour beyond... (table shifted to 0: add back)
    vrc = float(P.Pair().V0(P.RCUT)[0] * P.EH_K)                # the table is shifted by -v(rc) per pair
    shift = 0.5 * rho * 4 / 3 * np.pi * rc ** 3 * vrc
    tail = tail + shift
    dH = -uc - tail + T
    Re = float(np.mean(Res)); N = len(tr[0])
    gam = (np.mean(Us) - N * uc) / (4 * np.pi * Re ** 2)
    return dict(rho_core_A3=rho, rho_core_gcm3=gcm3, g_peak_A=float(np.mean(peaks)), u_core=uc, tail_K=float(tail),
                dHvap_K=float(dH), dHvap_kJmol=float(dH * 8.314462618e-3), R_e=Re, gamma_K_A2=float(gam),
                gamma_mNm=float(gam * 1.380649), atm_per_K=float(np.mean(atms)) if atms else None,
                nframes=len(tr))


def atm_energy(x, rcut=10.0):
    C9K = P.Pair().C9 * P.EH_K * P.A0_A ** 9
    tr = cKDTree(x)
    nb = tr.query_ball_point(x, rcut)
    E = 0.0
    for i in range(len(x)):
        js = [j for j in nb[i] if j > i]
        if len(js) < 2:
            continue
        js = np.array(js)
        a, b = np.triu_indices(len(js), 1)
        j, k = js[a], js[b]
        ok = np.linalg.norm(x[j] - x[k], axis=1) < rcut
        j, k = j[ok], k[ok]
        r12 = np.linalg.norm(x[i] - x[j], axis=1); r13 = np.linalg.norm(x[i] - x[k], axis=1)
        r23 = np.linalg.norm(x[j] - x[k], axis=1)
        E += float(np.sum(P.atm(C9K, r12, r13, r23)))
    return E


def job_coarse(T):
    res = load()
    L = res["learn"]
    btab = blob_table(res)
    out = {}
    for name, x0, nt in (("drop1000", cube1000(), 165.0), ("cluster147", np.load(os.path.join(OUT, "aggmol_frat_mackay147.npy")), 110.0)):
        k = int(round(len(x0) / BLOB))
        rng = np.random.default_rng(int(T) + k)
        cen, lab = blobs(x0, k, rng=rng)
        m = np.bincount(lab, minlength=k) * M
        S = D.System(cen.copy(), m.astype(float), btab, skin=2.0)
        S.v = rng.normal(size=S.x.shape) * np.sqrt(T / m[:, None])
        dtc = 0.01; nst = int(round(nt / TPS / dtc))
        rgs, ds = [], []
        for s in range(nst // 50):
            S.baoab(dtc, T, L["gamma_B"], rng, 50)
            if s * 50 * dtc * TPS > nt / 3:
                rgs.append(D.gyration(S.x)[0])
                ds.append(np.linalg.norm(S.x[:, None] - S.x[None], axis=-1)[np.triu_indices(k, 1)])
        out[name] = dict(k=k, rg_blob=float(np.mean(rgs)), g_peak=float(first_peak(np.concatenate(ds), k, len(ds))),
                         Tkin=float(2 * S.ekin() / (3 * k)))
        print("coarse", T, name, out[name], flush=True)
    res = load(); res[f"coarse_{T:.0f}K"] = out; save(res)


def first_peak(d, k, nf, rmax=25.0, nb=100):
    h, e = np.histogram(d, np.linspace(0, rmax, nb + 1))
    r = 0.5 * (e[1:] + e[:-1])
    g = h / (r ** 2)
    g = np.convolve(g, np.ones(3) / 3, "same")
    return r[np.argmax(g)]


def fine_blob_stats(T):
    """k-means blobs of the fine drop frames (same k as the coarse) -> blob-centre Rg and first blob-blob peak"""
    tr = np.load(os.path.join(OUT, f"aggmol_frat_drop{T:.0f}_traj.npy"))
    tt = np.arange(1, len(tr) + 1) * 100 * 10 * DT * TPS
    k = int(round(1000 / BLOB)); rng = np.random.default_rng(3)
    rgs, ds = [], []
    cen = None
    for x in tr[tt > 55].astype(float):
        big = largest_cluster(x)
        cen, lab = blobs(x[big], k, cen if cen is not None else None, rng=rng)
        C, _, _, m = blob_data(x[big], np.zeros_like(x[big]), np.zeros_like(x[big]), lab, k)
        rgs.append(D.gyration(C)[0])
        ds.append(np.linalg.norm(C[:, None] - C[None], axis=-1)[np.triu_indices(k, 1)])
    return dict(rg_blob=float(np.mean(rgs)), g_peak=float(first_peak(np.concatenate(ds), k, len(ds))))


def job_analyze():
    res = load()
    s1 = json.load(open(os.path.join(OUT, "aggmol_frat_step1.json")))
    V = {}
    dim = {round(d["R"], 2): d for d in s1["dimer"]}
    rel = [abs(dim[R]["quantum_J6"] - dim[R]["V0"]) / abs(dim[R]["V0"]) for R in (6.0, 6.5, 7.0)]
    V["S1"] = dict(rel=rel, PASS=bool(max(rel) <= 0.02))
    tri = {round(t["R"], 2): t for t in s1["trimer"]}
    rel = [abs(tri[R]["quantum_J4"] - tri[R]["ball"]) / abs(tri[R]["ball"]) for R in (6.5, 7.0)]
    V["S2"] = dict(rel=rel, PASS=bool(max(rel) <= 0.02))
    d20 = res["drop_20K"]
    V["D0"] = dict(big=d20["big"], asph=d20["asph"], PASS=bool(d20["big"] >= 0.95 and d20["asph"] <= 0.02))
    V["D1"] = dict(rho_gcm3=d20["rho_core_gcm3"], ratio=d20["rho_core_gcm3"] / 0.0708,
                   PASS=bool(1.10 <= d20["rho_core_gcm3"] / 0.0708 <= 1.50))
    V["D2"] = dict(peak_A=d20["g_peak_A"], PASS=bool(3.15 <= d20["g_peak_A"] <= 3.90))
    V["D3"] = dict(dHvap_kJmol=d20["dHvap_kJmol"], ratio=d20["dHvap_K"] / 108.3, PASS=bool(1.4 <= d20["dHvap_K"] / 108.3 <= 2.6))
    V["D4"] = dict(gamma_mNm=d20["gamma_mNm"], ratio=d20["gamma_mNm"] / 2.0, PASS=bool(1.2 <= d20["gamma_mNm"] / 2.0 <= 4.0))
    for T in (20, 15):
        if f"coarse_{T}K" in res and f"drop_{T}K" in res:
            fb = fine_blob_stats(T); res[f"fine_blobs_{T}K"] = fb
            c = res[f"coarse_{T}K"]["drop1000"]
            e1 = c["rg_blob"] / fb["rg_blob"] - 1; e2 = c["g_peak"] / fb["g_peak"] - 1
            if T == 20:
                V["L1"] = dict(coarse=c["rg_blob"], fine=fb["rg_blob"], rel=e1, PASS=bool(abs(e1) <= 0.10))
                V["L2"] = dict(coarse=c["g_peak"], fine=fb["g_peak"], rel=e2, PASS=bool(abs(e2) <= 0.10))
            else:
                V["L3_15K"] = dict(rg_rel=e1, peak_rel=e2)
    res["verdict"] = V; save(res)
    print(json.dumps(V, indent=1, default=float))


if __name__ == "__main__":
    j = sys.argv[1]
    if j in ("drop", "coarse"):
        dict(drop=job_drop, coarse=job_coarse)[j](float(sys.argv[2]))
    else:
        dict(step1=job_step1, clusters=job_clusters, learn=job_learn, analyze=job_analyze)[j]()

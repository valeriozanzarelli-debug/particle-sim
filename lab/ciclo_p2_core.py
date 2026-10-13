"""ciclo_p2_core.py - Phase 2 of CICLO_LIVELLI: the level-to-level loop in 1-D with AUTOMATIC coarse variables and a
choice of BLOCKING scheme: periodic b = 2, periodic b = 3, quasi-periodic Fibonacci (blocks of F_{j+1} and F_j cells:
mean factor phi = 1.618 per level).

One loop, identical at every level:
  1. slow set     : the current level's eigenmodes below a frequency cut omega_c (the fine dispersion at
                    k_c = pi / (2 B_next): half of the next level's Nyquist); positive band (the electron's sector)
  2. variables    : per block TYPE, the block density matrix of the slow set  rho = sum_m |v_m|_block><v_m|_block|;
                    keep the C leading eigenvectors, C = smallest with captured weight >= 1 - EPS_C   (automatic)
  3. rule         : coarse Hermitian generator H, nonzero only within w cells, entries tied by the local PATTERN of
                    block types (periodic: one pattern = a stencil); least squares H Q = Q Omega on the
                    Loewdin-orthonormalised coarse images Q of the slow modes (des Cloizeaux)
  4. check        : lowest eigenvalues of H against the level below (budget EPS_LEVEL, relative to the kinetic range
                    of the new level's quarter-Nyquist modes) + tame (no spurious state below the band bottom);
                    w = smallest in 1..WMAX that passes; none passes -> the level is not made (flag)
  5. repeat       : the next level uses the eigenmodes of the LEARNED H (never the fine rule again)
Potential (soft Coulomb): Galerkin restriction V_{n+1} = W V_n W^dag of the per-tick phase (no choice).
Frequencies are per fine tick (hbar = c = m = 1; tick = theta).
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import time
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

PHI = (1 + 5 ** 0.5) / 2
FRAC_CHECK = 1 / 12    # check: modes with >= 24 coarse cells per wavelength (the project's HD rule, CONTESTO s.6.1)
FRAC_SLOW = 1 / 6      # slow set (weight 1 in the fit, variables): 12 coarse cells per wavelength
NMIN_CELLS = 48        # stop when the next level would have fewer cells (box too small to judge)
MIN_CHECK = 9          # at least this many slowest modes in the check (small coarse lattices)
REL_FLOOR = 1e-2       # relative least squares: each mode's error counted against its own kinetic energy (+ floor)
W_FAST = 0.03          # weight of the non-slow modes of the sector in the fit (as agg_zz_core: keeps the rule tame)


# ------------------------------------------------------------------ fine 1-D walk
def fine_bloch_omega(k, th):
    """positive-branch frequency per tick: cos w = cos th cos k"""
    return np.arccos(np.cos(th) * np.cos(k))


def fine_modes_free(N, th, omega_c):
    """positive-band plane-wave eigenmodes of the free fine walk with omega <= omega_c: (2N, n) columns, omegas.
    Component order: index = 2 x + c (c = 0: R, 1: L)."""
    n = np.arange(N)
    kk = 2 * np.pi * np.fft.fftfreq(N)
    om = fine_bloch_omega(kk, th)
    sel = np.where(om <= omega_c + 1e-15)[0]
    sel = sel[np.argsort(om[sel], kind="stable")]
    c, s = np.cos(th), np.sin(th)
    V = np.zeros((2 * N, len(sel)), complex)
    for i, q in enumerate(sel):
        k = kk[q]
        M = np.array([[c, -1j * s], [-1j * s, c]]) @ np.diag([np.exp(-1j * k), np.exp(1j * k)])
        lam, U = np.linalg.eig(M)
        w = -np.angle(lam)
        u = U[:, np.argmax(w)]                      # positive branch (w in (0, pi))
        u = u / np.linalg.norm(u)
        ph = np.exp(1j * k * n) / np.sqrt(N)
        V[0::2, i] = u[0] * ph
        V[1::2, i] = u[1] * ph
    return V, om[sel]


def fine_U(N, th, Vtick=None):
    """sparse one-tick unitary: shift (R +1, L -1), coin, then the potential phase e^{-i V a} (Vtick = V a)"""
    c, s = np.cos(th), np.sin(th)
    rows, cols, vals = [], [], []
    for x in range(N):
        xr, xl = (x - 1) % N, (x + 1) % N           # R'(x) = R(x-1), L'(x) = L(x+1)
        ph = 1.0 if Vtick is None else np.exp(-1j * Vtick[x])
        rows += [2 * x, 2 * x, 2 * x + 1, 2 * x + 1]
        cols += [2 * xr, 2 * xl + 1, 2 * xr, 2 * xl + 1]
        vals += [c * ph, -1j * s * ph, -1j * s * ph, c * ph]
    return sps.csc_matrix((vals, (rows, cols)), shape=(2 * N, 2 * N))


def fine_bound_states(N, th, Vtick, nb, rest):
    """bound states of the fine walk with a potential: eigenvalues of U nearest e^{-i (rest - 0.6 * depth)}
    (shift-invert on the sparse unitary); returns omegas below the rest frequency, sorted"""
    U = fine_U(N, th, Vtick)
    depth = -Vtick.min() / th
    sig = np.exp(-1j * th * (1 - 0.5 * depth))
    lam = spla.eigs(U, k=4 * nb + 6, sigma=sig, return_eigenvectors=False, tol=1e-13, maxiter=20000)
    om = np.sort(-np.angle(lam))
    om = om[om < rest - 1e-9]
    return om[:nb]


# ------------------------------------------------------------------ blocking schemes
def fib_word(nsub):
    """cyclic Fibonacci word after nsub substitutions A->AB, B->A (start 'A'): list of 0 (A) / 1 (B)"""
    w = [0]
    for _ in range(nsub):
        out = []
        for c in w:
            out += [0, 1] if c == 0 else [0]
        w = out
    return w


def scheme_levels(kind, nlev, N0=None, nsub=None):
    """list of levels; level j has `types` (block types of its cells) and `groups` (indices of level j-1 cells forming
    each level-j cell); level 0 = fine cells. Returns (N0, levels)."""
    if kind in (2, 3):
        b = kind
        N = N0
        levels = [dict(types=np.zeros(N, int), groups=None, ncell=N)]
        n = N
        for j in range(nlev):
            if n % b:
                break
            m = n // b
            levels.append(dict(types=np.zeros(m, int), groups=[list(range(i * b, i * b + b)) for i in range(m)], ncell=m))
            n = m
        return N, levels
    # Fibonacci: level j cells = letters of W_{nsub - j}; deflation: coarse A <- fine 'A B', coarse B <- fine 'A'
    words = [fib_word(nsub - j) for j in range(min(nlev, nsub - 2) + 1)]
    levels = [dict(types=np.zeros(len(words[0]), int), groups=None, ncell=len(words[0]))]
    for j in range(1, len(words)):
        fine_w, coarse_w = words[j - 1], words[j]
        groups, p = [], 0
        for c in coarse_w:
            g = [p, p + 1] if c == 0 else [p]
            # check the substitution is consistent
            assert fine_w[p] == 0 and (c == 1 or fine_w[p + 1] == 1)
            groups.append(g)
            p += len(g)
        assert p == len(fine_w)
        levels.append(dict(types=np.array(coarse_w), groups=groups, ncell=len(coarse_w)))
    return len(words[0]), levels


# ------------------------------------------------------------------ automatic variables
def choose_variables(Vs, groups, types, Cprev, eps_c, cmax=3, C_force=None):
    """block density matrix per block type; returns W (n_coarse_dof x n_prev_dof), C, captured weights"""
    ntype = int(types.max()) + 1
    nb = len(groups)
    rho = {}
    for t in range(ntype):
        idxs = [g for g, ty in zip(groups, types) if ty == t]
        dim = len(idxs[0]) * Cprev
        R = np.zeros((dim, dim), complex)
        for g in idxs:
            ii = np.concatenate([np.arange(c * Cprev, c * Cprev + Cprev) for c in g])
            B = Vs[ii]
            R += B @ B.conj().T
        rho[t] = R / R.trace().real
    # one C for all types: smallest with captured >= 1 - eps_c for every type
    evs = {t: np.linalg.eigh(rho[t]) for t in rho}
    C = None
    for c in range(1, cmax + 1):
        if all(np.sum(np.sort(evs[t][0])[::-1][:c]) >= 1 - eps_c for t in evs if len(evs[t][0]) >= c):
            C = c
            break
    if C is None:
        C = cmax
    if C_force is not None:
        C = C_force
    caps = {t: float(np.sum(np.sort(evs[t][0])[::-1][:min(C, len(evs[t][0]))])) for t in evs}
    spectra = {t: np.sort(evs[t][0])[::-1][:4].tolist() for t in evs}
    vecs = {}
    for t, (e, U) in evs.items():
        o = np.argsort(e)[::-1][:C]
        Ut = U[:, o]
        # fix the gauge: make the summed amplitude of each kept vector real positive (cosmetic)
        for q in range(Ut.shape[1]):
            sm = Ut[:, q].sum()
            if abs(sm) > 1e-12:
                Ut[:, q] *= np.conj(sm) / abs(sm)
        vecs[t] = Ut
    nprev = Vs.shape[0]
    rows, cols, vals = [], [], []
    for I, (g, t) in enumerate(zip(groups, types)):
        ii = np.concatenate([np.arange(c * Cprev, c * Cprev + Cprev) for c in g])
        Ut = vecs[t]
        for q in range(min(C, Ut.shape[1])):
            rows += [I * C + q] * len(ii)
            cols += list(ii)
            vals += list(np.conj(Ut[:, q]))
    W = sps.csr_matrix((vals, (rows, cols)), shape=(nb * C, nprev))
    return W, C, caps, spectra, vecs


def hand_overlap(vecs, groups, types, Cprev, level_index):
    """fraction of the automatic block subspace inside the HAND descriptor subspace:
    level 1 with C = 1: block mean of the large component (R + L)/sqrt 2; otherwise: block means of every component
    of the level below (span of Cprev 'uniform' vectors). Returns tr(P_auto P_hand) / C per block type."""
    out = {}
    for t, Ut in vecs.items():
        g = [g for g, ty in zip(groups, types) if ty == t][0]
        n = len(g) * Cprev
        if level_index == 1 and Ut.shape[1] == 1:
            Hh = np.zeros((n, 1), complex); Hh[:, 0] = 1
        else:
            Hh = np.zeros((n, Cprev), complex)
            for c in range(Cprev):
                Hh[c::Cprev, c] = 1
        Hh, _ = np.linalg.qr(Hh)
        out[int(t)] = float(np.real(np.trace(Ut.conj().T @ Hh @ Hh.conj().T @ Ut)) / Ut.shape[1])
    return out


# ------------------------------------------------------------------ pattern-tied local generator
def patterns(types, w):
    n = len(types)
    pats = {}
    pid = np.zeros(n, int)
    for I in range(n):
        key = tuple(types[(I + d) % n] for d in range(-w, w + 1))
        if key not in pats:
            pats[key] = len(pats)
        pid[I] = pats[key]
    return pid, len(pats)


def fit_generator(Q, Om, types, C, w, ridge=1e-12, weights=None):
    """least squares for H = (M + M^dag)/2, M_{I,I+j} (j = 0..w) = tied by pattern(I): H Q = Q Omega.
    Q (n C, ns) complex, Om (ns,). Returns sparse H, number of real parameters."""
    n = len(types)
    pid, npat = patterns(types, w)
    ns = Q.shape[1]
    wt = np.ones(ns) if weights is None else weights
    Qw = Q * wt[None, :]
    QB = Qw.reshape(n, C, ns)
    nparc = npat * (w + 1) * C * C                    # complex params
    P = 2 * nparc
    AtA = np.zeros((P, P))
    Atb = np.zeros(P)
    Tgt = (QB * Om[None, None, :])

    def pidx(p, j, a, b):
        return ((p * (w + 1) + j) * C + a) * C + b
    # build A row block per cell I: rows (a, m) complex -> stacked real/imag
    for I in range(n):
        cols = []
        blocks = []                                    # each: (C, ns) complex contribution per unit real param
        for j in range(w + 1):
            p = pid[I]
            for a in range(C):
                for b in range(C):
                    base = pidx(p, j, a, b)
                    v = np.zeros((C, ns), complex); v[a] = 0.5 * QB[(I + j) % n, b]
                    cols.append(2 * base); blocks.append(v)
                    cols.append(2 * base + 1); blocks.append(1j * v)
            p2 = pid[(I - j) % n]
            for a in range(C):
                for b in range(C):
                    base = pidx(p2, j, a, b)
                    v = np.zeros((C, ns), complex); v[b] = 0.5 * QB[(I - j) % n, a]
                    cols.append(2 * base); blocks.append(v)              # conj(re) = re
                    cols.append(2 * base + 1); blocks.append(-1j * v)    # conj(i im) = -i im
        A = np.stack([bl.ravel() for bl in blocks], 1)          # (C ns, ncols) complex
        Ar = np.concatenate([A.real, A.imag], 0)
        y = Tgt[I].ravel()
        yr = np.concatenate([y.real, y.imag])
        cols = np.array(cols)
        # merge duplicate columns (j = 0 appears twice)
        uc, inv = np.unique(cols, return_inverse=True)
        Am = np.zeros((Ar.shape[0], len(uc)))
        np.add.at(Am.T, inv, Ar.T)
        AtA[np.ix_(uc, uc)] += Am.T @ Am
        Atb[uc] += Am.T @ yr
    lam = ridge * np.trace(AtA) / P
    x = np.linalg.solve(AtA + lam * np.eye(P), Atb)
    Mc = x[0::2] + 1j * x[1::2]
    rows, colsH, vals = [], [], []
    for I in range(n):
        p = pid[I]
        for j in range(w + 1):
            Jc = (I + j) % n
            for a in range(C):
                for b in range(C):
                    m = Mc[pidx(p, j, a, b)]
                    rows.append(I * C + a); colsH.append(Jc * C + b); vals.append(0.5 * m)
                    rows.append(Jc * C + b); colsH.append(I * C + a); vals.append(0.5 * np.conj(m))
    H = sps.csr_matrix((vals, (rows, colsH)), shape=(n * C, n * C))
    nreal = int(np.count_nonzero(np.abs(x) > 0))
    return H, P, npat


def lowdin(Q):
    S = Q.conj().T @ Q
    e, U = np.linalg.eigh(S)
    return Q @ (U * (1 / np.sqrt(np.maximum(e, 1e-300)))[None, :]) @ U.conj().T


# ------------------------------------------------------------------ the loop
def make_base(th, B0, wb=8, nk=20001, w_fast=0.03):
    """BASE level shared by all schemes: cells of B0 fine cells (Pauli regime, C = 1, the electron's positive band),
    kinetic generator per fine tick G(kappa) ~ eps(kappa / B0) = omega(k) - theta, learned as a real even stencil
    |j| <= wb by weighted least squares on a dense kappa grid (weight 1 for |kappa| <= pi/2, w_fast beyond: the same
    recipe as agg_zz_core). Returns stencil g (wb+1,), max relative error for |kappa| <= pi/2."""
    kap = np.linspace(-np.pi, np.pi, nk)
    tgt = fine_bloch_omega(kap / B0, th) - th
    wt = np.where(np.abs(kap) <= np.pi / 2, 1.0, w_fast)
    A = np.stack([np.ones_like(kap)] + [2 * np.cos(j * kap) for j in range(1, wb + 1)], 1)
    # transported exactly (phase-1 lesson): G(0) = 0 (rest phase removed) and the kinetic mass c2 of the fine rule
    d = 1e-4
    c2t = (fine_bloch_omega(d / B0, th) - th) / d ** 2
    Ceq = np.stack([np.concatenate([[1.0], 2 * np.ones(wb)]),                       # G(0)
                    np.concatenate([[0.0], -np.arange(1, wb + 1) ** 2.0])])          # kappa^2 coefficient
    deq = np.array([0.0, c2t])
    Aw = A * wt[:, None]
    K = np.block([[Aw.T @ Aw, Ceq.T], [Ceq, np.zeros((2, 2))]])
    rhs = np.concatenate([Aw.T @ (tgt * wt), deq])
    g = np.linalg.solve(K, rhs)[:wb + 1]
    fit = A @ g
    m = np.abs(kap) <= np.pi / 2
    err = float(np.max(np.abs(fit[m] - tgt[m])) / tgt[m].max())
    return g, err


def circulant(g, n):
    wb = len(g) - 1
    rows, cols, vals = [], [], []
    for I in range(n):
        for j in range(-wb, wb + 1):
            rows.append(I); cols.append((I + j) % n); vals.append(g[abs(j)])
    return sps.csr_matrix((vals, (rows, cols)), shape=(n, n))


def run_chain(scheme, H0, nlev=8, eps_level=1e-3, wmax=6, nobs=17, log=print, V0=None, nb=3, Bstop=80, cmax=2,
              eps_c=2e-2, nsub=None):
    """the loop (see module doc) on top of a base lattice with generator H0 (sparse Hermitian, C = 1, kinetic only:
    rest frequency removed, so slowness = the eigenvalue). Level sizes in base cells."""
    n0 = H0.shape[0]
    N, levels = scheme_levels(scheme, nlev, N0=n0, nsub=nsub)
    assert N == n0
    t0 = time.time()
    Bmean = [N / lv["ncell"] for lv in levels]
    e0, U0 = np.linalg.eigh(H0.toarray())
    om_exact = e0[:nobs]
    out = dict(scheme=str(scheme), N=N, levels=[], om_exact=om_exact.tolist())
    # kinetic scale of the base band at fraction f of the base Nyquist (from the base spectrum, by counting modes)
    def kin_at(Bcells, frac):
        # number of base modes with |kappa_base| <= frac * pi / Bcells  (two per kappa except 0)
        cnt = int(2 * np.floor(frac * N / (2 * Bcells)) + 1)
        return e0[min(cnt, N) - 1]
    if V0 is not None:
        eb0 = np.linalg.eigvalsh(H0.toarray() + sps.diags(V0).toarray())
        out["bound_base"] = np.sort(eb0)[:nb].tolist()
    Vs, oms = U0, e0
    Hprev = e0
    Cprev = 1
    Vpot = None if V0 is None else sps.diags(V0).tocsr()
    for j in range(1, len(levels)):
        lv = levels[j]
        tl = time.time()
        dchk = max(kin_at(Bmean[j], FRAC_CHECK), np.sort(Hprev)[MIN_CHECK - 1])   # >= MIN_CHECK modes in the check
        dhalf = max(kin_at(Bmean[j], FRAC_SLOW), dchk)
        ref = np.sort(Hprev[Hprev <= dchk * (1 + 1e-12)])
        best, tried = None, []
        chk = oms <= dchk * (1 + 1e-12)
        slow = oms <= dhalf * (1 + 1e-12)
        for Cc in range(1, cmax + 1):           # fewest numbers per block, then the most local rule, that pass
            nsec = lv["ncell"] * Cc
            if nsec > Vs.shape[1]:
                break
            W, C, caps, spectra, vecs = choose_variables(Vs[:, slow], lv["groups"], lv["types"], Cprev, eps_c, C_force=Cc)
            ho = hand_overlap(vecs, lv["groups"], lv["types"], Cprev, j)
            WV = W @ Vs
            capt = np.sum(np.abs(WV) ** 2, axis=0)
            cap_chk = float(np.mean(capt[chk]))
            T = np.sort(np.argsort(-capt, kind="stable")[:nsec])
            Q = lowdin(WV[:, T])
            wts = np.where(slow[T], 1.0, W_FAST) / (np.abs(oms[T]) + REL_FLOOR * dchk)   # relative least squares (as learn_h)
            missing = int(chk.sum() - chk[T].sum())
            for w in range(1, wmax + 1):
                if 2 * w + 1 > lv["ncell"]:
                    break
                H, P, npat = fit_generator(Q, oms[T], lv["types"], C, w, weights=wts)
                ev = np.linalg.eigvalsh(H.toarray())
                got = np.sort(ev)[:len(ref)]                 # the lowest modes, matched in order
                err = errc = float(np.max(np.abs(got - ref)) / dchk) if (not missing and len(ref)) else float("inf")
                ok = err <= eps_level and cap_chk >= 1 - eps_c
                tried.append(dict(C=C, w=w, err=err, err_common=errc, n_got=int(len(got)), n_ref=int(len(ref)),
                                  missing=missing, captured_check=cap_chk, params=P, patterns=npat, ok=ok))
                log(f"   [{scheme}] level {j} B={Bmean[j]:.3f} C={C} w={w}: err {err:.2e} ({len(got)}/{len(ref)}, miss {missing}, "
                    f"common {errc:.1e}, capt {cap_chk:.4f}) params {P} pat {npat} {'OK' if ok else ''}")
                if ok:
                    best = (w, H, ev, P, npat)
                    break
                if missing or cap_chk < 1 - eps_c:
                    break
            if best is not None:
                break
        rec = dict(level=j, B=Bmean[j], ncell=lv["ncell"], C=C, captured=caps, rho_top=spectra, hand_overlap=ho,
                   tried=tried)
        if best is None:
            rec["flag"] = "no admissible (C, w)"
            out["levels"].append(rec)
            log(f"   [{scheme}] level {j}: no admissible candidate -> stop")
            break
        w, H, ev, P, npat = best
        rec.update(w=w, params=P, patterns=npat, nnz=int(H.nnz), Hnorm=float(np.abs(ev).max()))
        mo = min(nobs, len(ev))
        rec["acc_err"] = float(np.max(np.abs(np.sort(ev)[:mo] - om_exact[:mo])) / (om_exact.max() - om_exact.min()))
        rec["leak"] = float(1 - np.mean(capt[slow]))
        if Vpot is not None:
            Vpot = (W @ Vpot @ W.conj().T).tocsr()
            evb = np.linalg.eigvalsh(H.toarray() + Vpot.toarray())
            eb = np.sort(evb)[:nb]
            rec["bound"] = eb.tolist()
            fb = np.array(out["bound_base"])
            rec["bound_rel_err"] = (np.abs(eb - fb) / np.abs(fb)).tolist()
        rec["sec"] = round(time.time() - tl, 1)
        out["levels"].append(rec)
        log(f"   [{scheme}] level {j} B={Bmean[j]:.3f}: C={C} w={w} acc {rec['acc_err']:.2e} leak {rec['leak']:.1e} "
            + (f"bound {np.round(rec['bound_rel_err'], 5)}" if Vpot is not None else "") + f" ({rec['sec']} s)")
        if j + 1 < len(levels) and levels[j + 1]["ncell"] >= NMIN_CELLS:
            e, U = np.linalg.eigh(H.toarray())
            Vs, oms, Hprev, Cprev = U, e, e, C
        else:
            break
    out["seconds"] = time.time() - t0
    out["Bmean"] = Bmean
    return out

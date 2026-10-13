"""aggmol_h2h2.py - phase 3 of study aggmol_: the force between two H2 molecules (4 electrons), the first step towards
liquid hydrogen (docs/AGGLOMERATO_DUE_ATOMI.md §13, §15). Pre-registration: lab/out/aggmol_h2h2_prereg.json.

Same operator as phase 2 (the main thread's rule as it is: p^2/2m_kin, m_kin = v3 m0; instantaneous Coulomb; singlet).
Rigid monomers at r = 1.4 a0 (Born-Oppenheimer + rigid rotor: a debt), centres 0 and R on z.

Method (a scale-down from a full 4-electron variational calculation, declared in the pre-registration):
  E_int = E1_HL + E2_ind + E2_disp
  E1_HL   Heitler-London energy of the antisymmetrised product of the two ECG monomers: a true 4-electron ECG
          calculation (all 24 permutations with the singlet x singlet spin factor); = electrostatics + first-order
          exchange (Pauli repulsion), to all orders in the overlap
  E2_ind  second-order induction, E2_disp second-order dispersion (no exchange), with the full Coulomb operator (no
          multipole expansion) from monomer pseudo-states (2-electron ECG, ground primitives + polarisation primitives),
          through 3-D Gaussian transition densities
  missing: exchange-dispersion, exchange-induction (both repulsive, small at the well), third and higher order.
C6 from the same pseudo-states (Casimir-Polder sum), cross-checked by the tail of E2_disp.
Jobs: mono | geom ORIENT | analyze ;   orientations H (parallel), X (crossed), T, L (collinear)."""
import os, sys, json, time, itertools
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_ecg as E

OUT = E.OUT
JS = os.path.join(OUT, "aggmol_h2h2.json")
NPZ = os.path.join(OUT, "aggmol_h2h2_mono.npz")       # not committed
RB = 1.4
K0 = int(os.environ.get("AGGMOL_K0", 18))
KAPPAS = [float(x) for x in os.environ.get("AGGMOL_KAPPAS", "0.6").split(",")]
WIDE = [float(x) for x in os.environ.get("AGGMOL_WIDE", "1").split(",")]
POLMODE = os.environ.get("AGGMOL_POLMODE", "e1")        # e1: shift electron 1; both: also both electrons together
RGRID = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0]
ORIENT = {"H": ((1, 0, 0), (1, 0, 0)), "X": ((1, 0, 0), (0, 1, 0)), "T": ((1, 0, 0), (0, 0, 1)),
          "L": ((0, 0, 1), (0, 0, 1))}
P12 = np.array([[0., 1.], [1., 0.]])


def load():
    return json.load(open(JS)) if os.path.exists(JS) else {}


def save(res):
    json.dump(res, open(JS, "w"), indent=1, default=float)


# ----------------------------------------------------------------------------------------------- 4-electron elements
def elems_n(A1, s1, A2, s2, nuclei, m=E.M_KIN):
    """as aggmol_ecg.elems, any number of electrons (all e-e pairs)"""
    N = A1.shape[-1]
    C = A1 + A2
    Ci = np.linalg.inv(C)
    det = np.linalg.det(C)
    v = A1 @ s1 + A2 @ s2
    R0 = Ci @ v
    c0 = np.einsum("...ix,...ij,...jx->...", s1, A1, s1) + np.einsum("...ix,...ij,...jx->...", s2, A2, s2)
    ex = np.einsum("...ix,...ij,...jx->...", v, Ci, v) - c0
    S = (np.pi ** N / det) ** 1.5 * np.exp(ex)
    Sig = 0.5 * Ci
    M = A1 @ A2
    d1, d2 = R0 - s1, R0 - s2
    Tk = 4 * (3 * np.einsum("...ij,...ji->...", M, Sig) + np.einsum("...ix,...ij,...jx->...", d1, M, d2)) / (2 * m)
    V = 0.0
    for i in range(N):
        for X in nuclei:
            V = V - E.coul(R0[..., i, :] - X, Sig[..., i, i])
        for j in range(i + 1, N):
            V = V + E.coul(R0[..., i, :] - R0[..., j, :], Sig[..., i, i] + Sig[..., j, j] - 2 * Sig[..., i, j])
    return S, Tk * S, V * S


def spin_factors():
    """c_P = sign(P) <chi|P chi>, chi = singlet(12) x singlet(34); P as a permutation matrix acting on electron labels"""
    up, dn = np.array([1., 0.]), np.array([0., 1.])
    th = (np.kron(up, dn) - np.kron(dn, up)) / np.sqrt(2)
    chi = np.kron(th, th).reshape((2,) * 4)
    out = []
    for p in itertools.permutations(range(4)):
        Pm = np.eye(4)[list(p)]
        sign = np.linalg.det(Pm)
        c = sign * np.sum(chi * np.transpose(chi, p))
        if abs(c) > 1e-12:
            out.append((p, Pm, float(c)))
    return out


# ----------------------------------------------------------------------------------------------- monomer
def mono_mats(A, s, nuc, A2=None, s2=None):
    """S, H with the (1 + P12) symmetriser (singlet), no inversion projection"""
    A2 = A if A2 is None else A2
    s2 = s if s2 is None else s2
    S = 0.0; H = 0.0
    for P in (np.eye(2), P12):
        Ap = np.einsum("ji,kjl,lm->kim", P, A2, P)
        sp = np.einsum("ji,kjx->kix", P, s2)
        a, t, v = E.elems(A[:, None], s[:, None], Ap[None], sp[None], nuc)
        S = S + a; H = H + t + v
    return S, H


OPS = {"X": [0], "Z": [2], "XZ": [0, 2], "ZZ": [2], "XY": [0, 1], "ZZZ": [2], "XZZ": [0, 2], "XYZ": [0, 1, 2],
       "X3": [0, 1]}
ROTS = {"Y": ("X", 90), "YZ": ("XZ", 90), "X2Y2": ("XY", 45), "YZZ": ("XZZ", 90), "ZX2Y2": ("XYZ", 45),
        "Y3": ("X3", 30)}
KPOL = int(os.environ.get("AGGMOL_KPOL", 15))


def op_value(g, op):
    c, v = g["c"], g["v"]
    return {"X": c[:, 0], "Y": c[:, 1], "Z": c[:, 2], "XZ": c[:, 0] * c[:, 2], "YZ": c[:, 1] * c[:, 2],
            "XY": c[:, 0] * c[:, 1], "ZZ": c[:, 2] ** 2 - 0.5 * (c[:, 0] ** 2 + c[:, 1] ** 2),
            "X2Y2": 0.5 * (c[:, 0] ** 2 - c[:, 1] ** 2),
            # octupoles (Gaussian moments: E[z^3] = cz^3 + 3 v cz, E[z x^2] = cz cx^2 + v cz, ...)
            "ZZZ": c[:, 2] ** 3 + 3 * v * c[:, 2] - 1.5 * (c[:, 2] * (c[:, 0] ** 2 + c[:, 1] ** 2) + 2 * v * c[:, 2]),
            "XZZ": c[:, 0] * (4 * (c[:, 2] ** 2 + v) - (c[:, 0] ** 2 + 3 * v) - (c[:, 1] ** 2 + v)),
            "XYZ": c[:, 0] * c[:, 1] * c[:, 2],
            "X3": c[:, 0] ** 3 + 3 * v * c[:, 0] - 3 * c[:, 0] * (c[:, 1] ** 2 + v)}[op]


def pseudo(S, H, c0, EA):
    """pseudo-states in the space orthogonal to psi0 (psi0 fixed); returns coefficients D and excitation energies"""
    Q = np.eye(len(S)) - np.outer(c0, S @ c0)
    Sq = Q.T @ S @ Q; Hq = Q.T @ H @ Q
    ws, Vs = np.linalg.eigh(Sq)
    keep = ws > 1e-10 * ws.max()
    X = Vs[:, keep] / np.sqrt(ws[keep])
    em, Wm = np.linalg.eigh(X.T @ Hq @ X)
    return Q @ X @ Wm, em - EA


def rotz(s, ang):
    c, si = np.cos(ang), np.sin(ang)
    Rm = np.array([[c, -si, 0], [si, c, 0], [0, 0, 1]])
    return s @ Rm.T


def job_mono():
    """ground state (K0 ECG, with inversion images) + polarisation functions chosen by the stochastic variational
    method on the Hylleraas functional of each operator (dipole x, z; quadrupole xz, zz, xy); y, yz, x^2-y^2 by rotation
    about the molecular axis. The pseudo-states of the whole set give the dynamic polarisabilities, C6, and E2."""
    t0 = time.time()
    res = load()
    rng = np.random.default_rng(int(os.environ.get("AGGMOL_SEED", 3)))
    nuc = [np.array((0, 0, RB / 2)), np.array((0, 0, -RB / 2))]
    e, B, t = E.build(2, RB, K0, ntrial=40, nref=int(os.environ.get("AGGMOL_NREF", 6)), budget=240)
    Ag = np.concatenate([B.A, B.A]); sg = np.concatenate([B.s, -B.s])
    ng = len(Ag)
    S, H = mono_mats(Ag, sg, nuc)
    w, U = sla.eigh(H, S)
    c = U[:, 0] / np.sqrt(U[:, 0] @ S @ U[:, 0])
    EA = float(c @ H @ c)
    gen = E.make_gen(2, RB)
    sets = {}
    alph = {}
    for op, axes in OPS.items():
        A_, s_ = Ag.copy(), sg.copy()
        S_, H_ = S.copy(), H.copy()
        best_a = 0.0
        for step in range(KPOL):
            At, st = gen(rng, 40)
            for k in range(40):
                if k % 2 == 0:                       # half: a ground primitive, widened, with shifted electrons
                    j = rng.integers(ng)
                    At[k] = Ag[j] * np.exp(rng.uniform(np.log(0.15), np.log(1.5))); st[k] = sg[j]
                a = np.diag(At[k])
                L = np.exp(rng.uniform(np.log(0.2), np.log(2.5)))
                for ee in range(2):
                    if ee == 0 or rng.random() < 0.5:
                        for ax in axes:
                            st[k, ee, ax] += rng.normal() * L / np.sqrt(a[ee])
            Sr, Hr = mono_mats(At, st, nuc, A_, s_)
            Sd, Hd = mono_mats(At, st, nuc)
            bi, bS, bH = -1, None, None
            for k in range(40):
                n = len(A_)
                S2 = np.zeros((n + 1, n + 1)); H2 = np.zeros((n + 1, n + 1))
                S2[:n, :n] = S_; H2[:n, :n] = H_
                S2[n, :n] = S2[:n, n] = Sr[k]; H2[n, :n] = H2[:n, n] = Hr[k]
                S2[n, n] = Sd[k, k]; H2[n, n] = Hd[k, k]
                nrm = np.sqrt(np.diag(S2))
                if np.max(np.abs(S2[n, :n]) / (nrm[n] * nrm[:n])) > 0.99:
                    continue
                c0 = np.concatenate([c, np.zeros(n + 1 - ng)])
                try:
                    D, dE = pseudo(S2, H2, c0, EA)
                except np.linalg.LinAlgError:
                    continue
                if dE.min() < 0.3:
                    continue
                g = transition_gaussians(np.concatenate([A_, At[k:k + 1]]), np.concatenate([s_, st[k:k + 1]]), c0, D, ng, nuc)
                b = g["T"].T @ op_value(g, op)
                al = 2 * np.sum(b ** 2 / dE)
                if al > best_a:
                    best_a, bi, bS, bH = al, k, S2, H2
            if bi < 0:
                continue
            A_ = np.concatenate([A_, At[bi:bi + 1]]); s_ = np.concatenate([s_, st[bi:bi + 1]])
            S_, H_ = bS, bH
        sets[op] = (A_[ng:], s_[ng:])
        alph[op] = float(best_a)
        print(f"  {op}: {len(A_) - ng} functions, alpha_{op} = {best_a:.5f} ({time.time() - t0:.0f} s)", flush=True)
    # y, yz, x^2 - y^2 by rotation about z
    for k, (src, ang) in ROTS.items():
        sets[k] = (sets[src][0], rotz(sets[src][1], np.deg2rad(ang)))
    Aall = np.concatenate([Ag] + [sets[k][0] for k in sets]); sall = np.concatenate([sg] + [sets[k][1] for k in sets])
    S, H = mono_mats(Aall, sall, nuc)
    c0 = np.concatenate([c, np.zeros(len(Aall) - ng)])
    D, dE = pseudo(S, H, c0, EA)
    np.savez(NPZ, A=Aall, s=sall, c0=c0, D=D, dE=dE, EA=EA, nground=ng)
    g = transition_gaussians(Aall, sall, c0, np.concatenate([c0[:, None], D], 1), ng, nuc)
    mom = {op: (g["T"].T @ op_value(g, op))[1:] for op in ("X", "Y", "Z", "XZ", "ZZ", "XY")}
    alpha = {op: float(2 * np.sum(mom[op] ** 2 / dE)) for op in mom}
    F = (mom["X"] ** 2 + mom["Y"] ** 2 + mom["Z"] ** 2) / 3
    C6 = float(6 * np.sum(F[:, None] * F[None] / (dE[:, None] + dE[None])))
    # anisotropic C6 for the orientations: C6(orient) = (1/pi) int sum over the dipole tensor (computed in analyze)
    res["mono"] = dict(K0=K0, KPOL=KPOL, E_A=EA, E_A_build=e, n_prim=len(Aall), n_states=len(dE),
                       dE_min=float(dE.min()), alpha_opt=alph, alpha=alpha,
                       alpha_iso=(alpha["X"] + alpha["Y"] + alpha["Z"]) / 3, C6_CP=C6,
                       ground_quadrupole=float(np.sum(g["T"][:, 0] * op_value(g, "ZZ"))) , sec=time.time() - t0)
    save(res)
    print(res["mono"], flush=True)


def transition_gaussians(A, s, c0, Cm, nground, nuc):
    """one-body transition densities rho_0m(r) = <Psi0| sum_e delta(r - r_e) |Psi_m> as sums of normalised 3-D
    Gaussians: centres c (G, 3), variances v (G,), weights T (G, M) (state coefficients folded in)"""
    Ai, si = A[:nground], s[:nground]
    cs, vs, Ts = [], [], []
    for P in (np.eye(2), P12):
        Ap = np.einsum("ji,kjl,lm->kim", P, A, P)
        sp = np.einsum("ji,kjx->kix", P, s)
        C = Ai[:, None] + Ap[None]
        Ci = np.linalg.inv(C); det = np.linalg.det(C)
        v = Ai[:, None] @ si[:, None] + Ap[None] @ sp[None]
        R0 = Ci @ v
        c0_ = np.einsum("...ix,...ij,...jx->...", si[:, None], Ai[:, None], si[:, None]) + \
            np.einsum("...ix,...ij,...jx->...", sp[None], Ap[None], sp[None])
        S = (np.pi ** 2 / det) ** 1.5 * np.exp(np.einsum("...ix,...ij,...jx->...", v, Ci, v) - c0_)   # (ng, K)
        W = (c0[:nground, None] * S)[:, :, None] * Cm[None]                                          # (ng, K, M)
        for e in range(2):
            cs.append(R0[:, :, e, :].reshape(-1, 3)); vs.append(0.5 * Ci[:, :, e, e].reshape(-1))
            Ts.append(W.reshape(-1, Cm.shape[1]))
    return dict(c=np.concatenate(cs), v=np.concatenate(vs), T=np.concatenate(Ts))


def fast_coul(ca, cb, va, vb):
    """Coulomb energy between unit 3-D Gaussian charges (centres ca, cb; variances per axis va, vb): all pairs"""
    from scipy.special import erf
    d2 = np.sum(ca ** 2, 1)[:, None] + np.sum(cb ** 2, 1)[None] - 2 * ca @ cb.T
    np.maximum(d2, 0, out=d2)
    s2 = va[:, None] + vb[None]
    x2 = d2 / (2 * s2)
    m = np.sqrt(d2)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = erf(np.sqrt(x2)) / m
    small = x2 < 1e-8
    if np.any(small):
        f[small] = (np.sqrt(2 / np.pi) / np.sqrt(s2) * (1 - x2 / 3))[small]
    return f


def place(g, axis, centre):
    """rotate the monomer frame (axis z) onto `axis` and translate"""
    z = np.array(axis, float); z /= np.linalg.norm(z)
    a = np.array([1., 0, 0]) if abs(z[0]) < 0.9 else np.array([0., 1, 0])
    x = a - z * (a @ z); x /= np.linalg.norm(x); y = np.cross(z, x)
    Rm = np.stack([x, y, z], 1)                       # columns: images of the frame axes
    return Rm, np.array(centre, float)


# ----------------------------------------------------------------------------------------------- one geometry
DEMAX = float(os.environ.get("AGGMOL_DEMAX", 1e9))


def geometry(orient, R, mono, spins, chunk=300, do_exch=True, check_direct=False, dEmax=None):
    A, s, c0, D, dE, EA, ng = (mono[k] for k in ("A", "s", "c0", "D", "dE", "EA", "nground"))
    EA = float(EA); ng = int(ng)
    dEmax = DEMAX if dEmax is None else dEmax
    keep = dE < dEmax                                    # pseudo-states kept (declared cut)
    D, dE = D[:, keep], dE[keep]
    ua, ub = ORIENT[orient]
    RA, cA = place(None, ua, (0, 0, 0)); RBm, cB = place(None, ub, (0, 0, R))
    nucA = [RA @ np.array((0, 0, sg * RB / 2)) + cA for sg in (1, -1)]
    nucB = [RBm @ np.array((0, 0, sg * RB / 2)) + cB for sg in (1, -1)]
    nuc_mono = [np.array((0, 0, RB / 2)), np.array((0, 0, -RB / 2))]
    g = transition_gaussians(A, s, c0, np.concatenate([c0[:, None], D], 1), ng, nuc_mono)
    gA = dict(c=g["c"] @ RA.T + cA, v=g["v"], T=g["T"]); gB = dict(c=g["c"] @ RBm.T + cB, v=g["v"], T=g["T"])
    # Coulomb between the transition densities: W = TA^T G TB, built in row chunks
    M = g["T"].shape[1]
    W = np.zeros((M, M))
    for i in range(0, len(gA["v"]), chunk):
        G = fast_coul(gA["c"][i:i + chunk], gB["c"], gA["v"][i:i + chunk], gB["v"])
        W += gA["T"][i:i + chunk].T @ (G @ gB["T"])
    # nuclear potentials on the transition densities
    nA = np.zeros(M); nB = np.zeros(M)
    for X in nucB:
        nA -= gA["T"].T @ E.coul(gA["c"] - X, gA["v"])      # <0|V_nucB|m> on A
    for X in nucA:
        nB -= gB["T"].T @ E.coul(gB["c"] - X, gB["v"])
    NN = sum(1 / np.linalg.norm(a - b) for a in nucA for b in nucB)
    Eelst = W[0, 0] + nA[0] + nB[0] + NN
    # induction: field of B (electrons + nuclei) on the excitations of A, and vice versa
    wA = W[1:, 0] + nA[1:]; wB = W[0, 1:] + nB[1:]
    Eind = -float(np.sum(wA ** 2 / dE) + np.sum(wB ** 2 / dE))
    Wd = W[1:, 1:]
    Edisp = -float(np.sum(Wd ** 2 / (dE[:, None] + dE[None])))
    out = dict(orient=orient, R=R, E_elst=float(Eelst), E_ind=Eind, E_disp=Edisp)
    if do_exch:
        out.update(heitler_london(A[:ng], s[:ng], c0[:ng], EA, RA, cA, RBm, cB, nucA + nucB, spins, Eelst,
                                  check_direct=check_direct))
    return out


def heitler_london(Ag, sg, c, EA, RA, cA, RBm, cB, nuc, spins, Eelst, chunk=40000, check_direct=False):
    """E1_HL - 2 E_A = [E_elst + X_H - 2E_A X_S] / (1 + X_S) (X_H with the nuclear repulsion), X over the permutations that move electrons between
    the molecules (4-electron ECG). Sum over all pairs (a, b) of product primitives; P and P^-1 give the same sum
    (counted once, twice the weight); for an involution only a <= b is computed. Optional check of the direct part
    against E_elst + 2 E_A."""
    n = len(Ag)
    sa = sg @ RA.T + cA; sb = sg @ RBm.T + cB
    I, J = np.meshgrid(np.arange(n), np.arange(n), indexing="ij"); I = I.ravel(); J = J.ravel()
    Ap = np.zeros((n * n, 4, 4)); Ap[:, :2, :2] = Ag[I]; Ap[:, 2:, 2:] = Ag[J]
    sp = np.concatenate([sa[I], sb[J]], 1)
    cp = c[I] * c[J]
    K = n * n
    direct = {(0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2)}
    done = set()
    tri_a, tri_b = np.triu_indices(K)
    out = dict(XH=0.0, XS=0.0, DH=0.0, DS=0.0)
    for p, Pm, cP in spins:
        if p in done or (p in direct and not check_direct):
            continue
        inv = tuple(int(i) for i in np.argsort(p))
        done.add(p); done.add(inv)
        A2 = np.einsum("ji,kjl,lm->kim", Pm, Ap, Pm)
        s2 = np.einsum("ji,kjx->kix", Pm, sp)
        h = o = 0.0
        if inv == p:                                   # involution: symmetric in (a, b)
            for i0 in range(0, len(tri_a), chunk):
                a, b = tri_a[i0:i0 + chunk], tri_b[i0:i0 + chunk]
                S_, T_, V_ = elems_n(Ap[a], sp[a], A2[b], s2[b], nuc)
                w = cp[a] * cp[b] * np.where(a == b, 1.0, 2.0)
                h += np.sum(w * (T_ + V_ - (0 if p in direct else 2 * EA) * S_)); o += np.sum(w * S_)
            mult = 1.0
        else:
            step = max(1, chunk // K)
            for i0 in range(0, K, step):
                i1 = min(K, i0 + step)
                S_, T_, V_ = elems_n(Ap[i0:i1, None], sp[i0:i1, None], A2[None], s2[None], nuc)
                w = cp[i0:i1, None] * cp[None]
                h += np.sum(w * (T_ + V_ - 2 * EA * S_)); o += np.sum(w * S_)
            mult = 2.0
        if p in direct:
            out["DH"] += cP * h; out["DS"] += cP * o
        else:
            out["XH"] += mult * cP * h; out["XS"] += mult * cP * o
    XH, XS = out["XH"], out["XS"]              # XH already holds <u|(H_el - 2 E_A) P|u>
    NN = sum(1 / np.linalg.norm(a - b) for a in nuc[:2] for b in nuc[2:])   # inter-molecular nuclear repulsion
    E1 = (Eelst + XH + NN * XS) / (1 + XS)
    r = dict(E1_HL=float(E1), E_exch1=float(E1 - Eelst), X_S=float(XS))
    if check_direct:
        r.update(direct_check=float(out["DH"] + NN * out["DS"] - 2 * EA * out["DS"] - Eelst), direct_norm=float(out["DS"]))
    return r


def job_geom(orient, budget=360):
    t0 = time.time()
    res = load()
    mono = dict(np.load(NPZ))
    spins = spin_factors()
    cur = res.setdefault("geom", {})
    for R in RGRID:
        key = f"{orient}_{R:.1f}"
        if key in cur:
            continue
        if time.time() - t0 > budget:
            print("checkpoint: run again", flush=True); return
        t1 = time.time()
        out = geometry(orient, R, mono, spins, check_direct=(orient == "H" and R == 6.5))
        out["sec"] = time.time() - t1
        out["E_int"] = out["E1_HL"] + out["E_ind"] + out["E_disp"]
        res = load(); res.setdefault("geom", {})[key] = out; save(res); cur = res["geom"]
        print(key, {k: (f"{v:.4e}" if isinstance(v, float) else v) for k, v in out.items()}, flush=True)
    print(f"ORIENT {orient} DONE", flush=True)


def job_analyze():
    from scipy.interpolate import CubicSpline
    res = load()
    G = res["geom"]
    Rs = np.array(RGRID)
    comp = ("E_elst", "E_exch1", "E_ind", "E_disp", "E_int")
    V = {o: {c: np.array([G[f"{o}_{R:.1f}"][c] for R in Rs]) for c in comp} for o in ORIENT}
    iso = {c: (2 * V["H"][c] + V["L"][c] + 4 * V["T"][c] + 2 * V["X"][c]) / 9 for c in comp}
    x = np.linspace(Rs[0], Rs[-1], 7001)

    def well(y):
        cs = CubicSpline(Rs, y); v = cs(x); i = int(np.argmin(v))
        return float(x[i]), float(v[i])
    wells = {o: well(V[o]["E_int"]) for o in ORIENT}
    Rm, Vm = well(iso["E_int"])
    m = res["mono"]
    # tail cross-check: -E_disp R^6 = C6 + C8 / R^2 on R = 10, 12
    y = -iso["E_disp"] * Rs ** 6
    sel = Rs >= 10
    C8, C6t = np.polyfit(Rs[sel] ** -2, y[sel], 1)
    pos = all(G[k]["E_exch1"] > 0 for k in G if G[k]["R"] <= 10)
    M0 = pos and abs(G["H_6.5"].get("direct_check", 1.0)) <= 1e-10
    verdict = {
        "M0": dict(direct_check=G["H_6.5"].get("direct_check"), exch_positive_R_le_10=pos, PASS=bool(M0)),
        "W1": dict(C6=m["C6_CP"], rel=m["C6_CP"] / 12.0 - 1, PASS=bool(abs(m["C6_CP"] / 12.0 - 1) <= 0.10)),
        "W2": dict(depth_mEh=-1e3 * Vm, depth_K=-Vm * 315775.02480407, rel=-Vm / 1.10e-4 - 1,
                   PASS=bool(abs(-Vm / 1.10e-4 - 1) <= 0.20)),
        "W3": dict(R_min=Rm, R_min_A=Rm * 0.529177210903, rel=Rm / 6.45 - 1, PASS=bool(abs(Rm / 6.45 - 1) <= 0.05)),
        "W4": {o: dict(R_min=wells[o][0], depth_mEh=-1e3 * wells[o][1]) for o in ORIENT},
    }
    res["iso"] = {c: iso[c].tolist() for c in comp}
    res["tail"] = dict(C6_fit=float(C6t), C8_fit=float(C8))
    res["verdict"] = verdict
    save(res)
    print(json.dumps(verdict, indent=1, default=float)); print(res["tail"])


if __name__ == "__main__":
    job = sys.argv[1]
    if job == "mono":
        job_mono()
    elif job == "geom":
        job_geom(sys.argv[2])
    elif job == "analyze":
        job_analyze()

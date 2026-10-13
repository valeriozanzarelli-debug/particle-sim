"""vuoto_f3_core.py - phase 3 of 'Vuoto che polarizza' (docs/VUOTO_POLARIZZA.md, section 14).

Question: the dressed electron of the rule has an O(alpha) velocity shift that does not go away with the lattice
spacing (phase 2, G3/G4).  Is it the Collins-Perez-Sudarsky-Urrutia-Vucetich effect (a Lorentz-breaking cutoff feeds
O(alpha) Lorentz violation to low-energy particles through loops)?  Which part is anisotropic, which isotropic, and
does a cyclic axis order (PROPOSAL, not the rule) remove the anisotropy?

Engine: the one-loop self-energy computed EXACTLY in the discrete time of the walk (no effective Hamiltonian):
the period is a list of operations (hop along axis j for the electron AND the light, coin for the electron only);
the light couples to the electron at each hop by the Peierls phase of the link (first order: -e dT_j/dk_j at the
midpoint momentum), the photon is the light walk itself (positive-frequency Riemann-Silberstein mode, helicity +
from F = E + iB at q, helicity - from G = E - iB, i.e. the walk at -q), its A-polarization at hop o is the walk's
state there (interpolation tau between before/after its own hop: DECLARED choice, tau = 1/2; tau = 0, 1 = systematic).
Second-order eigenvalue perturbation of the period operator F = O_N ... O_1 (exact for a unitary Floquet operator):
   delta mu_s = F2_ss + sum_m F1_sm F1_ms / (mu_s - mu_m),   delta phi_s = i delta mu_s / mu_s,
with hole theory (sea electrons b < 0 at p - q emit a photon -q and would land in s: that term is subtracted).
The seagull (A^2) term of the hop is the same for every electron state (dT^2/dk^2 = -T): dropped (cancels in dm, dv).
Coulomb part: instantaneous kernel 4 pi alpha/q^2 (DEBT, input as in phases 1-2), applied where the rule applies
the potential (after the coin, at the end of every tick): first-order exchange with the symmetric sea at each of
those times.  Species weight w(|q|) (K1 = 0.6, K2 = 1.05) = UV cutoff, as in phase 2 (declared).

Rules: 'R' = the frozen rule (ticks xyz, zyx; period 2).  'C' = PROPOSAL (ticks xyz, zyx, yzx, xzy, zxy, yxz;
period 6: the three palindromic blocks xyzzyx, yzxxzy, zxyyxz).
Units: lattice a = 1, tick = 1, m = theta.  Energies per tick.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from zz_core import ALPHA, BASIS_HALF, BASIS_ONE
from vuoto_core import ueig, species_weight, PI
from vuoto_f2_core import sphere_rule

RULES = {"R": [(0, 1, 2), (2, 1, 0)],
         "C": [(0, 1, 2), (2, 1, 0), (1, 2, 0), (0, 2, 1), (2, 0, 1), (1, 0, 2)]}
EYE4 = np.eye(4, dtype=complex)
SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex), np.diag([1.0 + 0j, -1])]


ROT = {"n": 0}          # test only: start the period after ROT["n"] operations (must not change any quasi-energy)


def op_list(rule):
    """[('h', j) ... ('c', None) ...] for one period"""
    out = []
    for order in RULES[rule]:
        for j in order:
            out.append(("h", j))
        out.append(("c", None))
    r = ROT["n"] % len(out)
    return out[r:] + out[:r]


def e_hop(K, j):
    """electron hop along j at momenta K (N,3): blockdiag(exp(+i k_j sigma_j) [L, h=-1], exp(-i k_j sigma_j) [R]).
    returns T (N,4,4) and dT/dk_j (N,4,4)"""
    k = K[:, j]
    c, s = np.cos(k), np.sin(k)
    S = SIG[j]
    TL = c[:, None, None] * np.eye(2) + 1j * s[:, None, None] * S
    TR = c[:, None, None] * np.eye(2) - 1j * s[:, None, None] * S
    T = np.zeros((len(K), 4, 4), complex)
    T[:, :2, :2], T[:, 2:, 2:] = TL, TR
    D = np.zeros_like(T)
    D[:, :2, :2] = 1j * S @ TL
    D[:, 2:, 2:] = -1j * S @ TR
    return T, D


def coin(th):
    ct, st = np.cos(th), np.sin(th)
    return np.block([[ct * np.eye(2), -1j * st * np.eye(2)], [-1j * st * np.eye(2), ct * np.eye(2)]]).astype(complex)


def g_hop(Q, j):
    """light hop along j: exp(-i q_j S_j) (spin 1), (N,3,3) real"""
    lam, P = BASIS_ONE[j]
    ph = np.exp(-1j * lam[None, :] * Q[:, j:j + 1])
    return np.einsum("ab,nb,cb->nac", P, ph, P.conj())


def W_e(K, th, rule, start=0):
    """electron period operator starting after op 'start' (cyclic representation)"""
    ops = op_list(rule)
    ops = ops[start:] + ops[:start]
    W = np.broadcast_to(EYE4, (len(K), 4, 4)).copy()
    C = coin(th)
    for kind, j in ops:
        W = (e_hop(K, j)[0] if kind == "h" else C[None]) @ W
    return W


def W_g(Q, rule):
    W = np.broadcast_to(np.eye(3, dtype=complex), (len(Q), 3, 3)).copy()
    for kind, j in op_list(rule):
        if kind == "h":
            W = g_hop(Q, j) @ W
    return W


def e_bands(K, th, rule):
    """Omega per PERIOD in (-pi, pi), eigenvectors (N,4,4)"""
    Om, lam, V = ueig(W_e(K, th, rule))
    return Om, V


def n_ticks(rule):
    return len(RULES[rule])


# ================================================================== exact discrete-time one-loop self-energy
def _photon_chain(Q, rule, nt):
    """positive-frequency light mode at Q along the period: g (N, nops+1, 3) chained through the light hops (coins
    = identity for light), period phase Omega (N,) in (0, pi), frequency per tick omega (N,) = Omega / nt"""
    W = W_g(Q, rule)
    lam, V = np.linalg.eig(W)
    Om = -np.angle(lam)                                             # e^{-i Om}
    i = np.argmax(Om, axis=1)                                       # largest positive Omega = helicity + of F at Q
    f = np.take_along_axis(V, i[:, None, None], 2)[:, :, 0]
    f /= np.linalg.norm(f, axis=1)[:, None]
    Omh = np.take_along_axis(Om, i[:, None], 1)[:, 0]
    ops = op_list(rule)
    g = [f]
    for kind, j in ops:
        g.append(np.einsum("nab,nb->na", g_hop(Q, j), g[-1]) if kind == "h" else g[-1])
    return np.stack(g, 1), Omh, Omh / nt


def _e_chain(K, th, rule, V0):
    """electron vectors V0 (N,4,M) chained through the period: (N, nops+1, 4, M)"""
    C = coin(th)
    out = [V0]
    for kind, j in op_list(rule):
        out.append((e_hop(K, j)[0] if kind == "h" else C[None]) @ out[-1])
    return np.stack(out, 1)


def self_energy_exact(p, th, rule="R", K12=(0.6, 1.05), nr=80, nct=16, nph=24, kmin_fac=1e-4, chunk=6000,
                      tau=0.5, coulomb=True, hole=True):
    """quasi-energy shifts per tick of the 4 electron states at momentum p (3,), exact discrete time, one loop.
    Returns dict with E0 (4,), dT (4,) transverse, dC (4,) Coulomb, imag (4,) and the band order (E0 ascending)."""
    p = np.asarray(p, float)
    nt = n_ticks(rule)
    ops = op_list(rule)
    w = species_weight(*K12)
    Qg, WW = sphere_rule(nr, nct, nph, kmin_fac * th, K12[1])
    WW = WW * w(Qg) / (2 * PI) ** 3
    Om_s, Vs = e_bands(p[None], th, rule)
    o = np.argsort(Om_s[0]); Om_s = Om_s[0][o]; Vs = Vs[0][:, o]
    vs = _e_chain(p[None], th, rule, Vs[None])[0]                    # (nops+1, 4, 4)
    mu_s = np.exp(-1j * Om_s)
    sgn_s = np.sign(Om_s)
    accT = np.zeros(4, complex)
    accC = np.zeros(4)
    hops = [(o_ + 1, j) for o_, (kind, j) in enumerate(ops) if kind == "h"]
    tick_end = [o_ + 1 for o_, (kind, j) in enumerate(ops) if kind == "c"]
    e2 = 4 * PI * ALPHA
    hole_flag = hole
    for c0 in range(0, len(Qg), chunk):
        Q, Wc = Qg[c0:c0 + chunk], WW[c0:c0 + chunk]
        N = len(Q)
        K = p[None] - Q                                              # intermediate electron
        Kmid = p[None] - 0.5 * Q
        Om_b, Vb = e_bands(K, th, rule)
        sb = np.sign(Om_b)
        Vc = _e_chain(K, th, rule, Vb)                               # (N, nops+1, 4, 4)
        gp, Omp, omp = _photon_chain(Q, rule, nt)                    # photon at +q (helicity +)
        gm, Omm, omm = _photon_chain(-Q, rule, nt)                   # walk at -q (helicity - of the photon at +q)
        # emission s -> b (at hop o):  e[o] = <V_{o,b}| D_o |v_{s,o-1}>,  absorption b -> s: a[o] = <v_{s,o}| D_o |V_{o-1,b}>
        e_list, a_list, js = [], [], []
        for (oi, j) in hops:
            D = e_hop(Kmid, j)[1]
            e_list.append(np.einsum("nab,nac,cs->nbs", Vc[:, oi].conj(), D, vs[oi - 1]))
            a_list.append(np.einsum("as,nac,ncb->nbs", vs[oi].conj(), D, Vc[:, oi - 1]))
            js.append((oi, j))
        e_arr = np.stack(e_list, 1); a_arr = np.stack(a_list, 1)     # (N, nh, 4b, 4s)

        def gt(g):                                                   # (N, nh) photon amplitude component at the hop
            return np.stack([(1 - tau) * g[:, oi - 1, j] + tau * g[:, oi, j] for (oi, j) in js], 1)
        # photon at +q: helicities from gp (+) and gm (-); photon at -q: helicities from gm (+) and gp (-)
        for (gA, OmA, omA), (gB, OmB, omB) in (((gp, Omp, omp), (gm, Omm, omm)), ((gm, Omm, omm), (gp, Omp, omp))):
            # gA: mode of the photon emitted by s (momentum +q) ; gB: the same helicity seen from -q (used by the sea)
            cA, cB = gt(gA), gt(gB)
            norm = e2 / (2 * omA)                                     # omA = omB
            # ---- single particle: s -> b + photon(+q, mode A) -> s
            Eamp = np.einsum("nhbs,nh->nbs", e_arr, cA.conj())
            Aamp = np.einsum("nhbs,nh->nbs", a_arr, cA)
            cum = np.cumsum(e_arr * cA.conj()[:, :, None, None], 1)  # sum_{o <= o'} emission
            cum = np.concatenate([np.zeros_like(cum[:, :1]), cum[:, :-1]], 1)   # strictly o < o'
            F2 = np.einsum("nhbs,nh,nhbs->nbs", a_arr, cA, cum)
            mu_m = np.exp(-1j * (Om_b + OmA[:, None]))               # (N, 4b)
            T = Aamp * Eamp * (mu_m[:, :, None] / (mu_s[None, None, :] - mu_m[:, :, None]))
            single = 1j * (F2 + T)                                   # (N, b, s)
            # ---- hole: sea b at p - q emits photon(-q, mode from gB) and would land in s
            Ep = np.einsum("nhbs,nh->nbs", a_arr, cB.conj())         # emission b -> s
            Ap = np.einsum("nhbs,nh->nbs", e_arr, cB)                # absorption s -> b
            cum2 = np.cumsum(a_arr * cB.conj()[:, :, None, None], 1)
            cum2 = np.concatenate([np.zeros_like(cum2[:, :1]), cum2[:, :-1]], 1)
            F2b = np.einsum("nhbs,nh,nhbs->nbs", e_arr, cB, cum2)
            mu_b = np.exp(-1j * Om_b)
            mu_mp = np.exp(-1j * (Om_s[None, :] + OmB[:, None]))     # (N, s)
            Tb = Ap * Ep * (mu_mp[:, None, :] / (mu_b[:, :, None] - mu_mp[:, None, :]))
            hole = -1j * (F2b + Tb)
            hole_on = hole_flag
            contrib = np.where((sb > 0)[:, :, None], single, hole) if hole_on else single
            accT += np.einsum("n,nbs->s", Wc * norm, contrib)
            # seagull (A^2 of the Peierls phase, dT^2/dk^2 = -T): delta mu/mu = -1/2 e^2 <A_j A_j> per hop; it is
            # state independent and purely imaginary in delta phi (it restores |mu| = 1: check of the normalization)
            accT += 1j * np.sum(Wc * norm * (-0.5) * np.sum(np.abs(cA) ** 2, 1))
        if coulomb:
            VC = 4 * PI * ALPHA / np.einsum("ni,ni->n", Q, Q)
            for oi in tick_end:
                ov = np.einsum("nab,as->nbs", Vc[:, oi].conj(), vs[oi])  # <b,t|s,t>
                accC += np.einsum("n,nb,nbs->s", Wc * VC * 0.5, sb, np.abs(ov) ** 2)
    res = {"E0": Om_s / nt, "dT": np.real(accT) / nt, "dC": accC / nt, "imag": np.imag(accT) / nt, "sgn": sgn_s}
    return res


def dv_exact(th, rule="R", K12=(0.6, 1.05), pd=0.05, axes=(0, 1, 2), **kw):
    """velocity renormalization per axis from the dressed dispersion (spectral, representation free):
    dE_sym = (dE_+ - dE_-)/2 (charge-conjugation symmetric, cancels state-independent constants);
    dv_j / v = [E0(p) dE_sym(p) - m dE_sym(0)] / (E0(p)^2 - m^2),  p = pd * th along axis j.  In units of alpha.
    Also returns dm/m (transverse, Coulomb) at rest."""
    def sym(r, key):
        return 0.5 * (r[key][2:].mean() - r[key][:2].mean())
    r0 = self_energy_exact(np.zeros(3), th, rule, K12, **kw)
    out = {"theta": th, "rule": rule, "K12": list(K12), "pd": pd,
           "dmT_over_m": sym(r0, "dT") / th, "dmC_over_m": sym(r0, "dC") / th, "E0_rest": float(r0["E0"][2])}
    m = r0["E0"][2:].mean()
    for part in ("T", "C"):
        out["dv%s_over_alpha" % part] = [None] * 3
    for j in axes:
        p = np.zeros(3); p[j] = pd * th
        r = self_energy_exact(p, th, rule, K12, **kw)
        E0 = r["E0"][2:].mean()
        for part in ("T", "C"):
            dv = (E0 * sym(r, "d" + part) - m * sym(r0, "d" + part)) / (E0**2 - m**2)
            out["dv%s_over_alpha" % part][j] = float(dv / ALPHA)
    out["dv_over_alpha"] = [None if out["dvT_over_alpha"][j] is None else out["dvT_over_alpha"][j] + out["dvC_over_alpha"][j]
                            for j in range(3)]
    return out

"""vuoto_f4_core.py - phase 4 of 'Vuoto che polarizza' (docs/VUOTO_POLARIZZA.md section 15).

(1) g - 2 with the EXACT vertex of the rule in discrete time (Feynman form, sub-step resolved; no effective Hamiltonian,
    hence no dependence on the time slice by construction).  (2) the dressed light: photon self-energy from the sea.

Conventions (same physics as lab/vuoto_f3_core.py, written in Feynman form so that 2- and 3-vertex loops are treated
alike).  Period = list of operations; a vertex at hop o replaces the operation o (maps slice o-1 to slice o).
Times t = (n, o).  Electron Feynman propagator from (n_a, o_a) to (n_b, o_b) at momentum k, band c:
    G = + |v_{c,o_b-1}><v_{c,o_a}| mu_c^(n_b - n_a)   if c > 0 and t_b > t_a  (strictly)
        - |v_{c,o_b-1}><v_{c,o_a}| mu_c^(n_b - n_a)   if c < 0 and t_b <= t_a (equal time = the sea density:
                                                      two DIFFERENT particles at the same hop, a real two-body term
                                                      of the discrete-time many-body walk)
Photon (light walk, positive frequency): from the earlier vertex to the later one; equal time: either formula.
Insertions in the electron line carry the weight lambda^-(n_out - n_in) and the prefactor mu_{out} ("mu units"):
the dressed quasi-energies solve  det[lambda - F0 - Sigma(lambda)] = 0.
The seagull (A^2) and the A_ext A^2 contact multiply every hop by the same real number: they rescale the whole period
operator and drop out of quasi-energy differences (checked in phase 3 through unitarity); not included.
The contact vertex A_ext A_gamma at the same hop (second-order Peierls term, d^2T/dk^2 = -T) IS included.
Coulomb: instantaneous 4 pi alpha/q^2 after each coin (input, as before), first-order exchange with the sea.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from zz_core import ALPHA
from vuoto_core import ueig, species_weight, PI
from vuoto_f2_core import sphere_rule
import vuoto_f3_core as F3
from vuoto_f3_core import op_list, e_hop, coin, n_ticks, e_bands, _e_chain, _photon_chain

E2 = 4 * PI * ALPHA
SW = {"T": 1.0, "C": 1.0}      # test switches (transverse, Coulomb)


# ------------------------------------------------------------------ geometric sums (Abel-regularized, |x| = 1)
def geo_ge(x, a):
    return x ** a / (1 - x)


def geo_le(x, b):
    return x ** b / (1 - 1 / x)


def dsum(x, y, t1, p1, t3, p3):
    """sum of x^n1 y^n3 over n1 in I1, n3 in I3 split by sign(n3 - n1): returns (gt, eq, lt).
    I = ('ge', alpha) or ('le', beta) with alpha in {0, 1}, beta in {0, -1}."""
    if t1 == "le" and t3 == "ge":
        T = geo_le(x, p1) * geo_ge(y, p3)
        eq = 1.0 if (p1 == 0 and p3 == 0) else 0.0
        return T - eq, eq + 0 * T, 0 * T
    if t1 == "ge" and t3 == "le":
        T = geo_ge(x, p1) * geo_le(y, p3)
        eq = 1.0 if (p1 == 0 and p3 == 0) else 0.0
        return 0 * T, eq + 0 * T, T - eq
    if t1 == "ge" and t3 == "ge":
        xy = x * y
        gt = (y / (1 - y)) * xy ** p1 / (1 - xy)
        lt = (x / (1 - x)) * xy ** p3 / (1 - xy)
        eq = xy ** max(p1, p3) / (1 - xy)
        return gt, eq, lt
    # both 'le': n = -m
    gt_s, eq_s, lt_s = dsum(1 / x, 1 / y, "ge", -p1, "ge", -p3)
    return lt_s, eq_s, gt_s


# ------------------------------------------------------------------ per-chunk ingredients
class Chunk:
    """everything that depends on the photon momenta Q (N,3) for external momenta p_in, p_out"""

    def __init__(self, Q, th, rule, p_in, p_out, ext_axis=None, wq=None):
        self.rule, self.th = rule, th
        self.Q = Q
        self.wq = wq
        self.ops = op_list(rule)
        self.L = len(self.ops)
        self.nt = n_ticks(rule)
        self.hops = [(o + 1, j) for o, (kind, j) in enumerate(self.ops) if kind == "h"]
        self.tick_end = [o + 1 for o, (kind, j) in enumerate(self.ops) if kind == "c"]
        N = len(Q)
        self.N = N
        self.k1 = p_in[None] - Q
        self.k2 = p_out[None] - Q
        self.Om1, V1 = e_bands(self.k1, th, rule)
        self.ch1 = _e_chain(self.k1, th, rule, V1)
        if np.array_equal(p_in, p_out):
            self.Om2, self.ch2 = self.Om1, self.ch1
        else:
            self.Om2, V2 = e_bands(self.k2, th, rule)
            self.ch2 = _e_chain(self.k2, th, rule, V2)
        gp, Omg, om = _photon_chain(Q, rule, self.nt)
        gm, _, _ = _photon_chain(-Q, rule, self.nt)
        self.gam = np.exp(-1j * Omg)                                  # photon eigenvalue per period
        self.norm = E2 / (2 * om)
        js = [j for (_, j) in self.hops]
        oi = [o for (o, _) in self.hops]

        def comp(g):                                                   # (N, nh) component j_o at hop o (after the hop;
            return np.stack([g[:, o, j] for (o, j) in self.hops], 1)   #  equal to before: S_j leaves e_j unchanged)
        # modes of the photon at +q: (gp, gm); at -q: (gm, gp)
        cq = [comp(gp), comp(gm)]
        cmq = [comp(gm), comp(gp)]
        self.Pf = sum(np.einsum("nb,na->nab", c, c.conj()) for c in cq)      # Pf[o1, o3] = sum_h g_o3 conj(g_o1)
        self.Pb = sum(np.einsum("na,nb->nab", c, c.conj()) for c in cmq)     # Pb[o1, o3] = sum_h g_o1 conj(g_o3)
        self.D1 = [e_hop(0.5 * (p_in[None] + self.k1), j)[1] for (_, j) in self.hops]    # photon vertex on the in side
        self.D3 = [e_hop(0.5 * (p_out[None] + self.k2), j)[1] for (_, j) in self.hops]   # on the out side
        self.ext_axis = ext_axis
        if ext_axis is not None:
            kx = 0.5 * (self.k1 + self.k2)
            self.Dx = {o: e_hop(kx, j)[1] for (o, j) in self.hops if j == ext_axis}
            kc1 = 0.5 * (p_in[None] + self.k2)                         # contact at the in vertex (p_in -> k2)
            kc3 = 0.5 * (self.k1 + p_out[None])                        # contact at the out vertex (k1 -> p_out)
            self.C1 = {o: -e_hop(kc1, j)[0] for (o, j) in self.hops if j == ext_axis}
            self.C3 = {o: -e_hop(kc3, j)[0] for (o, j) in self.hops if j == ext_axis}


def ext_states(p, th, rule):
    Om, V = e_bands(p[None], th, rule)
    o = np.argsort(Om[0])
    Om, V = Om[0][o], V[0][:, o]
    return Om, _e_chain(p[None], th, rule, V[None])[0]               # (4,), (L+1, 4, 4)


# ------------------------------------------------------------------ self-energy matrix (2 vertices), mu units
def sigma_T(ch, Om_s, vs, lam, side="in", eqsea=True, kets=None, bras=None):
    """transverse Sigma_{c'c}(lam) at the external momentum (in-side chunk: intermediate at k1), without the mu_{c'}
    prefactor.  kets / bras (optional, for the overlap diagrams): {o: (4, M)} replacing vs[o-1] at the first vertex /
    vs[o] at the second vertex (hops missing from the dict are skipped)."""
    hops = ch.hops
    Vb = ch.ch1 if side == "in" else ch.ch2
    Omb = ch.Om1 if side == "in" else ch.Om2
    Dv = ch.D1 if side == "in" else ch.D3
    mub = np.exp(-1j * Omb)                                            # (N,4)
    sb = Omb > 0
    out = 0
    # amplitudes: A1[o1] = <v_{b,o1}|D|v_{c,o1-1}>  (N,4b,4c);  A3[o3] = <v_{c',o3}|D|v_{b,o3-1}>  (N,4c',4b)
    A1, A3 = {}, {}
    for i, (o, j) in enumerate(hops):
        k_ = vs[o - 1] if kets is None else kets.get(o)
        if k_ is not None:
            A1[i] = np.einsum("nab,nac,cs->nbs", Vb[:, o].conj(), Dv[i], k_)
        b_ = vs[o] if bras is None else bras.get(o)
        if b_ is not None:
            A3[i] = np.einsum("as,nac,ncb->nsb", b_.conj(), Dv[i], Vb[:, o - 1])
    xf = mub * ch.gam[:, None] / lam                                   # forward: (mu_b gamma / lam)^N, N >= ...
    xb = lam * ch.gam[:, None] / mub                                   # backward: (lam gamma / mu_b)^m, m = -N >= ...
    for i1, (o1, j1) in enumerate(hops):
        if i1 not in A1:
            continue
        for i3, (o3, j3) in enumerate(hops):
            if i3 not in A3:
                continue
            # b > 0: N >= 1, or N = 0 with o3 > o1
            Sf = geo_ge(xf, 0 if o3 > o1 else 1)
            # b < 0: N <= -1, or N = 0 with o3 <= o1   (m = -N >= 0 or 1)
            Sb = geo_ge(xb, 0 if o3 <= o1 else 1)
            if o3 == o1 and not eqsea:                                 # test only: drop the equal-time sea term
                Sb = Sb - 1.0
            fac = np.where(sb, Sf * ch.Pf[:, i1, i3][:, None], -Sb * ch.Pb[:, i1, i3][:, None])   # (N, 4b)
            out = out + np.einsum("n,nsb,nb,nbc->sc", ch.wq * ch.norm, A3[i3], fac, A1[i1])
    return out * SW["T"]


def sigma_C(ch, vs, side="in", kets=None, bras=None):
    """Coulomb exchange with the symmetric sea at every tick end, without the mu_{c'} prefactor:
    returns sum_t <v_{c',t}| 1/2 sum_q V_q sum_b sgn_b |b_t><b_t| |v_{c,t}>   (phase units).
    kets / bras: {t: (4, M)} replacing vs[t] on the right / left (overlap diagrams)."""
    Vb = ch.ch1 if side == "in" else ch.ch2
    Omb = ch.Om1 if side == "in" else ch.Om2
    VC = 4 * PI * ALPHA / np.einsum("ni,ni->n", ch.Q, ch.Q)
    out = 0
    for t in ch.tick_end:
        kt = vs[t] if kets is None else kets.get(t)
        bt = vs[t] if bras is None else bras.get(t)
        if kt is None or bt is None:
            continue
        ovk = np.einsum("nab,ac->nbc", Vb[:, t].conj(), kt)             # <b_t|ket>
        ovb = np.einsum("nab,ac->nbc", Vb[:, t].conj(), bt)             # <b_t|bra>
        out = out + np.einsum("n,nb,nbs,nbc->sc", ch.wq * VC * 0.5, np.sign(Omb), ovb.conj(), ovk)
    return out * SW["C"]


def free_prop(p, th, rule, o_from, o_to, v):
    """free evolution at momentum p of the vectors v (4, M) from slice o_from to slice o_to >= o_from"""
    ops = op_list(rule)
    C = coin(th)
    for o in range(o_from + 1, o_to + 1):
        kind, j = ops[o - 1]
        v = (e_hop(p[None], j)[0][0] if kind == "h" else C) @ v
    return v


def overlap_vectors(p_in, p_out, th, rule, ext_axis, vin, vout):
    """kets / bras for the diagrams where the external vertex sits on an external line inside the same period as a
    loop vertex (period bookkeeping of the discrete time): see docs section 15"""
    ops = op_list(rule)
    L = len(ops)
    km = 0.5 * (p_in + p_out)[None]
    xo = [o_ + 1 for o_, (kind, j) in enumerate(ops) if kind == "h" and j == ext_axis]
    Dx = {o: e_hop(km, ext_axis)[1][0] for o in xo}
    kin, bout = {}, {}
    for o in range(1, L + 2):               # slice o-1 ket for a loop vertex at hop o (after an ext vertex o2 < o)
        acc = 0
        for o2 in xo:
            if o2 < o:
                acc = acc + free_prop(p_out, th, rule, o2, o - 1, Dx[o2] @ vin[o2 - 1])
        if not np.isscalar(acc):
            kin[o] = acc
    for o in range(0, L + 1):               # slice-o bra (as a column) before an ext vertex o2 > o
        acc = 0
        for o2 in xo:
            if o2 > o:
                # <vout_{o2}| Dx U(o2-1 <- o): column = U^dag Dx^dag vout_{o2}
                U = free_prop(p_in, th, rule, o, o2 - 1, np.eye(4, dtype=complex))
                acc = acc + U.conj().T @ Dx[o2].conj().T @ vout[o2]
        if not np.isscalar(acc):
            bout[o] = acc
    return kin, bout


# ------------------------------------------------------------------ chunks of the photon grid
def chunks(th, rule, p_in, p_out, ext_axis=None, K12=(0.6, 1.05), nr=64, nct=12, nph=16, kmin_fac=1e-3, chunk=4000,
           kmax=None):
    w = species_weight(*K12)
    Qg, WW = sphere_rule(nr, nct, nph, kmin_fac * th, K12[1] if kmax is None else kmax)
    WW = WW * w(Qg) / (2 * PI) ** 3
    for c0 in range(0, len(Qg), chunk):
        yield Chunk(Qg[c0:c0 + chunk], th, rule, np.asarray(p_in, float), np.asarray(p_out, float), ext_axis,
                    WW[c0:c0 + chunk])


# ------------------------------------------------------------------ vertex (3 vertices + 2 contacts), mu units
def _amps(ch, vin, vout):
    """A1 (N,nh,4b1,4c) at the in side, A3 (N,nh,4c',4b2) at the out side, Ax (N,nx,4b2,4b1) external"""
    V1, V2 = ch.ch1, ch.ch2
    A1 = np.stack([np.einsum("nab,nac,cs->nbs", V1[:, o].conj(), ch.D1[i], vin[o - 1])
                   for i, (o, j) in enumerate(ch.hops)], 1)
    A3 = np.stack([np.einsum("as,nac,ncb->nsb", vout[o].conj(), ch.D3[i], V2[:, o - 1])
                   for i, (o, j) in enumerate(ch.hops)], 1)
    xo = [o for (o, j) in ch.hops if j == ch.ext_axis]
    Ax = np.stack([np.einsum("nab,nac,ncd->nbd", V2[:, o].conj(), ch.Dx[o], V1[:, o - 1]) for o in xo], 1)
    return A1, A3, Ax, xo


def lam_vertex(ch, vin, vout, lam):
    """one-loop transverse vertex Lambda_{c'c}(lam) for the external field along ext_axis (static, inserted once per
    period at n2 = 0), without the mu_{c'} prefactor.  Includes the two contact diagrams."""
    A1, A3, Ax, xo = _amps(ch, vin, vout)
    hidx = {o: i for i, (o, j) in enumerate(ch.hops)}
    mu1, mu2 = np.exp(-1j * ch.Om1), np.exp(-1j * ch.Om2)
    s1, s2 = ch.Om1 > 0, ch.Om2 > 0                                    # (N,4)
    a = (lam / mu1)[:, None, :]                                         # (N,1,4b1)
    b = (mu2 / lam)[:, :, None]                                         # (N,4b2,1)
    g = ch.gam[:, None, None]
    S1 = s1[:, None, :]; S2 = s2[:, :, None]
    sign = np.where(S1, 1.0, -1.0) * np.where(S2, 1.0, -1.0)
    # factors for the 8 keys (o1 < o2, o3 > o2, o3 >= o1): forward part (x = a/g, y = b g), backward (x = a g, y = b/g)
    F = {}
    for k12 in (True, False):
        for k32 in (True, False):
            t1p = ("le", 0 if k12 else -1)      # c1 > 0
            t1m = ("ge", 1 if k12 else 0)       # c1 < 0 (n1 >= 0 if o1 >= o2)
            t3p = ("ge", 0 if k32 else 1)       # c2 > 0
            t3m = ("le", -1 if k32 else 0)      # c2 < 0 (n3 <= 0 if o3 <= o2)
            res = {}
            for direc, (x, y) in (("f", (a / g, b * g)), ("b", (a * g, b / g))):
                out = {}
                parts = [None] * 3
                acc = [0, 0, 0]
                for m1, T1 in ((S1, t1p), (~S1, t1m)):
                    for m2, T3 in ((S2, t3p), (~S2, t3m)):
                        with np.errstate(all="ignore"):
                            gt, eq, lt = dsum(x, y, T1[0], T1[1], T3[0], T3[1])
                        msk = m1 & m2
                        for i_, v in enumerate((gt, eq, lt)):
                            acc[i_] = acc[i_] + np.where(msk, v, 0)
                res[direc] = acc
            for k31 in (True, False):
                F[(k12, k32, k31)] = (sign * (res["f"][0] + (res["f"][1] if k31 else 0)),
                                      sign * (res["b"][2] + (0 if k31 else res["b"][1])))
    W = ch.wq * ch.norm
    out = np.zeros((4, 4), complex)
    for i1, (o1, j1) in enumerate(ch.hops):
        for i3, (o3, j3) in enumerate(ch.hops):
            pf, pb = ch.Pf[:, i1, i3][:, None, None], ch.Pb[:, i1, i3][:, None, None]
            G = 0
            for i2, o2 in enumerate(xo):
                Ff, Fb = F[(o1 < o2, o3 > o2, o3 >= o1)]
                G = G + (pf * Ff + pb * Fb) * Ax[:, i2]
            out += np.einsum("n,nsb,nbd,ndc->sc", W, A3[:, i3], G, A1[:, i1], optimize=True)
    # contact at the in vertex (photon + external at the same hop o2, p_in -> k2), then line 2 to the photon at o3
    V2 = ch.ch2; V1 = ch.ch1
    bb = (mu2 / lam)[:, :, None][:, :, 0]
    gg = ch.gam[:, None]
    for o2 in xo:
        AC1 = np.einsum("nab,nac,cs->nbs", V2[:, o2].conj(), ch.C1[o2], vin[o2 - 1])
        i2 = hidx[o2]
        for i3, (o3, j3) in enumerate(ch.hops):
            fp = ch.Pf[:, i2, i3][:, None] * geo_ge(bb * gg, 0 if o3 > o2 else 1)
            with np.errstate(all="ignore"):
                fm = -ch.Pb[:, i2, i3][:, None] * geo_le(bb / gg, 0 if o3 <= o2 else -1)
            fac = np.where(s2, fp, fm)
            out += np.einsum("n,nsb,nb,nbc->sc", W, A3[:, i3], fac, AC1)
    # contact at the out vertex (k1 -> p_out at o2), line 1 from the photon at o1
    aa = lam / mu1
    for o2 in xo:
        AC3 = np.einsum("as,nac,ncb->nsb", vout[o2].conj(), ch.C3[o2], V1[:, o2 - 1])
        i2 = hidx[o2]
        for i1, (o1, j1) in enumerate(ch.hops):
            with np.errstate(all="ignore"):
                fp = ch.Pf[:, i1, i2][:, None] * geo_le(aa / gg, 0 if o1 < o2 else -1)
            fm = -ch.Pb[:, i1, i2][:, None] * geo_ge(aa * gg, 0 if o1 >= o2 else 1)
            fac = np.where(s1, fp, fm)
            out += np.einsum("n,nsb,nb,nbc->sc", W, AC3, fac, A1[:, i1])
    return out * SW["T"]


def lam_coulomb(ch, vin, vout):
    """Coulomb vertex: -V_q delta P_minus at every tick end, delta P_minus = -delta G(t,t) from the external insertion;
    returns the phase-unit matrix (the effective coupling is -i mu_{c'} times it)"""
    V1, V2 = ch.ch1, ch.ch2
    mu1, mu2 = np.exp(-1j * ch.Om1), np.exp(-1j * ch.Om2)
    s1, s2 = (ch.Om1 > 0)[:, None, :], (ch.Om2 > 0)[:, :, None]
    X = mu1[:, None, :] / mu2[:, :, None]                               # (N,4b2,4b1)
    VC = 4 * PI * ALPHA / np.einsum("ni,ni->n", ch.Q, ch.Q)
    xo = [o for (o, j) in ch.hops if j == ch.ext_axis]
    Ax = {o: np.einsum("nab,nac,ncd->nbd", V2[:, o].conj(), ch.Dx[o], V1[:, o - 1]) for o in xo}
    out = np.zeros((4, 4), complex)
    with np.errstate(all="ignore"):
        for ts in ch.tick_end:
            dG = 0
            for o2 in xo:
                f_pm = -geo_ge(X, 0 if o2 > ts else 1)                  # b1 > 0, b2 < 0: vertex after t
                f_mp = -geo_le(X, 0 if o2 <= ts else -1)                # b1 < 0, b2 > 0: vertex before t
                fac = np.where(s1 & ~s2, f_pm, 0) + np.where(~s1 & s2, f_mp, 0)
                dG = dG + fac * Ax[o2]
            L = np.einsum("as,nab->nsb", vout[ts].conj(), V2[:, ts])         # <c'|b2>
            R = np.einsum("nab,ac->nbc", V1[:, ts].conj(), vin[ts])          # <b1|c>
            out += np.einsum("n,nsb,nbd,ndc->sc", ch.wq * VC, L, dG, R)
    return out * SW["C"]


def op_list_ticks(rule):
    return [o_ + 1 for o_, (kind, j) in enumerate(op_list(rule)) if kind == "c"]


def tree(p_in, p_out, th, rule, ext_axis, vin, vout, Om_out):
    ops = op_list(rule)
    km = 0.5 * (np.asarray(p_in, float) + np.asarray(p_out, float))[None]
    out = np.zeros((4, 4), complex)
    for o_, (kind, j) in enumerate(ops):
        if kind == "h" and j == ext_axis:
            o = o_ + 1
            out += vout[o].conj().T @ e_hop(km, j)[1][0] @ vin[o - 1]
    return np.exp(-1j * Om_out)[:, None] * out


# ------------------------------------------------------------------ assembly
def vertex_all(th, rule, p_in, p_out, ext_axis, lam_shift=0.0, dlam=1e-6, grid=None, parts=False):
    """effective first-order coupling of the static external field between the dressed positive states at p_in and
    p_out (2x2 on P), one loop, plus the bare tree; mu units"""
    grid = grid or {}
    p_in, p_out = np.asarray(p_in, float), np.asarray(p_out, float)
    Omi, vin = ext_states(p_in, th, rule)
    Omo, vout = ext_states(p_out, th, rule)
    P = [2, 3]; Nn = [0, 1]
    lam = np.exp(-1j * 0.5 * (Omi[P].mean() + Omo[P].mean())) * np.exp(-1j * lam_shift)
    F1 = tree(p_in, p_out, th, rule, ext_axis, vin, vout, Omo)
    S_in = np.zeros((4, 4), complex); S_out = np.zeros((4, 4), complex)
    dS_in = np.zeros((4, 4), complex); dS_out = np.zeros((4, 4), complex)
    C_in = np.zeros((4, 4), complex); C_out = np.zeros((4, 4), complex)
    Lam = np.zeros((4, 4), complex); LC = np.zeros((4, 4), complex)
    lp, lm = lam * np.exp(-1j * dlam), lam * np.exp(1j * dlam)
    kin, bout = overlap_vectors(p_in, p_out, th, rule, ext_axis, vin, vout)
    ticks = op_list_ticks(rule)
    # Coulomb at slice t: ext vertex before (o2 <= t: ket at slice t) or after (o2 > t: bra at slice t)
    kin_c = {t: kin[t + 1] for t in ticks if (t + 1) in kin}
    bout_c = {t: bout[t] for t in ticks if t in bout}
    for ch in chunks(th, rule, p_in, p_out, ext_axis, **grid):
        S_in += sigma_T(ch, Omi, vin, lam, "in")
        S_out += sigma_T(ch, Omo, vout, lam, "out")
        dS_in += (sigma_T(ch, Omi, vin, lp, "in") - sigma_T(ch, Omi, vin, lm, "in")) / (lp - lm)
        dS_out += (sigma_T(ch, Omo, vout, lp, "out") - sigma_T(ch, Omo, vout, lm, "out")) / (lp - lm)
        C_in += sigma_C(ch, vin, "in")
        C_out += sigma_C(ch, vout, "out")
        Lam += lam_vertex(ch, vin, vout, lam)
        LC += lam_coulomb(ch, vin, vout)
        # overlap diagrams (external vertex on an external line in the same period as a loop vertex)
        Lam += sigma_T(ch, Omo, vout, lam, "out", kets=kin) + sigma_T(ch, Omi, vin, lam, "in", bras=bout)
        LC += sigma_C(ch, vout, "out", kets={t: kin_c[t] for t in kin_c}, bras={t: vout[t] for t in kin_c}) \
            + sigma_C(ch, vin, "in", kets={t: vin[t] for t in bout_c}, bras=bout_c)
    mui, muo = np.exp(-1j * Omi), np.exp(-1j * Omo)
    Sig_in = mui[:, None] * (S_in - 1j * C_in)
    Sig_out = muo[:, None] * (S_out - 1j * C_out)
    z_in = mui[:, None] * dS_in
    z_out = muo[:, None] * dS_out
    LamT = muo[:, None] * Lam
    LamC = -1j * muo[:, None] * LC
    ix = np.ix_
    Vz = 0.5 * (z_out[ix(P, P)] @ F1[ix(P, P)] + F1[ix(P, P)] @ z_in[ix(P, P)])
    Rin = 1.0 / (lam - mui[Nn]); Rout = 1.0 / (lam - muo[Nn])
    leg = F1[ix(P, Nn)] @ (Rin[:, None] * Sig_in[ix(Nn, P)]) + Sig_out[ix(P, Nn)] @ (Rout[:, None] * F1[ix(Nn, P)])
    V = F1[ix(P, P)] + Vz + LamT[ix(P, P)] + LamC[ix(P, P)] + leg
    out = {"F1": F1[ix(P, P)], "V": V, "lam": lam, "Omi": Omi, "Omo": Omo}
    if parts:
        out.update({"Vz": Vz, "LamT": LamT[ix(P, P)], "LamC": LamC[ix(P, P)], "leg": leg,
                    "Sig_in": Sig_in, "Sig_out": Sig_out})
    return out


def dmu_onshell(th, rule, p, grid=None, eqsea=True):
    """on-shell first-order shift delta mu of the 4 states at p (diagonal Sigma_ss(mu_s)), Feynman form; also the
    phase-3-comparable delta Omega per period = i delta mu / mu (transverse, Coulomb)"""
    grid = grid or {}
    p = np.asarray(p, float)
    Om, v = ext_states(p, th, rule)
    mu = np.exp(-1j * Om)
    dT = np.zeros(4, complex); dC = np.zeros(4, complex)
    for ch in chunks(th, rule, p, p, None, **grid):
        for s in range(4):
            dT[s] += sigma_T(ch, Om, v, mu[s], "in", eqsea)[s, s]
        dC += np.diag(sigma_C(ch, v, "in"))
    return {"Om": Om, "dmu": mu * (dT - 1j * dC), "dOmT": 1j * dT, "dOmC": dC}


# ------------------------------------------------------------------ g - 2 and the dressed velocity (same engine)
def dmv(th, rule, grid=None, pd=0.05, axes=(0, 1, 2)):
    """dressed mass and velocity per axis (phase-3 definitions, Feynman-form engine): dm/m and dv_j (relative)"""
    def sym(r):
        d = np.real(r["dOmT"] + r["dOmC"])
        return 0.5 * (d[2:].mean() - d[:2].mean())
    r0 = dmu_onshell(th, rule, np.zeros(3), grid)
    nt = n_ticks(rule)
    m = r0["Om"][2:].mean() / nt
    dm = sym(r0) / nt
    out = {"dm_over_m": float(dm / m), "dv": [None] * 3,
           "rel_P": complex(np.mean(r0["dmu"][2:] / np.exp(-1j * r0["Om"][2:])))}
    for j in axes:
        p = np.zeros(3); p[j] = pd * th
        r = dmu_onshell(th, rule, p, grid)
        E0 = r["Om"][2:].mean() / nt
        out["dv"][j] = float((E0 * sym(r) / nt - m * dm) / (E0**2 - m**2))
    return out


def g2_pair(th, rule, qax, jax, grid=None, Qext=None):
    """Breit frame p = -Q/2 e_q, p' = +Q/2 e_q, static external A along jax: returns the dressed/bare ratio of the
    magnetic coupling in phase units r = tr(T^dag V)/tr(T^dag T) - delta mu_P / mu_P (seagull-free combination)"""
    grid = grid or {}
    kmin = grid.get("kmin_fac", 1e-3) * th
    Q = 1e-2 * kmin if Qext is None else Qext
    e = np.zeros(3); e[qax] = 1.0
    r = vertex_all(th, rule, -0.5 * Q * e, 0.5 * Q * e, jax, grid=grid, parts=True)
    T, V = r["F1"], r["V"]
    nT = np.trace(T.conj().T @ T)
    c = np.trace(T.conj().T @ V) / nT
    resid = np.linalg.norm(V - c * T) / np.linalg.norm(V)
    # on-shell delta mu / mu of the P pair at p (same grid), from the Sigma matrices already computed
    Si = r["Sig_in"]
    rel = np.mean(np.diag(Si)[2:] / np.exp(-1j * r["Omi"][2:]))
    out = {"c": complex(c), "rel": complex(rel), "r": complex(c - rel), "resid": float(resid), "Q": Q,
           "tree_norm": float(np.sqrt(nT.real))}
    for k in ("Vz", "LamT", "LamC", "leg"):
        out[k] = complex(np.trace(T.conj().T @ r[k]) / nT)
    return out


# ------------------------------------------------------------------ the dressed light: vacuum polarization of the sea
def kgrid(mode, th, K12=(0.6, 1.05), nr=64, nct=12, nph=16, M=48, kmin_fac=1e-2):
    """electron loop momenta: 'species' = sphere rule around k = 0 up to K2 with weight w (one-body weight returned
    separately), 'bz' = uniform midpoint grid on the whole cube [-pi, pi]^3 (all doublers, no weight)"""
    if mode == "bz":
        x = -PI + (np.arange(M) + 0.5) * 2 * PI / M
        K = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1).reshape(-1, 3)
        return K, np.full(len(K), 1.0 / M**3)
    Kg, WW = sphere_rule(nr, nct, nph, kmin_fac * th, K12[1])
    return Kg, WW / (2 * PI) ** 3


def photon_pol(th, rule, q, mode="species", K12=(0.6, 1.05), chunk=4000, eqsea=True, seagull=True, wmode="occ", **kg):
    """photon self-energy in discrete time, Feynman form (fermion loop, sign -1, Feynman electron propagators; the
    equal-time sea-sea term and the A^2 seagull included).  Returns the 2x2 matrix Sigma_gamma (delta gamma / gamma)
    in the basis of the two photon modes at q (helicity +: light walk at q, helicity -: light walk at -q), on shell,
    and the bare photon period phase.  delta Omega = i Sigma_gamma (2x2), per period."""
    q = np.asarray(q, float)
    nt = n_ticks(rule)
    gp, Omg, om = _photon_chain(q[None], rule, nt)
    gm, Omg2, _ = _photon_chain(-q[None], rule, nt)
    ops = op_list(rule)
    hops = [(o_ + 1, j) for o_, (kind, j) in enumerate(ops) if kind == "h"]
    G = np.array([[gp[0, o, j] for (o, j) in hops], [gm[0, o, j] for (o, j) in hops]])     # (2 modes, nh)
    lam = np.exp(-1j * Omg[0])
    norm = E2 / (2 * om[0])
    Kall, Wall = kgrid(mode, th, K12, **kg)
    w = species_weight(*K12)
    out = np.zeros((2, 2), complex)
    nsea = 0.0
    for c0 in range(0, len(Kall), chunk):
        K = Kall[c0:c0 + chunk]
        Wt = Wall[c0:c0 + chunk]
        Kq = K + q[None]
        if mode == "species":
            wk, wkq = w(K), w(Kq)
            nsea += np.sum(Wt * wk) * 2
        else:
            wk, wkq = np.ones(len(K)), np.ones(len(K))
            nsea += np.sum(Wt) * 2
        # sea occupation n(k) = w(k) (declared): pair weights by the filled state(s)
        if wmode == "occ":
            w_pm, w_mp, w_mm = wk, wkq, wk * wkq
        elif wmode == "occ_mean":
            w_pm, w_mp, w_mm = wk, wkq, 0.5 * (wk + wkq)
        else:                                               # phase-1 symmetric pair weight
            w_pm = w_mp = w_mm = np.sqrt(wk * wkq)
        Omc, Vc0 = e_bands(K, th, rule)
        Omb, Vb0 = e_bands(Kq, th, rule)
        Vc = _e_chain(K, th, rule, Vc0)
        Vb = _e_chain(Kq, th, rule, Vb0)
        sc, sb = Omc > 0, Omb > 0
        x = (np.exp(-1j * Omb)[:, :, None] / np.exp(-1j * Omc)[:, None, :]) / lam        # (N, b, c)
        Sb, Sc = sb[:, :, None], sc[:, None, :]
        km = K + 0.5 * q[None]
        Aup, Adn = [], []
        for (o, j) in hops:
            D = e_hop(km, j)[1]
            Aup.append(np.einsum("nab,nac,ncd->nbd", Vb[:, o].conj(), D, Vc[:, o - 1]))   # c(k) -> b(k+q)
            Adn.append(np.einsum("nab,nac,ncd->nbd", Vc[:, o].conj(), D, Vb[:, o - 1]))   # b(k+q) -> c(k)
        with np.errstate(all="ignore"):
            for i1, (o1, j1) in enumerate(hops):
                for i3, (o3, j3) in enumerate(hops):
                    f_pm = geo_ge(x, 0 if o3 > o1 else 1)
                    f_mp = geo_le(x, 0 if o3 < o1 else -1)
                    fac = (np.where(Sb & ~Sc, f_pm, 0) * w_pm[:, None, None]
                           + np.where(~Sb & Sc, f_mp, 0) * w_mp[:, None, None])
                    if o3 == o1 and eqsea:
                        fac = fac + np.where(~Sb & ~Sc, -1.0, 0) * w_mm[:, None, None]
                    loop = np.einsum("n,ncb,nbc,nbc->", Wt, Adn[i3], Aup[i1], fac)
                    out += loop * np.outer(G[:, i3].conj(), G[:, i1])
    out *= norm
    if seagull:
        out += -norm * nsea * (G.conj() @ G.T)
    return {"Sigma": out, "Om": float(Omg[0]), "Om_minus": float(Omg2[0]), "om_tick": float(om[0]), "nt": nt}


def photon_dc(th, rule, axis, mode="species", qfs=(0.05, 0.1), **kw):
    """dressed light speed along a direction (axis index or unit vector): delta c_gamma / c (in units of alpha) for
    each of the two dressed polarizations, after removing the q-independent photon mass of the species cut:
    omega^2 = mu^2 + c^2 q^2  ->  delta c / c = [W2 dW2 - W1 dW1] / (W2^2 - W1^2)  (W = period phase).
    Polarization eigenvalues are matched between the two q by overlap of the 2x2 eigenvectors."""
    e = np.zeros(3)
    if np.ndim(axis) == 0:
        e[int(axis)] = 1.0
    else:
        e = np.asarray(axis, float) / np.linalg.norm(axis)
    rows = []
    for qf in qfs:
        r = photon_pol(th, rule, qf * th * e, mode, **kw)
        dOm = 1j * r["Sigma"]
        H = 0.5 * (dOm + dOm.conj().T)                 # Hermitian part (the anti-Hermitian part = unitarity residue)
        ev, U = np.linalg.eigh(H)
        rows.append({"qf": qf, "Om": r["Om"], "Om_minus": r["Om_minus"], "dOm": dOm, "ev": ev, "U": U,
                     "imag_rel": float(np.abs(np.linalg.eigvals(0.5 * (dOm - dOm.conj().T) / 1j)).max())})
    a, b = rows[0], rows[1]
    ov = np.abs(a["U"].conj().T @ b["U"])
    perm = [int(np.argmax(ov[i])) for i in range(2)]
    W1, W2 = a["Om"], b["Om"]
    dc = [float((W2 * b["ev"][perm[i]] - W1 * a["ev"][i]) / (W2**2 - W1**2) / ALPHA) for i in range(2)]
    mass2 = float((W1 * W2) ** 2 * 0 + (W2**2 * 2 * W1 * a["ev"].mean() - W1**2 * 2 * W2 * b["ev"].mean())
                  / (W2**2 - W1**2) / ALPHA)
    return {"axis": e.tolist(), "dc_over_alpha": dc, "dc_mean": float(np.mean(dc)), "biref": float(dc[1] - dc[0]),
            "pol_vectors": [[complex(x) for x in a["U"][:, i]] for i in range(2)],
            "mass2_over_alpha": mass2, "imag_rel_over_alpha": [r["imag_rel"] / ALPHA / r["Om"] for r in rows],
            "Om": [W1, W2], "ev": [r["ev"].tolist() for r in rows], "offdiag": [complex(r["dOm"][0, 1]) for r in rows]}

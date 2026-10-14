"""tavg_core.py - il "blocco del nucleo" (studio tavg_): nucleo + elettroni interni agglomerati in una pallina.

L'agglomeratore fa sempre la stessa cosa: blocca, impara la regola del blocco dal livello sotto, controlla. Qui:
  * livello sotto  = l'atomo intero sullo zoom radiale (tavg_radial: HF ristretto del termine di Hund, regola non
                     relativistica con m_kin);
  * blocco         = nucleo + guscio interno (1s^2 per B..Ne; [Ne] per Na..Ar; ...);
  * regola imparata= come il blocco agisce sugli elettroni esterni: un potenziale per ogni l (s, p, d), ricavato
                     INVERTENDO l'equazione di HF dell'elettrone esterno, con l'orbitale esterno reso liscio dentro r_c
                     (forma di Troullier-Martins, stessa carica dentro r_c: "conservazione della norma") e uguale
                     all'orbitale vero fuori; poi forma separabile (Kleinman-Bylander) per la griglia 3-D;
  * controllo      = lo stesso blocco messo in configurazioni diverse (catione, eccitazioni) contro l'atomo intero.
Nessuna tabella: tutto viene dal livello sotto. In linguaggio standard e' uno "pseudopotenziale a conservazione di
norma generato da HF tutti-elettroni": lo dichiariamo come tale (debito: forma TM e r_c).

Unita' atomiche.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import math
import numpy as np
from scipy.optimize import brentq
import tavg_radial as tr

L_LETTER = tr.L_LETTER
R_TAIL_ADD = 2.75          # a0: beyond max(rc) + this, the block is a point charge Zv (dev: C eps to 2e-5 Eh)


# ------------------------------------------------------------------------------------------------ TM pseudization
def _derivs_at(g, P, l, rc, deg=8):
    """p(r) = ln(P/r^(l+1)) and its first four r-derivatives at rc, from a local polynomial fit on the log grid"""
    r = g.r
    sel = (r > rc * 0.82) & (r < rc * 1.22)
    x = r[sel] - rc
    y = np.log(np.abs(P[sel]) / r[sel] ** (l + 1))
    c = np.polyfit(x, y, deg)
    return [float(np.polyval(np.polyder(c, n), 0.0)) if n else float(np.polyval(c, 0.0)) for n in range(5)]


def tm_pseudize(g, P, l, rc):
    """Troullier-Martins nodeless orbital r^(l+1) exp(p(r)) inside rc, = P outside; returns phi_t, coefficients"""
    r = g.r
    pd = _derivs_at(g, P, l, rc)
    target = g.integ(np.where(r < rc, P * P, 0.0))
    pw = np.array([0, 2, 4, 6, 8, 10, 12])

    def coeffs(c2):
        c4 = -c2 * c2 / (2 * l + 5)
        # unknowns c0, c6, c8, c10, c12 ; equations: n-th derivative of p at rc, n = 0..4
        A = np.zeros((5, 5)); b = np.zeros(5)
        idx = [0, 6, 8, 10, 12]
        for n in range(5):
            for j, p in enumerate(idx):
                A[n, j] = _dpow(p, n, rc)
            b[n] = pd[n] - c2 * _dpow(2, n, rc) - c4 * _dpow(4, n, rc)
        sol = np.linalg.solve(A, b)
        return np.array([sol[0], c2, c4, sol[1], sol[2], sol[3], sol[4]])

    def pfun(c, rr):
        return sum(ci * rr ** p for ci, p in zip(c, pw))

    rin = r[r < rc]

    def norm_res(c2):
        c = coeffs(c2)
        ph = rin ** (l + 1) * np.exp(pfun(c, rin))
        return float(np.sum(ph * ph * rin) * g.h) - target

    grid = np.linspace(-30, 30, 1201)
    vals = np.array([norm_res(x) for x in grid])
    roots = []
    for i in range(len(grid) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            roots.append(brentq(norm_res, grid[i], grid[i + 1], xtol=1e-14))
    if not roots:
        raise RuntimeError(f"TM: no norm-conserving solution for l={l}, rc={rc}")
    c2 = min(roots, key=abs)
    c = coeffs(c2)
    sgn = float(np.sign(P[np.searchsorted(r, rc)]))
    phi = np.where(r < rc, r ** (l + 1) * np.exp(pfun(c, np.minimum(r, rc))), sgn * P)   # nodeless, positive
    return phi, c


def _dpow(p, n, r):
    """n-th derivative of r^p at r"""
    if p < n:
        return 0.0
    f = 1.0
    for i in range(n):
        f *= (p - i)
    return f * r ** (p - n)


def kinetic_tail(g, phi, l, kcut):
    """kinetic energy (Eh) of the orbital phi (P-like, on the log grid) carried by |k| > kcut"""
    from scipy.special import spherical_jn
    r = g.r
    ru = np.linspace(1e-4, 20.0, 6001); dr = ru[1] - ru[0]
    k = np.linspace(0, 40, 1201)
    f = np.interp(ru, r, phi)
    pk = math.sqrt(2 / math.pi) * (spherical_jn(l, np.outer(k, ru)) * (f * ru)[None, :]).sum(1) * dr
    dens = pk ** 2 * k ** 2
    return float(np.trapezoid(np.where(k > kcut, 0.5 * k ** 2 * dens, 0.0), k)), float(np.trapezoid(dens, k))


def choose_rc(g, P, l, kcut, budget, step=0.025, rmax=3.5):
    """agglomerator rule for the block size: the SMALLEST cut radius whose smooth orbital the coarse level can carry,
    i.e. kinetic energy beyond kcut <= budget. Scan from just outside the outermost node."""
    r = g.r
    sgn = np.sign(P)
    nodes = r[1:][(sgn[1:] * sgn[:-1] < 0) & (r[1:] < 10)]
    r0 = (nodes.max() * 1.2 if len(nodes) else 0.2)
    rc = math.ceil(r0 / step) * step
    hist = []
    while rc <= rmax:
        try:
            phi, _ = tm_pseudize(g, P, l, rc)
            t, nrm = kinetic_tail(g, phi, l, kcut)
            hist.append((rc, t))
            if t <= budget:
                return rc, hist
        except RuntimeError:
            hist.append((rc, None))
        rc = round(rc + step, 6)
    raise RuntimeError(f"no rc <= {rmax} meets the tail budget for l={l}")


# ------------------------------------------------------------------------------------------------ block learning
def relabel(config_val):
    """valence shells relabelled to the lowest n for each l (the block removes the inner ones)"""
    return [(l + 1 + sum(1 for (n2, l2, q2) in config_val[:i] if l2 == l), l, q) for i, (n, l, q) in enumerate(config_val)]


class Block:
    """learned core block. Fields: Z, Zv, ls, rc{l}, V{l} (semilocal, on the AE grid), Vloc, lloc, kb{l: (chi, E)},
    phi{l}, eps{l}, E_core (constant to add to valence energies), ae (the AE reference)"""

    def __init__(self, Z, ref_config, ncore_shells, rc, mkin=1.0, lloc=None, extra_refs=None, r_tail=None, h=0.03,
                 rmax=45.0, ae=None):
        """ref_config: AE config string (Hund ground term or a reference ionic/excited state); ncore_shells: number of
        leading shells in the block; rc: dict l -> cut radius; extra_refs: dict l -> AE config string used to get that
        channel when it is not occupied in ref_config (e.g. p for Li from 1s2 2p1)."""
        self.Z, self.mkin = Z, mkin
        cfg = tr.config_from_string(ref_config) if isinstance(ref_config, str) else ref_config
        self.ae = ae if ae is not None else tr.Atom(Z, cfg, mkin=mkin, h=h, rmax=rmax)
        if ae is None:
            self.ae.scf()
        g = self.g = self.ae.g
        r = g.r
        self.core = cfg[:ncore_shells]
        self.val = cfg[ncore_shells:]
        self.Zv = Z - sum(q for (n, l, q) in self.core)
        self.rc = dict(rc)
        self.val_ps = relabel(self.val)
        # pseudo orbitals for each valence shell (one shell per l allowed in this study)
        ls_val = [l for (n, l, q) in self.val]
        assert len(set(ls_val)) == len(ls_val), "one valence shell per l"
        self.phi, self.eps, self.V = {}, {}, {}
        # reference pseudo-atom (valence only) with the pseudo orbitals, to build the valence two-electron part
        refs = {l: (self.ae, ncore_shells + i) for i, l in enumerate(ls_val)}
        for l in rc:
            if l not in ls_val:
                cfg2 = tr.config_from_string(extra_refs[l])
                ae2 = tr.Atom(Z, cfg2, mkin=mkin, grid=g); ae2.scf()
                a2 = [i for i, c in enumerate(cfg2) if c[1] == l][-1]
                refs[l] = (ae2, a2)
        self.V = {}
        for l, (aeX, a) in refs.items():
            P = aeX.P[a]
            phi, c = tm_pseudize(g, P, l, self.rc[l])
            u = P / np.sqrt(r)
            eps = float(u @ aeX.fock(a) @ u) * g.h
            # valence part of that reference: shells after the core
            valX = aeX.config[ncore_shells:]
            psX = tr.Atom(self.Zv, relabel(valX), mkin=mkin, grid=g, nuc=np.zeros_like(r))
            psX.P = {}
            for i, (n2, l2, q2) in enumerate(valX):
                PP = aeX.P[ncore_shells + i]
                if l2 == l and ncore_shells + i == a:
                    psX.P[i] = phi
                else:
                    psX.P[i] = tm_pseudize(g, PP, l2, self.rc.get(l2, self.rc[l]))[0]
            ia = a - ncore_shells
            Fw = psX.fock(ia)                       # kinetic + centrifugal + valence two-electron (nuc = 0)
            ut = phi / np.sqrt(r)
            Vl = eps - (Fw @ ut) / (r ** 1.5 * phi)
            self.phi[l], self.eps[l] = phi, eps
            self.V[l] = Vl
        # tail: beyond r_tail the block is a point charge Zv (core charge and core exchange decay exponentially)
        rt = r_tail if r_tail is not None else max(self.rc.values()) + R_TAIL_ADD
        self.r_tail = rt
        w = 0.5 * (1 + np.tanh((r - rt) / 0.15))
        for l in self.V:
            self.V[l] = (1 - w) * self.V[l] + w * (-self.Zv / r)
        self.lloc = lloc if lloc is not None else max(self.V)
        self.Vloc = self.V[self.lloc]
        self.kb = {}
        for l in self.V:
            if l == self.lloc:
                continue
            dV = self.V[l] - self.Vloc
            chi = dV * self.phi[l]
            E = g.integ(self.phi[l] * chi)
            self.kb[l] = (chi, E)
        # core energy constant: E_AE(ref) - E_pseudo(ref)
        ps = self.pseudo_atom(self.val_ps)
        self.ps_ref = ps.scf()
        self.E_core = self.ae.E - ps.E

    def pseudo_atom(self, config_ps, semilocal=False):
        r = self.g.r
        if semilocal:
            return tr.Atom(self.Zv, config_ps, mkin=self.mkin, grid=self.g, nuc=self.Vloc,
                           extV={l: self.V[l] - self.Vloc for l in self.V if l != self.lloc})
        return tr.Atom(self.Zv, config_ps, mkin=self.mkin, grid=self.g, nuc=self.Vloc,
                       extNL={l: [self.kb[l]] for l in self.kb})

    def momentum_tail(self, kcuts):
        """kinetic energy (Eh) carried by |k| > kcut in each pseudo orbital (radial Bessel transform on a uniform r grid)"""
        from scipy.special import spherical_jn
        r = self.g.r
        ru = np.linspace(1e-4, 20.0, 8001); dr = ru[1] - ru[0]
        k = np.linspace(0, 30, 1501)
        out = {}
        for l, phi in self.phi.items():
            f = np.interp(ru, r, phi)
            J = spherical_jn(l, np.outer(k, ru))
            pk = math.sqrt(2 / math.pi) * (J * (f * ru)[None, :]).sum(1) * dr
            dens = pk ** 2 * k ** 2
            T = np.trapezoid(0.5 * k ** 2 * dens, k)
            out[l] = dict(norm=float(np.trapezoid(dens, k)), T=float(T),
                          tail={float(kc): float(np.trapezoid(np.where(k > kc, 0.5 * k ** 2 * dens, 0.0), k)) for kc in kcuts})
        return out


if __name__ == "__main__":
    import time
    t0 = time.time()
    b = Block(6, "1s2 2s2 2p2", 1, {0: 1.2, 1: 1.3})
    print("E_AE", b.ae.E, "E_core", b.E_core, "eps", b.eps, f"{time.time() - t0:.1f}s")
    print("pseudo ref", b.ps_ref)
    print("tail", b.momentum_tail([6.0]))
    for cs in ["2s2 2p1", "2s1 2p3", "2s2 2p3", "2s2"]:
        pass

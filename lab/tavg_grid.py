"""tavg_grid.py - Hartree-Fock sulla griglia agglomerata 3-D con i blocchi-nucleo imparati (studio tavg_).

UNA funzione pubblica per gli altri studi (molecole: LiH, N2, O2, H2O, NH3, CH4, CO2, ...):

    from tavg_grid import hf_griglia
    res = hf_griglia(nuclei=[(8, (0, 0, 0)), (1, (0, 1.43, 1.11)), (1, (0, -1.43, 1.11))],   # (Z, posizione in a0)
                     carica=0, spin=0,           # spin = N_su - N_giu (2S); il numero di elettroni viene da Z e carica
                     h=0.25, L=16.0)             # cella (a0) e lato della scatola (a0)
    res["E"]          energia totale (Eh): blocchi (costante imparata) + elettroni esterni + repulsione fra blocchi
    res["orbitali"]   lista di array N^3 (reali, normalizzati con sum(phi^2) h^3 = 1), res["occ"] = [(n_su, n_giu)]
    res["eps"]        energie degli orbitali (diagonale di Fock, Eh)
    res["densita"]    densita' degli elettroni esterni (elettroni / a0^3) sulla griglia
    res["griglia"]    coordinate 1-D (a0) della griglia (stessa per x, y, z; l'origine e' il punto N//2)

Cosa fa: minimizza l'energia di UN determinante (orbitali reali, occupazioni intere) con il gradiente proiettato
sulla varieta' degli orbitali ortonormali (precondizionato col cinetico), in tre varianti: 'R' (ristretto: stessi
orbitali per su e giu', anche a guscio aperto, come l'ROHF), 'U' (non ristretto: orbitali diversi per spin).
I nuclei sono BLOCCHI IMPARATI (tavg_core.Block): carica esterna Zv, potenziale locale + proiettori separabili;
idrogeno compreso (blocco senza elettroni interni). Le posizioni sono libere (non agganciate alla griglia): i blocchi
sono costruiti nello spazio k con la fase e^{-i k.R}.
Cinetico: 'regola' = banda della regola k^2/(2 m_kin) spettrale (default); 'imparato' = stencil imparato a0/4 o a0/8 di
aggmol (solo h = a0/4 o a0/8; vale fino a k ~ 6/a0, vedi documento).
Coulomb isolato: densita' (e coppie per lo scambio) sulla griglia N^3, imbottita a (2N)^3, nucleo 4 pi (1 - cos qR)/q^2
con R = L (esatto se la nube sta in una sfera di diametro L).
Unita' atomiche."""
import os, sys, time, json, math
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
import scipy.fft as sfft
from scipy.special import spherical_jn, erf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "out")
BLOCKDIR = os.path.join(OUT, "tavg_blocchi")


# ------------------------------------------------------------------------------------------------ blocks (from JSON)
def load_block(Z):
    """learned block for nuclear charge Z, as written by tavg_blocks.py (JSON: radial grid r, Vloc, kb list)"""
    fn = os.path.join(BLOCKDIR, f"Z{Z:02d}.json")
    if not os.path.exists(fn):
        raise FileNotFoundError(f"blocco per Z={Z} non imparato ({fn}): eseguire tavg_blocks.py")
    d = json.load(open(fn))
    d["r"] = np.array(d["r"]); d["Vloc"] = np.array(d["Vloc"])
    d["kb"] = [dict(l=k["l"], chi=np.array(k["chi"]), E=k["E"]) for k in d["kb"]]
    return d


def ylm_real(l, m, x, y, z, r):
    """real solid-harmonic angular parts (normalised on the sphere); r>0 assumed where used"""
    with np.errstate(invalid="ignore", divide="ignore"):
        if l == 0:
            return np.full(np.broadcast(x, y, z).shape, 0.5 / math.sqrt(math.pi))
        ux, uy, uz = x / r, y / r, z / r
        if l == 1:
            c = math.sqrt(3 / (4 * math.pi))
            return c * {-1: uy, 0: uz, 1: ux}[m]
        if l == 2:
            c = 0.5 * math.sqrt(15 / math.pi)
            return {-2: c * ux * uy, -1: c * uy * uz, 0: 0.25 * math.sqrt(5 / math.pi) * (3 * uz * uz - 1),
                    1: c * ux * uz, 2: 0.5 * c * (ux * ux - uy * uy)}[m]
    raise ValueError(l)


def radial_ft(r, f, l, k):
    """F(k) = 4 pi int j_l(kr) f(r) r^2 dr, f given on the log grid r (trapezoid in x = ln r)"""
    x = np.log(r); dx = x[1] - x[0]
    w = f * r ** 3 * dx
    out = np.empty_like(k)
    for i0 in range(0, len(k), 256):
        kk = k[i0:i0 + 256]
        out[i0:i0 + 256] = 4 * math.pi * (spherical_jn(l, np.outer(kk, r)) * w[None, :]).sum(1)
    return out


# ------------------------------------------------------------------------------------------------ the grid machine
class Machine:
    def __init__(self, nuclei, h=0.25, L=16.0, kinetic="regola", mkin=1.0, rg=0.6):
        self.h = h
        self.N = N = int(round(L / h))
        self.L = N * h
        self.dV = h ** 3
        self.x = (np.arange(N) - N // 2) * h
        self.mkin = mkin
        k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
        k1r = 2 * np.pi * np.fft.rfftfreq(N, d=h)
        KX, KY, KZ = np.meshgrid(k1, k1, k1r, indexing="ij", sparse=True)
        self.k2 = KX ** 2 + KY ** 2 + KZ ** 2
        if kinetic == "regola":
            self.Tk = self.k2 / (2 * mkin)
        elif kinetic == "imparato":
            self.Tk = self._learned_T()
        else:
            raise ValueError(kinetic)
        self.pre = 1.0 / (self.Tk + 1.0)
        # Coulomb kernel on the (2N)^3 padded grid, spherical cutoff R = L
        M = 2 * N
        q1 = 2 * np.pi * np.fft.fftfreq(M, d=h); q1r = 2 * np.pi * np.fft.rfftfreq(M, d=h)
        QX, QY, QZ = np.meshgrid(q1, q1, q1r, indexing="ij", sparse=True)
        q2 = QX ** 2 + QY ** 2 + QZ ** 2
        Rc = self.L
        with np.errstate(divide="ignore", invalid="ignore"):
            self.Kq = np.where(q2 > 0, 4 * np.pi * (1 - np.cos(np.sqrt(q2) * Rc)) / q2, 2 * np.pi * Rc ** 2)
        del q2
        self.M = M
        # nuclei / blocks
        self.nuclei = [(int(Z), np.array(p, float)) for Z, p in nuclei]
        self.blocks = {Z: load_block(Z) for Z in set(Z for Z, _ in self.nuclei)}
        self.Zv = [self.blocks[Z]["Zv"] for Z, _ in self.nuclei]
        self._build_ext(KX, KY, KZ, rg)

    def _learned_T(self):
        import aggmol_core as mc
        B = {0.25: 1024, 0.125: 512}[round(self.h, 6)]
        rec = mc.level_rec(B)
        T = mc.kinetic_grid(rec, self.N, B * mc.TH) / mc.ALPHA ** 2
        return np.ascontiguousarray(T[:, :, : self.N // 2 + 1])

    def _build_ext(self, KX, KY, KZ, rg):
        """local potential (erf long-range part in real space + short-range rest from k space) and KB projectors"""
        N, h = self.N, self.h
        X, Y, Zg = np.meshgrid(self.x, self.x, self.x, indexing="ij", sparse=True)
        k = np.sqrt(self.k2)
        kmax = float(k.max())
        ktab = np.linspace(0, kmax * 1.0001, 1200)
        Vloc = np.zeros((N, N, N))
        self.proj = []                                     # list of (beta array N^3, 1/E)
        vol = self.L ** 3
        # grid origin: index 0 is at x = -N//2 h; phase for a function centred at R relative to the FFT origin
        x0 = self.x[0]
        for (Z, R) in self.nuclei:
            b = self.blocks[Z]
            r = b["r"]; Zv = b["Zv"]
            rr = np.sqrt((X - R[0]) ** 2 + (Y - R[1]) ** 2 + (Zg - R[2]) ** 2)
            with np.errstate(divide="ignore", invalid="ignore"):
                lr = np.where(rr > 1e-12, -Zv * erf(rr / rg) / rr, -Zv * 2 / (math.sqrt(math.pi) * rg))
            Vloc += lr
            # short-range rest dV(r) = Vloc(r) + Zv erf(r/rg)/r -> k space
            dv = b["Vloc"] + Zv * erf(r / rg) / r
            ft = radial_ft(r, dv, 0, ktab)
            Fk = np.interp(k, ktab, ft)
            ph = np.exp(-1j * (KX * (R[0] - x0) + KY * (R[1] - x0) + KZ * (R[2] - x0)))
            Vloc += self._kspace_to_grid(Fk * ph)
            for kb in b["kb"]:
                l = kb["l"]
                ftl = radial_ft(r, kb["chi"] / r, l, ktab)          # chi is P-like (r * radial function)
                Fl = np.interp(k, ktab, ftl)
                for m in range(-l, l + 1):
                    with np.errstate(invalid="ignore", divide="ignore"):
                        kk = np.where(k > 0, k, 1.0)
                        Y = ylm_real(l, m, KX + 0 * k, KY + 0 * k, KZ + 0 * k, kk)
                        if l > 0:
                            Y = np.where(k > 0, Y, 0.0)
                    beta = self._kspace_to_grid(((-1j) ** l) * Fl * Y * ph, real_part=True)
                    self.proj.append((beta, 1.0 / kb["E"], Z))
        self.Vloc = Vloc
        # block-block repulsion (point charges Zv)
        Enn = 0.0
        for i in range(len(self.nuclei)):
            for j in range(i):
                Enn += self.Zv[i] * self.Zv[j] / np.linalg.norm(self.nuclei[i][1] - self.nuclei[j][1])
        self.Enn = Enn
        self.Ecore = sum(self.blocks[Z]["E_core"] for Z, _ in self.nuclei)

    def _kspace_to_grid(self, F, real_part=True):
        """f(r) = (1/V) sum_k F(k) e^{ik.r} for F given on the rfft half-space of the N^3 grid (Hermitian assumed)"""
        N = self.N
        return sfft.irfftn(F, s=(N, N, N)) * N ** 3 / self.L ** 3

    # -------------------------------------------------------------- Coulomb
    def coulomb(self, n):
        """isolated potential of a charge density n (per a0^3) on the grid"""
        N, M = self.N, self.M
        big = np.zeros((M, M, M))
        big[:N, :N, :N] = n
        V = sfft.irfftn(sfft.rfftn(big) * self.Kq, s=(M, M, M))[:N, :N, :N]
        return V

    def T(self, phi):
        return sfft.irfftn(self.Tk * sfft.rfftn(phi), s=phi.shape)

    def precond(self, g):
        return sfft.irfftn(self.pre * sfft.rfftn(g), s=g.shape)

    def Vnl(self, phi):
        out = np.zeros_like(phi)
        for beta, invE, _ in self.proj:
            c = float(np.sum(beta * phi)) * self.dV
            out += c * invE * beta
        return out


class State:
    """orbitals: list of arrays; occ: list of (nu, nd); sets: list of index lists that must stay mutually orthonormal"""

    def __init__(self, orbs, occ, sets):
        self.orbs, self.occ, self.sets = orbs, occ, sets


def energy_grad(mc, st, want_grad=True):
    dV = mc.dV
    orbs, occ = st.orbs, st.occ
    n = len(orbs)
    rho = np.zeros_like(orbs[0]); rs = [np.zeros_like(orbs[0]), np.zeros_like(orbs[0])]
    for p, (nu, nd) in zip(orbs, occ):
        rho += (nu + nd) * p * p
        rs[0] += nu * p * p; rs[1] += nd * p * p
    VH = mc.coulomb(rho)
    E_one = 0.0; E_kin = 0.0; E_loc = 0.0; E_nl = 0.0
    Hp = []
    for p, (nu, nd) in zip(orbs, occ):
        tp = mc.T(p); nlp = mc.Vnl(p)
        ek = float(np.sum(p * tp)) * dV; en = float(np.sum(p * nlp)) * dV
        el = float(np.sum(p * p * mc.Vloc)) * dV
        E_kin += (nu + nd) * ek; E_nl += (nu + nd) * en; E_loc += (nu + nd) * el
        Hp.append(tp + nlp + mc.Vloc * p)
    E_H = 0.5 * float(np.sum(rho * VH)) * dV
    # exchange: same-spin pairs
    E_x = 0.0
    Kp = [[np.zeros_like(orbs[0]) for _ in range(n)] for _ in range(2)] if want_grad else None
    for i in range(n):
        for j in range(i, n):
            wu = occ[i][0] * occ[j][0]; wd = occ[i][1] * occ[j][1]
            if wu == 0 and wd == 0:
                continue
            pij = orbs[i] * orbs[j]
            v = mc.coulomb(pij)
            e = float(np.sum(pij * v)) * dV
            f = 1.0 if i == j else 2.0
            E_x -= 0.5 * f * (wu + wd) * e
            if want_grad:
                for s, w in ((0, wu), (1, wd)):
                    if w:
                        Kp[s][i] += v * orbs[j]
                        if j != i:
                            Kp[s][j] += v * orbs[i]
    E = E_kin + E_loc + E_nl + E_H + E_x
    parts = dict(kin=E_kin, loc=E_loc, nl=E_nl, H=E_H, x=E_x)
    if not want_grad:
        return E, parts, None, None
    G = []; F = []
    for i, (p, (nu, nd)) in enumerate(zip(orbs, occ)):
        base = Hp[i] + VH * p
        Fu = base - Kp[0][i]; Fd = base - Kp[1][i]
        G.append(2 * (nu * Fu + nd * Fd))
        F.append((nu * Fu + nd * Fd) / max(nu + nd, 1))
    return E, parts, G, F


def _orthonormalise(mc, st):
    for s in st.sets:
        A = np.array([st.orbs[i].ravel() for i in s])
        S = A @ A.T * mc.dV
        w, U = np.linalg.eigh(S)
        X = U @ np.diag(w ** -0.5) @ U.T
        B = X @ A
        for k, i in enumerate(s):
            st.orbs[i] = B[k].reshape(st.orbs[i].shape)


def _tangent(mc, st, D):
    """remove from directions D the part that leaves the orthonormality manifold (keeps rotations between orbitals of
    different occupation): D_i -= sum_j phi_j sym(<phi_j|D_i>)"""
    out = [d.copy() for d in D]
    for s in st.sets:
        A = np.array([st.orbs[i].ravel() for i in s])
        Dm = np.array([D[i].ravel() for i in s])
        C = Dm @ A.T * mc.dV                       # C[i, j] = <phi_j | D_i>
        Cs = 0.5 * (C + C.T)
        R = Cs @ A
        for k, i in enumerate(s):
            out[i] = D[i] - R[k].reshape(D[i].shape)
    return out


def minimise(mc, st, tol=1e-7, maxiter=400, log=print, budget=None, t0=None):
    t0 = t0 or time.time()
    E, parts, G, F = energy_grad(mc, st)
    d_old = None; gp_old = None; g_old = None
    step = 0.5
    hist = [E]
    small = 0
    conv = False
    for it in range(maxiter):
        g = _tangent(mc, st, G)
        gp = _tangent(mc, st, [mc.precond(x) for x in g])
        gg = sum(float(np.sum(a * b)) for a, b in zip(g, gp)) * mc.dV
        if d_old is not None:
            num = sum(float(np.sum(a * (b - c))) for a, b, c in zip(g, gp, gp_old)) * mc.dV
            beta = max(0.0, num / g_old)
            d = [-a + beta * b for a, b in zip(gp, d_old)]
        else:
            d = [-a for a in gp]
        d = _tangent(mc, st, d)
        slope = sum(float(np.sum(a * b)) for a, b in zip(G, d)) * mc.dV
        if slope >= 0:
            d = [-a for a in gp]; slope = -gg
        orbs0 = [p.copy() for p in st.orbs]

        def at(t):
            st.orbs = [p + t * x for p, x in zip(orbs0, d)]
            _orthonormalise(mc, st)
            return energy_grad(mc, st, want_grad=False)[0]

        E1 = at(step)
        curv = (E1 - E - slope * step) / step ** 2
        tstar = -slope / (2 * curv) if curv > 0 else 2 * step
        tstar = float(np.clip(tstar, 0.1 * step, 4 * step))
        st.orbs = [p + tstar * x for p, x in zip(orbs0, d)]
        _orthonormalise(mc, st)
        En, parts, Gn, F = energy_grad(mc, st)
        if En > E and E1 < E:
            st.orbs = [p + step * x for p, x in zip(orbs0, d)]; _orthonormalise(mc, st)
            En, parts, Gn, F = energy_grad(mc, st); tstar = step
        elif En > E:
            st.orbs = orbs0; step *= 0.3; d_old = None
            log(f"   it {it}: risalita, passo ridotto a {step:.3g}")
            continue
        step = 0.5 * step + 0.5 * tstar
        dE = En - E
        E, G = En, Gn
        d_old, gp_old, g_old = d, gp, gg
        hist.append(E)
        if it % 10 == 0:
            log(f"   it {it}: E_val = {E:.9f}  dE = {dE:.2e}  |g| = {math.sqrt(abs(gg)):.2e}  passo {step:.3g}  ({time.time() - t0:.0f} s)")
        small = small + 1 if abs(dE) < tol else 0
        if small >= 3 and math.sqrt(abs(gg)) < 1e-4:
            conv = True
            break
        if budget and time.time() - t0 > budget:
            break
    return E, parts, F, hist, conv


def occupations(n_up, n_dn, restricted=True):
    """spatial orbitals with occupations: restricted -> doubly occupied first, then singly (spin up)"""
    if restricted:
        nd = min(n_up, n_dn); ns = abs(n_up - n_dn)
        occ = [(1, 1)] * nd + ([(1, 0)] if n_up > n_dn else [(0, 1)]) * ns
        return occ, [list(range(len(occ)))]
    occ = [(1, 0)] * n_up + [(0, 1)] * n_dn
    return occ, [list(range(n_up)), list(range(n_up, n_up + n_dn))]


def atomic_guess(mc, norb, restricted, n_up, n_dn, order=None):
    """starting orbitals from the blocks' own pseudo orbitals (s, then p x, y, z, then d), superposed over the nuclei;
    no shape is kept: the minimiser relaxes everything"""
    X, Y, Z = np.meshgrid(mc.x, mc.x, mc.x, indexing="ij", sparse=True)
    funcs = []
    for (Zn, R) in mc.nuclei:
        b = mc.blocks[Zn]
        rr = np.sqrt((X - R[0]) ** 2 + (Y - R[1]) ** 2 + (Z - R[2]) ** 2) + 1e-12
        for l, ms in ((0, [0]), (1, [1, -1, 0]), (2, [2, -2, 1, -1, 0])):
            if str(l) not in b["phi"]:
                continue
            ph = np.interp(rr, b["r"], np.array(b["phi"][str(l)])) / rr
            for m in ms:
                funcs.append(ph * ylm_real(l, m, X - R[0], Y - R[1], Z - R[2], rr))
    if order is not None:
        funcs = [funcs[i] for i in order]
    rng = np.random.default_rng(7)
    while len(funcs) < norb:
        funcs.append(rng.standard_normal(funcs[0].shape) * np.exp(-0.05 * (X ** 2 + Y ** 2 + Z ** 2)))
    return funcs


def hf_griglia(nuclei, carica=0, spin=None, h=0.25, L=16.0, restricted=True, kinetic="regola", mkin=1.0, tol=1e-7,
               maxiter=400, guess=None, order=None, occ=None, budget=None, log=print, machine=None):
    """Hartree-Fock on the agglomerated grid for arbitrary nuclei (Z, position in a0). See the module docstring."""
    mc = machine or Machine(nuclei, h=h, L=L, kinetic=kinetic, mkin=mkin)
    nel = int(round(sum(mc.Zv))) - carica
    if spin is None:
        spin = nel % 2
    n_up = (nel + spin) // 2; n_dn = nel - n_up
    if occ is None:
        occ, sets = occupations(n_up, n_dn, restricted)
    else:
        sets = [list(range(len(occ)))] if restricted else None
    norb = len(occ)
    if guess is None:
        if restricted:
            g = atomic_guess(mc, norb, restricted, n_up, n_dn, order)[:norb]
        else:
            g0 = atomic_guess(mc, max(n_up, n_dn), restricted, n_up, n_dn, order)
            g = g0[:n_up] + g0[:n_dn]
    else:
        g = [np.array(x, float) for x in guess]
    st = State([x.copy() for x in g], occ, sets)
    _orthonormalise(mc, st)
    t0 = time.time()
    Ev, parts, F, hist, conv = minimise(mc, st, tol=tol, maxiter=maxiter, log=log, budget=budget, t0=t0)
    eps = [float(np.sum(p * f)) * mc.dV for p, f in zip(st.orbs, F)]
    rho = sum((nu + nd) * p * p for p, (nu, nd) in zip(st.orbs, occ))
    E = Ev + mc.Ecore + mc.Enn
    return dict(E=E, E_val=Ev, E_blocchi=mc.Ecore, E_nn=mc.Enn, parti=parts, orbitali=st.orbs, occ=occ, eps=eps,
                densita=rho, griglia=mc.x, convergito=conv, storia=hist, secondi=time.time() - t0, h=mc.h, L=mc.L,
                N=mc.N, cinetico=kinetic, machine=mc)

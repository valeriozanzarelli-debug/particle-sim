"""vuoto_levels.py - production (pre-registered): r-space Uehling potential (V3) and level shifts (V4) from the
3-D lattice polarization table lab/out/vuoto_3d.json. Writes lab/out/vuoto_levels.json.
Run: timeout 590 python lab/vuoto_levels.py [theta_vac]   (default 0.0125)"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.special import gamma as Gfun
from zz_core import RadialWalk, dirac_formula, ALPHA, U_CHI
from vuoto_core import pitR_uehling, uehling_V, MHz, PI, C2_UEHLING

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
T0 = time.time()
TH_VAC = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0125


# ------------------------------------------------------------------ the lattice polarization (xyz average)
def load_pitR(th, K=(0.6, 1.05), n=64):
    D = json.load(open(os.path.join(OUT, "vuoto_3d.json")))
    beta = 6.0 if th < 0.02 else 5.0

    def g(dn, qm):
        k = "proj|%g|%s|%.8g|%d|%g|%g|%g" % (th, dn, qm, n, beta, K[0], K[1])
        return D[k][0]
    qs = sorted({float(k.split("|")[3]) for k in D if k.startswith("proj|%g|x|" % th) and k.endswith("|%d|%g|%g|%g" % (n, beta, K[0], K[1]))})
    qs = [q for q in qs if q not in (0.05, 0.025)]
    tab = []
    p0 = {}
    for dn in ("x", "y", "z"):
        p0[dn] = (4 * g(dn, 0.025) - g(dn, 0.05)) / 3
    for q in qs:
        tab.append(np.mean([p0[dn] - g(dn, q) for dn in ("x", "y", "z")]))
    return np.array(qs), np.array(tab), p0


def make_interp(qs, tab):
    """pitR_lat(q): spline of log(pitR) vs log q inside the table; q^2 law (lattice's own first point) below"""
    lq, lp = np.log(qs), np.log(tab)
    sp = CubicSpline(lq, lp)
    c2 = tab[0] / qs[0] ** 2

    def f(q):
        q = np.asarray(q, float)
        out = np.where(q < qs[0], c2 * q * q, 0.0)
        inside = (q >= qs[0]) & (q <= qs[-1])
        out = np.where(inside, np.exp(sp(np.log(np.clip(q, qs[0], qs[-1])))), out)
        return out          # 0 beyond the table (truncation, declared)
    return f


# ------------------------------------------------------------------ radial walk with an arbitrary potential
class RadialWalkV(RadialWalk):
    """zz_core.RadialWalk with V(r) = -Z alpha/r + dV(r) (same rule, extra potential phase)"""

    def __init__(self, kappa, a, rmax, Zalpha, dV):
        super().__init__(kappa, a, rmax, Zalpha=Zalpha)
        r = self.r
        V = -Zalpha / r + dV(r)
        hx, hz = kappa / r, np.ones_like(r)
        nrm = np.sqrt(hx**2 + hz**2)
        for attr, dt in (("M", a), ("Mh", 0.5 * a)):
            phi = dt * nrm
            M = np.zeros((self.N, 2, 2), complex)
            M[:, 0, 0] = np.cos(phi) - 1j * np.sin(phi) * hz / nrm
            M[:, 1, 1] = np.cos(phi) + 1j * np.sin(phi) * hz / nrm
            M[:, 0, 1] = -1j * np.sin(phi) * hx / nrm
            M[:, 1, 0] = -1j * np.sin(phi) * hx / nrm
            M *= np.exp(-1j * V * dt)[:, None, None]
            setattr(self, attr, np.einsum("ab,nbc,cd->nad", U_CHI.conj().T, M, U_CHI))


def state(kappa, n, th, Z, rmax, dV=None):
    Za = Z * ALPHA
    rw = RadialWalk(kappa, th, rmax, Zalpha=Za) if dV is None else RadialWalkV(kappa, th, rmax, Za, dV)
    E0 = dirac_formula(n, kappa, Za)
    om, u, res = rw.eigen(E0, iters=10)
    g, f = rw.gf(u)
    return om, rw.r, (g**2 + f**2) * th           # weights summing to 1


def formfactor(r, wgt, qs):
    out = np.empty(len(qs))
    for i, q in enumerate(qs):
        x = q * r
        out[i] = np.sum(wgt * np.where(x > 1e-8, np.sin(x) / np.maximum(x, 1e-300), 1.0))
    return out


def shift_q(pitR, r, wgt, Z, qmax):
    """first-order shift <V_U> = -(2 Z alpha/pi) int_0^qmax Pi_R(q) rho~(q) dq, Pi_R = 4 pi alpha pitR"""
    qs = np.concatenate([np.linspace(0, 0.02, 41)[1:], np.exp(np.linspace(np.log(0.02), np.log(qmax), 1500))[1:]])
    rho = formfactor(r, wgt, qs)
    integ = 4 * PI * ALPHA * pitR(qs) * rho
    return -(2 * Z * ALPHA / PI) * np.trapezoid(integ, qs)


# ------------------------------------------------------------------ analytic references (EXAM ONLY)
def ref_shift(P, rmax):
    """<V_U> = int P(r) V_U(r) dr, P = radial probability density (normalized), point nucleus Z in P's closure"""
    return quad(lambda r: P(r) * uehling_V(r, Za=ZREF[0] * ALPHA), 0, rmax, limit=400, points=[1e-3, 0.1, 1, 10], epsrel=1e-9)[0]


ZREF = [1]


def P_1s(lam):
    return lambda r: 4 * lam**3 * r * r * np.exp(-2 * lam * r)


def P_2s(lam):
    return lambda r: 0.5 * lam**3 * (1 - lam * r / 2) ** 2 * r * r * np.exp(-lam * r)


def P_dirac1s(Za):
    gm = np.sqrt(1 - Za * Za)
    lam = Za
    nrm = Gfun(2 * gm + 1) / (2 * lam) ** (2 * gm + 1)
    return lambda r: r ** (2 * gm) * np.exp(-2 * lam * r) / nrm


def main():
    R = {"theta_vac": TH_VAC}
    qs, tab, p0 = load_pitR(TH_VAC)
    pit = make_interp(qs, tab)
    qmax = qs[-1]
    R["qmax_over_m"] = qmax
    R["pit0_xyz"] = p0
    # ---- V3: r-space potential, nucleus smeared by exp(-(q/qc)^2) on both sides
    qc = 0.15 / TH_VAC          # AMENDED (lab/out/vuoto_amend.json): was 0.3/theta
    qg = np.unique(np.concatenate([np.exp(np.linspace(np.log(1e-3), np.log(0.5), 400)), np.arange(0.5, qmax, 0.005)]))
    damp = np.exp(-(qg / qc) ** 2)
    lat = pit(qg) * damp
    ref = np.array([pitR_uehling(q) for q in qg]) * damp
    rr = np.array([0.05, 0.1, 0.2, 0.5, 1, 2, 3])
    Vl = -(2 / (PI * rr)) * np.trapezoid(np.sin(np.outer(rr, qg)) * (4 * PI * ALPHA * lat / qg)[None, :], qg, axis=1)
    Vr = -(2 / (PI * rr)) * np.trapezoid(np.sin(np.outer(rr, qg)) * (4 * PI * ALPHA * ref / qg)[None, :], qg, axis=1)
    Vpt = np.array([uehling_V(r, Za=1.0) for r in rr])
    R["V3_r"] = rr.tolist(); R["V3_ratio_lat_ref"] = (Vl / Vr).tolist(); R["V3_smeared_ref_over_point"] = (Vr / Vpt).tolist()
    err = np.abs(Vl / Vr - 1)
    R["V3_PASS"] = bool(np.all(err[:-1] <= 0.02) and err[-1] <= 0.03)
    print("V3", np.round(Vl / Vr - 1, 5), "smeared/point", np.round(Vr / Vpt, 4), flush=True)
    # ---- V4: H levels
    out = {}
    ZREF[0] = 1
    lam = ALPHA
    refs = {"1s": ref_shift(P_1s(lam), 60 / lam), "2s": ref_shift(P_2s(lam), 120 / lam)}
    lead = {n: -4 * ALPHA * ALPHA**4 / (15 * PI * n**3) for n in (1, 2)}
    out["H_ref_MHz"] = {k: MHz(v) for k, v in refs.items()}
    out["H_leading_formula_MHz"] = {str(n): MHz(v) for n, v in lead.items()}
    for name, kappa, n, rmax in (("1s", -1, 1, 40 / ALPHA), ("2s", -1, 2, 80 / ALPHA), ("2p1/2", 1, 2, 80 / ALPHA)):
        vals = {}
        for th in (0.05, 0.025):
            om, r, wgt = state(kappa, n, th, 1, rmax)
            vals[th] = shift_q(pit, r, wgt, 1, qmax)
            print(name, th, MHz(vals[th]), "MHz", "%.0f s" % (time.time() - T0), flush=True)
        ext = (4 * vals[0.025] - vals[0.05]) / 3
        out[name] = {"theta_atom": {str(k): MHz(v) for k, v in vals.items()}, "extrap_MHz": MHz(ext)}
        if name in refs:
            out[name]["ratio_to_ref"] = ext / refs[name]
    R["H"] = out
    R["V4a_PASS"] = bool(abs(out["2s"]["ratio_to_ref"] - 1) <= 0.01)
    R["V4b_PASS"] = bool(abs(out["1s"]["ratio_to_ref"] - 1) <= 0.01)
    R["V4c_PASS"] = bool(abs(out["2p1/2"]["extrap_MHz"]) < 0.01)
    # ---- Nd59+ 1s
    Z = 60
    ZREF[0] = Z
    Za = Z * ALPHA
    refNd = ref_shift(P_dirac1s(Za), 40.0)
    nd = {}
    for th in (0.05, 0.025, 0.0125):
        om, r, wgt = state(-1, 1, th, Z, 40.0)
        nd[th] = shift_q(pit, r, wgt, Z, qmax)
        print("Nd 1s", th, nd[th] * 510998.95, "eV", flush=True)
    gm2 = 2 * np.sqrt(1 - Za * Za)
    s = 2 ** gm2
    ext2 = (s * nd[0.0125] - nd[0.025]) / (s - 1)
    # three-point: E0 + c th^(2 gamma) + d th^2
    A = np.array([[1, t ** gm2, t ** 2] for t in (0.05, 0.025, 0.0125)])
    ext3 = np.linalg.solve(A, np.array([nd[t] for t in (0.05, 0.025, 0.0125)]))[0]
    R["Nd"] = {"theta_atom_eV": {str(k): v * 510998.95 for k, v in nd.items()}, "extrap2_eV": ext2 * 510998.95,
               "extrap3_eV": ext3 * 510998.95, "ref_eV": refNd * 510998.95, "ratio2": ext2 / refNd, "ratio3": ext3 / refNd,
               "leading_formula_eV": -4 * ALPHA * Za**4 / (15 * PI) * 510998.95}
    R["V4d_PASS"] = bool(abs(ext2 / refNd - 1) <= 0.02)
    print("Nd ref", refNd * 510998.95, "eV; ratio", ext2 / refNd, ext3 / refNd, flush=True)
    # ---- V4e nonperturbative: put the smeared lattice V_U into the walk
    qcN = 0.15 / TH_VAC         # AMENDED, as qc
    qg2 = np.unique(np.concatenate([np.exp(np.linspace(np.log(1e-3), np.log(0.5), 300)), np.arange(0.5, qmax, 0.01)]))
    PiR = 4 * PI * ALPHA * pit(qg2) * np.exp(-(qg2 / qcN) ** 2)

    def dV(r):
        out = np.empty_like(r)
        for i0 in range(0, len(r), 2000):
            rs = r[i0:i0 + 2000]
            out[i0:i0 + 2000] = -(2 * Za / (PI * rs)) * np.trapezoid(np.sin(np.outer(rs, qg2)) * (PiR / qg2)[None, :], qg2, axis=1)
        return out
    th = 0.025
    om0, r, w0 = state(-1, 1, th, Z, 40.0)
    om1, _, _ = state(-1, 1, th, Z, 40.0, dV=dV)
    # first-order with the same smeared potential
    exp1 = np.sum(w0 * dV(r))
    R["V4e"] = {"dE_eigen_eV": (om1 - om0) * 510998.95, "dE_first_order_eV": exp1 * 510998.95,
                "ratio": (om1 - om0) / exp1}
    R["V4e_PASS"] = bool(abs((om1 - om0) / exp1 - 1) <= 0.01)
    print("V4e", R["V4e"], flush=True)
    R["time_s"] = time.time() - T0
    json.dump(R, open(os.path.join(OUT, "vuoto_levels%s.json" % ("" if TH_VAC == 0.0125 else "_th%g" % TH_VAC)), "w"), indent=1)
    print(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()

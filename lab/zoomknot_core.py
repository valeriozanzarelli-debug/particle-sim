"""zoomknot_core.py - multi-resolution ("zoom") Skyrme knot: nested cubic grids, finest in the core.

Levels are listed finest first. Level l: spacing h_l, points i*h_l for i = -N_l..N_l (origin is a grid point),
owned half-width R_l (|x|_inf <= R_l), R_l a multiple of h_{l+1}; the coarsest level owns everything (Dirichlet
vacuum outside, as knot_core). Composite field:
  * free (owned) points of level l: |x|_inf <= R_l - h_{l+1}; those inside R_{l-1} - h_l are INJECTED from level l-1;
  * points of level l with |x|_inf > R_l - h_{l+1} (one coarse cell of band + ghost layer) are PROLONGED from level l+1 by tensor-product cubic
    Lagrange interpolation, then renormalised to the 3-sphere.
Energy (conservative partition of volume): E = sum_l h_l^3 sum_x W_l(x) [k2_l e2 + k4_l e4](x),
  W_l = T_l - T_{l-1},  T_l(x) = prod_axes trapezoid indicator of the box R_l (1 inside, 1/2 on the face, 0 outside),
  so sum_l W_l integrates every region exactly once with a 2nd-order rule. e2, e4 = knot_core.static_density split.
k2_l, k4_l = couplings (1 = bare lattice; learned values = "agglomerate" rule). All levels relax jointly (L-BFGS on
u_l with psi_l = u_l * sqrt(h_ref/h_l)... implemented as psi = u / sqrt(h): equal diagonal Hessian scale per level)."""
import math, numpy as np, torch
import knot_core as kc

torch.set_num_threads(1)
DT = torch.float64


def ax(N, h):
    return torch.arange(-N, N + 1, dtype=DT) * h


def trap1(x, R, tol=1e-9):
    a = x.abs()
    one = torch.ones_like(x)
    return torch.where(a < R - tol, one, torch.where(a < R + tol, 0.5 * one, 0 * one))


def box_trap(N, h, R):
    if R is None:
        return torch.ones((2 * N + 1,) * 3, dtype=DT)
    t = trap1(ax(N, h), R)
    return t[:, None, None] * t[None, :, None] * t[None, None, :]


def box_mask(N, h, R):
    t = (ax(N, h).abs() <= R + 1e-9).to(DT)
    return t[:, None, None] * t[None, :, None] * t[None, None, :]


def lagrange_matrix(Nf, hf, Nc, hc):
    """P (nf x nc): cubic Lagrange interpolation from the coarse axis to the fine axis."""
    xf = ax(Nf, hf).numpy(); P = np.zeros((2 * Nf + 1, 2 * Nc + 1))
    for i, x in enumerate(xf):
        s = x / hc; j0 = int(math.floor(s + 1e-9)); t = s - j0
        if abs(t) < 1e-9:
            P[i, j0 + Nc] = 1.0; continue
        nodes = [j0 - 1, j0, j0 + 1, j0 + 2]
        for a, ja in enumerate(nodes):
            w = 1.0
            for b, jb in enumerate(nodes):
                if b != a:
                    w *= (s - jb) / (ja - jb)
            if 0 <= ja + Nc <= 2 * Nc:
                P[i, ja + Nc] = w
    return torch.tensor(P, dtype=DT)


def prolong(phi_c, P):
    x = torch.einsum("ai,dijk->dajk", P, phi_c)
    x = torch.einsum("bj,dajk->dabk", P, x)
    return torch.einsum("ck,dabk->dabc", P, x)


def densities(phi, h):
    """e2, e4 exactly as knot_core.static_density (forward/backward differences, 2nd order)."""
    gp, gm = kc.derivs_fb(phi, h)
    n2p = [(x * x).sum(0) for x in gp]; n2m = [(x * x).sum(0) for x in gm]
    e2 = 0.5 * sum(n2p[i] + n2m[i] for i in range(3))
    e4 = torch.zeros_like(e2)
    for i in range(3):
        for j in range(i + 1, 3):
            for gi, ni in ((gp[i], n2p[i]), (gm[i], n2m[i])):
                for gj, nj in ((gp[j], n2p[j]), (gm[j], n2m[j])):
                    e4 = e4 + 0.25 * (ni * nj - (gi * gj).sum(0) ** 2)
    return e2, e4, (gp, gm)


class Hier:
    def __init__(self, levels, k2=None, k4=None):
        """levels: list of dict(h, N, R) finest first; coarsest has R = None."""
        self.L = levels
        n = len(levels)
        self.k2 = list(k2) if k2 is not None else [1.0] * n
        self.k4 = list(k4) if k4 is not None else [1.0] * n
        self.W, self.M, self.P, self.inj, self.Mi = [], [], [], [], []
        for l, lv in enumerate(levels):
            h, N, R = lv["h"], lv["N"], lv["R"]
            T = box_trap(N, h, R)
            Tf = torch.zeros_like(T) if l == 0 else box_trap(N, h, levels[l - 1]["R"])
            self.W.append(T - Tf)
            # conforming band (fix, see ZOOM_NODO.md): fine points in R-w < |x|_inf <= R (w = coarser h) are NOT free,
            # they are the cubic interpolation of the coarse field, so the fine field cannot kink at the interface.
            wl = levels[l + 1]["h"] if l + 1 < n else 0.0
            self.M.append(box_mask(N, h, R - wl) if R is not None else None)
            self.Mi.append(None if l == 0 else box_mask(N, h, levels[l - 1]["R"] - h))
            if l + 1 < n:
                c = levels[l + 1]
                r = c["h"] / h; assert abs(r - round(r)) < 1e-9; r = int(round(r))
                assert abs(R / c["h"] - round(R / c["h"])) < 1e-9, "R must be a multiple of the coarser h"
                assert abs(R / h - round(R / h)) < 1e-9 and N * h >= R + 2 * h - 1e-9, "need >= 2 ghost points"
                need = N * h + 2 * c["h"]
                Rc = c["R"] if c["R"] is not None else c["N"] * c["h"]
                assert need <= Rc + 1e-9, f"prolongation buffer: level {l} needs coarse owned to {need}, has {Rc}"
                self.P.append(lagrange_matrix(N, h, c["N"], c["h"]))
                self.inj.append((r, int(round(R / c["h"])) - 1))
            else:
                self.P.append(None); self.inj.append(None)
        self.sites = sum((2 * lv["N"] + 1) ** 3 for lv in levels)

    def grids(self, l):
        lv = self.L[l]
        return kc.grid((2 * lv["N"] + 1,) * 3, lv["h"])

    def init_u(self, prof):
        return [kc.hedgehog(*self.grids(l), prof) * math.sqrt(lv["h"]) for l, lv in enumerate(self.L)]

    def fields(self, us):
        n = len(self.L)
        own = []
        for l in range(n):
            phi = kc.normalize(us[l])           # normalisation is scale-free: psi = u / sqrt(h) not needed here
            if l > 0:
                r, m = self.inj[l - 1]
                Nf = self.L[l - 1]["N"]; Nc = self.L[l]["N"]
                s = slice(Nf - m * r, Nf + m * r + 1, r)
                blk = own[l - 1][:, s, s, s]
                p = Nc - m
                blkp = torch.nn.functional.pad(blk, (p, p, p, p, p, p))
                phi = self.Mi[l] * blkp + (1 - self.Mi[l]) * phi
            own.append(phi)
        full = [None] * n
        full[n - 1] = own[n - 1]
        for l in range(n - 2, -1, -1):
            pr = kc.normalize(prolong(full[l + 1], self.P[l]))
            full[l] = self.M[l] * own[l] + (1 - self.M[l]) * pr
        return full

    def energy_parts(self, full):
        E2, E4 = [], []
        for l, lv in enumerate(self.L):
            e2, e4, _ = densities(full[l], lv["h"])
            w = self.W[l] * lv["h"] ** 3
            E2.append((w * e2).sum()); E4.append((w * e4).sum())
        return E2, E4

    def energy(self, us):
        E2, E4 = self.energy_parts(self.fields(us))
        return sum(self.k2[l] * E2[l] + self.k4[l] * E4[l] for l in range(len(self.L)))


def relax(H, us, maxit=100000, tol=1e-11, wall=None, hist=None, log=None, every=5):
    """L-BFGS on u_l (phi_l = u_l/|u_l|). Variables are stored with |u_l| = sqrt(h_l): since the energy is
    0-homogeneous in u, d E/d u scales like 1/|u|, which equalises the per-level Hessian scale (~h) seen by L-BFGS."""
    import time
    t0 = time.time()
    us = [u.detach().clone().requires_grad_(True) for u in us]
    opt = torch.optim.LBFGS(us, lr=1, max_iter=20, history_size=30, line_search_fn="strong_wolfe",
                            tolerance_grad=1e-14, tolerance_change=1e-16)
    hist = list(hist or [])

    def clo():
        opt.zero_grad(); E = H.energy(us); E.backward(); return E
    done = False
    while len(hist) < maxit:
        opt.step(clo)
        with torch.no_grad():
            for l, u in enumerate(us):
                u.data = kc.normalize(u.data) * math.sqrt(H.L[l]["h"])
            hist.append(float(H.energy(us)))
        if log and len(hist) % every == 0:
            log(f"    lbfgs {len(hist):5d}  E/12pi^2 = {hist[-1]/(12*math.pi**2):.9f}  t {time.time()-t0:.0f}s")
        if len(hist) > 10 and abs(hist[-6] - hist[-1]) / hist[-1] < 5 * tol:
            done = True; break
        if wall and time.time() - t0 > wall:
            break
    return [u.detach() for u in us], hist, done


# ------------------------------------------------------------------------------------------------ measurements
def stress_central(phi, h):
    """stress tensor T_ij from 4th-order central derivatives (as press_lattice.py)."""
    g = kc.derivs(phi, h)
    Dm = torch.stack([torch.stack([(g[i] * g[j]).sum(0) for j in range(3)]) for i in range(3)])
    trD = Dm[0, 0] + Dm[1, 1] + Dm[2, 2]
    D2 = torch.einsum("ik...,kj...->ij...", Dm, Dm)
    trD2 = D2[0, 0] + D2[1, 1] + D2[2, 2]
    e2 = trD; e4 = 0.5 * (trD**2 - trD2)
    T = 2 * Dm + 2 * trD * Dm - 2 * D2
    for i in range(3):
        T[i, i] -= e2 + e4
    return T, e2, e4


@torch.no_grad()
def measure(H, us, ode_npz, rbin=0.1, rmax=4.7):
    full = H.fields(us)
    O = np.load(ode_npz)
    p0_ode = float(O["p"][0]); r0_ode = float(O["r0"])
    out = {"levels": H.L, "k2": H.k2, "k4": H.k4, "sites": H.sites}
    E2s, E4s = H.energy_parts(full)
    nL = len(H.L)
    E2 = sum(H.k2[l] * float(E2s[l]) for l in range(nL)); E4 = sum(H.k4[l] * float(E4s[l]) for l in range(nL))
    out.update(E=E2 + E4, E2=E2, E4=E4, E_bare=float(sum(E2s) + sum(E4s)),
               E2_level=[float(x) for x in E2s], E4_level=[float(x) for x in E4s])
    B = 0.0; Lam = 0.0; Ip_fb = Ipabs_fb = Ip_c = Ipabs_c = 0.0
    R_all, pc_all, pfb_all, e_all, w_all = [], [], [], [], []
    for l, lv in enumerate(H.L):
        h = lv["h"]; phi = full[l]; w = H.W[l]
        X, Y, Z = H.grids(l)
        rr = torch.sqrt(X**2 + Y**2 + Z**2)
        e2, e4, (gp, gm) = densities(phi, h)
        e2k, e4k = H.k2[l] * e2, H.k4[l] * e4
        pfb = (e4k - e2k) / 3
        T, _, _ = stress_central(phi, h)
        pc = (T[0, 0] + T[1, 1] + T[2, 2]) / 3
        B += float((w * kc.baryon_density(phi, h)).sum()) * h**3
        v = torch.stack([torch.zeros_like(phi[0]), -phi[2], phi[1], torch.zeros_like(phi[0])])
        gc = [0.5 * (gp[i] + gm[i]) for i in range(3)]
        Tk = (v * v).sum(0) * (1 + e2) - sum((v * gc[i]).sum(0) ** 2 for i in range(3))
        Lam += 2 * float((w * Tk).sum()) * h**3      # isorotation about axis 3, omega = 1 (bare couplings)
        Ip_fb += float((w * pfb).sum()) * h**3; Ipabs_fb += float((w * pfb.abs()).sum()) * h**3
        Ip_c += float((w * pc).sum()) * h**3; Ipabs_c += float((w * pc.abs()).sum()) * h**3
        sel = w > 0
        R_all.append(rr[sel].numpy()); pc_all.append(pc[sel].numpy()); pfb_all.append(pfb[sel].numpy())
        e_all.append((e2k + e4k)[sel].numpy()); w_all.append((w[sel] * h**3).numpy())
        if l == 0:
            N = lv["N"]
            out["p0_central"] = float(pc[N, N, N]); out["p0_fb"] = float(pfb[N, N, N])
            out["T00_0"] = float((e2k + e4k)[N, N, N])
            f = torch.arccos(phi[0].clamp(-1, 1))
            m = (rr > 0) & (rr < 0.45)
            out["a_slope"] = float(np.polyfit(rr[m].numpy(), np.pi - f[m].numpy(), 1)[0]) if int(m.sum()) > 3 else None
    out.update(B=B, Lam=Lam, vonLaue_fb=Ip_fb / Ipabs_fb, vonLaue_central=Ip_c / Ipabs_c,
               virial=(E4 - E2) / (E2 + E4))
    R = np.concatenate(R_all); pc = np.concatenate(pc_all); pfb = np.concatenate(pfb_all); ee = np.concatenate(e_all)
    wv = np.concatenate(w_all)
    edges = np.arange(0, rmax + 1e-9, rbin); idx = np.digitize(R, edges) - 1; nb = len(edges) - 1
    sel = (idx >= 0) & (idx < nb)
    cnt = np.bincount(idx[sel], weights=wv[sel], minlength=nb)
    avg = lambda v: np.bincount(idx[sel], weights=(v * wv)[sel], minlength=nb) / np.maximum(cnt, 1e-300)
    rb = avg(R); pcb = avg(pc); pfbb = avg(pfb); eb = avg(ee)
    pOb = avg(np.interp(R, O["r"], O["p"])); eOb = avg(np.interp(R, O["r"], O["T00"]))
    ok = cnt > 0
    core = ok & (rb <= r0_ode / 2)
    out["core_max_dev_central"] = float((np.abs(pcb - pOb)[core] / p0_ode).max())
    out["core_max_dev_fb"] = float((np.abs(pfbb - pOb)[core] / p0_ode).max())
    sg = np.sign(pcb[ok]); rk = rb[ok]; j = int(np.argmax(sg < 0))
    out["r0_lattice"] = float(np.interp(0, [pcb[ok][j], pcb[ok][j - 1]], [rk[j], rk[j - 1]]))
    out["profile"] = dict(r=rb[ok].tolist(), p_central=pcb[ok].tolist(), p_fb=pfbb[ok].tolist(),
                          p_ode=pOb[ok].tolist(), T00=eb[ok].tolist(), T00_ode=eOb[ok].tolist())
    return out

"""nucleo1s_run.py - production runs of the study NUCLEO_1S, judged against lab/out/nucleo1s_prereg.json (commit 2e72499).
Usage: timeout 590 python lab/nucleo1s_run.py <part>    part in: law decomp zoom agglo z1 hybrid
Each part writes its block into lab/out/nucleo1s_results.json (checkpoint: a finished part is not redone)."""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_core import *
from nucleo1s_theory import inner_beta, contact_K, bulk_B
from zz_core import H_EVS

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
RES = os.path.join(OUT, "nucleo1s_results.json")
PRE = json.load(open(os.path.join(OUT, "nucleo1s_prereg.json")))
T0 = time.time()


def load():
    return json.load(open(RES)) if os.path.exists(RES) else {}


def save(part, block):
    import fcntl
    with open(RES + ".lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        r = load()
        r[part] = block
        tmp = RES + ".tmp"
        json.dump(r, open(tmp, "w"), indent=1, default=float)
        os.replace(tmp, RES)


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


def rmax_of(n, Za):
    return 30 + 9 * n * n / Za


# ------------------------------------------------------------------------------------------------- law
def part_law():
    out = load().get("law", {})
    for Z in (20, 40, 60, 80, 92):
        Za = Z * ALPHA
        pz = PRE["predictions_law"]["values"][str(Z)]
        K, B, g = pz["K"], pz["B"], pz["gamma"]
        blk = out.get(str(Z), {})
        for n, kappa in ((1, -1), (2, -2), (3, -3)):
            key = f"{n},{kappa}"
            ref = dirac_binding(n, kappa, Za)
            rmax = rmax_of(n, Za)
            rec = blk.get(key, dict(a=[], err=[], dirac=ref))
            gr = 1 + rec["err"][-1] if rec["err"] else (1.3 if kappa == -1 else 1.05)
            for k in range(len(rec["a"]), 12):
                a = 0.5 / 2**k
                b, info = level(Mesh.uniform(a, rmax), n, kappa, Za, guess_rel=gr)
                gr = b / ref
                rec["a"].append(a); rec["err"].append(b / ref - 1)
                blk[key] = rec; out[str(Z)] = blk; save("law", out)          # checkpoint after every lattice
            a_s, e_s = np.array(rec["a"]), np.array(rec["err"])
            if kappa == -1:
                R = (e_s[-1] - K * a_s[-1]**(2 * g)) / (B * a_s[-1]**2)
                d = e_s - B * a_s**2
                p = np.log(d[-2] / d[-1]) / np.log(2)
                rec.update(R=R, p_contact=p, two_gamma=2 * g, L1=bool(abs(R - 1) <= 0.02), L2=bool(abs(p - 2 * g) <= 0.02),
                           law_at_lab=dict(a=0.5, K_term=K * 0.5**(2 * g), B_term=B * 0.25, measured=e_s[0]))
            else:
                p = np.log(e_s[-2] / e_s[-1]) / np.log(2)
                rec.update(p=p, L3=bool(abs(p - 2) <= 0.01))
            blk[key] = rec; out[str(Z)] = blk; save("law", out)
            log(f"Z={Z} n={n} kappa={kappa}: err(0.5)={e_s[0]:+.4e} err(finest)={e_s[-1]:+.3e} " +
                (f"R={rec['R']:.5f} p={rec['p_contact']:.4f} (2g={2 * g:.4f})" if kappa == -1 else f"p={rec['p']:.5f}"))
    zs = [z for z in out if z.isdigit()]
    out["L1_pass"] = all(out[z]["1,-1"]["L1"] for z in zs)
    out["L2_pass"] = all(out[z]["1,-1"]["L2"] for z in zs)
    out["L3_pass"] = all(out[z][k]["L3"] for z in zs for k in ("2,-2", "3,-3"))
    save("law", out)


# ------------------------------------------------------------------------------------------------- decomposition
def part_decomp():
    Za = 60 * ALPHA
    ref = dirac_binding(1, -1, Za)
    rmax = rmax_of(1, Za)
    a = 0.5
    out = {}
    b0, info0 = level(Mesh.uniform(a, rmax), 1, -1, Za)
    e0 = b0 / ref - 1
    out["lab"] = e0
    # (1) the asymptotic law at a = 0.5
    pz = PRE["predictions_law"]["values"]["60"]
    out["law"] = dict(K_term=pz["K"] * a**pz["two_gamma"], B_term=pz["B"] * a * a,
                      rest=e0 - pz["K"] * a**pz["two_gamma"] - pz["B"] * a * a)
    # bulk split by ingredient (BCH, analytic continuation): only mass / kappa / V kept in A
    out["B_parts"] = {p: bulk_B(Za, parts=(p,)) for p in ("m", "k", "V")}
    out["B_all"] = bulk_B(Za)
    # (2) regions: cells inside R refined x256 (exact for that region), the rest 0.5
    reg = {}
    for R in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, rmax + 1):
        m = Mesh.graded(a, rmax, [(R, 256)])
        b, _ = level(m, 1, -1, Za, guess_rel=1.1)
        reg[str(R)] = b / ref - 1
    out["regions"] = reg
    # (3) the first cell as a mirror: reflection phase seen from r = 0.5, coarse cell vs fine cells (at w_1s Dirac)
    w1 = 1 - ref
    ph = {}
    for h in (0.5, 0.5 / 2**12):
        ein = np.arange(int(round(a / h)) + 1) * h
        rho = lambda w: inner_reflection(ein, -1, Za, w, a)
        ph[str(h)] = dict(phase=float(np.angle(rho(w1))),
                          delay=float(np.angle(rho(w1 + 1e-6) / rho(w1 - 1e-6)) / 2e-6))
    out["first_cell_mirror"] = ph
    dphi = ph["0.5"]["phase"] - ph[str(0.5 / 2**12)]["phase"]
    # first-order estimate of the level shift from the mirror phase error: dw = -dphi |x_in|^2 / (sum |x|^2 <l>)
    gw = info0["gw"]
    x = info0["x"]
    wgt = np.abs(x[gw.src_p[0]])**2 / np.sum(np.abs(x)**2) / gw.link_weight(x)
    out["first_cell_estimate"] = dict(dphi=dphi, dw=-dphi * wgt, rel_binding=dphi * wgt / ref)
    # (4) ingredients of the first cell (non additive; for the record)
    def run(M0=None, refl0=-1.0):
        gw2 = GradedWalk(Mesh.uniform(a, rmax), -1, Za, refl0=refl0)
        if M0 is not None:
            gw2.M[0] = M0
        w, _, _ = gw2.eigen(1 - 1.3 * ref)
        return (1 - w) / ref - 1

    def rot1(kr, V, m=1.0):
        hx, hz = kr, m
        nrm = np.hypot(hx, hz)
        phi = a * nrm
        c, s = np.cos(phi), np.sin(phi)
        M = np.array([[c - 1j * s * hz / nrm, -1j * s * hx / nrm], [-1j * s * hx / nrm, c + 1j * s * hz / nrm]]) * np.exp(-1j * V * a)
        return U_CHI.conj().T @ M @ U_CHI
    out["first_cell_variants"] = {
        "point (rule)": run(),
        "Coulomb off in cell 0": run(rot1(-1 / 0.25, 0.0)),
        "Coulomb shell-averaged (1.5/a)": run(rot1(-1 / 0.25, -Za * 1.5 / a)),
        "kappa/r shell-averaged": run(rot1(-1.5 / a, -Za / 0.25)),
        "both shell-averaged": run(rot1(-1.5 / a, -Za * 1.5 / a)),
        "other boundary at r=0 (f(0)=0)": run(refl0=+1.0)}
    # (5) D1 doubling: alias of the f(0)=0 level at w + pi/a in the g(0)=0 walk
    gwp = GradedWalk(Mesh.uniform(a, rmax), -1, Za, refl0=+1.0)
    wp, xp, _ = gwp.eigen(1 - 1.3 * ref)
    gwm = GradedWalk(Mesh.uniform(a, rmax), -1, Za, refl0=-1.0)
    wa, xa, _ = gwm.eigen(wp + np.pi / a + 1e-3)
    stag = np.repeat((-1.0)**np.arange(gwm.m.N), 2)
    ov = abs(np.vdot(xa, stag * xp)) / np.linalg.norm(xa) / np.linalg.norm(xp)
    d1 = abs((wa - np.pi / a) - wp) / abs(wp)
    out["D1"] = dict(w_f0=wp, w_alias_minus_pi_over_a=wa - np.pi / a, rel_diff=d1, overlap_staggered=ov,
                     w_1s=1 - b0, PASS=bool(d1 <= 1e-10))
    log("decomp:", json.dumps(out, default=float)[:1500])
    save("decomp", out)


# ------------------------------------------------------------------------------------------------- zoom
STATES = [(1, -1), (2, -1), (2, -2), (3, -1), (3, -2), (3, -3), (4, -1), (4, -2), (4, -3), (4, -4)]
NAMES = {(1, -1): "1s1/2", (2, -1): "2s1/2", (2, -2): "2p3/2", (3, -1): "3s1/2", (3, -2): "3p3/2", (3, -3): "3d5/2",
         (4, -1): "4s1/2", (4, -2): "4p3/2", (4, -3): "4d5/2", (4, -4): "4f7/2"}


def extend(edges, a, rmax):
    e = list(edges)
    while e[-1] < rmax - 1e-12:
        e.append(e[-1] + a)
    return np.array(e)


def nth_level(edges, n, kappa, Za):
    """level n of partial wave kappa: listen at the Dirac place and, if the rule's level found has the wrong
    number of nodes in g, listen again (the n - l - 1 count is the identification, as in orbitali.rs)"""
    ref = dirac_binding(n, kappa, Za)
    for gr in (1.3 if kappa == -1 else 1.05, 1.0, 1.15, 1.5, 0.9):
        b, info = level(Mesh(edges), n, kappa, Za, guess_rel=gr)
        g, f = info["gw"].gf(info["x"])
        gg = g.real
        big = np.abs(gg) > 1e-3 * np.abs(gg).max()   # (was 1e-6: tail noise counted as nodes, fixed after run 1)
        s = np.sign(gg[big])
        nodes = int(np.sum(s[1:] != s[:-1]))
        l = -kappa - 1 if kappa < 0 else kappa
        if nodes == n - l - 1:
            return b, nodes
    return b, nodes


def uniform_cells_for(err_target, Za, rmax):
    """smallest number of equal cells whose 1s error is <= err_target (bisection on log a)"""
    ref = dirac_binding(1, -1, Za)
    lo, hi = np.log(1e-4), np.log(0.5)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        a = np.exp(mid)
        b, _ = level(Mesh.uniform(a, rmax), 1, -1, Za, guess_rel=1.0 + err_target)
        if abs(b / ref - 1) <= err_target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-3:
            break
    a = np.exp(lo)
    return int(np.ceil(rmax / a)), a


def part_zoom():
    Za = 60 * ALPHA
    rmax1 = rmax_of(1, Za)
    res = load().get("zoom", {})
    ref = {s: dirac_binding(s[0], s[1], Za) for s in STATES}
    if "uniform_lab" not in res:
        u = {}
        for s in STATES:
            b, nd = nth_level(Mesh.uniform(0.5, rmax_of(s[0], Za)).e, s[0], s[1], Za)
            u[NAMES[s]] = b / ref[s] - 1
        res["uniform_lab"] = u
        save("zoom", res)
        log("uniform lab levels", u)
    for budget in (0.005, 0.001, 0.0002):
        key = f"budget_{budget}"
        if key in res:
            continue
        z = AutoZoom(0.5, rmax1, [(1, -1)], Za, budget, log=log).run()
        e1 = float(z.E[0] / ref[(1, -1)] - 1)
        est = z.history[-1]["est_err"][0]
        cells = len(z.edges) - 1
        nu, au = uniform_cells_for(abs(e1), Za, rmax1)
        lv = {}
        for s in STATES:
            ed = extend(z.edges, 0.5, rmax_of(s[0], Za))
            b, nd = nth_level(ed, s[0], s[1], Za)
            lv[NAMES[s]] = dict(err=b / ref[s] - 1, nodes=nd)
        U = res["uniform_lab"]
        z3 = {k: bool(lv[k]["err"] <= U[k] + 5e-4) for k in ("2p3/2", "3d5/2", "4f7/2")}
        res[key] = dict(err_1s=e1, est=est, ratio_est_true=est / abs(e1), cells=cells, h_min=float(np.diff(z.edges).min()),
                        uniform_cells_same_err=nu, uniform_a=au, cost_ratio=cells / nu, levels=lv,
                        history=z.history, edges=z.edges.tolist(),
                        Z1=bool(abs(e1) <= 0.01), Z2=bool(cells <= nu / 3), Z3=z3, Z3_all=all(z3.values()),
                        Z4=bool(0.5 <= est / abs(e1) <= 2))
        log(f"budget {budget}: 1s {100 * e1:+.4f} %, est {100 * est:.4f} %, cells {cells} vs uniform {nu} (a={au:.5f})",
            {k: f"{100 * v['err']:+.3f}%" for k, v in lv.items()})
        save("zoom", res)


# ------------------------------------------------------------------------------------------------- agglomerated nucleus
def part_agglo():
    Za = 60 * ALPHA
    zr = load()["zoom"]["budget_0.001"]
    edges = np.array(zr["edges"])
    R = 0.5
    ein = edges[edges <= R + 1e-12]
    eout = edges[edges >= R - 1e-12]
    a_first = eout[1] - eout[0]
    out = {"R": R, "cells_inside": len(ein) - 1, "cells_outside_1s": len(eout) - 1}
    learned = {}
    for kappa in (-1, -2, -3, -4):
        rho = lambda w, k=kappa: inner_reflection(ein, k, Za, w, a_first)
        ws = np.linspace(0.85, 1.0, 16)
        ph = np.unwrap(np.angle([rho(v) for v in ws]))
        tau, d0 = np.polyfit(ws, ph, 1)
        learned[kappa] = (d0, tau, rho, float(np.max(np.abs(ph - (d0 + tau * ws)))))
    out["learned"] = {str(k): dict(d0=v[0], tau=v[1], fit_resid_rad=v[3]) for k, v in learned.items()}
    lv = {}
    allok = True
    for s in STATES:
        n, kappa = s
        if kappa not in learned:
            continue
        ref = dirac_binding(n, kappa, Za)
        eo = extend(eout, 0.5, rmax_of(n, Za))
        full = np.concatenate([ein, eo[1:]])
        bf, _ = nth_level(full, n, kappa, Za)
        d0, tau, rho, _ = learned[kappa]
        ow = OuterWalk(Mesh(eo), kappa, Za, lambda w: np.exp(1j * (d0 + tau * w)))
        w2, _, _ = ow.eigen(1 - bf)
        ow1 = OuterWalk(Mesh(eo), kappa, Za, rho)
        w1, _, _ = ow1.eigen(1 - bf)
        d = abs((1 - w2) / bf - 1)
        ok = bool(d <= 5e-4)
        allok &= ok
        lv[NAMES[s]] = dict(with_cells=bf / ref - 1, exact_reflection=(1 - w1) / ref - 1, two_numbers=(1 - w2) / ref - 1,
                            rel_diff=d, N1=ok)
        log(NAMES[s], lv[NAMES[s]])
    out["levels"] = lv
    out["N1_pass"] = allok
    save("agglo", out)


# ------------------------------------------------------------------------------------------------- Z = 1
def part_z1():
    Za = ALPHA
    A0 = 1 / ALPHA
    rmax = 45 * A0
    MU = 1 / (1 + 1 / MP_ME)
    res = load().get("z1", {"levels": {}})
    for th in (0.05, 0.025, 0.0125):
        if str(th) in res["levels"]:
            continue
        m = mesh_zoom(th, rmax, 14, 16)
        lv = {}
        for n in (1, 2, 3):
            ref = dirac_binding(n, -1, Za)
            b, info = level(m, n, -1, Za, guess_rel=1 + th * th / 3)
            lv[n] = b
            log(f"theta {th} {n}s: rel err {b / ref - 1:+.4e} cells {m.N}")
        res["levels"][str(th)] = {str(k): v for k, v in lv.items()}
        res["cells"] = res.get("cells", {}); res["cells"][str(th)] = m.N
        save("z1", res)
    th = np.array([0.05, 0.025, 0.0125])
    ext = {}
    for n in (1, 2, 3):
        y = np.array([res["levels"][str(t)][str(n)] for t in th])
        c = np.polyfit(th**2, y, 2)
        ext[n] = c[-1]
    rel = {f"{n}s": ext[n] / dirac_binding(n, -1, Za) - 1 for n in (1, 2, 3)}
    res["extrapolated_rel_err"] = rel
    res["H1z"] = all(abs(v) <= 1e-5 for v in rel.values())
    res["H1z_1e-7"] = all(abs(v) <= 1e-7 for v in rel.values())
    f_rule = (ext[1] - ext[2]) * MC2_EV * MU / H_EVS
    f_dirac = (dirac_binding(1, -1, Za) - dirac_binding(2, -1, Za)) * MC2_EV * MU / H_EVS
    f_meas = 2466061413187035.0
    res["1S2S"] = dict(rule_Hz=f_rule, dirac_reduced_Hz=f_dirac, lattice_err_GHz=(f_rule - f_dirac) / 1e9,
                       vs_measured_GHz=(f_rule - f_meas) / 1e9, dirac_vs_measured_GHz=(f_dirac - f_meas) / 1e9)
    res["H2z"] = bool(abs(f_rule - f_dirac) <= 1e9)
    res["E1s_reduced_eV"] = ext[1] * MC2_EV * MU
    log("z1:", rel, res["1S2S"])
    save("z1", res)


# ------------------------------------------------------------------------------------------------- hybrid
def part_hybrid():
    from nucleo1s_hybrid import CoarseRadial, Hybrid, coarse_only, A0
    MU = 1 / (1 + 1 / MP_ME)
    out = {}

    def solve(c, Z, n, J0, rmax_a0):
        Za = Z * ALPHA
        Eref = -dirac_binding(n, -1, Za)
        hy = Hybrid(c, rmax_a0 * A0 / Z, J0=J0, Za=Za)
        Es = np.linspace(Eref * 1.06, Eref * 0.94, 25)
        rs = [hy.residual(E) for E in Es]
        br = [i for i in range(24) if np.sign(rs[i]) != np.sign(rs[i + 1])]
        if not br:
            return None, hy
        # bracket closest to the Dirac place (the rule decides inside it)
        i = min(br, key=lambda i: abs(0.5 * (Es[i] + Es[i + 1]) - Eref))
        E, u = hy.solve(Es[i], Es[i + 1])
        return E, hy

    for B, vmode in ((256, "point"), (512, "point"), (256, "shell")):
        c = CoarseRadial(B, vmode=vmode)
        key = f"B{B}_{vmode}"
        blk = dict(h=c.h, m_kin=c.m_kin, spin_part=c.spin_part, aniso_k_1_over_a0=c.aniso_1a0)
        for Z in ((1, 2) if (B == 256 and vmode == "point") else (1,)):
            for n in (1, 2):
                Za = Z * ALPHA
                Eref = -dirac_binding(n, -1, Za)
                rmax_a0 = 24 if n == 1 else 40
                Ec = coarse_only(c, rmax_a0 * A0 / Z, n, Za)
                rec = dict(dirac_eV=Eref * MC2_EV, coarse_only=Ec / Eref - 1)
                for J0 in (2, 3, 4):
                    E, hy = solve(c, Z, n, J0, rmax_a0)
                    rec[f"J0={J0}"] = None if E is None else E / Eref - 1
                    rec[f"J0={J0}_E_eV"] = None if E is None else E * MC2_EV
                    rec["cells_fine"], rec["cells_coarse"] = hy.cells_fine, hy.N
                blk[f"Z{Z}_{n}s"] = rec
                log(key, f"Z{Z} {n}s", {k: (f"{100 * v:+.4f}%" if isinstance(v, float) and abs(v) < 1 else v) for k, v in rec.items()})
        out[key] = blk
    p = out["B256_point"]
    e1 = p["Z1_1s"]["J0=3_E_eV"] * MU
    Ered = -dirac_binding(1, -1, ALPHA) * MC2_EV * MU
    out["HY1"] = dict(E1s_reduced_eV=e1, dirac_reduced_eV=Ered, rel=e1 / Ered - 1, vs_13_6057=e1 / -13.6057 - 1,
                      PASS=bool(abs(e1 / Ered - 1) <= 1e-3))
    e2 = out["B512_point"]["Z1_2s"]["J0=3"]
    out["HY2"] = dict(rel=e2, PASS=bool(e2 is not None and abs(e2) <= 1e-2))
    h3 = [p["Z2_1s"]["J0=3"], p["Z2_2s"]["J0=3"]]
    out["HY3"] = dict(rel_1s=h3[0], rel_2s=h3[1], PASS=bool(all(v is not None and abs(v) <= 5e-3 for v in h3)))
    log("HY", out["HY1"], out["HY2"], out["HY3"])
    save("hybrid", out)


if __name__ == "__main__":
    part = sys.argv[1]
    {"law": part_law, "decomp": part_decomp, "zoom": part_zoom, "agglo": part_agglo, "z1": part_z1,
     "hybrid": part_hybrid}[part]()
    log("done", part)

"""motore_f2_run.py - phase 2 runs (pre-registration lab/out/motore_f2_prereg.json).
Usage: python lab/motore_f2_run.py atom intel|ultra   (checkpointed: re-run until it prints DONE)
       python lab/motore_f2_run.py lot | pos | clip | judge
Writes lab/out/motore_f2_<part>.json (vectors of the atom in lab/out/motore_f2_vec_*.npy, not committed)."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_f2_core as f2
import motore_core as mc

PRE = json.load(open(os.path.join(f2.OUT, "motore_f2_prereg.json")))
T0 = time.time()
J = lambda p: os.path.join(f2.OUT, "motore_f2_%s.json" % p)


def load(p):
    return json.load(open(J(p))) if os.path.exists(J(p)) else {}


def part_atom(tier):
    res = load("atom")
    R = res.setdefault(tier, {})
    B = 256 if tier == "ultra" else 512
    n = PRE["tiers"][tier]["octant"]
    M = PRE["configs"]["M_" + tier]
    q = PRE["configs"]["q_a0"]
    lev = f2.load_level3(B)
    R.update(B=B, octant=n, box_a0=n * lev["h"] / f2.A0, cell_a0=lev["h"] / f2.A0, stencil_asym=lev["asym"])
    ats = {"even": f2.OctAtom(lev, n, "even"), "odd": f2.OctAtom(lev, n, "odd")}
    vec = {}
    for kind, sec, Yk in (("1s", "even", None), ("2s", "even", "1s"), ("2p", "odd", None)):
        fn = os.path.join(f2.OUT, f"motore_f2_vec_{tier}_{kind}.npy")
        if os.path.exists(fn) and kind in R:
            vec[kind] = np.load(fn)
            continue
        if time.time() - T0 > 400:
            json.dump(res, open(J("atom"), "w"), indent=1); print("again"); return
        t0 = time.time()
        at = ats[sec]
        Y = vec[Yk][:, None] if Yk else None
        v, lres, it = f2.solve(at, kind, Y=Y, maxiter=80)
        E, rres, lam = at.energy(v)
        np.save(fn, v)
        vec[kind] = v
        R[kind] = dict(E_eV=E * f2.MC2, rule_resid=rres, lobpcg_resid=lres, iterations=it, seconds=time.time() - t0)
        print(tier, kind, R[kind], flush=True)
        json.dump(res, open(J("atom"), "w"), indent=1)
    # kick
    X, Y_, Z = ats["even"].grid()
    qz = q * Z / f2.A0
    psi1 = vec["1s"].reshape((n,) * 3)
    ev = np.cos(qz) * psi1
    od = 1j * np.sin(qz) * psi1
    c1 = complex(np.vdot(vec["1s"], ev.ravel()))
    c2s = complex(np.vdot(vec["2s"], ev.ravel()))
    c2p = complex(np.vdot(vec["2p"], od.ravel()))
    norm0 = float(np.sum(np.abs(ev) ** 2) + np.sum(np.abs(od) ** 2))
    R["kick"] = dict(q_a0=q, norm=norm0, P1s=abs(c1) ** 2, P2s=abs(c2s) ** 2, P2p=abs(c2p) ** 2,
                     truncated=1 - (abs(c1) ** 2 + abs(c2s) ** 2 + abs(c2p) ** 2) / norm0)
    line = R["2p"]["E_eV"] - R["1s"]["E_eV"]
    R["line_eV"] = line
    print(tier, "kick", R["kick"], "line", line, flush=True)
    json.dump(res, open(J("atom"), "w"), indent=1)
    # evolution under the rule: conservation + cartolina check (G1c, G2)
    if "evolve" in R:
        print("DONE"); return
    lam = {k: complex(f2.OctAtom.energy(ats["even" if k != "2p" else "odd"], vec[k])[2]) for k in ("1s", "2s", "2p")}
    t0 = time.time()
    fe, fo = ev.astype(complex), od.astype(complex)
    for s in range(M):
        fe = ats["even"].W(fe)
        fo = ats["odd"].W(fo)
    t_grid = (time.time() - t0) / M
    normM = float(np.sum(np.abs(fe) ** 2) + np.sum(np.abs(fo) ** 2))
    d = {}
    for k, cc, f in (("1s", c1, fe), ("2s", c2s, fe), ("2p", c2p, fo)):
        d[k] = abs(complex(np.vdot(vec[k], f.ravel())) - cc * lam[k] ** M)
    rmax = max(R[k]["rule_resid"] for k in ("1s", "2s", "2p"))
    t1 = time.time()
    for s in range(100000):
        z = c1 * lam["1s"] ** s + c2s * lam["2s"] ** s + c2p * lam["2p"] ** s
    t_card = (time.time() - t1) / 100000
    R["evolve"] = dict(M=M, norm_change=normM / norm0 - 1, cartolina_diff=d, bound=M * rmax * 2 + 1e-12,
                       G2_pass=all(v <= M * rmax * 2 + 1e-12 for v in d.values()),
                       sec_per_step_grid=t_grid, sec_per_step_cartolina=t_card, cost_ratio=t_grid / t_card)
    print(tier, "evolve", R["evolve"], flush=True)
    json.dump(res, open(J("atom"), "w"), indent=1)
    print("DONE")


def part_lot():
    import reticolo_core as rc
    out = dict(light=[], locality=None, electron=[])
    L, e = PRE["configs"]["lot"]["L"], PRE["configs"]["lot"]["e"]
    for Nc in (32, 16, 12, 8):
        Nf = 2 * Nc
        nf = int(round(27 * Nf)); nf = int(np.ceil(nf / L) * L)
        Ff = rc.packet((nf, nf), 2, Nf, 30.0, [nf / 2, nf / 2], 3 * Nf)
        Ef = mc.energy(Ff)
        Fc = f2.lot_restrict(Ff, L, e, (0, 1))
        Fb = f2.lot_prolong(Fc, L, e, (0, 1))
        o = mc.overlap(Ff, Fb)
        row = dict(N_coarse=Nc, dE_restrict=4 * mc.energy(Fc) / Ef - 1, dE_roundtrip=mc.energy(Fb) / Ef - 1,
                   fidelity=float(abs(o) ** 2), centroid_jump_lambda=float(np.linalg.norm(mc.centroid(Fc) * 2 / Nf - mc.centroid(Ff) / Nf)),
                   wrong_band_coarse=1 - rc.helicity_fractions(Fc, 2)[0])
        out["light"].append(row)
        print("G3a", row, flush=True)
    # locality: perturb one fine block (2-D), look where the coarse field changes
    n = 256
    rng = np.random.default_rng(1)
    base = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    pert = base.copy()
    b0 = 4 * L
    pert[b0:b0 + L, b0:b0 + L] += rng.normal(size=(L, L))
    dC = np.abs(f2.lot_restrict(pert, L, e, (0, 1)) - f2.lot_restrict(base, L, e, (0, 1)))
    xs = np.arange(n // 2)
    lo, hi = (b0 - e) // 2, (b0 + L + e) // 2           # coarse index range that may change
    inside = (xs >= lo) & (xs < hi)
    outside_max = float(dC[~inside, :].max() if (~inside).any() else 0)
    outside_max = max(outside_max, float(dC[:, ~inside].max()))
    out["locality"] = dict(block_fine=[b0, b0 + L], allowed_coarse=[int(lo), int(hi)], max_change_outside=outside_max,
                           max_change_inside=float(dC.max()))
    print("G3b", out["locality"], flush=True)
    # electron charge: fine 1-D walk (theta = a0/4096), positive packet, 4 local halvings
    th, lv = mc.load_levels_1d()
    nf = 196608
    xf = (np.arange(nf) - nf / 2 + 0.5) * th
    k = 2 * np.pi * np.fft.fftfreq(nf)
    Ekin = np.arccos(np.cos(th) * np.cos(k)) / th - 1.0
    for k0a in (1.0, 2.0):
        psi = mc.positive_spinor_packet(xf, k0a / mc.A0, mc.A0 / 2, th)
        N0 = np.sum(np.abs(psi) ** 2)
        c = psi
        for _ in range(4):
            c = f2.lot_restrict(c, L, e, (1,))
        b = c
        for _ in range(4):
            b = f2.lot_prolong(b, L, e, (1,))
        fb, f0 = np.fft.fft(b, axis=-1), np.fft.fft(psi, axis=-1)
        Eb = np.sum(np.abs(fb) ** 2 * Ekin) / np.sum(np.abs(fb) ** 2)
        E0 = np.sum(np.abs(f0) ** 2 * Ekin) / np.sum(np.abs(f0) ** 2)
        row = dict(k0_a0=k0a, dN_coarse=float(16 * np.sum(np.abs(c) ** 2) / N0 - 1), dN_roundtrip=float(np.sum(np.abs(b) ** 2) / N0 - 1),
                   dEkin_roundtrip=float(Eb / E0 - 1), max_err=float(np.abs(b - psi).max() / np.abs(psi).max()))
        out["electron"].append(row)
        print("G3c", row, flush=True)
    sel = [r for r in out["light"] if r["N_coarse"] >= 12]
    out["G3a"] = "PASS" if all(abs(r["dE_restrict"]) <= 1e-6 and abs(r["dE_roundtrip"]) <= 1e-6 for r in sel) else "FAIL"
    out["G3b"] = "PASS" if outside_max <= 1e-14 else "FAIL"
    out["G3c"] = "PASS" if all(abs(r["dN_coarse"]) <= 1e-6 and abs(r["dN_roundtrip"]) <= 1e-6 and abs(r["dEkin_roundtrip"]) <= 1e-6 for r in out["electron"]) else "FAIL"
    print(out["G3a"], out["G3b"], out["G3c"])
    json.dump(out, open(J("lot"), "w"), indent=1)


def band_projector(rule, nc):
    """positive-band projector of a learned 1-D (R, L) rule on an nc grid: eigenvector with the larger weight on
    phi = (R + L)/sqrt 2 (as agg_zz_1d.dispersion)"""
    K = 2 * np.pi * np.fft.fftfreq(nc)
    M = rule.mult_at(K[:, None])
    lam, V = np.linalg.eig(M)
    wphi = np.abs(V[:, 0, :] + V[:, 1, :]) ** 2 / 2
    ip = np.argmax(wphi, axis=1)
    sel = np.zeros((nc, 2))
    sel[np.arange(nc), ip] = 1
    return np.einsum("nab,nb,nbc->nac", V, sel, np.linalg.inv(V))


def apply_k1(P, psi):
    f = np.fft.fft(psi, axis=-1).T
    return np.fft.ifft(np.einsum("nab,nb->na", P, f).T, axis=-1)


def part_pos():
    th, lv = mc.load_levels_1d()
    B = 192
    rule = lv[B]
    nf = 196608
    nc = nf // B
    xf = (np.arange(nf) - nf / 2 + 0.5) * th
    xc = (np.arange(nc) - nc / 2 + 0.5) * th * B
    tau = th * B
    P = band_projector(rule, nc)
    out = {}
    psi = mc.positive_spinor_packet(xf, 1 / mc.A0, mc.A0 / 2, th)
    c0 = mc.block1(psi, B)
    cp = apply_k1(P, c0)
    neg = float(np.sum(np.abs(c0 - cp) ** 2) / np.sum(np.abs(c0) ** 2))
    nst = 400
    alias = (2.0 * tau) % (2 * np.pi)
    alias = min(alias, 2 * np.pi - alias) / tau

    def run(c):
        PR, leak = [], 0.0
        s = c.copy()
        for i in range(nst):
            PR.append(float(np.sum(np.abs(s[0]) ** 2) / np.sum(np.abs(s) ** 2)))
            s2 = rule.step(s[None], 1)[0]
            leak = max(leak, float(np.sum(np.abs(s2 - apply_k1(P, s2)) ** 2) / np.sum(np.abs(s2) ** 2)))
            s = s2
        PR = np.array(PR)
        t = np.arange(nst) * tau
        amp = 2 * abs(np.mean((PR - PR.mean()) * np.exp(1j * alias * t)))
        return PR, amp, leak, s
    PR0, amp0, leak0, s0 = run(c0)
    PR1, amp1, leak1, s1 = run(cp)
    m0, w0, _ = mc.x_mean_width(s0, xc)
    m1, w1, _ = mc.x_mean_width(s1, xc)
    out["atomic_packet"] = dict(negative_weight_after_block=neg, alias_omega=alias, fake_beat_amp_unprojected=amp0,
                                fake_beat_amp_projected=amp1, leak_per_step_projected=leak1,
                                d_centroid_rel=abs(m1 - m0) / abs(m0) if m0 else abs(m1 - m0),
                                d_width_rel=abs(w1 / w0 - 1))
    # pure R at rest: must be flagged
    g = np.exp(-xf ** 2 / (4 * (mc.A0 / 2) ** 2))
    cz = mc.block1(np.stack([g, 0 * g]).astype(complex), B)
    negz = float(np.sum(np.abs(cz - apply_k1(P, cz)) ** 2) / np.sum(np.abs(cz) ** 2))
    PRz, ampz, _, _ = run(cz)
    # valley hypothesis (coordinator, BATTIMENTO_1S): is the slow beat the Zitterbewegung aliased by the coarse step,
    # or a mixing of momentum valleys? 1-D fine walk: valleys k = 0 (omega = +-theta per tick) and k = pi
    # (omega = +-(pi - theta)). Candidates folded into the coarse Nyquist band:
    def fold(om_per_time):
        ph = (om_per_time * tau) % (2 * np.pi)
        return min(ph, 2 * np.pi - ph) / tau
    zb = fold(2.0)
    valley = sorted(set(round(fold(((np.pi - th) - s1 * th) / th), 12) for s1 in (+1, -1)))
    t = np.arange(nst) * tau
    sp = np.abs(np.fft.rfft((PRz - PRz.mean()) * np.hanning(nst), n=16 * nst))
    fr = np.fft.rfftfreq(16 * nst, d=tau) * 2 * np.pi
    peak = float(fr[np.argmax(sp[1:]) + 1])
    kpi = float(np.sum(np.abs(np.fft.fft(cz, axis=-1)[:, nc // 2 - 2: nc // 2 + 3]) ** 2) / np.sum(np.abs(np.fft.fft(cz, axis=-1)) ** 2))
    out["pure_R"] = dict(negative_weight=negz, flagged=negz > 1e-3, fake_beat_amp_unprojected=ampz,
                         beat_peak_omega=peak, zb_alias_omega=zb, valley_alias_omegas=valley,
                         weight_near_coarse_pi=kpi, verdict_valley=("ZB aliasata" if abs(peak - zb) < 0.1 * zb else
                                                                    ("valli" if any(abs(peak - v) < 0.1 * v for v in valley) else "nessuna delle due")))
    a = out["atomic_packet"]
    out["G4"] = "PASS" if (a["fake_beat_amp_projected"] <= 1e-10 and a["leak_per_step_projected"] <= 1e-12 and
                           a["d_centroid_rel"] <= 1e-6 and a["d_width_rel"] <= 1e-6 and out["pure_R"]["flagged"]) else "FAIL"
    print(json.dumps(out, indent=1))
    json.dump(out, open(J("pos"), "w"), indent=1)


def part_clip():
    """shells of levels around the eye (3-D): each level k covers a cube of side S*2^k fine cells... with a fixed number
    of cells per side m; levels nested (the finer inside). Compare with a full fine world and with the camino of R4."""
    rows = []
    for m in (64, 96, 128):
        for K in (10, 20, 40):
            shell = m ** 3 - (m // 2) ** 3                      # level k minus the hole filled by level k-1
            cells = m ** 3 + (K - 1) * shell
            local = m ** 3 + sum(shell / 2 ** k for k in range(1, K))
            full = float(m) ** 3 * 8.0 ** (K - 1)
            rows.append(dict(m=m, levels=K, cells=cells, mem_GB_light=cells * 24 * 2 / 1e9, cost_per_fine_tick_local=local,
                             cost_global=cells, full_world_cells=full, saving_vs_full=full / cells,
                             world_side_over_fine_side=2.0 ** (K - 1)))
    json.dump(rows, open(J("clip"), "w"), indent=1)
    for r in rows:
        print(r)


def judge():
    A = load("atom"); Lo = load("lot"); Po = load("pos")
    ex = PRE["predictions"]
    out = {}
    for tier in ("ultra", "intel"):
        R = A[tier]; tol = PRE["tiers"][tier]["tol"]
        line = R["line_eV"]; kick = R["kick"]
        g1a = abs(line / ex["line_exact_eV"] - 1) <= tol
        g1b = abs(kick["P2p"] / ex["exact_P2p"] - 1) <= tol and abs(kick["P1s"] / ex["exact_P1s"] - 1) <= tol
        g1c = abs(R["evolve"]["norm_change"]) <= 1e-10
        omega_rel = line / ex["line_exact_eV"] - 1
        zeros_ok = all(abs((k * np.pi / (1 + omega_rel)) - k * np.pi) <= tol * k * np.pi for k in range(1, 5))
        g1d = (R["1s"]["E_eV"] < R["2p"]["E_eV"] and R["1s"]["E_eV"] < R["2s"]["E_eV"] and
               kick["P1s"] > kick["P2p"] > kick["P2s"] and zeros_ok)
        out[tier] = dict(line_eV=line, line_rel=omega_rel, P2p=kick["P2p"], P2p_rel=kick["P2p"] / ex["exact_P2p"] - 1,
                         P1s=kick["P1s"], P1s_rel=kick["P1s"] / ex["exact_P1s"] - 1, P2s=kick["P2s"],
                         truncated=kick["truncated"], norm_change=R["evolve"]["norm_change"],
                         G1a=g1a, G1b=g1b, G1c=g1c, G1d=g1d, G2=R["evolve"]["G2_pass"], cost_ratio=R["evolve"]["cost_ratio"])
    out["G1e"] = (abs(A["intel"]["line_eV"] / A["ultra"]["line_eV"] - 1) <= 0.10 and
                  abs(A["intel"]["kick"]["P2p"] / A["ultra"]["kick"]["P2p"] - 1) <= 0.10)
    out["G3"] = [Lo.get("G3a"), Lo.get("G3b"), Lo.get("G3c")]
    out["G4"] = Po.get("G4")
    print(json.dumps(out, indent=1, default=str))
    json.dump(out, open(J("judge"), "w"), indent=1, default=str)


if __name__ == "__main__":
    p = sys.argv[1]
    if p == "atom":
        part_atom(sys.argv[2])
    else:
        {"lot": part_lot, "pos": part_pos, "clip": part_clip, "judge": judge}[p]()

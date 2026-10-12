"""kz_k3.py - K3: born defects move.  Inertial (Hamiltonian) XY dynamics theta'' = -gamma theta' + sum_j sin(theta_j - theta_i).
 a) 2D vortex-antivortex pair at rest (d = 24, L = 128, gamma = 0): approach, annihilation, energy -> outgoing wave.
 b) same pair Lorentz-boosted at v = 0.5 perpendicular to the separation (report only).
 c) 2D KZ quench (L = 128, tauQ = 256) frozen at T = 0.2, then inertial dynamics gamma = 0.02: vortex number vs time.
 d) 3D K1 final configuration (L = 48, tauQ = 1024), cooled 0.5 t_hat, inertial gamma = 0.02: line length vs time.
-> lab/out/kz_k3.json, lab/out/kz_k3_tracks.npz
"""
import sys, os, json, time, math
import numpy as np, torch
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from kz_core import (force, inertial, vortex_map_2d, vortex_count, wrap, cool, quench, TWO_PI)

OUT = os.path.join(HERE, "out")
TBKT = 0.8929   # 2D XY BKT temperature (literature)


def energy_density(th, p):
    e = 0.5 * p * p
    for d in range(th.dim()):
        b = 1 - torch.cos(torch.roll(th, -1, d) - th)
        e = e + 0.5 * (b + torch.roll(b, 1, d))     # bond energy shared between its two sites
    return e


def pair_field(L, x1, x2, y0, gam=1.0, K=8):
    x = torch.arange(L, dtype=torch.float64).view(-1, 1); y = torch.arange(L, dtype=torch.float64).view(1, -1)
    th = torch.zeros(L, L, dtype=torch.float64)
    for nx in range(-K, K + 1):
        for ny in range(-K, K + 1):
            Y = gam * (y - y0 + ny * L)
            th += torch.atan2(Y, x - x1 + nx * L) - torch.atan2(Y, x - x2 + nx * L)
    return wrap(th)


def positions(th):
    w = vortex_map_2d(th)
    P = (w > 0).nonzero().double() + 0.5; Mn = (w < 0).nonzero().double() + 0.5
    return P.numpy(), Mn.numpy()


def part_a(L=128, d=24, dt=0.05, tmax=300.0):
    x0 = y0 = L / 2 + 0.5   # pair between lattice sites (plaquette centres)
    th = pair_field(L, x0 - d / 2, x0 + d / 2, y0 - 0.0)
    th = cool(th, 3.0, dt=0.1)
    p = torch.zeros_like(th)
    E0 = energy_density(th, p).sum().item()
    rec = {"t": [], "xp": [], "xm": [], "E": [], "Ekin": [], "nv": []}
    t = 0.0; t_ann = None; shells = []
    X = torch.arange(L, dtype=torch.float64).view(-1, 1); Y = torch.arange(L, dtype=torch.float64).view(1, -1)
    nstep = int(round(0.5 / dt))
    while t < tmax:
        P, Mn = positions(th)
        e = energy_density(th, p)
        rec["t"].append(t); rec["E"].append(e.sum().item()); rec["Ekin"].append(0.5 * (p * p).sum().item())
        rec["nv"].append(len(P) + len(Mn))
        rec["xp"].append(P[:, 0].mean() if len(P) else np.nan); rec["xm"].append(Mn[:, 0].mean() if len(Mn) else np.nan)
        if t_ann is None and len(P) + len(Mn) == 0 and t > 1:
            t_ann = t; xa = 0.5 * (rec["xp"][-2] + rec["xm"][-2])
        if t_ann is not None:
            r = torch.sqrt((X - xa) ** 2 + (Y - y0) ** 2)
            tot = e.sum().item()
            rs, idx = torch.sort(r.flatten()); cs = torch.cumsum(e.flatten()[idx], 0) / tot
            r50 = rs[(cs >= 0.5).nonzero()[0, 0]].item(); r90 = rs[(cs >= 0.9).nonzero()[0, 0]].item()
            out15 = 1 - cs[(rs <= 15).nonzero()[-1, 0]].item()
            shells.append((t - t_ann, r50, r90, out15))
            if t - t_ann > 60:
                break
        th, p = inertial(th, p, nstep, dt=dt, gamma=0.0)
        t += nstep * dt
    tt = np.array(rec["t"]); xp = np.array(rec["xp"]); xm = np.array(rec["xm"])
    sep = xm - xp
    # speeds from separation, smoothed over 5 time units (lattice-quantised positions)
    ok = ~np.isnan(sep)
    ts, ss = tt[ok], sep[ok]
    k = 10
    v_each = np.abs((ss[k:] - ss[:-k]) / (ts[k:] - ts[:-k])) / 2
    sh = np.array(shells) if shells else np.zeros((0, 4))
    win = (sh[:, 0] >= 10) & (sh[:, 0] <= 50)
    v90 = np.polyfit(sh[win, 0], sh[win, 2], 1)[0] if win.sum() > 3 else np.nan
    v50 = np.polyfit(sh[win, 0], sh[win, 1], 1)[0] if win.sum() > 3 else np.nan
    Erel = (np.max(rec["E"]) - np.min(rec["E"])) / E0
    out15 = sh[sh[:, 0] >= 40, 3][0] if (sh[:, 0] >= 40).any() else np.nan
    mv = float(v_each.max()) if len(v_each) else None
    res = {"d": d, "L": L, "E_pair": E0, "t_annihilation": t_ann, "energy_drift_rel": Erel,
           "max_vortex_speed_smoothed": mv, "separation_end": float(ss[-1]) if len(ss) else None,
           "speed_profile": {"t": ts[k:].tolist()[::4], "v": v_each.tolist()[::4]},
           "front_speed_r90": v90, "front_speed_r50": v50, "energy_fraction_outside_r15_at_40_after": out15,
           "Ekin_fraction_end": rec["Ekin"][-1] / rec["E"][-1]}
    res["checks"] = {"annihilated_before_300": t_ann is not None and t_ann < 300, "energy_conserved_1e-3": Erel < 1e-3,
                     "max_speed_below_c": mv is not None and mv < 1,
                     "energy_out_ge_0.9": bool(out15 >= 0.9), "front_speed_1pm0.1": bool(abs(v90 - 1) <= 0.1)}
    res["PASS"] = bool(all(res["checks"].values()))
    return res, rec, sh


def part_b(L=128, d=24, v=0.5, dt=0.05, tmax=40.0):
    x0 = y0 = L / 2 + 0.5
    g = 1 / math.sqrt(1 - v * v)
    th = pair_field(L, x0 - d / 2, x0 + d / 2, y0, gam=g)
    dy = wrap(torch.roll(th, -1, 1) - torch.roll(th, 1, 1)) / 2
    p = -v * dy
    tr = []
    t = 0.0; nstep = int(round(0.5 / dt))
    while t <= tmax:
        P, Mn = positions(th)
        if len(P) and len(Mn):
            tr.append((t, P[:, 0].mean(), P[:, 1].mean(), Mn[:, 0].mean(), Mn[:, 1].mean()))
        th, p = inertial(th, p, nstep, dt=dt)
        t += nstep * dt
    tr = np.array(tr)
    yc = 0.5 * (tr[:, 2] + tr[:, 4]); sep = tr[:, 3] - tr[:, 1]
    w = tr[:, 0] <= 20
    vy = np.polyfit(tr[w, 0], yc[w], 1)[0]
    vy_late = np.polyfit(tr[tr[:, 0] >= 20, 0], yc[tr[:, 0] >= 20], 1)[0] if (tr[:, 0] >= 20).sum() > 3 else None
    return {"v_boost": v, "centroid_speed_t0-20": vy, "centroid_speed_t20-40": vy_late,
            "separation_start_end": [float(sep[0]), float(sep[-1])], "t_last_tracked": float(tr[-1, 0])}, tr


def part_c(L=128, tauQ=256, Tstop=0.2, gamma=0.02, dt=0.05, tmax=200.0):
    tstop = (1 - Tstop / TBKT) * tauQ
    r, th = quench(L, 2, tauQ, [tstop], 777, Tc=TBKT, cool_times=(0.0,), dtype=torch.float64, keep_final=True)
    p = torch.zeros_like(th)
    rec = []
    t = 0.0; nstep = int(round(5.0 / dt))
    while t <= tmax:
        n, _ = vortex_count(th)
        rec.append((t, n, 0.5 * (p * p).sum().item(), energy_density(th, p).sum().item()))
        th, p = inertial(th, p, nstep, dt=dt, gamma=gamma)
        t += nstep * dt
    rec = np.array(rec)
    return {"t_stop_after_Tc": tstop, "n_vortices_start": int(rec[0, 1]), "n_vortices_t200": int(rec[-1, 1]),
            "drop_fraction": 1 - rec[-1, 1] / rec[0, 1], "Ekin_max": float(rec[:, 2].max()),
            "E_start": float(rec[0, 3]), "E_end": float(rec[-1, 3]),
            "PASS": bool(rec[-1, 1] <= 0.5 * rec[0, 1] and rec[:, 2].max() > 0)}, rec


def part_d(gamma=0.02, dt=0.05, tmax=100.0):
    fn = os.path.join(OUT, "kz_final_L48_t1024.pt")
    if not os.path.exists(fn):
        return {"missing": fn}, None
    th = torch.load(fn).double()
    th = cool(th, 0.5 * 1024 ** 0.573)
    p = torch.zeros_like(th)
    rec = []; t = 0.0; nstep = int(round(10.0 / dt))
    while t <= tmax:
        n, _ = vortex_count(th)
        rec.append((t, n / 48 ** 3, 0.5 * (p * p).sum().item()))
        th, p = inertial(th, p, nstep, dt=dt, gamma=gamma)
        t += nstep * dt
    rec = np.array(rec)
    return {"line_density_t": rec[:, 1].tolist(), "t": rec[:, 0].tolist(), "Ekin": rec[:, 2].tolist()}, rec


if __name__ == "__main__":
    T0 = time.time()
    ra, reca, sh = part_a(); print("a d=24 (pre-registered)", json.dumps(ra, default=float)[:900], "%.0fs" % (time.time() - T0), flush=True)
    # POST-HOC (not pre-registered): d = 24 stays pinned by the lattice (Peierls-Nabarro barrier); the same analysis
    # at d = 8, where the attraction 2pi/d beats the pinning (scan d = 6, 8, 12, 16 done after seeing d = 24 fail:
    # 6, 8, 12 annihilate at t ~ 10, 15, 85; 16 still pinned at t = 145)
    ra8, reca8, sh8 = part_a(d=8); print("a d=8 (post-hoc)", json.dumps(ra8, default=float)[:900], flush=True)
    rb, trb = part_b(); print("b", rb, "%.0fs" % (time.time() - T0), flush=True)
    rc, recc = part_c(); print("c", rc, "%.0fs" % (time.time() - T0), flush=True)
    rd, recd = part_d(); print("d", rd, "%.0fs" % (time.time() - T0), flush=True)
    with open(os.path.join(OUT, "kz_k3.json"), "w") as f:
        json.dump({"a_pair": ra, "a_pair_d8_posthoc": ra8, "b_kick": rb, "c_KZ_2D": rc, "d_KZ_3D": rd}, f, indent=1, default=float)
    np.savez(os.path.join(OUT, "kz_k3_tracks.npz"), a_t=reca["t"], a_xp=reca["xp"], a_xm=reca["xm"],
             a_E=reca["E"], a_Ekin=reca["Ekin"], a_shells=sh, a8_t=reca8["t"], a8_xp=reca8["xp"], a8_xm=reca8["xm"], a8_E=reca8["E"], a8_Ekin=reca8["Ekin"],
             a8_shells=sh8, b=trb, c=recc,
             d=recd if recd is not None else np.zeros(0))

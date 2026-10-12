"""knot_k3.py - K3: knot + anti-knot (and knot + knot) head-on collisions from the product ansatz.
usage: python knot_k3.py kind v T     kind in {SA (B=1 + B=-1), SS (B=1 + B=1 attractive channel), S1 (single, control)}
Resumable in chunks of <= 330 s (state in KNOT_SCRATCH).
"""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
kind = sys.argv[1]; v = float(sys.argv[2]); T = float(sys.argv[3])
h = 0.2; N = 64; DT = torch.float32; dt = 0.25 * h; RS = 4.5; snaps_t = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
# dt: prereg said 0.1 h; 0.25 h is a declared deviation (cost), energy drift checked in the S1 control
tag = f"{kind}_v{v:g}"
state_f = os.path.join(SCR, f"knot_k3_state_{tag}.pt")
X, Y, Z = kc.grid((N, N, N), h, DT)
Rr = torch.sqrt(X**2 + Y**2 + Z**2)
inside = Rr < RS
sponge = kc.sponge_profile((N, N, N), h, 1.4, 2.0, DT)
t_wall = time.time()

if os.path.exists(state_f):
    S = torch.load(state_f); phi, pi, step, rec, snaps = S["phi"], S["pi"], S["step"], S["rec"], S["snaps"]
else:
    F = np.load(os.path.join(SCR, "knot_field_h0.2_L64.npz"))
    rb, fb = F["rb"], F["fb"]
    prof = dict(r=np.concatenate([[0.0], rb, [rb[-1] + h, 1e3]]), f=np.concatenate([[np.pi], fb, [0.0, 0.0]]))
    g = 1 / math.sqrt(1 - v * v)
    tp = (3.5, 4.3)
    if kind == "S1":
        q = kc.hedgehog(X, Y, Z, prof, taper=tp)
        phi = q; vel = torch.zeros_like(q)
    else:
        rot = np.diag([-1.0, -1.0, 1.0]) if kind == "SS" else None
        sgn = -1 if kind == "SA" else 1
        q1 = kc.hedgehog(X, Y, Z, prof, center=(-2.0, 0, 0), gx=g, taper=tp)
        q2 = kc.hedgehog(X, Y, Z, prof, center=(2.0, 0, 0), gx=g, sign=sgn, iso=rot, taper=tp)
        d1 = kc.derivs(q1, h)[0]; d2 = kc.derivs(q2, h)[0]
        phi = kc.normalize(kc.qmul(q1, q2))
        vel = kc.qmul(-v * d1, q2) + kc.qmul(q1, v * d2)       # knot 1 moves +x, knot 2 moves -x
    vel = kc.tangent(phi, vel)
    pi = kc.momentum_from_velocity(phi, vel, h)
    step, rec, snaps = 0, [], {}

nsteps = int(round(T / dt)); every = max(1, int(round(0.25 / dt)))


def measure():
    H, kin, st, _, _ = kc.ham_terms(phi, pi, h, need_grad=False)
    b = kc.baryon_density(phi, h); e = kin + st
    W = e.sum()
    out = dict(t=step * dt, E=H, E_in=float(e[inside].sum() * h**3), B=float(b.sum() * h**3),
               absB=float(b.abs().sum() * h**3), B_in=float(b[inside].sum() * h**3),
               Ixx=float((e * X**2).sum() / W), Iyy=float((e * Y**2).sum() / W), Izz=float((e * Z**2).sum() / W),
               bpos_x=float((b.abs() * X.abs()).sum() / b.abs().sum()), bpos_y=float((b.abs() * Y.abs()).sum() / b.abs().sum()),
               bpos_z=float((b.abs() * Z.abs()).sum() / b.abs().sum()))
    return out, b, e


while step <= nsteps:
    if step % every == 0 and (not rec or rec[-1]["t"] < step * dt - 1e-9):
        r, b, e = measure(); rec.append(r)
        print(f"t={r['t']:.2f} E={r['E']:.2f} E_in={r['E_in']:.2f} B={r['B']:.4f} |B|={r['absB']:.3f} "
              f"Ixx,Iyy,Izz={r['Ixx']:.2f},{r['Iyy']:.2f},{r['Izz']:.2f} pos={r['bpos_x']:.2f},{r['bpos_y']:.2f},{r['bpos_z']:.2f}", flush=True)
        for ts in snaps_t:
            if abs(r["t"] - ts) < 1e-6:
                snaps[str(ts)] = dict(b_xy=b[:, :, N // 2].numpy().tolist(), e_xy=e[:, :, N // 2].numpy().tolist(),
                                      b_xz=b[:, N // 2, :].numpy().tolist())
    if step == nsteps:
        break
    phi, pi = kc.rk4(phi, pi, h, dt, sponge=sponge)
    step += 1
    if time.time() - t_wall > 330:
        torch.save(dict(phi=phi, pi=pi, step=step, rec=rec, snaps=snaps), state_f)
        print("CHECKPOINT", step, flush=True); sys.exit(3)

E0in = rec[0]["E_in"]
res = dict(kind=kind, v=v, h=h, N=N, R_sphere=RS, T=T, records=rec, snaps=snaps,
           frac_out=(E0in - rec[-1]["E_in"]) / E0in, absB_final=rec[-1]["absB"], absB_min=min(r["absB"] for r in rec),
           t_absB_below_0p2=next((r["t"] for r in rec if r["absB"] < 0.2), None),
           B_final=rec[-1]["B"], E0=rec[0]["E"])
print(json.dumps({k: v for k, v in res.items() if k not in ("records", "snaps")}, indent=1))
json.dump(res, open(os.path.join(OUT, f"knot_k3_{tag}.json"), "w"))
if os.path.exists(state_f):
    os.remove(state_f)

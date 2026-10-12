"""knot2_boost.py - S2: boosted knot with a selectable stencil (knot_k2.py with knot2_core).
usage: python knot2_boost.py scheme h v T   (resumable, chunks <= 330 s, state in KNOT_SCRATCH)
scheme c4 uses the c4-relaxed field (knot2_relax.py); c2/old use the old relaxed field (same static energy).
"""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot2_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
scheme = sys.argv[1]; kc.set_scheme(scheme)
h = float(sys.argv[2]); v = float(sys.argv[3]); T = float(sys.argv[4])
DT = torch.float32
NY = {0.2: 48, 0.15: 64}[h]; NX = {0.2: 112, 0.15: 150}[h]
tag = f"{scheme}_h{h:g}_v{v:g}"
state_f = os.path.join(SCR, f"knot2_boost_state_{tag}.pt")
dt = 0.25 * h          # prereg said 0.1 h; 0.25 h declared deviation (RK4 drift 2e-5 per time unit at v = 0.8, dev test)
g = 1 / math.sqrt(1 - v * v)
X, Y, Z = kc.grid((NX, NY, NY), h, DT)
t_wall = time.time()

if os.path.exists(state_f):
    S = torch.load(state_f)
    phi, pi, step, rec, meta = S["phi"], S["pi"], S["step"], S["rec"], S["meta"]
else:
    F = np.load(os.path.join(SCR, f"knot2_field_c4_h{h:g}_L{NY}.npz" if scheme == "c4" else f"knot_field_h{h:g}_L{NY}.npz"))
    rb, fb = F["rb"], F["fb"]
    prof = dict(r=np.concatenate([[0.0], rb, [rb[-1] + h, 1e3]]), f=np.concatenate([[np.pi], fb, [0.0, 0.0]]))
    x0 = -(NX - 1) / 2 * h + 3.6
    rest = kc.hedgehog(X, Y, Z, prof, center=(x0, 0, 0))
    E0 = float(kc.energy(rest, h))
    phi = kc.hedgehog(X, Y, Z, prof, center=(x0, 0, 0), gx=g)
    vel = kc.tangent(phi, -v * kc.derivs(phi, h)[0])
    pi = kc.momentum_from_velocity(phi, vel, h)
    step, rec = 0, []
    meta = dict(scheme=scheme, h=h, v=v, gamma=g, E0=E0, x0=x0, box=[NX, NY, NY], dt=dt,
                r_half_contracted_over_h=float(np.interp(np.pi / 2, fb[::-1], rb[::-1])) / g / h)

nsteps = int(round(T / dt)); every = max(1, int(round(0.2 / dt)))


def measure():
    H, kin, st, _, _ = kc.ham_terms(phi, pi, h, need_grad=False)
    b = kc.baryon_density(phi, h)
    e = kin + st
    out = dict(t=step * dt, E=H, B=float(b.sum() * h**3))
    for name, w in (("b", b.clamp(min=0)), ("e", e)):
        W = w.sum(); cx = (w * X).sum() / W
        out[name + "_cx"] = float(cx)
        out[name + "_sx"] = float(torch.sqrt((w * (X - cx) ** 2).sum() / W))
        out[name + "_sy"] = float(torch.sqrt(0.5 * (w * (Y**2 + Z**2)).sum() / W))
    return out


while step <= nsteps:
    if step % every == 0 and (not rec or rec[-1]["t"] < step * dt - 1e-9):
        rec.append(measure())
        r = rec[-1]
        print(f"t={r['t']:.2f} E/gE0={r['E']/(g*meta['E0']):.5f} B={r['B']:.4f} cx={r['b_cx']:.3f} "
              f"sx/sy(b)={r['b_sx']/r['b_sy']:.4f} sx/sy(e)={r['e_sx']/r['e_sy']:.4f}", flush=True)
    if step == nsteps:
        break
    phi, pi = kc.rk4(phi, pi, h, dt)
    step += 1
    if time.time() - t_wall > 330:
        torch.save(dict(phi=phi, pi=pi, step=step, rec=rec, meta=meta), state_f)
        print("CHECKPOINT", step, flush=True)
        sys.exit(3)

# analysis
t = np.array([r["t"] for r in rec]); half = t >= t[-1] / 2
ratio_b = np.array([r["b_sx"] / r["b_sy"] for r in rec])
ratio_e = np.array([r["e_sx"] / r["e_sy"] for r in rec])
vfit = np.polyfit(t[half], np.array([r["b_cx"] for r in rec])[half], 1)[0]
Es = np.array([r["E"] for r in rec])
res = dict(meta=meta, records=rec,
           E_over_gE0=float(Es[0] / (g * meta["E0"])), drift=float((Es.max() - Es.min()) / Es[0]),
           ratio_b_late=float(ratio_b[half].mean()), ratio_b_t0=float(ratio_b[0]), one_over_gamma=1 / g,
           ratio_e_late=float(ratio_e[half].mean()), ratio_e_t0=float(ratio_e[0]),
           e_sx_end_over_t0=float(rec[-1]["e_sx"] / rec[0]["e_sx"]), v_fit=float(vfit),
           B_min=float(min(r["B"] for r in rec)), B_max=float(max(r["B"] for r in rec)))
print(json.dumps({k: v for k, v in res.items() if k != "records"}, indent=1))
json.dump(res, open(os.path.join(OUT, f"knot2_boost_{tag}.json"), "w"), indent=1)
if os.path.exists(state_f):
    os.remove(state_f)

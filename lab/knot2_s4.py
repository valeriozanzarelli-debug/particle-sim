"""knot2_s4.py - S4: knot + anti-knot head-on (scheme c4, h = 0.2, box 80^3, sponge), long run, energy radiated out of
r < 4.5 and the spectrum of the outgoing waves (omega vs k in the shell 4.5 <= r <= 6.3).
usage: python knot2_s4.py v T   (resumable, chunks <= 330 s, state in KNOT_SCRATCH)  -> lab/out/knot2_s4_v<v>.json
"""
import sys, os, json, math, time, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot2_core as kc

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
SCR = os.environ.get("KNOT_SCRATCH", OUT)
v = float(sys.argv[1]); T = float(sys.argv[2])
kc.set_scheme("c4")
h = 0.2; N = 80; DT = torch.float32; dt = 0.25 * h; RS = 4.5
spec_t = [12.0, 16.0, 20.0, 24.0, 28.0]; snaps_t = [0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 30.0]
tag = f"v{v:g}"
state_f = os.path.join(SCR, f"knot2_s4_state_{tag}.pt")
X, Y, Z = kc.grid((N, N, N), h, DT)
Rr = torch.sqrt(X**2 + Y**2 + Z**2)
inside = Rr < RS
sponge = kc.sponge_profile((N, N, N), h, 1.4, 2.0, DT)
r1, r2, r3, r4 = 4.5, 4.8, 6.0, 6.3
mask = torch.where(Rr < r2, torch.sin(0.5 * math.pi * ((Rr - r1) / (r2 - r1)).clamp(0, 1)) ** 2,
                   torch.cos(0.5 * math.pi * ((Rr - r3) / (r4 - r3)).clamp(0, 1)) ** 2)
mask = torch.where((Rr < r1) | (Rr > r4), torch.zeros_like(mask), mask).numpy()
rbins = np.arange(0, 7.9 + 0.2, 0.2)
ridx = np.digitize(Rr.numpy().ravel(), rbins)
t_wall = time.time()

if os.path.exists(state_f):
    S = torch.load(state_f); phi, pi, step, rec, snaps, specs = S["phi"], S["pi"], S["step"], S["rec"], S["snaps"], S["specs"]
else:
    F = np.load(os.path.join(SCR, "knot2_field_c4_h0.2_L64.npz"))
    rb, fb = F["rb"], F["fb"]
    prof = dict(r=np.concatenate([[0.0], rb, [rb[-1] + h, 1e3]]), f=np.concatenate([[np.pi], fb, [0.0, 0.0]]))
    g = 1 / math.sqrt(1 - v * v); tp = (3.5, 4.3)
    q1 = kc.hedgehog(X, Y, Z, prof, center=(-2.0, 0, 0), gx=g, taper=tp)
    q2 = kc.hedgehog(X, Y, Z, prof, center=(2.0, 0, 0), gx=g, sign=-1, taper=tp)
    d1 = kc.derivs(q1, h)[0]; d2 = kc.derivs(q2, h)[0]
    phi = kc.normalize(kc.qmul(q1, q2))
    vel = kc.tangent(phi, kc.qmul(-v * d1, q2) + kc.qmul(q1, v * d2))
    pi = kc.momentum_from_velocity(phi, vel, h)
    step, rec, snaps, specs = 0, [], {}, {}

nsteps = int(round(T / dt)); every = max(1, int(round(0.25 / dt)))
K = 2 * math.pi * np.fft.fftfreq(N, d=h)
KX, KY, KZ = np.meshgrid(K, K, K, indexing="ij")
Kmag = np.sqrt(KX**2 + KY**2 + KZ**2)
Wlat = np.sqrt(sum(kc.omega_1d(np.abs(A) * h, "c4") ** 2 for A in (KX, KY, KZ))) / h
kb = np.arange(0, Kmag.max() + 0.4, 2 * math.pi / (N * h)); kidx = np.digitize(Kmag.ravel(), kb)


def spectrum(phi, phid):
    P = np.zeros(len(kb) + 1); Pd = np.zeros(len(kb) + 1); WL = np.zeros(len(kb) + 1)
    for a in (1, 2, 3):
        fa = np.abs(np.fft.fftn(phi[a].numpy() * mask)) ** 2
        fd = np.abs(np.fft.fftn(phid[a].numpy() * mask)) ** 2
        P += np.bincount(kidx, fa.ravel(), len(kb) + 1); Pd += np.bincount(kidx, fd.ravel(), len(kb) + 1)
        WL += np.bincount(kidx, (fd * Wlat).ravel(), len(kb) + 1)
    return dict(P=P.tolist(), Pd=Pd.tolist(), Wlat_w=WL.tolist())


def measure():
    H, kin, st, phid, _ = kc.ham_terms(phi, pi, h, need_grad=False)
    b = kc.baryon_density(phi, h); e = kin + st
    out = dict(t=step * dt, E=H, E_in=float(e[inside].sum() * h**3), B=float(b.sum() * h**3),
               absB=float(b.abs().sum() * h**3),
               e_r=(np.bincount(ridx, e.numpy().ravel(), len(rbins) + 1) * h**3).tolist())
    return out, b, e, phid


while step <= nsteps:
    if step % every == 0 and (not rec or rec[-1]["t"] < step * dt - 1e-9):
        r, b, e, phid = measure(); rec.append(r)
        print(f"t={r['t']:.2f} E={r['E']:.2f} E_in={r['E_in']:.2f} B={r['B']:.4f} |B|={r['absB']:.3f}", flush=True)
        for ts in snaps_t:
            if abs(r["t"] - ts) < 1e-6:
                snaps[str(ts)] = dict(e_xy=e[:, :, N // 2].numpy().tolist(), pi1_xy=phi[1, :, :, N // 2].numpy().tolist())
        for ts in spec_t:
            if abs(r["t"] - ts) < 1e-6:
                specs[str(ts)] = spectrum(phi, phid)
    if step == nsteps:
        break
    phi, pi = kc.rk4(phi, pi, h, dt, sponge=sponge)
    step += 1
    if time.time() - t_wall > 330:
        torch.save(dict(phi=phi, pi=pi, step=step, rec=rec, snaps=snaps, specs=specs), state_f)
        print("CHECKPOINT", step, flush=True); sys.exit(3)

E0in = rec[0]["E_in"]
fr = [(E0in - r["E_in"]) / E0in for r in rec]
res = dict(v=v, h=h, N=N, R_sphere=RS, T=T, kbins=kb.tolist(), records=rec, snaps=snaps, specs=specs, rbins=rbins.tolist(),
           frac_out_T=fr[-1], frac_out_24=next((f for r, f in zip(rec, fr) if abs(r["t"] - 24) < 1e-6), None),
           t_cross_70=next((r["t"] for r, f in zip(rec, fr) if f >= 0.7), None),
           absB_min=min(r["absB"] for r in rec), t_absB_below_0p2=next((r["t"] for r in rec if r["absB"] < 0.2), None),
           B_final=rec[-1]["B"], E0=rec[0]["E"])
print(json.dumps({k: v for k, v in res.items() if k not in ("records", "snaps", "specs", "kbins", "rbins")}, indent=1))
json.dump(res, open(os.path.join(OUT, f"knot2_s4_{tag}.json"), "w"))
if os.path.exists(state_f):
    os.remove(state_f)

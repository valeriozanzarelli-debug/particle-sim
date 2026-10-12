"""fp_run.py POL [budget_s] - production run of FOTONE A PALLINE (pre-registered in out/fp_prereg.json).
Runs in chunks with checkpoint (out/fp_ck_POL.npz); exit 3 = run me again. Per-tick metrics -> out/fp_run_POL.json,
ball frames (for the drawing) -> out/fp_frames_POL/tNNN.npz (not in git)."""
import sys, os, json, time
import numpy as np
from fp_core import LightK, pulse, to_k, to_x, energy, poynting, seed, advance, sample

PRE = json.load(open(os.path.join(os.path.dirname(__file__), 'out', 'fp_prereg.json')))
C = PRE['configurazioni']
C.update(json.load(open(os.path.join(os.path.dirname(__file__), 'out', 'fp_amend.json')))['configurazioni_nuove'])
POL = sys.argv[1]
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 540
LAM = C['celle_per_lambda']
SH = tuple(C['scatola_celle'])
NT = C['tick_volo']
NB = C['palline']
OUT = os.path.join(os.path.dirname(__file__), 'out')
CK = os.path.join(OUT, f'fp_ck_{POL}.npz')
FR = os.path.join(OUT, f'fp_frames_{POL}')
os.makedirs(FR, exist_ok=True)
AX = (SH[0] / 2, SH[1] / 2)              # beam axis (site coordinates)
BIN = 4


def proj_hist(pos):
    x = np.mod(pos[0], SH[0]); z = np.mod(pos[2], SH[2])
    h, _, _ = np.histogram2d(x, z, bins=(SH[0] // BIN, SH[2] // BIN), range=((-0.5, SH[0] - 0.5), (-0.5, SH[2] - 0.5)))
    return h / h.sum()


def proj_u(u):
    p = u.sum(1)
    p = p.reshape(SH[0] // BIN, BIN, SH[2] // BIN, BIN).sum((1, 3))
    return p / p.sum()


def f1(pos, u):
    h, p = proj_hist(pos), proj_u(u)
    return float(np.corrcoef(h.ravel(), p.ravel())[0, 1]), float(0.5 * np.abs(h - p).sum())


def field_metrics(F):
    u = energy(F)
    S = poynting(F)
    x = np.arange(SH[0])[:, None, None] - AX[0]
    y = np.arange(SH[1])[None, :, None] - AX[1]
    Jz = (x * S[1] - y * S[0]).sum()
    ang = 2 * np.pi * np.arange(SH[2]) / SH[2]
    uz = u.sum((0, 1))
    zc = np.angle((uz * np.exp(1j * ang)).sum()) * SH[2] / (2 * np.pi)
    return u, float(Jz / u.sum()), float(zc)


def ball_metrics(pos, v):
    rx, ry = pos[0] - AX[0], pos[1] - AX[1]
    return float(np.mean(rx * v[1] - ry * v[0])), float(np.mean(v[2]))


def save_frame(t, pos, F):
    E, B = sample(F, pos)
    np.savez(os.path.join(FR, f't{t:03d}.npz'), pos=pos.astype(np.float32), E=E.astype(np.float16), B=B.astype(np.float16))


t0 = time.time()
lk = LightK(SH)
if os.path.exists(CK):
    d = np.load(CK)
    Fk, pos, t = d['Fk'], d['pos'], int(d['t'])
    lk.tick_count = t
    rec = json.load(open(os.path.join(OUT, f'fp_run_{POL}.json')))
    F = to_x(Fk)
else:
    F = pulse(SH, LAM, C['w0_lambda'] * LAM, C['sigma_lungo_lambda'] * LAM, POL)
    Fk = to_k(F)
    rng = np.random.default_rng(C['seme'] + {'circ+': 0, 'circ-': 1, 'lin_x': 2}[POL])
    t = 0
    u, Jf, zc = field_metrics(F)
    pos = seed(u, NB, rng)
    c0, tv0 = f1(pos, u)
    rec = {'pol': POL, 'lam': LAM, 'shape': SH, 't': [], 'spin_balls': [], 'spin_field': [], 'vz_balls': [],
           'zc_field': [zc], 'F1_t0': [c0, tv0], 'stationary_frac': None}
    save_frame(0, pos, F)
    print(f'{POL} t=0 F1 corr {c0:.4f} TV {tv0:.4f} spin_field {Jf * 2 * np.pi / LAM:+.4f}', flush=True)
while t < NT:
    Fk1 = lk.tick(Fk)
    F1 = to_x(Fk1)
    pos, v = advance(pos, F, F1)
    Fk, F = Fk1, F1
    t += 1
    u, Jf, zc = field_metrics(F)
    sb, vz = ball_metrics(pos, v)
    rec['t'].append(t); rec['spin_balls'].append(sb * 2 * np.pi / LAM); rec['spin_field'].append(Jf * 2 * np.pi / LAM)
    rec['vz_balls'].append(vz); rec['zc_field'].append(zc)
    save_frame(t, pos, F)
    if t % 8 == 0:
        print(f'{POL} t={t} spin balls {rec["spin_balls"][-1]:+.4f} field {rec["spin_field"][-1]:+.4f} vz {vz:.4f} ({time.time() - t0:.0f}s)', flush=True)
    if t == NT:
        rec['F1_end'] = list(f1(pos, u))
        break
    if time.time() - t0 > BUDGET:
        break
json.dump(rec, open(os.path.join(OUT, f'fp_run_{POL}.json'), 'w'))
if t < NT:
    np.savez(CK, Fk=Fk, pos=pos, t=t)
    sys.exit(3)
if os.path.exists(CK):
    os.remove(CK)
print(POL, 'done', rec['F1_end'])

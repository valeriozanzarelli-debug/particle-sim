"""mg_r3.py - exam R3 (lab/out/mg_prereg.json): DRAWING levels of detail. Balls of the photon (FOTONE_PALLINE frames,
circ+ and lin_x, ticks 24..96 step 24) are aggregated by blocks of 2^k cells (k = 1, 2, 3): one drawn ball per
occupied block, weight = number of balls (= energy share), position = mean, velocity = mean.
Checks: total drawn energy, centre of energy, spin from the balls (J_z omega / U as in FOTONE_PALLINE F2)."""
import json, os
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
LAM, AX = 24, (96.0, 96.0)


def load(pol, t):
    d = np.load(os.path.join(OUT, f'fp_frames_{pol}', f't{t:03d}.npz'))
    p0 = np.load(os.path.join(OUT, f'fp_frames_{pol}', f't{t - 1:03d}.npz'))['pos'].astype(np.float64)
    pos = d['pos'].astype(np.float64)
    return pos, pos - p0          # velocity over the last tick (what the balls actually did)


def spin(w, pos, v):
    rx, ry = pos[0] - AX[0], pos[1] - AX[1]
    return float(np.sum(w * (rx * v[1] - ry * v[0])) / np.sum(w) * 2 * np.pi / LAM)


def aggregate(pos, v, b):
    key = np.floor(pos / b).astype(np.int64)
    k = (key[0] * 1000003 + key[1]) * 1000003 + key[2]
    u, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
    P = np.stack([np.bincount(inv, pos[c]) / cnt for c in range(3)])
    V = np.stack([np.bincount(inv, v[c]) / cnt for c in range(3)])
    return cnt.astype(float), P, V


res = {}
ok = True
for pol in ('circ+', 'lin_x'):
    for t in (24, 48, 72, 96):
        pos, v = load(pol, t)
        w0 = np.ones(pos.shape[1])
        s0, c0 = spin(w0, pos, v), pos.mean(1)
        for k in (1, 2, 3):
            w, P, V = aggregate(pos, v, 2 ** k)
            dE = abs(w.sum() - w0.sum()) / w0.sum()
            dc = float(np.abs((w * P).sum(1) / w.sum() - c0).max())
            s = spin(w, P, V)
            ds = abs(s - s0) / abs(s0) if pol != 'lin_x' else abs(s - s0)
            passed = dE <= 1e-6 and dc <= 0.05 and ds <= 0.01
            ok &= passed
            res[f'{pol}_t{t}_k{k}'] = {'palline': int(len(w)), 'dE': dE, 'd_baricentro': dc, 'spin0': s0, 'spin_k': s, 'scarto_spin': ds, 'pass': bool(passed)}
            print(pol, t, k, len(w), f'dE {dE:.1e} dc {dc:.2e} spin {s0:+.4f} -> {s:+.4f} ({ds:.2e})', passed)
res['R3_pass'] = bool(ok)
json.dump(res, open(os.path.join(OUT, 'mg_r3.json'), 'w'), indent=1)
print('R3', 'PASS' if ok else 'FAIL')

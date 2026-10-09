"""kk_cyclotron.py - test 3a: force from geometry. A uniform twist of the hidden circle in 2+1 D
(hop along y also slides the hidden coordinate by R b (x - x0)) -> cyclotron circles of opposite sense
for n = +1 / -1, straight line for n = 0. Deterministic. Writes lab/out/kk_cyclotron.json/.npz"""
import json, time
import numpy as np
import kk_core as K
from zz_core import bloch_batch

t0 = time.time()
out = {}

# ------------------------- (a) cyclotron in a uniform twist -------------------------
nx = ny = 128; N = 24; th = 0.2; p = 0.3; b = 0.015; T = 220
R = N / (2 * np.pi)
x0 = nx // 2
X = np.arange(nx)[:, None] * np.ones((1, ny))
D = {1: R * b * (X - x0)}            # hop along y also slides the hidden coordinate by R b (x - x0)
# lattice energy and group velocity at k = (0, p, 2 pi n / N) for the prediction
def Ev(n, k):
    KA = np.array([[0, k, 2 * np.pi * n / N]])
    om = -np.angle(np.linalg.eigvals(bloch_batch(KA, th)[0])) / 2
    return np.max(om)
cyc = {}; snaps = {}
for n in (+1, -1, 0):
    psi = K.packet(nx, ny, N, th, n, k0=(0.0, p), sigma=8.0, center=(x0, ny // 2))
    Wk = K.KKWalk(nx, ny, N, th, D=D)
    xs = np.arange(nx); ys = np.arange(ny)
    tr = []; pk = []; acc = np.zeros((nx, ny))
    for t in range((T if n != 0 else 120) + 1):
        if t % 2 == 0:
            rho = K.density(psi)
            # circular centroid (periodic box), unwrapped below
            tr.append((x0 + np.angle(np.sum(rho * np.exp(2j * np.pi * (xs[:, None] - x0) / nx))) * nx / (2 * np.pi),
                       ny // 2 + np.angle(np.sum(rho * np.exp(2j * np.pi * (ys[None, :] - ny // 2) / ny))) * ny / (2 * np.pi)))
            i_, j_ = np.unravel_index(np.argmax(rho), rho.shape)
            pk.append((i_, j_))
            if t % 20 == 0:
                acc += rho
        psi = Wk.tick(psi)
    tr = np.array(tr)
    tr[:, 0] = np.unwrap(tr[:, 0] * 2 * np.pi / nx) * nx / (2 * np.pi)
    tr[:, 1] = np.unwrap(tr[:, 1] * 2 * np.pi / ny) * ny / (2 * np.pi)
    E = Ev(n, p); v = (Ev(n, p + 1e-4) - Ev(n, p - 1e-4)) / 2e-4
    pk = np.array(pk, float)
    if n != 0:
        A2 = np.c_[2 * pk[:, 0], 2 * pk[:, 1], np.ones(len(pk))]
        px, py, p0 = np.linalg.lstsq(A2, (pk**2).sum(1), rcond=None)[0]
        r_peak = float(np.sqrt(p0 + px**2 + py**2))
        angp = np.unwrap(np.arctan2(pk[:, 1] - py, pk[:, 0] - px))
        om_peak = float(np.polyfit(2 * np.arange(len(pk)), angp, 1)[0])
        # fit a circle to the centroid track (algebraic least squares)
        A = np.c_[2 * tr[:, 0], 2 * tr[:, 1], np.ones(len(tr))]
        cx, cy, c0 = np.linalg.lstsq(A, (tr**2).sum(1), rcond=None)[0]
        r_fit = float(np.sqrt(c0 + cx**2 + cy**2))
        ang = np.unwrap(np.arctan2(tr[:, 1] - cy, tr[:, 0] - cx))
        omega = float(np.polyfit(2 * np.arange(len(tr)), ang, 1)[0])      # rad per tick, sign = sense
        r_pred = p / (abs(n) * b)
        T_pred = 2 * np.pi * p / (v * abs(n) * b)
        cyc[n] = dict(peak_r_fit=r_peak, peak_center=[float(px), float(py)], peak_period=float(2 * np.pi / abs(om_peak)),
                      center=[float(cx), float(cy)], r_fit=r_fit, r_pred_p_over_nb=r_pred,
                      sense="counter-clockwise" if omega > 0 else "clockwise", omega=omega,
                      period_fit=float(2 * np.pi / abs(omega)), period_pred=float(T_pred), E=float(E), v=float(v))
        print(f"(a) n={n:+d}: density PEAK: circle radius {r_peak:6.2f} (pred {p/(abs(n)*b):.2f}), centre x = {px:5.1f}, "
              f"period {2*np.pi/abs(om_peak):6.1f} ticks (pred {T_pred:.1f}), sense {'counter-clockwise' if om_peak > 0 else 'clockwise'}")
        print(f"(a) n={n:+d}: CENTROID radius fit {r_fit:6.2f} cells (pred p/(|n|b) = {r_pred:.2f}), centre x = {cx:6.1f} "
              f"(start {x0}), sense {cyc[n]['sense']}, period {2*np.pi/abs(omega):6.1f} ticks (pred {T_pred:.1f})")
    else:
        dev = float(np.max(np.abs(tr[:, 0] - x0)))
        dist = float(tr[-1, 1] - tr[0, 1])
        cyc[n] = dict(max_sideways_deviation=dev, distance_travelled=dist, v=float(v))
        print(f"(a) n= 0: max sideways deviation {dev:.2e} cells over {dist:.1f} cells travelled (straight)")
    snaps[n] = (acc, tr if n == 0 else pk)
out["cyclotron"] = {str(k): v for k, v in cyc.items()}
out["cyclotron_params"] = dict(nx=nx, N=N, theta=th, p=p, b=b, ticks=T)

out["time_s"] = time.time() - t0
json.dump(out, open("out/kk_cyclotron.json", "w"), indent=1)
np.savez_compressed("out/kk_cyclotron_tmp.npz", acc_p=snaps[1][0], tr_p=snaps[1][1], acc_m=snaps[-1][0],
                    tr_m=snaps[-1][1], acc_0=snaps[0][0], tr_0=snaps[0][1])
print("time %.1f s" % out["time_s"])

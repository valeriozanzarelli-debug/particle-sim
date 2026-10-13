"""emiss_levels.py - preparation (not a test): the levels of the 3-D zig-zag atom (Nd59+, a = 0.5) by LISTENING.

Start state: a level of the radial form of the same rule put on the cubic lattice (zz.rs con_stato); then the 3-D
walk runs and the balls' amplitudes are summed with a Gaussian time window (std S_T) at the frequency w0 (only
even ticks: the walk's period is two ticks).  Levels further than ~3/S_T from w0 cancel out.  Two passes: the
second listens at the Rayleigh frequency of the first.  The level's frequency = -arg<phi|W2|phi>/(2a) and its
purity = |W2 phi - lambda phi|.

Usage: python lab/emiss_levels.py N STATE [S_T]    STATE in 1s+, 1s-, 2p1/2+, 2p3/2+, 2p3/2+3
Writes lab/out/emiss_ck/level_N_STATE.npy (big, not committed) and lab/out/emiss_levels_N.json.
"""
import sys, json, time, os
import numpy as np
from emiss_core import *

STATES = {"2p1/2-": (1, 2, -0.5), "2p3/2-": (-2, 2, -0.5), "1s+": (-1, 1, 0.5), "1s-": (-1, 1, -0.5), "2p1/2+": (1, 2, 0.5), "2p3/2+": (-2, 2, 0.5), "2p3/2+3": (-2, 2, 1.5)}


def fname(n, st):
    return os.path.join(CK, f"level_{n}_{st.replace('/', '_')}.npy")


def rayleigh(D, phi):
    q = phi.copy()
    D.tick(q, 0)
    D.tick(q, 1)
    lam = np.vdot(phi, q) / np.vdot(phi, phi)
    om = -np.angle(lam) / (2 * D.a)
    res = np.sqrt(norm2(q - lam * phi) / norm2(phi))
    return om, res


def listen(D, psi0, w0, S_T, t_budget=520, ckname=None):
    """Gaussian-window lock-in (even ticks), resumable: raises SystemExit(3) after t_budget s with a checkpoint."""
    a = D.a
    nt = int(round(7 * S_T / a))
    nt += nt % 2
    c0 = nt / 2
    t_start = 0
    if ckname and os.path.exists(ckname):
        z = np.load(ckname)
        acc, psi, t_start = z["acc"], z["psi"], int(z["t"])
        print(f"listening resumed at tick {t_start}/{nt}")
    else:
        acc = np.zeros_like(psi0)
        psi = psi0.copy()
    t0 = time.time()
    for t in range(t_start, nt + 1):
        if t % 2 == 0:
            w = np.exp(-0.5 * ((t - c0) * a / S_T)**2)
            acc += (w * np.exp(1j * w0 * t * a)) * psi
            if ckname and time.time() - t0 > t_budget and t < nt:
                np.savez(ckname, acc=acc, psi=psi, t=t + 1)
                print(f"listening checkpoint at tick {t + 1}/{nt}; run again")
                raise SystemExit(3)
        if t < nt:
            D.tick(psi, t)
    if ckname and os.path.exists(ckname):
        os.remove(ckname)
    return acc / np.sqrt(norm2(acc))


if __name__ == "__main__":
    n = int(sys.argv[1]); st = sys.argv[2]; S_T = float(sys.argv[3]) if len(sys.argv) > 3 else 150.0
    a = 0.5
    Za = 60 * ALPHA
    D = Dirac3D(n, a)
    kap, nr, mj = STATES[st]
    t0 = time.time()
    if os.path.exists(fname(n, st)):
        psi = np.load(fname(n, st))
        print("refining the saved state")
    else:
        ar = 0.05 if kap == -1 else 0.5
        om_r, rg, g, f = radial_state(kap, nr, ar, Za)
        psi = dirac_3d_from_radial(rg, g, f, kap, mj, D.x)
    w0, r0 = rayleigh(D, psi)
    print(f"start: Rayleigh omega {w0:.6f}, residual {r0:.2e}")
    # first pass at the strongest level near the start (1s: the 3-D 1s sits at ~0.8995, the start's Rayleigh
    # value mixes it with the lattice levels 0.886 and 0.840): listen at the guess, then at the refined value
    # pre-pass: spectrum of <psi0|psi(t)> over 1200 ticks; the strongest peak in the state's window
    lo, hi = (0.89, 0.91) if kap == -1 else (0.955, 0.99)
    import glob
    old_ck = glob.glob(fname(n, st) + ".listen_*.npz")
    if old_ck:
        guess = float(old_ck[0].split("listen_")[1][:-4])
        phi = listen(D, psi, guess, S_T, ckname=old_ck[0])
        lo = None
    C, p = [], psi.copy()
    for t in range(1200 if lo is not None else 0):
        if t % 2 == 0:
            C.append(ovl(psi, p))
        D.tick(p, t)
    if lo is not None:
        C = np.array(C)
        tt = np.arange(len(C)) * 2 * a
        h = np.sin(np.pi * (np.arange(len(C)) + 0.5) / len(C))**2
        om = np.linspace(lo, hi, 2001)
        S = np.abs(np.exp(1j * np.outer(om, tt)) @ (h * C))
        guess = om[np.argmax(S)] if os.environ.get("GUESS") is None else float(os.environ["GUESS"])
        print(f"pre-pass: strongest level in [{lo}, {hi}] at {guess:.5f}")
        phi = listen(D, psi, guess, S_T, ckname=fname(n, st) + f'.listen_{guess:.5f}.npz')
    w1, r1 = rayleigh(D, phi)
    print(f"pass 1 (S_T {S_T}): omega {w1:.6f}, residual {r1:.2e}  ({time.time() - t0:.0f} s)")
    np.save(fname(n, st), phi)
    jf = os.path.join(OUT, f"emiss_levels_{n}.json")
    js = json.load(open(jf)) if os.path.exists(jf) else {}
    js[st] = dict(omega=w1, residual=r1, S_T=S_T, start_omega=w0, start_residual=r0, seconds=time.time() - t0,
                  dirac=dirac_formula(nr, kap, Za))
    json.dump(js, open(jf, "w"), indent=1)

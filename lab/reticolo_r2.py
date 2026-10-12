"""reticolo_r2.py - R2b/R2c (+ R2a on each side): a pulse crossing a level interface, fine->coarse and
coarse->fine, normal and 30 deg incidence, 2-D (and 3-D with --d3).  Usage:
   python reticolo_r2.py N [N ...] [--ang 30] [--d3]
Writes lab/out/reticolo_r2_cases/<tag>.json."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reticolo_core import (Level, World, packet, helicity_fractions, restrict_dense, OUT)

CASES = os.path.join(OUT, "reticolo_r2_cases")
os.makedirs(CASES, exist_ok=True)


def blocks_range(x0, x1, B, ny, d):
    """block coords covering fine-cell range [x0,x1) in x (x0,x1 multiples of B), all ny blocks in y (and z)"""
    xs = np.arange(x0 // B, x1 // B)
    grids = np.meshgrid(xs, *[np.arange(ny)] * (d - 1), indexing="ij")
    return np.stack(grids, -1).reshape(-1, d)


def run(N, ang, direction, d=2):
    t0 = time.time()
    B = 16 if d == 2 else 8
    CB = 2 * B                            # coarse block in fine cells
    lam_f = 2 * N
    sig = 3.0 * lam_f
    S = int(np.ceil(4 * sig))
    if ang == 0:
        Ly = CB
    else:
        lam_y = lam_f / np.sin(np.radians(ang))
        Ly = CB
        while abs(Ly / lam_y - round(Ly / lam_y)) > 1e-9:
            Ly += CB
    Xi = int(np.ceil((2 * S + 64) / CB) * CB)
    cosa = np.cos(np.radians(ang))
    D = 2 * S + 96
    T = int(np.ceil(D / cosa / 2.0)) * 2           # fine ticks (even)
    Xtot = int(np.ceil((Xi + D + S + 96) / CB) * CB)
    shape_f = (Xtot,) + (Ly,) * (d - 1)
    shape_c = (Xtot // 2,) + (Ly // 2,) * (d - 1)
    per_f = [0] + [Ly // B] * (d - 1)
    per_c = [0] + [Ly // CB] * (d - 1)
    pa = tuple(range(1, d))
    x0 = Xi - S - 32
    if direction == "f2c":
        Fd = packet(shape_f, d, lam_f, ang, [x0] + [Ly / 2] * (d - 1), sig, periodic_axes=pa)
        fine = Level(0, d, B, blocks_range(0, Xi, B, Ly // B, d), per_f)
        coarse = Level(1, d, B, blocks_range(0, Xtot // 2, B, Ly // CB, d), per_c)
        fine.set_dense([0] * d, Fd)
        coarse.set_dense([0] * d, restrict_dense(Fd, d, pa))
        ref = Level(0, d, B, blocks_range(0, Xtot, B, Ly // B, d), per_f)
        ref.set_dense([0] * d, Fd)
        ref_ticks, ref_shape, ref_is_fine = T, shape_f, True
    else:
        Fc = packet(shape_c, d, N, ang, [x0 / 2] + [Ly / 4] * (d - 1), sig / 2, periodic_axes=pa)
        fine = Level(0, d, B, blocks_range(Xi, Xtot, B, Ly // B, d), per_f)
        coarse = Level(1, d, B, blocks_range(0, Xtot // 2, B, Ly // CB, d), per_c)
        coarse.set_dense([0] * d, Fc)
        ref = Level(1, d, B, blocks_range(0, Xtot // 2, B, Ly // CB, d), per_c)
        ref.set_dense([0] * d, Fc)
        ref_ticks, ref_shape, ref_is_fine = T // 2, shape_c, False
    w = World([fine, coarse])
    w.sync_all_up()
    if direction == "c2f":
        w.pairs[0].prolong_all()            # fine region initialised from coarse (zero there anyway)
    E0 = w.energy()
    t_init = time.time() - t0
    cent = []
    for n in range(T // 2):
        w.step_top(1)
        if n % 4 == 0 or n == T // 2 - 1:
            E, M = w.moments()
            cent.append((2 * (n + 1), M[0], E))
    for _ in range(ref_ticks):
        ref.tick()
    R = w.energy(lambda X: X[:, 0] < Xi - 16) / E0
    Tt = w.energy(lambda X: X[:, 0] > Xi + 16) / E0
    stuck = w.energy(lambda X: (X[:, 0] >= Xi - 16) & (X[:, 0] <= Xi + 16)) / E0
    Etot = w.energy() / E0
    # compare on the coarse grid in the transmitted region
    Fc_out = coarse.get_dense([0] * d, shape_c)
    Rd = ref.get_dense([0] * d, ref_shape)
    Rc = restrict_dense(Rd, d, pa) if ref_is_fine else Rd
    xc = (np.arange(shape_c[0]) + 0.5) * 2
    msk = (xc > Xi + 48).reshape((-1,) + (1,) * (d - 1))
    if direction == "c2f":
        msk = msk & (xc < Xtot - 48).reshape((-1,) + (1,) * (d - 1))
    A = Fc_out * msk[..., None]
    Bv = Rc * msk[..., None]
    ov = np.vdot(Bv, A)
    fid = abs(ov) / np.sqrt(np.vdot(A, A).real * np.vdot(Bv, Bv).real)
    phase = float(np.angle(ov))
    amp_ratio = float(np.sqrt(np.vdot(A, A).real / np.vdot(Bv, Bv).real))
    hel = helicity_fractions(A, d)
    # velocities on each side (centroid slope where the packet is entirely on one side)
    c = np.array(cent)
    side1 = c[np.abs(c[:, 1] - Xi) > S + 32]
    v_before = v_after = None
    pre = side1[side1[:, 1] < Xi]
    post = side1[side1[:, 1] > Xi]
    if len(pre) > 2:
        v_before = float(np.polyfit(pre[:, 0], pre[:, 1], 1)[0] / cosa)
    if len(post) > 2:
        v_after = float(np.polyfit(post[:, 0], post[:, 1], 1)[0] / cosa)
    res = dict(N=N, ang=ang, direction=direction, d=d, lam_fine=lam_f, sigma_fine=sig, Ly=Ly, Xi=Xi, Xtot=Xtot,
               fine_ticks=T, T=Tt, R=R, stuck=stuck, E_total=Etot, lost=1 - Tt - R - stuck,
               fidelity=float(fid), phase=phase, amp_ratio=amp_ratio,
               hel_plus=hel[0], hel_minus=hel[1], hel_static=hel[2],
               v_before=v_before, v_after=v_after, cells=w.cells(), t_init=t_init, wall=time.time() - t0)
    return res


if __name__ == "__main__":
    args = sys.argv[1:]
    ang = 0
    d = 2
    dirs = ["f2c", "c2f"]
    if "--ang" in args:
        i = args.index("--ang"); ang = float(args[i + 1]); del args[i:i + 2]
    if "--d3" in args:
        args.remove("--d3"); d = 3
    if "--dir" in args:
        i = args.index("--dir"); dirs = [args[i + 1]]; del args[i:i + 2]
    for N in [int(a) for a in args]:
        for dr in dirs:
            r = run(N, ang, dr, d)
            tag = "d%d_N%d_a%d_%s" % (d, N, ang, dr)
            json.dump(r, open(os.path.join(CASES, tag + ".json"), "w"), indent=1)
            print(tag, "T=%.6f R=%.2e stuck=%.2e lost=%.2e fid=%.6f ph=%.2e amp=%.5f hel-=%.1e hel0=%.1e v=%s/%s %.0fs"
                  % (r["T"], r["R"], r["stuck"], r["lost"], r["fidelity"], r["phase"], r["amp_ratio"], r["hel_minus"],
                     r["hel_static"], r["v_before"], r["v_after"], r["wall"]), flush=True)

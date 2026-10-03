"""Checks of the machinery of the base rule (not experiments):
  1. Rust field sums == numpy reference of the same rule (2-D with blocking, 3-D vector kernel)
  2. 2-D: the balls' field behind the double slit vs the exact Rayleigh–Sommerfeld integral, near and far
  3. 3-D: the vector field is divergence-free and E_x reduces to the scalar RS integral; a line of 3-D
     origins along y reproduces the 2-D Hankel kernel (3-D → 2-D reduction of the kernel)
  4. speed
Run: python lab/base_test.py
"""
import time
import numpy as np
from base_rule import *
from optics import rs_field, sources, inc_amp, SLITS, A, W0, L


def double_slit_scene(rule=RULE, zs=10.0, jones=(1, 0), filters=None):
    sc = Scene2D(rule)
    sc.add(Obj2("source", "source", (0, -zs), (1, 0), 40.0, src=dict(w0=W0, zw=zs, jones=jones)))
    ints = [(c - A / 2, c + A / 2) for c in SLITS]
    sc.add(Obj2("wall", "aperture", (0, 0), (1, 0), 400.0, ints=ints))
    return sc


def main():
    log("== 1. Rust vs numpy reference")
    sc = double_slit_scene()
    rng = np.random.default_rng(1)
    P = np.column_stack([rng.uniform(-6, 6, 40), np.concatenate([rng.uniform(0.02, 0.5, 20), rng.uniform(0.5, 30, 20)])])
    t0 = time.time(); U, Ux, Uz = sc.field(P); t1 = time.time()
    Ur, Uxr, Uzr = field2_numpy(sc, P)
    e = max(np.abs(U - Ur).max(), np.abs(Ux - Uxr).max() / K, np.abs(Uz - Uzr).max() / K) / np.abs(Ur).max()
    log(f"  2-D: max rel diff Rust/numpy {e:.2e}   (Rust {1e3*(t1-t0):.1f} ms for 40 points)")
    # points on the source side too (they hear the source fan, not the slit fan; rays through the wall blocked)
    P2 = np.column_stack([rng.uniform(-6, 6, 10), rng.uniform(-9, -0.05, 10)])
    a, _, _ = sc.field(P2); b, _, _ = field2_numpy(sc, P2)
    log(f"  2-D source side: max rel diff {np.abs(a-b).max()/np.abs(b).max():.2e}")

    s3 = Scene3D()
    s3.add(Obj3("src", "source", 0.0, mask=("circ_in", 0, 0, 1.5), src=dict(w0=np.inf, jones=(1, 0.3j))))
    P3 = np.column_stack([rng.uniform(-2, 2, 12), rng.uniform(-2, 2, 12), np.concatenate([rng.uniform(0.02, 0.3, 6), rng.uniform(1, 8, 6)])])
    E, dE = s3.field(P3, 0)
    Er, dEr = field3_numpy(s3, P3, 0)
    log(f"  3-D: max rel diff Rust/numpy E {np.abs(E-Er).max()/np.abs(Er).max():.2e}  dE {np.abs(dE-dEr).max()/np.abs(dEr).max():.2e}")
    div = np.abs(dE[:, 0, 0] + dE[:, 1, 1] + dE[:, 2, 2]).max() / (K * np.abs(E).max())
    log(f"  3-D: |div E| / (k|E|) max {div:.2e}  (transversality of the vector RS field)")

    log("== 2. 2-D double slit: balls' field vs exact RS (incident = Gaussian waist, exact)")
    s0, w0 = sources(2000)
    for z in (0.02, 0.1, 0.5, 2, 10, 160):
        x = np.linspace(-6, 6, 301) if z < 20 else np.linspace(-60, 60, 601)
        P = np.column_stack([x, np.full(x.size, z)])
        U, Ux, Uz = sc.field(P)
        Ur, Uxr, Uzr = [q[:, 0] for q in rs_field(x, z, s0, w0)]
        jz = np.imag(np.conj(U[:, 0]) * Uz[:, 0]); jr = np.imag(np.conj(Ur) * Uzr)
        log(f"  z={z:6}: max|U-U_RS|/max|U_RS| {np.abs(U[:,0]-Ur).max()/np.abs(Ur).max():.2e}   rms(j_z err)/max j_z {np.sqrt(np.mean((jz-jr)**2))/jr.max():.2e}")
    log("  ball count:", balls_summary(sc))

    log("== 3. 3-D -> 2-D reduction of the kernel: a y-line of 3-D origins vs the 2-D Hankel wavelet")
    # long strip |x|<0.5 (Gaussian taper in y, wide) emitting a plane wave; field at y=0 vs 2-D strip
    wy = 40.0
    s3b = Scene3D()
    s3b.add(Obj3("strip", "source", 0.0, mask=("rects", [(-0.5, 0.5, -200, 200)]), src=dict(w0=np.inf)))
    s2b = Scene2D()
    s2b.add(Obj2("strip", "source", (0, 0), (1, 0), 0.5, src=dict(w0=np.inf)))
    # NOTE: the 3-D strip uses a plane wave over |y| < 50 (box of the plane-wave source); compare near the centre
    zz = np.array([0.3, 1.0, 3.0])
    P = np.column_stack([np.zeros(3), np.zeros(3), zz])
    E, _ = s3b.field(P, 0)
    U2, _, _ = s2b.field(np.column_stack([np.zeros(3), zz]))
    for z, a, b in zip(zz, E[:, 0], U2[:, 0]):
        log(f"  z={z}: 3-D strip (|y|<50) E_x = {a:.4f}   2-D strip U = {b:.4f}   rel diff {abs(a-b)/abs(b):.2e}")

    log("== 4. speed")
    P = np.column_stack([rng.uniform(-30, 30, 2000), rng.uniform(5, 150, 2000)])
    t0 = time.time(); sc.field(P); dt = time.time() - t0
    n0 = sum(len(sc.origins(i)[0]) for i in range(len(sc.fans)))
    log(f"  2-D: 2000 points x {n0} level-0 origins: {dt*1e3:.0f} ms  ({2000*n0/dt/1e6:.0f} M pairs/s, {NTHREADS} threads)")
    P = np.column_stack([rng.uniform(-1, 1, 2000), rng.uniform(-1, 1, 2000), rng.uniform(5, 50, 2000)])
    n0 = len(s3.origins(0)[0])
    t0 = time.time(); s3.field(P, 0); dt = time.time() - t0
    log(f"  3-D: 2000 points x {n0} level-0 origins: {dt*1e3:.0f} ms  ({2000*n0/dt/1e6:.0f} M pairs/s)")
    P[:, 2] = 0.03
    t0 = time.time(); s3.field(P, 0); dt = time.time() - t0
    log(f"  3-D near field (z=0.03, all levels): {dt*1e3/2000:.3f} ms per point")


if __name__ == "__main__":
    main()

"""aggmol_costs.py - task (1) of study aggmol_: the COST LAW of a box of two atoms at every level (calculation + timing, not
an exam). Writes lab/out/aggmol_costs.json.

For each level (cell a0/16, a0/8, a0/4, a0/2) and box (+-8 a0: the two-atom ground state; +-12 a0: excited states and the
cloud of 5-10 a0 with margin):
  cells N^3, memory of the one-electron machine (eigen-solver float64, real-time float32), the native time step
  tau = cell / c (the coarse step of the rule), steps per atomic unit of time;
  MEASURED here (this container, numpy/scipy FFT, ONE thread): time of one application of the coarse Hamiltonian
  (kinetic + 8-lantern Galerkin field) and of one free step (kinetic only);
  ESTIMATED on hardware with a bandwidth model (catalogue numbers of docs/MOTORE_OTTIMIZZATO.md section 8, NOT measured:
  Intel UHD ~25 GB/s; 16 GB card 288 GB/s (mid) and 717 GB/s (high)): bytes moved per cell per step = 96 B per 3-D FFT of a
  2-component complex64 spinor (3 axis passes x read+write x 16 B) x number of FFTs + 16 B per pointwise pass;
  a free step = 2 FFT + 1 multiply; a step with the Galerkin field = 2 + 16 FFT + 9 multiplies;
  the two-electron exact problem (configuration space, M = N^3 one-electron cells): M^2 numbers (singlet: M(M+1)/2).
Usage: timeout 590 python lab/aggmol_costs.py"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_core as mc

JS = os.path.join(mc.OUT, "aggmol_costs.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
BW = {"Intel UHD (~25 GB/s, catalogue)": 25e9, "16 GB mid (288 GB/s)": 288e9, "16 GB high (717 GB/s)": 717e9}
B_FFT, B_MUL = 96.0, 16.0
AU_TIME = 1.0 / mc.ALPHA ** 2            # one atomic unit of time in hbar / m c^2


def measure(N, h, reps=3):
    T = mc.kinetic_exact(N, h)
    t0 = time.time()
    f = mc.MultiField(N, h, [(0, 0, 0.7 * mc.A0 / h), (0, 0, -0.7 * mc.A0 / h)], dtype=np.float32 if N > 128 else np.float64)
    tfield = time.time() - t0
    m = mc.Mol(T, f, sym=(0, 0, 0))
    v = np.random.default_rng(1).normal(size=(N,) * 3)
    m.H(v)
    t0 = time.time()
    for _ in range(reps):
        m.H(v)
    tH = (time.time() - t0) / reps
    psi = (v + 1j * v).astype(np.complex128)
    t0 = time.time()
    for _ in range(reps):
        mc.ifftn(np.exp(-1j * T) * mc.fftn(psi))
    tfree = (time.time() - t0) / reps
    return dict(t_field_build_s=tfield, t_H_apply_s=tH, t_free_step_s=tfree)


if __name__ == "__main__":
    t_start = time.time()
    rows = res.get("rows", [])
    done = {(r["frac"], r["box_a0"]) for r in rows}
    for frac in (2, 4, 8, 16):
        for box in (8, 12):
            if (frac, box) in done:
                continue
            h = mc.A0 / frac
            N = 2 * box * frac
            M = N ** 3
            row = dict(frac=frac, box_a0=box, N=N, cells=M, cell_lambdaC=h, tau_hbar_over_mc2=h,
                       steps_per_au=AU_TIME / h,
                       mem_eigen_GB=M * (8 * 12 + 8 * 8 + 8) / 1e9,       # ~12 float64 vectors (LOBPCG) + 8 lanterns + T
                       mem_dynamics_GB=M * (16 + 32 + 4 + 32) / 1e9)      # spinor c64 + 8 lanterns f32 + T f32 + 2 work spinors
            byt_free = M * (2 * B_FFT + B_MUL)
            byt_gal = M * (18 * B_FFT + 9 * B_MUL)
            row["bytes_per_step_free_GB"] = byt_free / 1e9
            row["bytes_per_step_field_GB"] = byt_gal / 1e9
            row["est_ms_per_step_field"] = {k: 1e3 * byt_gal / bw for k, bw in BW.items()}
            row["two_electron_exact_numbers"] = float(M) ** 2
            row["two_electron_singlet_GB_c64"] = M * (M + 1) / 2 * 8 / 1e9
            if N <= 192:
                row.update(measure(N, h))
            rows.append(row)
            res["rows"] = rows
            json.dump(res, open(JS, "w"), indent=1, default=float)
            print(f"a0/{frac} box +-{box}: N={N} cells {M:.3g} mem eigen {row['mem_eigen_GB']:.3f} GB dyn {row['mem_dynamics_GB']:.3f} GB; "
                  f"H apply {row.get('t_H_apply_s', float('nan')):.3f} s, free step {row.get('t_free_step_s', float('nan')):.4f} s "
                  f"(CPU 1 thread); Intel est {row['est_ms_per_step_field']['Intel UHD (~25 GB/s, catalogue)']:.1f} ms/step; "
                  f"2e exact {row['two_electron_singlet_GB_c64']:.3g} GB", flush=True)
            if time.time() - t_start > 480:
                print("checkpoint: run again"); sys.exit(0)
    print("COSTS DONE")

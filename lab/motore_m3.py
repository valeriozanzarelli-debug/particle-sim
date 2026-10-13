"""motore_m3.py - transitions between levels: D1 (light, N -> N/2 -> N: spectral vs cubic vs block mean),
D2 (electron, fine -> learned level B -> fine: agglomerator block + spectral prolongation = deconvolution of the
block, vs piecewise-constant). Pre-registration: lab/out/motore_prereg.json.
Usage: python lab/motore_m3.py  -> lab/out/motore_m3.json"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_core as mc
import reticolo_core as rc
from motore_m1 import resample

WP = np.array([-0.0546875, 0.8203125, 0.2734375, -0.0390625])   # cubic Lagrange at +1/4 on nodes I-1..I+2


def cubic_restrict(F):
    return rc.restrict_dense(F, 2, periodic_axes=(0, 1))


def cubic_prolong(Fc):
    out = Fc
    for ax in (0, 1):
        r = lambda s: np.roll(out, -s, axis=ax)
        plus = WP[0] * r(-1) + WP[1] * r(0) + WP[2] * r(1) + WP[3] * r(2)           # fine 2I+1 (offset +1/4)
        minus = WP[3] * r(-2) + WP[2] * r(-1) + WP[1] * r(0) + WP[0] * r(1)         # fine 2I   (offset -1/4)
        sh = list(out.shape); sh[ax] *= 2
        new = np.empty(sh, complex)
        idx = [slice(None)] * out.ndim
        idx[ax] = slice(0, None, 2); new[tuple(idx)] = minus
        idx[ax] = slice(1, None, 2); new[tuple(idx)] = plus
        out = new
    return out


def block_restrict(F):
    n = F.shape[0]
    return F.reshape(n // 2, 2, n // 2, 2, 3).mean(axis=(1, 3))


def block_prolong(Fc):
    return np.repeat(np.repeat(Fc, 2, axis=0), 2, axis=1)


res = dict(D1=[], D2=[])
t0 = time.time()
for Nc in (32, 16, 12, 8):
    Nf = 2 * Nc
    nf = int(round(27 * Nf)); nf += nf % 4
    Ff = rc.packet((nf, nf), 2, Nf, 30.0, [nf / 2, nf / 2], 3 * Nf)          # on the fine lattice + band
    Ef = mc.energy(Ff)
    cf = mc.centroid(Ff) / Nf
    wrong_f = 1 - rc.helicity_fractions(Ff, 2)[0]
    for meth, R, P in (("spettrale", lambda F: resample(F, F.shape[0] // 2), lambda G: resample(G, 2 * G.shape[0])),
                       ("cubico", cubic_restrict, cubic_prolong),
                       ("media_blocco", block_restrict, block_prolong)):
        Fc = R(Ff)
        Ec = 4 * mc.energy(Fc)
        cc = mc.centroid(Fc) * 2 / Nf
        hf = rc.helicity_fractions(Fc, 2)
        Fb = P(Fc)
        o = mc.overlap(Ff, Fb)
        row = dict(N_coarse=Nc, method=meth, dE_restrict=Ec / Ef - 1, centroid_jump_lambda=float(np.linalg.norm(cc - cf)),
                   wrong_band_coarse=1 - hf[0], wrong_band_fine_before=wrong_f, dE_roundtrip=mc.energy(Fb) / Ef - 1,
                   fidelity_roundtrip=float(abs(o) ** 2), phase_roundtrip=float(np.angle(o)))
        res["D1"].append(row)
        print("D1 Nc=%2d %-12s dE=%+.2e jump=%.2e wrong=%.2e  roundtrip dE=%+.2e fid=%.8f" % (
            Nc, meth, row["dE_restrict"], row["centroid_jump_lambda"], row["wrong_band_coarse"], row["dE_roundtrip"],
            row["fidelity_roundtrip"]), flush=True)
sp = [r for r in res["D1"] if r["method"] == "spettrale" and r["N_coarse"] >= 12]
res["D1_verdict"] = "PASS" if all(abs(r["dE_restrict"]) <= 1e-6 and r["centroid_jump_lambda"] <= 1e-3 and r["wrong_band_coarse"] <= 1e-3 for r in sp) else "FAIL"
print("D1", res["D1_verdict"], flush=True)

# ---------------- D2 electron ----------------
th, lv = mc.load_levels_1d()
nf = 196608
xf = (np.arange(nf) - nf / 2 + 0.5) * th
k = 2 * np.pi * np.fft.fftfreq(nf)
Efine = np.arccos(np.cos(th) * np.cos(k)) / th - 1.0         # positive-band kinetic energy per unit time (mc^2 = 1)
for k0a in (1.0, 2.0):
    psi = mc.positive_spinor_packet(xf, k0a / mc.A0, mc.A0 / 2, th)
    N0 = np.sum(np.abs(psi) ** 2)
    f0 = np.fft.fft(psi, axis=-1)
    Ek0 = float(np.sum(np.abs(f0) ** 2 * Efine) / np.sum(np.abs(f0) ** 2))
    for B in (16, 192, 384):
        pc = mc.block1(psi, B)
        nc = nf // B
        Ncoarse = np.sum(np.abs(pc) ** 2)
        # spectral prolongation: invert the block in Fourier (Dirichlet kernel), zero above the coarse Nyquist
        kc = 2 * np.pi * np.fft.fftfreq(nc)                   # coarse wave numbers (1/coarse cell)
        kf_of = kc / B
        S = np.sum(np.exp(1j * np.outer(np.arange(B), kf_of)), axis=0) / B
        # c_I = sqrt(B)/B sum_j psi(IB+j)  ->  c_hat(K) = sqrt(B) S(k) psi_hat(k) / B, k = K/B (band-limited content)
        Fc = np.fft.fft(pc, axis=-1)
        Ff = np.zeros((2, nf), complex)
        mf = np.fft.fftfreq(nc) * nc
        idx = (mf.astype(int)) % nf
        Ff[:, idx] = Fc * np.sqrt(B) / S
        back_sp = np.fft.ifft(Ff, axis=-1)
        back_pc = np.repeat(pc, B, axis=-1) / np.sqrt(B)
        Q_dressed = float(np.sum(np.abs(Fc) ** 2 / np.abs(S) ** 2) / nc)
        out = {}
        for name, bk in (("spettrale", back_sp), ("gradini", back_pc)):
            fb = np.fft.fft(bk, axis=-1)
            Ekb = float(np.sum(np.abs(fb) ** 2 * Efine) / np.sum(np.abs(fb) ** 2))
            out[name] = dict(norm_roundtrip=float(np.sum(np.abs(bk) ** 2) / N0 - 1), Ekin_roundtrip=Ekb / Ek0 - 1,
                             max_err=float(np.abs(bk - psi).max() / np.abs(psi).max()))
        row = dict(k0_a0=k0a, B=B, descriptor_norm_jump=float(Ncoarse / N0 - 1), dressed_charge_jump=Q_dressed / N0 - 1,
                   spectral=out["spettrale"], steps=out["gradini"])
        res["D2"].append(row)
        print("D2 k0=%g B=%3d desc norm %+.2e dressed %+.2e | spectral norm %+.2e Ek %+.2e err %.1e | steps norm %+.2e Ek %+.2e" % (
            k0a, B, row["descriptor_norm_jump"], row["dressed_charge_jump"], out["spettrale"]["norm_roundtrip"],
            out["spettrale"]["Ekin_roundtrip"], out["spettrale"]["max_err"], out["gradini"]["norm_roundtrip"],
            out["gradini"]["Ekin_roundtrip"]), flush=True)
res["D2_verdict"] = "PASS" if all(abs(r["spectral"]["norm_roundtrip"]) <= 1e-10 for r in res["D2"]) else "FAIL"
print("D2", res["D2_verdict"], "wall", time.time() - t0)
json.dump(res, open(os.path.join(mc.OUT, "motore_m3.json"), "w"), indent=1)

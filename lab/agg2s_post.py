"""agg2s_post.py - extra observables of the converged agg2s_ states (term-by-term budget of the 2s), from the saved vectors:
  band_err_eV  = <T_learned - T_D>      (T_D = sqrt(m^2 + k^2) - m, m = the learned kinetic mass: the band-shape error of
                                        the learned kinetic beyond mass and the Dirac form)
  krem_D_eV    = <T_D - k^2/2m>         (relativistic kinetic part with the exact Dirac band, on the coarse state)
  tail_w       = weight of |psi(k)|^2 with |k_i| h > 0.8 pi for some i (outer shell of the Brillouin cube), and the same weight
                 predicted by the Born tail 64 pi^2 rho_B / k^8 (s states): the solver's high-k tail diagnostic.
Also recomputes <k^4>_c and the kinetic remainder on the TWISTED Fourier grid (u = psi e^{-i theta x/L}); the first
version of agg2s_run.py took the FFT of psi itself for twisted jobs (wrong grid: ~2 % on <k^4> of the 1s); these values
replace the run-time ones in the analysis.
Usage: timeout 590 python lab/agg2s_post.py   (resumes; lab/out/agg2s_post.json)"""
import os, sys, json, time
import numpy as np
import agg2s_core as c
import agg2s_run as r

JS = os.path.join(c.OUT, "agg2s_post.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
at_js = json.load(open(os.path.join(c.OUT, "agg2s_atom.json")))
born = c.born_constants()
t0 = time.time()
for job, info in at_js.items():
    if not info.get("converged") or job in res or time.time() - t0 > 450:
        continue
    kin, B, N, kind, tws, dep, warm = r.JOBS[job]
    tw = tuple(int(x) for x in tws)
    h = B * c.TH
    if kin == "S":
        T = c.kin_schroedinger(N, h, tw)
    else:
        T, _ = c.kin_learned(r.level(kin, B), N, h, tw)
    m = info["m_kin"]
    v = np.load(r.vecfile(job)).reshape((N,) * 3)
    n = np.arange(N) - N // 2
    for a in range(3):                               # u = psi e^{-i theta x/L}: the twisted Fourier grid
        if tw[a]:
            sh = [1, 1, 1]
            sh[a] = N
            v = v * np.exp(-1j * np.pi * n / N).reshape(sh)
    F = c.sfft.fftn(v, workers=1)
    del v
    w = np.empty(F.shape)
    np.multiply(F.real, F.real, out=w)
    w += F.imag ** 2
    del F
    ka = c.kgrid_axes(N, h, tw)                      # 1/lambda_C
    kyz2 = ka[1][:, None] ** 2 + ka[2][None, :] ** 2
    kmax = np.pi / h
    oy = np.abs(ka[1])[:, None] > 0.8 * kmax
    oz = np.abs(ka[2])[None, :] > 0.8 * kmax
    sw = sbe = skd = stail = sborn = sk4 = srem = 0.0
    rho_B = info["rho_c_a0"] / (1 - 8 * np.pi * born["I_coef"] * m * info["h_a0"]) ** 2
    rho_l = rho_B * (c.ALPHA ** 3)                   # lambda_C^-3
    vol = (N * h) ** 3
    for i in range(N):
        k2 = ka[0][i] ** 2 + kyz2
        TD = np.sqrt(m * m + k2) - m
        wi = w[i]
        sw += wi.sum()
        sbe += np.sum(wi * (T[i] - TD))
        sk4 += np.sum(wi * k2 ** 2)
        srem += np.sum(wi * (T[i] - k2 / (2 * m)))
        skd += np.sum(wi * (TD - k2 / (2 * m)))
        outer = (abs(ka[0][i]) > 0.8 * kmax) | oy | oz
        stail += np.sum(wi[outer])
        with np.errstate(divide="ignore"):
            pb = 64 * np.pi ** 2 * rho_l * (m * c.ALPHA) ** 2 / np.where(k2 > 0, k2, np.inf) ** 4   # |psi(k)|^2 (continuum norm)
        sborn += np.sum(pb[outer]) / vol
    del w
    res[job] = dict(k4_c_a0=float(sk4 / sw) * c.A0 ** 4, kin_remainder_eV=float(srem / sw) * c.MC2 if kin != "S" else 0.0,
                    band_err_eV=float(sbe / sw) * c.MC2, krem_D_eV=float(skd / sw) * c.MC2, tail_w=float(stail / sw),
                    tail_born=float(sborn) if kind != "2p" else None)
    json.dump(res, open(JS, "w"), indent=1)
    print(job, res[job], flush=True)

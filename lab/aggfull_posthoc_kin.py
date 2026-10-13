"""aggfull_posthoc_kin.py - POST HOC (written after the F3 numbers, not pre-registered): where does the large kinetic
remainder <T - k^2/2m_kin> of the s states come from?  For the coarse 1s/2s states (vectors of aggfull_run.py) sample
k ~ |psi(k)|^2 and compare, at the same k:
   learned:  T_learned(k) - k^2/(2 m_kin)             (the agglomerated generator of chain A)
   fine:     eps_fine(k)  - k^2/(2 m_kin,fine)         (the fine zig-zag walk's own positive band, Bloch, exact)
   Dirac:    sqrt(m^2 + k^2) - m - k^2/2m              (continuum; -k^4/8 at small k)
Writes lab/out/aggfull_posthoc_kin.json"""
import json, os
import numpy as np
import aggfull_core as ac
import aggfull_atom as aa
import zz_core as zz

OUT = ac.OUT
rng = np.random.default_rng(7)
res = {}
for job, B, N in (("e8_1s", 512, 256), ("e8_2s", 512, 256), ("e16_1s", 256, 256)):
    fn = os.path.join(OUT, f"aggfull_vec_{job}.npy")
    if not os.path.exists(fn):
        continue
    h = B * ac.TH_A
    v = np.load(fn).reshape((N,) * 3)
    w = np.abs(ac.fftn(v)) ** 2
    w = w.ravel() / w.sum()
    idx = rng.choice(w.size, size=60000, p=w)
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
    ii = np.unravel_index(idx, (N,) * 3)
    K = np.stack([k1[ii[0]], k1[ii[1]], k1[ii[2]]], 1)
    k2 = (K ** 2).sum(1)
    rec = [l for l in ac.load_chain_A() if l["B"] == B][0]
    mk, _ = aa.kinetic_mass(rec, h)
    r = ac.genrule(rec)
    g0 = np.real(0.5 * np.trace(r.gen_at(np.zeros((1, 3)))[0]))
    gp = np.real(0.5 * np.trace(r.gen_at(K * h), axis1=1, axis2=2))
    gm = np.real(0.5 * np.trace(r.gen_at(-K * h), axis1=1, axis2=2))
    T_l = (0.5 * (gp + gm) - g0) / h
    # fine walk positive band energy (two-tick Bloch, eigenphase / 2 ticks), relative to k = 0
    th = ac.TH_A
    W = zz.bloch_batch(K * th, th)
    lam = np.linalg.eigvals(W)
    om = np.sort(-np.angle(lam) / (2 * th), axis=1)[:, -1]
    W0 = zz.bloch_batch(np.zeros((1, 3)), th)
    om0 = np.sort(-np.angle(np.linalg.eigvals(W0)) / (2 * th), axis=1)[:, -1][0]
    eps_f = om - om0
    mf = 1 + th ** 2 / 3
    R_l = T_l - k2 / (2 * mk)
    # amendment (aggfull_amend.json): the re-learned kinetic (aggfull_relearn.pkl), same samples
    RK = None
    if os.path.exists(os.path.join(OUT, "aggfull_relearn.pkl")):
        import pickle
        recK = [l for l in pickle.load(open(os.path.join(OUT, "aggfull_relearn.pkl"), "rb")) if l["B"] == B][0]
        mkK, _ = aa.kinetic_mass(recK, h)
        rK = ac.genrule(recK)
        g0K = np.real(0.5 * np.trace(rK.gen_at(np.zeros((1, 3)))[0]))
        TK = (0.5 * (np.real(0.5 * np.trace(rK.gen_at(K * h), axis1=1, axis2=2)) +
                     np.real(0.5 * np.trace(rK.gen_at(-K * h), axis1=1, axis2=2))) - g0K) / h
        RK = TK - k2 / (2 * mkK)
    R_f = eps_f - k2 / (2 * mf)
    R_d = np.sqrt(1 + k2) - 1 - k2 / 2
    big = k2 > np.percentile(k2, 99)
    res[job] = dict(m_kin_learned=mk, R_learned_eV=float(R_l.mean() * ac.MC2), R_fine_eV=float(R_f.mean() * ac.MC2),
                    R_dirac_eV=float(R_d.mean() * ac.MC2), R_learned_top1pct_share=float(R_l[big].sum() / R_l.sum()),
                    learned_minus_fine_rel_at_top1pct=float(np.mean((T_l[big] - eps_f[big]) / eps_f[big])))
    if RK is not None:
        res[job].update(m_kin_relearned=mkK, R_relearned_eV=float(RK.mean() * ac.MC2),
                        relearned_minus_fine_rel_at_top1pct=float(np.mean((TK[big] - eps_f[big]) / eps_f[big])))
    print(job, {k: (f"{v:.4e}" if isinstance(v, float) else v) for k, v in res[job].items()}, flush=True)
json.dump(res, open(os.path.join(OUT, "aggfull_posthoc_kin.json"), "w"), indent=1)

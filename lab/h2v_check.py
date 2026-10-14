"""h2v_check.py - VIDEO H2: exams read back from the frames and from the rendered images (V0 summary, V5, V6, V7).
    python3 lab/h2v_check.py            -> adds to lab/out/h2v_esami.json
"""
import os, sys, json
import numpy as np
from scipy import ndimage
from scipy.optimize import curve_fit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h2v_core as C
import h2v_sim as S
import h2v_render as RN

OUT = C.OUT
JE = os.path.join(OUT, "h2v_esami.json")
es = json.load(open(JE))


def side_cam(w, h):
    return RN.camera(np.array([0.0, RN.CAM_D, 0.0]), (0, 0, 0), RN.FOV, w, h)


def proton_sep(f, w=960, h=540):
    """render only the proton glows from frame f with an exact side camera; distance of the two brightness centroids"""
    RN.SOLO_P = True
    cam = side_cam(w, h)
    out, meta, _, _ = RN.scena(f, w, h, cam)
    RN.SOLO_P = False
    lum = out.sum(0)
    lab, n = ndimage.label(lum > 0.5 * lum.max())
    cs = ndimage.center_of_mass(lum, lab, range(1, n + 1))
    m = np.array(ndimage.sum(lum, lab, range(1, n + 1)))
    i = np.argsort(m)[-2:]
    (y1, x1), (y2, x2) = cs[i[0]], cs[i[1]]
    fpx = 0.5 * h / np.tan(np.radians(RN.FOV) / 2)
    return float(np.hypot(x1 - x2, y1 - y2) / fpx * RN.CAM_D), meta


def main():
    wf = json.load(open(os.path.join(OUT, "h2v_wf.json")))["E"]
    d = {k: v["diff"] for k, v in wf.items()}
    es["V0"] = dict(diff_per_R=d, max_abs=float(max(abs(x) for x in d.values())), n_R=len(d),
                    PASS=bool(len(d) == 11 and max(abs(x) for x in d.values()) <= 1e-6))
    # V5
    Rimg, meta = proton_sep(0)
    es["V5"] = dict(R_image=Rimg, R_true=1.4, rel=Rimg / 1.4 - 1, PASS=bool(abs(Rimg / 1.4 - 1) <= 0.02))
    print("V0", es["V0"]["max_abs"], "V5", es["V5"], flush=True)
    # V6: proton separation measured on rendered images of part C (every 2nd frame)
    fs = list(range(S.NA + S.NB, S.NF, 2))
    tt, RR, Rmeta = [], [], []
    for f in fs:
        r, m = proton_sep(f, 640, 360)
        tt.append((f - S.NA - S.NB) / S.FPS); RR.append(r); Rmeta.append(m["R"])
    tt, RR = np.array(tt), np.array(RR)
    fit = lambda t, a, b, w, p: a + b * np.cos(w * t + p)
    pp, _ = curve_fit(fit, tt, RR, p0=[RR.mean(), (RR.max() - RR.min()) / 2, 2 * np.pi / S.SLOW_PERIOD_S, 0])
    T_s = 2 * np.pi / abs(pp[2])
    V = C.vib()
    Tau = 2 * np.pi / (V["E"][1] - V["E"][0])
    T_fs = T_s / S.SLOW_PERIOD_S * Tau * C.AU_FS
    nu01 = json.load(open(os.path.join(OUT, "aggmol_h2.json")))["vib"]["nu01_cm"]
    T_ref = 1e15 / (2.99792458e10 * nu01)
    es["V6"] = dict(T_fs_frames=T_fs, T_fs_from_nu01=T_ref, rel=T_fs / T_ref - 1, R_mean=float(pp[0]), R_amp=float(abs(pp[1])),
                    slowdown=S.SLOW_PERIOD_S / (T_ref * 1e-15), PASS=bool(abs(T_fs / T_ref - 1) <= 0.01),
                    note="consistenza della catena (livelli -> fotogrammi -> immagine), non fisica nuova")
    np.save(os.path.join(S.SCR, "h2v_V6.npy"), np.stack([tt, RR, np.array(Rmeta)]))
    print("V6", es["V6"], flush=True)
    # V7: <|z|> of the electron balls at the extremes of <R>(t) vs analytic mixture
    frames = list(range(S.NA + S.NB, S.NF))
    Rm = np.array([json.loads(str(np.load(os.path.join(S.FR, f"f_{f:04d}.npz"))["meta"]))["R"] for f in frames])
    absz = {R: C.WF(R).abs_z() for R in C.RGRID}
    res = {}
    for name, sel in (("max", np.argsort(Rm)[-16:]), ("min", np.argsort(Rm)[:16])):
        bz, az = [], []
        for i in sel:
            z = np.load(os.path.join(S.FR, f"f_{frames[i]:04d}.npz"))
            m = json.loads(str(z["meta"]))
            P = z["p"][z["c"] == 0]
            bz.append(np.mean(np.abs(P[:, 2])))
            az.append(sum(w * absz[R] for w, R in zip(m["w"], C.RGRID)))
        res[name] = dict(R=float(Rm[sel].mean()), balls=float(np.mean(bz)), analytic=float(np.mean(az)))
    ok = all(abs(v["balls"] - v["analytic"]) <= 0.03 for v in res.values())
    ok &= np.sign(res["max"]["balls"] - res["min"]["balls"]) == np.sign(res["max"]["analytic"] - res["min"]["analytic"])
    es["V7"] = dict(res, abs_z_per_R=absz, PASS=bool(ok))
    print("V7", es["V7"], flush=True)
    json.dump(es, open(JE, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()

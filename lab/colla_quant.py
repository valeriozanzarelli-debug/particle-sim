"""colla_quant.py - C2a: rigid-rotor spin/isospin quantisation of the AN (m_pi = 138 MeV) Skyrmions B = 1..4.
Inertia tensors U (iso), V (spin), W (mixed) from the relaxed lattice fields (colla_ckpt_an_B*.pt), then
  nucleon  : M_N = E_1 + 3 hbar^2/(8 Lambda),   Lambda = mean of diag(U), diag(V), -diag(W)
  deuteron : J=1, I=0 (Braaten-Carson ground state)  -> E_rot = hbar^2 / V_perp  (needs W_perp ~ 0, checked)
  3He/3H   : J=I=1/2, K+L=0 (pre-registered) -> E_rot = (3/4)(u+v+2w)/(2(uv-w^2)) hbar^2, w = sign(det W)|W|
  4He      : J=I=0 -> E_rot = 0
Writes lab/out/colla_quant.json."""
import json, os, math, numpy as np, torch
import colla_core as cc
torch.set_num_threads(1)
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
U_ = json.load(open(os.path.join(OUT, "colla_units.json")))
AN = U_["AN"]; Eu, hb = AN["E_unit"], AN["hbar"]
H = 0.2
res = dict(units=AN, B={})
for B in (1, 2, 3, 4):
    st = torch.load(os.path.join(OUT, f"colla_ckpt_an_B{B}.pt"), weights_only=False)
    phi = st["phi"]
    E = float(cc.energy(phi, H, m=AN["m"]))
    U, V, W, c = cc.inertia(phi, H)
    eu, ev = np.linalg.eigvalsh(U), np.linalg.eigvalsh(V)
    sw = np.linalg.svd(W, compute_uv=False)
    r = dict(E=E, U=U.tolist(), V=V.tolist(), W=W.tolist(), U_eig=eu.tolist(), V_eig=ev.tolist(), W_sv=sw.tolist(),
             detW=float(np.linalg.det(W)), centroid=c)
    if B == 1:
        Lam = float(np.mean(np.concatenate([np.diag(U), np.diag(V), -np.diag(W)])))
        r.update(Lam=Lam, Lam_1D=AN["Lam"], Erot=0.75 * hb**2 / (2 * Lam), Erot_Delta=3.75 * hb**2 / (2 * Lam))
    elif B == 2:
        # body frame of V: the torus axis is the eigenvector whose eigenvalue is far from the other two
        w_, Q = np.linalg.eigh(V)
        gaps = [abs(w_[0] - w_[1]) + abs(w_[0] - w_[2]), 0, abs(w_[2] - w_[0]) + abs(w_[2] - w_[1])]
        ia = int(np.argmax(gaps)); perp = [i for i in range(3) if i != ia]
        wu, Qu = np.linalg.eigh(U)
        gu = [abs(wu[0] - wu[1]) + abs(wu[0] - wu[2]), 0, abs(wu[2] - wu[0]) + abs(wu[2] - wu[1])]
        iu = int(np.argmax(gu)); pu = [i for i in range(3) if i != iu]
        Wb = Qu.T @ W @ Q          # iso body frame x spin body frame
        V_perp = float(np.mean(w_[perp])); V_ax = float(w_[ia]); U_perp = float(np.mean(wu[pu])); U_ax = float(wu[iu])
        r.update(V_perp=V_perp, V_axis=V_ax, U_perp=U_perp, U_axis=U_ax, W_body=Wb.tolist(),
                 W_axis=float(Wb[iu, ia]), W_perp_max=float(np.abs(Wb[np.ix_(pu, perp)]).max()),
                 checks=dict(U11_U22=float(abs(wu[pu[0]] - wu[pu[1]]) / U_perp), V11_V22=float(abs(w_[perp[0]] - w_[perp[1]]) / V_perp),
                             V33_over_U33=V_ax / U_ax, W33_over_U33=abs(float(Wb[iu, ia])) / U_ax,
                             W11_over_U11=float(np.abs(Wb[np.ix_(pu, perp)]).max()) / U_perp),
                 Erot=hb**2 / V_perp)
    elif B == 3:
        # Schur: W intertwines the spatial and the (conjugated) iso action of T, so W = c R with R a proper rotation;
        # c = sign(det W) |W| (the hedgehog gives c = -Lambda, det W < 0). In the frame rotated by R the symmetry is
        # 'equal angles' and the allowed J=I=1/2 state has K+L=0 (correction of the first analysis, which assumed c < 0)
        u, v, w = float(eu.mean()), float(ev.mean()), float(np.sign(np.linalg.det(W))) * float(sw.mean())
        r["Erot_assumed_c_negative"] = 0.75 * (u + v - 2 * abs(w)) / (2 * (u * v - w * w)) * hb**2
        r.update(u=u, v=v, w=w, aniso=dict(U=float(np.ptp(eu) / u), V=float(np.ptp(ev) / v), W=float(np.ptp(sw) / abs(w))),
                 Erot=0.75 * (u + v + 2 * w) / (2 * (u * v - w * w)) * hb**2,
                 Erot_alt_M1=None)
        # M=1 (forbidden: the vector rep of T has no 1-dim invariant subspace), shown only for comparison: K.L = 1/4
        r["Erot_alt_M1"] = (v * 0.75 + u * 0.75 - 2 * w * 0.25) / (2 * (u * v - w * w)) * hb**2
    else:
        r.update(Erot=0.0)
    res["B"][B] = r
    print(B, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k not in ("U", "V", "W", "W_body")}, flush=True)

b = res["B"]
E1, Lam = b[1]["E"], b[1]["Lam"]
MN = E1 + b[1]["Erot"]
out = dict(M_N_MeV=MN * Eu, M_Delta_MeV=(E1 + b[1]["Erot_Delta"]) * Eu, E1_MeV=E1 * Eu, Erot_N_MeV=b[1]["Erot"] * Eu)
for B, name in ((2, "d"), (3, "3He"), (4, "4He")):
    cl = (B * E1 - b[B]["E"]) * Eu
    q = (B * MN - b[B]["E"] - b[B]["Erot"]) * Eu
    out[name] = dict(classical=cl, Erot_MeV=b[B]["Erot"] * Eu, quantised=q, mass_MeV=(b[B]["E"] + b[B]["Erot"]) * Eu)
    if B == 3:
        out[name]["quantised_alt_M1"] = (B * MN - b[B]["E"] - b[B]["Erot_alt_M1"]) * Eu
res["binding"] = out
print(json.dumps(out, indent=1))
json.dump(res, open(os.path.join(OUT, "colla_quant.json"), "w"), indent=1)

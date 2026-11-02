"""qg_export_check.py - reads back assets/scene/quark_gluoni/ (as the program would, from the raw bytes) and checks it
against the source data within the quantization error. Output: lab/out/qg_export_check.json.

    python3 lab/qg_export_check.py
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_run import SCR  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DST = os.path.join(REPO, "assets", "scene", "quark_gluoni")
OUT = os.path.join(REPO, "lab", "out", "qg_export_check.json")


def check_quark():
    import qg_render1 as Q
    d = os.path.join(DST, "quark_Y")
    m = json.load(open(os.path.join(d, "manifest.json")))
    raw = np.fromfile(os.path.join(d, "palline.f32"), dtype="<f4")
    assert raw.size == m["n"] * 5, "size"
    P = raw.reshape(m["n"], 5)
    Q.NB = m["n"]
    sc = Q.carica("prod", m["d"], m["T"], seed=m["seme"])
    ref = np.stack([sc["pos"][0], sc["pos"][1], sc["pos"][2], sc["quota"], sc["livello"]], 1)
    err_p = float(np.abs(P - ref).max())
    C = np.fromfile(os.path.join(d, "campo.f32"), dtype="<f4").reshape(m["campo"]["forma"])
    F, fe, _ = Q.dati("prod", m["d"], m["T"])
    err_c = float(np.abs(C - np.stack([F[2], fe, F[0], F[1]])).max())
    tol_p = float(np.abs(ref).max()) * 2 ** -23
    tol_c = float(np.abs(F).max()) * 2 ** -23
    return dict(n=m["n"], byte=os.path.getsize(os.path.join(d, "palline.f32")), max_err_palline=err_p,
                tolleranza_f32=tol_p, max_err_campo=err_c, tolleranza_campo=tol_c,
                range_quota=[float(P[:, 3].min()), float(P[:, 3].max())],
                range_livello=[float(P[:, 4].min()), float(P[:, 4].max())],
                ok=bool(err_p <= tol_p and err_c <= tol_c))


def check_gluoni(run):
    d = os.path.join(DST, f"gluoni_{run}")
    m = json.load(open(os.path.join(d, "manifest.json")))
    N, nf = m["n"], m["n_fotogrammi"]
    raw = np.fromfile(os.path.join(d, "palline.bin"), dtype="u1")
    assert raw.size == nf * 8 * N, "size"
    fs = sorted(glob.glob(os.path.join(SCR, f"qg_ym_{run}_f[0-9]*.npz")))
    shape = np.array(m["forma"], float)
    e_q = e_rb = e_lv = e_pos = 0.0
    e_t = 0.0
    for j, i in enumerate(m["fotogrammi_tenuti"]):
        b = raw[j * 8 * N:(j + 1) * 8 * N]
        q = b[:6 * N].view("<u2").reshape(N, 3).T
        rb = b[6 * N:7 * N].astype(float) / 255
        lv = b[7 * N:].astype(float)
        z = np.load(fs[i])
        e_q = max(e_q, float(np.abs(q.astype(int) - z["q"][:, :N].astype(int)).max()))
        pos = q.astype(float) / 65535 * shape[:, None]
        pos_src = z["q"][:, :N].astype(float) / 65535 * shape[:, None]
        e_pos = max(e_pos, float(np.abs(pos - pos_src).max()))
        e_rb = max(e_rb, float(np.abs(rb - np.clip(z["rb"][:N].astype(float), 0, 1)).max()))
        x = np.clip(z["eb"][:N].astype(float) / m["eps_max"], 0, 1) ** 0.35 * 255
        e_lv = max(e_lv, float(np.abs(lv - x).max()))
        e_t = max(e_t, abs(m["t"][j] - float(z["t"])))
    return dict(n=N, fotogrammi=nf, byte=os.path.getsize(os.path.join(d, "palline.bin")),
                max_err_q_u16=e_q, max_err_rb=e_rb, tolleranza_rb=0.5 / 255, max_err_lv_su_255=e_lv,
                tolleranza_lv=0.5, max_err_t=e_t, G7_min=m["G7_min"], S8_min=m["S8_min"],
                ok=bool(e_q == 0 and e_rb <= 0.5 / 255 + 1e-9 and e_lv <= 0.5 + 1e-9 and e_t == 0))


if __name__ == "__main__":
    out = {"quark_Y": check_quark()}
    for d in sorted(glob.glob(os.path.join(DST, "gluoni_*"))):
        run = os.path.basename(d)[len("gluoni_"):]
        out[f"gluoni_{run}"] = check_gluoni(run)
    out["tutto_ok"] = all(v["ok"] for v in out.values() if isinstance(v, dict))
    tot = sum(os.path.getsize(f) for f in glob.glob(os.path.join(DST, "**", "*"), recursive=True) if os.path.isfile(f))
    out["byte_totali"] = tot
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out, indent=1))

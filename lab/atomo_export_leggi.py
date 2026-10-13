"""Reader for the export of the coupled light + atom (test `accoppiato_esporta`, src/atomo/tests.rs).

The export is a directory with `manifest.json`, raw little-endian float16 files (stored = physical x scale),
`fotogrammi.csv` (scalars per frame) and `serie.csv` (scalars per light tick). See docs/ATOMO.md, "Export".

    python3 lab/atomo_export_leggi.py [DIR] [--npz] [--figura]

    from atomo_export_leggi import carica
    d = carica(DIR)               # memory-mapped, converted lazily
    rho = d.rho(k)                # (48,48,48) float32, |psi|^2 per lambda_C^3, frame k
    j = d.j(k)                    # (3,48,48,48) probability current per lambda_C^3 (units c), averaged over the frame
    E = d.E(k)                    # (3,64,64,64) e*E of the light, m c^2 / lambda_C
    d.fotogrammi["p_exc"]         # scalars per frame (numpy structured array)
"""
import csv
import json
import os
import sys

import numpy as np

DEFAULT = "/mnt/project-files/ripresa/dati/atomo_accoppiato"
QUI = os.path.dirname(os.path.abspath(__file__))


def _csv(path):
    with open(path) as f:
        r = csv.reader(f)
        head = next(r)
        rows = np.array([[float(v) for v in row] for row in r])
    return {h: rows[:, i] for i, h in enumerate(head)}


class Export:
    def __init__(self, cartella):
        self.dir = cartella
        with open(os.path.join(cartella, "manifest.json")) as f:
            self.man = json.load(f)
        self.raw = {}
        for nome, info in self.man["file"].items():
            if not nome.endswith(".f16"):
                continue
            self.raw[nome[:-4]] = (np.memmap(os.path.join(cartella, nome), dtype="<f2", mode="r",
                                             shape=tuple(info["forma"])), float(info["scala"]))
        self.fotogrammi = _csv(os.path.join(cartella, "fotogrammi.csv"))
        self.serie = _csv(os.path.join(cartella, "serie.csv"))
        self.n = self.man["fotogrammi"]

    def _get(self, nome, k):
        a, s = self.raw[nome]
        return np.asarray(a[k], dtype=np.float32) / s

    def rho(self, k): return self._get("rho", k)
    def j(self, k): return self._get("j", k)
    def fase(self, k): return self._get("fase", k)
    def E(self, k): return self._get("E", k)

    def assi(self):
        """Coordinates (lambda_C) of the atom grid and of the (downsampled) light grid."""
        ga, gl = self.man["griglia_atomo"], self.man["griglia_luce_originale"]
        n, a = ga["n"], ga["passo"]
        xa = (np.arange(n) - n / 2 + 0.5) * a
        nl, b = gl["n"], gl["passo"]
        xl = (2 * np.arange(nl // 2) + 0.5 - nl / 2) * b
        return xa, xl

    def a_npz(self, uscita=None):
        """All frames in one compressed .npz (float16, stored values: divide by the scales in `scale`)."""
        uscita = uscita or os.path.join(self.dir, "atomo_accoppiato.npz")
        arr = {k: np.asarray(v[0]) for k, v in self.raw.items()}
        arr["scale"] = np.array([self.raw[k][1] for k in self.raw], dtype=np.float64)
        arr["scale_nomi"] = np.array(list(self.raw))
        for k, v in self.fotogrammi.items():
            arr["fotogrammi_" + k] = v
        for k, v in self.serie.items():
            arr["serie_" + k] = v
        arr["manifest"] = np.array(json.dumps(self.man))
        np.savez_compressed(uscita, **arr)
        return uscita


def carica(cartella=DEFAULT):
    return Export(cartella)


def figura(d, uscite):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fr, se = d.fotogrammi, d.serie
    t = fr["t_dal_sparo"]
    # three moments: before (frame 0), absorption (strongest field at the atom), re-emission (pulse gone:
    # light energy in the box < 1e-3 of the pulse, field at the atom only the atom's own)
    pmax = fr["p_exc"][t < 400].max()
    k_ass = int(np.argmax(fr["p_exc"] > 0.5 * pmax))   # half-way up the absorption rise, pulse on the atom
    e_box = fr["energia_luce_box"]
    dopo = np.where(e_box < 1e-3 * e_box[0])[0]
    k_rie = int(dopo[np.argmin(np.abs(t[dopo] - 800.0))]) if len(dopo) else d.n - 1
    scelti = [(0, "sparo: l'impulso (lungo +x) sta arrivando"), (k_ass, "assorbimento"), (k_rie, "riemissione (impulso uscito)")]

    xa, xl = d.assi()
    na, nl = len(xa), len(xl)
    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.15, 1.15], hspace=0.38, wspace=0.25)
    ax = fig.add_subplot(gs[0, :])
    ax.plot(se["t_dal_sparo"], se["p_exc"], color="#c0392b", lw=1.2, label="p_exc = 1 − |⟨1s|ψ⟩|² (ogni tick di luce)")
    ax.plot(t, fr["p2p"] * 50, color="#7f8c8d", lw=1, label="p(2p½, m=½) × 50")
    ax2 = ax.twinx()
    ax2.plot(t, fr["eE_atomo_z"] * 1e3, color="#2e86c1", lw=0.8, alpha=0.8)
    ax2.set_ylabel("e·E_z sull'atomo (10⁻³ mc²/λ_C)", color="#2e86c1")
    for k, nome in scelti:
        ax.axvline(t[k], color="k", ls=":", lw=1)
        ax.text(t[k], 1.02, nome.split(":")[0], rotation=0, fontsize=8, ha="left", va="bottom")
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0, 1)
    ax.set_xlabel("t dallo sparo (λ_C/c)")
    ax.set_ylabel("popolazione")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"Luce + atomo accoppiati (Z = {d.man['Z']:.0f}, K = {d.man['K_ACC']:.0f}): {d.n} fotogrammi ogni "
                 f"{d.man['dt_fotogramma']} λ_C/c; ω21 = {d.man['omega_21']:.4f}", fontsize=10, pad=16)

    ext_a = [xa[0], xa[-1], xa[0], xa[-1]]
    ext_l = [xl[0], xl[-1], xl[0], xl[-1]]
    ia, il = na // 2, nl // 2   # y = 0 lies between the two middle cells (atom) / near them (light: I = 31.75)
    for c, (k, nome) in enumerate(scelti):
        rho = d.rho(k)
        sl = 0.5 * (rho[:, ia - 1, :] + rho[:, ia, :])   # y = 0 plane (between the two middle cells)
        a1 = fig.add_subplot(gs[1, c])
        a1.imshow(np.sqrt(sl.T), origin="lower", extent=ext_a, cmap="magma")
        a1.set_title(f"{nome}\nt = {t[k]:.0f}, p_exc = {fr['p_exc'][k]:.2f}, ⟨z⟩ = {fr['dipolo_z'][k]:+.2f}", fontsize=9)
        a1.set_xlabel("x (λ_C)"); a1.set_ylabel("z (λ_C)")
        a1.text(0.02, 0.95, "√ρ, piano y=0", transform=a1.transAxes, color="w", fontsize=8, va="top")
        ez = d.E(k)[2]
        sz = 0.5 * (ez[:, il - 1, :] + ez[:, il, :])
        a2 = fig.add_subplot(gs[2, c])
        # colour range: the pulse in the first two panels; the atom's own (much weaker) light in the third
        lim = np.percentile(np.abs(sz), 99.5) if c < 2 else np.percentile(np.abs(sz[np.hypot(*np.meshgrid(xl, xl, indexing="ij")) > 40]), 99.5)
        a2.imshow(sz.T, origin="lower", extent=ext_l, cmap="RdBu_r", vmin=-lim, vmax=lim)
        a2.add_patch(plt.Rectangle((xa[0], xa[0]), xa[-1] - xa[0], xa[-1] - xa[0], fill=False, ec="k", lw=0.6))
        a2.set_title(f"e·E_z della luce, piano y=0 (scala ±{lim:.1e})", fontsize=9)
        a2.set_xlabel("x (λ_C)"); a2.set_ylabel("z (λ_C)")
    for u in uscite:
        os.makedirs(os.path.dirname(u), exist_ok=True)
        fig.savefig(u, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return scelti


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cart = args[0] if args else DEFAULT
    d = carica(cart)
    m = d.man
    print(f"{cart}: {d.n} frames, every {m['tick_per_fotogramma']} atom ticks ({m['dt_fotogramma']} λ_C/c)")
    for k, (a, s) in d.raw.items():
        print(f"  {k:5s} {a.shape} float16, scale {s:g}")
    fr = d.fotogrammi
    kp = int(np.argmax(d.serie["p_exc"]))
    print(f"  p_exc max {d.serie['p_exc'][kp]:.3f} at t = {d.serie['t_dal_sparo'][kp]:.0f}; "
          f"norm check sum(rho)*a^3 frame 0 = {d.rho(0).sum() * m['tick_atomo'] ** 3:.4f}, last = {d.rho(d.n - 1).sum() * m['tick_atomo'] ** 3:.4f}")
    if "--npz" in sys.argv:
        print("npz:", d.a_npz())
    if "--figura" in sys.argv:
        sc = figura(d, [os.path.join(cart, "controllo.png"), os.path.join(QUI, "out", "atomo_export_controllo.png")])
        print("figure frames:", [(k, round(float(fr["t_dal_sparo"][k]), 1), n) for k, n in sc])

"""Figura dello studio SOTTORETICOLO (docs/SOTTORETICOLO.md): dai risultati di sottoreticolo_esame."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("lab/out/sottoreticolo_results.json"))
s = d["serie_E2_per_tick_luce"]
t = [8 * 0.5 * (i + 1) for i in range(len(s["p_exc_intero"]))]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
ax[0].plot(t, s["p_exc_intero"], label="cubo intero (Lab di prima)", lw=1.6)
ax[0].plot(t, s["p_exc_colore"], "--", label="un colore (00)", lw=1.2)
ax[0].set(xlabel="t (λ_C/c)", ylabel="p_exc", title="E2: p_exc, cubo intero vs un colore")
ax[0].legend()
diff = [a - b for a, b in zip(s["p_exc_colore"], s["p_exc_intero"])]
ax[1].plot(t, diff, color="C3")
ax[1].set(xlabel="t (λ_C/c)", ylabel="Δ p_exc", title=f"differenza (max {max(map(abs, diff)):.1e}): un colore è anisotropo")
S1, S2, S3 = d["S1_avanti"]["rapporto"], d["S2_passo_accoppiato"]["rapporto"], d["S3_prepara_s"]["rapporto"]
ax[2].bar(["regola\n(avanti)", "passo accoppiato\n(+ luce 128³)", "preparazione\n(1s esatto)"], [S1, S2, S3], color=["C2", "C0", "C1"])
for i, v in enumerate([S1, S2, S3]):
    ax[2].text(i, v + 0.05, f"×{v:.2f}", ha="center")
ax[2].axhline(4, color="gray", ls=":")
ax[2].set(ylabel="velocità: intero / un colore", title="guadagno nel cloud (2 thread)")
e = d["esiti"]
fig.suptitle("SOTTORETICOLO — " + ", ".join(f"{k} {v}" for k, v in e.items()), fontsize=10)
fig.tight_layout()
fig.savefig("lab/out/sottoreticolo_figure.png", dpi=110)
print("ok")

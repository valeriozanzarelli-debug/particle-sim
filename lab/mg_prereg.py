"""mg_prereg.py - pre-registration of the study MOTORE GRAFICO A PALLINE (rendering engine of the Laboratorio
Atomo: GPU instanced balls, rule on GPU, quality presets up to Ultra, display levels of detail).
Written and committed BEFORE any production measurement. Writes lab/out/mg_prereg.json.
Physics multi-level (which agglomerate level a region uses, local time step) belongs to docs/MOTORE_OTTIMIZZATO.md:
here only the DRAWING levels."""
import json, os

PRE = {
    "studio": "motore_grafico",
    "data": "2026-10-13",
    "hardware_di_prova": {
        "cloud": "CPU 4 core + Vulkan software (llvmpipe) sotto Xvfb: correttezza e foto, non tempi",
        "portatile": "Intel UHD (la GTX 1650 Ti e' rotta): tempi veri via ponte (foto con --misure del tempo)",
        "obiettivo": "GPU da 16 GB (Ultra): solo stime, nessun esame finche' non c'e' la scheda",
    },
    "preset": {
        "bassa":  {"celle_per_lambda": 12, "n": 96,  "palline": 20000,   "scia": 2, "regola": "cpu"},
        "media":  {"celle_per_lambda": 16, "n": 128, "palline": 60000,   "scia": 3, "regola": "cpu"},
        "alta":   {"celle_per_lambda": 24, "n": 192, "palline": 250000,  "scia": 3, "regola": "gpu"},
        "ultra":  {"celle_per_lambda": 32, "n": 384, "palline": 4000000, "scia": 4, "regola": "gpu"},
    },
    "esami": {
        "R1_regola_gpu": "la regola della luce in WGSL contro fotone3d.rs (CPU, f32) sullo stesso impulso, 96 tick, "
                         "n = 96: max |F_gpu - F_cpu| / max |F_cpu| <= 1e-4 (arrotondamenti f32 diversi)",
        "R2_istanze": "disegno a istanze (un buffer di 32 byte per pallina, niente mesh ricostruita) contro il "
                      "disegno a mesh di oggi, stesso fotogramma deterministico: PSNR >= 35 dB",
        "R3_livelli_disegno": "palline aggregate per blocchi (livello k: blocchi 2^k celle) per la vista da lontano: "
                              "energia totale disegnata conservata entro 1e-6, baricentro entro 0.05 celle, spin "
                              "letto dalle palline aggregate (come FOTONE_PALLINE F2) entro 1% del livello 0, per k = 1, 2, 3",
        "R4_tempi_intel": "sul portatile (Intel UHD) al preset bassa: >= 20 fotogrammi/s e <= 50 ms per tick nella "
                          "scena fotone; media: solo riportato",
        "R5_niente_geometria": "ogni preset disegna solo palline (lettura del codice): nessuna linea, tubo, superficie",
    },
    "debiti_noti": [
        "Ultra non provato su una scheda vera (nessuna GPU da 16 GB disponibile)",
        "i tempi del cloud (llvmpipe) non dicono niente sulla GPU",
    ],
}

if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "mg_prereg.json")
    json.dump(PRE, open(out, "w"), indent=1, ensure_ascii=False)
    print("scritto", out)

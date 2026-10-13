# Problemi aperti attaccati con la regola zig-zag (notte 13/10/2026)

Cinque studi, uno per problema aperto di fisica atomica/quantistica. Ogni studio ha fatto ricerca,
pre-registrazione prima delle corse, corse con convergenza in a → 0, documento con debiti.
Parametri liberi della regola: 0 in tutti. Ramo: `claude/problemi-aperti-56mnej` (non su master).

| Studio | Documento | Validazione | Nuovo (piccolo, controllato) | Limite trovato |
|---|---|---|---|---|
| Riflusso quantistico | `RIFLUSSO.md` | c_bm 0,0384498 vs 0,0384517; Dirac vs Ashfaque 2019 entro 1e-3 | limite Δ ≤ mc²T/(4πħ) (Cauchy–Schwarz), saturato nell'ultrarelativistico (0,99998 a ε=1e4); piccola correzione alla tabella di Ashfaque a ε=0,5 (+0,24 %) | la regola a passo finito non supera mai c_bm (congettura); Q1 FAIL; stati legati non ben posti |
| Tempo di tunnel | `TEMPO_TUNNEL.md` | Hartman = Dirac a 4e-6; Larmor = dwell a 2e-7; Büttiker a 3e-5 | contare le rotazioni dentro la barriera = orologio di Larmor sulla massa = tempo proprio; vale L/c esatto a E = V0 (Büttiker–Landauer dà 0), va all'indietro per V0 > E | tempo di arrivo dipende dal rivelatore (E13 FAIL); Ramos 2020 in ms non confrontabile |
| Atomo supercritico | `SUPERCRITICO.md` | Z_cr 177,463 vs 177,46 (Godunov 2017); sfera 169,706; Γ entro 0,2–1,8 % | primo film del tuffo dell'1s e del positrone che esce fatto con un walk di Dirac (nuovo come metodo, non come fisica) | nucleo puntiforme: il reticolo fa tuffare l'1s a Z* = 131,37 < 137 per ogni passo a → servono nuclei estesi |
| Schwinger assistito | `SCHWINGER_ASSISTITO.md` | Schwinger 0,1 %; Keldysh 0,998–1,0001; γ_crit 1,606 vs 1,592 | nessuna scoperta (QED 1+1 omogenea è risolta); numeri esatti del modello dove WKB sbaglia di ×2–3e4 | 1+1D; laser reali lontani (XFEL 10 keV → e^−321); A1 FAIL a ε<0 per cancellazioni. Coppie contate solo nella valle vera (lontano da k = π) |
| Reticolo vs Lorentz | `LIMITI_LORENTZ.md` | dispersione di ZIGZAG.md, formule chiuse a 1e-15 | niente termine lineare (lo uccide il tick palindromico), niente birifrangenza a tutti gli ordini; a ≤ 4,8e-31 m ≈ 3e4 l_P dai fotoni PeV di LHAASO | reticolo del Laboratorio escluso di ~17 ordini: il reticolo del calcolo non può essere quello fisico; rischio radiativo (Collins 2004) non calcolato |

Figure: `lab/out/rifl_figura.png`, `tunnel_figura.png`, `supercrit_figura.png`, `assist_figura.png`, `lv_figure.png`.

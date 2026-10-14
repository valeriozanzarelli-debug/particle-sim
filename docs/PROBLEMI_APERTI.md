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

## Giri 2 e 3 (13/10/2026)

| Studio | Documento | Validazione | Nuovo (piccolo, controllato) | Limite trovato |
|---|---|---|---|---|
| Protocollo per il tempo proprio nel tunnel | `TUNNEL_PROTOCOLLO.md` | lettura con ancella = previsione a 3e-6; Dirac continuo a 6,5e-4 | formule esatte del reticolo (pendenza 1/cos θ, rapporto tan(εa)/tan θ); cammino in fibra: ~1200 trasmessi per punto per 5σ | l'idea di orologio di massa non è nostra (Schach–Giese 2024); il rumore (≤ 0,003 rad/passo) conta più della statistica; 7 FAIL su 21 |
| Causa del tuffo finto a Z* = 131,4 | `ZSTAR.md` | ritrova β(Z=60) = 0,10301 di nucleo1s_ | formula senza parametri Z* = 131,3614 (regola 131,3612); tuffi finti previsti e confermati anche per p3/2 (206,0) e d3/2 (42,05) | Z* è un artefatto della condizione al centro (varianti da 95 a 130); rimedio: celle a (i+1)a o (i+√3/2)a, poi Z_cr → 137 solo logaritmicamente; film Z=200 esportato in `lab/out/zstar_film.json` |
| Tempi di arrivo spin 1/2 | `TEMPI_ARRIVO.md` | flusso = esatto; ABR di Jozani–Tumulka 0,965 vs 96 % | in 2-D il termine di spin nasce da solo (taglio 7,215 vs 7,223); la cella-rivelatore della regola = ABCD ideale di Tumulka (κ = 2mc/ħ) | nessun rivelatore lineare e simmetrico vede lo spin; τ_max di Das–Dürr −3,2 % (E1 FAIL); niente 3-D cilindrico per doppioni e sottoreticoli |
| Dilatazione quantistica del tempo | `DILATAZIONE_QUANTISTICA.md` | Smith–Ahmadi entro 1,8 % (lento), continuo esatto 4e-4 | ritmo dell'orologio √(1−v²)·cos(ka); correzione nulla a p = mc/√2 | per particella libera la regola non può differire dalla letteratura; nessuna scoperta |
| Violazione di Lorentz dai loop | `LORENTZ_LOOP.md` | quadro di Collins et al. 2004 | simmetria vera della regola = D4h con asse y (16/48), ordine a giro = T_h | con la fase 3 di Debiti: regola congelata esclusa come reticolo dell'universo per qualunque a (anisotropia 1e-4 vs 1e-21); ordine a giro escluso con riserva (gap +0,044 α vs luce nuda) |
| Velocità protetta (problema di Collins) | `VELOCITA_PROTETTA.md` | ω = 2ω_W(k/2) di Bisio et al. | la luce della regola è esattamente una coppia di walker di Weyl a k/2 (1,3e-15); c vestita dell'elettrone indipendente dall'impulso | nessun meccanismo noto protegge c_e = c_γ nella regola; il fotone composto eredita c solo al livello delle autoenergie (gap ≤ 2,5e-4 α, ne serve 1e-18) |

## Giro 4: doppioni e FW1 (13/10/2026)

| Studio | Documento | Risultato |
|---|---|---|
| Regola senza doppioni | `SENZA_DOPPIONI.md` | teorema: nessuna regola fatta di pezzi a una direzione evita i gemelli (W₃ = 0, anche BCC/FCC); proposta FW1 = congelata + un Wilson ritorto (h = 0,4, debito) sui tick pari: gemelli ≥ 0,37/a, coppie 0,99 di una specie senza proiettore; T rotta |
| FW1, vuoto senza proiettore | `FW1_VUOTO.md` | mare intero: 2α/3π a −0,14 %, Uehling −0,2/−0,3 %; congelata a mare intero: 4 specie meno 4 coni (validazione del problema); luce vestita +0,043 α e birifrangente (Collins non curato); FW1 rompe x↔z |
| FW1, due righe, battimento, atomo 3-D | `FW1_ATOMO.md` | le due righe a θ = 0,5 sono la coppia di spin spaccata dalla T rotta (∝ θ^3,3); una sola riga per θ ≤ 1/3; battimento congelata ritrovato (0,84037/0,88627/0,89930/0,89960); atomo 3-D senza toppe: FW1 un solo 1s al ≥ 99,5 % nella valle vera, errore tra −2,5 % e +2,7 % (scatola); costo ×2,1–2,4 |

## Pausa (14/10/2026, 01:00 UTC) — dove si riprende

Progetto in pausa per decisione di Valerio. Nessun agente attivo. Tutto è sul ramo `claude/problemi-aperti-56mnej` (non su master: manca l'ok scritto di Valerio).

FW1 è fermo qui: è una variante IN PROVA accanto alla congelata (che resta il riferimento); h = 0,4 è un debito. Funziona per θ ≤ 1/3 (vuoto senza proiettore, un solo 1s nella valle vera). Al passo del Laboratorio (θ = 0,5) l'1s è una coppia di spin spaccata dalla T rotta (∝ θ^3,3).

Prossimi passi, in ordine:
1. Una variante di FW1 che conservi T esatta (es. Wilson ritorto simmetrico sui due tick, o fase che si compensa nel periodo), con lo stesso test: coppie senza proiettore, 1s a θ = 0,5 una riga sola.
2. c dell'elettrone vestito con FW1 (manca): poi c_e − c_luce per FW1 contro la congelata.
3. Atomo 3-D con FW1: scatola più grande per portare l'errore sotto l'1 % (oggi tra −2,5 % e +2,7 %).
4. Portare su master (serve l'ok di Valerio) e passare FW1 al Laboratorio Atomo (thread Ripresa) e al renderer (Modello standard).

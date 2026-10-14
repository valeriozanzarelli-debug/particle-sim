# Tavola periodica, parte 2: dal boro in su sulla griglia agglomerata — PARZIALE

Studio `tavg_`. Stato al 14/10/2026: **PARZIALE** (fermato su richiesta, a un checkpoint). Pre-registrazione della
fase 1 committata prima delle corse: `lab/out/tavg_prereg.json` (commit 69414b4). Fase 2 (correlazione) non iniziata.
Tier di questi risultati: **"approssimato"** (Hartree–Fock + blocchi imparati; non esatto come H, He, H₂, Li, Be).

## Cosa si fa, detto semplice

L'ECG esatto non va oltre 4–5 elettroni. Per B…Kr la strada che scala è a **tre livelli dell'agglomeratore**:

1. **Livello fine = zoom radiale** (`lab/tavg_radial.py`): l'atomo intero, tutti gli elettroni, con l'operatore della
   regola (p²/2m_kin − Z/r + 1/r_ij), Hartree–Fock ristretto del termine di Hund. Nessun orbitale da tabella: si parte
   dagli orbitali del nucleo nudo e il campo medio si costruisce da solo.
2. **Blocco-nucleo** (`lab/tavg_core.py`, `lab/tavg_blocks.py`): nucleo + gusci interni diventano UNA pallina. La
   sua regola (come spinge gli elettroni esterni, un potenziale per ogni l) si impara invertendo l'equazione del livello
   fine; la grandezza del blocco r_c la sceglie la macchina (il più piccolo blocco che la griglia grossa riesce a
   portare), e la cella h la sceglie la macchina (la più grossa fra a₀/4, a₀/6, a₀/8 che trasferisce il catione entro
   0,8 mEh). In linguaggio standard è uno pseudopotenziale a conservazione di norma generato da noi: lo dichiaro.
3. **Griglia 3-D** (`lab/tavg_grid.py`): gli elettroni esterni su celle cubiche, con i blocchi in posizioni qualsiasi.

## Numeri che ci sono (produzione, contro la pre-registrazione)

| esame | cosa | risultato | esito |
|---|---|---|---|
| F1 (B, C) | ROHF fine contro limite HF | B −24,5290607 (chiave −24,529061, +3·10⁻⁷); C −37,6886190 (chiave −37,688619, ~0) | **PASS** (2 su 6) |
| F2 (B⁺, C⁺) | catione | B⁺ −24,2375752 (−24,237575); C⁺ −37,2922238 (−37,292224) | **PASS** (2 su 6) |
| F4 (B, C) | aufbau: la configurazione più bassa fra i candidati è quella vera | B 2s²2p, B⁺ 2s², C 2s²2p², C⁺ 2s²2p | **PASS** (4 casi su 40) |
| IP HF (fine) | ΔSCF | B 7,932 eV, C 10,786 eV (sperimentali 8,298 e 11,260: manca la correlazione, fase 2) | solo riportato |
| REL | stima relativistica su IP | B −0,0023 eV | trascurabile |
| B2 (B) | blocco → catione | +0,37 mEh a a₀/4 | riportato |
| B3 (B) | promozione 2s→2p | −1,51 mEh (soglia 1,6) | **PASS** |
| G1–G7 | griglia 3-D | **non ancora misurati in produzione** | — |

Sviluppo (dichiarato, prima della pre-registrazione): livello fine Ne +4·10⁻¹⁰ Eh, Kr +9·10⁻⁹, Fe +2·10⁻⁶ dai limiti HF;
griglia 3-D contro il radiale con lo stesso blocco: B +0,001 mEh, C −0,11, Ne +0,08, O −3,0 mEh (a a₀/4 e a₀/6: è
rottura di simmetria vera, il p doppiamente occupato prende un raggio diverso); una corsa di griglia dura 5–20 s.
Ne ha bisogno della cella a₀/8 (catione: 24 mEh a a₀/4, 4,1 a a₀/6, 0,74 a a₀/8). C: la promozione 2s→2p sbaglia di
~4,5 mEh qualunque r_c (integrali valenza–valenza dentro il blocco): FAIL atteso per B3 su C, N, O.
Lo stencil imparato a₀/4 di aggmol vale solo fino a k ≈ 6/a₀ (−2,4 % a 7,9, −8 % a 9,4): per C…Ne la griglia usa la
banda della regola k²/2m spettrale (debito: ri-imparare il livello con l'osservabile dei gusci 2p).

## Funzione per gli altri studi (molecole)

`lab/tavg_grid.py`: `hf_griglia(nuclei=[(Z, (x, y, z)), ...], carica=0, spin=None, h=None, L=16.0, restricted=True)`
→ dict con `E` (totale, Eh), `E_val`, `orbitali`, `occ`, `eps`, `densita`, `griglia`. h = None prende la cella più
fine scelta dai blocchi presenti. Richiede i blocchi in `lab/out/tavg_blocchi/Z??.json` (oggi solo **Z05**, boro;
gli altri si rifanno con `./tavg_drive.sh tavg_blocks_prod.py <Z...>`; H, Li, Be, C, N, O sono nella lista).
Collaudo pre-registrato G7 (H₂, N₂) non ancora corso.

## Debiti

Forma TM del blocco e coda a carica puntiforme oltre r_c + 2,75 a₀ (manopole: bilancio di coda 2·10⁻⁴ Eh, 0,8 π/h,
0,8 mEh per scegliere h, canale locale p); nucleo congelato; Born–Oppenheimer; determinanti a orbitali reali (per d², d³,
d⁷, d⁸ non sono termini puri); cinetico della regola spettrale invece dello stencil imparato oltre k ≈ 6; relatività
solo stimata. Parametri liberi di fisica: zero (ingressi Z e costanti).

## Come riprendere

```
cd lab
./tavg_drive.sh tavg_fine.py 5 36                       # livello fine (riprende dal checkpoint, B e C fatti)
./tavg_drive.sh tavg_blocks_prod.py 6 7 8 9 10 1 3 4 11 ... 36
./tavg_drive.sh tavg_run_grid.py atom 5 6 7 8 9 10       # poi conv, learned, shift, mol H2 N2
```
Poi: verdetti G1–G7, export in `/mnt/project-files/debiti/dati/<nome>/` (tier "approssimato"), figura, fase 2
(correlazione imparata dagli esatti He-like/Li-like/Be-like della regola, `lab/tavg_ecg.py` copia di multi_e_ecg).

## Dove l'intuizione di Valerio può aiutare / Idee geometriche e visive

- **Il blocco-nucleo come pallina visibile**: a schermo il nucleo con i gusci interni è una sfera morbida di raggio r_c
  (0,8–1,5 a₀) che "sa" respingere gli elettroni esterni; zoomando dentro, la sfera si apre nei gusci veri del livello
  fine. Costo zero (i dati ci sono), rispetta le regole (è un livello dell'agglomeratore). Domanda: ti va che lo zoom
  sull'atomo mostri proprio questo passaggio blocco → gusci?
- **La rottura di simmetria dell'ossigeno**: il p doppiamente occupato è più largo (⟨r⟩ 1,275 contro 1,204 a₀): a
  schermo una nube leggermente schiacciata lungo un asse. È fisica vera di HF, non un difetto.

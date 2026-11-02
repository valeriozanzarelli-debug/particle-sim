# Tavola periodica, parte 2: dal boro in su sulla griglia agglomerata — PARZIALE (ossigeno: HF sì, correlazione no)

Studio `tavg_`. Stato al 2/11/2026: **PARZIALE**. Pre-registrazione della fase 1: `lab/out/tavg_prereg.json` (commit
69414b4). Emendamento fase 2 (correlazione, bersaglio ossigeno): `lab/out/tavg_prereg2.json` (commit ec4a4ba, prima
di ogni fit). Verdetti di O (con N, C, H): `lab/out/tavg_verdetti_ossigeno.json` (`python lab/tavg_verdetti_o.py`).

## Ossigeno, in breve (2/11)

- **Hartree–Fock dell'ossigeno funziona sui tre livelli.** Fine: O −74,8093985 (chiave −74,809398, 5·10⁻⁷), O⁺
  −74,3726057 (3·10⁻⁷). Griglia 3-D a a₀/4 col blocco imparato (r_c 1,2/1,175 a₀): O −74,812490, O⁺ −74,372471.
- **La rottura di simmetria di −3,1 mEh è fisica di HF, non la griglia** (prova `lab/tavg_simm.py`): con i tre p
  costretti ad avere la stessa forma radiale (quella del livello fine) la griglia ridà il radiale entro +0,11 mEh; se li
  lascio liberi scende di −3,09 mEh a a₀/4 e −3,05 a a₀/6, uguale. Il p doppiamente occupato si allarga (⟨r⟩ 1,276
  contro 1,204 a₀). È HF ristretto con orbitali non equivalenti (noto in letteratura), non UHF (che dà −10,4 mEh).
- Costo: il −3,1 mEh sta solo nel neutro, così l'IP HF della griglia (11,974 eV) supera quello del livello fine
  (11,886) di 0,088 eV: **G4 FAIL** (soglia 0,05). È la stessa fisica, ma la soglia è quella e resta.
- **La correlazione imparata NON regge sull'ossigeno.** Imparata solo dagli esatti della regola (H⁻, He, He 2³S, Li⁺,
  Li, Li 2P, Be⁺, ioni tipo He Z = 6–10; 12 punti), forma LYP a 4 manopole (scelta da sola per validazione incrociata:
  1,6 mEh contro 2,4 della forma a 2), prevista alla cieca: E_c(O) = −172 mEh contro −258 (**−33 %, C1 FAIL**);
  IP(O) = 13,204 eV contro 13,618 (**−3,04 %, C2 FAIL** per un soffio; con l'HF del livello fine 13,116, −3,7 %).
  Va invece bene la riga: IP(N) 14,139 contro 14,534 (−2,7 %, C5 PASS), caduta N→O 0,935 eV contro 0,916 (C3, C4 PASS).
- **Perché**: nell'allenamento ci sono solo coppie del guscio interno (1s², 1s2s, 1s2p) e lì E_c è quasi costante
  (42–47 mEh). Nell'ossigeno ~¾ della correlazione sta nelle coppie del guscio 2 (2s², 2s2p, 2p²), che nessun dato
  ha mai mostrato. Prova (post-hoc, non contata): le costanti LYP pubblicate (prese da He soltanto) danno −265 mEh
  su O ma sbagliano i NOSTRI esatti (Be⁺ −61 invece di −47): la forma a densità non distingue i tipi di coppia, il
  risultato giusto su O della letteratura è in parte fortuna. Togliere gli ioni He-like o cambiare i pesi non cambia
  niente (−171/−174 mEh).
- **Export di O: NON fatto** (G4 e C1, C2 falliscono). Il blocco O (`lab/out/tavg_blocchi/Z08.json`, HF) è usabile
  con `hf_griglia` con errore dichiarato (−3,1 mEh sul totale, +0,09 eV sull'IP HF).
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
| F1, F2 (N, O) | fine contro chiave | N −2·10⁻⁷, N⁺ −1·10⁻⁸, O −5·10⁻⁷, O⁺ +3·10⁻⁷ | **PASS** |
| F3, F4 (N, O) | gusci, Kato, aufbau | 2 massimi; Kato 6,9993 / 7,9992; aufbau giusto per N, N⁺, O, O⁺ | **PASS** |
| REL (N, O) | stima relativistica su IP | N −0,007 eV, O −0,010 eV | sotto la bandiera |
| B3 (N, O) | promozione 2s→2p | N −3,05, O −3,40 mEh | **FAIL** (atteso) |
| G1 (N, O) | griglia contro radiale, stesso blocco | N −0,02, N⁺ −0,08; O −3,09 (finestra non sferica), O⁺ +0,07 mEh | **PASS** |
| G2 (O) | a₀/4 contro a₀/6, L = 12 | 0,06 mEh | **PASS** |
| G3 (N, O) | totale griglia contro limite HF | N −0,02, N⁺ +0,24; O −3,09, O⁺ +0,14 mEh | **PASS** |
| G4 (N, O) | IP HF griglia contro fine | N +0,007 eV; O +0,088 eV | N **PASS**, O **FAIL** |
| G5 (O) | stencil imparato a₀/4 | −5,5 mEh | **FAIL** (atteso) |
| G6 (C, O) | uova: atomo spostato di (0,37, 0,21, 0,11) celle | C 0,098 mEh; O 0,165 mEh (a a₀/6: 6·10⁻⁷) | C **PASS**, O **FAIL** |
| G7 | molecole: H₂, N₂ | H₂ −0,014 mEh; N₂ +3,62 mEh (soglia 3) | H₂ **PASS**, N₂ **FAIL** |
| C1–C5 | correlazione imparata, cieca | vedi sopra | C1, C2 **FAIL**; C3, C4, C5 **PASS** |

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
oggi Z01 H, Z05 B, Z06 C, Z07 N, Z08 O; gli altri con `./tavg_drive.sh tavg_blocks_prod.py <Z...>`).
**Bug trovato e corretto il 2/11**: in `_build_ext` il nome `Y` (griglia) veniva riusato per le armoniche sferiche:
con due nuclei e proiettori KB (N₂, e quindi H₂O) la funzione si rompeva. Un atomo solo e H₂ non erano toccati.
Per l'acqua: la O–H è a 1,81 a₀ e il blocco O ha coda fino a 1,2 + 2,75 a₀; N₂ (2,07 a₀) sbaglia di +3,6 mEh: aspettarsi
un errore di qualche mEh da sovrapposizione dei blocchi finché non si misura su H₂O stessa.

## Debiti

Forma TM del blocco e coda a carica puntiforme oltre r_c + 2,75 a₀ (manopole: bilancio di coda 2·10⁻⁴ Eh, 0,8 π/h,
0,8 mEh per scegliere h, canale locale p); nucleo congelato; Born–Oppenheimer; determinanti a orbitali reali (per d², d³,
d⁷, d⁸ non sono termini puri); cinetico della regola spettrale invece dello stencil imparato oltre k ≈ 6; relatività
solo stimata. Parametri liberi di fisica in HF: zero (ingressi Z e costanti). **Correlazione**: forma LYP presa dalla
letteratura (debito), 4 manopole fittate sui nostri 12 esatti (a = 0,0248, ab = 6,9·10⁻⁴, c = 0,30, d = 0,05 — d al
bordo della griglia pre-registrata: i dati non lo fissano). Esclusi dall'allenamento: Be (studio non chiuso), H₂ e LiH
(manca la loro densità HF sullo stesso livello).

## Come stringere la regola (proposta, invece di abbassare le soglie)

Il problema non è la forma: è che l'agglomeratore non ha mai visto una coppia del guscio 2. Regola proposta,
"correlazione per coppia, imparata nel campo vero":
1. **Coppie di valenza esatte sul blocco**: due elettroni nel campo del blocco O già imparato (1s² congelato) —
   cioè O⁶⁺ con la valenza 2s², 2s2p (³P, ¹P), 2p² (³P, ¹D, ¹S) — con l'ECG della regola (2 elettroni: costo piccolo).
   Pezzo nuovo da scrivere: il potenziale separabile del blocco dentro l'ECG. Uguale per C⁴⁺ e N⁵⁺ per vedere come la
   coppia cambia con la densità.
2. **Be quando chiude** (2s² con la quasi-degenerazione 2s/2p, ~45 mEh): primo esatto di una coppia del guscio 2 a
   densità "atomo neutro".
3. **Legge per tipo di coppia**: E_c = Σ_coppie ε_tipo(densità della coppia), con le manopole contate per tipo
   (s², sp, p² singoletto/tripletto, interguscio), poi un controllo a 3 elettroni (O⁵⁺ 2s²2p) per vedere se le coppie
   si sommano. Solo dopo, di nuovo O alla cieca con le stesse soglie (C1 5 %, C2 3 %).
Diagnosi gratis per la regola di h: anche le "uova" (G6) devono stare nel bilancio; O a a₀/6 passa (6·10⁻⁷).

## Come riprendere

```
cd lab
./tavg_drive.sh tavg_fine.py 5 36                       # livello fine (riprende dal checkpoint, B e C fatti)
./tavg_drive.sh tavg_blocks_prod.py 6 7 8 9 10 1 3 4 11 ... 36
./tavg_drive.sh tavg_run_grid.py atom 5 6 7 8 9 10       # poi conv, learned, shift, mol H2 N2
./tavg_drive.sh tavg_simm.py                             # prova della rottura di simmetria di O
python tavg_corr.py train && python tavg_corr.py predict # fase 2 (dopo tavg_prereg2.json)
python tavg_verdetti_o.py && python tavg_figure.py       # verdetti e figura lab/out/tavg_figure.png
```
Poi: verdetti G1–G7, export in `/mnt/project-files/debiti/dati/<nome>/` (tier "approssimato"), figura, fase 2
(correlazione imparata dagli esatti He-like/Li-like/Be-like della regola, `lab/tavg_ecg.py` copia di multi_e_ecg).

## Dove l'intuizione di Valerio può aiutare / Idee geometriche e visive

- **Il blocco-nucleo come pallina visibile**: a schermo il nucleo con i gusci interni è una sfera morbida di raggio r_c
  (0,8–1,5 a₀) che "sa" respingere gli elettroni esterni; zoomando dentro, la sfera si apre nei gusci veri del livello
  fine. Costo zero (i dati ci sono), rispetta le regole (è un livello dell'agglomeratore). Domanda: ti va che lo zoom
  sull'atomo mostri proprio questo passaggio blocco → gusci?
- **Le coppie come palline doppie**: la proposta sopra è vicina all'idea dei livelli: due elettroni che si evitano
  sono un "mini-agglomerato" con la sua regola imparata. Domanda per te: ti torna che ogni coppia di elettroni sia un
  oggetto con la sua correzione, come una molla che si impara una volta e si riusa nell'acqua?
- **La rottura di simmetria dell'ossigeno**: il p doppiamente occupato è più largo (⟨r⟩ 1,275 contro 1,204 a₀): a
  schermo una nube leggermente schiacciata lungo un asse. È fisica vera di HF, non un difetto.

Figura: `lab/out/tavg_figure.png` (a: correlazione imparata, allenamento e cieco; b: IP di N e O; c: la prova della
rottura di simmetria).
